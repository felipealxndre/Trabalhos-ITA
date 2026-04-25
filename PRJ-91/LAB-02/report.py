"""
Saida: tabelas (console + CSV) e graficos de potencia e polar de velocidades.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from config import (
    Aircraft, Mission, Atm, G,
    mps_to_kt, kg_to_lb, W_to_hp, ft2_to_m2, kt_to_mps, isa, ftpm_to_mps,
)
from aero import (
    compute_power_curve, find_char_speeds, available_power,
    PowerCurve, CharSpeeds,
)
from simulation import MissionResult


# =================================================================
# DataFrames
# =================================================================

def _results_df(res: MissionResult) -> pd.DataFrame:
    rows = []
    for p in res.phases:
        v_str = (
            f"{mps_to_kt(p.V_mps):.1f} ({p.V_label})"
            if p.V_mps is not None and p.V_mps > 0
            else "--"
        )
        rows.append({
            "Fase": p.name,
            "t [min]": p.duration_s / 60.0,
            "d [NM]": p.distance_m / 1852.0,
            "mi [lb]": kg_to_lb(p.mass_i_kg),
            "mf [lb]": kg_to_lb(p.mass_f_kg),
            "Comb [lb]": p.fuel_lb,
            "V [kt]": v_str,
            "Pi [kW]": p.power.induced / 1e3,
            "Pp [kW]": p.power.profile / 1e3,
            "Pf [kW]": p.power.parasitic / 1e3,
            "Pmisc [kW]": p.power.misc / 1e3,
            "Pclimb [kW]": p.power.climb / 1e3,
            "Ptot [kW]": p.power.total / 1e3,
        })
    return pd.DataFrame(rows)


def _speeds_df(res: MissionResult) -> pd.DataFrame:
    rows = []
    for p in res.phases:
        if p.chars is None:
            continue
        c = p.chars
        rows.append({
            "Fase": p.name,
            "Vy [kt]": f"{mps_to_kt(c.Vy):.1f}",
            "VAM [kt]": f"{mps_to_kt(c.VAM):.1f}",
            "VDM [kt]": f"{mps_to_kt(c.VDM):.1f}",
            "VMR [kt]": f"{mps_to_kt(c.VAM):.1f}",
            "VH [kt]": f"{mps_to_kt(c.VH):.1f}",
            "Rz_max [ft/min]": f"{c.Vz_max:.0f}",
        })
    return pd.DataFrame(rows)


# =================================================================
# Console
# =================================================================

def print_summary(res: MissionResult, fuel_cap_lb: float, PU_max_W: float):
    sep = "=" * 95
    print(f"\n{sep}\n  {res.case_name}\n{sep}")

    df = _results_df(res)
    print("\n--- Consumo por Fase ---")
    print(df.to_string(index=False, float_format=lambda x: f"{x:.2f}"))

    df_s = _speeds_df(res)
    if not df_s.empty:
        print("\n--- Velocidades Caracteristicas (fases 2-5) ---")
        print("  (VMR = VAM -- menor razao de descida em autorrotacao)")
        print(df_s.to_string(index=False))

    print("\n--- Resumo de Combustivel ---")
    print(f"  Massa inicial:     {kg_to_lb(res.mass_i):>10.1f} lb")
    print(f"  Massa final:       {kg_to_lb(res.mass_f):>10.1f} lb")
    print(f"  Comb. total:       {res.fuel_lb:>10.1f} lb")
    print(f"  Capacidade tanque: {fuel_cap_lb:>10.1f} lb")
    margin = fuel_cap_lb - res.fuel_lb
    status = "OK" if margin >= 0 else "*** INSUFICIENTE ***"
    print(f"  Margem:            {margin:>10.1f} lb  ({status})")

    print("\n--- Verificacao de Potencia ---")
    for p in res.phases:
        pu_w = p.PU_at_V if p.PU_at_V else PU_max_W
        pn_w = p.power.total
        ok = "OK" if pn_w <= pu_w else "*** EXCEDE PU ***"
        print(
            f"  {p.name:<30s}  PN={pn_w/1e3:>7.1f} kW  "
            f"PU={pu_w/1e3:>7.1f} kW  [{ok}]"
        )
    print()


# =================================================================
# CSV
# =================================================================

def save_csv(res: MissionResult, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    name = res.case_name.replace(" ", "_").replace("/", "-")
    _results_df(res).to_csv(
        os.path.join(out_dir, f"{name}_consumo.csv"),
        index=False, float_format="%.3f",
    )
    df_s = _speeds_df(res)
    if not df_s.empty:
        df_s.to_csv(
            os.path.join(out_dir, f"{name}_velocidades.csv"), index=False,
        )


# =================================================================
# Plot do balanco de potencia (PN x V0)
# =================================================================

def plot_power_balance(
    ac: Aircraft,
    par: Mission,
    mass_kg: float,
    Zp_ft: float,
    Vz: float,
    case_name: str,
    phase_label: str,
    out_dir: str,
    headwind_mps: float = 0.0,
    neglect_Vz_on_inflow: bool = False,
    show_level_comparison: bool = False,
):
    """
    Balanco de potencia PN(V0) para qualquer fase.

    Se show_level_comparison=True (util na subida), sobrepoe a curva
    de voo nivelado (Vz=0) para evidenciar o efeito da correcao de
    lambda_c sobre a potencia induzida.
    """
    atm = isa(Zp_ft, par.dT)
    f_m2 = ft2_to_m2(par.f_ft2)

    curve = compute_power_curve(
        mass_kg, atm.rho, ac.rotor, par.Cd0, par.ki, f_m2, par.eta_m,
        Vz=Vz, neglect_Vz_on_inflow=neglect_Vz_on_inflow,
    )
    pu_fn = available_power(ac.P_max, atm)
    chars = find_char_speeds(
        curve, ac.P_max, mass_kg * G, par.eta_m, pu_fn, headwind_mps,
    )

    V_kt = curve.V0 / 0.514444
    fig, ax = plt.subplots(figsize=(12, 7))

    ax.plot(V_kt, curve.induced / 1e3, "--", alpha=0.7, label="Induzida")
    ax.plot(V_kt, curve.profile / 1e3, "--", alpha=0.7, label="Perfil")
    ax.plot(V_kt, curve.parasitic / 1e3, "--", alpha=0.7, label="Parasita")
    ax.plot(V_kt, curve.misc / 1e3, "--", alpha=0.7, label="Miscelanea")
    if np.any(curve.climb != 0):
        label_clb = "Subida" if Vz > 0 else "Descida"
        ax.plot(V_kt, curve.climb / 1e3, "--", color="brown", alpha=0.7,
                label=label_clb)
    ax.plot(V_kt, curve.total / 1e3, "k-", lw=2.5, label="PN total")

    if show_level_comparison and Vz != 0:
        curve_lv = compute_power_curve(
            mass_kg, atm.rho, ac.rotor, par.Cd0, par.ki, f_m2, par.eta_m,
        )
        ax.plot(V_kt, curve_lv.total / 1e3, "k--", lw=1.5, alpha=0.4,
                label="PN nivelado (sem correcao lc)")

    PU_arr = np.array([pu_fn(v) for v in curve.V0]) / 1e3
    ax.plot(V_kt, PU_arr, "r:", lw=2, label="PU (com Ram)")

    markers = [
        (chars.Vy, chars.PN_Vy, "D", "orange", "Vy"),
        (chars.VAM, chars.PN_VAM, "o", "green", "VAM"),
        (chars.VDM, chars.PN_VDM, "s", "blue", "VDM"),
    ]
    for spd, pn, mk, color, lbl in markers:
        v_kt = mps_to_kt(spd)
        ax.scatter(
            [v_kt], [pn / 1e3], s=120, marker=mk, color=color,
            edgecolors="k", zorder=5, label=f"{lbl} = {v_kt:.1f} kt",
        )

    if not np.isnan(chars.VH):
        vh_kt = mps_to_kt(chars.VH)
        idx = int(np.argmin(np.abs(V_kt - vh_kt)))
        ax.scatter(
            [vh_kt], [curve.total[idx] / 1e3], s=120, marker="^",
            color="purple", edgecolors="k", zorder=5,
            label=f"VH = {vh_kt:.1f} kt",
        )

    hw_kt = mps_to_kt(headwind_mps) if headwind_mps > 0 else 0.0
    vdm_kt = mps_to_kt(chars.VDM)
    origin_kt = hw_kt
    if vdm_kt > origin_kt:
        slope = (chars.PN_VDM / 1e3) / (vdm_kt - origin_kt)
        V_line = np.linspace(origin_kt, max(V_kt), 50)
        lbl_tang = "Tangente (VDM)"
        if hw_kt > 0:
            lbl_tang += f" - origem Vw={hw_kt:.0f} kt"
            ax.scatter(
                [origin_kt], [0], s=60, marker="x", color="red", zorder=5,
            )
        ax.plot(V_line, slope * (V_line - origin_kt), "b:", alpha=0.4,
                label=lbl_tang)

    vz_info = ""
    if Vz > 0:
        vz_info = f", Vz=+{Vz/0.3048*60:.0f} ft/min"
    elif Vz < 0:
        vz_info = f", Vz={Vz/0.3048*60:.0f} ft/min"

    ax.set_xlabel("V0 [kt]", fontsize=12)
    ax.set_ylabel("Potencia [kW]", fontsize=12)
    ax.set_title(
        f"{case_name} -- {phase_label}\n"
        f"Balanco de potencia -- Zp={Zp_ft:.0f} ft, ISA+{par.dT:.0f}, "
        f"m={mass_kg:.0f} kg{vz_info}",
        fontsize=13,
    )
    ax.legend(loc="upper center", ncol=3, fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max(V_kt))
    pu_max = max(PU_arr)
    ymax = max(pu_max * 1.15, max(curve.total / 1e3) * 1.1)
    ymin = min(0, np.min(curve.climb / 1e3) * 1.2) if np.any(curve.climb < 0) else 0
    ax.set_ylim(ymin, ymax)
    plt.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    safe = case_name.replace(" ", "_").replace("/", "-")
    ph = phase_label.replace(" ", "_")
    path = os.path.join(out_dir, f"{safe}_{ph}_balanco.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot salvo: {path}")


# =================================================================
# Polar de velocidades (Vz x V0)
# =================================================================

def plot_velocity_polar(
    ac: Aircraft,
    par: Mission,
    mass_kg: float,
    Zp_ft: float,
    case_name: str,
    phase_label: str,
    out_dir: str,
    Vz_for_curve: float = 0.0,
    neglect_Vz_on_inflow: bool = False,
):
    """
    Polar de velocidades Vz(V0) = (PU - PN_aero) * eta_m / W.

    Mostra VvM (= Vy, max Rz), Vvm (min |Rz| em descida) e VH.
    Para a subida, Vz_for_curve > 0 inclui a correcao de lambda_c
    no calculo da potencia induzida.
    """
    atm = isa(Zp_ft, par.dT)
    f_m2 = ft2_to_m2(par.f_ft2)

    curve = compute_power_curve(
        mass_kg, atm.rho, ac.rotor, par.Cd0, par.ki, f_m2, par.eta_m,
        Vz=Vz_for_curve, neglect_Vz_on_inflow=neglect_Vz_on_inflow,
    )
    pu_fn = available_power(ac.P_max, atm)
    W = mass_kg * G

    PN_aero = curve.total - curve.climb
    PU_arr = np.array([pu_fn(v) for v in curve.V0])
    Vz_arr_mps = (PU_arr - PN_aero) * par.eta_m / W
    Vz_arr_ftpm = Vz_arr_mps / 0.3048 * 60.0

    V_kt = curve.V0 / 0.514444

    i_vmax = int(np.argmax(Vz_arr_ftpm))
    Vy_kt = V_kt[i_vmax]
    Vz_max_ftpm = Vz_arr_ftpm[i_vmax]

    neg_mask = Vz_arr_ftpm < 0
    if neg_mask.any():
        neg_idx = np.where(neg_mask)[0]
        i_vmin_abs = neg_idx[int(np.argmax(Vz_arr_ftpm[neg_mask]))]
        Vvm_kt = V_kt[i_vmin_abs]
        Vz_min_ftpm = Vz_arr_ftpm[i_vmin_abs]
    else:
        Vvm_kt = None
        Vz_min_ftpm = None

    zero_cross = np.where(np.diff(np.sign(Vz_arr_ftpm)))[0]
    VH_kt = V_kt[zero_cross[-1]] if len(zero_cross) > 0 else None

    fig, ax = plt.subplots(figsize=(12, 6))

    ax.plot(V_kt, Vz_arr_ftpm, "k-", lw=2, label="Vz(V0)")
    ax.axhline(0, color="gray", lw=0.8, ls="--")

    ax.scatter([Vy_kt], [Vz_max_ftpm], s=120, marker="D", color="orange",
               edgecolors="k", zorder=5,
               label=f"VvM (Vy) = {Vy_kt:.1f} kt, Rz_max = {Vz_max_ftpm:.0f} ft/min")

    if Vvm_kt is not None:
        ax.scatter([Vvm_kt], [Vz_min_ftpm], s=120, marker="v", color="green",
                   edgecolors="k", zorder=5,
                   label=f"Vvm = {Vvm_kt:.1f} kt, Rz_min = {Vz_min_ftpm:.0f} ft/min")

    if VH_kt is not None:
        ax.scatter([VH_kt], [0], s=100, marker="^", color="purple",
                   edgecolors="k", zorder=5,
                   label=f"VH = {VH_kt:.1f} kt")

    vz_info = ""
    if Vz_for_curve > 0:
        vz_info = f" (com correcao lc, Vz=+{Vz_for_curve/0.3048*60:.0f} ft/min)"
    elif Vz_for_curve < 0:
        vz_info = f" (Vz={Vz_for_curve/0.3048*60:.0f} ft/min)"

    ax.set_xlabel("V0 [kt]", fontsize=12)
    ax.set_ylabel("Vz [ft/min]", fontsize=12)
    ax.set_title(
        f"{case_name} -- {phase_label}\n"
        f"Polar de velocidades -- Zp={Zp_ft:.0f} ft, ISA+{par.dT:.0f}, "
        f"m={mass_kg:.0f} kg{vz_info}",
        fontsize=13,
    )
    ax.legend(loc="lower left", fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max(V_kt))
    plt.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    safe = case_name.replace(" ", "_").replace("/", "-")
    ph = phase_label.replace(" ", "_")
    path = os.path.join(out_dir, f"{safe}_{ph}_polar.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot salvo: {path}")


# =================================================================
# Gera todos os plots para fases 2-5
# =================================================================

def plot_all_phases(ac: Aircraft, par: Mission, res: MissionResult,
                    out_dir: str):
    """Gera balanco de potencia + polar de velocidades para fases 2 a 5."""

    Vz_climb = ftpm_to_mps(par.climb_rate_ftpm)
    Vz_desc = -ftpm_to_mps(par.desc_rate_ftpm)
    Zp_climb_avg = 0.5 * (par.Zp_to_ft + par.Zp_cr_ft)
    Zp_desc_avg = Zp_climb_avg
    hw_mps = kt_to_mps(par.headwind_kt) if par.headwind_kt else 0.0

    phase_defs = [
        {
            "idx": 1, "label": "Fase 2 - Subida",
            "Zp": Zp_climb_avg, "Vz": Vz_climb,
            "neglect": False, "hw": 0.0, "show_lv": True,
        },
        {
            "idx": 2, "label": "Fase 3 - Cruzeiro VDM",
            "Zp": par.Zp_cr_ft, "Vz": 0.0,
            "neglect": False, "hw": hw_mps, "show_lv": False,
        },
        {
            "idx": 3, "label": "Fase 4 - Reserva VAM",
            "Zp": par.Zp_cr_ft, "Vz": 0.0,
            "neglect": False, "hw": 0.0, "show_lv": False,
        },
        {
            "idx": 4, "label": "Fase 5 - Descida",
            "Zp": Zp_desc_avg, "Vz": Vz_desc,
            "neglect": True, "hw": 0.0, "show_lv": False,
        },
    ]

    for pd_ in phase_defs:
        ph = res.phases[pd_["idx"]]
        m_avg = 0.5 * (ph.mass_i_kg + ph.mass_f_kg)

        plot_power_balance(
            ac, par, m_avg, pd_["Zp"], pd_["Vz"],
            res.case_name, pd_["label"], out_dir,
            headwind_mps=pd_["hw"],
            neglect_Vz_on_inflow=pd_["neglect"],
            show_level_comparison=pd_["show_lv"],
        )

        Vz_polar = pd_["Vz"] if pd_["Vz"] > 0 else 0.0
        plot_velocity_polar(
            ac, par, m_avg, pd_["Zp"],
            res.case_name, pd_["label"], out_dir,
            Vz_for_curve=Vz_polar,
        )
