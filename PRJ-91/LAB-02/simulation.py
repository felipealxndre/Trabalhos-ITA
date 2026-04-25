"""
Fases da missao e simulacao completa do AH-1S.

Sequencia de 6 fases:
  1) Pairado IGE (5 min, h=6 ft, SL)
  2) Subida em Vy a Rz prescrita ate FL050
  3) Cruzeiro nivelado em VDM (FL050, distancia prescrita)
  4) Reserva em VAM (30 min, FL050)
  5) Descida em Vy a Rz prescrita ate SL
  6) Pairado IGE (5 min, h=6 ft, SL)
"""

from dataclasses import dataclass
from typing import List, Optional

from config import (
    Aircraft, Mission, Atm, G,
    isa, ft_to_m, ft2_to_m2, ftpm_to_mps, kt_to_mps,
    mps_to_kt, NM_to_m, sfc_to_SI, kg_to_lb,
)
from aero import (
    Power, PowerCurve, CharSpeeds,
    power_required, ground_effect, available_power,
    compute_power_curve, find_char_speeds,
)


@dataclass
class PhaseResult:
    name: str
    duration_s: float
    distance_m: float
    mass_i_kg: float
    mass_f_kg: float
    fuel_kg: float
    power: Power
    V_mps: Optional[float] = None
    V_label: Optional[str] = None
    chars: Optional[CharSpeeds] = None
    atm: Optional[Atm] = None
    PU_at_V: Optional[float] = None
    notes: str = ""

    @property
    def fuel_lb(self) -> float:
        return kg_to_lb(self.fuel_kg)


def _fuel(P_W: float, t_s: float, sfc_SI: float) -> float:
    """Combustivel [kg] = SFC [kg/(W*s)] * P [W] * t [s]."""
    return sfc_SI * P_W * t_s


def _avg_weight_iter(m0, t_s, calc_pwr, sfc_SI, n=3):
    """Itera peso medio para fases longas (refinamento do consumo)."""
    fuel = _fuel(calc_pwr(m0).total, t_s, sfc_SI)
    for _ in range(n):
        m_avg = m0 - 0.5 * fuel
        pwr = calc_pwr(m_avg)
        fuel = _fuel(pwr.total, t_s, sfc_SI)
    return fuel, pwr, m_avg


def _curve_and_chars(ac, par, atm, mass, headwind_mps=0.0, Vz=0.0,
                     neglect_Vz_on_inflow=False):
    """Helper: gera curva PN e extrai velocidades com efeito Ram."""
    f_m2 = ft2_to_m2(par.f_ft2)
    curve = compute_power_curve(
        mass, atm.rho, ac.rotor, par.Cd0, par.ki, f_m2, par.eta_m,
        Vz=Vz, neglect_Vz_on_inflow=neglect_Vz_on_inflow,
    )
    pu_fn = available_power(ac.P_max, atm)
    chars = find_char_speeds(
        curve, ac.P_max, mass * G, par.eta_m, pu_fn, headwind_mps
    )
    return curve, chars, f_m2, pu_fn


# Fase 1 / 6: Pairado IGE

def _phase_hover_ige(mass, dur_min, h_ft, Zp_ft, ac, par, name):
    atm = isa(Zp_ft, par.dT)
    t = dur_min * 60.0
    sfc_SI = sfc_to_SI(par.sfc)
    f_m2 = ft2_to_m2(par.f_ft2)
    ge = ground_effect(ft_to_m(h_ft) / (2.0 * ac.rotor.radius))

    pwr = power_required(
        mass, 0.0, 0.0, atm.rho, ac.rotor,
        par.Cd0, par.ki, f_m2, par.eta_m, ge,
    )
    fuel = _fuel(pwr.total, t, sfc_SI)
    pu0 = ac.P_max

    return PhaseResult(
        name, t, 0.0, mass, mass - fuel, fuel, pwr,
        0.0, "hover", None, atm, pu0,
        f"IGE h={h_ft:.0f} ft, z/D={ft_to_m(h_ft)/(2*ac.rotor.radius):.3f}",
    )


# Fase 2: Subida em Vy

def _phase_climb(mass, Zp_i, Zp_f, ac, par):
    Zp_avg = 0.5 * (Zp_i + Zp_f)
    atm = isa(Zp_avg, par.dT)
    Vz = ftpm_to_mps(par.climb_rate_ftpm)
    t = (Zp_f - Zp_i) / par.climb_rate_ftpm * 60.0
    sfc_SI = sfc_to_SI(par.sfc)

    _, chars, f_m2, pu_fn = _curve_and_chars(ac, par, atm, mass, Vz=Vz)
    Vy = chars.Vy

    calc = lambda m: power_required(
        m, Vy, Vz, atm.rho, ac.rotor,
        par.Cd0, par.ki, f_m2, par.eta_m,
    )
    fuel, pwr, _ = _avg_weight_iter(mass, t, calc, sfc_SI)

    return PhaseResult(
        "2 - Subida (Vy)", t, Vy * t, mass, mass - fuel, fuel, pwr,
        Vy, "Vy", chars, atm, pu_fn(Vy),
        f"Vz={par.climb_rate_ftpm:.0f} ft/min, Zp_avg={Zp_avg:.0f} ft",
    )


# Fase 3: Cruzeiro em VDM

def _phase_cruise(mass, ac, par):
    atm = isa(par.Zp_cr_ft, par.dT)
    sfc_SI = sfc_to_SI(par.sfc)
    hw_mps = kt_to_mps(par.headwind_kt)

    _, chars, f_m2, pu_fn = _curve_and_chars(ac, par, atm, mass, hw_mps)
    VDM = chars.VDM
    Vg = VDM - hw_mps
    if Vg <= 0:
        raise ValueError(
            f"Vento ({par.headwind_kt} kt) >= VDM ({mps_to_kt(VDM):.1f} kt)"
        )

    dist = NM_to_m(par.distance_NM)
    t = dist / Vg

    calc = lambda m: power_required(
        m, VDM, 0.0, atm.rho, ac.rotor,
        par.Cd0, par.ki, f_m2, par.eta_m,
    )
    fuel, pwr, m_avg = _avg_weight_iter(mass, t, calc, sfc_SI)

    _, chars_avg, _, _ = _curve_and_chars(ac, par, atm, m_avg, hw_mps)

    return PhaseResult(
        "3 - Cruzeiro (VDM)", t, dist, mass, mass - fuel, fuel, pwr,
        VDM, "VDM", chars_avg, atm, pu_fn(VDM),
        f"{par.distance_NM:.0f} NM, vento={par.headwind_kt:.0f} kt, "
        f"Vg={mps_to_kt(Vg):.1f} kt",
    )


# Fase 4: Reserva em VAM

def _phase_reserve(mass, ac, par):
    atm = isa(par.Zp_cr_ft, par.dT)
    sfc_SI = sfc_to_SI(par.sfc)
    t = par.reserve_min * 60.0

    _, chars, f_m2, pu_fn = _curve_and_chars(ac, par, atm, mass)
    VAM = chars.VAM

    calc = lambda m: power_required(
        m, VAM, 0.0, atm.rho, ac.rotor,
        par.Cd0, par.ki, f_m2, par.eta_m,
    )
    fuel, pwr, _ = _avg_weight_iter(mass, t, calc, sfc_SI)

    return PhaseResult(
        "4 - Reserva (VAM)", t, VAM * t, mass, mass - fuel, fuel, pwr,
        VAM, "VAM", chars, atm, pu_fn(VAM),
        f"{par.reserve_min:.0f} min",
    )


# Fase 5: Descida em Vy
# Pi calculada como em voo nivelado (lc=0 no Glauert).
# lc*CT mantido para contabilizar ganho de energia potencial.

def _phase_descent(mass, Zp_i, Zp_f, ac, par):
    Zp_avg = 0.5 * (Zp_i + Zp_f)
    atm = isa(Zp_avg, par.dT)
    Vz = -ftpm_to_mps(par.desc_rate_ftpm)
    t = (Zp_i - Zp_f) / par.desc_rate_ftpm * 60.0
    sfc_SI = sfc_to_SI(par.sfc)

    _, chars, f_m2, pu_fn = _curve_and_chars(ac, par, atm, mass)
    Vy = chars.Vy

    calc = lambda m: power_required(
        m, Vy, Vz, atm.rho, ac.rotor,
        par.Cd0, par.ki, f_m2, par.eta_m,
        neglect_Vz_on_inflow=True,
    )
    fuel, pwr, _ = _avg_weight_iter(mass, t, calc, sfc_SI)

    return PhaseResult(
        "5 - Descida (Vy)", t, Vy * t, mass, mass - fuel, fuel, pwr,
        Vy, "Vy", chars, atm, pu_fn(Vy),
        f"Vz=-{par.desc_rate_ftpm:.0f} ft/min, Pi com lc=0 (lab)",
    )


@dataclass
class MissionResult:
    case_name: str
    phases: List[PhaseResult]
    mass_i: float
    mass_f: float

    @property
    def fuel_kg(self) -> float:
        return sum(p.fuel_kg for p in self.phases)

    @property
    def fuel_lb(self) -> float:
        return kg_to_lb(self.fuel_kg)


def simulate(ac: Aircraft, par: Mission) -> MissionResult:
    """Roda as 6 fases da missao e retorna resultado consolidado."""
    m = ac.max_takeoff
    m0 = m
    phases: List[PhaseResult] = []

    p = _phase_hover_ige(m, par.hover_min, par.hover_h_ft, par.Zp_to_ft,
                         ac, par, "1 - Pairado IGE (inicio)")
    phases.append(p); m = p.mass_f_kg

    p = _phase_climb(m, par.Zp_to_ft, par.Zp_cr_ft, ac, par)
    phases.append(p); m = p.mass_f_kg

    p = _phase_cruise(m, ac, par)
    phases.append(p); m = p.mass_f_kg

    p = _phase_reserve(m, ac, par)
    phases.append(p); m = p.mass_f_kg

    p = _phase_descent(m, par.Zp_cr_ft, par.Zp_to_ft, ac, par)
    phases.append(p); m = p.mass_f_kg

    p = _phase_hover_ige(m, par.hover_min, par.hover_h_ft, par.Zp_to_ft,
                         ac, par, "6 - Pairado IGE (final)")
    phases.append(p); m = p.mass_f_kg

    return MissionResult(par.name, phases, m0, m)
