"""
Aerodinamica do rotor: inflow de Glauert (iterativo), potencia necessaria,
efeito solo, curva de desempenho PN(V0), velocidades caracteristicas.
"""

import math
import numpy as np
from dataclasses import dataclass
from typing import Callable, Optional
from scipy.optimize import curve_fit

from config import Rotor, Atm, G, ft2_to_m2

# Efeito solo -- ajuste a Fig. 3 do Prouty (curva tracejada, ensaio em voo)

_GE_DATA = np.array([
    [0.10, 0.55], [0.15, 0.65], [0.20, 0.72], [0.25, 0.77],
    [0.30, 0.82], [0.40, 0.87], [0.50, 0.92], [0.60, 0.95],
    [0.75, 0.97], [1.00, 1.00], [1.20, 1.00],
])
_ge_model = lambda z, a, b: 1.0 - a * np.exp(-b * z)
_A_GE, _B_GE = curve_fit(
    _ge_model, _GE_DATA[:, 0], _GE_DATA[:, 1], p0=[0.5, 3.0]
)[0]


def ground_effect(z_over_D: float) -> float:
    """Vi_IGE / Vi_OGE -- fator multiplicativo sobre a potencia induzida."""
    if z_over_D >= 1.0:
        return 1.0
    if z_over_D <= 0.0:
        return 1.0 - _A_GE
    return float(_ge_model(z_over_D, _A_GE, _B_GE))


def solve_inflow(CT: float, mu: float, lc: float = 0.0) -> float:
    """Resolve lambda_i implicito por iteracao de ponto fixo."""
    li = math.sqrt(max(CT, 0.0) / 2.0)
    for _ in range(200):
        d = 2.0 * math.sqrt(mu ** 2 + (lc + li) ** 2)
        if d < 1e-12:
            return 0.0
        f = CT / d
        li_new = li + 0.5 * (f - li)
        if abs(li_new - li) < 1e-8:
            return li_new
        li = li_new
    raise RuntimeError(f"Inflow nao convergiu (CT={CT:.5f}, mu={mu:.4f})")


@dataclass
class Power:
    """Decomposicao de potencia [W]."""
    induced: float
    profile: float
    parasitic: float
    misc: float
    climb: float

    @property
    def shaft(self) -> float:
        """Potencia no eixo do rotor (ind + perfil + parasita)."""
        return self.induced + self.profile + self.parasitic

    @property
    def total(self) -> float:
        """Potencia total no eixo do motor."""
        return self.shaft + self.misc + self.climb


def power_required(
    mass_kg: float,
    V0: float,
    Vz: float,
    rho: float,
    rotor: Rotor,
    Cd0: float,
    ki: float,
    f_m2: float,
    eta_m: float,
    ge_factor: float = 1.0,
    neglect_Vz_on_inflow: bool = False,
) -> Power:
    """
    Potencia total no eixo do motor [W].

    Se neglect_Vz_on_inflow=True e Vz<0, a potencia induzida e calculada
    como em voo nivelado (lc=0 no Glauert), conforme instrucao do lab
    para a fase de descida. A parcela lc*CT continua usando o Vz real.
    """
    CT = mass_kg * G / (rho * rotor.area * rotor.tip_speed ** 2)
    mu = V0 / rotor.tip_speed
    lc = Vz / rotor.tip_speed

    lc_inflow = 0.0 if (neglect_Vz_on_inflow and Vz < 0) else lc
    li = solve_inflow(CT, mu, lc_inflow)

    scale = rho * rotor.area * rotor.tip_speed ** 3
    denom = 2.0 * math.sqrt(mu ** 2 + (lc_inflow + li) ** 2)

    Cp_i = ki * CT ** 2 / denom * ge_factor
    Cp_p = (rotor.solidity * Cd0 / 8.0) * (1.0 + 4.65 * mu ** 2)
    Cp_f = 0.5 * (f_m2 / rotor.area) * mu ** 3
    Cp_eR = Cp_i + Cp_p + Cp_f
    Cp_misc = (1.0 / eta_m - 1.0) * Cp_eR
    Cp_climb = lc * CT

    return Power(
        Cp_i * scale, Cp_p * scale, Cp_f * scale,
        Cp_misc * scale, Cp_climb * scale,
    )


def available_power(PU0: float, atm: Atm, k_ram: float = 0.7) -> Callable:
    """Retorna PU(V0) com efeito Ram simplificado."""
    c = k_ram * 0.5 * atm.rho / atm.P
    return lambda V0: PU0 * (1.0 + c * V0 ** 2)


@dataclass
class PowerCurve:
    V0: np.ndarray       # m/s
    total: np.ndarray    # W
    induced: np.ndarray
    profile: np.ndarray
    parasitic: np.ndarray
    misc: np.ndarray
    climb: np.ndarray    # W (positivo subida, negativo descida, zero nivelado)


def compute_power_curve(
    mass_kg: float,
    rho: float,
    rotor: Rotor,
    Cd0: float,
    ki: float,
    f_m2: float,
    eta_m: float,
    Vz: float = 0.0,
    neglect_Vz_on_inflow: bool = False,
    V_max: float = 100.0,
    n: int = 2000,
) -> PowerCurve:
    """
    Gera curva PN(V0) com 2000 pontos (~0.1 kt).

    Se Vz != 0, inclui efeito de lambda_c no inflow (subida) e parcela
    lambda_c*CT. Para descida (Vz<0), usar neglect_Vz_on_inflow=True.
    """
    V = np.linspace(0.5, V_max, n)
    tot = np.zeros(n)
    ind = np.zeros(n)
    pro = np.zeros(n)
    par = np.zeros(n)
    mis = np.zeros(n)
    clb = np.zeros(n)
    for j, v in enumerate(V):
        pw = power_required(
            mass_kg, v, Vz, rho, rotor, Cd0, ki, f_m2, eta_m,
            neglect_Vz_on_inflow=neglect_Vz_on_inflow,
        )
        tot[j] = pw.total
        ind[j] = pw.induced
        pro[j] = pw.profile
        par[j] = pw.parasitic
        mis[j] = pw.misc
        clb[j] = pw.climb
    return PowerCurve(V, tot, ind, pro, par, mis, clb)


@dataclass
class CharSpeeds:
    Vy: float
    VAM: float
    VDM: float
    VH: float
    PN_Vy: float
    PN_VAM: float
    PN_VDM: float
    Vz_max: float  # ft/min


def find_char_speeds(
    curve: PowerCurve,
    PU0_W: float,
    weight_N: float,
    eta_m: float,
    PU_func: Optional[Callable] = None,
    headwind_mps: float = 0.0,
) -> CharSpeeds:
    """
    Extrai velocidades caracteristicas da curva PN(V0).
    VDM com vento: minimiza PN/(V_TAS - V_headwind).
    """
    V, PN = curve.V0, curve.total
    PU = np.array([PU_func(v) for v in V]) if PU_func else np.full_like(V, PU0_W)

    i_am = int(np.argmin(PN))

    Vg = V - headwind_mps
    ratio = np.where(Vg > 1.0, PN / Vg, np.inf)
    i_dm = int(np.argmin(ratio))

    PN_aero = PN - curve.climb
    excess = PU - PN_aero
    i_vy = int(np.argmax(excess))
    Vz_max_mps = excess[i_vy] * eta_m / weight_N

    ok = PN <= PU
    VH = float(V[ok][-1]) if ok.any() else float("nan")

    return CharSpeeds(
        Vy=V[i_vy], VAM=V[i_am], VDM=V[i_dm], VH=VH,
        PN_Vy=PN[i_vy], PN_VAM=PN[i_am], PN_VDM=PN[i_dm],
        Vz_max=Vz_max_mps / 0.3048 * 60.0,
    )
