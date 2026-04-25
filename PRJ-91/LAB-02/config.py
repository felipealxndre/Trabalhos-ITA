"""
Dados do AH-1S Cobra, parametros de missao, atmosfera ISA e conversoes.
Fonte: Prouty -- Helicopter Performance, Stability and Control.
"""

import math
from dataclasses import dataclass, replace
from typing import List

G = 9.80665  # m/s^2

# Conversoes de unidades

def ft_to_m(x):     return x * 0.3048
def ft2_to_m2(x):   return x * 0.09290304
def kt_to_mps(x):   return x * 0.514444
def mps_to_kt(x):   return x / 0.514444
def ftps_to_mps(x): return x * 0.3048
def ftpm_to_mps(x): return x * 0.3048 / 60.0
def lb_to_kg(x):    return x * 0.45359237
def kg_to_lb(x):    return x / 0.45359237
def hp_to_W(x):     return x * 745.7
def W_to_hp(x):     return x / 745.7
def NM_to_m(x):     return x * 1852.0

def sfc_to_SI(sfc_lb_per_hphr):
    """lb/(hp*hr) -> kg/(W*s)."""
    return sfc_lb_per_hphr * 0.45359237 / (745.7 * 3600.0)

# Atmosfera ISA

_T0 = 288.15       # K   (temperatura padrao ao nivel do mar)
_P0 = 101325.0     # Pa  (pressao padrao)
_LAPSE = 0.0065    # K/m (gradiente termico troposferico)
_R = 287.05287     # J/(kg*K)


@dataclass(frozen=True)
class Atm:
    """Estado atmosferico."""
    T: float    # K
    P: float    # Pa
    rho: float  # kg/m^3


def isa(Zp_ft: float, dT: float = 0.0) -> Atm:
    """
    ISA com offset de temperatura.
    Pressao usa T padrao; densidade usa T efetiva (convencao Zp + ISA+dT).
    """
    h = Zp_ft * 0.3048
    T_std = _T0 - _LAPSE * h
    T_eff = T_std + dT
    P = _P0 * (T_std / _T0) ** (G / (_R * _LAPSE))
    return Atm(T_eff, P, P / (_R * T_eff))

# Rotor e aeronave

@dataclass(frozen=True)
class Rotor:
    radius: float       # m
    solidity: float     # sigma
    tip_speed: float    # Omega*R [m/s]

    @property
    def area(self) -> float:
        return math.pi * self.radius ** 2


@dataclass(frozen=True)
class Aircraft:
    max_takeoff: float      # kg
    fuel_capacity: float    # kg
    P_max: float            # W  (potencia max. decolagem)
    P_cont: float           # W  (potencia max. continua)
    rotor: Rotor


AH1S = Aircraft(
    max_takeoff=lb_to_kg(10000),
    fuel_capacity=lb_to_kg(1684),
    P_max=hp_to_W(1800),
    P_cont=hp_to_W(1290),
    rotor=Rotor(
        radius=ft_to_m(22.0),
        solidity=0.065,
        tip_speed=ftps_to_mps(746.0),
    ),
)

# Parametros de missao

@dataclass(frozen=True)
class Mission:
    name: str
    distance_NM: float
    climb_rate_ftpm: float
    headwind_kt: float
    Cd0: float = 0.012
    ki: float = 1.15
    eta_m: float = 0.85
    f_ft2: float = 19.45            # Area placa plana equivalente [ft^2]
    sfc: float = 0.458              # lb/(hp*hr)
    dT: float = 20.0                # ISA+20
    Zp_to_ft: float = 0.0          # altitude partida/chegada
    Zp_cr_ft: float = 5000.0       # altitude cruzeiro
    hover_h_ft: float = 6.0        # altura pairado IGE
    hover_min: float = 5.0
    reserve_min: float = 30.0
    desc_rate_ftpm: float = 1000.0


CASE1 = Mission("Caso 1 - Baseline", 400, 1000, 0)
CASE2 = replace(CASE1, name="Caso 2 - Vento 15 kt", headwind_kt=15)
CASE3 = replace(CASE1, name="Caso 3 - Subida 2000 ft/min", climb_rate_ftpm=2000)
CASE4 = replace(CASE1, name="Caso 4 - Distancia 440 NM", distance_NM=440)
ALL_CASES: List[Mission] = [CASE1, CASE2, CASE3, CASE4]
