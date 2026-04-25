"""
PRJ-91 / Laboratorio no. 3 -- Requisitos de Missao do AH-1S Cobra.

Roda os 4 cenarios, gera tabelas de consumo e plots da curva PN x V0.

Uso:
    python main.py
"""

import os

from config import AH1S, ALL_CASES, kg_to_lb
from simulation import simulate
from report import print_summary, save_csv, plot_all_phases

RESULTS_DIR = os.path.join(os.path.dirname(__file__), "output", "results")
PLOTS_DIR = os.path.join(os.path.dirname(__file__), "output", "plots")


def main():
    fuel_cap_lb = kg_to_lb(AH1S.fuel_capacity)

    print("\n" + "=" * 95)
    print("  PRJ-91 / Laboratorio no. 3 -- AH-1S Cobra")
    print("=" * 95)
    print(f"  Massa max. decolagem:   {kg_to_lb(AH1S.max_takeoff):.0f} lb")
    print(f"  Capacidade combustivel: {fuel_cap_lb:.0f} lb")
    print(f"  Potencia max. decol.:   {AH1S.P_max/745.7:.0f} hp")
    print(f"  Potencia max. contin.:  {AH1S.P_cont/745.7:.0f} hp")

    for case in ALL_CASES:
        res = simulate(AH1S, case)

        print_summary(res, fuel_cap_lb, AH1S.P_max)
        save_csv(res, RESULTS_DIR)
        plot_all_phases(AH1S, case, res, PLOTS_DIR)

    print("\n" + "=" * 95)
    print("  Simulacao concluida.")
    print(f"  CSVs em:  {RESULTS_DIR}")
    print(f"  Plots em: {PLOTS_DIR}")
    print("=" * 95 + "\n")


if __name__ == "__main__":
    main()
