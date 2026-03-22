"""Ponto de entrada: python -m plots"""

import os
import matplotlib
matplotlib.use('Agg')
import pandas as pd

from . import PLOTS


results_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')


def main():
    os.makedirs(results_dir, exist_ok=True)

    link = 'https://docs.google.com/spreadsheets/d/e/2PACX-1vR4a2Stgdml-pOp_Bu9sjHo8ptlAghpXB-kO386j-A9_MXgTGcAVnDBJ5m4B-sZRJWoHFGsmI2h2cI2/pub?output=xlsx'
    df = pd.read_excel(link, sheet_name='Tabela Final')

    print(df.head(10))
    civil = df[df['Unidades'] > 0].copy()

    for name, plot_fn in PLOTS:
        plot_fn(civil, results_dir)

if __name__ == '__main__':
    main()
