"""Gráfico de pizza: distribuição da frota civil brasileira por OEM."""

import os
from .styles import pizza


def plot(civil, results_dir):
    oem = civil.groupby('Fabricante Atual')['Unidades'].sum().sort_values(ascending=False)

    top = oem.head(6)
    outros = oem.iloc[6:].sum()
    if outros > 0:
        top['Outros'] = outros

    out = os.path.join(results_dir, 'distribuicao_frota')
    os.makedirs(out, exist_ok=True)
    pizza(dict(top), os.path.join(out, 'participacao_oem.png'))
