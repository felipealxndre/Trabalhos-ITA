"""Gráfico de barras: distribuição da frota a turbina por categoria de peso."""

import os
from .styles import barras

CAT_ORDER = ['Leve', 'Médio', 'Médio-pesado', 'Pesado']


def plot(civil, results_dir):
    turbina = civil[civil['Tipo de Motor'] == 'Turbina'].copy()
    totals = turbina.groupby('Categoria de Peso')['Unidades'].sum().reindex(CAT_ORDER)

    out = os.path.join(results_dir, 'categoria_peso')
    os.makedirs(out, exist_ok=True)
    barras(
        dict(totals),
        'Categoria de Peso', 'Nº de Aeronaves',
        os.path.join(out, 'frota_por_categoria_peso.png'),
    )
