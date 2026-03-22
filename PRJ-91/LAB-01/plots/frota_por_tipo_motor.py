"""Gráfico de barras: distribuição da frota por categoria de motor."""

import os
from .styles import barras

CAT_ORDER = ['Pistão', 'Monomotor\nTurbina', 'Bimotor\nTurbina']


def _classificar_motor(row):
    if row['Tipo de Motor'] == 'Pistão':
        return 'Pistão'
    elif row['Nº Motores'] == 1:
        return 'Monomotor\nTurbina'
    else:
        return 'Bimotor\nTurbina'


def plot(civil, results_dir):
    civil = civil.copy()
    civil['Motor Cat'] = civil.apply(_classificar_motor, axis=1)
    totals = civil.groupby('Motor Cat')['Unidades'].sum().reindex(CAT_ORDER)

    out = os.path.join(results_dir, 'categoria_motor')
    os.makedirs(out, exist_ok=True)
    barras(
        dict(totals),
        'Categoria', 'Nº de Aeronaves',
        os.path.join(out, 'frota_por_tipo_motor.png'),
    )
