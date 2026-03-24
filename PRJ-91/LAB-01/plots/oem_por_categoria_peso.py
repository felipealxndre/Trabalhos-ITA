"""Gráficos de pizza individuais: OEM por categoria de peso (turbina)."""

import os
import re
import unicodedata
from .styles import pizza

CAT_ORDER = ['Leve', 'Médio', 'Médio-pesado', 'Pesado']


def _top_oem(sub, n=4):
    oem = sub.groupby('Fabricante Atual')['Unidades'].sum().sort_values(ascending=False)
    top = oem.head(n)
    outros = oem.iloc[n:].sum()
    if outros > 0:
        top['Outros'] = outros
    return dict(top)


def _slug(name):
    s = unicodedata.normalize('NFKD', name).encode('ascii', 'ignore').decode()
    return re.sub(r'[^a-z0-9]+', '_', s.lower()).strip('_')


def plot(civil, results_dir):
    turbina = civil[civil['Tipo de Motor'] == 'Turbina'].copy()

    out = os.path.join(results_dir, 'categoria_peso')
    os.makedirs(out, exist_ok=True)
    for cat in CAT_ORDER:
        sub = turbina[turbina['Categoria de Peso'] == cat]
        data = _top_oem(sub)
        pizza(data, os.path.join(out, f'oem_peso_{_slug(cat)}.png'),
              callout=['Outros'])
