"""Gráficos de pizza individuais: OEM por categoria de motor."""

import os
import re
import unicodedata
from .styles import pizza

CAT_ORDER = ['Pistão', 'Monomotor Turbina', 'Bimotor Turbina']


def _classificar_motor(row):
    if row['Tipo de Motor'] == 'Pistão':
        return 'Pistão'
    elif row['Nº Motores'] == 1:
        return 'Monomotor Turbina'
    else:
        return 'Bimotor Turbina'


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
    civil = civil.copy()
    civil['Motor Cat'] = civil.apply(_classificar_motor, axis=1)

    out = os.path.join(results_dir, 'categoria_motor')
    os.makedirs(out, exist_ok=True)
    callout_map = {
        'Pistão': ['Outros', 'Enstrom'],
        'Monomotor Turbina': ['Outros'],
        'Bimotor Turbina': ['Outros'],
    }

    for cat in CAT_ORDER:
        sub = civil[civil['Motor Cat'] == cat]
        data = _top_oem(sub)
        pizza(data, os.path.join(out, f'oem_motor_{_slug(cat)}.png'),
              callout=callout_map.get(cat, []))
