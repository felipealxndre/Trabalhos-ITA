from .participacao_oem import plot as plot_participacao_oem
from .frota_por_tipo_motor import plot as plot_frota_por_tipo_motor
from .frota_por_categoria_peso import plot as plot_frota_por_categoria_peso
from .oem_por_categoria_peso import plot as plot_oem_por_categoria_peso
from .oem_por_tipo_motor import plot as plot_oem_por_tipo_motor

PLOTS = [
    ('participacao_oem', plot_participacao_oem),
    ('frota_por_tipo_motor', plot_frota_por_tipo_motor),
    ('frota_por_categoria_peso', plot_frota_por_categoria_peso),
    ('oem_por_categoria_peso', plot_oem_por_categoria_peso),
    ('oem_por_tipo_motor', plot_oem_por_tipo_motor),
]
