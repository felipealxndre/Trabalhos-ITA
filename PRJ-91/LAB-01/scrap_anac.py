import re
import pandas as pd
import requests
from bs4 import BeautifulSoup
from difflib import get_close_matches

df = pd.read_csv('https://docs.google.com/spreadsheets/d/1Bo60Tu0unJoOLFab4Mrt0wk93rJD0KBDaH0HuiL7Ds4/export?format=csv')

BASE_URL = "https://aeronaves.anac.gov.br/aeronaves"

resp = requests.get(f"{BASE_URL}/get_dados_dominio_filtros.asp?tipo=tipo_icao")
icao_anac = {item["id"]: item["nome"] for item in resp.json()}
icao_ids = list(icao_anac.keys())

def encontrar_icao(type_designator: str) -> str | None:
    td = type_designator.strip().upper()
    if td in icao_anac:
        return td
    matches = get_close_matches(td, icao_ids, n=1, cutoff=0.6)
    return matches[0] if matches else None

def get_num_aeronaves(icao_code: str) -> int:
    url = (
        f"{BASE_URL}/consulta-rab-1-iframe-resposta-4.asp"
        f"?textMarca=&selectHabilitacao=&selectIcao={icao_code}"
        f"&selectModelo=&selectFabricante=&textNumeroSerie="
    )
    resp = requests.get(url)
    resp.encoding = "iso-8859-1"
    soup = BeautifulSoup(resp.text, "html.parser")

    for tag in soup.find_all(string=True):
        if "Total de registros encontrados" in tag:
            match = re.search(r'(\d+)', tag)
            if match:
                return int(match.group(1))
    return 0

# Para cada Type Designator único, acha o ICAO correspondente e busca o total
icao_map = {} 
total_map = {}

for td in df["Type Designator"].dropna().unique():
    icao = encontrar_icao(td)
    icao_map[td] = icao
    total_map[td] = get_num_aeronaves(icao) if icao else None

df["icao_match"] = df["Type Designator"].map(icao_map)
df["num_aeronaves"] = df["Type Designator"].map(total_map)

print(df[["Type Designator", "Model", "icao_match", "num_aeronaves"]])

df.to_excel("aeronaves_anac.xlsx", index=False)