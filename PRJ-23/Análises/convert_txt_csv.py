import pandas as pd
import io

def convert_forces_txt_to_csv(txt_filename, output_prefix="strip_forces"):
    
    try:
        with open(txt_filename, 'r') as f:
            file_content = f.read()
    except FileNotFoundError:
        print(f"Erro: Arquivo '{txt_filename}' não encontrado.")
        return

    # O cabeçalho da tabela de Strip Forces é o mesmo para todas as superfícies
    HEADER_LINE = "j      Yle    Chord     Area     c cl      ai      cl_norm  cl       cd       cdv    cm_c/4    cm_LE  C.P.x/c"
    
    # Colunas a serem usadas
    column_names = ['j', 'Yle', 'Chord', 'Area', 'c cl', 'ai', 'cl_norm', 'cl', 'cd', 'cdv', 'cm_c/4', 'cm_LE', 'C.P.x/c']
    
    # Divisor principal para separar as seções de superfície
    surfaces = file_content.split("---------------------------------------------------------------")
    
    # Ignorar o primeiro elemento vazio antes da primeira linha de separação
    surfaces = surfaces[1:] 
    
    for i, surface_block in enumerate(surfaces):
        # 1. Identificar o nome da superfície
        # Procurar pela linha que começa com "Surface # X     [Nome da Superfície]"
        surface_name = "Unknown"
        if "Surface #" in surface_block:
            # Encontrar a linha e extrair o nome após o número
            name_line = [line for line in surface_block.split('\n') if "Surface #" in line][0]
            # Ex: "Surface # 1     Asa Principal" -> "Asa Principal"
            surface_name = name_line.split("Surface #")[1].split("     ")[-1].strip()
            # Limpar caracteres de parênteses e vírgulas para nomear o arquivo
            safe_name = surface_name.replace(" (YDUP)", "_YDUP").replace(" ", "_").replace(".", "").replace(",", "")
        
        # 2. Encontrar o bloco de dados de Strip Forces
        start_index = surface_block.find(HEADER_LINE)
        
        if start_index == -1:
            # A superfície não contém dados de Strip Forces
            continue

        data_block_start = surface_block[start_index:]
        
    
        data_io = io.StringIO(data_block_start)
        
        df = pd.read_csv(
            data_io, 
            sep='\s+',  # Usa regex para separar por um ou mais espaços (mais robusto)
            names=column_names, 
            skiprows=[0, 1],
            engine='python' # O uso de 'sep' com regex exige o engine Python
        )
        
        # Remover linhas potencialmente vazias ou de rodapé
        df = df.dropna(subset=['Yle'])
        
        # 4. Salvar em um arquivo CSV
        output_filename = f"PRJ-23\\Resultados\\Estabilidade\\{output_prefix}_{i+1}_{safe_name}.csv"
        df.to_csv(output_filename, index=False)


convert_forces_txt_to_csv('PRJ-23\\fs_aft_nottrim.txt')