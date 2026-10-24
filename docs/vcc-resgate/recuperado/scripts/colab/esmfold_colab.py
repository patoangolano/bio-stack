# bio-stack — ESMFold na T4 do Colab, via mcp__colab-proxy-mcp__run_code_cell.
#
# Por que este arquivo existe: o endpoint `esmfold` da NVIDIA respondeu 410 em
# 2026-10-21 e saiu do catalogo, sem slug novo. E o ESMFold nao cabe na RTX
# 5070 Laptop: o pico medido e ~14 GiB em entrada curta e ~41 GiB a 1000
# residuos, contra 8151 MiB de VRAM. A T4 do Colab tem 16 GB e resolve a
# faixa curta e media.
#
# As duas linhas que fazem caber, e que sao o ponto do arquivo:
#   model.esm = model.esm.half()     -> fp16 no tronco de linguagem, que e a
#                                       maior parte dos parametros
#   model.trunk.set_chunk_size(64)   -> troca memoria por tempo no Evoformer
#
# Sem as duas, a T4 estoura em sequencia de tamanho util. Com elas, monomero
# de algumas centenas de residuos passa. O teste do NIH mediu que o
# --cpu-offload NAO salva as sequencias que falham no default, entao nao
# conte com ele: reduza o chunk ou corte a sequencia.
#
# Para AlphaFold2 (endpoint `alphafold2` tambem morto), nao reimplemente:
# use o notebook oficial ColabFold, que traz o MSA pelo servidor MMseqs2.
#   https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb

SEQUENCIA = (
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHS"
    "LAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATE"
)

CELULA_INSTALACAO = r"""
!pip -q install "transformers>=4.44" accelerate
import torch
print("GPU:", torch.cuda.get_device_name(0))
print("VRAM total (GiB):", round(torch.cuda.get_device_properties(0).total_memory / 1024**3, 2))
"""

CELULA_DOBRAMENTO = r"""
import torch
from transformers import AutoTokenizer, EsmForProteinFolding

tokenizador = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
modelo = EsmForProteinFolding.from_pretrained(
    "facebook/esmfold_v1", low_cpu_mem_usage=True
)
modelo = modelo.cuda()

# As duas alavancas de memoria. A ordem importa: half() antes do chunk.
modelo.esm = modelo.esm.half()
modelo.trunk.set_chunk_size(64)

sequencia = {sequencia!r}
entradas = tokenizador([sequencia], return_tensors="pt", add_special_tokens=False)
entradas = {{chave: valor.cuda() for chave, valor in entradas.items()}}

with torch.no_grad():
    saida = modelo(**entradas)

pdb = modelo.output_to_pdb(saida)[0]
with open("predicao.pdb", "w") as arquivo:
    arquivo.write(pdb)

# pLDDT medio: o numero que decide se a predicao vale uso. Acima de 70 o
# backbone e confiavel; abaixo de 50 a regiao provavelmente e desordenada.
plddt = saida["plddt"][0, :, 1].mean().item()
print("residuos:", len(sequencia))
print("pLDDT medio:", round(plddt, 2))
print("pico de VRAM (GiB):", round(torch.cuda.max_memory_allocated() / 1024**3, 2))
"""


def celulas(sequencia: str = SEQUENCIA) -> list[str]:
    """Devolve as duas celulas prontas para `run_code_cell`, na ordem."""
    return [
        CELULA_INSTALACAO,
        CELULA_DOBRAMENTO.format(sequencia=sequencia),
    ]


if __name__ == "__main__":
    for indice, celula in enumerate(celulas(), start=1):
        print(f"{'=' * 30} celula {indice} {'=' * 30}")
        print(celula)
