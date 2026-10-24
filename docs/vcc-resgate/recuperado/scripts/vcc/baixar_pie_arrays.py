"""
Baixa os arrays de resposta dos repositorios PIE que alimentam o transfer_fit
(etapa 3 de docs/VCC-2026-ROTA.md).

Cada repo preprocessado da Arc guarda uma matriz (contexto, perturbacao) x gene
por tipo de evidencia. Baixamos cinco das seis:

    fold_changes   log2FC por gene              float32   alvo de magnitude
    delta_p        delta discretizado do PIE    float32   alvo auxiliar
    fdr            q-valor por gene             float32   alvo de DE
    tested         mascara de gene testado      bool      universo por linha
    ctrl_means     expressao media do controle  float32   features de contexto

`lfc_true` fica de fora de proposito: 2,36 GB em float64 so no Orion, e o alvo
de magnitude ja sai de fold_changes em float32 pela metade do tamanho. Os
parquet de `de/` tambem ficam de fora: sao a mesma informacao por linha, em
formato de tabela longa, e somam 11,4 GB so no Orion.

Tahoe-100M nao entra: e perturbacao por droga e dose, nao knockout genico, e o
eixo de perturbacao nao cruza com os 300 alvos de CRISPRi do desafio.

Retomavel. `hf_hub_download` continua um arquivo parcial e nao rebaixa o que ja
esta completo, entao rodar de novo depois de uma queda de rede e barato.

Rodar da raiz do projeto:
    ./.venv/Scripts/python.exe -I scripts/vcc/baixar_pie_arrays.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

from huggingface_hub import HfApi, hf_hub_download
from huggingface_hub.errors import HfHubHTTPError

DESTINO = Path("data/raw/pie/arrays")
REPOS = (
    "PIE_x_atlas_orion",
    "PIE_replogle_nadig_essential",
    "PIE_arc_vcc_25",
    "PIE_jiang",
)
ARQUIVOS = (
    "preprocessed/meta.json",
    "preprocessed/ctx_ids.npy",
    "preprocessed/pert_ids.npy",
    "preprocessed/ctrl_means.npy",
    "preprocessed/tested.npy",
    "preprocessed/fold_changes.npy",
    "preprocessed/fdr.npy",
    "preprocessed/delta_p.npy",
)


def tamanhos(api: HfApi, repo_id: str) -> dict[str, int]:
    """Tamanho em bytes de cada arquivo do repo, pelo metadado do Hub.

    `dataset_info` em vez de `repo_info`: o segundo devolve uma uniao que
    inclui KernelInfo, que nao tem `siblings`.
    """
    info = api.dataset_info(repo_id, files_metadata=True)
    return {s.rfilename: (s.size or 0) for s in (info.siblings or [])}


def main() -> None:
    api = HfApi()
    DESTINO.mkdir(parents=True, exist_ok=True)
    plano: list[tuple[str, str, int]] = []
    ausentes: list[tuple[str, str]] = []

    for repo in REPOS:
        repo_id = f"arcinstitute/{repo}"
        try:
            disponiveis = tamanhos(api, repo_id)
        except HfHubHTTPError as erro:
            print(f"!! {repo}: nao foi possivel listar ({erro})", file=sys.stderr)
            raise
        for arquivo in ARQUIVOS:
            if arquivo in disponiveis:
                plano.append((repo, arquivo, disponiveis[arquivo]))
            else:
                # arc_vcc_25 nao publica ctrl_means: o repo da competicao traz
                # so as linhas perturbadas. Ausencia e fato do repo, nao erro.
                ausentes.append((repo, arquivo))

    total = sum(tam for _, _, tam in plano)
    print(f"== plano: {len(plano)} arquivos, {total / 1e9:.2f} GB ==")
    for repo, arquivo, tam in plano:
        print(f"   {tam / 1e6:9.1f} MB  {repo}/{arquivo}")
    if ausentes:
        print("\n== nao publicados pelo repo (seguimos sem eles) ==")
        for repo, arquivo in ausentes:
            print(f"   {repo}/{arquivo}")

    print(f"\n== baixando para {DESTINO} ==", flush=True)
    baixados = 0
    for indice, (repo, arquivo, tam) in enumerate(plano, start=1):
        alvo = DESTINO / repo
        alvo.mkdir(parents=True, exist_ok=True)
        inicio = time.monotonic()
        caminho = hf_hub_download(
            repo_id=f"arcinstitute/{repo}",
            filename=arquivo,
            repo_type="dataset",
            local_dir=str(alvo),
        )
        decorrido = time.monotonic() - inicio
        real = Path(caminho).stat().st_size
        if real != tam:
            raise OSError(f"{repo}/{arquivo}: {real} bytes no disco, {tam} anunciados pelo Hub")
        baixados += real
        taxa = (real / 1e6 / decorrido) if decorrido > 0.05 else float("inf")
        print(
            f"[{indice:2d}/{len(plano)}] ok {real / 1e6:9.1f} MB  {taxa:6.1f} MB/s  {repo}/{arquivo}",
            flush=True,
        )

    registro = {
        "destino": str(DESTINO),
        "repos": list(REPOS),
        "arquivos_por_repo": list(ARQUIVOS),
        "nao_publicados": [f"{r}/{a}" for r, a in ausentes],
        "bytes_baixados": baixados,
        "lfc_true_excluido": "float64, 2,36 GB no Orion; magnitude vem de fold_changes",
        "tahoe_excluido": "droga-dose, nao knockout genico",
    }
    (DESTINO / "registro_download.json").write_text(
        json.dumps(registro, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"\n== concluido: {baixados / 1e9:.2f} GB ==")


if __name__ == "__main__":
    main()
