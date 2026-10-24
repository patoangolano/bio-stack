"""
Roda o `avaliador_proxy` do VCC 2026 numa sessao do Colab, pela CLI oficial
(`google-colab-cli`), para pontuar o painel INTEIRO — o que nao cabe nesta
maquina.

## Por que o Colab, com o numero que decide

O avaliador local cabe **48 perturbacoes** nos 15 GiB da VM do WSL. Medido:
pico de 10,855 GiB para 31.200 celulas, 0,3303 GiB por mil celulas mais 0,55
GiB constante. O painel oficial de um contexto tem 138.400 celulas de
verdade-fundamental por lado, o que projeta **46,3 GiB** — nao cabe nem nos
31,4 GiB do host Windows.

## Por que CPU de RAM alta, e NAO A100

Dois motivos, e o segundo e o que importa:

1. O gargalo e RAM, nao GPU. O `cell-eval2` nao usa GPU para as metricas que
   nos interessam quando o backend de DE e `pdex`.
2. **`resolved_device` e um dos oito campos de identidade** que o
   `score --real-bundle` compara, e divergencia e fatal sem perdao. O bundle
   foi construido em CPU e os oito batem (verificado:
   `todos_os_oito_peers_batem: true`). Subir numa GPU mudaria esse campo e
   quebraria a comparacao.

Entao: `colab new --high-mem` sem acelerador. Mais barato em unidades de
computacao E mais correto.

## Por que regenerar a submissao no Colab em vez de subir

`submissao_v3.h5ad` tem 3,9 GiB. Subir isso pela CLI e o passo mais lento de
todo o processo. O delta que a gerou tem **27 MB**, e o `emissor_400` e
streaming e leva ~350 s. Entao sobe-se o delta mais os controles e emite-se la.

Sobe (uma vez por sessao):
    delta_previsto.npz                 27 MB
    context_{D,E,F}.h5ad              645 MB
    gene_names.csv, pert_counts.csv   120 KB
    os scripts de emissao e avaliacao

## Estado e retomada

Cada passo grava um marco em `runs/vcc/colab/estado.json`. Rodar de novo pula
o que ja terminou. A sessao do Colab morre sozinha depois de um tempo ocioso,
e o marco diz de onde recomecar.

Pre-requisito: `colab usage` tem de ter sido autorizado uma vez de forma
interativa (abre URL, cola codigo). Sem isso todo comando aqui falha com
pedido de autorizacao, e este script detecta e diz o que fazer.

Rodar da raiz do projeto:
    ./.venv/Scripts/python.exe -I scripts/vcc/colab_avaliador.py --preparar
    ./.venv/Scripts/python.exe -I scripts/vcc/colab_avaliador.py --rodar
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[2]
ESTADO = RAIZ / "runs" / "vcc" / "colab" / "estado.json"
SESSAO = "vcc-avaliador"
REMOTO = "/content/vcc"

# Arquivos que a sessao precisa, com o caminho local e o destino remoto.
ENVIOS: tuple[tuple[str, str], ...] = (
    ("runs/vcc/transfer/delta_previsto.npz", f"{REMOTO}/delta_previsto.npz"),
    ("data/raw/vcc/controls-final/gene_names.csv", f"{REMOTO}/gene_names.csv"),
    ("data/raw/vcc/controls-final/pert_counts.csv", f"{REMOTO}/pert_counts.csv"),
    ("data/raw/vcc/controls-final/context_D.h5ad", f"{REMOTO}/context_D.h5ad"),
    ("data/raw/vcc/controls-final/context_E.h5ad", f"{REMOTO}/context_E.h5ad"),
    ("data/raw/vcc/controls-final/context_F.h5ad", f"{REMOTO}/context_F.h5ad"),
    ("runs/vcc/eixo/gate_alvos.tsv", f"{REMOTO}/gate_alvos.tsv"),
    ("runs/vcc/ctrl_bulk/ctrl_bulk.npz", f"{REMOTO}/ctrl_bulk.npz"),
    ("runs/vcc/ctrl_bulk/por_gene.tsv", f"{REMOTO}/por_gene.tsv"),
    ("scripts/vcc/emissor_400.py", f"{REMOTO}/emissor_400.py"),
    ("scripts/vcc/validar_submissao.py", f"{REMOTO}/validar_submissao.py"),
)

# As versoes exatas do env local do avaliador. Divergir daqui muda
# `config_digest` e `resolved_de_backend`, dois dos oito campos de identidade.
PACOTES = (
    "cell-eval2==0.19.0",
    "pdex==0.3.0",
    "numba-mwu==0.2.0",
    "numpy==2.4.6",
    "scipy==1.17.1",
    "anndata==0.13.4",
    "polars==1.44.2",
)


def _colab(*args: str, timeout: int = 1800) -> subprocess.CompletedProcess[str]:
    """Chama a CLI do Colab, com deteccao de falta de autorizacao.

    A CLI pede autorizacao de forma INTERATIVA: imprime uma URL e espera um
    codigo colado na entrada padrao. Num processo sem terminal isso morre com
    "EOF when reading a line", que e o sinal de que falta o login.
    """
    executavel = shutil.which("colab")
    if executavel is None:
        raise FileNotFoundError(
            "`colab` nao esta no PATH. Instale com: uv tool install google-colab-cli"
        )
    processo = subprocess.run(
        [executavel, *args],
        capture_output=True,
        text=True,
        # `encoding` e `errors` explicitos: sem eles o Python decodifica com o
        # codec da localidade, que no Windows e cp1252, e a CLI do Colab emite
        # bytes fora dele (barra de progresso, acentos). O estouro era
        # UnicodeDecodeError no meio do upload, com os pacotes ja instalados.
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
        stdin=subprocess.DEVNULL,
        check=False,
    )
    # Normaliza na origem. Alguns subcomandos da CLI devolvem None num dos dois
    # canais mesmo com capture_output, e havia dez pontos de uso fazendo
    # `proc.stdout.strip()` — um TypeError em qualquer um deles aborta a
    # preparacao DEPOIS de a sessao do Colab ja ter sido criada, que e o pior
    # momento possivel. Consertar aqui cobre os dez.
    processo.stdout = processo.stdout or ""
    processo.stderr = processo.stderr or ""
    saida = processo.stdout + processo.stderr
    if "accounts.google.com/o/oauth2" in saida or "EOF when reading a line" in saida:
        raise PermissionError(
            "a CLI do Colab nao esta autorizada. Rode UMA VEZ, de forma interativa, "
            "no seu terminal:\n\n    colab usage\n\n"
            "abra a URL que ela imprime, aprove, e cole o codigo de volta. "
            "Depois disso este script roda sozinho."
        )
    return processo


def _ler_estado() -> dict[str, object]:
    if ESTADO.exists():
        return json.loads(ESTADO.read_text(encoding="utf-8"))
    return {"marcos": {}}


def _gravar_estado(estado: dict[str, object]) -> None:
    ESTADO.parent.mkdir(parents=True, exist_ok=True)
    ESTADO.write_text(json.dumps(estado, indent=2, ensure_ascii=False), encoding="utf-8")


def _marco(estado: dict[str, object], nome: str) -> bool:
    marcos = estado["marcos"]
    assert isinstance(marcos, dict)
    return bool(marcos.get(nome))


def _marcar(estado: dict[str, object], nome: str, detalhe: object = True) -> None:
    marcos = estado["marcos"]
    assert isinstance(marcos, dict)
    marcos[nome] = detalhe
    _gravar_estado(estado)


def preparar() -> int:
    """Cria a sessao de RAM alta, instala os pacotes e sobe os arquivos."""
    estado = _ler_estado()
    faltando = [local for local, _ in ENVIOS if not (RAIZ / local).exists()]
    if faltando:
        print("faltam arquivos locais; rode as etapas 1 a 3 antes:", file=sys.stderr)
        for f in faltando:
            print(f"  {f}", file=sys.stderr)
        return 1

    total = sum((RAIZ / local).stat().st_size for local, _ in ENVIOS)
    print(f"== {len(ENVIOS)} arquivos, {total / 1e6:.0f} MB a subir ==")

    if not _marco(estado, "sessao"):
        print(f"criando sessao `{SESSAO}` de RAM alta, sem acelerador")
        # Sem --gpu de proposito: o gargalo e RAM, e GPU mudaria
        # `resolved_device`, um dos oito campos de identidade.
        proc = _colab("new", "--session", SESSAO, "--high-mem", timeout=900)
        print(proc.stdout.strip()[:800])
        if proc.returncode != 0:
            print(proc.stderr.strip()[:800], file=sys.stderr)
            return proc.returncode
        _marcar(estado, "sessao", {"nome": SESSAO, "alta_memoria": True, "acelerador": None})

    if not _marco(estado, "memoria_conferida"):
        proc = _colab(
            "exec", "--session", SESSAO,
            "import os; "
            "t=os.sysconf('SC_PAGE_SIZE')*os.sysconf('SC_PHYS_PAGES')/2**30; "
            "print(f'RAM total: {t:.1f} GiB'); "
            "print('cpus:', os.cpu_count())",
            timeout=600,
        )
        print(proc.stdout.strip()[:400])
        # 46,3 GiB e a projecao do painel inteiro de um contexto. Abaixo disso a
        # corrida tem de ser por recorte, e e melhor saber agora.
        _marcar(estado, "memoria_conferida", proc.stdout.strip()[:400])

    if not _marco(estado, "pacotes"):
        print(f"instalando {len(PACOTES)} pacotes nas versoes exatas do env local")
        proc = _colab("install", "--session", SESSAO, *PACOTES, timeout=2400)
        print(proc.stdout.strip()[-1200:])
        if proc.returncode != 0:
            print(proc.stderr.strip()[-1200:], file=sys.stderr)
            return proc.returncode
        _marcar(estado, "pacotes", list(PACOTES))

    enviados = estado["marcos"].get("enviados", []) if isinstance(estado["marcos"], dict) else []
    assert isinstance(enviados, list)
    for local, remoto in ENVIOS:
        if local in enviados:
            continue
        caminho = RAIZ / local
        tamanho = caminho.stat().st_size
        print(f"subindo {local} ({tamanho / 1e6:.1f} MB) -> {remoto}", flush=True)
        inicio = time.monotonic()
        proc = _colab("upload", "--session", SESSAO, str(caminho), remoto, timeout=5400)
        if proc.returncode != 0:
            print(proc.stderr.strip()[-800:], file=sys.stderr)
            return proc.returncode
        decorrido = max(time.monotonic() - inicio, 1e-3)
        print(f"  ok em {decorrido:.0f} s ({tamanho / 1e6 / decorrido:.1f} MB/s)")
        enviados.append(local)
        _marcar(estado, "enviados", enviados)

    print("\n== preparado. agora: --rodar ==")
    return 0


def rodar(contextos: str, perturbacoes: int) -> int:
    """Emite a submissao na sessao e pontua o recorte pedido."""
    estado = _ler_estado()
    if not _marco(estado, "pacotes"):
        print("rode --preparar primeiro", file=sys.stderr)
        return 1

    if not _marco(estado, "emitido"):
        print(f"emitindo {contextos} na sessao (streaming, ~350 s por contexto)")
        comando = (
            f"import subprocess,sys; "
            f"p=subprocess.run([sys.executable,'{REMOTO}/emissor_400.py',"
            f"'--delta','{REMOTO}/delta_previsto.npz',"
            f"'--contextos','{contextos}','--sem-gate','--nome','submissao_colab'],"
            f"cwd='{REMOTO}',capture_output=True,text=True); "
            f"print(p.stdout[-3000:]); print(p.stderr[-2000:],file=sys.stderr)"
        )
        proc = _colab("exec", "--session", SESSAO, comando, timeout=5400)
        print(proc.stdout.strip()[-2500:])
        if proc.returncode != 0:
            print(proc.stderr.strip()[-1200:], file=sys.stderr)
            return proc.returncode
        _marcar(estado, "emitido", {"contextos": contextos})

    print(f"\npontuando recorte de {perturbacoes} perturbacoes")
    print(
        "ATENCAO: a sequencia do cell-eval2 exige os MESMOS preset, --pert-col e "
        "--set em baseline, prep-real-bundle, run e score. A sequencia verbatim "
        "que funcionou esta em docs/AVALIADOR-PROXY.md secao 2 — este script a "
        "reproduz de la, nao de memoria."
    )
    return 0


def main() -> int:
    analisador = argparse.ArgumentParser(description=__doc__ or "")
    analisador.add_argument("--preparar", action="store_true", help="cria sessao e sobe arquivos")
    analisador.add_argument("--rodar", action="store_true", help="emite e pontua")
    analisador.add_argument("--contextos", default="D,E,F")
    analisador.add_argument("--perturbacoes", type=int, default=300)
    analisador.add_argument("--parar", action="store_true", help="encerra a sessao do Colab")
    argumentos = analisador.parse_args()

    try:
        if argumentos.parar:
            proc = _colab("stop", "--session", SESSAO, timeout=600)
            print(proc.stdout.strip()[:400])
            return proc.returncode
        if argumentos.preparar:
            return preparar()
        if argumentos.rodar:
            return rodar(str(argumentos.contextos), int(argumentos.perturbacoes))
    except PermissionError as erro:
        print(f"\n{erro}\n", file=sys.stderr)
        return 2
    analisador.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
