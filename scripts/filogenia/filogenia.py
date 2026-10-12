"""
Roda o pipeline de filogenia (pipelines/filogenia) a partir do Windows.

    uv run python scripts/filogenia/filogenia.py --input pipelines/filogenia/assets/test/locos --outgroup Homo_sapiens

O Nextflow so roda no WSL, e chamar o WSL de fora tem tres armadilhas medidas
nesta maquina (CLAUDE.md): o wsl.exe engole `$var` da linha de comando, o Git
Bash reescreve caminho POSIX, e o PATH do nextflow so existe em shell de login.
Este script contorna duas delas: escreve o comando num arquivo
(runs/filogenia/<id>/executar.sh, que fica como registro exato do que rodou),
entao nenhum `$var` passa pela linha do wsl.exe, e executa com `bash -l`. A do
Git Bash nao alcanca a chamada que ele faz ao wsl.exe (quem chama e o Python,
nao o MSYS), mas alcanca os argumentos dele: rodado do Git Bash, caminho do
WSL em --input ou nos extras (/home/..., /mnt/c/...) chega como
C:/Program Files/Git/... Nesse caso exporte MSYS_NO_PATHCONV=1 antes, e ponha
~ entre aspas: sem aspas o bash troca o til pelo HOME do Windows.

Layout de uma execucao:
  ~/nf-work/filogenia/<id>/      diretorio de lancamento (.nextflow/, logs) — no WSL
  ~/nf-work/filogenia/work/      work dir compartilhado (cache do -resume) — no WSL
  runs/filogenia/<id>/           resultados, figuras, executar.sh — no projeto

Ao terminar, desenha as arvores (desenhar_arvore.py) e grava uma linha de
proveniencia em runs/provenance.jsonl com parametros, versoes e saidas.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
import re
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

RAIZ = Path(__file__).resolve().parents[2]
PIPELINE = RAIZ / "pipelines" / "filogenia"
RUNS = RAIZ / "runs" / "filogenia"
DISTRO = "Ubuntu-24.04"


_UNC_WSL = re.compile(r"^[\\/]{2}wsl(?:\.localhost|\$)[\\/]([^\\/]+)(.*)$", re.IGNORECASE)
_GIT_BASH = re.compile(r"^/([A-Za-z])(/.*)?$")


def para_wsl(caminho: str | Path) -> str:
    r"""
    Caminho do Windows -> caminho do WSL. Aceita C:\..., /c/... (Git Bash com
    MSYS_NO_PATHCONV=1) e \\wsl.localhost\<distro>\... (Explorer); caminho que
    ja e do WSL (/home/..., /mnt/c/..., ~/...) passa intacto.
    """
    s = str(caminho)
    if m := _UNC_WSL.match(s):
        if m[1].lower() != DISTRO.lower():
            raise SystemExit(f"[filogenia] {s} esta na distro {m[1]}; o pipeline roda em {DISTRO}")
        return m[2].replace("\\", "/") or "/"
    if m := _GIT_BASH.match(s):
        return f"/mnt/{m[1].lower()}{m[2] or ''}"
    if s.startswith(("/", "~")):
        return s
    s = str(Path(s).resolve()).replace("\\", "/")
    if len(s) > 1 and s[1] == ":":
        return f"/mnt/{s[0].lower()}{s[2:]}"
    return s


def _aspas(s: str) -> str:
    """
    Aspas simples de shell, para caminhos com espaco ou caractere especial. O til
    inicial sai das aspas como "$HOME": entre aspas o bash nao o expande, e o
    file() do Nextflow tambem nao (medido: file('~/x') vira <lancamento>/~/x).
    """
    if s == "~" or s.startswith("~/"):
        return '"$HOME"' + (_aspas(s[1:]) if len(s) > 1 else "")
    return "'" + s.replace("'", "'\"'\"'") + "'"


def preparar_entrada(entrada: str, saida: Path) -> str:
    """
    CSV loco,fasta com caminho do Windows nao abre no WSL. Reescreve o CSV em
    runs/filogenia/<id>/entrada.csv com cada fasta ja como caminho do WSL — o
    relativo resolvido contra o diretorio do CSV, como o pipeline faria. Outra
    entrada (diretorio, FASTA, glob) so passa por para_wsl.
    """
    origem = Path(entrada)
    if not entrada.lower().endswith(".csv") or not origem.is_file():
        return para_wsl(entrada)
    with origem.open(encoding="utf-8-sig", newline="") as fh:
        leitor = csv.DictReader(fh)
        campos = [c.strip() for c in leitor.fieldnames or []]
        linhas = [{k.strip(): v for k, v in linha.items() if k} for linha in leitor]
    if "fasta" not in campos:
        return para_wsl(entrada)  # o pipeline recusa com a mensagem certa
    for linha in linhas:
        fasta = (linha.get("fasta") or "").strip()
        if not fasta:
            continue
        if fasta.startswith(("/", "~")) or Path(fasta).is_absolute():
            linha["fasta"] = para_wsl(fasta)
        else:
            linha["fasta"] = para_wsl(origem.parent / fasta)
    destino = saida / "entrada.csv"
    with destino.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, lineterminator="\n")
        w.writeheader()
        w.writerows(linhas)
    return para_wsl(destino)


def saneia_outgroup(outgroup: str) -> str:
    """
    O pipeline saneia o outgroup como saneia os taxons ('Homo sapiens' vira
    'Homo_sapiens' na arvore). A figura tem de procurar o mesmo nome, senao
    enraiza no ponto medio. Usa a funcao do proprio pipeline, nao uma copia.
    """
    spec = importlib.util.spec_from_file_location(
        "filogenia_util", PIPELINE / "bin" / "filogenia_util.py"
    )
    if spec is None or spec.loader is None:
        return outgroup
    util = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(util)
    nomes = (util.sanear_nome(o.strip(), "completo") for o in outgroup.split(","))
    return ",".join(n for n in nomes if n)


def montar_script(a: argparse.Namespace, run_id: str, saida: Path) -> str:
    entrada = preparar_entrada(a.input, saida)
    params = {
        "input": entrada,
        "outdir": para_wsl(saida),
        "alinhador": a.alinhador,
        "limpeza": a.limpeza,
        "outgroup": a.outgroup,
        "tipo_seq": a.tipo_seq,
        "nome_taxon": a.nome_taxon,
        "cpus_concatenado": a.cpus_concatenado,
        "seed": a.seed,
    }
    args = " \\\n    ".join(
        f"--{k} {_aspas(str(v))}" for k, v in params.items() if v not in (None, "")
    )
    extras = " ".join(_aspas(x) for x in a.extras)
    lancamento = f"~/nf-work/filogenia/{run_id}"
    return (
        "#!/usr/bin/env bash\n"
        f"# bio-stack/filogenia — execucao {run_id}\n"
        f"# gerado por scripts/filogenia/filogenia.py em {datetime.now(UTC).isoformat(timespec='seconds')}\n"
        "set -euo pipefail\n"
        "export NXF_ANSI_LOG=false\n"
        f"mkdir -p {lancamento} && cd {lancamento}\n"
        f"nextflow run {_aspas(para_wsl(PIPELINE))} \\\n"
        f"    -profile {a.perfil} \\\n"
        "    -w ~/nf-work/filogenia/work \\\n"
        "    -resume \\\n"
        f"    {args}" + (f" \\\n    {extras}" if extras else "") + "\n"
    )


def registrar_proveniencia(
    run_id: str, saida: Path, script: str, ok: bool, segundos: float
) -> None:
    resumo: dict = {}
    if (saida / "resumo.json").is_file():
        resumo = json.loads((saida / "resumo.json").read_text(encoding="utf-8"))
    entrada = {
        "v": 3,
        "ts": datetime.now(UTC).isoformat(timespec="seconds"),
        "server": "filogenia",
        "kind": "pipeline-run",
        "tool": "pipelines/filogenia",
        "ok": ok,
        "run_id": run_id,
        "segundos": round(segundos),
        "comando": script,
        "params": resumo.get("parametros", {}),
        "versoes": resumo.get("versoes", {}),
        # as_posix: o log e lido no WSL e por /metodos; barra invertida do Windows
        # viraria escape em qualquer ferramenta que trate o campo como caminho.
        "outputs": [p.relative_to(RAIZ).as_posix() for p in sorted(saida.glob("*.*"))],
    }
    log = Path(os.environ.get("BIO_STACK_PROVENANCE", RAIZ / "runs" / "provenance.jsonl"))
    try:
        with log.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(entrada, ensure_ascii=False) + "\n")
    except OSError as exc:
        print(f"[filogenia] AVISO: proveniencia nao gravada em {log}: {exc}", file=sys.stderr)


def desenhar(saida: Path, outgroup: str) -> list[Path]:
    """Figuras com o mesmo Python do projeto (Bio.Phylo + matplotlib)."""
    from desenhar_arvore import ROTULOS_PADRAO, desenhar_arvore  # mesmo diretorio deste script

    # A ordem dos campos do rotulo depende de quais suportes a execucao pediu;
    # o resumo.json registra. Sem ele, os padroes do pipeline.
    rotulos = ROTULOS_PADRAO
    if (saida / "resumo.json").is_file():
        gravados = json.loads((saida / "resumo.json").read_text(encoding="utf-8")).get(
            "rotulos_suporte"
        )
        if gravados:
            rotulos = tuple(gravados)

    figuras = saida / "figuras"
    outgroup = saneia_outgroup(outgroup)
    feitas: list[Path] = []
    # (arvore, titulo, unidade do eixo, legenda). O wASTRAL hibrido escreve
    # comprimento so onde consegue, e esse numero nao e unidade coalescente
    # (mistura comprimento de gene com o do ASTRAL) — a figura sai sem escala.
    alvos = [
        (
            saida / "05_concatenado" / "concatenado.treefile",
            "Supermatriz particionada (ML)",
            "substituições por sítio",
            "supermatriz",
        ),
        (
            saida / "07_arvore_especies" / "especies.tre",
            "Árvore de espécies (wASTRAL)",
            "comprimento do wASTRAL (não é unidade coalescente)",
            "especies",
        ),
    ]
    alvos += [
        (
            arq,
            # "loco_C" ja diz loco; "Loco loco_C" repetiria a palavra.
            arq.stem if arq.stem.lower().startswith("loco") else f"Loco {arq.stem}",
            "substituições por sítio",
            "genes",
        )
        for arq in sorted((saida / "04_arvores_genes").glob("*.treefile"))
    ]
    for arvore, titulo, unidade, tipo in alvos:
        if arvore.is_file():
            feitas += desenhar_arvore(
                arvore,
                figuras / arvore.stem,
                titulo=titulo,
                outgroup=outgroup,
                unidade=unidade,
                rotulos=rotulos,
                tipo=tipo,
            )
    return feitas


def imprimir_resumo(saida: Path) -> None:
    caminho = saida / "resumo.json"
    if not caminho.is_file():
        return
    r = json.loads(caminho.read_text(encoding="utf-8"))
    print(f"\nLocos: {r['n_locos_inferidos']} inferidos de {r['n_locos_entrada']}", end="")
    print(f" (fora: {r['locos_fora']})" if r.get("locos_fora") else "")
    for loco in r.get("locos", []):
        print(
            f"  {loco['loco']:<24} {loco.get('n_seqs', '?'):>4} seqs  {loco.get('n_sitios', '?'):>6} sitios  "
            f"{loco.get('modelo', '?'):<14} UFBoot>=95: {loco.get('frac_ufboot_95', '?')}"
        )
    if c := r.get("concatenado"):
        particoes = f"{c.get('n_particoes_entrada')} particoes"
        if c.get("n_particoes_finais") is not None:
            particoes += f" -> {c['n_particoes_finais']} apos a fusao (MERGE)"
        print(f"Supermatriz: {c.get('n_sitios')} sitios, {particoes}; lnL {c.get('lnL')}")
        # Com muitos sitios o UFBoot da supermatriz satura em 100 mesmo em ramo
        # que metade dos genes contradiz; quem julga o ramo e a concordancia.
        print(
            f"  suporte satura em supermatriz (julgue pelo gCF/sCFL abaixo); "
            f"ramos com SH-aLRT>=80 e UFBoot>=95: {c.get('frac_sh80_ufb95')}"
        )
    if cf := r.get("concordancia"):
        print(
            f"Concordancia: gCF medio {cf.get('gCF_medio')} (min {cf.get('gCF_min')}), sCFL medio {cf.get('sCF_medio')}"
        )
    # Aviso do preparo (stop interno removido, registro vazio descartado) muda
    # o dado; tem de aparecer aqui, nao so no stderr da tarefa.
    for loco, avisos in (r.get("avisos_preparo") or {}).items():
        for aviso in avisos:
            print(f"  AVISO {loco}: {aviso}")
    print(f"Metodos: {saida / 'metodos.md'}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--input",
        required=True,
        help="diretorio (um FASTA por loco), FASTA unico, glob ou CSV loco,fasta",
    )
    ap.add_argument("--outgroup", default="", help="taxon(s) para enraizar, separados por virgula")
    ap.add_argument("--nome", default="", help="rotulo da execucao (padrao: nome da entrada)")
    ap.add_argument("--alinhador", default="mafft", choices=["mafft", "famsa", "nenhum"])
    ap.add_argument("--limpeza", default="clipkit", choices=["clipkit", "trimal", "nenhuma"])
    ap.add_argument("--tipo-seq", default="", choices=["", "DNA", "AA"])
    ap.add_argument(
        "--nome-taxon", default="primeiro", choices=["primeiro", "completo", "colchete"]
    )
    ap.add_argument(
        "--cpus-concatenado",
        type=int,
        default=1,
        help=">1 acelera a supermatriz, mas perde reproducao bit a bit",
    )
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument(
        "--perfil", default="docker", help="perfil do Nextflow (docker; 'test,docker' roda o teste)"
    )
    ap.add_argument(
        "--retomar", default="", help="id de uma execucao anterior para continuar com -resume"
    )
    ap.add_argument("--sem-figura", action="store_true")
    ap.add_argument(
        "extras", nargs="*", help="argumentos extras repassados ao nextflow (depois de --)"
    )
    a = ap.parse_args(argv)

    if a.retomar:
        run_id = a.retomar
    else:
        base = a.nome or Path(str(a.input).rstrip("/\\")).stem or "filogenia"
        rotulo = re.sub(r"[^A-Za-z0-9_-]+", "_", base).strip("_")[:40]
        run_id = f"{datetime.now().strftime('%Y%m%d-%H%M%S')}-{rotulo}"
    saida = RUNS / run_id
    saida.mkdir(parents=True, exist_ok=True)

    script = montar_script(a, run_id, saida)
    arq_script = saida / "executar.sh"
    arq_script.write_text(script, encoding="utf-8", newline="\n")
    print(
        f"[filogenia] execucao {run_id}\n[filogenia] comando gravado em {arq_script.relative_to(RAIZ)}"
    )

    inicio = datetime.now(UTC)
    with (saida / "nextflow_console.log").open("w", encoding="utf-8") as console:
        proc = subprocess.Popen(
            ["wsl", "-d", DISTRO, "--", "bash", "-l", para_wsl(arq_script)],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            env={**os.environ, "MSYS_NO_PATHCONV": "1"},
        )
        assert proc.stdout is not None
        for linha in proc.stdout:
            sys.stdout.write(linha)
            console.write(linha)
        codigo = proc.wait()
    segundos = (datetime.now(UTC) - inicio).total_seconds()

    if codigo == 0 and not a.sem_figura:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        try:
            for fig in desenhar(saida, a.outgroup):
                print(f"[filogenia] figura: {fig.relative_to(RAIZ)}")
        except ImportError as exc:
            print(f"[filogenia] figuras puladas: {exc} (rode com `uv run`)", file=sys.stderr)

    registrar_proveniencia(run_id, saida, script, codigo == 0, segundos)
    if codigo == 0:
        imprimir_resumo(saida)
    else:
        print(
            f"[filogenia] falhou (exit {codigo}). Log: {saida / 'nextflow_console.log'}",
            file=sys.stderr,
        )
    return codigo


if __name__ == "__main__":
    sys.exit(main())
