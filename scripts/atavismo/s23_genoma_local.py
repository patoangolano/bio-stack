"""tBLASTn genoma-inteiro de Gallus gallus, LOCAL e controlado.
Baixa GRCg7b uma vez e roda diamond blastx. O FASTA nunca e lido, so alimenta o diamond."""

from __future__ import annotations

import subprocess
import zipfile
from pathlib import Path

import requests

SCR = Path(
    r"C:\Users\Manec\AppData\Local\Temp\claude\C--Users-Manec-bio-stack\e2209347-cc77-4f1d-8a19-fcf7abdb21d6\scratchpad"
)
SCR.mkdir(parents=True, exist_ok=True)
SEQ = Path(r"C:\Users\Manec\bio-stack\runs\atavismo\seqs")
DMD = Path(r"C:\Users\Manec\bio-stack\.venv\Lib\site-packages\gget\bins\Windows\diamond.exe")
zp = SCR / "gallus.zip"
fa = SCR / "gallus_GRCg7b.fna"
if not fa.exists():
    if not zp.exists() or zp.stat().st_size < 1e8:
        url = (
            "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
            "GCF_016699485.2/download?include_annotation_type=GENOME_FASTA"
        )
        print("baixando GRCg7b...", flush=True)
        with requests.get(url, stream=True, timeout=3600) as r:
            r.raise_for_status()
            n = 0
            with open(zp, "wb") as f:
                for ch in r.iter_content(1 << 22):
                    f.write(ch)
                    n += len(ch)
                    if n % (1 << 28) < (1 << 22):
                        print(f"  {n / 1e6:.0f} MB", flush=True)
        print(f"zip: {zp.stat().st_size / 1e6:.0f} MB", flush=True)
    with zipfile.ZipFile(zp) as z:
        names = [n for n in z.namelist() if n.endswith(("_genomic.fna", ".fna"))]
        print("extraindo:", names, flush=True)
        with open(fa, "wb") as out:
            for n in names:
                with z.open(n) as src:
                    while True:
                        b = src.read(1 << 22)
                        if not b:
                            break
                        out.write(b)
print(f"FASTA: {fa.stat().st_size / 1e6:.0f} MB", flush=True)
cmd = [
    str(DMD),
    "blastx",
    "-q",
    str(fa),
    "-d",
    str(SEQ / "ref"),
    "--sensitive",
    "--masking",
    "0",
    "--comp-based-stats",
    "0",
    "--frameshift",
    "15",
    "--evalue",
    "1e-3",
    "--max-target-seqs",
    "200",
    "--threads",
    "16",
    "--block-size",
    "4",
    "--quiet",
    "--outfmt",
    "6",
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qframe",
    "-o",
    str(SEQ / "hits_GENOMA_INTEIRO.tsv"),
]
print("rodando diamond blastx no genoma inteiro...", flush=True)
subprocess.run(cmd, check=True)
print("diamond OK", flush=True)
