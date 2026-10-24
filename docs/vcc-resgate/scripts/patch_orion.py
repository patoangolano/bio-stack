import pathlib

p = pathlib.Path("C:/Users/Manec/bio-stack/scripts/vcc/avaliador_proxy.py")
s = p.read_text(encoding="utf-8")
PARES = []

PARES.append((
    '''RAIZ = Path(__file__).resolve().parents[2]
CONTROLES = RAIZ / "data" / "raw" / "vcc" / "controls-final"
PIE_ARC = RAIZ / "data" / "raw" / "pie" / "arrays" / "PIE_arc_vcc_25" / "preprocessed"

#: O rotulo e a coluna vem do `manifest.json` do painel, nao de convencao nossa.
COLUNA_PERT = "target_gene"
ROTULO_CONTROLE = "non-targeting"
#: Contexto de `arc_vcc_25` com os 190 nocautes medidos; `ARC_H1_VAL` e o fold de
#: validacao e nao traz linha propria em `fold_changes`.
CONTEXTO_PIE = "ARC_H1"''',
    '''RAIZ = Path(__file__).resolve().parents[2]
CONTROLES = RAIZ / "data" / "raw" / "vcc" / "controls-final"
PIE = RAIZ / "data" / "raw" / "pie" / "arrays"

#: O rotulo e a coluna vem do `manifest.json` do painel, nao de convencao nossa.
COLUNA_PERT = "target_gene"
ROTULO_CONTROLE = "non-targeting"

#: De onde vem a razao de dobra medida, e qual contexto de cada fonte usar.
#: A escolha decide quantas perturbacoes do PAINEL o proxy consegue rotular, e a
#: diferenca e grande — medido em 2026-10-22:
#:
#:   arc_vcc_25 / ARC_H1   146 nocautes,  12 deles entre os 300 alvos do VCC
#:   orion      / hct116  8 987 linhas, 148 deles entre os 300 alvos do VCC
#:
#: Os 12 do arc_vcc_25 sao o mesmo 12 que o memorando da rota tirou do
#: `train.json` do PIE; os 148 do orion sao o mesmo 148. Para um proxy com
#: rotulos REAIS do painel e mais de uma duzia de perturbacoes, orion e a unica
#: fonte que serve.
FONTES: dict[str, tuple[str, str]] = {
    # `ARC_H1_VAL` existe mas e o fold de validacao, e nao e o default.
    "arc_vcc_25": ("PIE_arc_vcc_25", "ARC_H1"),
    "orion": ("PIE_x_atlas_orion", "hct116"),
    "jiang": ("PIE_jiang", ""),
    "replogle": ("PIE_replogle_nadig_essential", ""),
}'''))

PARES.append((
    '''def _carregar_efeitos_pie(eixo_genes: np.ndarray, alvos_vcc: set[str], n_alvos: int,
                          *, min_genes_moventes: int, limiar_movente: float,
                          so_alvos_vcc: bool) -> Efeitos:''',
    '''def _carregar_efeitos_pie(eixo_genes: np.ndarray, alvos_vcc: set[str], n_alvos: int,
                          *, min_genes_moventes: int, limiar_movente: float,
                          so_alvos_vcc: bool, fonte_efeitos: str,
                          contexto_pie: str | None) -> Efeitos:'''))

PARES.append((
    '''    meta = json.loads((PIE_ARC / "meta.json").read_text(encoding="utf-8"))
    genes_pie = np.asarray(meta["genes"], dtype=str)
    pert_para_id = meta["pert_to_id"]
    ctx_id = meta["context_to_id"][CONTEXTO_PIE]

    fc = np.load(PIE_ARC / "fold_changes.npy")
    tested = np.load(PIE_ARC / "tested.npy")
    pert_ids = np.load(PIE_ARC / "pert_ids.npy")
    ctx_ids = np.load(PIE_ARC / "ctx_ids.npy")''',
    '''    diretorio, ctx_default = FONTES[fonte_efeitos]
    base = PIE / diretorio / "preprocessed"
    meta = json.loads((base / "meta.json").read_text(encoding="utf-8"))
    genes_pie = np.asarray(meta["genes"], dtype=str)
    pert_para_id = meta["pert_to_id"]
    nome_ctx = contexto_pie or ctx_default
    if nome_ctx not in meta["context_to_id"]:
        raise SystemExit(
            f"contexto {nome_ctx!r} nao existe em {fonte_efeitos}; disponiveis: "
            f"{sorted(meta['context_to_id'])}"
        )
    ctx_id = meta["context_to_id"][nome_ctx]

    # mmap, nao leitura inteira: o `fold_changes.npy` do orion tem 1,18 GB e o
    # `tested.npy` 295 MB, e so as linhas dos candidatos sao lidas (148 de 8 987
    # no caso que importa). Ler os tres inteiros custaria 2,7 GB para usar 13 MB.
    fc = np.load(base / "fold_changes.npy", mmap_mode="r")
    tested = np.load(base / "tested.npy", mmap_mode="r")
    pert_ids = np.load(base / "pert_ids.npy")
    ctx_ids = np.load(base / "ctx_ids.npy")'''))

PARES.append((
    '''    if not candidatos:
        raise SystemExit(
            f"nenhum nocaute de arc_vcc_25 em {CONTEXTO_PIE} atende o filtro — confira "
            "pert_counts.csv e PIE_arc_vcc_25/preprocessed/meta.json"
        )''',
    '''    if not candidatos:
        raise SystemExit(
            f"nenhum nocaute de {fonte_efeitos} em {nome_ctx} atende o filtro — confira "
            f"pert_counts.csv e {diretorio}/preprocessed/meta.json"
        )'''))

PARES.append((
    '''    return Efeitos(
        alvos=[nome for _p, _n, nome, _l in escolhidos],
        no_painel_vcc=[bool(c[0]) for c in escolhidos],
        fold_changes=out,
        n_genes_cobertos=int(col_nossa.size),
        n_pares_efetivos=efetivos,
    )''',
    '''    return Efeitos(
        alvos=[nome for _p, _n, nome, _l in escolhidos],
        no_painel_vcc=[bool(c[0]) for c in escolhidos],
        fold_changes=out,
        n_genes_cobertos=int(col_nossa.size),
        n_pares_efetivos=efetivos,
        fonte=fonte_efeitos,
        contexto_fonte=nome_ctx,
    )'''))

PARES.append((
    '''    fold_changes: np.ndarray  # (n_alvos, n_genes) float64, 1.0 = sem efeito
    n_genes_cobertos: int
    n_pares_efetivos: int''',
    '''    fold_changes: np.ndarray  # (n_alvos, n_genes) float64, 1.0 = sem efeito
    n_genes_cobertos: int
    n_pares_efetivos: int
    fonte: str
    contexto_fonte: str'''))

PARES.append((
    '''    p.add_argument("--contexto", default="D", choices=["D", "E", "F"])''',
    '''    p.add_argument("--contexto", default="D", choices=["D", "E", "F"])
    p.add_argument("--fonte-efeitos", default="arc_vcc_25", choices=sorted(FONTES),
                   help="de onde vem a razao de dobra medida. `orion` cobre 148 dos 300 "
                        "alvos do painel contra 12 de `arc_vcc_25` [default: %(default)s]")
    p.add_argument("--contexto-pie", default=None,
                   help="contexto dentro da fonte de efeitos; o default de cada fonte esta "
                        "em FONTES")'''))

PARES.append((
    '''    efeitos = _carregar_efeitos_pie(eixo_genes, alvos_vcc, args.n_perts,
                                    min_genes_moventes=args.min_genes_moventes,
                                    limiar_movente=args.limiar_movente,
                                    so_alvos_vcc=args.so_alvos_vcc)''',
    '''    efeitos = _carregar_efeitos_pie(eixo_genes, alvos_vcc, args.n_perts,
                                    min_genes_moventes=args.min_genes_moventes,
                                    limiar_movente=args.limiar_movente,
                                    so_alvos_vcc=args.so_alvos_vcc,
                                    fonte_efeitos=args.fonte_efeitos,
                                    contexto_pie=args.contexto_pie)'''))

PARES.append((
    '''        "fonte_efeitos": str((PIE_ARC / "fold_changes.npy").relative_to(RAIZ)).replace("\\\\", "/"),
        "contexto_pie": CONTEXTO_PIE,''',
    '''        "fonte_efeitos": efeitos.fonte,
        "contexto_fonte_efeitos": efeitos.contexto_fonte,'''))

for antigo, novo in PARES:
    if antigo not in s:
        raise SystemExit("trecho nao encontrado:\\n" + antigo[:200])
    s = s.replace(antigo, novo)

p.write_text(s, encoding="utf-8")
print("ok")
