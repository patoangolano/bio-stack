import pathlib

p = pathlib.Path("C:/Users/Manec/bio-stack/scripts/vcc/avaliador_proxy.py")
s = p.read_text(encoding="utf-8")

PARES = []

PARES.append((
    '''def _carregar_efeitos_pie(eixo_genes: np.ndarray, alvos_vcc: set[str], n_alvos: int,
                          *, min_genes_moventes: int) -> Efeitos:
    """Razoes de dobra medidas por alvo, projetadas no nosso eixo de genes.

    Escolhe os alvos com MAIS genes moventes primeiro. O gate do preset exige que
    a referencia tenha >= 10 genes significativos por perturbacao (`min_gate_size`);
    um alvo com efeito fraco atravessa o avaliador como NaN e sai do agregado, o
    que esvazia a verificacao em vez de testa-la.
    """''',
    '''def _carregar_efeitos_pie(eixo_genes: np.ndarray, alvos_vcc: set[str], n_alvos: int,
                          *, min_genes_moventes: int, limiar_movente: float,
                          so_alvos_vcc: bool) -> Efeitos:
    """Razoes de dobra medidas por alvo, projetadas no nosso eixo de genes.

    Escolhe os alvos com MAIS genes moventes primeiro, e entre eles da preferencia
    aos que estao no painel VCC. O gate do preset exige que a referencia tenha
    >= 10 genes significativos por perturbacao (`min_gate_size`); um alvo com
    efeito fraco atravessa o avaliador como NaN e sai do agregado, o que esvazia a
    verificacao em vez de testa-la.

    ⚠️ O painel VCC limita a 12, nao a 32. Medido em 2026-10-22: `arc_vcc_25` tem
    146 nocautes no contexto `ARC_H1`, e exatamente **12** deles estao entre os
    300 alvos do VCC 2026 — o mesmo 12 que o memorando da rota tirou do
    `train.json` do PIE. (Contando `ARC_H1_VAL`, que e o fold de validacao, sao
    17.) Para um proxy com 32 ou 64 perturbacoes e preciso usar os outros 134
    nocautes, que nao sao alvos do painel. Isso nao prejudica o avaliador: o que
    ele mede e a FORMA do dado — numero de perturbacoes, profundidade, estrutura
    de DE —, nao a identidade dos genes. Com `--so-alvos-vcc` o teto volta a 12.
    """'''))

PARES.append((
    '''    id_para_pert = {v: k for k, v in pert_para_id.items()}
    candidatos: list[tuple[int, str, int]] = []
    for linha in range(fc.shape[0]):
        if int(ctx_ids[linha]) != ctx_id:
            continue
        nome = id_para_pert.get(int(pert_ids[linha]))
        if nome is None or nome not in alvos_vcc:
            continue
        # "Movente" = testado, finito, positivo e a >= 20% de 1. Abaixo disso o
        # arredondamento para inteiro apaga o efeito em celula de baixa contagem.
        v = fc[linha, col_pie]
        t = tested[linha, col_pie]
        movente = t & np.isfinite(v) & (v > 0.0) & (np.abs(np.log2(np.where(v > 0, v, 1.0))) > 0.263)
        candidatos.append((int(movente.sum()), nome, linha))

    if not candidatos:
        raise SystemExit(
            "nenhum dos 190 nocautes de arc_vcc_25 em ARC_H1 esta entre os 300 alvos do "
            "painel VCC 2026 — confira pert_counts.csv e PIE_arc_vcc_25/preprocessed/meta.json"
        )
    candidatos.sort(key=lambda c: -c[0])
    escolhidos = [c for c in candidatos if c[0] >= min_genes_moventes][:n_alvos]
    if len(escolhidos) < n_alvos:
        raise SystemExit(
            f"so {len(escolhidos)} alvos tem >= {min_genes_moventes} genes moventes; "
            f"pedidos {n_alvos}. Baixe --min-genes-moventes ou --n-perts."
        )

    out = np.ones((len(escolhidos), eixo_genes.size), dtype=np.float64)
    efetivos = 0
    for k, (_, _nome, linha) in enumerate(escolhidos):
        v = fc[linha, col_pie].astype(np.float64)
        t = tested[linha, col_pie]
        valido = t & np.isfinite(v) & (v > 0.0)
        out[k, col_nossa[valido]] = v[valido]
        efetivos += int(valido.sum())
    return Efeitos(
        alvos=[nome for _, nome, _ in escolhidos],
        fold_changes=out,
        n_genes_cobertos=int(col_nossa.size),
        n_pares_efetivos=efetivos,
    )''',
    '''    id_para_pert = {v: k for k, v in pert_para_id.items()}
    candidatos: list[tuple[bool, int, str, int]] = []
    for linha in range(fc.shape[0]):
        if int(ctx_ids[linha]) != ctx_id:
            continue
        nome = id_para_pert.get(int(pert_ids[linha]))
        if nome is None:
            continue
        no_painel = nome in alvos_vcc
        if so_alvos_vcc and not no_painel:
            continue
        # "Movente" = testado, finito, positivo e a mais de `limiar_movente` em
        # |log2fc| de 1. O default 0,1 (~7%) nao e generoso: a referencia tem
        # 18 400 celulas contra 400 perturbadas, e nesse desequilibrio o Wilcoxon
        # detecta deslocamento bem abaixo de 20%. O que o limiar protege e contra
        # contar como efeito um gene cujo deslocamento o arredondamento para
        # inteiro apagaria.
        v = fc[linha, col_pie]
        t = tested[linha, col_pie]
        seguro = np.where(v > 0.0, v, 1.0)
        movente = (t & np.isfinite(v) & (v > 0.0)
                   & (np.abs(np.log2(seguro)) > limiar_movente))
        candidatos.append((no_painel, int(movente.sum()), nome, linha))

    if not candidatos:
        raise SystemExit(
            f"nenhum nocaute de arc_vcc_25 em {CONTEXTO_PIE} atende o filtro — confira "
            "pert_counts.csv e PIE_arc_vcc_25/preprocessed/meta.json"
        )
    # Alvo do painel primeiro, depois o efeito mais forte: se o proxy couber nos
    # 12 do painel, usa so eles, e o relatorio diz quantos sao.
    candidatos.sort(key=lambda c: (not c[0], -c[1]))
    escolhidos = [c for c in candidatos if c[1] >= min_genes_moventes][:n_alvos]
    if len(escolhidos) < n_alvos:
        raise SystemExit(
            f"so {len(escolhidos)} nocautes tem >= {min_genes_moventes} genes moventes a "
            f"|log2fc| > {limiar_movente}; pedidos {n_alvos}. Baixe --min-genes-moventes, "
            "baixe --n-perts, ou tire --so-alvos-vcc (que limita a 12)."
        )

    out = np.ones((len(escolhidos), eixo_genes.size), dtype=np.float64)
    efetivos = 0
    for k, (_painel, _n, _nome, linha) in enumerate(escolhidos):
        v = fc[linha, col_pie].astype(np.float64)
        t = tested[linha, col_pie]
        valido = t & np.isfinite(v) & (v > 0.0)
        out[k, col_nossa[valido]] = v[valido]
        efetivos += int(valido.sum())
    return Efeitos(
        alvos=[nome for _p, _n, nome, _l in escolhidos],
        no_painel_vcc=[bool(c[0]) for c in escolhidos],
        fold_changes=out,
        n_genes_cobertos=int(col_nossa.size),
        n_pares_efetivos=efetivos,
    )'''))

PARES.append((
    '''    alvos: list[str]
    fold_changes: np.ndarray  # (n_alvos, n_genes) float64, 1.0 = sem efeito''',
    '''    alvos: list[str]
    no_painel_vcc: list[bool]
    fold_changes: np.ndarray  # (n_alvos, n_genes) float64, 1.0 = sem efeito'''))

PARES.append((
    '''    p.add_argument("--min-genes-moventes", type=int, default=200)''',
    '''    p.add_argument("--min-genes-moventes", type=int, default=200,
                   help="piso de genes com efeito por alvo [default: %(default)s]")
    p.add_argument("--limiar-movente", type=float, default=0.1,
                   help="|log2fc| acima do qual um gene conta como movente [default: %(default)s]")
    p.add_argument("--so-alvos-vcc", action="store_true",
                   help="usa somente nocautes que sao alvos do painel VCC 2026. Teto "
                        "medido: 12 em ARC_H1")'''))

PARES.append((
    '''    efeitos = _carregar_efeitos_pie(eixo_genes, alvos_vcc, args.n_perts,
                                    min_genes_moventes=args.min_genes_moventes)''',
    '''    efeitos = _carregar_efeitos_pie(eixo_genes, alvos_vcc, args.n_perts,
                                    min_genes_moventes=args.min_genes_moventes,
                                    limiar_movente=args.limiar_movente,
                                    so_alvos_vcc=args.so_alvos_vcc)'''))

PARES.append((
    '''        "alvos": efeitos.alvos,''',
    '''        "alvos": efeitos.alvos,
        "n_alvos_no_painel_vcc": sum(efeitos.no_painel_vcc),
        "so_alvos_vcc": args.so_alvos_vcc,
        "limiar_movente": args.limiar_movente,
        "min_genes_moventes": args.min_genes_moventes,'''))

for antigo, novo in PARES:
    if antigo not in s:
        raise SystemExit("trecho nao encontrado:\n" + antigo[:120])
    s = s.replace(antigo, novo)

p.write_text(s, encoding="utf-8")
print("ok")
