# Pesquisa e expansão da bio-stack — outubro de 2026

Levantamento feito em 2026-10-11 por um workflow de agentes: nove domínios pesquisados em paralelo; cada candidato conferido ao vivo por um verificador cético independente (versão, licença, rota de instalação e encaixe nesta máquina); uma síntese; e um crítico de completude e correção, cujas lacunas foram pesquisadas e verificadas na mesma rodada e cujos erros — todos medidos — foram corrigidos na síntese final. A primeira seção diz o que de fato entrou na stack; o relatório que vem depois é o mapa completo, com receita para o que ficou de fora.

## Estado no fim da sessão de 2026-10-11

### Aplicado e conferido

Tudo abaixo está instalado, testado nesta máquina e documentado. Nada foi
commitado: os arquivos novos aparecem como não rastreados no `git status`.

| Camada | O que entrou | Como conferir |
|---|---|---|
| Pipeline de filogenia | `pipelines/filogenia`: Nextflow DSL2 com módulos nf-core, IQ-TREE 3.1.4 (patch registrado no módulo), MAFFT ou FAMSA, ClipKIT ou trimAl, ModelFinder + SH-aLRT + aBayes + UFBoot2 por loco, supermatriz MFP+MERGE, gCF e sCFL, wASTRAL, RO-Crate e BCO pelo nf-prov, e `metodos.md` com versões lidas dos containers e DOIs conferidos. Invólucro `scripts/filogenia/filogenia.py`, figuras com suporte reatribuído por bipartição, agente `phylo-analyst`, comando `/filogenia`, guia [`FILOGENIA.md`](FILOGENIA.md) | `uv run --script tests/check_filogenia.py` (RF 0/18 contra a árvore simulada, 66 s) e `tests/check_filogenia_figura.py` |
| CLIs no WSL | pixi 0.81.0 e 37 ambientes declarados em `bootstrap/pixi-global.toml` (OrthoFinder 3.1.5, subread 2.1.1 e StringTie 3.0.3 inclusos, os três vindos da pesquisa), snapshot com build e sha256 de 45 pacotes, base do TaxonKit (taxdump de 2026-10-12T02:13Z, md5 conferido) | `uv run --script tests/check_wsl_ferramentas.py` |
| Python (Windows, uv) | grupos `evolucao`, `estatistica`, `pesquisa` e `notebooks` (16 pacotes diretos, 57 no fecho), teto `polars<2`; nenhuma versão existente mudou | `uv run python -c "import phykit, pingouin, decoupler, pyalex, marimo"` |
| R (Windows, renv) | `r/`: R 4.6.1 + Bioconductor 3.23 — filogenia comparativa (ape, phangorn, phytools, castor, geiger, corHMM, TreeDist, phylter), figura (ggtree, treeio, ggtreeExtra), DE de referência (DESeq2, edgeR, limma, sva, apeglm, ashr, IHW, qvalue, variancePartition, tximport); 237 pacotes em `r/renv.lock` | `cd r; Rscript -e "renv::status()"` |
| Citação | estilo CSL ABNT fixado por commit em `docs/csl/` | `docs/csl/README.md` |
| Proveniência | RO-Crate (WRROC) e BioCompute Object por execução; trace com `cpu_model`, `attempt` e `workdir`; linha `pipeline-run` em `runs/provenance.jsonl` | `runs/filogenia/<id>/` |
| Revisão adversarial | duas rodadas de agentes: o pipeline, o invólucro, os testes e a integração revistos em cinco lentes, com cada achado reproduzido e conferido por um cético; e 555 afirmações factuais da documentação conferidas contra código, configuração e ferramentas instaladas (41 corrigidas). Os dados do disco também: o `df` do WSL mostra ~1 TB virtual, mas o limite real é o C: (98 GB livres em 2026-10-12) — CLAUDE.md e AGENTS.md dizem isso agora | `tests/check_filogenia_figura.py` trava as regressões achadas |
| Correção de dado | em `PROJETO-ARCHOSAURIA.md`, *Chrysemys picta* passa a `GCF_011386835.1`; o `GCF_007399415.2` que estava lá é *Gopherus evgoodei*. Nenhum script tinha usado o accession errado | `datasets summary genome accession GCF_007399415.2` |

### Depende de decisão sua

Cada item tem receita pronta no relatório abaixo; nenhum foi aplicado porque
muda o repositório público, instala programa no Windows, exige conta externa
ou remove algo que análises antigas citam.

| Item | Por quê | Custo de não fazer |
|---|---|---|
| `LICENSE` + `CITATION.cff` | repositório público sem licença equivale a todos os direitos reservados | Zenodo e WorkflowHub bloqueados; o RO-Crate de cada execução sai inválido sem `manifest.license` |
| Primeiro commit do que entrou | `pipelines/`, `r/`, `bootstrap/pixi-global.*`, `docs/csl/`, `scripts/filogenia/`, os testes novos, o agente e o comando estão fora do git | o campo `commit` do provenance aponta para um código que não contém o pipeline que rodou |
| Varredura de segredos antes desse commit | `runs/provenance.jsonl` guarda milhares de linhas de comando | uma chave que tenha escapado do `redact()` vira pública |
| Nextflow 26.04.7 | correção de dependências vulneráveis | nenhum hoje; o pipeline foi validado no 26.04.6, então revalide depois de subir |
| `sudo apt remove samtools bcftools tabix && sudo apt autoremove` | o htslib 1.19 do apt tem CVEs no decodificador CRAM e responde fora do shell de login; é o `autoremove` que tira o `libhts3t64` | o 1.19 continua respondendo em `bash -c`, cron e systemd. Nada registrado no `provenance.jsonl` depende dele |
| `rm ~/.local/bin/seqkit` | o 2.13.0 solto fica sombreado pelo 2.14.0 do pixi | nenhum em shell de login |
| Python 3.12 | destrava scvi-tools 1.5, pertpy, PyMC 6, anndata 0.13 | nenhum até um projeto pedir; a migração exige `pandas<3` e `anndata<0.13` enquanto `scanpy==1.11.5` |
| Quarto 1.10 e Zotero 10 + Better BibTeX | `metodos.md` → DOCX/PDF com o CSL ABNT; gerenciador de referências com chave estável | o `/metodos` continua saindo só em Markdown |
| Chave da API do OpenAlex | desde 2026-02-13 o OpenAlex exige chave e conta créditos por uso, e o polite pool por e-mail acabou (anúncio de 2026-01-13 na lista openalex-users) | pyalex e o tooluniverse batem no limite sem chave |
| air + jarl no `lint_gate` | o ramo R do hook está inerte (o lintr não existe) | código R entra sem lint |
| Containers por digest | tag de container pode ser republicada | a reprodução depende de a tag não mudar |
| Manifesto sha256 dos dados | dados públicos re-obteníveis | o Métodos não prova qual arquivo entrou |

---

> As marcas *feito* e *falta* do relatório abaixo refletem a máquina no momento da síntese. O estado final da sessão é o da seção anterior: o que ela lista como aplicado está aplicado, mesmo que o relatório ainda diga *falta*.

São 13 domínios, cada candidato com o veredito de um verificador cético. Entram só os confirmados ou corrigidos, já com as correções aplicadas. Versões conferidas ao vivo em 2026-10-11.

## 1. Sumário executivo

- **O núcleo está montado:**
  - 37 ambientes pixi no WSL, com manifesto e snapshot;
  - grupos uv `evolucao`, `estatistica`, `pesquisa` e `notebooks`;
  - renv em `r/` (R 4.6.1, Bioc 3.23, 237 pacotes);
  - `pipelines/filogenia` validado (RF 0/18).
- **O pipeline IQ-TREE 3 reproduz byte a byte.** Usa 3.1.4 por patch, `-T 1`, `--seed`, aBayes como peso do wASTRAL e MAFFT `--threadit 0`. Ainda não está no git (`git ls-files pipelines` = 0).
- **P0 agora:**
  - commit (lista corrigida), LICENSE e CITATION.cff;
  - containers por digest, `NXF_VER` e ponte Nextflow → `provenance.jsonl`;
  - air e jarl, Quarto e Zotero;
  - RepeatMasker no pixi (subread e StringTie já entraram, depois desta síntese);
  - seis downloads pequenos de dados.
- **Regulação de galinha nativa em GRCg7b.** Pelo Ensembl Regulation 116, ENSGALRZ_9RW4D (cromatina aberta) cobre o início de ψ-ENAM, e a janela de ψ-AMBN fica vazia. Nenhuma fonte tem tecido craniofacial.
- **AlphaGenome e GPN-STAR não pontuam ave.** Um só conhece humano e camundongo; o outro tem modelos hg38 com MSA de mamíferos. Ficam como controle em mamífero (P2).
- **Calibração fóssil.** O fossilcalibrations.org morreu; o dump está no Figshare (CC0). Crocodylia vem de Davies 2026. Os limites precisam ser recalculados pelo ICS v2026/06.
- **Manifesto pixi.** Toda instalação pixi nova vai para `bootstrap/pixi-global.toml`: o `pixi global sync` do bootstrap apaga o ambiente que faltar no arquivo (medido).
- **Python 3.11 por ora.** A migração para 3.12 (`pandas<3`, `anndata<0.13`) é P1.

## 2. O que entra agora (P0)

**Já feito e conferido:**

- **pixi.** 37 ambientes com IQ-TREE, ASTER, HyPhy, PAML, MACSE, PRANK, OrthoFinder, subread, StringTie, BLAST+, DIAMOND, miniprot, MMseqs2, HMMER, hts 1.24, seqkit, taxonkit (taxdump baixado), datasets, kent, AGAT, jq, DuckDB, parallel e nf-test.
- **uv.** phykit 2.8.0, decoupler 2.2.0, pingouin 0.7.0, scikit-posthocs 0.17.1, pyhmmer, pygenomeviz, pyalex, habanero, marimo, jupytext e nbstripout; teto `polars<2`; transformers 5.19.0 com `esmc` nativo.
- **R.** ape, phangorn, phytools, castor, geiger, corHMM, TreeDist, phylter, ggtree, ggtreeExtra, DESeq2, edgeR, limma, variancePartition, tximport, sva, apeglm, IHW e qvalue.
- **Pipeline.** Trace com `cpu_model`, nf-prov 1.7.0, `docker -u`, publish `copy`.
- **Documentos.** PROJETO-ARCHOSAURIA §2 corrigido; `docs/csl` com o CSL ABNT; STACK.md já desaconselha o bio-mcp.

| # | Ação | Rota | Versão | Por quê |
|---|---|---|---|---|
| 1 | Commit do pipeline e anexos | git | — | O commit do provenance não contém o código que rodou |
| 2 | LICENSE + CITATION.cff + `manifest.license` | repo | CFF 1.2.0 | Repositório público sem licença; o crate não valida sem ela |
| 3 | Containers por digest | config | — | Só a imagem python está fixada |
| 4 | `nextflowVersion '!>=26.04.6'` + `NXF_VER=26.04.7` | config | 26.04.7 | Patch de segurança |
| 5 | Ponte `pipeline-run` / `env-snapshot` | script | — | O `/metodos` deixa de escrever [VERSAO A CONFIRMAR] |
| 6 | `data/manifesto_raw.tsv` (sha256) | script | — | Dado público re-obtível; DVC desnecessário |
| 7 | air + jarl; `lint_gate` R → `jarl check` | uv tool | 0.12.0 / 0.6.0 | Ramo R inerte (sem lintr) |
| 8 | Quarto + `QUARTO_PYTHON` | winget | 1.10.19 | Métodos e TCC em DOCX/PDF sem LaTeX |
| 9 | Zotero + Better BibTeX (auto-export CSL-JSON) | winget | 10.0.6 / 9.0.71 | Sem gerenciador de referências |
| 10 | subread / StringTie | px | 2.1.1 / 3.0.3 | Transcrição sobre a janela não anotada — **feito** depois da síntese |
| 11 | RepeatMasker + raiz do Dfam 4.0 | px + download | 4.2.4 | TE na janela, com Kimura-CpG |
| 12 | Ensembl Regulation GRCg7b (anotação, atividade, 49 .bb ATAC) | download | 116 | Cis nativo; 42 dos 49 ATAC são de desenvolvimento |
| 13 | Chains galGal6↔GRCg7b + chromAlias | download | 2026-03-11 | Ponte para ChickenGTEx e Pan 2023 |
| 14 | FCD + Davies 2026 + Pan 2023 | download | Figshare | Calibrações; estados ChromHMM |
| 15 | nf-core/rnaseq `-profile test,docker` | Nextflow | 3.27.0 | Medir o teto de memória |
| 16 | Apagar `~/.local/bin/seqkit`; PATH sem duplicata | WSL | — | 2.13.0 sombreado; `~/.local/bin` aparece 3× |
| 17 | ESM C pelo transformers | — | ESMC-600M | Sem instalar; o pacote `esm` exige ≥3.12 e colide com o fair-esm |

No WSL, todo comando passa por script em arquivo, com `MSYS_NO_PATHCONV=1` e `bash -l`:

```bash
git add pipelines r scripts/filogenia tests/check_filogenia.py tests/check_filogenia_figura.py \
  tests/check_wsl_ferramentas.py docs/FILOGENIA.md docs/csl .claude/commands/filogenia.md \
  .claude/agents/phylo-analyst.md .codex/agents/phylo-analyst.toml bootstrap/05_wsl_ferramentas.sh \
  bootstrap/pixi-global.toml bootstrap/pixi-global.snapshot.tsv   # dry-run: .nextflow e renv/library ficam fora
pixi global install -c conda-forge -c bioconda --environment repeatmasker "repeatmasker==4.2.4"
cp ~/.pixi/manifests/pixi-global.toml /mnt/c/Users/Manec/bio-stack/bootstrap/pixi-global.toml
bash /mnt/c/Users/Manec/bio-stack/bootstrap/05_wsl_ferramentas.sh      # sync + snapshot
B=https://ftp.ensembl.org/pub/release-116/regulation/gallus_gallus/GRCg7b
curl --remote-name-all $B/annotation/Gallus_gallus.GRCg7b.regulatory_features.v116.gff3.gz $B/annotation/Gallus_gallus.GRCg7b.regulatory_activity.v116.tsv.gz
curl -L -o fcd.sql.zip https://ndownloader.figshare.com/files/69543690   # md5 651c4ed9570f905f79ca1fae1fa5e185
curl -s -H 'Content-Type:application/json' 'https://rest.ensembl.org/overlap/region/gallus_gallus/Z:55140000-55190000?feature=regulatory'
```

Na lista do `git add` entram também `.agents/skills/source-command-filogenia/`
e o próprio `docs/PESQUISA-STACK-2026-10.md`; os arquivos já rastreados que
mudaram (README, CLAUDE.md, AGENTS.md, docs/STACK.md, pyproject.toml, uv.lock
e outros) entram com `git add -u`. Revise o `git status` antes: há mudanças
anteriores a esta sessão no mesmo conjunto.

Cuidados ao rodar:

- **RepeatMasker.** Expor `famdb.py`, `trf` e os utilitários (`calcDivergenceFromAlign.pl`, `buildRMLibFromEMBL.pl`), conferindo os nomes em `~/.pixi/envs/repeatmasker/bin`.
- **Dfam.** Raiz + curados, ~89 MB, em `~/dfam`, ligados em `share/famdb-3.0.0/Libraries/famdb`. O não curado de um clado vem pela API (`include_raw=true&format=embl`), sem a partição de 4,13 GB.
- **Downloads.** Usar `curl --remote-name-all` com URLs explícitas: chaves `{a,b}` no bash quebram o `-O`.
- **CITATION.cff.** Validar com `uvx --python 3.12 cffconvert --validate`. Sem `.zenodo.json`, que anula o CFF.

## 3. Estratégia de linguagens e ambientes

A divisão é por plataforma e por lock:

- **Windows.** Python (uv), R (renv), Quarto, lint e LSP.
- **WSL.** CLIs bioconda (pixi global, uma por ambiente), Nextflow + Docker e Python sem wheel Windows.
- **Locks que viram Métodos:** `uv.lock`, `r/renv.lock`, `pixi-global.snapshot.tsv` e o tópico `versions` dos containers.

| Tarefa | Ambiente | Por quê |
|---|---|---|
| Pipeline | Nextflow 26.04 DSL2 (WSL), parser estrito, params tipados, `-profile docker`; em nf-core, sempre `-params-file` | Versão vem do container |
| Testes de pipeline | nf-test 0.9.5 + nft-utils 1.2.1 | Snapshot com seed fixa |
| Exploração | Python 3.11 (uv), `.py`, marimo, jupytext | Lock com hash; lint cobre `.py` |
| BAM/VCF em Python | venv no WSL (pysam, cyvcf2, pybedtools manylinux) | Sem wheel Windows; WSL sem gcc |
| Estatística | Python travado; R (DESeq2, edgeR, limma, dream) como referência | Versão citável |
| Filogenia | IQ-TREE 3.1.4: pixi para uso interativo, container no pipeline | Mesma versão nos dois |
| Comparativa, ASR de caráter | R: phytools, castor, corHMM, geiger | CRAN com binário Windows |
| Seleção, ASR de sequência | IQ-TREE `--ancestral`, PAML, HyPhy | Posterior por sítio |
| Figuras de árvore | ggtree / toytree | Instalados |
| scRNA | scanpy CPU; rapids-singlecell no WSL (P1) | Já travado |
| CLIs | pixi global; apt só para a base | apt tem htslib 1.19 com CVEs |
| Tabelas grandes | Parquet + DuckDB 1.5.6 | Tipos preservados |
| Documentos | Quarto 1.10.19 + CSL ABNT + Typst | Sem LaTeX |
| Lint | ruff, basedpyright, air, jarl; `nextflow lint` via WSL | Binários rápidos |
| Redes filogenéticas | Julia + SNaQ, só com reticulação (P2) | Método de referência |

**Ambientes separados:**

- pacote `esm`: ≥3.12, transformers <5;
- gpn: 3.13 e `transformers==5.15.0`;
- Geneformer: `transformers==4.46`;
- clipkit no Windows: `uv tool --python 3.11`;
- roc-validator e cffconvert: `uvx`;
- TOGA2: ≥3.13 no WSL;
- LiftOn, OpenSpliceAI e rapids-singlecell: venv no WSL.

**Não adotar:**

- Snakemake; conda ou mamba; apt para bio; R pelo conda-forge (sem build R 4.6, sem corHMM); DVC.
- Julia sem pergunta; código novo em Go ou Rust; Seurat como via principal; `.ipynb` sem par; `curl | sh`.
- Homônimos e defasados: `aster` e `pixi` do PyPI, `uvx air`, pal2nal 14.1 do bioconda (Perl de 2011, bugs silenciosos), decoupler e pygenomeviz do bioconda, ETE3, MCMCtreeR.

## 4. Por domínio

Rota: px = pixi no WSL (entra no manifesto); uv = Windows; R = renv; ct = container; web = sem instalação; ok = já instalado.

### 4.1 Filogenia

| Ferramenta | O que faz | Por que aqui | Rota | Versão | Prio |
|---|---|---|---|---|---|
| IQ-TREE 3 | ML, UFBoot2/SH-aLRT/aBayes, gCF/sCFL, AU, ASR, IQ2MC | Núcleo do pipeline | px ok; ct por patch | 3.1.4 | P0 |
| ASTER (wASTRAL) | Árvore de espécies | Contraponto coalescente; CASTER-site aceita MSA | px ok; ct | 1.25 (exe Windows v1.25.4.8) | P0 |
| MAFFT, ClipKIT, trimAl | Alinhar e limpar | Etapas do pipeline | px 7.525/2.14.0/1.5.1; ct 7.520/2.11.4/1.5.0 | — | P0 |
| PhyKIT, gotree/goalign, OrthoFinder | Saturação; árvores; ortogrupos | QC e Projeto 1 | uv/px ok | 2.8.0 / 0.5.2 / 3.1.5 | P0 |
| MUSCLE5, FAMSA2, RAxML-NG | Ensembles; proteína em escala; 2ª opinião | Incerteza de alinhamento | px ok | 5.3 / 2.4.1 / 2.0.3 | P1 |
| TreeShrink, IQ2MC, SplitsTree, phyloplace | Outliers; datação; redes; posicionar relíquias | Amostragem maior; aves × crocodilianos; DSPP ~50 pb | px / ok / .msi / ct | 1.4.0 / — / 6.9.7 / 2.3.0 | P1 |

Pendentes (P1):

- **Teste AU, corrigido.** Sem `-s` o 3.1.4 aborta (`File not found`), porque o `.best_model.nex` não traz a sequência:
  `iqtree3 -s 05_concatenado/supermatriz.fa -p 05_concatenado/concatenado.best_model.nex -n 0 -z candidatas.nwk -zb 10000 -au --seed 12345 -T 1 --prefix au`
- **IQ2MC, corrigido.**
  - No `.ctl`: apagar `ckpfile` e `hessianfile`; pôr o flag no `BDparas` (ex. `1 1 0.1 C`), sem o qual o mcmctree 4.10.10 aborta.
  - Sementes distintas por cadeia; renomear `.hessian` para `in.BV`.
  - Rodar `usedata` 0, depois 2, em ≥2 cadeias, e juntar com `mcmctree --combine`.

Armadilhas:

- O modo `-S` não grava UFBoot.
- wASTRAL e gotree leem só o 1º número de `a/b`.
- `-T >1` muda o suporte em 1–2 pontos.
- `--mem` e `--symtest-remove-bad` falham com partições.
- `--ancestral` preenche clado só de gap.
- O ClipKIT mudou padrões na 2.12 e na 2.14; o pipeline fica no 2.11.4 até haver motivo medido.

### 4.2 Evolução molecular e perda gênica

| Ferramenta | O que faz | Por que aqui | Rota | Versão | Prio |
|---|---|---|---|---|---|
| HyPhy (+ hyphy-analyses) | RELAX, BUSTED-E, aBSREL, Contrast-FEL | Relaxamento em Aves+Testudines | px ok; ct 2.5.103 | 2.5.103 | P0 |
| PAML, MACSE, IQ-TREE `--ancestral` | codeml de ramo; códon com frameshift; posterior por sítio | Itens 1–4 do atavismo | px ok | 4.10.10 / 2.07 / 3.1.4 | P0 |
| Datação por ω em ramo misto | Tf = T·(ωp−ωm)/(ωp−ωf) | Item 1; T pelo TimeTree | web + codeml | — | P0 |
| PRANK, CESAR, CAPHEINE | 2º alinhador; projeção de éxon; módulos relax/contrastfel | Regra 4; lesão real; o nf-core não tem esses módulos | px ok / ct / copiar | 251117 / 1.03 / commit | P1 |
| PHAST + PhyloAcc, CAFE5 | CNEE acelerado; famílias | Regra 3; Projeto 1 | ct | 1.9.9 / 2.5.0 / 5.1.0 | P1 |

Subworkflow de seleção (P1), com flags conferidas no `--help` da 2.5.103:

1. MACSE `-seq`/`-seq_lr`; `exportAlignment` com NNN, guardando a cópia sem máscara.
2. IQ-TREE `--ancestral`.
3. HyPhy:
   - RELAX e BUSTED: `--srv Yes --multiple-hits Double+Triple --error-sink Yes`;
   - aBSREL: `--srv Yes --multiple-hits Double+Triple` (sem error-sink);
   - Contrast-FEL: `--srv Yes` com `--branch-set`, `--permutations` e `--q-value` (sem multiple-hits nem error-sink).
4. codeml.
5. PhyKIT.
6. Tabela com nível de evidência e FDR.

Para reprodutibilidade bit a bit: `RANDOM_SEED` e CPU=1.

Armadilhas:

- Stop interno trava HyPhy e codeml, e mascarar apaga a lesão.
- dS satura a 243 Ma.
- Os módulos nf-core `hyphy/*` fixam 2.5.101.
- Proteína RefSeq "LOW QUALITY" vem corrigida: reconfira no genoma.
- EPO e multiz não respeitam códon: realinhe antes de declarar lesão compartilhada.
- Nada de TreeShrink em conjunto com relíquias.

### 4.3 Genômica comparativa

| Ferramenta | O que faz | Por que aqui | Rota | Versão | Prio |
|---|---|---|---|---|---|
| datasets, miniprot, kent, taxonkit, gffread, AGAT | Dados; lesões; MAF/liftOver; anotação | Receita das relíquias | px ok | 18.38.0 / 0.18 / 482 / 0.20.0 / 0.12.9 / 1.7.0 | P0 |
| compleasm / BUSCO 6 | Completude odb12 | Etapa 2 do Projeto 1 | ct (`compleasm run`) | 0.2.10 / 6.1.0 | P1 |
| TOGA2 (+ portal) | Classes de perda | ~950 aves já classificadas | web; WSL Python ≥3.13 | 2.0.10 | P1 |
| make_lastz_chains, pairgenomealign, genomeqc | Cadeias; dotplot; QC | Ave × crocodiliano | ct | 3.1.7 / 3.0.4 / 1.0.1 | P1 |
| LiftOn, eggNOG-mapper, fetchngs | Anotar; função; leituras | rCroNil1/C. siamensis; Projeto 1 | venv WSL / px / ct | 1.0.14 / 2.1.15 / 1.13.0 | P1 |

Montagens:

| Espécie | Montagem | Observação |
|---|---|---|
| Gallus gallus | GCF_016699485.2 | — |
| Alligator mississippiensis | GCF_030867095.1 | — |
| Crocodylus niloticus | GCA_054825085.1 | Embargo VGP de 1 ano |
| Crocodylus siamensis | GCA_052110815.1 | — |
| Gavialis gangeticus | GCA_030020295.1 | Não o "reference" de 2016 |
| Chrysemys picta bellii | GCF_011386835.1 | — |
| Struthio camelus | GCF_040807025.1 | — |

Armadilhas:

- O miniprot chega a 15–22 GB em genoma humano: indexe por cromossomo ou use `-M 2`.
- Ensembl 116 congelado: registre `release=116`.
- 9 scripts ainda chamam `/datasets/v2alpha/`.
- Estar no NCBI Orthologs não prova gene íntegro.
- GENESPACE exige OrthoFinder 2.5.
- OrthoDB devolve 0 em nível inexistente.

### 4.4 CLI no WSL

Já instalados: pixi 0.81.0, jq, DuckDB, parallel, csvtk, hts 1.24, BLAST+, DIAMOND, HMMER, MMseqs2 e minimap2.

P1 (cada um entra no manifesto):

- sra-tools 3.4.1, com `vdb-config` uma vez;
- go-yq 4.54.1;
- hyperfine 2.0.0, ripgrep 15.2.0 e fd 10.5.0;
- Foldseek 10-941cd33 + ProstT5, para MMP20 e KLK4. As EMPs são desordenadas.

Varredura de DSPP/AMTN com as correções medidas. Rodar em `~/`, medindo a RAM com `/usr/bin/time -v`:

```bash
seqkit sliding -g -W 2000000 -s 1900000 GRCg7b.fna -o janelas.fa.gz   # sem -g somem scaffolds < 2 Mb e a ultima janela
seqkit split2 -p 40 -O blocos janelas.fa.gz
parallel --joblog dmd.log --resume-failed --memfree 3G -j 4 diamond blastx -q {} -d ref -F 15 --sensitive --max-hsps 0 --threads 4 --outfmt 6 -o {}.tsv ::: blocos/*
# sem --max-hsps 0: 1 HSP por par e reliquia multi-exon vira fragmento (medido, 2.2.8)
```

Cruzar com `tblastn -seg no -comp_based_stats 0` e com `miniprot -j2 -N 50 --outs=0.5`, usando controle positivo (SULT1E1), e unir no DuckDB.

Armadilhas:

- O `latest_version` do anaconda ordena como texto.
- O DIAMOND não tem tblastn.
- O MMseqs2 GPU funcionou no sm_120 via PTX.
- O Docker não tem nvidia-container-toolkit.

### 4.5 Linguagens

- **P0.** air, jarl e Quarto.
- **P1.**
  - Python 3.12, com `pandas<3` e `anndata<0.13` enquanto `scanpy==1.11.5`, que chama `AnnData.concatenate`.
  - Nextflow LS v26.04.4 por `.lsp.json`.
  - `nbstripout --install`.
- **P2.** piqtree e Julia/SNaQ.
- **Evitar 3.13:** o rpy2-rinterface não tem wheel cp313 para Windows.

### 4.6 Ômica

| Ferramenta | O que faz | Por que aqui | Rota | Versão | Prio |
|---|---|---|---|---|---|
| decoupler, pingouin, scikit-posthocs | Atividade de TF; estatística | Releases que corrigiram resultados errados | uv ok | 2.2.0 / 0.7.0 / 0.17.1 | P0 |
| DE do Bioconductor | DESeq2, edgeR, limma/dream | Referência; efeito aleatório | R ok | Bioc 3.23 | P0 |
| scvi-tools, pertpy, PyMC | Integração; pseudobulk; GLMM | Entre espécies; poucos indivíduos | uv (1.4.2 no 3.11; resto após 3.12) | 1.5.1 / 1.4.0 / 6.3.2 | P1 |
| rapids-singlecell | scanpy em GPU | Atlas acima de 50 mil células | venv WSL cu13 | 0.18.0 | P1 |
| dreamlet, muscat, clusterProfiler, BgeeCall | scRNA multiamostra; enriquecimento; presente/ausente | Galinha; transcrição da relíquia | R (BgeeCall com `download_kallisto = TRUE`) | Bioc 3.23 | P1 |
| differentialabundance | DE em container | Validação cruzada | ct, `-params-file` | 2.0.0 | P1 |

Armadilhas:

- No 3.11 o uv instala versão antiga sem avisar.
- O HCOP para galinha é hipótese.
- Bgee usa GRCg6a.
- Em scRNA, a unidade é a amostra.

### 4.7 Pesquisa e literatura

| Ferramenta | O que faz | Por que aqui | Rota | Versão | Prio |
|---|---|---|---|---|---|
| Zotero + Better BibTeX; Quarto + CSL ABNT | Biblioteca; Métodos e TCC | Sem gerenciador; `docs/csl` pronto | winget | 10.0.6 / 9.0.71; 1.10.19 | P0 |
| OpenAlex (pyalex), Crossref (habanero) | Filtros e contagens; DOI → CSL-JSON; retratação | Núcleo do `/metodos` | uv ok; chave P1 | 0.21 / 2.9.2 | P1 |
| pyzotero MCP, bibverify, LitSense 2.0, ASReview | Biblioteca local; checar .bib; busca por sentença; triagem | Lit-scout reprodutível | uv tool / web | 1.15.2 / 0.5.0 / — / 3.0.8 | P1 |
| SciELO, BDTD; template ABNT (Quarto `abnt`/`abntyp`) | Literatura brasileira; TCC | Senac | web / GitHub | — | P1 |

Armadilhas:

- OpenAlex é pago desde 2026-02-13; a chave vai no header, e o `CHAVES.md` está desatualizado.
- A ABNT do doi.org vem inconsistente: formatar localmente.
- O módulo nf-core `iqtree` cita o IQ-TREE 2; cite doi:10.1093/molbev/msag117.
- O Senac está fora da CAFe.

### 4.8 IA para biologia

| Ferramenta | O que faz | Por que aqui | Rota | Versão | Prio |
|---|---|---|---|---|---|
| ESM C 300M/600M | pLM | Eixo independente de dN/dS | transformers ok | biohub/ESMC-600M | P0 |
| pyhmmer, MMseqs2 | Perfil; busca com GPU | Relíquias | ok | 0.12.3 / 18.8cc5c | P0 |
| Foldseek + ProstT5; Biohub API (ESMFold2, ESMC-6B) | Estrutura sem dobrar; 2º preditor | Globulares; não cabe em 8 GB | px; token | 10-941cd33 | P1 |
| OpenSpliceAI, NT v2, SaProt, Geneformer, bitsandbytes | Splice; linhas de base; célula única; 4 bits | Gallus; comparação | venvs próprios | 0.0.7 / — / 650M / V2 / 0.50.2 | P1 |
| AlphaGenome; GPN-STAR | Regulação; DNA LM com MSA | Só humano/camundongo; hg38 com mamíferos: não pontuam ave | MCP / uv 3.13 | 0.10.0 / 0.9.0 | P2 |

Armadilhas:

- flash-attn não tem sm_120.
- Os pins de transformers são exclusivos entre si.
- ProstT5 dá só 3Di (sem TM-score).
- 4 bits desloca o logit de variante.

### 4.9 Reprodutibilidade

| Ferramenta | O que faz | Por que aqui | Rota | Versão | Prio |
|---|---|---|---|---|---|
| Ponte de proveniência, digests, LICENSE/CFF, manifesto de dados, `NXF_VER` | Métodos e citação | Itens 1–6 do P0 | script/config | — | P0 |
| nf-prov | WRROC + BCO | Feito; BCO falha no schema (#58/#59) | plugin | 1.7.0 | ok |
| lineage, nf-test local, CI | LIDs; snapshot; Actions | Após o commit | config / px ok | 26.04 / 0.9.5 | P1 |
| Zenodo, Software Heritage, WorkflowHub, roc-validator | DOI e SWHID | TCC | web / uvx | 0.12.2 | P1 |
| Higiene do WSL | `appendWindowsPath=false`, `memory=20GB`, zstd, WSLENV | Buscas de comando: 1557 ms → 9 ms | config | WSL 2.7.14 | P1 |
| Apptainer + Slurm, Galaxy, restic/rclone | Cluster; jobs avulsos; backup | Além do laptop | .deb / uvx / bin | 1.5.4 / 1.11.1 / 0.19.1 | P1 |

Armadilhas:

- Fixar a versão de plugin.
- O nf-prov copia params para o crate: nada de segredo em param.
- O `.lineage/` precisa ir para o `.gitignore`.
- `image prune --filter until` olha a data de build: proibido.
- Compactar o VHDX só à mão.
- Shen 2020: 18% das árvores do IQ-TREE não se reproduzem.

### 4.10 Lacuna — de FASTQ a contagem, e transcrição sobre a lesão

- **Matriz de contagem:** nf-core/rnaseq.
- **"A relíquia é transcrita?":** alinhamento ao genoma + contagem em intervalos próprios. Salmon e kallisto dão zero por construção: a janela não é anotada.

| Ferramenta | O que faz | Por que aqui | Rota | Versão | Prio |
|---|---|---|---|---|---|
| nf-core/rnaseq | BAM, `SJ.out.tab`, contagens; RSeQC classifica junções | Fecha fetchngs → DE | ct, `-params-file`, limite 13 GB | 3.27.0 | P0 |
| subread, StringTie | `featureCounts -F SAF -J -G`; montagem na janela | Relíquia contra controles intergênicos | px | 2.1.1 / 3.0.3 | P0 |
| gffcompare, GenMap, igv-reports | Classe `u`; mapeabilidade; figura por locus | Fechar o argumento; afastar parálogo SCPP | px / uv tool WSL | 0.12.10 / 1.3.0 / 1.17.0 | P1 |
| nf-core/scrnaseq | Single-cell | Gallus tem 7626 corridas; Crocodylia, zero | ct (simpleaf, sem CellBender) | 4.2.0 | P1 |

Para o grupo externo de 2–2,6 Gb, o índice STAR não cabe: use `aligner: hisat2`. Nunca `bowtie2_salmon`, que é a rota procariótica.

Armadilhas:

- No Nextflow 26.04, booleano na linha de comando vira String: sempre `-params-file`.
- Labels de 72–200 GB: definir `resourceLimits`.
- Anotação RefSeq: `skip_biotype_qc`, sem `--gencode`.
- O salmon 2.7 muda o TPM baixo.
- Multimapeamento: relatar com e sem filtro.
- No ENA, Crocodylia tem 292 corridas RNA-Seq; falta metadado de tecido.

### 4.11 Lacuna — repetições e TE

| Ferramenta | O que faz | Por que aqui | Rota | Versão | Prio |
|---|---|---|---|---|---|
| UCSC GenArk + API | rmsk com %div, repeatModeler, .2bit | Pronto para rCroNil1 (Dfam 3.9, 41,2%), gavial, frango e jacaré | web; `twoBitToFa` remoto | 2026-10-06 | P0 |
| Dfam 4.0 + RepeatMasker | Famílias (CC0); `.align` e Kimura-CpG | Sem família curada de Archosauria/Aves: usar `-uncurated` ou API por clado | px + download | 4.0 / 4.2.4 | P0 |
| RepeatModeler, ULTRA + longdust | De novo (>32 GB); tandem | C. siamensis; delimitar o DPP | px (`ltr_retriever==2.9.0`) | 2.0.9 / 1.2.1 | P1 |
| TOGA2 × RepeatMasker | BIG_INS/INTRON_GAIN contra família | A inserção é TE? | bedtools | 2.0.10 | P1 |

Receita na janela: API rmsk → `twoBitToFa -noMask` → `RepeatMasker -lib (customLib + DR do clado) -a` → `calcDivergenceFromAlign.pl` → `bedtools intersect` com as lesões.

Armadilhas:

- O .2bit do GenArk é WindowMasker no gavial e no frango.
- A máscara apaga o DPP: busque sem máscara.
- O soft-mask quebra a semente do LASTZ.
- A inserção só conta se for ortóloga e compartilhada.

### 4.12 Lacuna — calibração fóssil

| Recurso | Conteúdo | Por que aqui | Rota | Versão | Prio |
|---|---|---|---|---|---|
| FCD (Figshare, CC0) | 270 calibrações | Site morto; filtrar `PublicationStatus=4`; `MaxAge=0` = sem máximo | download | v1 | P0 |
| Benton 2015; Davies 2026 | Archosauria 247,1–260,2; Neornithes 66–86,8; Crocodylia 74,997–95 | Esqueleto do Projeto 2 | PDF / Figshare | GTS2012 / v1 | P0 |
| PBDB, MCMCtree | Estratigrafia; L/U/B/SN/ST, `usedata=0` | Conferir o espécime; prior efetivo | web / px ok | data1.2 / 4.10.10 | P0 |
| mcmc3r, tbea, CladeDate; B10K (ERDA) | Densidades; c da L; ST; strings de calibração de aves | Prior e Neoaves | R / download | f6a8f78bc5 / 1.8.0 / 1.2 | P1 |

Árvore-esboço (unidade de 100 Ma), já com os limites a recalcular:

```
9 1
(Anolis_carolinensis, (Chrysemys_picta, ((Alligator_mississippiensis, Crocodylus_porosus) 'B(0.74997, 0.95, 1e-300, 0.025)', (Struthio_camelus, ((Gallus_gallus, Meleagris_gallopavo), (Taeniopygia_guttata, Melopsittacus_undulatus) 'L(0.535, 0.1, 0.0374, 0.025)') 'L(0.66, 0.1, 0.5, 1e-300)') 'U(0.868, 0.025)') 'B(2.471, 2.602, 1e-300, 0.025)')) 'B(2.559, 2.959, 1e-300, 0.025)';
```

Ajustes antes de usar:

- Recalcular pelo ICS v2026/06: Archosauria ~247,0–259,94.
- Escolher o mínimo de Neognathae: 66,0, 66,7 ou 68,4 (Vegavis 2025).
- Archelosauria fica sem calibração.

Armadilhas:

- A L padrão põe o 97,5% em 24×tL: ajustar com `tbea::c_truncauchy`.
- Não repetir B em nós aninhados.
- O mínimo duro é `pL=1e-300`.
- Portugalosuchus não serve de mínimo.

### 4.13 Lacuna — atlas regulatório de galinha

| Recurso | Conteúdo | Por que aqui | Rota | Versão | Prio |
|---|---|---|---|---|---|
| Ensembl Regulation GRCg7b | Features + matriz ACTIVE/INACTIVE em 49 epigenomas; 120 .bb (~203 MB) | Nativo; 9RW4D sobre ψ-ENAM; ψ-AMBN vazio | FTP + REST; `.bw` por URL | 116 | P0 |
| Chains + chromAlias | galGal6↔GRCg7b; Z = chrZ = NC_052572.1 | Janela testada ida e volta | download | 2026-03-11 | P0 |
| Pan 2023 | 15 estados ChromHMM, pares enhancer→gene (galGal6) | Fonte única | Figshare (CC BY 4.0) | v3 | P0 |
| ChickenGTEx | 257 bigWigs; molQTL; expressão em 28 tecidos (galGal6, adultos) | Cis e transcrição residual | API / download.cncb.ac.cn | v0 | P1 |
| FAANG API | Metadados; 8 Hi-C | Procedência; contato | web (barra final) | — | P1 |

Armadilhas:

- Sem passar pelo chromAlias, o `bedtools` devolve zero sem erro.
- O GFF3 traz gene_id da release 113.
- ChickenGTEx não declara licença.
- Nenhuma fonte tem mandíbula ou germe dentário: "INACTIVE" só vale para os tecidos amostrados.

## 5. Recursos sem instalação

- **Ensembl REST 116.** EPO de 27 saurópsidos e Regulation GRCg7b. Congelado.
- **UCSC.** API, Table Browser (MAF), 77-way, phyloP, GenArk e chains.
- **TOGA/TOGA2 e coleções de genomas.** VGP/GenomeArk, B10K, DNA Zoo.
- **Ortologia e leituras.** NCBI Datasets v2, OrthoDB v12.2, OMA, ENA Portal.
- **Tempo e fósseis.**
  - TimeTree 5: Gallus–Struthio 107,6 Ma; Gallus–Alligator 243,3 Ma.
  - Open Tree v3, PBDB e FCD/Davies no Figshare.
  - Calibrações do B10K no ERDA.
- **Expressão e regulação.** Bgee 15.2, Expression Atlas/SCEA, GEISHA, ChickenGTEx, Pan 2023, FAANG e g:Profiler e114.
- **Repetições e serviços.** Dfam 4.0 API, Datamonkey (só dado público), Galaxy público, Seqera Containers.
- **Literatura.**
  - LitSense 2.0; OpenAlex (DOI grátis); Crossref/DataCite; SciELO; BDTD.
  - No tooluniverse: PubTator3, Unpaywall, retratação, DataCite, SemanticScholar, ORCID.
- **Publicação.** Zenodo (sandbox), Software Heritage e WorkflowHub.

## 6. Watchlist (P2)

- AlphaGenome e GPN-STAR: controle em mamífero.
- PhyloNetworks/SNaQ: com reticulação.
- BEAST 2.7/X; piqtree.
- GENESPACE (OrthoFinder 2.5), ntSynt, SyRI, Cactus.
- GRASP, AleRax, RERconverge, PseudoChecker.
- SAMap, CAGEE, Seurat, celltypist, squidpy, limpa.
- regtools e mosdepth (`samtools bedcov -j` já cobre); scnanoseq.
- Earl Grey (≥48 GB), TE-Aid.
- bigWig em Python (pyBigWig, deeptools, CrossMap).
- triton-windows 3.6.0.post26; scGPT, UCE, Caduceus, ProtT5/Ankh.
- DataLad, Seqera Platform, r2u; chave Semantic Scholar.

**Descartados:**

- NIM AlphaFold2: 410 desde 2026-09-24.
- Phyloformer: sem licença.
- PubTator3, Unpaywall e MCPs OpenAlex/arXiv/scite: redundantes ou pagos.
- softmask, palaeoverse, ChromHMM, Track Hub Registry, micromamba, bio-mcp.

## 7. Armadilhas transversais

- **Manifesto.** O bootstrap roda `pixi global sync`: ambiente fora de `bootstrap/pixi-global.toml` é apagado. Copie o manifesto de volta depois de cada instalação.
- **anaconda.org.** O `latest_version` ordena como texto: errou mmseqs2, foldseek, prank e r-phytools (há 2.5_2 no conda-forge).
- **WSL.**
  - Script em arquivo, `MSYS_NO_PATHCONV=1` e `bash -l`.
  - O guard bloqueia `curl|sh`, diskpart e escrita em `data/raw`.
- **Nextflow 26.04.** Parser estrito, params tipados, `-params-file`; nada de `def` no config; `env('HOME')`.
- **Montagens.**
  - GRCg7b: relíquias, Ensembl, TOGA2, g:Profiler.
  - galGal6: 77-way, ChickenGTEx, Pan 2023, Bgee.
  - Converta pelo chromAlias e teste o liftOver na janela.
- **Versão citada é a executada.** Tópico `versions` e snapshots. Os pares ClipKIT 2.11.4/2.14.0 e IQ-TREE 3.1.3/3.1.4/3.1.2 divergem.
- **Reprodutibilidade numérica.**
  - IQ-TREE `-T 1`; MAFFT `--threadit 0`; HyPhy CPU=1.
  - O `-t` do wASTRAL muda o arquivo.
  - Sementes distintas por cadeia no MCMCtree.
- **Ausência de evidência não é evidência de ausência.**
  - Zeros silenciosos: BLAST web, OrthoDB, `bedtools` com nome errado, Salmon fora da anotação.
  - Use dois ou três métodos com controle positivo.
- **Memória.** O WSL tem 15 GB e não há `.wslconfig`. Medir antes de paralelizar:
  - miniprot: 15–22 GB;
  - STAR: ≥32 GB;
  - GenMap: 10n;
  - LAST MAM8: ~50 GB;
  - RepeatModeler: >32 GB.
- **Compilação.** O WSL não tem gcc e o Windows não tem Rtools: prefira binários.
- **Python 3.11.** Pacotes que exigem ≥3.12 são resolvidos para versão antiga sem aviso.
- **Licenças.**
  - AGPL (ASTER, RAxML-NG): só como serviço.
  - CC-BY-NC-SA: NT v2, Ankh.
  - GPL-3: OpenSpliceAI.
  - ChickenGTEx: sem licença.
  - VGP: embargo de 1 ano.
  - `.over.chain`: só uso não comercial.
- **Serviços congelados ou pagos.** Ensembl 116, OpenAlex pago, FCD morto, OSF sem projeto novo após 2026-11-16.

## 8. Fontes

- **Filogenia e evolução:**
  - https://github.com/iqtree/iqtree3/releases/tag/v3.1.4
  - https://doi.org/10.1093/molbev/msag117
  - https://github.com/chaoszhang/ASTER
  - https://github.com/veg/hyphy-analyses
  - https://github.com/veg/CAPHEINE
  - https://github.com/abacus-gene/paml/releases/tag/v4.10.9
  - https://www.agap-ge2pop.org/macsee-pipelines/
- **Genômica e CLI:**
  - https://github.com/lh3/miniprot
  - https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/multiz77way/
  - https://genome.senckenberg.de/download/TOGA2/
  - https://pixi.prefix.dev/latest/global_tools/manifest/
  - https://github.com/bbuchfink/diamond/releases/tag/v2.2.8
- **Ambientes e reprodutibilidade:**
  - https://docs.seqera.io/nextflow/migrations/26-04
  - https://nf-co.re/blog/2026/parameter-types
  - https://github.com/quarto-dev/quarto-cli/releases/tag/v1.10.19
  - https://github.com/etiennebacher/jarl
  - https://github.com/nextflow-io/nf-prov
  - https://doi.org/10.1038/s41467-020-20005-6
- **Ômica, IA e pesquisa:**
  - https://github.com/scverse/decoupler
  - https://huggingface.co/biohub/ESMC-600M
  - https://pypi.org/project/alphagenome/
  - https://help.openalex.org/api/authentication/
  - https://retorque.re/zotero-better-bibtex/
- **Lacunas:**
  - https://nf-co.re/rnaseq/3.27.0
  - https://www.dfam.org/releases/Dfam_4.0/families/FamDB/
  - https://doi.org/10.6084/m9.figshare.34046970
  - https://doi.org/10.26879/424
  - https://doi.org/10.1098/rspb.2026.0635
  - https://stratigraphy.org/news/156
  - https://ftp.ensembl.org/pub/release-116/regulation/gallus_gallus/GRCg7b/annotation/
  - https://doi.org/10.1126/sciadv.ade1204
  - https://doi.org/10.1093/nar/gkaf731

## Correções que a crítica mediu nas receitas

Um crítico de completude e correção leu a primeira síntese e testou as receitas nesta máquina. Os erros abaixo foram medidos por ele e entraram corrigidos no relatório acima; ficam listados para quem reaproveitar uma receita antiga.

- Comando do teste AU (P1) não roda: `iqtree3 -p concat.best_model.nex -n 0 -z candidatas.nwk -zb 10000 -au --seed 12345 -T 1 --prefix au`. Medido no WSL com o `concatenado.best_model.nex` do run 20261011-212604: o IQ-TREE 3.1.4 lê os intervalos do charset como nomes de arquivo e aborta — `Reading alignment file 1-938 939-2203 2204-2954 2955-4055 ... ERROR: File not found`. O `.best_model.nex` só carrega modelo e partição, não a sequência; com `-s supermatriz.fa` acrescentado, o mesmo comando completa e escreve a seção USER TREES com p-AU. Corrigir o comando para incluir `-s <supermatriz>`.
- Receita do Projeto 2, passo 5: `seqkit sliding -W 2000000 -s 1900000 GRCg7b.fna -o janelas.fa.gz` sem `-g/--greedy` perde parte do genoma. Medido com seqkit 2.14.0: sem `-g`, a última janela menor que W é descartada e uma sequência mais curta que W não gera janela nenhuma (teste: registro de 25 pb com -W 10 -s 9 deu 2 janelas sem `-g` e 3 com; registro de 6 pb deu zero sem `-g`). Numa varredura irrestrita de DSPP/AMTN isso apaga todo scaffold menor que 2 Mb — justamente os não-placed, onde relíquia costuma cair — e o trecho terminal de cada cromossomo. Acrescentar `-g`.
- Mesma receita, passo 5: o DIAMOND blastx como está reporta um HSP só por par consulta-alvo, então relíquia multi-éxon aparece como fragmento único. Medido com diamond 2.2.8 sobre uma janela sintética com 3 éxons da mesma isca separados por 1,5 kb: `-F 15 --sensitive --outfmt 6` devolveu apenas o primeiro HSP (501-680 / 1-60), e `--long-reads` também; só com `--max-hsps 0` saíram os três (501-680, 2181-2360, 3861-4040). Sem `--max-hsps 0` a contagem de éxons e as coordenadas do bloco ficam erradas, e o cruzamento com miniprot/tBLASTn perde o alinhamento.
- AlphaGenome não cobre galinha, logo o item P1 `uv add alphagenome==0.10.0 — predições regulatórias em lote` não atende 'o lado perdido por regulação do /atavismo'. Medido no wheel 0.10.0: `alphagenome/models/v1/dna_model.py` define `class Organism(enum.Enum)` apenas com HOMO_SAPIENS e MUS_MUSCULUS, e o README trata só de genoma humano (e GENCODE humano/camundongo). Serve como controle em mamífero, não em GRCg7b, Crocodylus, kakapo, tartaruga ou Anolis.
- Item P1 do subworkflow de seleção aplica o mesmo conjunto de flags aos quatro métodos do HyPhy — `--srv Yes --multiple-hits Double+Triple --error-sink Yes` para RELAX/BUSTED-E/aBSREL/Contrast-FEL. Medido com hyphy 2.5.103 (`hyphy <método> --help`): RELAX e BUSTED expõem srv, multiple-hits e error-sink; aBSREL expõe srv e multiple-hits, mas não error-sink; Contrast-FEL expõe srv e não expõe multiple-hits nem error-sink (tem branch-set, permutations, p-value, q-value). O texto da armadilha em §4.2 está correto; é o item do plano que generaliza demais.
- As receitas P1 que instalam com `pixi global install` (treeshrink, foldseek, sra-tools, mosdepth, go-yq, hyperfine/ripgrep/fd, eggnog-mapper) são desfeitas pelo próprio P0 do plano se o ambiente não for escrito de volta em `bootstrap/pixi-global.toml`. Medido em PIXI_HOME isolado: instalei dois ambientes, removi um do manifesto e rodei `pixi global sync` — a saída foi `(removed) t_tree` e o ambiente e o binário exposto desapareceram. Como `bootstrap/05_wsl_ferramentas.sh` copia o manifesto do repositório para `~/.pixi/manifests/` antes do sync, qualquer ambiente instalado à mão e ausente do arquivo versionado é apagado na próxima execução do bootstrap (que o plano manda rodar e que `tests/check_wsl_ferramentas.py` exercita). Cada item P1 precisa do passo 'acrescentar [envs.<nome>] ao bootstrap/pixi-global.toml', como o item do OrthoFinder já faz.
- A lista de `git add` do P0 não cumpre o objetivo declarado ('o commit que o provenance grava não contém o código que rodou'): ela nomeia pipelines, bootstrap/05_wsl_ferramentas.sh, bootstrap/pixi-global.toml, o snapshot, r, docs/csl e tests/check_wsl_ferramentas.py, e deixa de fora arquivos não rastreados que o próprio plano cita como parte do pipeline — `scripts/filogenia/filogenia.py` (o lançador que gera o script com NXF_VER e que o P0 do nextflow.config manda editar), `tests/check_filogenia.py` (o teste de RF=0, também citado nesse mesmo P0), `docs/FILOGENIA.md` (a decisão sobre ClipKIT/trimAl que a tabela §2.5 referencia), `.claude/commands/filogenia.md` e `.claude/agents/phylo-analyst.md`. Verificado com `git status`/`git ls-files`: os cinco estão untracked e `git ls-files pipelines` retorna 0.

---

## Anexo — catálogo completo levantado e verificado

Cada linha é um candidato levantado por um pesquisador e conferido por um verificador cético independente, com fonte consultada ao vivo em 2026-10-11. *Prioridade* é a final (depois da verificação): P0 instalar agora, P1 receita pronta para quando o caso aparecer, P2 observar. Refutados ficam listados para ninguém propor de novo. As seções marcadas *lacuna* foram pesquisadas depois, a partir do que o crítico de completude apontou como faltando.


### IA para biologia (local e hospedada)

| Candidato | O que é | Rota | Versão | Prioridade | Veredito |
|---|---|---|---|---|---|
| ESM C 300M / 600M (pacote esm, Biohub) | Modelo de linguagem de proteína de ponta da EvolutionaryScale, agora publicado pelo Chan Zuckerberg Biohub. As variantes 300M e 600M são so… | Rota barata, sem instalar nada: os pesos carregam com o transformers 5.19.0 JA instalado + torch 2.11.0+cu128, no venv … | esm 3.4.1.post1 (PyPI 2026-09-16; GitHu… | P0 | corrigido |
| Foldseek | Busca e clusterização de estruturas de proteína por alfabeto estrutural 3Di. Inclui Foldseek-Multimer para complexos, aceleração por GPU no… | WSL: `conda/pixi`-free nao impede — o bioconda JA tem foldseek 10.941cd33 linux-64, mesma versao do release. Use o tarb… | Release 10-941cd33 (2025-01-19); master… | P0 | corrigido |
| ProstT5 dentro do Foldseek (createdb --prostt5-model) | Modelo de tradução sequência→3Di. Permite montar banco estrutural a partir de FASTA e buscar estruturalmente sem nenhuma estrutura, por vol… | `foldseek databases ProstT5 weights tmp` e depois `foldseek createdb db.fasta db --prostt5-model weights [--gpu 1]` — o… | Comandos verificados no README do maste… | P0 | corrigido |
| MMseqs2 (com build de GPU) | Busca de sequência e perfil em escala, clusterização, taxonomia e workflow de download de bancos (UniRef, NR, Swiss-Prot, Pfam). O build de… | WSL: bioconda ja tem 18.8cc5c linux-64 (paridade com o release). Para GPU use mmseqs-linux-gpu.tar.gz (148 MB, glibc >=… | Release 18-8cc5c (2025-07-27); master c… | P0 | corrigido |
| pyhmmer | Ligações Cython para o motor HMMER, com busca de perfil (hmmsearch, hmmscan, phmmer, nhmmer) executada dentro do processo Python e paraleli… | uv pip install pyhmmer — há wheel win_amd64 publicado no release atual, além de manylinux e macOS. | 0.12.3 (2026-08-19), igual ao declarado | P0 | confirmado |
| Biohub Platform API (ESMFold2, ESMC-6B, SAE hospedados) | API hospedada que serve os modelos ESM sem setup local. Modelos nomeados no README: esmc-600m-2024-12, esmc-6b-2024-12 e esmfold2-fast-2026… | Token no developer console de biohub.ai (cliente dentro do pacote esm, que exige ambiente proprio com Python >=3.12). A… | Modelos esmfold2-fast-2026-05, esmc-6b-… | P1 | corrigido |
| NIM AlphaFold2 e AlphaFold2-Multimer | Os dois únicos cartões da coleção de biologia do build.nvidia.com que o servidor nvidia-bio desta stack não cobre. O NIM expõe três rotas: … | Acrescentar as chaves ao catálogo em servers/nvidia_bio_mcp/server.py. O slug hospedado não foi verificado: descubra co… | Cartões alphafold2 e alphafold2-multime… | descartar | **refutado** |
| AlphaGenome (cliente Python oficial) | Cliente da API do AlphaGenome para predição de tracks regulatórias (expressão, acessibilidade, splicing, contato) a partir de sequência, em… | uv pip install alphagenome (wheel py3-none-any, requires-python >=3.10) + chave de API. | 0.10.0 (2026-10-08), igual ao declarado | P1 | confirmado |
| GPN-STAR e GPN-MSA (songlab) | Modelos de linguagem genômica treinados sobre alinhamentos multiespécie. O repositório publica a família gpn-star-hg38 em 25M, 85M e 200M p… | uv pip install gpn (wheel py3-none-any) + pesos do HF. Exige ambiente uv separado: requires-python >=3.13,<3.14 e pin t… | gpn 0.9.0 (PyPI e release v0.9.0 em 202… | P1 | confirmado |
| OpenSpliceAI | Reimplementação em PyTorch do SpliceAI, com treino, fine-tuning para outras espécies, predição e análise de variante. O SpliceAI original (… | No WSL: uv pip install openspliceai. Depende de pysam e mappy, que não têm wheel para Windows; os classificadores do pa… | 0.0.7 (2026-06-24); repo com push em 20… | P1 | confirmado |
| Nucleotide Transformer v2 (InstaDeep) | Família de DNA LM em 50M a 500M parâmetros; a variante multi-species de 500M tem 498 346 364 parâmetros em F32 no card. | transformers + pesos do HF (InstaDeepAI/nucleotide-transformer-v2-500m-multi-species). | Card de 500m modificado em 2025-10-22; … | P1 | confirmado |
| SaProt 650M | Modelo de linguagem de proteína com vocabulário estrutura-ciente: cada token combina aminoácido e token 3Di do Foldseek. Variantes AF2 e PD… | transformers + pesos westlake-repl/SaProt_650M_AF2; o 3Di vem do Foldseek (ProstT5 ou estrutura real). | Pesos westlake-repl/SaProt_650M_AF2 atu… | P1 | confirmado |
| Geneformer V2 | Modelo de fundação de transcriptoma de célula única treinado por ranking de expressão. O repositório HF publica Geneformer-V1-10M e Genefor… | pip install do repo HF (git) em ambiente uv SEPARADO: requirements.txt pinna transformers==4.46, incompativel com o tra… | Repo HF ctheodoris/Geneformer modificad… | P1 | corrigido |
| scvi-tools | Modelos variacionais probabilísticos para célula única: scVI, scANVI, totalVI, destVI, integração com correção de lote e quantificação de i… | `uv pip install scvi-tools` no ambiente atual NAO entrega a 1.5.1: o venv e Python 3.11 (.venv/pyvenv.cfg version_info … | 1.5.1 (2026-09-10) exige Python >=3.12;… | P1 | corrigido |
| CellTypist | Anotação automática de tipo celular por regressão logística com modelos pré-treinados, mais rotulagem hierárquica e por maioria sobre clust… | uv pip install celltypist (wheel py3-none-any; dependências leves: scanpy, scikit-learn, leidenalg — todas já presentes… | 1.7.1 (2025-06-25), igual ao declarado | P1 | confirmado |
| triton-windows + bitsandbytes | Dois habilitadores de kernel no Windows: triton-windows traz wheels win_amd64 do Triton (o Triton oficial publica só manylinux), e bitsandb… | bitsandbytes: `uv pip install bitsandbytes` no grupo gpu (wheel win_amd64, MIT, torch<3,>=2.4) — pronto para uso. trito… | triton-windows 3.8.0.post29 (2026-09-28… | P1 | corrigido |
| scGPT e UCE (modelos de fundação de célula única) | scGPT é transformer generativo para célula única (MIT); UCE é o Universal Cell Embedding de Stanford, que embeda célula sem precisar de rót… | scGPT: nao existe checkpoint oficial no HF — os pesos sao pastas do Google Drive linkadas no README (whole-human, conti… | scgpt 0.2.4 (2025-03-31), repo com push… | P2 | corrigido |
| Caduceus (via mamba-ssm) | DNA LM bidirecional equivariante a complemento reverso, baseado em espaço de estados (Mamba), com janela de 131 kb. O checkpoint ps_seqlen-… | Nao e preciso compilar. O release v2.3.2.post1 no GitHub traz 80 wheels pre-construidos linux (cu11/cu12/cu13 x torch 2… | Card do Caduceus modificado em 2025-10-… | P2 | corrigido |
| ProtT5-XL (encoder, half) e Ankh | Geração anterior de pLM ainda muito usada: Rostlab/prot_t5_xl_half_uniref50-enc (619 577 downloads em 30 dias) e a família Ankh (base e lar… | transformers + pesos do HF. Ankh exige atenção à licença. | prot_t5_xl_half_uniref50-enc modificado… | P2 | corrigido |
| Phyloformer | Estimativa de filogenia por rede de autoatenção: prediz distâncias evolutivas a partir do alinhamento, sem busca em espaço de árvores. | Só pelo GitHub (clone). Não instalar agora. | Repositório lucanest/Phyloformer com pu… | descartar | confirmado |
| BixBench, LAB-Bench, Biomni e HeurekaBench | Conjunto de referência para avaliar agentes de ciência: BixBench (mais de 50 cenários reais de análise e ~300 perguntas abertas), LAB-Bench… | Só leitura por ora. O pacote biomni existe no PyPI (0.0.8, 2025-10-27, wheel py3-none-any) se quiser inspecionar o ambi… | BixBench arXiv 2503.00096 v1 2025-02-28… | P2 | corrigido |

Armadilhas registradas:

- flash-attn está fora desta máquina, e isso é estrutural. O PyPI 2.8.3.post1 (2026-06-11) publica apenas sdist, zero wheels, com classificador 'Operating System :: Unix'. O README lista como suportadas Ampere, Ada e Hopper — Blackwell de consumo (sm_120) não aparece; o FlashAttention-4 é escrito em CuTeDSL para Hopper e Blackwell de datacenter (B200), não para RTX 5070. E sobre Windows o próprio R…
- Três pins de transformers mutuamente exclusivos. O pacote esm exige transformers<5.0.0,>=4.57.6; o gpn pinna transformers==5.15.0; o PyPI já está em transformers 5.19.0 (2026-10-06). Instalar os dois no mesmo ambiente resolve para algo que quebra um dos dois em silêncio — o sintoma é erro de atributo em tempo de carga de modelo, não erro de instalação. São três ambientes uv separados: um para esm…
- O gpn 0.9.0 exige requires-python >=3.13,<3.14 — faixa de uma única versão menor. Não compartilha interpretador com um ambiente de 3.12 nem de 3.14, e o Python do sistema nesta máquina é 3.14. Precisa de `uv python install 3.13` e ambiente dedicado.
- scgpt 0.2.4 não resolve contra esta stack: pinna scvi-tools<1.0 (o atual é 1.5.1), datasets<3.0.0, orbax<0.1.8 e depende de torchtext, que foi descontinuado e não tem build para torch 2.11. Tentar instalar no ambiente principal rebaixa scanpy/scvi-tools e pode arrastar torch.
- bioconda atrasa em relação ao upstream nas ferramentas do Steinegger: foldseek está em 8.ef4e960 no bioconda contra release 10-941cd33 no GitHub (2025-01-19) e master com push em 2026-09-09; mmseqs2 está em 17.b804f no bioconda contra release 18-8cc5c. E o bioconda não tem build de GPU nem build win-64 do foldseek — para GPU é o tarball oficial.
- GPU do Foldseek e do MMseqs2 só acelera busca de sequência e de perfil, exige driver >=525.60.13 (tem 617.14, ok) e glibc >=2.29 para o mmseqs / >=2.17 para o foldseek (Ubuntu 24.04 tem 2.39, ok). Os READMEs exigem Turing como piso e Ampere ou mais novo para velocidade plena, e avisam que os binários pré-compilados não rodam em Volta/Pascal. sm_120 não é nomeado em nenhum dos dois: faça teste de …
- Banco estrutural construído com ProstT5 contém só sequência 3Di, sem coordenada de Cα. Isso desabilita --alignment-type 1, TM-score e saída de LDDT, e restringe a busca a monômero. Reportar TM-score a partir de um banco ProstT5 é erro de método, não de ferramenta — o README avisa explicitamente.
- Os bancos do Foldseek escalam feio: Alphafold/UniProt e Alphafold/UniProt50 cobrem o universo inteiro do AFDB (tamanho não verificado nesta sessão; a tentativa de resolver content-length nas URLs de download falhou). Comece por Alphafold/Swiss-Prot ou PDB. E baixe dentro de ~/ no WSL: /mnt/c atravessa a fronteira de sistema de arquivos e é ordens de grandeza mais lento, exatamente no padrão de mu…
- Licenças não comerciais em dois modelos recomendados: Nucleotide Transformer v2 e a família Ankh são CC-BY-NC-SA-4.0. Uso acadêmico cabe, mas a cláusula ShareAlike contamina derivados e a restrição precisa constar na seção de Métodos gerada por /metodos. OpenSpliceAI é GPL-3.0-only: não copie trechos do código para dentro de scripts/ do repositório sem decidir a licença do repositório primeiro.
- openspliceai, tiledbsoma e as dependências pysam/mappy não têm wheel win_amd64 — os classificadores do openspliceai listam só MacOS e Linux. Isso define a fronteira: tudo que toca pysam é WSL. O cellxgene-census é puro Python, mas depende de tiledbsoma, então vai para o WSL junto.
- Os classificadores do PyPI mentem por desatualização nos dois sentidos: pyhmmer 0.12.3 declara 'Operating System :: POSIX' e no entanto publica wheel win_amd64 no release atual. Decida pela lista de arquivos do release, não pelo classificador.
- O caminho hospedado do NIM não é o caminho documentado do NIM local. A documentação do AlphaFold2 NIM descreve protein-structure/alphafold2/predict-structure-from-sequence, predict-msa-from-sequence e predict-structure-from-msa; o servidor nvidia_bio_mcp usa a forma biology/<vendor>/<modelo>/<ação>. Descubra o slug com nim_probe_endpoints antes de acrescentar o modelo ao catálogo, como já foi fei…
- Windows não tem wheel do Triton oficial — só manylinux. Qualquer uso de torch.compile ou de kernel Triton do lado Windows depende do fork triton-windows (3.8.0.post29). O pareamento de versão entre Triton 3.8 e torch 2.11 não foi verificado nesta sessão; valide antes de contar com compilação em produção.
- ESMC-6B em fp16 pede ~12 GB e não cabe nos 8 GB — e a tentação é carregar em 4-bit e seguir. Para pontuação de variante por verossimilhança, quantização agressiva desloca os logits, que é exatamente a grandeza que você está medindo. Use a API hospedada do Biohub para 6B e deixe o 4-bit para tarefas de embedding, onde o erro é absorvido pelo classificador downstream.

### Reprodutibilidade, dados e execução

| Candidato | O que é | Rota | Versão | Prioridade | Veredito |
|---|---|---|---|---|---|
| pixi (workspace linux-64 com pixi.lock para as CLIs do WSL) | Gerenciador de pacotes conda escrito em Rust. Usa manifesto pixi.toml e lockfile pixi.lock multiplataforma, resolve PyPI com o uv embutido … | No WSL, sempre por script em arquivo com 'bash -l': 1. 'curl -fsSL -o ~/pixi-install.sh https://pixi.prefix.dev/install… | 0.81.0 (GitHub v0.81.0 de 2026-09-15; c… | P0 | corrigido |
| Ponte Nextflow/WSL -> runs/provenance.jsonl (linhas pipeline-run e en… | Dois scripts pequenos, sem dependencia nova: 1. nf_para_jsonl.py (Python do Windows): depois de cada run, le pipeline_info/ (trace e versoe… | 1. Commitar pipelines/filogenia antes dos runs. Se nao der, gravar na linha pipeline-run o sha256 do conteudo do pipeli… | n/a (schema 3 do provenance.jsonl; topi… | P0 | corrigido |
| nf-prov 1.7.0 (Workflow Run RO-Crate e BioCompute Object por run) | Plugin da nextflow-io que grava proveniencia padronizada no fim do run. Formatos: - wrroc: RO-Crate 1.1 com os perfis Process/Workflow/Prov… | Nao precisa instalar; e configuracao no nextflow.config: - plugins { id 'nf-prov@1.7.0' } - prov { enabled = true; form… | 1.7.0 (2026-01-05; exige Nextflow >=25.… | P0 | confirmado |
| Relatorios nativos do Nextflow com campos de reprodutibilidade (trace… | Relatorios que o Nextflow ja gera: - trace TSV por tarefa: 42 campos disponiveis, incluindo hash, container, cpus, realtime, peak_rss, cpu_… | Mexer so no que falta em pipelines/filogenia/nextflow.config: - params.trace_report_suffix = new java.util.Date().forma… | Nextflow 26.04.6 instalado; 26.04.7 lan… | P0 | corrigido |
| nextflow lineage (LIDs por run, tarefa e arquivo) | Metadados de linhagem nativos do Nextflow. Entraram na 25.04 e ainda estao marcados 'Experimental: may change in a future release'. Cada ru… | No nextflow.config: lineage { enabled = true; store.location = "${env('HOME')}/nf-lineage/filogenia" }. Acrescentar .li… | Nextflow 26.04.6 (recurso existe desde … | P1 | corrigido |
| nf-test 0.9.5 | Framework de testes para processos, subworkflows e pipelines Nextflow, com snapshots. - 0.9.3: --stop-on-first-failure e --smart-testing; -… | No WSL, por script em arquivo: 1. 'curl -fsSL -o ~/get-nftest.sh https://get.nf-test.com'; LER o arquivo. 2. 'bash ~/ge… | 0.9.5 (2026-03-19). Precisa ser >=0.9.4… | P1 | corrigido |
| Containers por digest: biocontainers para ferramenta unica, Seqera Co… | Referenciar cada imagem por @sha256: em vez de tag. - biocontainers (quay.io/biocontainers): uma imagem por pacote bioconda; serve para fer… | Nada a instalar. Fluxo: 1. 'nextflow inspect main.nf -profile docker -format config > conf/containers.config'; 2. troca… | docker-ce 29.9.0 + buildx v0.38.0 (WSL)… | P0 | confirmado |
| Higiene do WSL para Nextflow e Docker (work dir, limpeza, VHDX, backu… | Seis regras: - work dir e lineage no ext4 (~/); - publicacao em modo copy para /mnt/c; - container rodando com o UID do usuario; - limpeza … | 1. Ja feito: work dir ~/nf-work/filogenia (comando /filogenia), publishDir mode 'copy' e docker.runOptions '-u $(id -u)… | WSL 2.7.14.0 instalado (kernel 6.18.33.… | P1 | corrigido |
| LICENSE + CITATION.cff (validado com cffconvert) | Arquivo de licenca e metadados de citacao no formato CFF 1.2.0, na raiz do repo. O GitHub mostra 'Cite this repository'. O Zenodo le um sub… | 1. LICENSE na raiz: MIT ou Apache-2.0 para o codigo. Declarar no README que docs e dados derivados ficam sob CC-BY-4.0.… | CFF 1.2.0. cffconvert 2.0.0 e de 2021-0… | P0 | corrigido |
| Manifesto de dados brutos com sha256 (opcao 'nada' para versionar dad… | Em vez de DVC ou DataLad, versionar so data/raw/MANIFESTO.tsv e o script de download. O manifesto guarda accession, URL, release da fonte, … | Um script, por exemplo scripts/dados/manifesto.py (hashlib.sha256 em streaming), gera e confere data/manifesto_raw.tsv,… | n/a (ncbi-datasets-cli 18.38.0 no conda… | P0 | corrigido |
| CI no GitHub Actions (nf-test + pytest + pre-commit) | Workflows em .github/workflows: - job em ubuntu-latest com nf-core/setup-nextflow e nf-core/setup-nf-test, rodando 'nf-test test --profile … | Criar .github/workflows/ci.yml com setup-nextflow v3.0.1, setup-nf-test v2.1.0, setup-uv v10.3.0 e setup-pixi v0.11.0. … | setup-nextflow v3.0.1; setup-nf-test v2… | P1 | confirmado |
| Zenodo (integracao GitHub) + Software Heritage: DOI e SWHID | Zenodo: cada release do GitHub vira um registro com DOI (um DOI de conceito e um DOI por versao), usando metadados do .zenodo.json ou do CI… | - Zenodo: entrar no zenodo.org com o GitHub, ligar o repo e publicar uma release (testar antes no sandbox.zenodo.org). … | servico (consultado em 2026-10-11) | P1 | corrigido |
| WorkflowHub (registro do pipelines/filogenia) | Registro FAIR de workflows, independente do motor (Nextflow e um dos mais usados). Usa Workflow RO-Crate e cunha DOI DataCite por versao co… | 1. Criar conta com GitHub ou LS Login. 2. Criar ou entrar num Team dentro de um Space. 3. Importar o workflow do git. 4… | servico (consultado em 2026-10-11) | P1 | confirmado |
| ro-crate-py (rocrate) + roc-validator | rocrate e uma biblioteca e CLI para criar e ler RO-Crates: RO-Crate 1.3 por padrao, 1.0 a 1.2 suportados; comandos 'rocrate init', 'rocrate… | No Windows, sob demanda e isolado, sem tocar no uv.lock: - 'uvx --from roc-validator rocrate-validator validate <dir>';… | rocrate 0.16.0; roc-validator 0.12.2 | P1 | corrigido |
| Apptainer + perfil Slurm proprio (cluster da universidade ou SDumont/… | Runtime de containers para HPC sem daemon, sucessor do Singularity. O Nextflow tem o escopo apptainer: enabled, autoMounts (padrao true), c… | No WSL: - instalar o .deb oficial: 'sudo apt install ./apptainer_1.5.4_amd64.deb' (release do GitHub); ou apptainer 1.5… | Apptainer 1.5.4 | P1 | confirmado |
| Galaxy publico (usegalaxy.eu/.org) via galaxy-mcp ou BioBlend | Servidores Galaxy gratuitos, com ferramentas filogeneticas instaladas e exportacao das invocacoes como Workflow Run RO-Crate. O galaxy-mcp … | Criar conta e API key no usegalaxy.eu ou no .org. - MCP: 'uvx galaxy-mcp' com GALAXY_URL e GALAXY_API_KEY no ambiente (… | galaxy-mcp 1.11.1; bioblend 1.9.0 | P1 | confirmado |
| restic (+ rclone) para backup | Backup deduplicado e criptografado, com snapshots verificaveis por 'restic check'. Binario unico para Windows e Linux. Via rclone, pode usa… | - Windows: binario oficial restic_0.19.1_windows_amd64.zip (release do GitHub), ou restic e rclone do conda-forge (win-… | restic 0.19.1; rclone 1.75.2 | P1 | confirmado |
| renv (lockfile de pacotes R no Windows) | Lockfile de pacotes R por projeto: o renv.lock registra versao e origem de cada pacote (CRAN, Bioconductor ou GitHub). | No Windows: 'install.packages("renv")' (CRAN 1.3.1), depois 'renv::init()' em analises/R/ e 'renv::snapshot()' a cada i… | 1.3.1 | P1 | confirmado |
| DataLad + git-annex (com comparacao com DVC) | Versionamento de dados sobre git + git-annex: o conteudo fica fora do git e so ponteiros entram no repositorio. 'datalad run' grava comando… | No WSL: - 'uv tool install datalad'; - git-annex pelo wheel do PyPI (10.20260901.post1) ou pelo conda-forge (10.2026071… | datalad 1.7.1 (PyPI; o conda-forge aind… | P2 | confirmado |
| Seqera Platform Cloud (Basic ou Academic Program) + tower-cli | Painel web para acompanhar e lancar runs Nextflow. - Cloud Basic, gratis: ate 3 usuarios, 3 runs simultaneos e historico de 250 runs. - Aca… | 1. Criar conta em cloud.seqera.io. 2. Pôr TOWER_ACCESS_TOKEN no ambiente do WSL. 3. Rodar 'nextflow run ... -with-tower… | tower-cli 0.42.0 (2026-10-08); no bioco… | P2 | corrigido |

Armadilhas registradas:

- pixi: o instalador poe o PATH so no ~/.bashrc (comentario do script: 'Default to bashrc'). Com 'bash -l arquivo.sh', o ~/.bashrc do Ubuntu sai cedo em shell nao interativa, entao pixi e as CLIs expostas dao 'command not found'. E a mesma falha parcial do nextflow e do seqkit. Correcao: por o PATH no bloco 'bio-stack PATH' do ~/.profile.
- pixi com o manifesto em /mnt/c cria .pixi/envs no NTFS, com milhares de arquivos pequenos (lento). Usar detached-environments = true em ~/.pixi/config.toml. Nao adicionar [tool.pixi] ao pyproject.toml: o uv e o pixi passariam a disputar o mesmo arquivo.
- O modulo nf-core IQTREE fixa '-nt AUTO -ntmax $task.cpus' e nao passa --seed; a seed padrao do IQ-TREE vem do relogio. O '-nt AUTO' vem depois de $args, entao um -T em ext.args nao e confiavel. Usar 'nf-core modules patch iqtree' ou um modulo local com --seed e -T fixos. Shen et al. 2020: 18% das arvores do IQ-TREE nao se reproduzem com seed e 2 threads iguais, e 3 ou mais threads pioram. Para co…
- Versoes divergentes no mesmo pipeline: - modulo nf-core IQTREE: iqtree 3.1.3; modulo local: 3.1.4; - trimal: 1.5.0 no modulo, 1.5.1 no bioconda; - clipkit: 2.11.4 no modulo, 2.14.0 no bioconda; - mafft: 7.520 no modulo, 7.525 no bioconda. O Metodos deve citar o que foi executado (versoes do topic), nao a ultima versao do catalogo.
- nf-prov WRROC (conferido no WrrocRenderer.groovy): - grava TODOS os params do run; - COPIA para o diretorio do crate todo arquivo local passado direto como param. Nunca passar chave de API como param (usar 'nextflow secrets'). Passar genomas e alinhamentos grandes por samplesheet, nao como param direto. O crate declara RO-Crate 1.1 e perfis WRROC 0.1.
- Plugin declarado sem versao baixa a ultima no momento do run: o run deixa de ser reproduzivel. Fixar 'nf-prov@1.7.0'. O cache local tem nf-schema 2.7.2, mas a ultima e 3.0.0 (salto de versao major).
- Containers do biocontainers rodam como root: no Docker nativo do WSL deixam arquivos de root em work/ e no outdir, e o 'nextflow clean' falha. Usar docker.runOptions = '-u $(id -u):$(id -g)', como no template nf-core.
- O publishDir publica por symlink por padrao. Links de /mnt/c para o work dir no ext4 nao abrem pelo Windows e quebram depois do 'nextflow clean'. Usar mode 'copy' em toda saida que o Windows vai ler.
- O lineage grava em ./.lineage no diretorio de lancamento. Lancado de /mnt/c, vai para o NTFS e aparece no git status, porque o .gitignore do bio-stack nao tem .lineage/. pipelines/filogenia ja tem um .nextflow/ criado em /mnt/c.
- Repo publico: nao incluir o campo 'env' em trace.fields (pode vazar variaveis de ambiente). Qualquer gravador novo no provenance.jsonl deve reutilizar o redact() do provenance.py.
- O classify() do hook marca qualquer Bash com 'python' ou 'nextflow' como analysis-run. As sondagens desta pesquisa (curl e scripts no scratchpad) entraram no log como 'analise' e vao inflar as contagens do /metodos.
- VHDX esparso: na linha 2.x, 'wsl --manage --set-sparse true' esta bloqueado por risco de corrupcao (issue #13075, exige --allow-unsafe). So a 3.0.2 pre-release reabilita, ainda como experimental. O Windows 11 Home nao tem Optimize-VHD (conferido). Compactar so com backup feito.
- 'wsl --export' leva junto o /var/lib/docker (hoje 3,9 GB de imagens + 7,0 GB de build cache). Rodar o prune antes do export.
- 'nextflow clean -f' sem -k apaga os metadados que o 'nextflow log' usa. Se o /metodos depender do log, publicar trace e report no outdir antes de limpar.
- Limites anonimos do Wave: 25 builds por dia e 100 pulls por hora; cada manifest request conta como um pull, as camadas nao. Nao construir container sob demanda a cada run. Congelar a imagem no Seqera Containers e referenciar por digest.
- Zenodo: - o .zenodo.json anula o CITATION.cff por inteiro; - so um subconjunto do CFF e lido; - a integracao com o GitHub so funciona em repo publico; - limite de 50 GB e 100 arquivos por registro.
- Repo publico sem LICENSE equivale legalmente a 'todos os direitos reservados'. Arquivar no Zenodo ou no WorkflowHub sem licenca deixa o reuso juridicamente bloqueado.
- Galaxy: IQ-TREE 3.1.3 (e nao 3.1.4), ASTRAL 5.7.8 no lugar do ASTER, e sem FAMSA. Misturar runs do Galaxy e locais exige citar as duas versoes.
- CI: o runner publico tem 14 GB de SSD. Snapshot do nf-test com md5 de .treefile so e estavel com seed e threads fixos. Preferir checar a topologia ou ignorar arquivos instaveis com o nft-utils.
- Versoes do apt no noble estao velhas: git-annex 10.20240129, restic 0.16.4, rclone 1.60.1, jq 1.7.1, skopeo 1.13.3; nao ha pacote datalad. Instalar pelo pixi/conda-forge ou pelos binarios oficiais.
- DVC: a manutencao passou da Iterative para a lakeFS em 2025-11, e o ultimo release e de 2026-03-31. Risco de ritmo mais lento para um projeto de varios anos.
- Clusters: em varios deles o apptainer so existe nos nos de computacao, entao o Nextflow precisa ser lancado de um no alocado. No SDumont, o suporte a Apptainer/Singularity nao foi verificado.
- Os creditos de pesquisa do Google Cloud so valem para docentes, doutorandos e pos-docs, nao para graduacao ou mestrado; expiram 365 dias apos o resgate.

### Infraestrutura de pesquisa e literatura

| Candidato | O que é | Rota | Versão | Prioridade | Veredito |
|---|---|---|---|---|---|
| OpenAlex API + pyalex | Índice aberto (dados CC0) de obras, autores, instituições, fontes e financiadores. O pyalex é o cliente Python (MIT), mantido pelo mesmo au… | Windows: `uv add --group lit pyalex` (e `lit` em [tool.uv] default-groups), reservado a consultas em lote ou agregadas:… | pyalex 0.21 (2026-02-23); API com cobra… | P1 | corrigido |
| Resolução e checagem de DOI: Crossref/DataCite por content negotiatio… | Receita sem instalação. `https://doi.org/<DOI>` com `Accept: application/vnd.citationstyles.csl+json` devolve CSL-JSON, tanto para DOI Cros… | Checagem pontual sem instalar nada: Crossref_check_retraction, Crossref_get_work, DataCite_get_doi e OpenCitations_get_… | Crossref REST com limites novos desde 2… | P1 | corrigido |
| Zotero 10 + Better BibTeX 9 | Gerenciador de referências desktop (Zotero 10.0.6, AGPL-3.0) com o plugin Better BibTeX 9.0.71 (MIT). O plugin dá chaves de citação estávei… | Windows: `winget install DigitalScholar.Zotero` (10.0.6 no winget em 2026-10-11). Depois, em Tools > Plugins > Install … | Zotero 10.0.6; Better BibTeX 9.0.71 (20… | P0 | confirmado |
| Quarto 1.10 (com Pandoc 3.10 e Typst 0.15.1 embutidos) | Sistema de publicação científica em Markdown + código. O tipo de projeto 'manuscript' gera HTML, DOCX, PDF e o pacote MECA, com células de … | Windows: `winget install Posit.Quarto` (1.10.19) ou o .zip, que não exige admin. No WSL, o .deb do GitHub, se for rende… | 1.10.19 (2026-10-06); pré-release 1.11.5 | P0 | confirmado |
| Estilo CSL ABNT: associacao-brasileira-de-normas-tecnicas.csl | CSL oficial do repositório citation-style-language/styles para ABNT autor-data (default-locale pt-BR). O mesmo repositório tem variantes: a… | Baixar fixando o commit em https://raw.githubusercontent.com/citation-style-language/styles/c7de5bedd4caf9c2ce193e4ccc5… | último commit c7de5be (2025-05-15); o c… | P0 | confirmado |
| nf-prov 1.7.0 (Workflow Run RO-Crate / BioCompute) + WorkflowHub | Plugin oficial do Nextflow que grava a proveniência de cada execução em formatos padrão: `wrroc` (Workflow Run RO-Crate, com os perfis Proc… | Mesmo bloco em pipelines/filogenia/nextflow.config: `plugins { id 'nf-prov@1.7.0' }` e `prov { formats { wrroc { file =… | 1.7.0 (2026-01-05; o registro de plugin… | P0 | corrigido |
| pyzotero 1.15 com CLI e servidor MCP | Cliente Python da API do Zotero (web e local) com CLI e servidor MCP opcionais. O MCP começa somente leitura; `--enable-writes` e `--enable… | `uv tool install "pyzotero[cli,mcp]"` (fora do uv.lock do projeto) e, no .mcp.json, `"zotero": {"command": "pyzotero-mc… | 1.15.2 (2026-09-18) | P1 | confirmado |
| bibverify | Verificador de referências BibTeX (CLI, API Python e servidor MCP). Procura cada entrada em Crossref, OpenAlex, Semantic Scholar, PubMed, E… | `uv tool install bibverify` (Python >=3.11, wheel puro) ou `npx @hylouis233/bibverify`. Variáveis opcionais: BIBVERIFY_… | 0.5.0 (2026-10-09) | P1 | confirmado |
| PubTator3 API (+ Europe PMC Annotations API) | API do NCBI com anotação automática de genes, espécies, químicos, doenças, variantes e linhagens celulares em PubMed e PMC-OA. As menções s… | Nada a adicionar. Usar PubTator3_EntityAutocomplete, PubTator3_LiteratureSearch, PubTator3_GetEntityRelations e PubTato… | PubTator 3.0 (serviço, sem número de ve… | descartar | **refutado** |
| LitSense 2.0 API | Busca semântica do NCBI em nível de sentença e de parágrafo, sobre resumos do PubMed e textos completos do PMC-OA. Segundo resumo do NLM: 3… | Sem instalação nem chave: GET em https://www.ncbi.nlm.nih.gov/research/litsense2-api/api/sentences/?query=<pergunta> ou… | LitSense 2.0 (NAR 2025) | P1 | confirmado |
| Unpaywall API v2 | API que localiza a melhor cópia em acesso aberto de um DOI: URL do PDF ou da landing page, licença e versão. | Nada a adicionar. Usar Unpaywall_get_full_text_url ou Unpaywall_check_oa_status do tooluniverse, ou o best_oa_location … | API v2 (consultada 2026-10-11) | descartar | **refutado** |
| ASReview LAB 3 (+ PRISMA2020; Rayyan como comparação) | Triagem de títulos e resumos por aprendizado ativo, local e de código aberto: aplicação web local (`asreview lab`), simulação reprodutível … | `uv tool install asreview` (Python >=3.10, wheel puro) e depois `asreview lab`. A entrada é RIS/CSV exportado do Zotero… | 3.0.8 (2026-06-18) | P1 | confirmado |
| Zenodo (integração com GitHub) + CITATION.cff | Repositório do CERN que arquiva cada release do GitHub e emite um DOI DataCite por versão, mais um DOI de conceito que reúne todas as versõ… | Sem instalação: ativar o repositório em zenodo.org/account/settings/github e publicar um release. Escrever o CITATION.c… | serviço InvenioRDM (consultado 2026-10-… | P1 | confirmado |
| SciELO: ArticleMeta API + prefixo DOI 10.1590 via OpenAlex/Crossref | Três vias para a coleção SciELO Brasil: a API ArticleMeta (metadados e XML SciELO PS por coleção), os DOIs Crossref do prefixo 10.1590 (mem… | Sem instalação: OpenAlex (pyalex) com filtro de prefixo + busca, e ArticleMeta em http://articlemeta.scielo.org/api/v1/… | ArticleMeta API v1 (consultada 2026-10-… | P1 | confirmado |
| BDTD e Oasisbr (IBICT): API VuFind | A Biblioteca Digital Brasileira de Teses e Dissertações e o agregador Oasisbr, ambos sobre VuFind, com API REST de busca em JSON. | Sem instalação: GET https://bdtd.ibict.br/vufind/api/v1/search?lookfor=<termos>&type=AllFields&limit=20 (o mesmo caminh… | VuFind API v1 (consultada 2026-10-11) | P1 | confirmado |
| Modelo ABNT para TCC: formato Quarto `abnt` (abnTeX2) ou Typst `abnty… | Formato Quarto para teses e dissertações ABNT (danielvartan/abnt v1.0.1) que reescreve abnTeX2 + biblatex no padrão Pandoc. Sem LaTeX, as a… | Formato `abnt`: clonar github.com/danielvartan/abnt (é um Quarto book), rodar `renv::restore()` no R e instalar o TinyT… | abnt v1.0.1 (2026-01-28); abnTeX2 1.9.7… | P1 | confirmado |
| Semantic Scholar API | API do Allen Institute sobre o grafo acadêmico S2: Graph, Recommendations (SPECTER2) e Datasets. | Nada a instalar: usar as SemanticScholar_* do tooluniverse (search_papers, get_paper, get_paper_citations, get_paper_re… | Graph API v1; semanticscholar 0.12.0 (2… | P2 | corrigido |
| Portal de Periódicos CAPES (CAFe) e Plataforma Lattes/ORCID | O Portal federal oferece conteúdo assinado (por IP institucional ou remoto pela federação CAFe/RNP) e acordos transformativos que pagam APC… | Nada a instalar. Para o Lattes, `install.packages('getLattes')` no R 4.6.1 se for processar o XML. | serviços (consultados 2026-10-11); getL… | P2 | confirmado |
| Mapas de literatura: Connected Papers, ResearchRabbit/Litmaps, Local … | Ferramentas web de exploração por grafo de citações, mais as 'citações inteligentes' do scite, que classificam se a citação apoia ou contes… | Nada a instalar (navegador). Quando o resultado precisar entrar no Métodos, preferir o Local Citation Network ou um gra… | serviços web (consultados 2026-10-11) | P2 | confirmado |
| MCPs de pesquisa avaliados: OpenAlex (oksure), arXiv (blazickjp), sci… | Servidores MCP de terceiros: openalex-research-mcp 0.5.0 (npm, TypeScript) para OpenAlex, arxiv-mcp-server 0.8.2 (PyPI) para arXiv e o MCP … | Só se o caso aparecer: `npx openalex-research-mcp` (com OPENALEX_API_KEY) ou `uvx arxiv-mcp-server`. | openalex-research-mcp 0.5.0 (2026-06-22… | descartar | **refutado** |

Armadilhas registradas:

- OpenAlex: desde 2026-02-13 o polite pool e o mailto acabaram. Sem chave, o orçamento é US$0,10/dia: busca de texto custa US$0,001, inclusive via filtro title_and_abstract.search, e repetir a mesma chamada paga de novo. Com chave grátis, US$1/dia. Lookup por DOI/ID é grátis. A tabela 'O que é aberto' do docs/CHAVES.md ficou desatualizada nesse ponto.
- Mandar a chave do OpenAlex no cabeçalho `Authorization: Bearer`, como o pyalex já faz. Com `api_key=` na query string, a chave vai em claro para o log do servidor, repetindo o incidente do NCBI_API_KEY.
- Crossref sem mailto cai no pool público. Medido em 2026-10-11: 5 req/s com concorrência 1 para DOI único e 1 req/s com concorrência 1 para listas (x-api-pool: public-array). Não paralelizar consultas de lista. O polite pool (10/3 e 3/3) exige o e-mail real do usuário.
- A referência ABNT pronta que vem do doi.org (`Accept: text/x-bibliography; style=associacao-brasileira-de-normas-tecnicas`) sai inconsistente: no Crossref faltaram páginas e DOI; no DataCite vieram tags HTML (<i>, <b>), vírgula dupla e 'Disponível em: <...>'. Baixar CSL-JSON e formatar localmente com a CSL fixada.
- O próprio arquivo CSL ABNT declara 'NBR 10520:2023 e NBR 6023:2002' e ainda põe a URL entre < >, padrão de 2002 que a NBR 6023:2018 aboliu. Nenhuma ferramenta (CSL, biblatex-abnt, abnTeX2) garante a 6023:2018 completa: conferir com o manual de normalização do Senac. As normas ABNT não são públicas; o acesso costuma vir via GEDWeb pela biblioteca.
- O DOI dos módulos nf-core fica para trás da ferramenta. O meta.yml do iqtree cita o IQ-TREE 2 (10.1093/molbev/msaa015), mas o pipeline usa o IQ-TREE 3.1.4: citar Wong et al. 2026, MBE, 10.1093/molbev/msag117. O mafft/align cita o artigo de 2002 (10.1093/nar/gkf436); para a v7, o artigo é Katoh & Standley 2013 (10.1093/molbev/mst010).
- O bio.tools erra citações por colisão de nome. 'aster' aponta para um ASTER de scATAC-seq (Chen 2022, 10.1093/bioinformatics/btac842), não para o ASTER de filogenômica (Zhang, Nielsen & Mirarab 2025, MBE, 10.1093/molbev/msaf172). 'iqtree' aponta para o IQ-TREE 1 (msu300) e 'clipkit' para 'ClipKIT in the browser'. Nunca resolver DOI de software automaticamente pelo bio.tools.
- PubTator3: a busca em texto livre casa iniciais de autor ('AMBN AND chicken' achou 'Madeira AMBN'). Usar o ID do autocomplete (`@GENE_AMBN`). O limite de taxa oficial não foi encontrado; os wrappers citam 3 req/s.
- LitSense 2.0: a API não é documentada (os endpoints foram extraídos do bundle da página). Pergunta em linguagem natural funciona, mas uma lista de palavras-chave piora o resultado. Cada chamada leva ~5 s e devolve 100 resultados. Guardar PMIDs, não trechos (direito autoral e regra do lit-scout).
- O WebFetch resume o conteúdo com um modelo pequeno. Em APIs JSON (PubTator3, LitSense, OpenAlex, BDTD) isso perde IDs, e o lit-scout não tem Bash. Dar a ele um script/Bash ou expor as APIs por um MCP próprio.
- Unpaywall: sem e-mail a resposta é HTTP 422, mas o e-mail não é validado. Os campos `evidence` e `updated` voltam 'deprecated', e `url_for_pdf` veio nulo para um artigo CC BY. Manter fallback via PMC ou locations do OpenAlex.
- Better BibTeX 9.0.71 exige Zotero >=8.0.1 (manifesto do XPI: min 8.0.1, max 10.*); o Zotero 7 perdeu suporte. Escrita pela API local exige Zotero >=10 e autorização explícita. A API local precisa ser ligada em Settings > Advanced e só escuta em localhost: do WSL, só funciona com rede espelhada.
- zotero-mcp-server: o extra [semantic] baixa ChromaDB + sentence-transformers e duplica o nvidia-nemo. O toolset scite usa endpoints públicos não documentados, e o plugin que roda dentro do Zotero é AGPL.
- Um grupo de dependências novo (ex.: `lit`) precisa entrar em `default-groups`; senão o `uv sync` puro o desinstala, a mesma armadilha já medida no pyproject. CLIs como bibverify, asreview e pyzotero-mcp ficam melhor como `uv tool`, fora do uv.lock.
- OSF: a partir de 2026-11-16 não se criam projetos novos, e em 2027-02-19 tudo vira somente leitura; só Registries e Preprints continuam. O osfclient está parado desde 2021. Não começar projeto no OSF.
- Dryad cobra US$150 até 5 GB de autor sem instituição parceira (vigente desde 2025-05-06). O Brasil não tem isenção automática, e o pedido de isenção só vale até 10 GB.
- Zenodo: se existir .zenodo.json, o CITATION.cff é ignorado por inteiro. Cada release do GitHub gera uma versão com DOI novo, então testar antes no sandbox.zenodo.org. O limite é 50 GB e 100 arquivos por registro.
- Portal CAPES: o Senac não aparece no metadado da federação CAFe (382 entidades em 2026-10-11), então não dá para contar com acesso remoto a conteúdo assinado. Os acordos transformativos exigem autor correspondente de instituição participante.
- O Google Scholar não tem API. A ferramenta bc_search_google_scholar_publications do biocontext depende do scholarly 1.7.11 (raspagem, último release em 2023-01), sujeito a CAPTCHA e bloqueio. Não usar como fonte do Métodos.
- O Quarto no conda-forge está em 1.9.38, contra 1.10.19 no GitHub e no winget. O Quarto embute o Pandoc 3.10; instalar o pandoc 3.12.1 avulso deixa duas versões no PATH, e o Métodos pode citar a errada.
- O OAI-PMH da BDTD não está configurado (404) e o da SciELO deu 502 no teste. Usar as APIs REST (VuFind, ArticleMeta) e gravar as respostas em UTF-8.
- O nf-prov 1.7.0 exige Nextflow >=25.10 (o manifest do pipeline já pede !>=25.10.0) e baixa o plugin no primeiro run; o formato `legacy` foi removido na 1.5.0. Fixar `nf-prov@1.7.0` para o Métodos citar a versão.
- Contagens de citação divergem por fonte. O IQ-TREE 2 em 2026-10-11 tinha OpenAlex 16.866, OpenCitations 14.561 e Europe PMC 13.056. Citar sempre a fonte e a data.
- scite: desde outubro de 2026 contas grátis não têm créditos de MCP. O Rayyan Free não gera diagrama PRISMA, e os dados ficam na nuvem.

### Análise ômica e estatística

| Candidato | O que é | Rota | Versão | Prioridade | Veredito |
|---|---|---|---|---|---|
| Python 3.12 no venv uv da stack (requires-python >=3.12) | Trocar o interpretador do venv Windows de CPython 3.11.17 para 3.12 e subir requires-python para '>=3.12'. É o portão para o ecossistema sc… | uv python pin 3.12; requires-python = '>=3.12'; manter scanpy==1.11.5, pydeseq2==0.5.4 e statsmodels==0.15.0; acrescent… | alvo: CPython 3.12; venv atual: 3.11.17… | P0 | corrigido |
| decoupler 2.2.0 + recursos OmniPath (CollecTRI, PROGENy, DoRothEA, Ha… | Infere atividade de TF e de vias direto em AnnData com vários métodos: dc.mt.ulm, mlm, viper, aucell, gsea, gsva, ora, zscore, waggr, mdt, … | uv add 'decoupler>=2.2' (wheel py3-none-any, requires_python >=3.11). Não usar o bioconda, que está parado em 1.5.0. | 2.2.0 (PyPI, 2026-07-11) | P0 | confirmado |
| Núcleo Bioconductor 3.23 de expressão diferencial: DESeq2, edgeR 4, l… | A referência de DE em bulk, com as implementações originais: - DESeq2 1.52.0; - edgeR 4.10.5 (quasi-likelihood); - limma 3.68.5 (voom, dupl… | Num projeto R da stack: install.packages('renv'); renv::init(bioconductor = '3.23'); renv::install(c('bioc::DESeq2','bi… | Bioconductor 3.23 (lançado 2026-04-29, … | P0 | confirmado |
| pingouin 0.7.0 + scikit-posthocs 0.17.1 | pingouin faz estatística clássica com tamanho de efeito e IC: ANOVA de uma via, de medidas repetidas e mista; ANCOVA; correlações parciais … | uv add 'pingouin>=0.7.0' 'scikit-posthocs>=0.17.1'. Os dois são wheel py3-none-any; o pingouin exige Python >=3.11 e en… | pingouin 0.7.0 (2026-09-26); scikit-pos… | P0 | confirmado |
| g:Profiler e114 (g:GOSt, g:Orth, g:Convert) | Serviço web de enriquecimento (GO, Reactome, KEGG, WikiPathways, TRANSFAC, miRTarBase, HPA) e de conversão de IDs e ortólogos, com dados do… | Sem instalação: POST para https://biit.cs.ut.ee/gprofiler/api/gost/profile/ com httpx, que já está no projeto. Em análi… | dados e114_eg62_p19_27110d83; gprofiler… | P0 | confirmado |
| NCBI Datasets v2: ortólogos (API REST e datasets CLI 18.38.0) | Conjuntos de ortólogos do NCBI Gene, gerados pelo pipeline de anotação eucariótica, consultáveis por gene_id. Cobrem centenas de genomas an… | API: GET https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/<gene_id>/orthologs via httpx; uma chave NCBI aumenta o limit… | API v2; ncbi-datasets-cli 18.38.0 (cond… | P1 | confirmado |
| scvi-tools 1.5.1 (scVI, scANVI, SysVI) | Modelos generativos em PyTorch/Lightning para single-cell: scVI e scANVI (integração e transferência de rótulo), SysVI (lotes fortes, como … | Depois da migração para 3.12: uv add 'scvi-tools>=1.5.1' (wheel py3-none-any). Não usar os extras [cuda] e [cuda13]: el… | 1.5.1 (PyPI, 2026-09-10; conda-forge 1.… | P1 | confirmado |
| pertpy 1.4.0 | Framework scverse para experimentos com várias condições. Inclui: - DE pseudobulk unificado: PyDESeq2, EdgeR, Statsmodels, Wilcoxon, t e pe… | Depois da migração para 3.12: uv add 'pertpy[de]>=1.4' (puro Python; o extra [de] traz pydeseq2 e formulaic). O extra [… | 1.4.0 (PyPI, 2026-10-02) | P1 | confirmado |
| rapids-singlecell 0.18.0 (wheels cu13) no WSL | Versão em GPU do fluxo scanpy (QC, HVG, PCA, vizinhos, Leiden, UMAP, harmony, rank_genes_groups) sobre CuPy e RAPIDS (cuML, cuGraph, cuVS),… | No WSL, venv próprio: uv venv -p 3.12 ~/rsc; source ~/rsc/bin/activate; uv pip install 'rapids-singlecell-cu13[rapids]'… | rapids-singlecell 0.18.0 (2026-10-06; w… | P1 | confirmado |
| PyMC 6.3.2 + nutpie 0.16.11 + ArviZ 1.3 + Bambi 0.21.0 | Modelagem bayesiana: - PyMC 6: PyTensor 3 com backend Numba por padrão; - nutpie: NUTS em Rust, sampler padrão quando instalado; - ArviZ 1.… | uv add 'pymc>=6.3' nutpie 'bambi>=0.21' (sem o extra bambi[nutpie]). Stan só se for necessário: no WSL, rodar sudo apt … | pymc 6.3.2 (2026-09-08), pytensor 3.3.3… | P1 | corrigido |
| Bioconductor single-cell multi-amostra: dreamlet + muscat + glmGamPoi… | Cinco pacotes para scRNA com várias amostras: - dreamlet 1.10.0: DE pseudobulk com efeitos aleatórios, sobre o dream; - muscat 1.26.0: esta… | No mesmo projeto renv do núcleo de DE: renv::install(c('bioc::dreamlet','bioc::muscat','bioc::glmGamPoi','bioc::scDblFi… | dreamlet 1.10.0, muscat 1.26.0, glmGamP… | P1 | confirmado |
| Enriquecimento em R: clusterProfiler 4.20.2 + fgsea 1.38.0 + msigdbr … | ORA e GSEA com GO, KEGG, Reactome e MSigDB, usando anotação de galinha: - clusterProfiler: enrichGO, gseGO e compareCluster com OrgDb; - fg… | renv::install(c('bioc::clusterProfiler','bioc::fgsea','bioc::enrichplot','msigdbr','bioc::org.Gg.eg.db','bioc::GO.db','… | clusterProfiler 4.20.2, fgsea 1.38.0, e… | P1 | confirmado |
| Expressão de referência para galinha: Bgee 15.2 (BgeeDB, BgeeCall) + … | Quatro recursos complementares: - Bgee: expressão normalizada e chamadas presente/ausente por órgão e estágio em 52 espécies; - BgeeCall: g… | Windows nativo: renv::install(c('bioc::BgeeCall', 'bioc::BgeeDB')) com binários 3.23 e, no BgeeCall, KallistoMetadata(d… | Bgee 15.2 (/ftp/current; pastas bgee_v1… | P1 | corrigido |
| nf-core/differentialabundance 2.0.0 | Pipeline nf-core de abundância diferencial a partir de matriz de contagem: - DESeq2, limma/voom e variancePartition/dream; - contrastes por… | Script em arquivo, executado com MSYS_NO_PATHCONV=1 e bash -l, em ~/ no WSL (não em /mnt/c): nextflow run nf-core/diffe… | 2.0.0 (2026-06-23) | P1 | confirmado |
| SAMap 3.0.1 (PyPI sc-samap) | Mapeia atlas single-cell entre espécies distantes com um grafo de homologia gênica feito por BLAST recíproco (inclui parálogos) e o manifol… | No WSL: sudo apt install build-essential python3.12-dev; uv venv -p 3.12 ~/samap; uv pip install sc-samap. Alternativa … | 3.0.1 (PyPI, 2026-03-31); o bioconda 's… | P2 | corrigido |
| CAGEE 1.2 (evolução de expressão gênica na filogenia) | Infere a evolução da expressão gênica ao longo de uma árvore de espécies com movimento browniano limitado. Estima taxas genome-wide, estado… | Para BM limitado ou ilimitado: pixi global install -c conda-forge -c bioconda cagee (1.2, linux-64). Para FreeRate: sud… | 1.2 (2024-09-04); os modelos FreeRate e… | P2 | corrigido |
| Seurat 5.6.0 | Ecossistema R de single-cell: SCTransform v2 com glmGamPoi, integração CCA/RPCA/Harmony, BPCells para matriz em disco e mapeamento por refe… | No projeto renv: install.packages('Seurat', type = 'binary'), que fica na 5.5.1 até o CRAN publicar o .zip da 5.6.0 par… | 5.6.0 no CRAN (só fonte no Windows); o … | P2 | corrigido |
| celltypist 1.7.1 | Anotação automática de tipo celular por regressão logística, com modelos pré-treinados (61 modelos no registro). | uv add celltypist (puro Python, requires_python >=3.6), só se aparecer dado humano ou de camundongo. | 1.7.1 (PyPI, 2025-06-25); registro de m… | P2 | confirmado |
| squidpy 1.8.3 / muon 0.2.0 / mudata 0.4.1 (espacial e multimodal scve… | squidpy faz análise espacial (vizinhança, autocorrelação, imagem) sobre spatialdata. mudata e muon são o contêiner e as ferramentas multimo… | Depois da migração para 3.12: uv add squidpy muon (puro Python). O squidpy traz dask e spatialdata e fixa imagecodecs<2… | squidpy 1.8.3 (2026-06-29), muon 0.2.0 … | P2 | confirmado |
| Proteômica e metabolômica leves: limpa 1.4.2 (R) + pyOpenMS 3.6.0 + m… | Três pacotes: - limpa: quantifica proteína a partir de precursores DIA-NN/Spectronaut, trata valores ausentes por curva de probabilidade de… | renv::install('bioc::limpa'); uv add pyopenms matchms. O pyopenms tem wheel win_amd64 abi3; o matchms é puro Python e e… | limpa 1.4.2 (Bioc 3.23; devel 1.5.0), p… | P2 | confirmado |

Armadilhas registradas:

- Venv em Python 3.11: o uv resolve por marcador de versão e instala versões antigas sem aviso (pertpy 1.0.3 em vez de 1.4.0, scvi-tools 1.4.2, PyMC 5.28.5, ArviZ 0.23.4, squidpy 1.8.2). Confira com 'uv pip list' depois de cada 'uv add'.
- Migrar para Python 3.12 troca, já pelo uv.lock atual, anndata 0.12.19 por 0.13.3.post0 e pandas 2.3.3 por 3.0.5, mesmo com scanpy==1.11.5 fixo. O anndata 0.13 torna X copy-on-write, remove AnnData.concatenate e __setitem__ e exige zarr ≥3; o pandas 3 muda strings e cópia; o liana exige pandas<3. Rode a bateria de testes antes de aceitar a troca.
- Python 3.13 no Windows: o rpy2-rinterface 3.6.7 publica wheels win_amd64 para cp311, cp312 e cp314, mas não para cp313.
- pingouin abaixo de 0.7.0 dava resultados errados sem erro. Casos corrigidos na 0.7.0: - rcorr com padjust contava p=0 fictícios, deixando o FDR muito liberal; - ptests ignorava padjust; - mixed_anova tinha correção GG inconsistente; - linear_regression dava erro-padrão dependente da escala do preditor; - anova e ancova eram corrompidas por níveis de fator não usados. Fixe pingouin>=0.7.0 e scikit…
- PyMC 6 e ArviZ 1.0 mudaram o intervalo de credibilidade padrão de HDI 94% para ETI 89% e trocaram InferenceData por xarray.DataTree. O nutpie, padrão quando instalado, usa tune=400, contra 1000 do NUTS do PyMC. Declare tune e tipo de intervalo nos Métodos.
- scvi-tools 1.5.1 passou a usar Adam/AdamW 'fused' por padrão e a 1.5.0 removeu JAX, então os resultados não batem bit a bit com a 1.4.x. Não instale os extras [cuda] ou [cuda13]: eles puxam torchvision, torchaudio e cupy do PyPI, fora do índice cu128 explícito do uv.
- RAPIDS só roda no WSL, e 8 GB está abaixo dos 16 GB recomendados. Windows e WSL não fazem oversubscription de memória gerenciada (rmm managed_memory não estende a VRAM). Use wheels cu13, que trazem cubins Blackwell; o cu12 recompila PTX por JIT. A rota Docker com --gpus falha hoje porque o nvidia-container-toolkit não está instalado no WSL (verificado).
- JAX no Windows é só CPU, com wheels experimentais. O pertpy[jax] (scCODA, tascCODA, scGen) roda em CPU no Windows; GPU no WSL é experimental.
- Pacotes sem wheel Windows: - harmonypy 2.1.0; - hnswlib 0.8.0 (dependência do SAMap); - tiledbsoma (cellxgene-census); - inmoose 0.9.1 (só sdist, com Cython e C++17); - rapids-singlecell-cu12/cu13 (só manylinux). No Windows o uv tenta compilar do sdist e falha sem MSVC.
- O R 4.6.1 do Windows tem 30 pacotes base, nenhum Rtools e .libPaths() aponta só para C:/Program Files/R/R-4.6.1/library. Use renv (biblioteca do projeto) e só binários. IHW, OmnipathR e zellkonverter não têm .zip Windows no VIEWS 3.23, mas são NeedsCompilation=no e instalam do fonte sem Rtools.
- O zellkonverter usa basilisk e baixa um ambiente Python próprio na primeira chamada. Para a ponte h5ad–R, use o anndataR 1.2.2, que lê h5ad e zarr nativamente.
- O bioconda está defasado: decoupler 1.5.0 (PyPI tem 2.2.0), samap 1.0.15 (PyPI sc-samap 3.0.1) e bioconductor-deseq2 1.50.2 sem win-64 (Bioc 3.23 tem 1.52.0). Fixe canal e versão em containers e nos Métodos.
- O decoupler 2.x quebrou a API 1.x (run_ulm e get_collectri viraram dc.mt.ulm e dc.op.collectri). Código de tutorial antigo ou gerado por LLM falha.
- A tradução HCOP (decoupler, e babelgene no msigdbr) assume que a regulação humana se conserva em galinha e anole, então o resultado é hipótese. O HCOP não cobre crocodilianos, psitacídeos nem tartarugas.
- O MSigDB v2026.1 existe só para humano e camundongo. O msigdbr 26.1.1 baixa os dados do Zenodo no primeiro uso e projeta para outras espécies com o babelgene 22.9, um snapshot do HCOP de 2022 — diferente do HCOP atual que o decoupler baixa.
- As fontes usam montagens diferentes de galinha: o Bgee 15.2 usa GRCg6a (Ensembl 102), o g:Profiler e114 usa GRCg7b e o NCBI usa a anotação RefSeq corrente. Coordenadas do chrZ e IDs de gene do projeto de atavismo precisam de liftover ou conversão antes de cruzar as fontes.
- Ausência de ortólogo anotado no NCBI (por exemplo, AMBN sem registro em aves e em Chrysemys) é o nível 1, o mais fraco, da escala do agente evo-devo. Nunca vira 'gene perdido' sem tBLASTn com controle positivo.
- O scanpy 1.12.x muda o cálculo do rank_genes_groups 'vs_rest' (1.12.3), corrige normalize_total entre matriz esparsa e densa (1.12.4) e passa a avisar que o teste não serve para a maioria dos usos. Trocar o pin 1.11.5 invalida as frases verificadas da skill rigor-estatistico até que sejam reverificadas.
- Em scRNA, a unidade de comparação entre condições é a amostra: use pseudobulk via decoupler com pydeseq2, edgeR ou dreamlet. rank_genes_groups e o DE por célula do scVI/scANVI tratam célula como réplica.
- O kallisto 0.52.0 não publica binário Windows, e o BgeeCall exige kallisto no PATH: rode no WSL ou em Docker. Não verifiquei se o BgeeCall 1.28.4 é compatível com kallisto ≥0.50.
- O Bioconductor 3.24 está em devel e sai em breve. Trave a 3.23 no renv para que as versões citadas nos Métodos não mudem numa reinstalação.
- Licenças GPL (pingouin, inmoose, pytximport, harmonypy, anndata2ri, BgeeCall) são aceitáveis como dependência em uso acadêmico, mas não copie código desses pacotes para o repositório que será publicado.

### Linguagens e ambientes

| Candidato | O que é | Rota | Versão | Prioridade | Veredito |
|---|---|---|---|---|---|
| pixi (ambiente WSL: conda-forge + bioconda + PyPI num lock só) | Gerenciador de ambientes conda escrito em Rust (prefix.dev). Usa manifesto pixi.toml (ou [tool.pixi] no pyproject) e lockfile pixi.lock. Te… | WSL: curl -fsSL https://pixi.sh/install.sh -o ~/pixi-install.sh, ler o script e depois rodar bash ~/pixi-install.sh. O … | 0.81.0 (release GitHub 2026-09-15; cond… | P0 | confirmado |
| Matriz de colocação Python: uv no Windows × pixi no WSL | Regra decidida pela existência de wheel win_amd64 no PyPI, conferida em 2026-10-11. Pacote com wheel Windows, ou puro, vai para o uv do Win… | Windows: libs com uv add num grupo próprio (ex.: phylo), acrescentado a default-groups. CLIs Python puras com uv tool i… | versões PyPI listadas acima, conferidas… | P0 | confirmado |
| Python 3.12 como piso comum, com pins explícitos e corte por data | Subir a venv do Windows de 3.11 para 3.12 num commit dedicado. Antes da troca, fixar os pacotes que o lock troca de major nessa fronteira. … | Fixar pandas e anndata: manter 2.3.3/0.12.19 (py3.11) ou adotar pandas==3.0.6 e anndata==0.13.4 (ou 0.13.3.post0, a ver… | uv 0.13.0 (2026-10-09) confere. Versões… | P1 | corrigido |
| R 4.6.1 no Windows: rig + pak + renv + Bioconductor 3.23 (Rtools45 so… | O rig gerencia versões de R, biblioteca de usuário, Rtools, pak e repositórios P3M. O pak instala de CRAN, Bioconductor e GitHub em paralel… | winget install --id Posit.rig -e. winget install --id RProject.Rtools -e só quando for compilar (o Rtools45 serve R 4.5… | R 4.6.1 (2026-06-24 ucrt); Bioconductor… | P0 | confirmado |
| Pacotes R para filogenia comparativa, ASR de caráter e figura de árvo… | CRAN: ape 5.8-1, phangorn 2.12.1, phytools 2.5-2, castor 1.8.7, geiger 2.0.12, corHMM 2.8, phylolm 2.6.5, caper 1.0.4, TreeTools 2.4.1 e ti… | Igual à proposta (pak dentro do renv, binários Windows R 4.6), mas com o motivo corrigido. Opcional: instalar o corHMM … | CRAN, todos com binário Windows R 4.6: … | P0 | corrigido |
| r2u (CRAN como binários Ubuntu) para R dentro do WSL | Repositório apt que entrega o CRAN e parte do Bioconductor como pacotes .deb para Ubuntu. O bspm liga o install.packages() ao apt. Declara … | Seguir o setup do r2u para noble (chave, fonte apt, pin de prioridade 700 e bspm) com o R 4.6 do repositório apt do CRA… | página do r2u: noble suportado, R 4.6.0… | P2 | confirmado |
| air + jarl (formatador e linter de R em Rust) | O air (Posit) formata código R. O jarl é um linter de R construído sobre o air: tem mais de 85 regras, correção automática, saída concisa, … | Windows: uv tool install air-formatter==0.12.0 e uv tool install jarl-linter==0.6.0. Os wheels trazem air.exe e jarl.ex… | air 0.12.0 (2026-10-01; PyPI air-format… | P0 | confirmado |
| LSP de Nextflow e R no Claude Code (Nextflow Language Server no Java … | O Nextflow Language Server é um JAR autônomo que pede Java 17+; o release v26.04.4 traz o language-server-all.jar, de 15 MB. Plugins do Cla… | Baixar language-server-all.jar do release v26.04.4 para uma pasta fixa e criar um plugin local com .lsp.json: {"nextflo… | Nextflow LS v26.04.4 (2026-09-16); exte… | P1 | corrigido |
| Nextflow 26.04 com parser estrito, nextflow lint, records e registro … | Na 26.04, o parser estrito (v2) é o padrão de nextflow run. O nextflow lint tem -format, -o json, -project-dir e -files-from. Records (reco… | WSL: fixar export NXF_VER=26.04.7 no ~/.profile (o launcher baixa a versão) ou rodar nextflow self-update. Continuar in… | 26.04.7 stable (2026-10-08); 26.09.2-ed… | P0 | confirmado |
| nf-test 0.9.5 com nft-utils | Framework de teste para processos, subworkflows e workflows Nextflow, com snapshots. A 0.9.5 acrescentou teste de canais topic, que os módu… | WSL pelo pixi (bioconda nf-test==0.9.5, noarch; o Java 21 já está presente). O instalador oficial curl -fsSL https://ge… | 0.9.5 (GitHub 2026-03-19; bioconda 0.9.… | P0 | confirmado |
| Seqera Containers / Wave para passos sem biocontainer | Registro público gratuito (community.wave.seqera.io) que constrói imagens Docker ou Singularity a partir de pacotes conda (bioconda e conda… | Gerar a imagem em seqera.io/containers e fixar a URI resultante no container do módulo local. Habilitar o Wave no nextf… | Wave v1.40.0 (2026-10-05). CLI: wave-cl… | P1 | corrigido |
| nf-prov e data lineage do Nextflow para a seção de Métodos | O nf-prov 1.7.0 gera BioCompute Object, Workflow Run RO-Crate, DAG em HTML e GEXF ao fim do run, sem alterar o script. O data lineage (expe… | No nextflow.config: plugins { id 'nf-prov@1.7.0' } com um bloco prov { formats { wrroc { file = ... } bco { file = ... … | nf-prov 1.7.0 (2026-01-05); lineage exp… | P1 | confirmado |
| Snakemake 9 como alternativa (não adotar) | Motor de workflow em Python. A documentação recomenda instalar pelo pixi e, no Windows, rodar dentro do WSL. | Se um projeto externo exigir: pixi no WSL (bioconda snakemake 9.27.0). | 9.27.0 (GitHub 2026-09-11; PyPI e bioco… | P2 | confirmado |
| Quarto 1.10 para relatório, Métodos e manuscrito | Sistema de publicação que renderiza .qmd, .ipynb e .md para HTML, DOCX e PDF. Inclui o Pandoc e o CLI do Typst, o que permite PDF sem LaTeX… | winget install --id Posit.Quarto -e. Definir QUARTO_PYTHON apontando para .venv/Scripts/python.exe do repositório, porq… | 1.10.19 (release no GitHub e no winget … | P0 | confirmado |
| marimo (notebook reativo salvo como .py) | Notebook em que as células formam um grafo de dependências, sem estado oculto, salvo como .py puro que também roda como script. O marimo ch… | Windows: uv add --group notebooks marimo==0.25.1, com o grupo em default-groups. Notebooks avulsos: uvx marimo edit --s… | 0.25.1 (PyPI 2026-10-01; conda-forge 0.… | P1 | corrigido |
| jupytext + nbstripout para quem fica no Jupyter | O jupytext pareia cada .ipynb com um .py no formato percent, que vai para o git e para o lint. O nbstripout remove as saídas do .ipynb no c… | Windows: uv add --group notebooks jupytext==1.19.6 nbstripout==0.9.1. Em cada notebook, jupytext --set-formats ipynb,py… | jupytext 1.19.6 (2026-10-04); nbstripou… | P1 | confirmado |
| piqtree (IQ-TREE 3 como biblioteca Python) | Bindings Python oficiais do IQ-TREE, mantidos pela organização iqtree, com interface por objetos do cogent3. Permite reconstruir árvores e … | No Windows, só depois da venv em Python 3.12 (uv add piqtree==0.8.3). No WSL, em [pypi-dependencies] do pixi. Não está … | 0.8.3 (PyPI 2026-06-09); embute o IQ-TR… | P2 | confirmado |
| Julia + PhyloNetworks.jl / SNaQ.jl | Julia 1.13 instalada via juliaup. O PhyloNetworks 1.4.0 trata de redes filogenéticas. O SNaQ 1.2.1 infere redes por pseudo-verossimilhança … | Se o caso aparecer: winget install --id Julialang.Juliaup -e, depois juliaup add 1.13 e Pkg.add de PhyloNetworks, SNaQ … | Julia 1.13.1 (2026-09-26; LTS 1.10.12);… | P2 | confirmado |
| Go e Rust só como consumidores: binários estáticos no Windows | Não escrever código em Go nem em Rust; consumir binários prontos. seqkit 2.14.0, csvtk 0.38.0, gotree 0.5.2 e goalign 0.4.1 publicam .exe p… | Baixar seqkit_windows_amd64.exe.tar.gz e csvtk_windows_amd64.exe.tar.gz e conferir o .md5.txt publicado. gotree_v0.5.2_… | seqkit 2.14.0 e csvtk 0.38.0 (2026-09-2… | P1 | corrigido |
| DuckDB + Parquet como camada SQL e formato de troca entre Python e R | O DuckDB 1.5.6 tem wheel Python para Windows, pacote R com binário Windows e CLI no conda-forge. Ele consulta CSV, TSV, JSONL e Parquet dir… | Windows: uv add duckdb==1.5.6. R: pak::pkg_install(c('duckdb','arrow')) no renv. WSL: duckdb-cli==1.5.6 no pixi. | 1.5.6 (PyPI 2026-09-28; CRAN 2026-09-29… | P1 | confirmado |

Armadilhas registradas:

- O guard_data.py bloqueia curl/wget com pipe para bash/sh e iwr/Invoke-WebRequest com pipe para iex. Por isso os instaladores de uma linha de pixi, juliaup, nf-test (get.nf-test.com), uv e air-installer.ps1 falham quando o Claude os executa. Rota: baixar para arquivo, ler e executar; ou usar winget (Posit.Quarto 1.10.19, Posit.rig 0.11.0, RProject.Rtools, Julialang.Juliaup), scoop (pixi 0.81.0, ju…
- O instalador do pixi grava o PATH só no ~/.bashrc. A stack chama o WSL com bash -l arquivo.sh, que lê o ~/.profile, e o .bashrc do Ubuntu desiste em shell não interativa. O resultado é pixi: command not found com exit 127, a mesma falha parcial da armadilha 3 do CLAUDE.md. Pôr o export PATH=$HOME/.pixi/bin:$PATH no ~/.profile e conferir o mesmo em qualquer outro instalador.
- uv sync ou uv venv rodados no WSL dentro de /mnt/c/Users/Manec/bio-stack recriam a .venv com layout Linux e destroem a venv do Windows. A dica de RAPIDS impressa pelo bootstrap/02_python_env.ps1 (uv venv --python 3.11 && source .venv/bin/activate) faz exatamente isso se for executada na raiz do repositório. No WSL, trabalhar em ~/ ou exportar UV_PROJECT_ENVIRONMENT=$HOME/.venvs/bio-stack.
- Um env pixi ou conda dentro de /mnt/c significa milhares de arquivos pequenos atravessando a fronteira NTFS: fica lento e sujeito a erro de permissão e de symlink. Usar detached-environments = true, que põe o env no ext4, e deixar no repositório só o pixi.toml e o pixi.lock.
- O uv.lock já bifurca na fronteira do Python 3.12. Com 3.11 resolve pandas 2.3.3, anndata 0.12.19 e biotite 1.6.0; com 3.12 ou mais, pandas 3.0.5, anndata 0.13.3 e biotite 1.7.1. Trocar o Python da venv, ou criar uma venv 3.12 no WSL a partir do mesmo lock, muda o major do pandas sem nenhum aviso.
- O polars 2.0.0 saiu em 2026-10-06, e o pyproject aceita essa versão (polars>=1.44), então o próximo uv lock --upgrade puxa a 2.0. O polars-bio 0.36.1 exige polars<2 e pyarrow<25. Fixar polars>=1.44,<2 até o ecossistema acompanhar.
- O ruff 0.17.0 saiu em 2026-10-09, e a stack fixa ruff==0.16.7. Subir de propósito, porque regras e formatação mudam entre minors.
- Não têm wheel win_amd64: pysam 0.24.1, cyvcf2 0.34.0, pybedtools 0.12.1, pytrimal 0.8.5, pyBigWig 0.3.26, HTSeq 2.1.2 e ruranges/gtfreader (dependências do pyranges1). Só têm sdist: ete4 4.4.0 e mappy 2.31. No Windows o uv tenta compilar e falha sem MSVC. No WSL os que têm wheel manylinux instalam, mas ete4 e mappy falham porque não há gcc nem make. Bioconda/conda-forge pelo pixi resolve.
- Exigem Python 3.12 ou mais: piqtree 0.8.3, biotite 1.7.1 e pyranges1 1.4.5 (via ruranges). Nenhum deles entra na venv 3.11 do Windows.
- Há três versões do IQ-TREE ao alcance: o piqtree 0.8.3 embute a 3.1.2, o módulo nf-core iqtree usa o container 3.1.3 e o bioconda está na 3.1.4. Fixar uma versão por projeto (override do container do módulo ou pin no pixi) e registrar qual gerou cada árvore.
- O lint_gate.py está inerte para Nextflow e para R no Windows, sem aviso. No Nextflow, have('nextflow') é verdadeiro só porque o uv existe, mas o hook chama o nextflow direto; como ele não está no PATH do Windows, o OSError vira lista vazia. No R, o R do Windows não tem lintr, o erro vai para stderr e o hook lê só stdout.
- O lint_gate.py não cobre .ipynb (o NotebookEdit cai fora do DISPATCH), então os notebooks escapam de ruff e basedpyright. O marimo (.py) ou o par .py do jupytext resolvem.
- Colisão de nomes: no PyPI, air é um framework web Python (Air 0.48.1); no bucket main do scoop, air é o live-reload do Go (1.67.4). O formatador R da Posit é air-formatter no PyPI (executável air) e air no conda-forge. uvx air e scoop install air instalam a coisa errada, e o fallback tool_cmd do lint_gate cairia nisso.
- Canais conda atrasados: (1) r-phytools no conda-forge parado em 1.2_0 (2022), contra 2.5-2 no CRAN; (2) bioconductor-ggtree 4.0.5 e bioconductor-treeio 1.34.0 no bioconda, contra 4.2.0 e 1.36.1 no Bioconductor 3.23; (3) csvtk 0.31.0 no bioconda, contra 0.38.0 no GitHub; (4) quarto 1.9.38 no conda-forge, contra 1.10.19; (5) phykit 2.7.0 no bioconda, contra 2.8.0 no PyPI. O R de análise vem do CRAN…
- O bioconda só funciona com o conda-forge antes dele e prioridade estrita. Nunca listar o bioconda primeiro no pixi.toml.
- No Windows, o Quarto escolhe o Python pelo launcher py. Aqui existem C:/Python314 e um Python 3.12 de usuário, e nenhum deles tem o Jupyter da stack. Sem QUARTO_PYTHON apontando para .venv/Scripts/python.exe, o render com motor jupyter falha ou usa pacotes errados. Cada .qmd tem um único motor: R e Python no mesmo documento só via knitr + reticulate.
- Nextflow 26.04: o parser estrito é o padrão de nextflow run, e NXF_SYNTAX_PARSER=v1 reabilita o antigo. O hash de entradas map passou a incluir as chaves, então processos com meta map reexecutam no primeiro -resume depois de subir da 25.10. Tipos estáticos são preview (nextflow.enable.types = true em cada script), e a sintaxe mudou entre a 25.10 e a 26.04; ficam fora do pipeline da tese.
- O nf-schema 3.0.0 (2026-09-29) quebra a API: paramsHelp passa a usar paramsHelp(parameter: ...), paramsSummaryMap() e paramsSummaryLog() deixam de receber o objeto workflow, e failUnrecognisedParams e failUnrecognisedHeaders foram removidos. O template do nf-core tools 4.1.0 fixa nf-schema@2.5.1 e nextflowVersion !>=25.10.4; não subir o plugin num pipeline gerado pelo template.
- O Wave anônimo permite 25 builds por dia e 250 pulls por hora. As imagens do Seqera Containers são públicas, então nada de dado ou código privado dentro delas. O repositório seqeralabs/wave-cli foi arquivado; a CLI agora sai com tags cli-v* no seqeralabs/wave.
- O R 4.6 usa o Rtools45 (não existe Rtools46; no winget, RProject.Rtools 4.5.6768). Ele só é necessário para pacote compilado vindo do GitHub, como o RERconverge (Rcpp/RcppArmadillo). No CRAN, o binário Windows do renv ainda é o 1.2.4, enquanto o fonte já é o 1.3.1.
- Em modo implícito, o renv só grava no renv.lock os pacotes que aparecem em código R do projeto. Num repositório que ainda não tem .R, o renv::snapshot() sai vazio. Usar renv::snapshot(type = 'all') ou ajustar o snapshot.type.
- O Bioconductor 3.23 é o release (2026-04-29, R 4.6). O 3.24 está em devel sem data confirmada; o histórico aponta fim de outubro. Fixar o 3.23 no renv agora; a troca para o 3.24 não exige trocar o R.
- O ete4 renderiza imagem pelo PyQt6 (extra treeview) ou pelo selenium (smartview). No WSL sem display, isso exige QT_QPA_PLATFORM=offscreen e bibliotecas Qt do sistema. Para figura de publicação, usar ggtree ou toytree.
- O basedpyright do Windows marca import pysam e import cyvcf2 como ausentes em scripts escritos para o WSL. Com reportMissingImports=warning o lint_gate não bloqueia, mas o LSP perde os tipos. Separar os scripts do WSL numa pasta própria.

### Genômica comparativa e aquisição de dados

| Candidato | O que é | Rota | Versão | Prioridade | Veredito |
|---|---|---|---|---|---|
| NCBI Datasets CLI (datasets + dataformat) | CLI oficial do NCBI que baixa e resume genomas, genes, ortólogos (NCBI Orthologs) e taxonomia em pacotes zip com relatório JSONL. O datafor… | Windows: conda-forge ncbi-datasets-cli (tem win-64) ou os executáveis datasets.exe e dataformat.exe de https://ftp.ncbi… | 18.38.0 (2026-09-28) | P0 | confirmado |
| miniprot | Alinha proteína contra genoma com splicing e frameshift. O GFF3 de saída traz, por alinhamento, os atributos Frameshift (número de frameshi… | WSL: bioconda miniprot 0.18 via pixi, ou docker quay.io/biocontainers/miniprot:0.18--h577a1d6_0, ou make (só depende de… | 0.18 (2025-07-11) | P0 | corrigido |
| Utilitários UCSC (kent) — twoBitToFa, faToTwoBit, mafsInRegion, mafSp… | Utilitários de linha de comando do Genome Browser para 2bit, MAF, cadeias e nets, liftOver, bigWig e checagem de hubs. | WSL: binários pré-compilados v504 de https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/ (download individual ou rs… | kent v504 (binários de 2026-10-06); bio… | P0 | corrigido |
| Recursos UCSC hospedados: GenArk + API REST + multiz77way/phyloP77way… | Quatro recursos hospedados. (1) GenArk: um hub por accession GCA/GCF, com 2bit soft-masked, RepeatMasker, genes e chromAlias. (2) A API RES… | Nenhuma: é HTTP. Teste verificado: https://api.genome.ucsc.edu/getData/sequence?genome=GCF_016699485.2;chrom=NC_052572.… | GenArk dataTime 2026-10-06; multiz77way… | P0 | corrigido |
| TaxonKit + NCBI taxdump | Toolkit local para a taxonomia NCBI: lineage, reformat, list (descendentes), name2taxid, lca, filter e rastreio de mudanças de taxid. | Windows: binário taxonkit_windows_amd64.exe.tar.gz da release v0.20.0. WSL: bioconda taxonkit ou biocontainer taxonkit:… | 0.20.0 (2025-04-29) | P0 | confirmado |
| gffread + AGAT | O gffread extrai CDS, proteína e transcritos a partir de GFF3+FASTA e converte GFF↔GTF. O AGAT é uma suíte Perl que padroniza GFF3, guarda … | WSL: bioconda gffread e agat (noarch), ou os biocontainers gffread:0.12.9--hf426362_0 e agat:1.7.0--pl5321hdfd78af_0. E… | gffread 0.12.9 (2026-04-10); AGAT 1.7.0… | P0 | corrigido |
| Coleções de genomas: VGP/GenomeArk, B10K, DNA Zoo (tabela de montagen… | Fontes das montagens de referência. O VGP publica no GenomeArk, com bucket S3 público. O B10K (Bird 10K) tem genomas de aves e os alinhamen… | Nenhuma instalação; acesso por accession com o datasets CLI e pelos buckets S3 públicos. Antes de uma análise genome-wi… | Estado em 2026-10-11. VGP no NCBI (PRJN… | P0 | corrigido |
| compleasm + BUSCO 6 (odb12) + OMArk | Medem completude de montagem e anotação com ortólogos de cópia única do OrthoDB v12 (BUSCO). O compleasm reimplementa o BUSCO sobre minipro… | WSL: biocontainers compleasm:0.2.10--pyhdfd78af_0 e busco:6.1.0--pyhdfd78af_2, ou os módulos nf-core compleasm/run e bu… | compleasm 0.2.10 (2026-10-07, suporta o… | P0 | confirmado |
| TOGA2 (+ resultados TOGA pré-computados do Hiller lab) | Projeta a anotação e infere ortologia a partir de cadeias de alinhamento genômico, com CESAR2 e um classificador XGBoost. Cada projeção, tr… | WSL: 'git clone --recurse-submodules https://github.com/hillerlab/TOGA2' e depois 'make', dentro de um venv com Python … | TOGA2 v2.0.10 (2026-09-10); TOGA1 v1.1.… | P1 | confirmado |
| make_lastz_chains (LASTZ) + nf-core/pairgenomealign (LAST) | Dois pipelines Nextflow de alinhamento par a par de genomas inteiros. O make_lastz_chains, do Hiller lab, encadeia LASTZ, axtChain, RepeatF… | WSL. make_lastz_chains: git clone na tag v3.1.7 e 'nextflow run main.nf -params-file params.json -profile docker', com … | make_lastz_chains v3.1.7 (2026-08-09); … | P1 | corrigido |
| GENESPACE + MCScanX | Sintenia entre vários genomas ancorada em proteína (OrthoFinder + MCScanX), com gráficos riparian, acompanhamento de regiões de interesse e… | WSL: R via CRAN Ubuntu ou pixi r-base, depois devtools::install_github('jtlovell/GENESPACE'). OrthoFinder 2.5.5 em ambi… | GENESPACE 1.3.1 no DESCRIPTION do maste… | P2 | corrigido |
| nf-core/genomeqc | Pipeline nf-core que baixa genomas e anotações do NCBI e compara qualidade: BUSCO com ideograma, QUAST, AGAT, sobreposição de genes, isofor… | WSL: 'nextflow run nf-core/genomeqc -r 1.0.1 -profile docker --busco_lineage sauropsida_odb12', sem passar --gxdb. O li… | 1.0.1 (2026-09-16) | P1 | confirmado |
| LiftOn (com EGAPx e GeMoMa como referência) | Transfere anotação entre montagens combinando Liftoff (DNA) e miniprot (proteína) com um algoritmo de maximização de proteína. Desde a v1.0… | WSL: 'uv pip install lifton==1.0.14' num venv, com minimap2 e miniprot do bioconda via pixi, ou o biocontainer lifton:1… | LiftOn 1.0.14 (2026-09-25); Liftoff 1.6… | P1 | corrigido |
| Ortologia pré-computada: OrthoDB v12.2 (API) + OMA (API/omadb) + eggN… | Bancos de ortologia hierárquica consultáveis por REST. O OrthoDB tem grupos por nível taxonômico, perfil filético e FASTA. O OMA tem HOGs e… | Os REST não precisam de instalação; o httpx já está instalado. Windows: 'uv add omadb' (2.2.7, puro Python). eggNOG-map… | OrthoDB v12.2 (endpoint orthodb_release… | P1 | confirmado |
| TimeTree 5 (API) + Open Tree of Life (API v3, rotl/opentree) | O TimeTree dá tempos de divergência sintetizados de estudos publicados (pairwise, timeline, timetree). O Open Tree dá a árvore sintética, o… | TimeTree: 'curl https://timetree.temple.edu/api/pairwise/<taxidA>/<taxidB>', que devolve CSV. Open Tree: no R do Window… | TimeTree 5 (API viva); Open Tree síntes… | P1 | corrigido |
| nf-core/fetchngs + ENA Portal API | Pipeline nf-core que baixa metadados e FASTQ do SRA, ENA, DDBJ e GEO por accession e gera samplesheet para outros pipelines nf-core. A ENA … | WSL: 'nextflow run nf-core/fetchngs -r 1.13.0 -profile docker --input ids.csv'. Para descobrir runs: curl 'https://www.… | fetchngs 1.13.0 (2026-09-15); ENA Porta… | P1 | confirmado |
| pyGenomeViz (com clinker) | Visualiza genomas e loci comparativos em Python a partir de GenBank ou GFF, com links entre faixas, e exporta PNG, SVG, PDF ou HTML. Tem GU… | Windows: 'uv add pygenomeviz' (1.7.0, wheel py3-none-any; as dependências biopython, matplotlib e numpy já estão no amb… | 1.7.0 (2026-06-13); clinker 0.0.32 (202… | P1 | confirmado |
| Ensembl REST (Compara, release 116, congelado) + GraphQL da nova plat… | A API REST do Ensembl expõe o Compara: árvores de genes, homologias, sintenia e alinhamentos EPO e Pecan. Ela está congelada na release 116… | Nenhuma: HTTP com httpx. Registre 'Ensembl 116' no provenance.jsonl. Para dados novos, use https://www.ensembl.org/data… | REST release 116 (rest.ensembl.org/info… | P1 | confirmado |
| Progressive Cactus + HAL (alinhamentos B10K 363 aves / 605 vertebrado… | Alinhador múltiplo de genomas inteiros guiado por árvore, com saída em HAL e as ferramentas hal2maf, halLiftover, halSynteny e cactus-phast. | WSL: bioconda cactus 3.3.0 (linux-64), binários Linux ou imagem Docker da release v3.3.0, só para operações em HAL (hal… | Cactus v3.3.0 (GitHub 2026-08-11; bioco… | P2 | corrigido |
| ntSynt (+ ntSynt-viz) e SyRI/plotsr | Sintenia nucleotídica. O ntSynt encontra blocos macrossintênicos em vários genomas por grafos de minimizadores. O SyRI e o plotsr detectam … | WSL: bioconda ntsynt 1.0.9, ntsynt-viz 1.1.4, syri 1.8.2 e plotsr 1.2.0, ou o biocontainer ntsynt:1.0.9--py312h9b2995c_… | ntSynt 1.0.9 (2026-08-20); ntSynt-viz 1… | P2 | confirmado |

Armadilhas registradas:

- A montagem marcada como 'reference' no NCBI nem sempre é a melhor. Para Gavialis gangeticus, --reference (ou filters.reference_only) devolve o GavGan_comp1 (GCF_001723915.1, de 2016, contig N50 ~10 kb), mas existe o rGavGan2.hap2 do VGP (GCA_030020295.1, contig N50 84,2 Mb, anotação VGL de 2026-01). Crocodylus porosus não tem montagem melhor que a CroPor_comp1; use C. niloticus rCroNil1 (GCA_0548…
- Há um erro de accession em docs/PROJETO-ARCHOSAURIA.md (§2, grupos externos): GCF_007399415.2 é Gopherus evgoodei (rGopEvg1_v1.p), não Chrysemys picta. A tartaruga-pintada é GCF_011386835.1 (ASM1138683v2, RefSeq 2024-05). Conferido na API v2 em 2026-10-11.
- Nove scripts de scripts/atavismo e o §2 do PROJETO-ARCHOSAURIA chamam /datasets/v2alpha/. Hoje ainda responde 200, mas a versão estável é /datasets/v2/ desde 2024, e versões aposentadas passam a devolver HTTP 410. Migre antes que quebre.
- Estar no NCBI Orthologs ou no RefSeq não prova gene íntegro. ENAM aparece como PROTEIN_CODING em Chlamydotis macqueenii (XP_086627479.1, 263 aa) e em Chrysemys, Trachemys e Malaclemys (253–505 aa), contra 1092 aa no jacaré: são fragmentos. Antes de chamar de 'presente', valide com o miniprot (StopCodon/Frameshift) ou com o TOGA2 (classe L ou UL).
- O Ensembl 116 é a última release no modelo antigo: APIs, MySQL e FTP estão congelados, beta.ensembl.org já redireciona para www.ensembl.org, e o GraphQL novo não expõe homologias nem árvores. O Compara 116 usa montagens antigas (C. porosus CroPor_comp1, Chrysemys 3.0.3, Strigops bStrHab1_v1.p, avestruz ASM69896v1) e não tem Alligator mississippiensis. Coordenadas e IDs não batem com as montagens …
- O liftOver da UCSC foi feito para converter entre montagens do mesmo organismo; entre espécies, use cadeias e nets ou o halLiftover. O liftOver e os arquivos .chain só podem ser usados sem fins comerciais, o que serve para trabalho acadêmico.
- galGal6 (GRCg6a) não é GRCg7b. As coordenadas das relíquias (GRCg7b, NC_052572.1) precisam passar pelo GCF_016699485.2ToGalGal6.over.chain.gz antes de recortar o multiz77way ou o phyloP77way. Os nomes de cromossomo também diferem (NC_052572.1 contra chrZ); use o chromAlias da GenArk.
- Os pacotes ucsc-* do bioconda estão na versão 482 (2025), enquanto os binários oficiais estão na v504 (2026-10-06). Para comandos novos, use os binários do hgdownload.
- Pelo README, o GENESPACE só funciona com OrthoFinder 2.5.x; o OrthoFinder 3 é incompatível. Isso conflita com o OrthoFinder 3.1.5 previsto no pipeline de filogenia, então isole o GENESPACE em ambiente ou container próprio (orthofinder 2.5.5).
- O TOGA2 não tem pacote conda nem módulo nf-core. Exige Python ≥3.13 (o WSL tem 3.12.3), gcc, cargo e Nextflow. As cadeias precisam ter como alvo a referência. O parser chaintools quebra se faltar a linha em branco entre cadeias ou se houver linhas '#'. O score mínimo padrão (15.000) pode descartar cadeias de montagens fragmentadas.
- Os alinhamentos de códon do TOGA mascaram toda mutação inativadora como XXX. Servem para HyPhy e PAML, mas não para contar lesões; no atavismo, use geneInactivatingMutations.tsv e loss_summ_data.tsv. O conjunto pré-computado com referência galGal6 não tem crocodilianos nem genes ausentes no frango (ENAM, AMBN, AMELX...).
- O make_lastz_chains exige genomas soft-masked (minúsculas), nunca hard-masked. A GenArk entrega sequência soft-masked (medido: 30% de minúsculas em 100 kb do chr1 de GRCg7b). FASTA de outras fontes pode vir sem máscara e fazer o tempo do LASTZ explodir.
- O miniprot não é otimizado para homólogos distantes e só indexa ORFs acima de -L (padrão 30). Uma relíquia com stop a cada poucos códons pode perder as sementes. Localize primeiro com DIAMOND ou tBLASTn (mascaramento desligado só dentro da janela sintênica, regra que o projeto já aprendeu) e refine com o miniprot na janela.
- No nf-core/genomeqc, o busco_lineage padrão do config é hymenoptera_odb10: passe sauropsida_odb12. Não forneça --gxdb, porque o FCS-GX exige 470 GiB de disco e 512 GiB de RAM.
- BUSCO odb12: não existe linhagem de Crocodylia nem de Testudines, então use sauropsida_odb12 para todas as espécies. odb12 e odb12.2 diferem (segundo o compleasm, o odb12 não traz lengths_cutoff); registre a versão do dataset no provenance.jsonl.
- Sintenia por k-mer (ntSynt, minimap2, SyRI) não atravessa ave×crocodiliano (~243 Ma pelo TimeTree); o preset máximo do ntSynt é '>10%' de divergência. Use âncoras proteicas (GENESPACE/MCScanX) ou cadeias LASTZ/LAST.
- Tamanhos dos alinhamentos B10K: 363-avian-2020.hal tem 417 GB e só aves; 605-vertebrate-2020.hal tem 1,19 TB, sem crocodilianos nem tartarugas; o MAF single-copy tem 232,6 GB. O multiz77way do galGal6 (chrZ com 2,7 GB) já traz 2 jacarés, 4 tartarugas e 3 escamados.
- O OrthoDB só tem alguns níveis: Sauropsida (8457) e Aves (8782) existem; Archosauria, Sauria e Archelosauria não. Nesses níveis a busca devolve count 0 com status ok, o mesmo 'zero silencioso' do BLAST web. Confira se o nível existe antes de ler zero como ausência.
- O OMA All.May2026 não tem Alligator, Gavialis nem Struthio, e o CROPO usa a montagem antiga do C. porosus.
- O eggNOG-mapper v3 (beta, eggNOG 7) não é compatível com os bancos do eggNOG 5 e muda as colunas de saída. Para reprodutibilidade, fixe a 2.1.15.
- Rodam só em Linux: jcvi (sem wheel Windows), LiftOn (por causa do pysam), ete4 (só sdist; o conda-forge não tem win-64), além de miniprot, LASTZ e MCScanX (o bioconda não tem win-64). O pygenomeviz do bioconda está parado em 0.4.4; use a 1.7.0 do PyPI.
- GenomeArk: o domínio raiz genomeark.org não resolve. Use https://www.genomeark.org/ ou o bucket S3 'genomeark'.
- A API do TimeTree devolve CSV (text/plain), não JSON. Os tempos são síntese de estudos publicados e não substituem datação própria.

### Ferramentas de linha de comando no WSL

| Candidato | O que é | Rota | Versão | Prioridade | Veredito |
|---|---|---|---|---|---|
| pixi 0.81.0 (modo pixi global): gerenciador das CLIs de bio no WSL | Gerenciador de pacotes conda escrito em Rust (binário único, da prefix.dev). O modo `pixi global` instala cada ferramenta num ambiente isol… | Script em arquivo (armadilha 1), executado com MSYS_NO_PATHCONV=1 wsl -d Ubuntu-24.04 -- bash -l <script>. Nada de curl… | 0.81.0 | P0 | confirmado |
| jq 1.8.2 (+ yq do mikefarah, pacote go-yq 4.54.1, como par para YAML/… | Processador de JSON e JSON Lines na linha de comando, em C e sem dependências. | pixi global install -c conda-forge "jq==1.8.2" (ou a entrada [envs.jq] do manifesto). Para o P1: pixi global install -c… | 1.8.2 | P0 | confirmado |
| DuckDB CLI 1.5.6 (duckdb-cli; extensão duckhts como P2) | Banco analítico embutido, colunar e vetorizado. Roda SQL direto sobre CSV/TSV, Parquet e JSON(L) sem carga prévia, com spill em disco quand… | pixi global install -c conda-forge "duckdb-cli==1.5.6" (o executável é `duckdb`, conferido no teste do feedstock duckdb… | 1.5.6 | P0 | confirmado |
| seqkit 2.14.0 (migração do binário solto 2.13.0 para o pixi) | Canivete suíço de FASTA/FASTQ em Go: stats, grep, subseq, sliding, split2, fx2tab, translate e afins. | pixi global install -c conda-forge -c bioconda "seqkit==2.14.0" e depois rm ~/.local/bin/seqkit. Conferir com `type -a … | 2.14.0 | P0 | confirmado |
| samtools / bcftools / htslib 1.24 | Suíte htslib para SAM/BAM/CRAM, VCF/BCF, faidx, bgzip e tabix. | Ambiente [envs.hts] com samtools, bcftools e htslib ==1.24 (canais conda-forge e bioconda), expondo samtools, bcftools,… | 1.24 | P0 | confirmado |
| BLAST+ 2.17.0 (tblastn local) | Suíte BLAST do NCBI: makeblastdb, blastn, blastp, blastx, tblastn, tblastx e blastdbcmd. | pixi global install -c conda-forge -c bioconda "blast==2.17.0". Traz entrez-direct, ncbi-vdb e perl como dependências, … | 2.17.0 | P0 | confirmado |
| DIAMOND 2.2.8 | Alinhador de proteínas acelerado (blastp e blastx, com alinhamento sensível a frameshift via -F) e fluxos de clusterização (cluster, linclu… | pixi global install -c conda-forge -c bioconda "diamond==2.2.8" | 2.2.8 | P0 | confirmado |
| miniprot 0.18 (alinhamento proteína-genoma com splicing e frameshift) | Alinha proteínas a um genoma inteiro com íntrons e frameshifts. Por alinhamento, reporta o número de frameshifts (tag PAF fs:i, GFF Framesh… | pixi global install -c conda-forge -c bioconda "miniprot==0.18". Opcional: minisplice 0.4 (bioconda, 2025-07) para esco… | 0.18 | P0 | confirmado |
| ncbi-datasets-cli 18.38.0 (datasets + dataformat) | CLI oficial do NCBI Datasets (API v2). Baixa genoma, proteína, CDS, GFF e relatórios em pacotes, desidratados ou não, com `rehydrate`, e co… | pixi global install -c conda-forge "ncbi-datasets-cli==18.38.0" (expõe datasets e dataformat). Não está no bioconda. Ba… | 18.38.0 | P0 | confirmado |
| GNU parallel 20260922 | Executor de jobs em paralelo para a shell, com joblog, retomada, novas tentativas e limites de carga e memória. | pixi global install -c conda-forge "parallel==20260922" (noarch, Perl). O apt oferece 20231122. | 20260922 | P0 | confirmado |
| csvtk 0.38.0 (conda-forge) | Toolkit de CSV/TSV com cabeçalho para pipes: join, filter2, mutate, summary, pretty, csv2html. | pixi global install -c conda-forge "csvtk==0.38.0". Em ambiente com canais conda-forge e bioconda, o conda-forge vence. | 0.38.0 | P1 | confirmado |
| HMMER 3.4 (hmmbuild/hmmsearch/nhmmer) + pyhmmer 0.12.3 no Windows | Busca de homologia por HMM de perfil, em proteína (hmmsearch, hmmscan) e em DNA (nhmmer). | pixi global install -c conda-forge -c bioconda "hmmer==3.4", expondo hmmbuild, hmmsearch, hmmscan, hmmpress e nhmmer. O… | 3.4 (pyhmmer 0.12.3) | P1 | confirmado |
| MMseqs2 18-8cc5c (CPU e modo GPU) | Suíte de busca e clusterização de proteínas e nucleotídeos. O modo GPU usa `makepaddedseqdb` e depois `search --gpu 1` (ou gpuserver). | pixi global install -c conda-forge -c bioconda "mmseqs2==18.8cc5c" (expõe mmseqs; traz aria2 como dependência). Primeir… | 18-8cc5c (bioconda 18.8cc5c) | P1 | confirmado |
| minimap2 2.31 (+ bwa 0.7.19 como alternativa viável ao bwa-mem2) | Alinhador genérico: leituras curtas e longas, cDNA com splicing e genoma contra genoma (asm5/asm20), com paftools. | pixi global install -c conda-forge -c bioconda "minimap2==2.31" | 2.31 | P1 | confirmado |
| sra-tools 3.4.1 + mosdepth 0.3.14 (verificação em nível de leitura; f… | sra-tools: prefetch, fasterq-dump e vdb-config, para baixar leituras do SRA. mosdepth: cobertura por base ou por região a partir de BAM/CRA… | pixi global install -c conda-forge -c bioconda "sra-tools==3.4.1" e "mosdepth==0.3.14". O post-link do sra-tools (que c… | sra-tools 3.4.1; mosdepth 0.3.14 | P1 | confirmado |
| Higiene do WSL: PATH, apt mínimo, memória e variáveis | Quatro ajustes de configuração, todos medidos nesta máquina, sem instalar nenhuma ferramenta de bio. | 1. No /etc/wsl.conf, acrescentar [interop] appendWindowsPath=false e manter [boot] systemd=true e [user] default=manec.… | WSL 2.7.14.0; Ubuntu 24.04.5; kernel 6.… | P1 | corrigido |
| hyperfine 2.0.0 + ripgrep 15.2.0 + fd 10.5.0 (produtividade no WSL) | hyperfine mede tempo e memória de comandos; ripgrep busca texto; fd busca arquivos. | pixi global install -c conda-forge "hyperfine==2.0.0" "ripgrep==15.2.0" "fd-find==10.5.0" (expõe hyperfine, rg e fd; no… | hyperfine 2.0.0; ripgrep 15.2.0; fd 10.… | P1 | confirmado |
| micromamba 2.9.0 (alternativa avaliada e preterida) | Implementação C++ estática do conda/mamba: ambientes por YAML, activate e `micromamba run`. | Não instalar agora. Se for preciso: binário micromamba-linux-64 com .sha256 da release 2.9.0-0, MAMBA_ROOT_PREFIX em ~/… | 2.9.0 | P2 | confirmado |
| Foldseek 10-941cd33 | Busca estrutural rápida com alfabeto 3Di + aminoácido. Tem modo GPU desde a release 10 (easy-search e easy-multimersearch) e predição de 3D… | Quando houver caso: pixi global install -c conda-forge -c bioconda "foldseek==10.941cd33". Bancos (PDB, AFDB) em ~/; GP… | 10-941cd33 | P2 | confirmado |
| bio-mcp (organização GitHub: bio-mcp-blast, -samtools, -seqkit, -bwa,… | Servidores MCP finos que embrulham CLIs de bioinformática (BLAST, samtools, bcftools, seqkit, bwa, bedtools, fastqc, interpro, amber, evo2,… | Não instalar. Remover o exemplo bio-mcp de docs/STACK.md e chamar as CLIs instaladas pelo pixi diretamente. | 0.1.0 (pyproject; sem release) | P2 | confirmado |

Armadilhas registradas:

- O campo latest_version de api.anaconda.org/package/<canal>/<pkg> erra versões com hash: indicou mmseqs2 17.b804f e foldseek 8.ef4e960, mas o bioconda já publica 18.8cc5c (2025-07-27) e 10.941cd33 (2025-01-20). Ordene por files[].upload_time ou confira a release no GitHub.
- O csvtk saiu do bioconda (receita 404) e vive no conda-forge desde 2024-11. O bioconda parou em 0.31.0 (2024-11-02); o atual é 0.38.0 (2026-09-24). Pinar csvtk=0.31 congela quatro versões atrás.
- O instalador do pixi e o `micromamba shell init` gravam o PATH no ~/.bashrc. Neste WSL, o ~/.bashrc desiste em shell não-interativa, então `wsl -d Ubuntu-24.04 -- bash -l script.sh` não acha o que foi instalado. É a mesma falha parcial da armadilha 3. Use PIXI_NO_PATH_UPDATE=1 e exporte no fim do ~/.profile.
- O hook guard_data bloqueia curl\|sh e bloqueou até o download do install.sh do pixi feito só para leitura. Instale pixi e micromamba pelo tarball da release com verificação do .sha256.
- pixi global não tem lockfile (issue #1280, aberta desde 2024-04-26): o manifesto fixa versão, não build. Para o Métodos, pine `==versão` e registre `pixi global list` (ou ~/.pixi/envs/*/conda-meta/*.json) no provenance. Para reprodutibilidade estrita, use workspace pixi com pixi.lock e detached-environments fora de /mnt/c.
- O pixi não executa scripts post-link por padrão (run-post-link-scripts=false). O post-link.sh do sra-tools no bioconda, que grava /LIBS/GUID em ~/.ncbi/user-settings.mkfg, não roda: execute `vdb-config` uma vez antes do primeiro prefetch.
- Depois da migração haverá dois samtools no PATH: /usr/bin/samtools 1.19.2 (apt) e ~/.pixi/bin/samtools 1.24. Idem para o seqkit 2.13.0 em ~/.local/bin. Quem vence é a ordem do PATH: registre `samtools --version` no provenance e remova o binário solto. CLAUDE.md e WSL_TOOLS em session_context.py precisam refletir a nova origem.
- A htslib 1.19 do apt tem CVEs abertas no decodificador CRAM e no índice GZI (CVE-2026-31962 e CVE-2026-31964 a -31971, corrigidas em 1.23.1 e 1.24). No Ubuntu noble, o CVE-2026-31965 está 'needs-triage'. Não abra CRAM ou BCF de terceiros com a 1.19.
- samtools 1.24 mudou a semente padrão de `view --subsample`: agora vem do cabeçalho, não de 0. O mesmo comando sem `--subsample-seed` dá outra subamostra. Desde a 1.22, o CRAM de saída é 3.1 (htsjdk/GATK podem não ler) e a htslib não busca mais referência CRAM no EBI: passe --reference ou REF_PATH/REF_CACHE.
- O nome `yq` é ambíguo. O `yq` do apt (3.1.0) e do conda-forge (4.4.5) é o wrapper Python do kislyuk, com sintaxe de jq. O yq do mikefarah, que lê e escreve TOML, é o pacote `go-yq` 4.54.1.
- O jq 1.8.0 quebrou compatibilidade: indices/index/rindex passaram a contar por code point, e tonumber rejeita espaços. Não misture o jq 1.7.1 do apt com o 1.8.2 do pixi no mesmo projeto.
- O bwa-mem2 pede cerca de 28×N GB de RAM para indexar (README; N em Gb): ~30 GB para Gallus, que tem 1,07 Gb, e muito mais para crocodilianos. Não cabe nos 15 GB do WSL. Use minimap2 -ax sr ou bwa 0.7.19.
- O DIAMOND não tem modo tblastn, só blastp e blastx. O blastx com o genoma inteiro como consulta única já abortou no Windows ('Masking queries', 1,07 Gb, rc=127). No WSL, fatie o genoma (seqkit sliding/split2), paralelize e converta a coordenada da janela para a do cromossomo antes de juntar.
- O DIAMOND muda defaults dentro da série 2.2. A 2.2.3 trocou o mascaramento do clustering para SEG; a 2.2.6 alterou a sensibilidade do --very-sensitive e passou a exigir BMI/BMI2/LZCNT/F16C/MOVBE no caminho AVX2. Fixe a versão exata no Métodos.
- O miniprot só emite, por padrão, alinhamentos com score de pelo menos 0,99 do melhor (--outs) e até 30 secundários (-N), e não é otimizado para homólogos distantes (README). Uma relíquia divergente pode simplesmente não aparecer. Ausência no miniprot não é ausência: combine com tBLASTn e nhmmer, com controle positivo no mesmo genoma.
- MMseqs2 na GPU: o pré-filtro GPU é sem gaps, sempre na sensibilidade máxima, e pode diferir do pré-filtro k-mer da CPU (issue #1083, mascaramento no banco padded). A clusterização não usa GPU. Declare modo e banco no Métodos. No WSL2, a memória pinned é limitada (guia NVIDIA), o que pode afetar bancos maiores que os 8 GB de VRAM.
- O Docker do WSL não tem nvidia-container-toolkit (só o runtime runc, medido). `--gpus all` e o perfil gpu do Nextflow falham. Até instalar o toolkit, MMseqs2 e Foldseek só usam GPU nativos via pixi. Não instale driver Linux dentro do WSL.
- O PATH do Windows é anexado ao WSL: 90 entradas, 78 em /mnt/c. Vinte buscas de comando inexistente levaram 1557 ms contra 9 ms só com PATH Linux (medido), e npm/npx do Windows resolvem dentro do WSL.
- Sem .wslconfig, o WSL usa 15 GB dos 31,4 GB do host e 4 GB de swap (medido). Índice do miniprot e DIAMOND sobre genoma inteiro podem morrer por OOM (exit 137) se rodarem em paralelo sem --memfree.
- A NCBI_API_KEY foi removida em 2026-09-12, e variáveis do Windows não chegam ao WSL sem WSLENV. Sem chave, o NCBI Datasets limita a 5 req/s (10 com chave). Varreduras em 180 genomas sentem isso.
- Os scripts do atavismo (s23, 01_inventario) chamam /datasets/v2alpha. O NCBI promoveu a v2alpha a v2 em 2024 e retirou a v1; clientes antigos já recebem 410 Gone. O estado atual dos caminhos v2alpha não foi verificado.
- hyperfine 2.0.0 (2026-10-07) roda os comandos sem shell por padrão (pipes e redirecionamentos exigem -S) e mudou o formato do --export-json.
- fastp 1.3.x teve várias correções de travamento ligadas a threads e BGZF (1.3.2 a 1.4.0). Fixe 1.4.0 ou superior se usar fora do nf-core.
- O qsv no conda-forge está parado em 14.0.0 (2026-01-13), enquanto o GitHub está em 24.0.0 (2026-10-05).
- O apt usa nomes próprios: fd-find instala `fdfind`, não `fd`; o pacote `yq` é o do kislyuk.
- Bancos e genomas em /mnt/c atravessam o 9P. Bancos BLAST/DIAMOND/miniprot e FASTA devem ficar em ~/ (ext4); só os resultados finais voltam para runs/.
- O nf-core não tem módulo de nhmmer nem de ncbi datasets. Se a varredura virar pipeline, esses passos precisam de módulo local.
- O pacote 'pixi' do PyPI (1.0.1, de 2020) não é o gerenciador da prefix.dev. Nunca instale via pip.

### Evolução molecular, seleção e perda gênica

| Candidato | O que é | Rota | Versão | Prioridade | Veredito |
|---|---|---|---|---|---|
| HyPhy 2.5.103 (+ hyphy-analyses) | Plataforma de testes de seleção com modelos de códon. Testes embutidos: RELAX (relaxamento ou intensificação, parâmetro K; --srv, --multipl… | WSL/Nextflow: quay.io/biocontainers/hyphy:2.5.103--h526e2cb_0 com -profile docker. Uso interativo: pixi global install … | 2.5.103 (2026-10-06) | P0 | confirmado |
| PAML 4.10.10 (codeml/baseml/MCMCtree) | Máxima verossimilhança com modelos de códon: modelo de ramo (ω por categoria de ramo), ramo-sítio (Test 2), clade model C e modelos de síti… | Windows: paml-4.10.10-win-x86_64.tar.gz das releases do GitHub. WSL/Nextflow: quay.io/biocontainers/paml:4.10.10--h87e0… | 4.10.10 (2026-01-29) | P0 | confirmado |
| MACSE 2.07 (+ pipelines OMM_MACSE e MACSE_ALFIX) | Alinhador de sequências codificantes guiado pela tradução, que admite frameshift e códon de parada com custos próprios. Subprogramas: trimN… | Windows e WSL: java -jar macse_v2.07.jar (sha256 96873a14…6583963 confere com a bioconda). Nextflow: quay.io/biocontain… | MACSE 2.07 (março/2023); OMM_MACSE v12.… | P0 | corrigido |
| miniprot 0.18 | Alinhamento proteína→genoma com splicing e frameshift sobre o genoma inteiro. Gera GFF3/PAF (--gff, --aln) com frameshifts e stops, a parti… | WSL: curl -L https://github.com/lh3/miniprot/releases/download/v0.18/miniprot-0.18_x64-linux.tar.bz2 \| tar -jxf -, ou … | 0.18 (2025-07-11) | P0 | confirmado |
| PhyKIT 2.8.0 | Toolkit de filogenômica em Python. Comandos relevantes: codon_dnds (dN, dS e ω par a par por NG86, LWL85, YN00 ou ML; --stop-policy skip; N… | Windows: uv add phykit==2.8.0. A wheel é py3-none-any e as dependências (biopython, numpy, scipy, scikit-learn 1.9.1, u… | 2.8.0 (PyPI/GitHub, 2026-10-09); biocon… | P0 | confirmado |
| Ensembl Compara 116 — EPO-Extended de 27 saurópsidos (REST) | Alinhamento múltiplo de genoma inteiro (EPO-Extended) de 27 saurópsidos, com o frango em GRCg7b, mais GERP constrained elements de saurópsi… | Nenhuma instalação. GET https://rest.ensembl.org/alignment/region/gallus_gallus/Z:55155480-55157294:1?method=EPO_EXTEND… | Ensembl release 116; frango bGalGal1.ma… | P0 | confirmado |
| UCSC galGal6 — multiz77way, phyloP77way e phastCons77way (API) + lift… | Alinhamento multiz de 77 vertebrados com o frango galGal6 como referência, mais phyloP por base, elementos conservados (phastConsElements77… | Nenhuma instalação. Conservação por região: https://api.genome.ucsc.edu/getData/track?genome=galGal6;track=phyloP77way;… | phyloP77way dataTime 2018-12-10; phastC… | P0 | corrigido |
| IQ-TREE 3.1.4 — modo --ancestral | Reconstrução ancestral marginal por Bayes empírico (--ancestral/-asr, --asr-min; -asr-joint também é aceito) sobre árvore fixa (-te). Gera … | Mesmo binário do pipeline: iqtree 3.1.4 (bioconda/quay) no WSL ou iqtree-3.1.4-Windows.zip. Exemplo: iqtree3 -s aln.fa … | 3.1.4 (2026-09-10) | P0 | confirmado |
| Datação da inativação por ω em ramo misto + lesões compartilhadas (Me… | Protocolo em três passos. (1) Mutações inativadoras compartilhadas em coluna homóloga posicionam a perda num ramo. (2) codeml, modelo de ra… | Sem instalação. https://timetree.org/api/pairwise/9031/8801 (Gallus–Struthio camelus; 8799 é a família Struthionidae e … | TimeTree 5 (API ativa em 2026-10-11) | P0 | corrigido |
| PRANK 251117 | Alinhador consciente da filogenia, que posiciona os gaps pela árvore. Tem modelo de códon (-codon) e exporta ancestrais com indels (-showan… | WSL: quay.io/biocontainers/prank:251117--hd63eeec_0, bioconda prank=251117 (build de 2026-04-09) ou o binário prank.lin… | 251117 (release no GitHub em 2026-03-19… | P1 | confirmado |
| CESAR 2.0 (bioconda cesar 1.02; GitHub v1.03) | Realinhador de éxons ou genes inteiros por HMM, que respeita o quadro de leitura e os sítios de splice. O modo gene detecta deleção complet… | WSL: quay.io/biocontainers/cesar:1.02--hde5307d_0 ou bioconda cesar; o código-fonte v1.03 está em github.com/hillerlab/… | v1.03 (GitHub, 2026-07-18); bioconda 1.… | P1 | confirmado |
| TOGA2 2.0.10 + portal de dados pré-computados | Anotação por projeção de genoma que classifica ortologia e perda gênica (com as mutações inativadoras listadas) e marca retrogenes e pseudo… | Dados: https://genome.senckenberg.de/download/TOGA2/ (pastas TOGA2/, GenomeAlignments/ e Genomes/, mais assemblies_and_… | TOGA2 v2.0.10 (2026-09-10); portal com … | P1 | confirmado |
| PseudoChecker (servidor web) | Plataforma online de inferência de inativação gênica: compara a região genômica da espécie-alvo com o CDS e a estrutura de uma referência, … | Nenhuma: https://pseudochecker.ciimar.up.pt/pseudochecker/ (só formulário). Registrar parâmetros e resultados à mão em … | serviço ativo em 2026-10-11; site © 202… | P2 | corrigido |
| GRASP 2.0.0 (bnkit.jar) | Reconstrução ancestral por máxima verossimilhança com grafos de ordem parcial. A versão 2.0.0 infere indels ancestrais por programação inte… | java -jar bnkit.jar da release GRASP_v2.0.0 (Java 21 presente). Usar só em alinhamento proteico, com -Xmx dentro do lim… | GRASP_v2.0.0 (2026-08-14) | P2 | corrigido |
| PHAST 1.9.9 (phyloFit/phyloP/phastCons) + PhyloAcc 2.5.0 | PHAST: ajuste de modelo neutro (phyloFit), teste de conservação ou aceleração por base e por elemento (phyloP --mode CONACC com --subtree o… | WSL: quay.io/biocontainers/phast:1.9.9--h3ffd07f_0 e quay.io/biocontainers/phyloacc:2.5.0--py313h0aeee62_0, ou pixi -c … | PHAST 1.9.9 (2026-07-25); PhyloAcc 2.5.… | P1 | confirmado |
| R: phytools 2.5-2 + corHMM 2.8 (com ape 5.8-1, castor 1.8.7 e geiger … | Métodos comparativos: modelos Mk, mapeamento estocástico, ASR de caráter discreto, teste de evolução correlacionada (fitPagel) e modelos de… | Windows: install.packages(c('ape','phytools','corHMM','castor','geiger')). Existem binários CRAN para R 4.6, então não … | phytools 2.5-2; corHMM 2.8; ape 5.8-1; … | P1 | confirmado |
| CAPHEINE (veg) — módulos Nextflow de HyPhy para reaproveitar | Pipeline Nextflow DSL2 do grupo do HyPhy, construído no template nf-core. Tem módulos locais hyphy/relax, hyphy/contrastfel, hyphy/labeltre… | Copiar modules/local/hyphy/{relax,contrastfel,labeltree} de github.com/veg/CAPHEINE (exige Nextflow ≥25.10; o WSL tem 2… | sem release; último commit em 2026-06-15 | P1 | confirmado |
| CAFE5 5.1.0 | Modelo de nascimento-morte para o tamanho de famílias gênicas, com variação de taxa entre famílias (gamma) e entre linhagens. | Nextflow: módulo nf-core cafe (quay.io/biocontainers/cafe:5.1.0--h43eeafb_0). Na bioconda só existe build linux-64. | 5.1.0 (release de 2023-04-17) | P1 | confirmado |
| AleRax 1.4.0 (sucessor do GeneRax 2.1.3) | Reconciliação probabilística entre árvore gênica e árvore de espécies, com duplicação, transferência e perda, coestimando as árvores gênica… | WSL: quay.io/biocontainers/alerax:1.4.0--heeb98d5_0 ou bioconda alerax (linux-64/osx-64). GeneRax: quay.io/biocontainer… | 1.4.0 (2025-11-19) | P2 | confirmado |
| RERconverge 0.3.0 | Calcula taxas evolutivas relativas por gene e as associa a um fenótipo convergente em escala genômica, com significância por permulações. | WSL: quay.io/biocontainers/rerconverge:0.3.0--r44h7b50bb2_3 ou bioconda rerconverge. No Windows exige compilar Rcpp/Rcp… | 0.3.0 (release de 2021; repositório ati… | P2 | confirmado |

Armadilhas registradas:

- HyPhy recusa stop interno: o io.CheckAssertion exige alinhamento sem stop e com comprimento múltiplo de 3. O codeml também falha com stop. Exportar do MACSE com -codonForInternalStop NNN -codonForInternalFS NNN, mas guardar a versão sem máscara: mascarar apaga a própria lesão que é a evidência de nível 3/4.
- O MUSCLE, previsto no item 1 do PROJETO-ATAVISMO, não é consciente de códon nem modela frameshift; use MACSE com -seq_lr. O MACSE_ALFIX, feito para conjuntos com muitos frameshifts (pseudogenes), só existe como .sif (Singularity). O Docker do WSL não roda .sif: precisaria instalar Apptainer.
- dS satura em Aves × Crocodylia (243,3 Ma pelo TimeTree). Antes de interpretar ω, inspecionar o dS por ramo (JSON do HyPhy, free-ratio do codeml) e por par (PhyKIT codon_dnds --method YN00 --verbose). Para encurtar ramos longos, adicionar táxons intermediários: Gavialis, Crocodylus, paleognatas e várias tartarugas.
- Falso-positivo de seleção vem de três fontes conhecidas: erro de alinhamento (Fletcher & Yang 2010 doi:10.1093/molbev/msq115; Jordan & Goldman 2012 doi:10.1093/molbev/msr272; Markova-Raina & Petrov 2011 doi:10.1101/gr.115949.110), mutações multinucleotídicas (Venkat et al. 2018 doi:10.1038/s41559-018-0584-5) e variação de taxa sinônima (Wisotsky et al. 2020 doi:10.1093/molbev/msaa037). Rodar BUST…
- RELAX em pseudogene: o ramo de transição, onde a perda ocorreu, precisa de rótulo próprio, como a categoria 'misto' de Meredith 2009. Misturá-lo aos ramos de pseudogene completo dilui K, e ramos sem rótulo viram 'unclassified'.
- Os módulos nf-core hyphy/* fixam HyPhy 2.5.101, mas a 2.5.102/2.5.103 trouxe correções no BUSTED-PH e no GARD: sobrescreva o container para 2.5.103--h526e2cb_0. O nf-core também não tem hyphy/relax nem contrastfel.
- Os módulos nf-core miniprot/index e miniprot/align fixam miniprot 0.11; a versão atual é 0.18.
- O miniprot não é otimizado para homólogos distantes e tem uma heurística contra pseudogenes processados. Um negativo a ~245 Ma de divergência é evidência fraca (nível 2): use várias sondas (jacaré, gavial, Anolis, tartaruga, humano) e cruze com DIAMOND/tblastn rodando no WSL.
- DSPP: o domínio DPP é repetitivo e de baixa complexidade, vem mascarado e mal montado nos bancos (McKnight & Fisher 2009). Buscar pela região DSP N-terminal e pela sintenia com SPP1/DMP1/IBSP, não pelo DPP.
- Os modelos RefSeq de aves marcados 'LOW QUALITY PROTEIN' são corrigidos pelo Gnomon. No XP_046793377.1 de Gallus, por exemplo, a nota diz 'modified relative to its source genomic sequence ... deleted 1 base in 1 codon'. Um pseudogene pode aparecer como proteína íntegra: sempre reconferir na sequência genômica.
- Ausência em montagem antiga não é perda. Sire et al. 2008 declarou AMBN e ENAM 'desaparecidos' no galGal3, e Al-Hashimi et al. 2010 depois achou ψ-ENAM no frango. Regiões GC-ricas faltam em montagens de leitura curta (Hron 2015; Botero-Castro 2017). Reconferir em GRCg7b ou no T2T (Huang et al. 2023) e checar completude com BUSCO/compleasm.
- Ensembl EPO e UCSC multiz são alinhamentos de genoma inteiro, não conscientes de códon: um indel 'compartilhado' ali pode ser artefato. Antes de declarar lesão compartilhada (nível 4), reextrair os segmentos e realinhar com MACSE. Conferir também que o GRCg7b do Ensembl (GCA_016699485.1) bate, base a base, com o GCF_016699485.2 usado nas coordenadas do projeto.
- Em Psittaciformes o cariótipo é rearranjado. No EPO, o segmento ortólogo a ψ-ENAM, que está no chrZ do frango, cai no cromossomo '3' do Strigops. Diferença de cromossomo não refuta ortologia (item 4 em aberto): decidir pela sintenia local.
- galGal6 (GRCg6a, Red Jungle Fowl) não é o mesmo genoma que o GRCg7b (broiler). O liftOver pode falhar exatamente nas relíquias. Além disso, as outras aves do 77-way são montagens de leitura curta de 2014.
- O IQ-TREE --ancestral trata gap como dado ausente e pode preencher estado em posições que o ancestral talvez nem tivesse. A ASR também tende ao consenso (Regra 4). Guardar o .state com a posterior por sítio, testar ancestrais alternativos (Eick et al., MBE, doi:10.1093/molbev/msw223) e usar GRASP ou FastML para os indels.
- TreeShrink e PhyKIT spurious_sequence removem ramos longos, mas em conjuntos com relíquias o ramo longo do pseudogene é sinal, não erro: não aplicar nesses conjuntos.
- ClipKIT com --codon --remove_stop_codons all mascara os stops. É útil para gerar a entrada do HyPhy, mas apaga a evidência; manter uma cópia sem máscara.
- PAML 4.10.9 mudou o formato do arquivo de controle do MCMCtree: BDparas agora exige o flag C ou M, e .ctl antigos falham.
- TOGA2 exige Python ≥3.13 (o WSL tem 3.12.3), Cargo e Nextflow, e o contêiner oficial é Apptainer, não Docker. Downloads do portal anteriores a dez/2025 tinham chains incompletas para algumas espécies, depois recomputadas.
- O campo latest_version do api.anaconda.org ordena versões como texto: para o PRANK ele devolve v.170427, embora o build 251117 de 2026-04-09 exista. Scripts de verificação de versão devem ordenar por data de upload.
- Zoonomia 241-way, 447-way e 470-way são só de mamíferos, com hg38 como referência: servem ao controle de toolkit em Mammalia, não às coordenadas de Gallus. Selectome (v7) e OrthoMaM (só mamíferos) estão desatualizados ou fora do clado.
- godon: a última atualização no Bitbucket é de 2021-09-13 e não há pacote bioconda nem PyPI. Não adotar; ficar com HyPhy e PAML.
- Uma varredura genômica (Projeto 1) testa milhares de genes: corrigir por FDR (skill rigor-estatistico) e confirmar os acertos de aBSREL/BUSTED-E com o codeml ramo-sítio Test 2, depois de checar o dS.
- Datamonkey envia os dados para um servidor de terceiros, o que só é aceitável com dados públicos. A interface não lista BUSTED-E nem BUSTED-PH e não declara a versão do HyPhy.
- Licenças: CAFE5 tem licença própria (IU Open Source License) e build bioconda só para linux-64; GeneRax é AGPL-3.0; o repositório bnkit (GRASP) não tem LICENSE declarado no GitHub.

### Filogenética e filogenômica

| Candidato | O que é | Rota | Versão | Prioridade | Veredito |
|---|---|---|---|---|---|
| IQ-TREE 3 | Inferência por máxima verossimilhança. Inclui ModelFinder e PartitionFinder (com mAIC no 3.1.4), MixtureFinder (-m MIX+MFP), MAST (+T), PoM… | No pipeline: container quay.io/biocontainers/iqtree:3.1.4--h8471819_0 (o bioconda expõe iqtree3 e o link iqtree). Para … | 3.1.4 (GitHub 2026-09-10; bioconda 2026… | P0 | confirmado |
| ASTER (wASTRAL, ASTRAL-IV, ASTRAL-Pro3, CASTER) | Família de métodos de árvore de espécies baseados em quartetos: - wASTRAL: árvores de genes com suporte e comprimento de ramo, ponderadas, … | Pipeline: quay.io/biocontainers/aster:1.25--hc52dbad_1 (build da tag v1.25). Windows: a pasta exe/ do zip do branch Win… | v1.25 (GitHub 2026-04-20). O bioconda e… | P0 | corrigido |
| MAFFT 7.52x | Alinhador múltiplo. L-INS-i (--localpair --maxiterate 1000) é 'provavelmente o mais preciso; recomendado para menos de 200 sequências' segu… | Pipeline: manter o nf-core mafft/align com o container mulled oficial (mafft 7.520 + pigz 2.8). Se precisar da 7.525, g… | Upstream 7.526 (2024-04); bioconda/quay… | P0 | corrigido |
| MACSE v2 (+ PAL2NAL) | O MACSE alinha sequências codificantes em espaço de códon e tolera frameshifts e códons de parada, o que o torna a solução automática para … | MACSE: o módulo nf-core macse/refinealignment (quay.io/biocontainers/macse:2.07--hdfd78af_0), mais um módulo local com … | MACSE 2.07 (mar/2023; bioconda 2.07). P… | P0 | corrigido |
| ClipKIT 2.x + trimAl 1.5 (limpeza de alinhamento) | ClipKIT: 15 modos de aparar (smart-gap, que é o padrão; gappy; kpic*; entropy; gappyout; block-gappy; composition-bias; heterotachy), --cod… | Pipeline: módulos nf-core clipkit e trimal, com patch para quay.io/biocontainers/clipkit:2.14.0--pyhdfd78af_0 e trimal:… | ClipKIT 2.14.0 (2026-08-06; bioconda e … | P0 | corrigido |
| PhyKIT 2.8.0 | Kit de controle de qualidade filogenômico, em CLI e Python. Cobre saturação (saturation), DVMC, sítios parcimônia-informativos (pis), treen… | No Windows: uv add phykit==2.8.0 (wheel py3-none-any; exige Python ≥3.10, e o venv do projeto é 3.11.17). No pipeline/W… | 2.8.0 (PyPI e GitHub 2026-10-09); bioco… | P0 | confirmado |
| gotree 0.5.2 + goalign 0.4.1 | Binários Go únicos para manipular árvores e alinhamentos. O gotree reenraíza por grupo externo com vários táxons (--strict), colapsa ramos … | No Windows: baixar gotree_v0.5.2_amd64.exe e goalign_v0.4.1_amd64.exe dos releases para C:\Users\Manec\tools e pôr no P… | gotree v0.5.2 e goalign v0.4.1 (GitHub … | P0 | confirmado |
| Pilha R 4.6: ape, phangorn, phytools, TreeDist, castor + ggtree/treei… | Ecossistema R de filogenia: - ape: leitura e manipulação de árvores. - phytools e castor: métodos comparativos e reconstrução de caráter. -… | No subdiretório r/: Rscript -e "install.packages('renv'); renv::init()". Depois: install.packages(c('ape','phangorn','p… | ape 5.8-1, phangorn 2.12.1, phytools 2.… | P0 | corrigido |
| OrthoFinder 3 | Inferência filogenética de ortologia: ortogrupos, árvores de gene, ortólogos de cópia única e árvore de espécies enraizada. A versão 3 é vo… | No pipeline: módulo nf-core orthofinder (fixa 3.1.3). Troque a imagem por quay.io/biocontainers/orthofinder:3.1.5--hdfd… | 3.1.5 (GitHub 2026-05-13; bioconda 2026… | P0 | confirmado |
| MUSCLE 5.3 (ensembles) | Alinhador PPP (-align) e Super5 (-super5) com ensembles de alinhamentos: -stratified gera 16 MSAs, -disperse mede a dispersão, addconfseq/l… | No Windows: muscle-win64.v5.3.exe do release. No pipeline: quay.io/biocontainers/muscle:5.3--h9948957_3. O módulo nf-co… | 5.3 (GitHub 2024-11-11; bioconda 5.3) | P1 | confirmado |
| FAMSA2 | Alinhador progressivo de proteína com árvore-guia medoide, feito para famílias enormes. O README cita 3 milhões de transportadores ABC alin… | No pipeline: módulo nf-core famsa/align, que fixa 2.4.1 (quay.io/biocontainers/famsa:2.4.1--h9ee0642_0); a série 2.5.x … | 2.5.2 (GitHub 2026-03-23); bioconda e c… | P1 | confirmado |
| RAxML-NG 2.0.3 | Máxima verossimilhança com busca adaptativa à dificuldade do conjunto (Pythia) e modo --fast com parada antecipada. Faz seleção automática … | No pipeline/WSL: quay.io/biocontainers/raxml-ng:2.0.3--h870a6a7_0. O módulo nf-core injeta --seed=42 quando falta. Não … | 2.0.3 (GitHub 2026-09-03; bioconda 2026… | P1 | confirmado |
| compleasm 0.2.10 / BUSCO 6.1.0 com sauropsida_odb12 | Avalia a completude do genoma e extrai genes de cópia única (BUSCOs) alinhando proteínas ao genoma com miniprot. O compleasm aceita as linh… | Pipeline: o módulo nf-core compleasm/run, com patch opcional para quay.io/biocontainers/compleasm:0.2.10--pyhdfd78af_0 … | compleasm 0.2.10 (GitHub 2026-10-07; bi… | P1 | corrigido |
| TreeShrink 1.4.0 + PhylteR 0.9.12 (sequências outlier) | O TreeShrink remove táxons que inflam o diâmetro das árvores de genes, ou seja, ramos longos espúrios. O PhylteR (R) acha sequências outlie… | TreeShrink: Python puro desde a 1.4.0 (sem R) e fora do PyPI. No Windows: uv pip install git+https://github.com/uym2/Tr… | TreeShrink 1.4.0 (GitHub 2026-05-21; bi… | P1 | confirmado |
| IQ2MC: IQ-TREE 3 + MCMCtree (PAML) | Datação bayesiana por verossimilhança aproximada. O IQ-TREE 3 calcula gradiente e Hessiana, inclusive com partições e misturas, e gera o ar… | Passo do IQ-TREE no container 3.1.4: iqtree3 -s super.fa -p best_model.nex -te especies_enraizada.nwk --dating mcmctree… | --dating mcmctree desde o IQ-TREE 3.0.1… | P1 | confirmado |
| SplitsTree App 6.9.7 | Redes de splits (NeighborNet, redes de consenso de árvores de gene), árvores e figuras interativas. Instala por .msi no Windows; os release… | No Windows: SplitsTree-6.9.7.msi do release em husonlab/splitstree6 | 6.9.7 (GitHub 2026-09-29) | P1 | confirmado |
| nf-core/phyloplace (EPA-ng + gappa) | Posicionamento filogenético em quatro passos: 1. Alinha as sequências de consulta ao alinhamento de referência (HMMER, MAFFT ou Clustal Ome… | No WSL: nextflow run nf-core/phyloplace -r 2.3.0 -profile docker --input ... (Nextflow 26.04.6 e Docker já verificados … | pipeline 2.3.0 (2026-09-23); EPA-ng 0.3… | P1 | confirmado |
| PhyloNetworks.jl + SNaQ.jl (Julia) | Redes filogenéticas explícitas. O PhyloNetworks.jl é o núcleo: lê e manipula redes e calcula CFs de quartetos. O SNaQ.jl estima redes de ní… | Julia via juliaup. No Windows: winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore. No WSL: curl -fsSL https://… | PhyloNetworks.jl v1.4.0 (2026-08-23); S… | P2 | confirmado |
| BEAST 2.7 / BEAST X | Inferência bayesiana de árvores datadas, com relógios relaxados e modelos de árvore. O BEAST X (10.x) traz novos modelos e amostradores; o … | BEAST 2: instale BEAST.v2.7.7.Windows.zip (ou Linux.x86.tgz no WSL) e atualize BEAST.base e BEAST.app para 2.7.8 pelo P… | BEAST 2.7.8 (GitHub 2025-06-17; biocond… | P2 | corrigido |
| Bibliotecas Python de árvore: toytree 3, DendroPy 5, TreeSwift, ete4 | Manipulação e desenho de árvores em Python: - toytree: desenho SVG/HTML em notebook. - DendroPy: biblioteca madura de computação filogenéti… | Windows: uv add toytree==3.0.11 dendropy==5.1.0 treeswift==1.1.51. ete4 no WSL: sudo apt install build-essential e depo… | toytree 3.0.11 (2025-11-26); DendroPy 5… | P2 | corrigido |

Armadilhas registradas:

- Nomes de táxon (medido, IQ-TREE 3.1.4): o FASTA é lido só até o primeiro espaço ('Gallus gallus:Z(ENAM)' vira 'Gallus'), e os caracteres : , ( ) [ ] ; ' viram '_' com apenas um WARNING. Dois táxons do mesmo gênero colidem, e gCF e ASTER deixam de casar os nomes. Saneie para [A-Za-z0-9_.-] antes e mantenha uma tabela de mapeamento.
- Modo -S DIR do IQ-TREE 3.1.4 (medido): o -B 1000 roda, mas o loci.treefile sai sem UFBoot e sem .contree/.ufboot. Com --abayes, o rótulo vira '/0.738'. Para ter suporte por loco, rode uma tarefa por loco, como o pipeline já faz.
- wASTRAL v1.25.4.8 (medido): com rótulo 'SH-aLRT/UFBoot', usa o primeiro número sem aviso, ou seja, pondera pelo SH-aLRT. Com rótulo que começa com '/', aborta com std::invalid_argument (stod). O sed do módulo ASTER_WASTRAL, que mantém o UFBoot, funciona com dois campos, mas quebra se --abayes for adicionado e o rótulo virar 'SH-aLRT/aBayes/UFBoot' (medido).
- gotree collapse support (medido) também lê só o primeiro número de rótulos compostos. Com rótulos 'SH-aLRT/UFBoot', colapsaria pelo SH-aLRT. Separe o rótulo antes.
- Threads e reprodutibilidade (medido no build Windows do 3.1.4, example.phy 17×1998, -B 1000 --alrt 1000 --seed 7): -T 1 reproduziu os arquivos byte a byte; -T 2 repetido mudou um UFBoot de 96 para 97; -T 1, 2 e 4 mantiveram topologia e logL, mas SH-aLRT e UFBoot variaram 1 a 2 pontos. O -T fixo do patch atual do pipeline não garante identidade bit a bit; só -T 1 garante. O build Linux não foi tes…
- --symtest-remove-bad com partições abortou no 3.1.4 (assert em SuperAlignment::buildPattern, medido). Use --symtest-only, leia o .symtest.csv e remova as partições você mesmo no arquivo de partição.
- --mem é recusado com modelos particionados ('-mem option does not work with partition models yet', saída 2, medido). O módulo nf-core já omite --mem nesse caso. Acompanhe a linha 'NOTE: X MB RAM is required' do log.
- --ancestral não modela indels (medido): o nó cujo clado inteiro só tem gap recebeu estado com p = 0,875. Nas relíquias com deleção, mascare essas posições. O PR #223 (ASR com gaps) ainda está aberto; o GRASP é a referência para ASR com indels.
- Versões fixadas nos módulos nf-core estão atrás das atuais: iqtree 3.1.3, clipkit 2.11.4, trimal 1.5.0, mafft 7.520, muscle5 5.1, fasttree 2.1.10, famsa 2.4.1 e orthofinder 3.1.3. A 3.1.3 do IQ-TREE aborta em setRootNode quando o outgroup do -o falta numa partição (#203, corrigido no 3.1.4). Em Métodos, reporte a versão registrada no tópico versions, não a 'mais recente'.
- O ClipKIT mudou de comportamento entre versões. Na 2.12.0, -g sem -m passou a usar o modo gappy em vez de smart-gap. Na 2.14.0, --ambiguity_handling missing virou padrão e ambiguidades IUPAC passaram a contar como ausentes. Subir de 2.11.4 para 2.14.0 muda o que é aparado: fixe a versão e passe -m explícito.
- Os modos kpic* removem singletons (autapomorfias) e encurtam os ramos terminais, o que enviesa datação, saturação, dN/dS e TreeShrink. smart-gap e kpic-smart-gap não tiveram benchmark próprio: o artigo de 2020 avaliou gappy, kpic, kpi e combinações.
- Teste AU: com -zb 1000, o 3.1.4 avisa 'At least -zb 10000 for reliable results!' (medido). O tutorial avançado ainda mostra -zb 1000.
- --scfl tinha erro na detecção de gaps com heterogeneidade de taxas até a 3.1.0 (corrigido na 3.1.1). Não compare sCFL calculado em versões diferentes.
- O --root do ASTER aceita um táxon só (issue #57). Com dois grupos externos (Chrysemys, Anolis), enraíze por Anolis ou reenraíze depois com gotree reroot outgroup. O wASTRAL deixa os ramos terminais sem comprimento, e o ape e outras ferramentas podem rejeitar a árvore.
- UFBoot não é bootstrap padrão. O módulo ASTER_WASTRAL passa valores de UFBoot no modo -S (escala de bootstrap de 0 a 100), combinação que o artigo do wASTRAL não avaliou; o artigo recomenda o aBayes do IQ-TREE, que teve a melhor acurácia. No modo -S, o suporte padrão para nó sem rótulo (-d) é 0.
- Checkpoint do IQ-TREE com --runs: retomar falha com assert em Checkpoint::getArray (#216, aberta). No Nextflow, use --redo dentro da tarefa e deixe o cache para o -resume.
- Binários Windows: IQ-TREE 3.1.4, wASTRAL e gotree rodaram nativos. O trimal.exe 1.5.1 do release saiu com código 127 sem mensagem; ele importa libstdc++-6.dll e libgcc_s_seh-1.dll do MinGW. O RAxML-NG não tem binário Windows.
- Memória: o IQ-TREE reportou 'Host ... 31 GB RAM' no Windows, enquanto o WSL tem 15 GB. Uma corrida grande pode usar o binário Windows ou subir o memory= no .wslconfig, mas aí sai do Nextflow e do registro de proveniência.
- Nomes enganosos: o pacote 'aster' do PyPI (0.0.0a1) não é o ASTER filogenético. O r-phytools do conda-forge parou em 1.2_0 (2022), contra 2.5-2 no CRAN. davidemms/OrthoFinder (2.5.5) é o repositório antigo; o atual é OrthoFinder/OrthoFinder.
- Licença AGPL-3.0 (ASTER, RAxML-NG, EPA-ng, GeneRax): o uso acadêmico local é livre, mas expor essas ferramentas como serviço de rede a terceiros, por exemplo num servidor MCP compartilhado, obriga a oferecer o código-fonte.
- Python: o ete4 no PyPI só tem sdist com Cython, então exige compilador no Windows. O piqtree 0.8.3 exige Python ≥3.12, e o venv do projeto é 3.11.17. A skill phylogenetics usa ETE3, que quebra no Python 3.13 (import de cgi).
- O FAMSA2 foi desenhado e avaliado para proteína ('multiple-protein-sequence alignment'). Para nucleotídeo ou CDS, fique com MAFFT ou MACSE.
- Sequências idênticas: o IQ-TREE mantém duas e recoloca as demais no fim com ramo zero (medido). Esses ramos de comprimento zero entram no gCF e no ASTER como resoluções arbitrárias; deduplique por espécie antes.
- BEAGLE com CUDA em Blackwell sm_120 não foi verificado. Trate BEAST e MrBayes como só CPU.

### Lacuna — Atlas funcionais de galinha para o lado regulatório da Regra 3 (ChickenGTEx, atlas de elementos regulatórios de Pan 2023, hubs FAANG)

| Candidato | O que é | Rota | Versão | Prioridade | Veredito |
|---|---|---|---|---|---|
| Ensembl Regulation — Regulatory Build de Gallus gallus em GRCg7b (v11… | Anotação regulatória do Ensembl para galinha: features (promoter, enhancer, open_chromatin_region, CTCF_binding_site), EMARs (Epigeneticall… | Em bash (Git Bash/WSL): curl --remote-name-all https://ftp.ensembl.org/pub/release-116/regulation/gallus_gallus/GRCg7b/… | Ensembl 116, ainda a corrente: REST /in… | P0 | corrigido |
| Ensembl Regulation — peaks e signal por epigenoma, GRCg7b | Picos (bigBed) e sinal (bigWig) por epigenoma, remapeados pelo Ensembl para GRCg7b: 49 tracks ATAC-seq, 57 ChIP-seq de histona (h3k4me1 15,… | Picos: baixar os 120 .bb (~203 MB) de peaks/atac-seq, peaks/dnase-seq, peaks/chip-seq/tf e peaks/chip-seq/histone/{h3k4… | release 116. peaks/: 120 .bb, todos de … | P0 | corrigido |
| UCSC liftOver chain galGal6 → GRCg7b + tabelas chromAlias | Chain de conversao galGal6ToGCF_016699485.2.over.chain.gz (3,5 MB) e a variante chainBridge (3,65 MB, Last-Modified 2026-03-11), mais as ta… | Mantém os dois arquivos citados e acrescenta a chain inversa, necessária para consultar bigWigs galGal6 (ChickenGTEx, P… | chain galGal6ToGCF_016699485.2.over.cha… | P0 | confirmado |
| Pan et al. 2023 — atlas processado de elementos regulatorios de galin… | Deposito processado do atlas de ~1,57 M elementos em 15 estados de cromatina de 23 tecidos adultos. 16 arquivos, 568.137.448 bytes no total… | Listar: curl -s https://api.figshare.com/v2/articles/20032103 \| jq -r '.files[] \| [.name,.download_url] \| @tsv'. Bai… | v3, DOI 10.6084/m9.figshare.20032103.v3… | P0 | corrigido |
| ChickenGTEx — 257 bigWigs epigenomicos pela API de manifesto do IGV | Manifesto JSON dos 257 tracks epigenomicos do ChickenGTEx, servido por GET simples em /chickengtex/api/igv/list/Epigenetic_Marks. Cada entr… | curl -s https://ngdc.cncb.ac.cn/chickengtex/api/igv/list/Epigenetic_Marks > manifesto.json. Converter as regiões GRCg7b… | Sem versão. 257 tracks, todos com Last-… | P1 | corrigido |
| ChickenGTEx — downloads de molQTL e expressao, via api/download/list | molQTL cis em 28 tecidos (eQTL, sQTL, exQTL, lncQTL, 3'aQTL), matrizes de expressao, genotipos, covariaveis, molQTL fine-mapped e sumarios … | curl -s https://ngdc.cncb.ac.cn/chickengtex/api/download/list \| jq -r '.result[] \| [.[4], .[5]] \| @tsv'. Os registro… | chicken_gtex_v0, 28 entradas. Arquivos … | P1 | corrigido |
| FAANG Data Portal — API Elasticsearch (api.faang.org) | API de busca estilo Elasticsearch sobre os metadados do FAANG Data Coordination Centre do EMBL-EBI. Aceita payload de query/aggs por POST. … | curl -s 'https://api.faang.org/data/dataset/_search/?size=0' -H 'Content-Type: application/json' --data-raw '{"query":{… | Sem versionamento publico. Respondeu HT… | P1 | confirmado |
| Rota de leitura de bigWig: ucsc-bigwigsummary + pyBigWig/deeptools/Cr… | ucsc-bigwigsummary 482 (bioconda, linux-64) para media/max de sinal numa janela, aceitando URL remota. Em Python: pyBigWig 0.3.26 (MIT) par… | Instalar só sob demanda, na convenção do manifesto: pixi global install -c conda-forge -c bioconda --environment kent "… | ucsc-bigwigsummary 482 (bioconda linux-… | P2 | corrigido |
| ChromHMM 1.28 / MACS3 3.0.5 — re-segmentacao propria em GRCg7b | ChromHMM 1.28 (bioconda, linux-64/noarch/osx-64) para aprender e chamar estados de cromatina a partir de marcas; MACS3 3.0.5 para chamada d… | pixi global install -c bioconda chromhmm=1.28 macs3=3.0.5 | chromhmm 1.28 (bioconda, subdirs linux-… | descartar | **refutado** |
| Track Hub Registry — verificado e descartado para este projeto | Registro do EMBL-EBI de track hubs de terceiros, com API de busca. Consultado para achar hubs FAANG de galinha. | Nenhuma instalacao. Consulta: curl -s https://www.trackhubregistry.org/api/info/assemblies | API respondeu HTTP 200 em 2026-10-11; i… | descartar | confirmado |

Armadilhas registradas:

- DUAS MOEDAS DE MONTAGEM NO MESMO LOCO. ChickenGTEx declara a anotacao como Gallus_gallus.GRCg6a.102.gtf.gz e a pagina do IGV confirma galGal6; Pan 2023 e galGal6; as relicias do reliquias_coordenadas.tsv estao em GRCg7b. O Ensembl Regulation v116 publica diretorio para AS DUAS (release-116/regulation/gallus_gallus/{GRCg6a,GRCg7b}/) — e so o de GRCg7b tem a pasta annotation/. Pegue o de GRCg7b e n…
- TRES NOMENCLATURAS DE CROMOSSOMO, E A FALHA E SILENCIOSA. O GFF3 do Ensembl diz 'Z'. O chain e o hub do UCSC dizem 'chrZ'. As coordenadas das relicias dizem 'NC_052572.1'. A linha do chromAlias de GRCg7b e literalmente: NC_052572.1 \| CM028522 \| Z \| Z \| CM028522.1 \| Z \| chrZ. Sem converter, bedtools intersect devolve zero intervalos, exit 0, sem mensagem — a mesma classe de armadilha do $var…
- NAO EXISTE CHAIN DE VOLTA GRCg7b→galGal6 PUBLICADO. O diretorio liftOver/ sob o hub de GCF_016699485.2 so tem chains para GCF_002263795.3 e GCF_014805655.1. Para ir de GRCg7b para galGal6 (necessario se quiser consultar ChickenGTEx a partir de uma coordenada de relicia) e preciso inverter o chain com chainSwap ou usar a variante chainBridge. Converter no sentido errado e descobrir isso no meio do…
- A JANELA DAS RELICIAS E O PIOR CASO PARA liftOver. psi-ENAM e psi-AMBN caem num vao NAO ANOTADO do RefSeq entre NRG1 e LOC112530610, em sequencia de baixa complexidade — exatamente o perfil que cai em gap de chain. Teste o liftOver nessa janela especifica ANTES de montar qualquer pipeline em cima dele; se a janela nao mapear, o candidato 1 (Ensembl GRCg7b nativo) deixa de ser conveniencia e passa…
- molQTL.tar.gz DO ChickenGTEx TEM 621 GB, e exQTL nominal 22 GB. O espaco real e o livre em C: (98 GB em 2026-10-12; o `df` do WSL mostra ~930 GB, mas e o teto do disco virtual): o nominal inteiro nao cabe nem de longe. Os arquivos significant_cis_* (104–448 MB) respondem a pergunta deste projeto; o nominal so serve para colocalizacao.
- pyBigWig E pysam NAO TEM WHEEL win_amd64. pyBigWig 0.3.26 publica so manylinux_2_27/2_28 x86_64 e osx. Consequencia que engana: o wheel do CrossMap 0.7.4 E py3-none-any e instala no Windows sem reclamar, mas requires_dist lista pysam>=0.22, bx-python>=0.11 e pyBigWig>=0.3.22 — a instalacao quebra na dependencia, nao no pacote. Qualquer leitura de bigWig em Python e rota WSL.
- bigWig REMOTO DEPENDE DE HTTP RANGE, E ISSO FOI MEDIDO. curl -r 0-63 em ngdc.cncb.ac.cn/chickengtex/static/Epigenetic_marks/BW/ATAC_Bursa_CA.bw devolve HTTP 206, Accept-Ranges: bytes, Content-Range bytes 0-63/131116427. Se um proxy ou espelho engolir o Range, cada consulta baixa 131 MB inteiros e os 257 tracks viram ~30 GB. Verifique o 206 antes de colocar a URL num laco.
- AS LICENCAS DIVERGEM ENTRE AS QUATRO FONTES, E NENHUMA E UNIFORME. figshare Pan 2023: CC BY 4.0 (reuso e redistribuicao livres com citacao). Artigo Sci Adv 2023: licenca AAAS, NAO CC — nao redistribua figura nem texto. NAR 2025: CC BY 4.0. Nat Genet 2025: sem CC, apenas politica de text-and-data-mining da Springer Nature. Dados do ChickenGTEx: SEM licenca declarada, so 'Copyright © the FarmGTEx-C…
- O GFF3 DO ENSEMBL v116 CARREGA IDs DE GENE DA RELEASE 113. O cabecalho do arquivo diz textualmente: '# Gene IDs in promoter attributes are from Ensembl release 113'. Cruzar o campo gene_id dos promoters com anotacao da 116 pode bater em ID aposentado e produzir juncao silenciosamente incompleta.
- OS ENDPOINTS DE BUSCA POR GENE DO ChickenGTEx NAO SAO REST DOCUMENTADO. api/download/list e api/igv/list/Epigenetic_Marks respondem a GET simples e sao confiaveis. Ja api/qtls/search/paged e api/gene/specificity devolveram {"code":"500"} e {"code":"300","message":"No valid data!"} com parametros de query convencionais — sao endpoints DataTables que esperam o payload completo do navegador (draw/co…
- A ROTA api.faang.org EXIGE BARRA FINAL. /data/dataset/_search/?size=1 devolve 200 com 2.962 datasets; /data/dataset?size=1 devolve 404 do Django e / tambem 404. Quem testa a raiz conclui errado que a API nao existe. E o espelho data.faang.org/api/dataset/_search/ responde 200 mas serve um indice congelado de 2021 (_index: 2021-03-22_faang_build_3_dataset) — resultado plausivel e desatualizado, qu…
- NENHUM TECIDO CRANIOFACIAL EM NENHUMA DAS FONTES. Enumerei os 49 tracks ATAC do Ensembl GRCg7b: bone_marrow, cerebellum, heart, hindbrain, hindlimb_muscle, kidney, liver, lung, mesonephros, skin_of_back, small_intestine, testis, thymus, trachea, uterus. Os 23 tecidos de Pan 2023 e do ChickenGTEx sao todos adultos e nenhum e mandibula ou germe dentario. As histonas do Ensembl sao todas m_20_w (adu…

### Lacuna — Rota de FASTQ → matriz de contagem (nf-core/rnaseq, scrnaseq) e evidência de transcrição sobre a lesão

| Candidato | O que é | Rota | Versão | Prioridade | Veredito |
|---|---|---|---|---|---|
| nf-core/rnaseq | Pipeline Nextflow DSL2 de FASTQ a matriz de contagem: FastQC/fq lint, trimming (Trim Galore 2.3.0 ou fastp 1.3.6), remocao opcional de rRNA… | Nada a instalar. Fumaça: `wsl -d Ubuntu-24.04 -- bash -lc 'cd ~/rnaseq && nextflow run nf-core/rnaseq -r 3.27.0 -profil… | 3.27.0 (publicada em 2026-09-23; pina n… | P0 | corrigido |
| subread (featureCounts) | Sumarizacao de leituras alinhadas sobre intervalos genomicos. Aceita GTF/GFF ou SAF (`-F SAF`: cinco colunas GeneID/Chr/Start/End/Strand), … | `pixi global install -c conda-forge -c bioconda --environment subread "subread==2.1.1"` no WSL (o pixi 0.81.0 está pres… | 2.1.1 (bioconda linux-64, build h577a1d… | P0 | corrigido |
| regtools | `regtools junctions extract` le um BAM indexado e emite BED12 de junçoes de splice, com a coluna score = numero de leituras que sustentam a… | Se for instalar: `pixi global install -c conda-forge -c bioconda --environment regtools "regtools==1.0.0"` no WSL. Na 1… | 1.0.0 (release no GitHub em 2023-01-10,… | P2 | corrigido |
| StringTie | Montagem de transcritos a partir de alinhamento ao genoma, guiada por referencia (`-G`) ou de novo, com quantificaçao por FPKM/TPM. A serie… | `pixi global install -c conda-forge -c bioconda --environment stringtie "stringtie==3.0.3"` no WSL. Para a janela: `sam… | 3.0.3 | P0 | confirmado |
| gffcompare | Compara um GTF montado contra uma anotaçao de referencia e classifica cada transcrito com um codigo de classe (`=` identico, `j` novo isofo… | `pixi global install -c conda-forge -c bioconda --environment gffcompare "gffcompare==0.12.10"` no WSL. Para comparar u… | 0.12.10 (GitHub em 2025-07-12; bioconda… | P1 | corrigido |
| GenMap | Calcula mapeabilidade (k,e) de um genoma: para cada posiçao, quantas vezes o k-mero que comeca ali ocorre no genoma com ate `e` erros. Said… | `pixi global install -c conda-forge -c bioconda --environment genmap "genmap==1.3.0"` no WSL. Depois, `genmap index -F … | 1.3.0 (nenhum release desde 2020-06-17;… | P1 | corrigido |
| igv-reports | Gera um relatorio HTML autocontido com igv.js embutido: recebe um arquivo de sitios (BED, VCF, MAF ou TSV), um FASTA de referencia e faixas… | No WSL (nao no Python do Windows): `uv tool install igv-reports==1.17.0` ou `pixi global install -c conda-forge -c bioc… | 1.17.0 (PyPI) / 1.16.0 (bioconda) | P1 | confirmado |
| nf-core/scrnaseq | Pipeline Nextflow DSL2 de FASTQ single-cell a matriz de contagem em AnnData/MTX. Alinhadores: simpleaf (padrao, salmon/piscem + alevin-fry)… | Nada a instalar: `wsl -d Ubuntu-24.04 -- bash -lc 'cd ~/scrna && nextflow run nf-core/scrnaseq -r 4.2.0 -profile docker… | 4.2.0 (publicada em 2026-07-04; pina nf… | P1 | corrigido |
| mosdepth | Calculo rapido de profundidade por base, por janela (`--by <tamanho>`) ou por regiao de um BED (`--by arquivo.bed`), com saida em bed.gz e … | `pixi global install -c conda-forge -c bioconda --environment mosdepth "mosdepth==0.3.14"` no WSL. Em RNA-seq, NUNCA us… | 0.3.14 (GitHub em 2026-04-24; bioconda … | P2 | corrigido |
| nf-core/scnanoseq | Pipeline Nextflow DSL2 para single-cell/single-nuclei 10X em leitura longa Oxford Nanopore Q20+ (R10.4): QC com NanoPlot/NanoComp/ToulligQC… | Quando houver dado 10X com ONT: `nextflow run nf-core/scnanoseq -r 1.3.0 -profile docker -params-file params.yml -c lim… | 1.3.0 (publicada em 2026-06-26; pina nf… | P2 | corrigido |

Armadilhas registradas:

- ARMADILHA PRINCIPAL — Nextflow 26 quebra flag booleana nua. O WSL tem nextflow 26.04.6, e a partir da 26.04 o parser de sintaxe estrita e o padrao (reversivel com NXF_SYNTAX_PARSER=v1), o que removeu a conversao de tipo na linha de comando. Em pipeline preso a nf-schema < 2.7.2, `--save_reference` ou `--remove_ribo_rna` falham com `Value is [string] but should be [boolean]` — e `--save_reference …
- Os labels de recurso do pipeline nao cabem nesta maquina e a falha nao e obvia. `conf/base.config` do rnaseq 3.27.0 define process_high = 12 cpus / 72.GB x task.attempt, process_high_memory = 200.GB, e `conf/modules/prepare_genome.config` fixa HISAT2_BUILD em 200.GB x attempt. Com 15 GB de RAM no WSL, essas tarefas nunca sao agendadas. Obrigatorio passar um `-c` com `process.resourceLimits = [cpu…
- Indice STAR nao cabe para o grupo externo, e cabe mal para o frango. O README do STAR diz que genoma de mamifero exige no minimo 16 GB e idealmente 32 GB; o docs do nf-core/rnaseq diz ~38 GB para GRCh37; o parametro `limitGenomeGenerateRAM` tem default 31.000.000.000 e o modulo nf-core passa `--limitGenomeGenerateRAM task.memory-100MB`. Tamanhos de referencia conferidos hoje na API do NCBI Datase…
- `--aligner hisat2` nao produz contagem — e isso pode ser exatamente o que se quer, se souber. O docs/usage.md do rnaseq 3.27.0 afirma que a quantificacao nao e executada com HISAT2 por falta de opcao adequada (issue #822). Para a matriz de contagem isso e um beco; para o teste da reliquia e a rota barata: HISAT2 entrega o BAM genomico com junçoes, e a contagem sobre os intervalos feitos a mao sai…
- Quantificacao por transcrito responde zero por construcao, nao por biologia. Salmon e kallisto quantificam somente o que esta no FASTA/GTF de entrada; ψ-ENAM e ψ-AMBN estao, segundo o proprio docs/PROJETO-ATAVISMO.md, num vao nao anotado do chrZ. Rodar `--pseudo_aligner salmon` e concluir 'nao transcrito' repete, com outra ferramenta, o erro que o projeto ja registrou ter cometido e corrigido (le…
- Anotacao RefSeq/NCBI e o padrao destes projetos e quebra o QC de biotipo. O atributo de biotipo em GTF do NCBI aparece so nas linhas `gene`, nao nas linhas `exon` que o featureCounts le; a 3.27.0 passou a pular o QC de biotipo com aviso em vez de abortar (PR #1919), mas o certo e passar `--skip_biotype_qc`. E NAO passe `--gencode` com anotacao RefSeq: esse flag troca `--featurecounts_group_type` …
- Numero de contagem muda entre versoes do pipeline, e muda mais onde dói. O release 3.27.0 subiu salmon 1.10.3 -> 2.7.0 e stringtie 2.2.3 -> 3.0.3, e o proprio CHANGELOG avisa: oscilacoes de TPM/FPKM de dois digitos percentuais e inversoes de detectado/nao-detectado em genes de expressao baixa a moderada, com alinhamento STAR e inferencia de tipo de biblioteca inalterados. Uma reliquia transcrita …
- I/O sobre /mnt/c inviabiliza um run de Nextflow. O CLAUDE.md da stack ja registra que `/mnt/c` atravessa a fronteira de sistema de arquivos e e ordens de magnitude mais lento para muitos arquivos pequenos; o `work/` do Nextflow e o pior caso possivel desse padrao. Rode em `~` dentro do WSL e publique so o `--outdir` final para `/mnt/c/Users/Manec/bio-stack/runs/`. Chamada de fora do WSL continua …
- Multi-mapeamento e a explicacao alternativa que derruba o achado se nao for enfrentada de frente. A familia SCPP tem parálogos no proprio genoma do frango (AMTN, ODAM, SCPPQ1 e vizinhos, ja mapeados no projeto), e leitura de parálogo expresso atribuida a janela chrZ produz exatamente o sinal que se busca. featureCounts por padrao descarta multi-mapeadores (sem `-M`), e `-Q 10` e `--primary` apert…
- Realidade de dado publico por clado, conferida ao vivo hoje na API do ENA Portal (POST em /ena/portal/api/count, result=read_run). Crocodylia (tax_tree 1294634): 1799 read_run no total, mas somente 292 com library_strategy=RNA-Seq e 808 com library_source=TRANSCRIPTOMIC; zero single-cell; zero em plataforma ONT/PacBio. Psittaciformes (9223): 10 single-cell, 3 long-read (a contagem total de RNA-Se…
- Plataforma: nada disto vai para o Python do Windows. subread, regtools, mosdepth, genmap, stringtie e gffcompare nao tem build win-64 no bioconda (conferido: linux-64, linux-aarch64, osx-64, osx-arm64). `pysam` 0.24.1 nao publica wheel win_amd64 no PyPI (so manylinux e macOS), o que arrasta igv-reports e htseq para o WSL. `pytximport` 0.13.0 e noarch e puro Python, logo funcionaria no uv do Windo…
- Upstream quieto nao e upstream quebrado, mas precisa constar nos Metodos. regtools: ultimo release 1.0.0 em 2023-01-10, ultimo push 2024-08-21, build do bioconda refeito em 2025-05-03. GenMap: ultimo release 2020-06-17, ultimo push 2024-06-20, build refeito em 2025-02-19. Qualimap 2.3 (dentro do rnaseq) e de 2023-11-21. Divergencia de licenca a resolver antes de citar: o bioconda declara gffcompa…
- `--save_reference` e obrigatorio na primeira execucao e e justamente um dos flags que a armadilha de booleano atinge. Sem ele, o indice STAR/Salmon e reconstruido a cada run — e nesta maquina a construcao e a etapa que quase nao cabe. Passe `save_reference: true` no `-params-file`, nunca na linha de comando. Os indices saem em `<outdir>/genome/index` e devem ser reaproveitados depois via `--star_…

### Lacuna — Anotação de elementos repetitivos e TE na janela da relíquia (RepeatMasker/Dfam, RepeatModeler2, EarlGrey)

| Candidato | O que é | Rota | Versão | Prioridade | Veredito |
|---|---|---|---|---|---|
| UCSC GenArk + API REST da UCSC (trilhas repeatMasker/repeatModeler/si… | Hubs da UCSC para assemblies do NCBI. Trazem RepeatMasker (cada registro carrega a linha .out: SW, %div, %del, %ins, família, classe, coord… | Nada a instalar. - Anotação: GET https://api.genome.ucsc.edu/getData/track?genome=<GCA/GCF>;track=repeatMasker\|repeatM… | hubs lidos em 2026-10-11: rCroNil1.hap2… | P0 | confirmado |
| Dfam 4.0 + FamDB 3.0.0 (famdb.py) + API REST do Dfam | Banco aberto (CC0) de famílias de TE, com consensos, HMMs e alinhamentos-semente. O FamDB 3.0 particiona o banco por componente (curado/não… | O famdb.py 3.0.0 entra como dependência do repeatmasker 4.2.4 no bioconda. - Baixar da URL versionada https://www.dfam.… | Dfam 4.0 (release 2026-05; API 6.0.0; 7… | P0 | confirmado |
| RepeatMasker 4.2.4 (+ RMBlast 2.17, TRF 4.10, calcDivergenceFromAlign… | Anota e mascara TEs por homologia contra uma biblioteca (FamDB ou -lib). '-a' grava o alinhamento de cada cópia; util/calcDivergenceFromAli… | pixi (WSL): novo [envs.repeatmasker] em bootstrap/pixi-global.toml com repeatmasker==4.2.4. Ele puxa famdb 3.0.0, rmbla… | 4.2.4 (bioconda, 2026-06-22); GitHub v4… | P0 | confirmado |
| WindowMasker + SDust/DustMasker (BLAST+ 2.17, já no pixi) e o soft-ma… | Mascara repetições pela frequência de k-mers, sem biblioteca e sem classificar. O NCBI já entrega todo genomic.fna eucariótico com as repet… | Nada a instalar: já estão no env blast e no env minimap2 do pixi. Para C. siamensis e os demais, usar o genomic.fna do … | BLAST+ 2.17.0 (windowmasker 1.0.0), já … | P2 | corrigido |
| RepeatModeler 2.0.9 (-LTRStruct) | Descobre famílias de TE de novo (RepeatScout + RECON + LTRharvest/LTR_retriever, com extensão por RepeatAfterMe). Produz consensos classifi… | pixi (WSL): [envs.repeatmodeler] com repeatmodeler==2.0.9 e ltr_retriever==2.9.0 fixado. O RepModelConfig.pm descreve L… | 2.0.9 (GitHub 2026-06-11; bioconda 2026… | P1 | confirmado |
| hillerlab/softmask v0.0.2 (Nextflow) | Workflow Nextflow do Hiller Lab. Renomeia cabeçalhos, roda RepeatModeler (só quando não há biblioteca), aplica RepeatMasker em blocos de 10… | git clone https://github.com/hillerlab/softmask (branch master, tag v0.0.2) e rodar no WSL: 'nextflow run src/main.nf -… | v0.0.2 (2026-07-21); repositório criado… | descartar | **refutado** |
| ULTRA 1.2.1 (repetições em tandem) + longdust 1.4 (regiões longas de … | ULTRA rotula repetições em tandem, inclusive degeneradas, com um HMM de pontuação estável. Aceita período máximo (-p), FDR (--fdr) e sub-re… | pixi (WSL): [envs.ultra] com ultra==1.2.1 (bioconda linux-64); longdust==1.4 opcional no mesmo env. Não usar o módulo n… | ULTRA 1.2.1 (GitHub 2025-10-08; biocond… | P1 | confirmado |
| Classificador de lesão por TE: mutações BIG_INS/INTRON_GAIN do TOGA2 … | Receita que cruza as mutações inativadoras do TOGA2 com a anotação de repetições da mesma janela. Diz se a lesão é uma inserção de TE, de q… | bedtools 2.31.1 + RepeatMasker (P0), interseção na janela. Ler todas as mutações do TOGA2, inclusive as mascaradas: BIG… | TOGA2 v2.0.10 no GitHub; a imagem toga2… | P1 | corrigido |
| Earl Grey 7.3.1 | Pipeline automático de TE: RepeatModeler2, curadoria iterativa (TEstrainer), RepeatMasker e desfragmentação. Gera GFF com divergência Kimur… | Só via pixi/bioconda: earlgrey==7.3.1 (linux-64), mais as partições do Dfam 4.0 baixadas pelo configure_dfam40.sh gerad… | 7.3.1 (GitHub/bioconda); o Docker Hub p… | P2 | corrigido |
| Curadoria pontual da família implicada: TE-Aid 1.0.0 (alternativas: T… | O TE-Aid monta um painel de diagnóstico de um consenso (cobertura no genoma, autoalinhamento, ORFs, hits em proteínas de TE) para decidir l… | pixi (WSL): te-aid==1.0.0 (bioconda noarch; usa R + BLAST). TEtrimmer==1.7.2 (bioconda noarch; exige python <3.11, pfam… | TE-Aid 1.0.0 (2025-03-05); TEtrimmer 1.… | P2 | confirmado |

Armadilhas registradas:

- O Dfam padrão engana em arcossauros. Sem -uncurated, o RepeatMasker usa só famílias curadas, e o Dfam 4.0 não tem nenhuma atribuída a Archosauria, Crocodylia, Aves ou Gallus: sobram 197 genéricas. O genoma parece pouco repetitivo e um TE dentro do éxon passa despercebido.
- -species muda a biblioteca. Com -uncurated, 'Crocodylus siamensis' alcança 197 famílias, 'Crocodylus' 1.294 (inclui as 1.097 de C. porosus), 'Crocodylia' 6.178 e 'Gallus gallus' 374. Espécie sem biblioteca própria pede gênero ou ordem, ou -lib; foi o que a UCSC fez no rCroNil1 (-species Crocodylus).
- As bibliotecas DR do Dfam vêm de montagens antigas. C. porosus foi modelado pela EBI em 2020 sobre GCA_001723895.1 (contig N50 32,9 kb), a mesma montagem que o PROJETO-ARCHOSAURIA descarta; famílias longas ou jovens podem faltar. As do gavial (Gozashti et al., 2025) estão no estado 'Submitted', sem revisão de editor.
- O .2bit padrão do GenArk nem sempre é RepeatMasker. No gavial é WindowMasker+SDust (32,91%; o rmsk com Dfam 3.6 cobre só 12,81%). No frango também é WindowMasker (21,30%, contra 12,49% do rmsk). Para o make_lastz_chains no gavial, use GCA_030020295.1.repeatModeler.2bit (45,1%) e registre no provenance qual máscara entrou.
- A anotação rmsk do frango no GenArk usou RepBase 20181026 (crossmatch -s). Não dá para reproduzi-la localmente sem a assinatura RepBase do GIRI, e uma rodada local só com Dfam vai sair diferente. Compare as duas na janela e cite a origem.
- No FamDB 3.0, '-uncurated' depende de dfam40.uncurated.consensus.1.h5: a partição cobre Eukaryota inteiro (4,13 GB gz), sem recorte por táxon. O HMM de Archelosauria (uh.69, 2,28 GB gz) só serve ao engine nhmmer. Baixe pela URL versionada Dfam_4.0: o link 'current' já mudou de release e quebrou o instalador do Earl Grey (corrigido na 7.2.6).
- O pacote bioconda 'dfam' parou na 3.7 (2023); a 4.0 só existe em dfam.org. Partições baixadas para dentro do prefixo de um env pixi podem sumir quando o env é recriado: guarde-as fora e ligue por symlink.
- O RepeatModeler 2.0.9 espera LTR_retriever 2.9.0 (RepModelConfig.pm), mas o bioconda declara '>=2.9' e o mais novo é 3.0.5. Fixe ltr_retriever==2.9.0, senão o -LTRStruct pode quebrar.
- RAM: o RepeatModeler pede >32 GB e o Earl Grey ≥3 GB por thread; o WSL tem 15 GB (12 GB livres, medido). Uma rodada de novo num genoma de 2,1 a 2,3 Gb nesta máquina arrisca OOM e semanas de execução. Prefira GenArk e Dfam, e peça à UCSC o hub de C. siamensis.
- Máscara e busca de relíquia não se misturam. RepeatMasker, TRF, WindowMasker e SDust mascaram o arranjo DPP (Ser-Ser-Asp) do DSPP como repetição simples. Máscara serve para semear o alinhamento (LASTZ); a busca da relíquia roda em sequência sem máscara, ou com -nolow.
- O soft-mask bloqueia a semente do LASTZ. Um TE inserido no éxon fica em minúsculas, a cadeia quebra ali e o TOGA2 pode chamar o éxon de missing/deleted em vez de BIG_INS. A prova da inserção vem da anotação direta da janela sem máscara (±100 kb), não da ausência no TOGA2.
- O TOGA2 redime inserções que ganham sítios de splicing (INTRON_GAIN, via CESAR2+SpliceAI): um TE no éxon flanqueado por GT...AG pode ser lido como íntron novo, e não como lesão. Inserções menores que 50 pb entram como FS_INS, não como BIG_INS (BIG_INDEL_SIZE=50).
- Uma inserção de TE num pseudogene pode ser posterior à perda: sem restrição, as inserções se acumulam. Ela só conta como lesão inativadora se for ortóloga e compartilhada por todas as linhagens que perderam o gene, sem mutação compartilhada mais antiga. A divergência Kimura-CpG até o consenso dá um limite superior de idade; convertê-la em tempo exige a taxa neutra de cada linhagem.
- As bibliotecas de novo têm muito Unknown: 68% no gavial e 52% no frango (GenArk). Sem família atribuída, não há como datar pela paisagem da família.
- Colisão de nomes no nf-core: o módulo 'ultra' é o uLTRA (alinhador de RNA long-read), não o ULTRA de repetições em tandem. Os módulos repeatmasker e repeatmodeler usam containers 4.1.5 e 2.0.5, anteriores ao FamDB 3.0: use -lib ou troque o container.
- O hillerlab/softmask vem configurado para cluster: pede 24 GB × tentativa e 16 CPUs, com 8 retentativas, use_container=false (conda) e imagem ':latest'. Sem override, o executor local recusa a tarefa. Fixe a imagem por digest e confira se os cabeçalhos renomeados voltam aos accessions (CM.../NC_...).
- O repeatmasker.org estava em manutenção em 2026-10-11. A v4.2.5 (2026-09-28, com correções no ProcessRepeats que juntam fragmentos de LINE e DNA) existe só no GitHub. Registre a 4.2.4 e atualize quando chegar ao bioconda. O TRF é AGPL-3.0: declare na seção de Métodos.

### Lacuna — Fonte de calibração fóssil para a datação (Fossil Calibration Database, PBDB, conjunto de Benton 2015) e a sintaxe de prior do MCMCtree

| Candidato | O que é | Rota | Versão | Prioridade | Veredito |
|---|---|---|---|---|---|
| Fossil Calibration Database (dump CC0 no Figshare) | Base curada da série 'Fossil Calibrations' da Palaeontologia Electronica: Benton 2015, Ksepka & Clarke 2015, Joyce 2013, Head 2015/2016, Sm… | curl -L https://ndownloader.figshare.com/files/69543690 -o data/raw/fcd/FossilCalibration_2024-10-07_14-17-26.sql.zip; … | Figshare v1 (2026-10-01, CC0); dump Mar… | P0 | corrigido |
| Paleobiology Database — Data Service API data1.2 | API REST/JSON sem chave. Endpoints úteis: - taxa/single e taxa/list com show=app: primeiro e último aparecimento (fea/fla, tei). - occs/lis… | Nenhuma: httpx ou curl. Exemplos: https://paleobiodb.org/data1.2/occs/list.json?base_name=Vegavis&show=coords,strat,ref… | data1.2 (respondendo em 2026-10-12 UTC) | P0 | confirmado |
| Benton et al. 2015 — Constraints on the timescale of animal evolution… | Calibrações de nós de Metazoa, cada uma com mínimo (fóssil de grupo-coroa, horizonte e idade) e máximo suave justificados. Para Sauropsida:… | Só leitura: https://palaeo-electronica.org/content/pdfs/424.pdf para as justificativas. Os valores estruturados vêm do … | Publicado em fev/2015 (idades no GTS201… | P0 | corrigido |
| Davies et al. 2026 — calibrações de nós de Crocodylia (Proc R Soc B 2… | Seis calibrações justificadas por Brochu, Donoghue e coautores (docx 2.1_node_calibrations e node_prior.nex do MrBayes): - Coroa de Crocody… | Download de 1,1 MB: https://ndownloader.figshare.com/files/67122326 para data/raw/crocs_davies2026/ (sha256 98e433adf0d… | artigo publicado em 2026-10-07; dataset… | P0 | confirmado |
| MCMCtree (PAML 4.10.10) — sintaxe de calibração e checagem do prior e… | As calibrações vão no arquivo de árvore como rótulos entre aspas simples. A árvore precisa ser enraizada e não pode ter comprimentos de ram… | WSL (bash -l, script em arquivo): iqtree3 -s aln -m MODELO -te arvore_calibrada.nwk --dating mcmctree --prefix X; no X.… | PAML 4.10.10 (2026-01-29; bioconda h87e… | P0 | corrigido |
| mcmc3r (R, dosreislab) | Pacote do grupo de dos Reis e Álvarez-Carretero para o MCMCtree: - densidade, CDF e quantis exatos das calibrações do programa (dB/pB/qB, d… | Windows R 4.6.1: install.packages(c('remotes','ape','corpcor','coda')); remotes::install_github('dosreislab/mcmc3r@f6a8… | 0.6.0 (DESCRIPTION 2025-04-17; HEAD f6a… | P1 | corrigido |
| tbea (R, CRAN) | Pré e pós-processamento de análises bayesianas: - c_truncauchy: estima o c da L(tL,p,c,pL) do MCMCtree para o percentil escolhido (0,975 po… | install.packages('tbea') no R 4.6.1 do Windows (binário r-release 1.8.0) ou no WSL. | 1.8.0 (CRAN, 2026-09-06) | P1 | confirmado |
| CladeDate (R, GitHub) | Gera densidades empíricas de calibração a partir das idades das ocorrências fósseis do clado. Faz Monte Carlo sobre o intervalo mínimo–máxi… | No R: install.packages(c('remotes','phangorn','sn')); remotes::install_github('evolucionario/CladeDate@<commit>'). É R … | 1.2 no DESCRIPTION (o README diz 1.3); … | P1 | confirmado |
| Stiller et al. 2024 (B10K, Nature 629:851) — 34 calibrações ST de ave… | Datação da árvore de famílias de aves (363 espécies). Configuração: MCMCtree 4.9h, verossimilhança aproximada, GTR+G, relógio independente … | Baixar só os arquivos pequenos: curl -L 'https://erda.ku.dk/archives/341f72708302f1d0c461ad616e783b86/B10K/data_upload/… | Nature 629:851–860 (2024); arquivo ERDA… | P1 | corrigido |
| palaeoverse (R, CRAN) | Utilitários de paleobiologia: look_up e time_bins convertem nomes de intervalo em idades pela GTS2020 (padrão), pela GTS2012 ou por escalas… | install.packages('palaeoverse') no Windows (binário 1.5.0; depende de sf e h3jsr). | 1.5.0 (CRAN, 2026-06-11) | descartar | **refutado** |

Armadilhas registradas:

- fossilcalibrations.org devolve HTTP 200, mas é falso positivo. O FCD anunciou o fim da vida em 2024, em 2024-10-01 o domínio estava estacionado e hoje é um blog japonês de forex (WordPress). Não citar, não raspar, não linkar nos Métodos: use doi:10.6084/m9.figshare.34046970 e os artigos da série da PE. A REST NESCent/fcdb-api e os links do Fossilworks também estão mortos.
- O dump do FCD tem 48 'Private Draft' que nunca foram publicados no site. Exemplos: Joyce 2013 Pleurodira 111–165,2 e Cryptodira 124–200,2; Ksepka & Phillips 2015 Casuariiformes e Tinamiformes. Filtre PublicationStatus=4 ou confira no artigo. CalibrationQuality vale 'Current' em todas as 270, então o campo não informa nada.
- No dump, os campos de idade por fóssil (Link_CalibrationFossil.MinAge/MaxAge/AgeType) estão quase todos em 0 ou no padrão. Use calibrations.MinAge/MaxAge e os textos de justificativa; MaxAge=0 significa 'sem máximo'. A tabela L_timescales só lista GSA 1999/2009, então a escala real de cada limite está no texto (Benton 2015 usa GTS2012).
- As definições de nó do FCD estão presas a topologias e a uma taxonomia antigas. 'Reptilia' (ID 235) é Testudines × Sauria (NCBI:32561), o que é incompatível com Archelosauria; na topologia molecular, o nó equivalente é 'Diapsida' (ID 236, Lepidosauria×Archosauria), com os mesmos 255,9–295,9. 'Neoaves' (ID 242) é Columbea×Passerea (Jarvis 2014). O taxid 8493, que o FCD rotula 'Crocodylia', hoje é …
- Archelosauria (Testudines×Archosauria) não tem calibração no FCD nem em Benton 2015. Não invente um valor: deixe o nó sob o prior BD entre a raiz e Archosauria até achar fonte.
- Escala de tempo: Benton 2015 e o FCD usam GTS2012, a PBDB tem outra tabela e o ICS vigente é v2024/12. Três exemplos de divergência: - Máximo de Aves (base do Santoniano + incerteza): 86,8 (GTS2012), 85,7 (PBDB), 85,9 (ICS v2024/12: 85,7±0,2). - Limite Olenekiano/Anisiano, base do mínimo de Archosauria: 247,1, 247,0 e 246,7. - Base do Wuchiapingiano, máximo de Archosauria: 259,8±0,4 → 260,2; 259,…
- Vegavis iaai foi revisto depois de 2015. Torres et al. 2025 (Nature 638:146) dataram o fóssil em 69,2–68,4 Ma e reforçaram sua posição de anseriforme dentro das aves-coroa. Benton 2015 usa 66,0 e Ksepka & Clarke 2015 usam 66,5 para o mesmo fóssil. Asteriornis (66,8–66,7 Ma; Field et al. 2020) fica perto do ancestral de Galloanserae. É preciso decidir e justificar o mínimo de Neognathae: 66,0, 66,…
- Crocodylia: não use Portugalosuchus azenhae (Cenomaniano, ~95 Ma) como mínimo. Com dados moleculares ele cai fora da coroa (Darlim et al. 2022), e Davies et al. 2026 o usam só como máximo suave.
- Primeiro aparecimento na PBDB não é calibração de grupo-coroa. Medido hoje: 'Crocodylia' 237 Ma, Crocodylidae no Kimmeridgiano, Gavialidae no Bathoniano, Neornithes e Neognathae em 132,6 Ma, Testudines no Carniano, Psittaciformes em 61,66 Ma. Os intervalos também são grossos: Vegavis aparece como 72,2–66 na PBDB e 69,2–68,4 na publicação; Leidyosuchus como 83,6–72,2 contra ≥75,017 por U–Pb.
- A 'L(tL)' padrão (p=0,1, c=1) tem cauda pesadíssima: o 97,5% cai em 24,43×tL, e 'L(0.535)' leva o 97,5% a 1.307 Ma. Ajuste c (tbea::c_truncauchy) ou garanta um máximo herdado do ancestral, e sempre confira o prior efetivo.
- Não repita o mesmo B em nós aninhados (Benton 2015 dá 66–86,8 a Neornithes e a Neognathae). A regra ancestral mais velho que descendente desloca as marginais. Ponha o mínimo do fóssil no nó a que ele pertence e o máximo no nó mais antigo que ele restringe. Depois rode usedata=0 e compare as marginais (mcmc3r) antes de qualquer posterior.
- Unidade de tempo de 100 Ma: divida os limites por 100. Em ST e SN, divida location e scale e mantenha shape e df. Em G(α,β), multiplique β por 100. Em L só tL muda, porque p e c são relativos. rgene_gamma depende da unidade; sigma2_gamma só precisa mudar com clock=3.
- O parser aceita 'S2N', mas a tabela do manual escreve 'SN2'. Árvore com comprimentos de ramo é rejeitada; os rótulos vão entre aspas simples; a árvore precisa ser enraizada; mínimo duro é pL=1e-300, nunca 0.
- IQ2MC: o .ctl gerado usa hessianfile e ckpfile, que só existem no fork iqtree/paml (tag 4.10.8.iq2mc, 2025-05-22). O PAML 4.10.10 da camada pixi não reconhece essas linhas: remova-as e renomeie .mcmctree.hessian para in.BV. Use sempre o .rooted.nwk gerado, porque a Hessiana depende da ordem dos táxons, e confira se as calibrações sobreviveram nesse arquivo ou reinsira-as.
- O seed padrão do mcmctree.ctl é -1. Fixe a semente (regra da stack), rode pelo menos 2 cadeias independentes e confira convergência e ESS antes de calcular T. Calcule T amostra a amostra a partir do mcmc.txt, não pelas médias.
- MCMCtreeR foi removido do CRAN em 2026-05-08 (e-mail do mantenedor inválido; GitHub parado desde 2019): não adotar. O rootBBB (idade de clado por ponte browniana bayesiana; dsilvestro) não tem licença declarada.
- LSD2 no IQ-TREE (--date): datas antigas entram negativas (-66 = 66 Ma). A doc mostra faixas 'min:max' só para pontas e data pontual para nó interno. A sintaxe 'mrca(A,B) b(min,max)' é do LSD2 standalone. É uma checagem, não substitui o MCMCtree.
