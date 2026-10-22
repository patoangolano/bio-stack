# Modelos de IA para biologia — o que há, o que falta, o que é redundante

Auditoria de 2026-10-22. Critério: **nicho funcional**, não nome de modelo. Um
modelo só entra na lista de lacunas se nenhum outro na stack cobrir o que ele
faz. Modelo melhor no mesmo nicho de um que já está instalado é troca, não
ganho, e está na seção de redundâncias com o motivo.

## Os dois tetos que decidem tudo

| Recurso | Teto | O que limita de verdade |
|---|---|---|
| RTX 5070 Laptop | 8151 MiB, sm_120 | VRAM |
| Colab Pro, A100 80 GB | 81 920 MiB | **unidades de computação, não VRAM** |

O A100 de 80 GB cobre qualquer modelo desta página. A restrição passou a ser
orçamento: a 7,52 CU/h, as 100 CU/mês do Colab Pro dão **~13 h de A100 80 GB**.
A L4 (22,5 GB) a 1,71 CU/h dá ~58 h, e a T4 (15 GB usáveis) a 1,19 CU/h dá
~84 h. Colab não garante GPU específica — Pro+ ainda pode receber T4.

Consequência operacional, e é a regra mais importante aqui: **nunca gaste hora
de A100 no que o NIM hospedado faz de graça.** Dobramento de monômero e de
complexo sai no `openfold3` e no `boltz2` sem consumir CU. A100 é para o que
nada mais cobre — ver `MACE-OFF` abaixo.

## NIM da NVIDIA — superfície medida, não documentada

`health.api.nvidia.com` **não expõe `/v1/models`** (404). Por isso o servidor
em `servers/nvidia_bio_mcp/` crava os slugs. A superfície real foi medida por
sondagem: 404 = rota inexistente, 422 = rota viva, 410 = aposentado.

**Vivos (11).** `boltz2` (complexo + afinidade IC50), `openfold2` (precisa de
`msa_search` antes), `openfold3` (complexos biomoleculares), `msa_search`,
`proteinmpnn` (dobramento inverso), `rfdiffusion` (backbone de novo),
`evo2` + `evo2_generate` (LM de DNA), `genmol` (fragmento→molécula),
`molmim` (otimização de propriedade) e `diffdock`.

O `molmim` foi confirmado com payload real em 2026-10-22: HTTP 200 em 1,4 s,
três análogos de aspirina com QED. O 500 no corpo vazio é só validação, não
indisponibilidade — e a resposta confirma as duas peculiaridades do modelo:
`molecules` volta como JSON **dentro de uma string**, e o campo do SMILES se
chama `sample`, não `smiles`. O `diffdock` segue como único desconhecido do
catálogo: precisa de receptor real para sair do 500, o que se faz por
`mcp__nvidia-bio__dock_ligand` com um `.pdb` em disco.

**Mortos (3), sem slug novo.** `esmfold`, `alphafold2`, `alphafold2_multimer`.
Responderam 410 em 2026-10-21 e saíram do catálogo — `nim_list_models` devolve
11, não 14.

**Sondados e ausentes (11).** `esm2-650m`, `esm2/embedding`, `geneformer`,
`proteina`, `la-proteina`, `moflow`, `vista3d`, `maisi`, `rfdiffusion2`,
`openfold2` sob slug curto, `chai1`. Todos 404, com `boltz2` como controle em
422 na mesma rodada. A chave não destrava nada além dos 11 — não há camada
oculta a descobrir.

No `integrate.api.nvidia.com` a chave lista 83 modelos, praticamente todos LLM
de propósito geral. O único de domínio é `writer/palmyra-med-70b`.

## Lacunas que valem preencher

Ordenadas por quanto abrem de trabalho novo.

### 1. Campo de força por aprendizado de máquina — cobertura zero

A stack tem GROMACS 2026.3 e OpenMM 8.6.1 (campo de força clássico) e o BioEmu
1.4.1 (ensemble de equilíbrio, só backbone). **Nada em acurácia quântica.**

`MACE-OFF` fecha isso e tem integração OpenMM nativa — e o OpenMM já está no
pixi. Erro médio de barreira no TorsionNet-500: **0,25 kcal/mol**, contra 0,5
do AIMNet2, 1,25 do GFN2-xTB e 2,0 do campo clássico SAGE. Reproduz a
superfície de energia livre do tripeptídeo de alanina em água explícita com os
mesmos quatro mínimos do AMBER14SB/TIP3P. Vazão de 2,1 × 10⁶ passos/dia no
A100 80 GB.

Alternativas no mesmo nicho: `AIMNet2` (14 elementos, neutros e carregados,
termo de longo alcance físico) e `LiTEN-FF` (Nature Communications 2026,
~2,7× mais rápido que o MACE-OFF23(M) num teste de 6 869 átomos).

**É o melhor destino para as horas de A100.** Nenhuma outra coisa na stack
chega perto, e o custo por hora se justifica porque o resultado é dinâmica que
o campo clássico não consegue produzir.

### 2. LigandMPNN — o `proteinmpnn` é cego a ligante

O `proteinmpnn`, hospedado e local, ignora ligante, cofator e metal ao
desenhar sequência. O `LigandMPNN` condiciona neles. 2,62 M de parâmetros
contra 1,66 M, 0,9 s por 100 resíduos em CPU. Não disputa VRAM com nada.

Desenhar o entorno de um sítio ativo com cofator sem isso é desenhar no
escuro. Lacuna pequena de instalar e grande de consequência.

### 3. Docking clássico — só há o aprendido

`diffdock` é docking cego aprendido. Docking físico (`gnina`, `smina`, Vina)
falha de maneira **diferente**, e é justamente por isso que serve: rodar os
dois e comparar pose é como se pega alucinação de pose do modelo aprendido.
CPU, barato.

### 4. Função enzimática e homologia remota — cobertura zero

Consulta ao UniProt resolve sequência anotada. Para sequência sem anotação não
há nada: `CLEAN` (número EC) e `ProstT5` (homologia remota por estrutura
inferida) preenchem. Nicho inteiro descoberto.

### 5. AlphaMissense — sobreposição parcial, ganho modesto

Patogenicidade de missense. O `evo2` opera em DNA e o ESM-2 dá
verossimilhança zero-shot; os dois são aproximações de LM. O AlphaMissense é
feito para a tarefa. Entra por último porque o ganho sobre o que já há é
incremental, não categórico.

## Redundâncias — não instalar, e o motivo

### ESMFold

Roda no A100 80 GB (pico de ~41 GiB a 1000 resíduos), mas **é redundante**. O
Protenix-v1 é da família AlphaFold3 e o supera, e `openfold3` e `boltz2` estão
hospedados e de graça. A única vantagem do ESMFold é velocidade por dispensar
MSA. O 410 não é perda que mereça substituição — é um modelo que a stack já
superou por outro caminho.

Para quem ainda quiser: `scripts/colab/esmfold_colab.py` tem as duas alavancas
de memória (`model.esm.half()` e `trunk.set_chunk_size(64)`). O teste do NIH
mediu que `--cpu-offload` **não** salva as sequências que falham no default.

### Borzoi e Enformer

A stack já tem **AlphaGenome** por `mcp__biomcp__alphagenome_predictor`, e ele
ganha do Borzoi: **+25,5%** em sinal de eQTL, +14,7% em log-fold change de
expressão por gene, e recupera **41% dos eQTL do GTEx contra 19%** no limiar
de 90% de acurácia de sinal. Entrada de 1 Mb com resolução de base única —
5× o comprimento e 128× a resolução do Enformer, mais embeddings 2D para
contato de cromatina e splicing.

Ressalva honesta: um estudo estratificado por ancestralidade deu AUROC
ligeiramente maior ao Borzoi (0,837 contra 0,820), e há fatores regulatórios
em que o Borzoi leva vantagem modesta. Não justifica segunda instalação.

### scGPT e Geneformer

Os benchmarks são duros com eles. Um estudo intitula-se literalmente *"one PCA
still rules them all"*; outro mediu o Geneformer **pior que o baseline** em
identificação de célula; em predição de perturbação, baselines aditivos e
lineares simples igualam ou superam os modelos de fundação. Num benchmark da
Genome Biology sobre 22 tecidos, scCello/UCE/scVI/scFoundation/scGPT/Geneformer
venceram em 9/6/3/2/1/1 tecidos — nenhum domina.

`scanpy` + `pydeseq2` não é o elo fraco da stack. O `scVI` é defensável (3 de
22) e já existe o skill `scvi-tools`, mas o ganho não é categórico.

### Chai-1

404 na chave, e sobrepõe `openfold3` e `boltz2`.

## Onde cada coisa deve rodar

| Nicho | Rota |
|---|---|
| dobramento de monômero | `openfold3` hospedado; Protenix-v1 local até 500 tok (6,1 GB) |
| complexo + afinidade | `boltz2` hospedado |
| desenho de sequência | `proteinmpnn`; **LigandMPNN** local quando houver ligante |
| backbone de novo | `rfdiffusion` → encadeia `proteinmpnn` |
| MD clássica | GROMACS, OpenMM (pixi do WSL) |
| **MD em acurácia quântica** | **MACE-OFF no A100 80 GB** |
| ensemble de equilíbrio | BioEmu (`envs/bioemu`, backbone só) |
| LM de proteína | ESM-2 650M local; `esm3_sm_open_v1` 1,4 B (`envs/esm3`) |
| LM de DNA | `evo2` hospedado; DNABERT-2 local |
| efeito de variante regulatória | AlphaGenome via `mcp__biomcp` |
| docking | `diffdock` aprendido **+ gnina físico para conferir** |

## Fontes

Colab: [burn rates medidos em 03/2026](http://mccormickml.com/2024/04/23/colab-gpus-features-and-pricing/).
MACE-OFF: [JACS](https://pubs.acs.org/jacsat/article/147/21/17598/3736678/MACE-OFF-Short-Range-Transferable-Machine-Learning) ·
[arXiv](https://arxiv.org/abs/2312.15211) ·
[LiTEN-FF, Nat Commun 2026](https://www.nature.com/articles/s41467-026-70377-4).
LigandMPNN: [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.12.22.573103v1.full).
AlphaGenome: [Nature](https://nature.com/articles/s41586-025-10014-0) ·
[ressalva por ancestralidade](https://pmc.ncbi.nlm.nih.gov/articles/PMC13320881/).
Single-cell: ["one PCA still rules them all"](https://arxiv.org/pdf/2410.13956) ·
[limitações dos modelos de fundação](https://www.biorxiv.org/content/10.1101/2025.06.26.661767v1.full.pdf) ·
[benchmark Genome Biology](https://link.springer.com/article/10.1186/s13059-025-03781-6).
ESMFold: [ESM](https://github.com/facebookresearch/esm) ·
[medição no NIH Biowulf](https://hpc.nih.gov/apps/esm.html).
Protenix: [benchmark de inferência](https://github.com/bytedance/Protenix/blob/main/docs/training_inference_instructions.md).
