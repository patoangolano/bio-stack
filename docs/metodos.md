# Métodos — relíquias de genes do esmalte em amniotas sem dente

Gerado em 2026-10-17 pela rotina `/metodos`, a partir de `runs/provenance.jsonl`
(5.645 registros entre 2026-10-12T00:00Z e 2026-10-17T14:18Z: 2.316 comandos de
shell, 1.911 execuções de análise, 539 edições de código, 88 consultas a bases),
dos scripts `scripts/atavismo/s24`–`s36` e das versões lidas no ambiente. Duas
linhas antigas do log não são JSON válido e foram ignoradas na contagem.

## Registro e reprodutibilidade

Todas as operações foram registradas automaticamente, com parâmetros e horário
UTC, em `runs/provenance.jsonl`. Cada etapa corresponde a um script versionado no
repositório; a ordem de execução e o comando de cada um estão no cabeçalho do
próprio script. A aleatoriedade foi controlada por gerador local
(`numpy.random.default_rng`, semente 20261017). Não usamos estado global de
gerador.

## Espécies e montagens

Usamos 18 montagens RefSeq do NCBI, acessadas pela API NCBI Datasets v2alpha e
pelo EDirect (efetch 26.0) entre 2026-10-12 e 2026-10-17 (Tabela M1). O
*Melopsittacus undulatus* (GCF_012275295.1) entrou apenas na contagem de genes
anotados. *Casuarius*, *Aquila* e o próprio *Melopsittacus* não puderam ser
analisados por intervalo, porque NRG1 ou LPL não estão anotados nessas montagens.

**Tabela M1.** Montagens e sequências usadas (o acesso da sequência foi conferido
pelo título do registro).

| Grupo | Espécie | Montagem | Sequência do intervalo NRG1..LPL |
|---|---|---|---|
| Crocodylia (com dente) | *Alligator mississippiensis* | GCF_030867095.1 (rAllMis1) | NC_081826.1 |
| | *Gavialis gangeticus* | GCF_001723915.1 (GavGan_comp1) | NW_017729024.1 |
| | *Crocodylus porosus* | GCF_001723895.1 (CroPor_comp1) | NW_017728950.1 |
| Squamata (com dente) | *Anolis carolinensis* | GCF_035594765.1 (rAnoCar3.1.pri) | NC_085842.1 |
| | *Python bivittatus* | GCF_059951505.1 (Pbiv_v1) | NW_028277112.1 |
| | *Podarcis muralis* | GCF_964188315.1 (rPodMur119.hap1.1) | NC_135671.1 |
| Testudines (sem dente) | *Chrysemys picta* | GCF_011386835.1 (ASM1138683v2) | NC_088796.1 |
| | *Chelonia mydas* | GCF_015237465.2 (rCheMyd1.pri.v2) | NC_051245.2 |
| | *Gopherus evgoodei* | GCF_007399415.2 (rGopEvg1_v1.p) | NC_044327.1 |
| | *Pelodiscus sinensis* | GCF_049634645.1 (ASM4963464v1) | NC_134716.1 |
| Aves (sem dente) | *Gallus gallus* | GCF_016699485.2 (GRCg7b) | NC_052572.1 (chrZ) |
| | *Anas platyrhynchos* | GCF_047663525.1 (IASCAAS_PekinDuck_T2T) | NC_092621.1 |
| | *Strigops habroptila* | GCF_004027225.2 (bStrHab1.2.pri) | NC_044302.2 (chrZ) |
| | *Taeniopygia guttata* | GCF_048771995.1 (bTaeGut7.mat) | NC_133063.1 |
| | *Dromaius novaehollandiae* | GCF_036370855.1 (bDroNov1.hap1) | NC_088132.1 (só NRG1) |
| | *Struthio camelus* | GCF_040807025.1 (bStrCam1.hap1) | NC_090982.1 |
| | *Apteryx rowi* | GCF_003343035.1 (aptRow1) | NW_020447342.1 |
| | *Tinamus guttatus* | GCF_059628265.1 (ASM5962826v1) | NW_028288906.1 (só LPL) |

## Proteínas de consulta

As consultas principais foram as proteínas do jacaré: ENAM (QOX06214.1), AMBN
(XP_019346386.1), AMTN (KYO35092.1), o controle positivo LPL (XP_019346320.1) e
o controle negativo DSPP (XP_059578269.1). Para os genes do painel fora do
cluster usamos AMELX (XP_014463267.1), MMP20 (XP_006272555.1), DSPP e ACP4
(XP_019336727.1), além de todas as proteínas anotadas em ±1,5 Mb de cada alvo no
jacaré, que serviram de competidoras. A segunda consulta usou as proteínas RefSeq
de ENAM, AMBN, AMTN e LPL do *Anolis carolinensis*, a primeira isoforma da
tabela de features de cada gene [ACESSO A CONFIRMAR]. O controle de composição
usou a IBSP do *Anolis* (XP_003221343.2).

## Intervalos sintênicos

O cluster ENAM–AMBN–AMTN fica no chr3 do jacaré entre NRG1 e LPL, na ordem
NRG1 → ENAM → AMBN → AMTN → LPL → PSD3. Para cada espécie localizamos as
âncoras NRG1 e LPL pela anotação RefSeq e extraímos o intervalo entre elas com
60 kb de margem de cada lado (`s26_relicto_comparativo.py`,
`s30_especies_extra.py`). Para AMELX, MMP20, DSPP e ACP4, que ficam em outros
cromossomos do jacaré, cada gene recebeu um intervalo próprio
(`s27_painel_intervalos.py`): as âncoras foram o vizinho nomeado mais próximo de
cada lado no jacaré, exigidas nas espécies sem dente e no *Anolis*, com 60 kb de
margem e teto de 3 Mb. Quando o par de âncoras não estava no mesmo registro,
usamos uma única âncora com janela de 750 kb, e o índice dessas espécies deve
ser lido como limite inferior.

## Busca por tBLASTn

As proteínas de consulta foram buscadas nos intervalos com tBLASTn (BLAST+
2.17.0) com `-evalue 10 -seg yes -max_target_seqs 50` (`s26_alinhar.sh`,
`s27_alinhar.sh`). O limiar de e-value ficou frouxo de propósito: a decisão foi
tomada pelo bitscore, que não depende do tamanho do banco. Os intervalos vão de
170 kb a 1,8 Mb, e o mesmo alinhamento muda de e-value conforme o intervalo.
O piso de ruído foi medido com o controle negativo DSPP do jacaré contra o
intervalo do cluster (22–30 bits).

**Competição local** (`s27_indice.py`). Todas as proteínas do intervalo do
jacaré foram buscadas juntas. Um trecho do alvo foi descartado quando outra
proteína alinhou pelo menos 50% dele com bitscore maior.

**Melhor hit recíproco** (`s27_reciproco.sh`). Cada trecho que sobrou foi
traduzido e buscado com BLASTX contra o proteoma inteiro do jacaré
(GCF_030867095.1, 48.040 proteínas; `-evalue 10 -max_target_seqs 5`). O trecho
só contou se o melhor hit foi a proteína de consulta. Dos 69 trechos, 66
confirmaram; os 3 rejeitados eram parálogos de fora do intervalo.

No cluster do frango (`s24_alinhar.sh`, `s25_verificar_relicto.sh`), as regiões
candidatas passaram ainda por dois testes: número de cópias (dc-megablast da
região contra o próprio intervalo, e ≤ 1e-3, exigindo um único hit) e BLASTX
recíproco contra as cinco proteínas de consulta.

**Índices.** Calculamos dois índices: (i) a soma dos melhores bitscores de ENAM,
AMBN e AMTN dividida pelo melhor bitscore do LPL na mesma espécie (`s26`); e
(ii) a soma dos bitscores dos trechos recíprocos dividida pelo valor do próprio
jacaré (`s27`).

## Nulo por janela e segunda consulta

Para separar sinal de ruído (`s33_indice_nulo.py`, `s33_alinhar.sh`), geramos 20
versões embaralhadas de cada proteína de consulta, com a mesma composição e o
mesmo comprimento e a ordem dos resíduos destruída (`rng.permutation`). Elas
foram buscadas no mesmo intervalo de cada espécie com o mesmo tBLASTn. Para cada
gene e espécie, só o excesso do bitscore observado sobre o percentil 95 do nulo
entrou no índice, que é o índice (i) refeito com esses excessos. Repetimos tudo
com a segunda consulta, as proteínas do *Anolis*. Ela está à mesma distância
filogenética de crocodilianos, tartarugas e aves (Lepidosauria × Archelosauria,
cerca de 281 Ma; Kumar et al. 2022), o que elimina a vantagem dos crocodilianos
como parentes próximos da consulta.

## Controle de composição

A IBSP é uma sialoproteína óssea intrinsecamente desordenada como ENAM e AMBN, e
é funcional nas aves. A IBSP do *Anolis* foi buscada nos intervalos do DSPP do
`s27`, onde a IBSP é âncora, com o mesmo nulo do `s33` (`s34_composicao.py`,
`s34_alinhar.sh`). A queda do sinal da IBSP nas aves estima quanto da queda do
ENAM se explica só pela cegueira do tBLASTn para essa classe de proteína.

## Calibração com pseudogene conhecido

O Ψ-AMEL da galinha (Sire et al. 2008; GenBank EU340348) foi buscado por BLASTN
(e ≤ 1e-5) no intervalo do AMELX do `s27` no GRCg7b, para medir se o tBLASTn
com a proteína de crocodiliano detecta um pseudogene de ave reconhecido
(`s27_calibra_amelx.sh`). Repetimos o tBLASTn do AMELX do jacaré com `-seg no`
para excluir o filtro de baixa complexidade como causa.

## Lesões em quadro (nível 3)

As proteínas do jacaré foram alinhadas a cada intervalo com miniprot 0.18-r281,
que modela íntrons e marca frameshifts e códons de parada, com
`-j 2 -k 4 -M 0 -n 1 -l 4 --outc=0.01 --outs=0.5 -t 8` (`s29_lesoes.sh`).
Frameshifts foram lidos das operações F/G da string `cg` e os stops dos códons
`*taa`, `*tag` e `*tga` da string `cs` (`s29_lesoes.py`). Só contamos
alinhamentos que se sobrepõem a um trecho recíproco do `s27`. Antes de ler as
aves, calibramos a taxa de lesão do alinhador nas seis espécies com dente, cujos
genes são funcionais, incluindo três escamados com identidade ao jacaré
(0,36–0,39 no ENAM) na mesma faixa das relíquias.

## Lesões compartilhadas (nível 4)

A relíquia de ENAM do kakapo, recortada do alinhamento do miniprot, foi usada
como consulta em nucleotídeo (BLASTN `-task dc-megablast`, e ≤ 1e-5) contra os
intervalos das demais aves (`s31_aves_nt.sh`). Os trechos homólogos e o CDS de
três crocodilianos foram alinhados com MAFFT v7.525 (`--localpair --maxiterate
1000`, L-INS-i) e com MUSCLE 5.3 (`-align`), na fita codificante e no quadro do
jacaré (`s31_nivel4.py`). Polarizamos cada lesão pelos crocodilianos: deleção é
gap na ave onde a maioria dos crocodilianos tem base; inserção é o inverso; stop
é códon de parada no quadro do jacaré. Uma lesão só contou quando (i) apareceu
com o mesmo tipo, posição e tamanho em pelo menos duas aves, (ii) apareceu nos
dois alinhadores e (iii) não estava em homopolímero (corrida de 4 ou mais bases
iguais). Datamos cada lesão pelo menor clado que contém todas as aves que a têm,
exigindo que toda ave desse clado que cobre a posição também a tenha. Lesões que
violam essa condição foram marcadas como incongruentes com a árvore. A topologia
usada foi Palaeognathae (*Struthio*, (*Apteryx*, *Dromaius*)) e Neognathae
(*Anas*, (*Strigops*, *Taeniopygia*)).

## ODAPH e ODAM

O ortólogo de ODAPH no jacaré foi identificado por posição entre CDKL2 e THAP6
e por estrutura de éxons (Liang et al. 2021), e conferido com BLASTP contra o
Swiss-Prot humano e o proteoma do jacaré (`s32_odam_odaph.py/.sh`; consultas do
UniProt REST, acessado em 2026-10-17). Buscamos ODAPH e ODAM nos intervalos
sintênicos de todas as espécies com tBLASTn e miniprot (mesmos parâmetros acima)
e, depois, nos genomas inteiros de frango, kakapo e jacaré (`s35_genomas.sh`,
`s35_busca_genoma.sh/.py`; tBLASTn `-evalue 10 -seg yes -max_target_seqs 50`,
decisão por bitscore, recíproco por BLASTX contra o proteoma do jacaré). Cada
locus candidato foi comparado entre as duas aves pelos genes vizinhos em ±400 kb
(`s35_sintenia_odaph.py`). A sensibilidade da busca foi medida com uma régua: a
ODAPH do jacaré contra as ODAPH funcionais de *Anolis* e *Podarcis*
(`s35_regua_odaph.sh`).

## Taxa de recombinação

Usamos o mapa de ligação do frango de Elferink et al. (2010), Additional file 5
(bins de cerca de 500 kb, taxa de macho, porque o Z só recombina em ZZ), obtido
pela API REST do Europe PMC (`s28_recombinacao_z.py`). O intervalo NRG1..LPL foi
localizado em galGal3 pelo refGene do UCSC (API REST, 2026-10-17):
chrZ:52.836.615–53.408.327. A taxa da janela foi a média dos bins ponderada pela
sobreposição, comparada com a distribuição dos 122 bins do Z e dos 1.697 bins
autossômicos.

## Validação por alinhamento de genoma inteiro

Baixamos o MAF do chr4 do multiz100way do hg38 (UCSC; 4,6 GB; Blanchette et al.
2004) e o recortamos com `mafsInRegion` (kent 482) nos 29 éxons codificantes de
ENAM, AMBN e AMTN humanos (RefSeq curado NM_031889.3, NM_016519.6 e NM_212557.4),
o que deu 73 blocos (`s36_multiz.py/.sh`). Para cada gene e montagem, medimos a
fração do CDS humano com base alinhada.

## Tempos de divergência

Os tempos de divergência citados vieram da API do TimeTree 5 (Kumar et al. 2022),
acessada em 2026-10-17: *Alligator* × *Gallus* 243 Ma (IC 240–247);
*Alligator* × *Anolis* 281 Ma (277–285); *Alligator* × *Chrysemys* 257 Ma
(249–265); *Dromaius* × *Gallus* 108 Ma (102–114); *Gallus* × *Strigops* 96 Ma
(90–103).

## Estatística

Não calculamos teste de hipótese. As espécies não são independentes
(compartilham filogenia), e o número de linhagens de perda é dois (aves e
tartarugas). As comparações entre grupos são descritivas: faixas e razões, sem
p-valor.

## Software e ambiente

BLAST+ 2.17.0 (bioconda, build hb02a186_1); miniprot 0.18-r281; MAFFT v7.525;
MUSCLE 5.3; samtools/htslib 1.24; SeqKit 2.14.0; NCBI datasets 18.38.0; EDirect
(efetch) 26.0; kent `mafsInRegion` 482. Ferramentas de linha de comando no WSL2
(Ubuntu 24.04, Python 3.12.3), instaladas por pixi a partir de
`bootstrap/pixi-global.toml`. Análise em Python 3.11.17 (numpy 2.4.6, pandas
2.3.3, scipy 1.17.1, requests 2.34.2, matplotlib 3.11.2, reportlab 5.0.1),
gerenciado por uv. Hardware: notebook com 20 núcleos lógicos e 15 GB de RAM no
WSL2. A GPU (RTX 5070, 8 GB) não foi usada nesta análise.

## Tentativas que mudaram o rumo

- **lncRNAs como relíquias** (`s24`). A hipótese inicial era que os lncRNAs do
  Ensembl no intervalo do frango fossem os restos transcritos de ENAM, AMBN e
  AMTN. As relíquias que o tBLASTn detectou caíram nos vãos entre os lncRNAs,
  em sequência sem anotação. O único lncRNA com similaridade nucleotídica ao
  cluster (72,8% em 158 pb, dentro do AMTN do jacaré) não codifica proteína do
  cluster quando traduzido (BLASTX, melhor e = 1,2).
- **Índice em bits sem nulo** (`s26`). Com seis espécies, o índice separou com e
  sem dente. A crítica externa mostrou que três genes no piso de ruído (cerca de
  25 bits cada) divididos por cerca de 170 bits do LPL dão 0,4 sem relíquia
  nenhuma. Ao ampliar para 18 espécies com o nulo do `s33`, a separação
  desapareceu com a consulta do jacaré e só reapareceu com a do *Anolis*. A
  partir disso as conclusões passaram a se apoiar nas lesões (nível 3), e não no
  índice.
- **Afirmações de independência entre linhagens** (`s26`). A inversão AMBN/AMTN
  entre frango e kakapo foi retirada porque os dois valores (23,5 e 24,6 bits)
  estão no piso de ruído.
- **Ponto cego do tBLASTn.** O Ψ-AMEL conhecido não foi detectado pela proteína
  de crocodiliano (zero HSP; 25,8 bits sem SEG). Por isso adotamos a consulta em
  nucleotídeo a partir de uma relíquia de ave para o nível 4.
- **E-value no genoma inteiro** (`s35`). A 1e-2, só a proteína do próprio jacaré
  se encontrou num banco de 2,3 Gb. A decisão passou a ser por bitscore.
- **Table Browser do UCSC** bloqueou o acesso automatizado (Cloudflare
  Turnstile), e o alinhamento PECAN de amniotas do Ensembl só cobre primatas na
  região. Por isso baixamos o MAF do chr4 inteiro.
- **Frameshift pelo multiz** (`s36`). A contagem de corridas de gap fora de
  múltiplo de 3 reprovou no controle: o ENAM funcional do jacaré teve 54 contra o
  humano. A métrica foi descartada e só a cobertura foi usada.

## Uso de inteligência artificial

O desenho das análises, a escrita dos scripts, a execução e a redação deste
texto foram feitos com assistência do Claude Code (Anthropic, modelo Claude
Opus 5.5), sob direção e revisão do autor. A crítica externa do método (§7 de
`docs/METODO-RELIQUIAS.md`) foi produzida por um agente de IA especializado
(`evo-devo`) e conferida contra os dados. Nenhum dado biológico usado aqui foi
gerado por IA: todas as sequências são montagens públicas do NCBI, do UCSC e do
Ensembl.

## Referências citadas nesta seção

- Blanchette M et al. 2004. *Genome Res* 14(4):708–715. doi:10.1101/gr.1933104
- Camacho C et al. 2009. *BMC Bioinformatics* 10:421. doi:10.1186/1471-2105-10-421
- Edgar RC. 2022. *Nat Commun* 13(1):6968. doi:10.1038/s41467-022-34630-w
- Elferink MG et al. 2010. *BMC Genet* 11:11. doi:10.1186/1471-2156-11-11
- Katoh K, Standley DM. 2013. *Mol Biol Evol* 30(4):772–780. doi:10.1093/molbev/mst010
- Kumar S et al. 2022. *Mol Biol Evol* 39(8):msac174. doi:10.1093/molbev/msac174
- Li H. 2023. *Bioinformatics* 39(1):btad014. doi:10.1093/bioinformatics/btad014
- Liang T et al. 2021. *Sci Rep* 11(1):1132. doi:10.1038/s41598-020-80912-y
- Sire JY, Delgado SC, Girondot M. 2008. *BMC Evol Biol* 8:246. doi:10.1186/1471-2148-8-246

## Pendências

1. **Montagem do gavial.** O docstring do `s26_relicto_comparativo.py` diz
   GCA_030020295.1, mas a execução usou GCF_001723915.1 (GavGan_comp1,
   NW_017729024.1). Corrigir o docstring ou refazer com a montagem nova.
2. **Acessos das proteínas do *Anolis*** usadas como segunda consulta no `s33`:
   o cabeçalho do FASTA não guardou o acesso. Recuperar e citar.
3. **Proteínas de consulta GenBank.** ENAM (QOX06214.1) e AMTN (KYO35092.1) do
   jacaré são GenBank, não RefSeq. Conferir se coincidem com os XP atuais da
   GCF_030867095.1.
4. **Datas de acesso por espécie.** O log permite datar cada download; o texto
   dá o intervalo 2026-10-12 a 2026-10-17.
5. **Versão do multiz100way** (data do arquivo baixado do UCSC) e versão da
   montagem humana usada pelo RefSeq curado.
6. **Coordenadas em galGal3** no `s28`: a relíquia foi mapeada em GRCg7b e a
   taxa vem de galGal3. Não houve liftOver da relíquia; a janela é a do par de
   âncoras.
7. **Assistência de IA**: confirmar se a instituição ou a revista exige
   declaração em formato próprio.
8. **Repositório**: o paper cita github.com/patoangolano/bio-stack como
   disponibilidade de código; confirmar que ele está público (ou trocar por um
   depósito com DOI, como Zenodo) antes de circular o texto.
9. **Afiliação do autor** no cabeçalho do paper (`[AFILIAÇÃO A PREENCHER]`).
