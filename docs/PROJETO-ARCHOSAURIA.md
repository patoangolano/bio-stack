# Projeto 1 — Psittaciformes × Crocodylia: genômica comparativa de Archosauria

Plano de 2026-09-14. Deriva de uma pesquisa de 108 agentes com verificação
adversarial (130 alegações extraídas → 25 testadas → **7 sobreviveram**) mais
medição direta contra a API do NCBI Datasets.

Leia a §5 antes de tudo: ela diz o que a pergunta central **não** pode responder.

---

## 1. A pergunta, reformulada para ser respondível

Você escreveu: *"projetar quais são as mutações e mudanças que podem ter
ocorrido nos últimos 250 milhões de anos que os separam"*.

Como está, a pergunta não tem resposta — nem com computação infinita. O que é
inferível é outra coisa, e ainda assim é bastante:

| Pergunta | Respondível? |
|---|---|
| Quais genes são ortólogos entre periquito e jacaré | **Sim**, diretamente |
| Qual o estado *provável* de um sítio no ancestral comum | **Sim**, com probabilidade posterior por sítio |
| Quantas substituições ocorreram em cada ramo | **Sim**, como contagem inferida |
| Quais famílias gênicas expandiram ou contraíram em cada linhagem | **Sim**, com ressalvas |
| Em que **ordem** as mutações ocorreram dentro de um ramo | **Não** |
| Quais mutações *de fato* ocorreram (a trajetória histórica real) | **Não** |

Os próprios autores do Cactus escrevem que reconstruções de genoma ancestral são
*"only a hypothetical reconstruction of the true history"* — não um registro
histórico. Essa frase deve aparecer na sua Discussão, citada. Ela não enfraquece
o trabalho: delimitá-la é o que separa um TCC defensável de uma alegação
indefensável.

**Reformulação sugerida:** *"Quais diferenças moleculares entre Psittaciformes e
Crocodylia são derivadas em cada linhagem, e em quais genes o padrão de
substituição indica mudança de regime evolutivo desde o ancestral comum de
Archosauria?"*

---

## 2. Espécies e acessões

Medido na API do NCBI em 2026-09-14. **Não copie estes números para o texto sem
reconferir** — release de anotação muda por baixo.

### Conjunto principal

| Papel | Espécie | Accession | Nível | ctg N50 | Anotação |
|---|---|---|---|---|---|
| Psittaciformes | *Melopsittacus undulatus* | `GCF_012275295.1` | cromossomo | 16,6 Mb | RefSeq EGAP 2024-11 |
| Psittaciformes | *Strigops habroptila* | `GCF_004027225.2` | cromossomo | 9,5 Mb | RefSeq EGAP 2020-03 |
| Crocodylia | *Alligator mississippiensis* | `GCF_030867095.1` | cromossomo | 38,6 Mb | RefSeq EGAP 2023-09 |
| Crocodylia | *Gavialis gangeticus* | `GCA_030020295.1` | cromossomo | 84,2 Mb | VGL EGAPx 2026-01 |

Os dois crocodilianos cobrem Alligatoridae e Gavialidae — as duas raízes da
divergência dentro de Crocodylia. *Crocodylus niloticus* `GCA_054825085.1`
(ctg N50 111,5 Mb) fecharia o trio, mas **não tem anotação**: serve para
sintenia, não para análise de gene.

### Grupos externos — não são opcionais

| Papel | Espécie | Accession |
|---|---|---|
| Ave fora de Neoaves | *Gallus gallus* | `GCF_016699485.2` |
| Testudines | *Chrysemys picta* | `GCF_007399415.2` |
| Lepidosauria | *Anolis carolinensis* | `GCF_035594765.1` |

Sem grupo externo você não **polariza**. Se periquito tem A e jacaré tem G, só a
tartaruga e o lagarto dizem qual é o estado ancestral e qual é a mudança
derivada — que é exatamente a sua pergunta. Um trabalho sem grupo externo
descreve diferenças; com grupo externo, descreve **mudanças**.

### Armadilha confirmada: não use o gavial antigo

`GCF_001723915.1` (GavGan_comp1, 2016) anuncia scaffold N50 de 96 Mb. É
artefato: foi montado com **Ragout usando o próprio jacaré como referência**, o
contig N50 real é 14,8 kb, e 26,4% da sequência é padding de N. Compará-lo ao
jacaré é circular por construção. O mesmo vale, em menor grau, para
*Crocodylus porosus* `GCF_001723895.1` (contig N50 32,9 kb).

O `GCA_030020295.1` (rGavGan2.hap2, VGP) resolve isso: contig N50 de 84,2 Mb,
zero gaps, anotado.

### Download

Nenhuma credencial, nenhum CLI, nenhum root. Verificado: `HTTP 200`, ZIP válido.

```
https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_012275295.1/download?include_annotation_type=PROT_FASTA,CDS_FASTA,GENOME_FASTA,GENOME_GFF
```

Baixe para `~/` **dentro do WSL**, nunca para `/mnt/c` — a travessia de sistema
de arquivos é ordens de grandeza mais lenta para muitos arquivos pequenos.

---

## 3. A restrição de hardware que define o desenho

**Progressive Cactus está fora de alcance.** Verificado em três fontes
independentes: a documentação oficial provisiona 208–416 GiB por job de
alinhamento; o pipeline de saurópsidas da Harvard FAS exige cluster SLURM; o
alinhamento de 605 amniotas consumiu centenas de nós AWS de 32 vCPU. Você tem
15 GB. Não é questão de paciência — não cabe.

Pior: o consumo de memória do Cactus **escala com o número de núcleos alocados**
enquanto o tempo satura por volta de 24 núcleos. Jogar os 20 núcleos nele piora
o teto de memória em vez de ajudar.

Caminho viável para alinhamento de genoma inteiro:

- `minimap2 -x asm20` — pareado, rápido, cabe folgado
- LASTZ + chain/net — o padrão clássico do UCSC para pares distantes
- `nf-core/pairgenomealign` — embrulha isso em Nextflow, que você já tem rodando

**Estratégia central:** faça o grosso do trabalho em **espaço de proteína e
códon**, não em genoma inteiro. Os conjuntos de proteínas têm ~16–20 mil
sequências por espécie (dezenas de MB), e ortologia + alinhamento de códon +
teste de seleção cabem confortavelmente em 20 núcleos / 15 GB. O alinhamento de
genoma inteiro fica para a pergunta de sintenia, em cima de um par só.

---

## 4. Cadeia de análise

Ordem de dependência: cada etapa produz o insumo da seguinte.

1. **Baixar** os 7 genomas (proteína + CDS + GFF). ~100 MB de proteína/CDS.
2. **BUSCO** em cada conjunto, mesma linhagem para todos — mede se a anotação é
   comparável entre espécies **antes** de qualquer conclusão depender disso.
3. **Ortologia**: OrthoFinder sobre os conjuntos de proteína. Saída: grupos de
   ortólogos e ortólogos de cópia única.
4. **Alinhamento de códon**: alinhar proteína, retro-traduzir para códon
   (nunca alinhe DNA codificante diretamente). Filtrar colunas mal alinhadas.
5. **Árvore**: IQ-TREE sobre a matriz concatenada dos ortólogos de cópia única.
6. **Teste de seleção** por gene, com modelos de ramo e de ramo-sítio.
7. **Reconstrução ancestral** nos nós internos, guardando a probabilidade
   posterior por sítio — sem ela o resultado não é interpretável.
8. **Estrutura** (opcional): para os genes com sinal mais forte, dobrar a
   proteína atual e a ancestral reconstruída via NIM hospedado e comparar onde
   as substituições caem — superfície? sítio ativo? núcleo hidrofóbico?

O passo 8 é o que a sua stack faz e uma stack comum não faz. É também onde o TCC
fica visualmente forte.

---

## 5. Armadilhas estatísticas

Estas **não** foram verificadas adversarialmente nesta pesquisa — o ângulo de
seleção não produziu alegação sobrevivente. Trate como hipóteses a conferir na
literatura primária antes de escrever Métodos, não como fato estabelecido:

- **Saturação de sítios sinônimos.** A ~250 milhões de anos, dS provavelmente
  satura. Um dN/dS calculado sobre dS saturado não significa o que parece.
  Olhe a distribuição de dS antes de interpretar qualquer dN/dS.
- **Erro de alinhamento gera falso-positivo de seleção.** É o modo de falha mais
  documentado da área. A escolha do alinhador e do filtro muda a lista de genes.
  Rode com dois alinhadores diferentes e reporte a interseção.
- **Teste múltiplo.** Uma varredura genômica testa ~10 mil genes. Sem correção,
  centenas de falsos-positivos. A skill `rigor-estatistico` vale aqui — as regras
  de FDR são as mesmas do RNA-seq.
- **Confundidor de anotação.** Jacaré e periquito usam EGAP; o gavial usa
  **EGAPx**, pipeline diferente. Contagens de gene entre pipelines diferentes não
  são estritamente comparáveis. Para ortologia isso é tolerável; para
  "família X contraiu na linhagem Y", não é — ali o BUSCO da etapa 2 é o
  controle que torna a afirmação sustentável ou não.
- **Viés da reconstrução ancestral.** Há literatura de que ASR tende a produzir
  sequências tipo-consenso. Confira antes de afirmar propriedade da proteína
  ancestral.

---

## 6. Itens em aberto

Registrado como ausência de cobertura, não como resposta negativa.

1. **A data.** "250 milhões de anos" é a cifra comum para a divergência
   ave–crocodiliano, mas **nenhuma alegação sobre datação molecular sobreviveu à
   verificação**. Antes de pôr um número no texto, busque a estimativa e o
   intervalo de confiança em estudo calibrado por fóssil.
2. **Neosuchia × Crocodylia.** Nenhuma alegação sobre nomenclatura de clado
   sobreviveu. Provisoriamente: Neosuchia inclui formas extintas; para genomas
   sequenciados o termo operacional é Crocodylia. Confirme antes de titular.
3. **A pergunta biológica específica.** Vocalização, beta-queratinas, receptores
   olfativos, longevidade — nenhuma foi verificada quanto a precedente
   publicado. Escolher uma com literatura de apoio é o que ancora a comparação e
   dá com o que comparar o seu resultado.
4. **Versões e requisitos das ferramentas** da §4 não foram medidos nesta
   máquina. Próximo passo prático: instalar uma de cada vez e medir, em vez de
   confiar em requisito anunciado.
