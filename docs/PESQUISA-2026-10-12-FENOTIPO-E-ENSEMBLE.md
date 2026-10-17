# Programa de pesquisa — eixo fenotípico e eixo de ensemble

Levantamento de 2026-10-12. Produzido com dois agentes `lit-scout` em paralelo,
busca ranqueada no Consensus, pesquisa profunda do NotebookLM (43 fontes
encontradas, 31 importadas) e **leitura do código do RERconverge 0.3.0 instalado
nesta máquina**.

Framework de decisão: Fischbach & Walsh, *Problem choice and decision trees in
science and engineering*, Cell 187:1828–1833 (2024), via a skill
`bio-research:scientific-problem-selection`.

Caderno NotebookLM com as 31 fontes:
<https://notebook.google.com/notebook/d2264cfa-0a51-48d5-8f08-d92a7792eedf>

> **Leia a §1 antes de qualquer coisa.** Um go/no-go disparou e o resultado foi
> **não** para o Desenho 1 como estava escopado. O desenho sobrevive, mas com um
> parâmetro trocado, e a troca é barata.

---

## 1. O go/no-go que disparou: 7 genomas não bastam

Três medidas independentes, todas obtidas hoje, convergem no mesmo ponto.

**Medida 1 — o piso do próprio pacote.** Lido do RERconverge 0.3.0 instalado em
`r/renv/library/windows/R-4.6/x86_64-w64-mingw32/RERconverge`:

```r
args(getAllCor)
# function (RERmat, charP, method = "auto", min.sp = 10, min.pos = 2, ...)

args(correlateWithBinaryPhenotype)
# function (RERmat, charP, min.sp = 10, min.pos = 2, ...)
```

A documentação instalada define `min.sp` como *"Minimum number of species that
must be present for a gene"*. O padrão é **10**. Com 7 genomas, todo gene é
descartado antes de qualquer teste. Não é recomendação de estilo — é o filtro
que roda.

**Medida 2 — o tamanho de amostra do precedente.** Da meta-análise metodológica
sobre as fontes importadas:

| Estudo | Caráter | Pontas |
|---|---|---|
| Meredith et al. 2014 | perda de dente mineralizado em aves | **48 genomas de aves** + jacaré + grupos externos |
| Kowalczyk et al. 2022 | perda de pelo em mamíferos (RERconverge) | **62 genomas de mamíferos** |
| Springer et al. 2016 | inativação de C4orf26 em placentários sem dente | 18 sequências na análise de seleção (103 na filogenia de suporte) |

O menor número de pontas com resultado publicado numa análise quantitativa
dessas fontes é **18**, e aquela era dN/dS no PAML, não taxa evolutiva relativa.
Os dois estudos que usaram taxa relativa de fato operaram com 48 e 62 pontas.

**Medida 3 — onde o sinal vive.** Taxa evolutiva **relativa** é, por construção,
uma medida de desvio do gene contra o fundo da árvore. Com 7 pontas o fundo é
estimado sobre 7 ramos. O estimador existe, mas a distribuição nula contra a
qual o desvio seria julgado é rasa.

### O que fazer — trocar o parâmetro, não o desenho

Fischbach & Walsh, armadilha 4: *parâmetros fixos em excesso*. O Desenho 1 havia
fixado três coisas de uma vez — o método (RERconverge), o conjunto de táxons (os
7 genomas do Projeto 1) e o caráter (dentição). Fixe **um**.

**Fixe o caráter. Deixe o conjunto de táxons flutuar.**

Aves têm dezenas de genomas anotados publicados — Meredith usou 48 em 2014, e o
número hoje é muito maior. Expandir o conjunto é **download e anotação, não
método novo**. Os 7 genomas do `docs/PROJETO-ARCHOSAURIA.md` continuam sendo o
núcleo da comparação Psittaciformes × Crocodylia; o Desenho 1 usa um conjunto
maior porque a estatística dele exige, e isso não conflita com o Projeto 1 — são
dois recortes do mesmo material.

Antes de baixar nada: medir espaço em `C:` (`df -h /c`) e conferir o orçamento
de anotação. Genoma de ave anotado é da ordem de 1 GB por espécie em
proteína+CDS+GFF; 50 espécies é da ordem de 50 GB, e havia 95 GB livres.

---

## 2. A lacuna confirmada — e ela é o ativo do projeto

**Nenhum estudo nas 31 fontes aplicou RERconverge, ou análise equivalente de
taxa evolutiva relativa, à perda de dentição ou de esmalte.**

Os três trabalhos que existem sobre o tema usaram outra coisa:

- Meredith et al. 2014 — detecção direta de mutação inativadora (*frameshift*,
  códon de parada prematuro, deleção) + taxa de *frameshift* + dN/dS
- Springer et al. 2016 — Ka/Ks com modelos de ramo no CODEML/PAML
- Emerling et al. 2023 (xenartros) — pseudogenização compartilhada + relaxamento
  de seleção no PAML

Todos perguntam *o gene está quebrado?*. Nenhum pergunta *a taxa de evolução
deste gene desviou junto com a perda do caráter?*. São perguntas diferentes: a
primeira é binária e detecta o evento; a segunda é contínua e detecta regime.

Isto responde o **item aberto 3** do `docs/PROJETO-ARCHOSAURIA.md` — "a pergunta
biológica específica… nenhuma foi verificada quanto a precedente publicado". A
pergunta tem precedente de método (Kowalczyk 2022 em pelo, Partha 2017 em visão)
e não tem precedente de aplicação ao caráter. É o lugar certo para estar.

> **Armadilha 2 de Fischbach & Walsh, enfrentada de frente:** "aplicar ferramenta
> existente a sistema novo" é incrementalismo quando a ferramenta já é
> amplamente aplicada ali. Aqui ela não é — a área de perda de dente usa
> detecção de pseudogene, não taxa relativa. A contribuição é metodológica e
> biológica ao mesmo tempo. Mas o teste honesto é o da §4: a resposta precisa
> ser diferente da que a detecção de pseudogene já dá, senão é ornamento.

---

## 3. Dois itens abertos do Projeto 1, fechados com fonte

### Item 1 — a data

`docs/PROJETO-ARCHOSAURIA.md` registrou: *"'250 milhões de anos' é a cifra comum
para a divergência ave–crocodiliano, mas nenhuma alegação sobre datação molecular
sobreviveu à verificação"*. A divergência continua a conferir, mas a data que
este projeto precisa é outra — a da perda do caráter:

**A perda do revestimento de esmalte ocorreu há ~116 milhões de anos, IC 95%
105,7–128,3 Ma, no ancestral comum de Neornithes** (Meredith et al. 2014). Vem
com intervalo de confiança, que é exatamente o que o item pedia.

### O conjunto de genes é SEIS, não quatro

O trabalho que este projeto vinha seguindo rastreia 4 genes. O conjunto canônico
de Meredith et al. 2014 é **6**:

| Gene | Papel | Estado em 48 aves | Estado em *Alligator mississippiensis* |
|---|---|---|---|
| `DSPP` | dentina | inativado | intacto |
| `AMTN` | matriz de esmalte | inativado | intacto |
| `AMBN` | matriz de esmalte | inativado | intacto |
| `ENAM` | matriz de esmalte | inativado | intacto |
| **`AMELX`** | matriz de esmalte | inativado | intacto |
| `MMP20` | matriz de esmalte | inativado | intacto |

`AMELX` e `DSPP` estão fora do conjunto que o `scripts/atavismo/s18_chrZ.py`
mediu. Acrescentar os dois é a primeira coisa a fazer, e é barata.

Nota a conferir, não afirmação: `AMELX` é ligado ao X em mamíferos. Onde ele
cai no genoma de ave, e se cai no chrZ junto com os outros, é pergunta direta
para o inventário — e se cair, reforça a §5.

Springer et al. 2016 acrescenta um sétimo candidato, `C4orf26`, inativado em
placentários sem dente. Estado em aves: **não localizado nas fontes**. Vale
checar.

---

## 4. O regime estatístico, lido do código e confirmado na literatura

### FDR existe, mas não é onde eu disse primeiro

`correlateWithBinaryPhenotype` não contém `p.adjust` no próprio corpo. Ela
delega para `getAllCor`, e é lá que a correção acontece:

```r
corout$p.adj = p.adjust(corout$P, method = "BH")
```

Benjamini–Hochberg, sobre a coluna `P`. Então o `p.adj` existe.

### E o `p.adj` não é o número defensável

O artigo de permulações filogenéticas (**Saputra et al. 2021**, entre as fontes
importadas) estabelece o que torna isso insuficiente, e é preciso:

- p-valores paramétricos em teste de associação filogenética são **mal
  calibrados e anticonservativos**, porque as espécies que compartilham o caráter
  estão agrupadas por ancestralidade e o teste trata cada ponta como ocorrência
  independente;
- **aplicar FDR sobre p-valor descalibrado não conserta**: o FDR ajusta a
  multiplicidade, não a violação do pressuposto de independência que distorceu
  a distribuição nula;
- **permulações** resolvem gerando distribuição nula empírica que preserva a
  estrutura filogenética e a contagem do grupo focal.

Isto é o mesmo erro estrutural que o item 6 do checklist da skill
`rigor-estatistico` descreve em single-cell: *"com 5.000 células de 3
camundongos, o n não é 5.000 — é 3"*. Aqui: com 50 aves num clado, o n não é 50
— é o número de **origens independentes** da perda. E a perda de dente nas aves
teve **uma única origem** (Meredith et al. 2014). Um evento.

> **Consequência que molda o desenho todo.** Se o caráter tem uma origem só, não
> há convergência para detectar, e a potência de qualquer método de taxa
> relativa é mínima — por mais genomas que se empilhe. A força de Kowalczyk 2022
> em pelo e de Partha 2017 em visão vem de **múltiplas perdas independentes**.
>
> Saída, e é ela que torna o desenho viável: **sair de Neornithes**. Perda de
> dente ou de esmalte tem origens independentes em amniotas — xenartros
> (Emerling 2023), baleias misticetas, pangolim, tamanduá, placentários sem dente
> (Springer 2016), além das aves. Aí há réplica evolutiva, e a comparação
> Psittaciformes × Crocodylia entra como **um** dos eventos, com o jacaré de
> controle positivo — os 6 genes intactos nele são o único ponto de ancoragem
> funcional da árvore.

### O pacote tem o maquinário

Verificado nos exports do RERconverge 0.3.0 instalado — 22 funções de
permulação, entre elas:

```
getPermsBinary  getPermsBinaryExtantOnly  getPermsBinaryFudged
getPermsContinuous  permpvalcor  adaptivePermulation
plotPositivesFromPermulations  categoricalPermulations
```

`permpvalcor(realcor, permvals)` devolve, segundo a documentação instalada,
*"a data frame containing permulation p-values"*. É esse o número que vai para
a tabela do artigo. O `p.adj` do BH vai para o material suplementar, nomeado
como paramétrico.

---

## 5. Eixos novos — o que os pumps de intuição produziram

Três eixos que não estavam no plano, cada um com a literatura que o sustenta ou
o ameaça.

### E1. Por que a relíquia ainda é detectável? (*dogs that don't bark*)

116 milhões de anos de decaimento neutro deveriam ter apagado a homologia
reconhecível. As relíquias estão lá. Por quê?

Três explicações concorrentes, e as três são testáveis só com sequência:

**(a) Pleiotropia retarda a perda.** Zaraisky et al. 2024, *Biological Reviews*:
em muitos casos a perda de genes associados à regeneração em amniotas foi
**significativamente atrasada** em relação ao momento em que a capacidade foi
perdida, e a hipótese dos autores é pleiotropia — o gene só cai quando algo
compensa suas funções não relacionadas ao caráter. Perda de *enhancer*
altamente específico, por contraste, pode ser neutra e rápida.
→ **Predição testável:** genes de esmalte com função pleiotrópica documentada
(`MMP20` é metaloproteinase de matriz; `AMTN` tem papel relatado além do esmalte)
decaem **mais lentamente** que os de função exclusiva (`AMBN`, `ENAM`). Ordenar
os 6 genes por grau de decaimento e correlacionar com pleiotropia anotada.

**(b) Deleção de DNA, não pseudogenização, é o modo dominante.** Crombez et al.
2025, *Nature Communications*: simulação de pseudogene neutro mostra que, **na
ausência de deleção de DNA**, pseudogenes deveriam permanecer detectáveis por
muito mais tempo do que se observa. Observam enriquecimento de pseudogene em
regiões de **alta** recombinação, e Ka/Ks não-neutro em alguns pseudogenes.
→ **Predição testável:** se a relíquia persiste por baixa taxa de deleção, ela
deve estar em região de baixa recombinação — e não por seleção.

**(c) O chrZ não é ambiente neutro.** Charlesworth 2020, *Phil Trans R Soc B*: a
degeneração começa quando o crossing-over cessa, com taxa inicialmente rápida e
depois muito mais lenta; as duas quantidades que governam são o tamanho da região
ligada ao sexo e o momento da supressão de recombinação. Wu et al. 2015
(*J Mol Evol*, mamão) mostram que o efeito Hill-Robertson em região
pericentromérica reduz a eficácia seletiva **também** no lado homólogo — posição
cromossômica, não só ligação ao sexo.
→ **Confundidor que precisa de controle explícito.** As relíquias foram medidas
**no chrZ**. Afirmar seleção residual sem separar região recombinante de região
suprimida é afirmar sobre o mecanismo errado.

E1 é o eixo mais barato de todos: **nenhum fenótipo medido, nenhuma segmentação
de crânio, nenhum genoma novo além dos que o Desenho 1 já exige.** É sequência e
coordenada genômica.

### E2. Sinal de seleção sem dN/dS — e o alerta que quase o mata

A §5 do `PROJETO-ARCHOSAURIA.md` avisa que a 250 Myr o dS satura. Spielman &
Wilke 2015 (*MBE*) acrescenta um **segundo** mecanismo de viés, independente da
saturação: se os códons sinônimos diferem em aptidão, dN/dS pode assumir valores
arbitrariamente altos **com toda a seleção sendo purificadora**; os modelos
estilo Goldman-Yang, os mais usados, enviesam substancialmente; e as estimativas
menos enviesadas **nunca** vêm do modelo de melhor ajuste por AIC/BIC. Escolher
modelo por critério de ajuste é armadilha.

Alternativa tentadora: verossimilhança de ESM-2 650M, que cabe nos 8 GB desta
GPU, como medida de restrição evolutiva.

**Pugh et al. 2025 (bioRxiv) quase mata a ideia e entrega o conserto.**
Verossimilhança de sequência reflete não só restrição de aptidão mas também
**estrutura filogenética e viés de amostragem**, e o problema **piora** conforme
a capacidade do modelo aumenta. Usar verossimilhança crua como sinal de seleção
numa comparação **filogenética** é precisamente esse modo de falha. O conserto
proposto é *Likelihood-Fitness Bridging* — média de escores sobre sequências
sujeitas a pressão semelhante, que sob um modelo de **Ornstein-Uhlenbeck**
equivale a marginalizar a deriva. Aplica-se sem retreinar.

Gordon et al. 2024 (bioRxiv) acrescenta a verificação prévia obrigatória: tanto
sequência tipo-selvagem **super**-preferida quanto **sub**-preferida piora o
desempenho, com cauda de lei de potência devida a homologia.
→ **Antes de usar ESM-2 como medida, localizar onde cada proteína de ave e de
crocodiliano cai na distribuição de verossimilhança do modelo.** Se caírem em
regimes diferentes de preferência, a comparação entre linhagens mede o viés do
modelo, não a biologia. Isto é um controle, não uma formalidade.

O modelo OU que o LFB assume já está na stack — `geiger` e `phytools`, em
`r/renv.lock`.

### E3. Matriz sistemática de perda (*make it systematic*)

Sharma et al. 2018, *Nature Communications* (215 citações) é a referência
metodológica: detecção acurada de perda gênica associada a adaptação em
mamíferos aquáticos e voadores, mostrando que perda por relaxamento **após** a
adaptação informa sobre a biologia da espécie.

Lyubetsky et al. 2023, *Biology Direct* (método WINEGRET) é o irmão
metodológico do Desenho 1: busca genes perdidos **no mesmo estágio evolutivo** em
que o caráter foi perdido, usando homologia e sintenia local em vez de base de
ortólogos.
→ **Ressalva que o exclui por ora:** WINEGRET exige RNA-seq profundo do caráter
em desenvolvimento no organismo modelo. Esse dado não existe aqui. O caminho por
taxa evolutiva relativa continua sendo o indicado.

E a revisão que nomeia o campo: **"From trees to traits: a review of advances in
PhyloG2P methods and future directions"**. **`PhyloG2P`** é o termo de busca que
faltava — é o nome da classe de método, e procurar por ele abre a literatura
inteira de uma vez.

---

## 6. Desenho 2 — o que a controvérsia entrega

A questão **está aberta**, com duas posições nomeadas:

| Posição | Quem sustenta |
|---|---|
| A hiperestabilidade é propriedade do **ambiente ancestral** (Terra mais quente) | Gaucher et al. 2003/2008; Perez-Jimenez et al. 2011; Risso et al. 2013/2014; Akanuma et al. 2013; Benner |
| A hiperestabilidade é **artefato de consenso** do método de ASR | Trudeau, Kaltenbach & Tawfik 2016; Wheeler, Lim, Marqusee & Harms 2016/2021; Lim et al. 2016 |

O argumento do lado do artefato: ASR por máxima verossimilhança enviesa para
aminoácido de consenso em posição ambígua, e mutação em direção ao consenso
tende a ser **aditivamente estabilizadora** — o método fabrica superioridade
ancestral que não sempre se reproduz *in vivo*.

**Onde o Desenho 2 entra.** A controvérsia inteira foi travada em
**termoestabilidade medida** — temperatura de desnaturação, atividade residual.
Ninguém, nas fontes, a atacou por **ensemble conformacional**. O ensemble é uma
medida diferente e mais informativa: não "derrete a que temperatura" mas "quanta
heterogeneidade o estado nativo admite".

A predição é limpa e falsificável: **se o viés de consenso for a explicação, a
reconstrução ancestral deve ser mais rígida que TODAS as pontas extantes** — não
mais rígida que a média, mais rígida que o máximo. Rigidez maior que a média é
compatível com ancestral genuinamente estável; rigidez fora do envelope de toda
a diversidade extante é assinatura de construção.

Isso transforma a ressalva da §5 do `PROJETO-ARCHOSAURIA.md` — *"há literatura de
que ASR tende a produzir sequências tipo-consenso; confira antes de afirmar
propriedade da proteína ancestral"* — de citação obrigatória em **teste que se
roda**.

Fonte de métrica localizada: *"Ligand binding remodels protein side-chain
conformational heterogeneity"* (eLife) trata heterogeneidade conformacional
quantitativamente. **Pendente:** qual métrica é aceita em artigo revisado para
afirmar "mais rígido" — e se BioEmu, que prevê **só backbone**, sustenta a
afirmação sem dinâmica explícita, ou se o OpenMM é obrigatório. Os dois
`lit-scout` ainda estavam rodando quando este documento foi escrito.

---

## 7. Matriz de risco

Escala de risco 1–5 (5 = existencial). "Tempo até ler" é tempo até o resultado
que confirma ou derruba a premissa.

| # | Premissa | Tipo | Risco | Tempo até ler | Estado |
|---|---|---|---|---|---|
| 1 | 7 genomas bastam para taxa evolutiva relativa | Técnica | 5 | **já leu** | **DERRUBADA, mas não pelo motivo que eu disse.** Existe precedente com 7 pontas (Gaspar 2026, perda de voo em rails). O que derruba é **colinearidade**: com 7 pontas de Archosauria o caráter é perfeitamente confundido com a bipartição Aves/Crocodylia e o espaço de permulação é vazio. `min.sp=10` é o sintoma; a colinearidade é a causa |
| 2 | Há réplica evolutiva suficiente para taxa relativa | Biológica | 5 → **2** | **já leu** | **RESOLVIDA.** 12 origens independentes de perda de dente ou de esmalte em amniotas (Meredith 2014 + Springer 2016 + Emerling 2023). Condição: o recorte tem de ser amniota, não Neornithes |
| 3 | RERconverge em perda de dente é novo | Técnica | 3 | **já leu** | **CONFIRMADA.** Nenhum estudo nas 31 fontes o fez |
| 4 | A relíquia persiste por seleção, não por ambiente genômico | Biológica | 4 | 2–4 semanas | **ABERTA.** Três explicações concorrentes (§5 E1); confundidor do chrZ é real |
| 5 | Verossimilhança de ESM-2 serve como sinal de seleção filogenético | Técnica | 4 | 2 semanas | **CONDICIONADA.** Exige LFB (Pugh 2025) + verificação de preferência (Gordon 2024) |
| 6 | BioEmu sustenta afirmação de rigidez em artigo | Técnica | 3 | 1 semana | **PENDENTE** de resposta do levantamento |
| 7 | A dentição é mensurável como fenótipo nos táxons de interesse | Técnica | 3 | 1–2 meses | Gargalo é segmentação; TotalSegmentator é treinado em TC clínica humana |

### Passo 0 executado no mesmo dia — e o resultado virou o jogo

A premissa 2 era a que decidia o projeto. Foi lida hoje, não em uma semana.

**São 12 origens evolutivas independentes** de perda de dentição ou de esmalte em
amniotas, documentadas nas fontes. Meredith et al. 2014 abre dizendo que
*"edentulism… has evolved convergently among vertebrates, including birds,
turtles, and several lineages of mammals"*:

| Linhagem | Perdeu | Genomas nas fontes |
|---|---|---|
| Aves (Neornithes) | dente inteiro | 48 |
| **Testudines** | dente inteiro | 3 (*Chrysemys picta*, *Chelonia mydas*, *Pelodiscus sinensis*) |
| Pholidota (pangolim) | dente inteiro | 1 |
| Mysticeti | dente inteiro no adulto | 11 espécies / 3 genomas |
| Ancestral de Pilosa | só esmalte | maioria dos xenartros viventes |
| Vermilingua (tamanduás) | dente inteiro | maioria dos tamanduás viventes |
| Cingulata (tatu) | só esmalte | 1 (*Dasypus novemcinctus*) |
| Tubulidentata (aardvark) | só esmalte | 1 |
| Kogiidae | só esmalte | não quantificado |
| Monodontidae (narval) | só esmalte | não quantificado |
| Tachyglossidae (equidnas) | dente inteiro | não quantificado |
| Sirenia (*Hydrodamalis gigas*) | dente inteiro | não quantificado |

Dois ganhos que não estavam previstos.

**Primeiro: o caráter é graduado, não binário.** "Dente inteiro" e "só esmalte"
são dois estados de perda. Isso permite codificação **categórica**, que o
RERconverge suporta — a função `correlateWithCategoricalPhenotype` e o
`categoricalPermulations` estão nos exports do pacote instalado, e há artigo
específico no caderno (*RERconverge Expansion: Using Relative Evolutionary Rates
to Study Complex Categorical Trait Evolution*). Caráter graduado tem mais
informação que binário, e a §4 do levantamento era escrita supondo binário.

**Segundo, e é o melhor achado do dia: o conjunto de 7 genomas da §2 do
`PROJETO-ARCHOSAURIA.md` já está estruturado certo.** Conferido contra a tabela
de acessões:

| Papel no desenho | Espécies que o projeto já tem |
|---|---|
| Perda, evento 1 (Neornithes) | *Melopsittacus undulatus*, *Strigops habroptila*, *Gallus gallus* |
| **Perda, evento 2 — independente** | ***Chrysemys picta*** (Testudines) |
| **Controle positivo funcional** | *Alligator mississippiensis*, *Gavialis gangeticus* — os 6 genes **intactos** |
| Polarização | *Anolis carolinensis* |

A tartaruga não entrou ali como evento de perda — entrou como grupo externo para
polarizar. Mas pela Meredith 2014 ela **é** uma perda independente. O conjunto
tem, sem nenhum download novo, **duas origens de perda, três linhagens com dente
retido e o melhor controle positivo que existe** para esses 6 genes.

**Consequência para o escopo.** Archosauria deixa de ser o conjunto de dados e
passa a ser a **âncora**: fornece um evento de perda com data e IC, um segundo
evento independente no grupo-irmão, e o controle positivo que nenhum conjunto de
mamífero oferece. Os outros 10 eventos entram como réplica de primeiro plano, e
são todos mamíferos com genoma publicado. Duas origens não bastam para taxa
relativa; doze bastam. E a estrutura não muda — só cresce.

**Bandeiras vermelhas, pelo critério do framework.** Duas premissas em risco 5 é
o limite do que se chama "múltiplos milagres". A premissa 1 já caiu e tem
conserto barato. A premissa **2 é a que decide o projeto**, e ela não tem
conserto por esforço — tem conserto por **mudança de recorte taxonômico**.

**Luz verde.** O go/no-go de maior risco era testável em menos de um dia, e foi
testado hoje. É o oposto do padrão do avestruz que o framework descreve.

---

## 8. Ordem de execução

Reordenada pelo que o levantamento mostrou. A ordem anterior (1 depois 2) não
sobrevive: o Desenho 1 ganhou uma dependência de dado que o Desenho 2 não tem.

**Passo 0 — FEITO em 2026-10-12.** 12 origens independentes, acima do limiar de
~5 que eu havia fixado para seguir. O Desenho 1 tem potência se o recorte for
amniota. Ver o quadro na §7.

**Passo 1 — Desenho 2, que não depende do Passo 0.** Precisa de um gene, não de
um conjunto de táxons. Candidato natural: a proteína de esmalte com melhor
alinhamento e maior profundidade de amostragem entre as 6, reconstruída no nó de
Archosauria. Entrega o teste do viés de ASR, que é publicável por si e entra na
Discussão do Projeto 1 como validação do passo 7 da §4.

**Passo 2 — E1, persistência da relíquia.** Barato, independente do Passo 0, e
dá o controle que o Desenho 1 vai precisar de qualquer forma. Acrescentar `AMELX`
e `DSPP` ao inventário do `s18_chrZ.py`, checar `C4orf26`, mapear recombinação
no chrZ, ordenar os genes por decaimento, correlacionar com pleiotropia anotada.

**Passo 3 — Desenho 1, com o recorte corrigido.** Só depois do Passo 0, e só se
ele devolver réplica evolutiva suficiente. Conjunto de táxons expandido,
permulações obrigatórias, `p.adj` do BH no suplementar.

**Passo 4 — E2, ESM-2 com LFB.** Depende de ter o alinhamento do Passo 3.

**E3 e o eixo de mecânica/morfometria ficam fora deste programa por ora.** Não
porque sejam ruins, mas porque a §7 do catálogo mede o gargalo: segmentação de
crânio não humano é trabalho manual por espécime, e nada acima precisa dele.

---

## 9. Buscas a fazer — o que este levantamento deixou em aberto

Com os termos certos, que é o que a meta-análise entregou de mais útil.

1. **`PhyloG2P`** — o nome da classe de método. A revisão *From trees to traits*
   está no caderno; seguir as referências dela.
2. **Origens independentes de perda de dentição em amniotas**, por evento, com
   genoma publicado. É o Passo 0.
3. **`C4orf26` em aves** — estado não localizado nas fontes.
4. **Posição genômica de `AMELX` em aves** — chrZ ou não.
5. **Métrica de rigidez de ensemble aceita em revisão por pares**, e se BioEmu
   (backbone apenas) basta. Levantamento em curso.
6. **Potência de permulações em função do número de origens independentes** — há
   estudo de simulação? Isso responderia o Passo 0 quantitativamente em vez de
   por contagem.
7. **Mapa de recombinação do chrZ de *Gallus gallus*** — necessário para o
   confundidor da §5 E1.

---

## Procedência

Tudo nesta sessão está em `runs/provenance.jsonl`. As afirmações sobre o
RERconverge são reproduzíveis com:

```bash
cd r && Rscript -e 'suppressMessages(library(RERconverge)); print(args(getAllCor))'
```

As 31 fontes estão no caderno do NotebookLM, consultável por
`mcp__gemini-notebook-mcp__notebook_query` com o id
`d2264cfa-0a51-48d5-8f08-d92a7792eedf`.

Afirmações marcadas **não localizado nas fontes** são ausência de cobertura
deste levantamento, não resposta negativa.

---

# Adendo da mesma sessão — o que os dois levantamentos mudaram

Os dois agentes `lit-scout` voltaram depois que o corpo acima já estava escrito.
Trouxeram duas correções que derrubam premissas minhas e uma medida que vindica
o trabalho que já existe. Fica como adendo, não reescrito por cima, para o
histórico ficar legível.

## A1. A premissa do chrZ — alarme falso que rendeu a melhor medida do dia

O levantamento alertou: **nenhuma fonte publicada coloca relíquia de gene de
esmalte no cromossomo Z de *Gallus***. Por Sire, Delgado & Girondot 2008
(`10.1186/1471-2148-8-246`) a região candidata de AMBN/ENAM na galinha é o
**cromossomo 4**, flanqueada por SULT1E1 e IGJ. Busca dirigida ao chrZ: zero.

O alarme leu o título errado. O `scripts/atavismo/s18_chrZ.py` não afirma que os
genes estão no chrZ — a docstring dele diz, literalmente, *"Âncoras sintênicas
ARCOSSAURIAS (não mamíferas): o cluster de esmalte do jacaré é flanqueado por
LPL/PSD3. Em Gallus, LPL está no cromossomo Z."* É afirmação sobre **vizinhança
sintênica**, não sobre posição de gene. Duas âncoras publicadas diferentes,
dois cromossomos — e isso não é contradição, é consequência do rearranjo que
Al-Hashimi et al. 2010 (`10.1093/molbev/msq098`) já havia reportado: o
pseudogene de enamelina da galinha aparece **em região não esperada pela
sintenia de mamífero**, sugerindo translocação em ancestral sauropsídeo.

Então a âncora arcossauria é a certa para esta pergunta, e a de Sire é a
mamífera. Medido agora pela API do NCBI Datasets, a mesma que o `s18_chrZ.py`
usa:

**Jacaré — *Alligator mississippiensis* (rAllMis1), cromossomo 3, tudo em fita
negativa:**

| Gene | Intervalo | Tamanho |
|---|---|---|
| ENAM | 157.293.890 – 157.307.957 | 14.067 bp |
| AMBN | 157.326.619 – 157.341.756 | 15.137 bp |
| AMTN | 157.359.292 – 157.374.521 | 15.229 bp |
| LPL | 157.488.689 – 157.512.736 | 24.047 bp |
| PSD3 | 157.601.092 – 157.780.011 | 178.919 bp |

Ordem no cromossomo: **ENAM → AMBN → AMTN → [114 kb] → LPL → [88 kb] → PSD3**.

**Frango — *Gallus gallus* (GRCg7b), cromossomo Z:**

| Gene | Intervalo |
|---|---|
| LPL | 55.215.002 – 55.232.438 |
| PSD3 | 55.231.641 – 55.383.432 |

**O vão intergênico LPL↔PSD3 passou de 88.356 bp no jacaré para −797 bp no
frango — os dois genes se sobrepõem.** Encurtamento de 89.153 bp entre dois
vizinhos conservados.

Três consequências, e a terceira é o próximo passo.

1. **A escolha de âncora do projeto está verificada.** LPL e PSD3 são de fato
   vizinhos sintênicos da região do cluster de esmalte no jacaré, e caem no chrZ
   no frango. A frase "a vizinhança arcossauria do esmalte mapeia no chrZ de
   *Gallus*" é **sustentada por medida**, não por asserção.
2. **O colapso de 88 kb é sinal real** naquele loco — deleção grande na linhagem
   da ave, na vizinhança imediata.
3. **Mas o cluster de esmalte nunca esteve dentro daquele vão.** ENAM/AMBN/AMTN
   ficam 114–181 kb **além** de LPL, fora da janela LPL..PSD3. Logo o colapso
   medido é *adjacente* ao antigo endereço do esmalte, não coincidente com ele.
   Para afirmar que existe cicatriz de deleção do cluster, a janela precisa ser
   estendida para além de LPL e comparada ao que há imediatamente além de LPL no
   chrZ do frango. **É uma chamada de API, e é o próximo passo do E1.**

Quarta medida, que muda como descrever o achado: **nenhum dos oito genes
(ENAM, AMBN, AMTN, MMP20, DSPP, AMELX, ODAM, ODAPH) está anotado em nenhum
assembly de *Gallus gallus*.** Então "relíquia" aqui não pode significar
"pseudogene anotado" — significa homologia detectável abaixo do limiar de
anotação. A frase de Métodos tem de dizer isso, e dizer com que ferramenta e que
limiar a homologia foi detectada.

E confirma Sire 2008 por outro caminho: ele reporta que **AMBN e ENAM
desapareceram após rearranjo cromossômico** na galinha. Desapareceram, não
pseudogenizaram. Isso é deleção, que é precisamente o modo dominante de perda
que Crombez et al. 2025 identificou (§5 E1, explicação **b**).

## A2. O quase-caso do Desenho 1, e o diferencial que sobra

Existe um trabalho muito próximo, e ele não é RERconverge:

**Abramyan, Li & Khansa 2025**, *"Reduction of Tooth Replacement
Disproportionately Affects the Evolution of Enamel Matrix Proteins"*
(`10.1007/s00239-025-10258-4`). **53 espécies** (41 Iguania — 24 acrodontes, 17
pleurodontes — mais 12 mamíferos onívoros), **seis genes** (AMBN, AMEL, AMTN,
ACP4, ENAM, MMP20), método **dN/dS por ramo no CODEML/PAML 4.8** com modelos
*free-ratio* e de 1 a 4 razões, grupos definidos **a priori** pelo modo de
substituição dentária, mais DIVERGE, PCOC e baseml. Efeito mais forte em AMEL.

Então o Desenho 1 continua **não sendo replicação**, mas o diferencial encolheu
e precisa ser dito com precisão: Abramyan compara **modelos de ω em grupos
fixados pelo pesquisador**; o Desenho 1 testaria **taxa relativa — resíduo
normalizado pelo fundo genômico — com nulo empírico por permulação**. São
inferências diferentes sobre o mesmo material. Essa frase, exatamente assim,
é a que justifica o trabalho na introdução.

E o levantamento sustenta a afirmação de ausência com as buscas feitas: seis
consultas literais no Europe PMC, e em cinco delas nenhum registro aplica taxa
relativa a perda de dente ou de esmalte.

### Mais dois genes para o painel

- **ODAM** — Springer et al. 2019 (`10.1186/s12862-019-1359-6`): inativado em
  placentários sem dente e sem esmalte, e em odontocetos. Entra no painel.
- **ACP4** — usado por Abramyan 2025. Não estava em nenhuma lista anterior.

O painel vai de 4 para **nove** candidatos: ENAM, AMBN, AMTN, AMELX, MMP20,
DSPP, ODAM, ODAPH, ACP4.

### A distinção que o painel permite, e que vale ouro

Meredith 2014 mede, em mamíferos sem esmalte (preguiça, tatu, porco-da-terra):
**DSPP intacto, MMP20/ENAM/AMTN/AMELX decaídos**. Isto separa operacionalmente
**perda de esmalte** de **perda de dente inteiro** — DSPP é dentina, e sobrevive
quando só o esmalte vai. Com 12 origens independentes e dois estados de perda, o
painel de nove genes dá um **padrão esperado por estado**, não só um teste. Essa
é a estrutura de um resultado forte.

## A3. Desenho 2 — a premissa que eu te dei está errada

Eu escrevi, na §6: *"se o viés de consenso for a explicação, a reconstrução
ancestral deve ser mais rígida que TODAS as pontas extantes"*. Duas coisas
derrubam isso.

**Primeira: "consenso é rígido por construção" não é lei geral.**
**Nixon, Lim, Sternke, Barrick, Harms & Marqusee 2024**, *Protein Science*
33:e5011 (`10.1002/pro.5011`) — em RNase H, o consenso do **alinhamento
completo** não foi mais estável e **não se comportou como proteína bem
enovelada**; o consenso de um conjunto **restrito a um clado** foi
significativamente mais estável e cooperativo. A estabilidade do consenso
depende do conjunto de entrada. Note os autores: Barrick, que em 2019 publicou
o desenho por consenso como estratégia de hiperestabilização
(`10.1073/pnas.1816707116`), e Marqusee/Harms, do lado do artefato.

Antes dele, e no mesmo sentido: **Risso, Gavira, Gaucher & Sanchez-Ruiz 2014**,
*Proteins* 82:887–896 (`10.1002/prot.24575`) — três variantes consenso de
β-lactamase contra as ressurreições precambrianas; **só uma das três pôde ser
expressa e purificada**, e foi **menos** estável que as ressurreições.
Conclusão dos autores: consenso não é proxy de ancestral.

**Consequência de desenho:** o comparador não pode ser "ancestral vs. pontas".
Precisa de **três braços** — o ancestral reconstruído, as pontas extantes, e o
**consenso explícito do mesmo alinhamento**, em duas versões (alinhamento
completo e clado restrito). Sem o terceiro braço o resultado não distingue
artefato de história.

**Segunda: a direção não está decidida na literatura.** "Ancestral mais rígido"
não é expectativa segura:

| Trabalho | O que achou |
|---|---|
| Babková et al. 2020, *CSBJ* (`10.1016/j.csbj.2020.06.021`) | dehalogenases ancestrais: ΔTm até 24 °C, e MD mostra o ancestral confinado a **um** estado conformacional enquanto as modernas alternam entre **dois** |
| Amritkar, Cuevas-Zuviría & Kaçar 2024, *MBE* (`10.1093/molbev/msae268`) | RuBisCO: complexos ancestrais **mais flexíveis** (RMSD e RMSF regionais maiores) |
| Zou et al. 2015, *MBE* (`10.1093/molbev/msu281`) | a dinâmica ancestral é a que **permite** promiscuidade: mais flexível no sítio |
| Modi et al. 2021, *Nat Commun* (`10.1038/s41467-021-22089-0`) | hinge-shift: a rigidez se **redistribui**, não aumenta monotonicamente |

**Reformulação do teste, agora direção-agnóstica:** a assinatura de artefato não
é rigidez, é **estar fora do envelope**. A pergunta certa é *o ancestral
reconstruído cai fora da faixa de variação conformacional medida entre as pontas
extantes, em qualquer direção?* Fora do envelope de toda a diversidade extante é
assinatura de construção; dentro dele, por mais extremo que seja, é compatível
com história.

E há magnitude medida para calibrar: **Williams, Pollock, Blackburne & Goldstein
2006**, *PLoS Comput Biol* 2:e69 (`10.1371/journal.pcbi.0020069`) — simulação
com ancestral verdadeiro conhecido. Viés de estabilidade (ΔG reconstruído menos
verdadeiro): **ML ≈ −1,5 kcal/mol** (~6 °C de Tm), parcimônia ≈ −0,4,
**bayesiano amostrado ≈ −0,05, isto é, sem viés**.

Isso entrega o controle interno mais limpo do desenho: **reconstruir o mesmo nó
por ML e por amostragem do posterior bayesiano e comparar os dois ensembles.**
Se o viés existir, o ML desvia e o bayesiano amostrado não. A diferença entre os
dois métodos **é** o sinal, medida contra um valor esperado publicado. Não
precisa nem de proteína ressuscitada em bancada.

## A4. BioEmu não sustenta a afirmação sozinho

**Zha et al. 2026**, *Int J Mol Sci* 27:2896 (`10.3390/ijms27062896`) — avaliação
externa do BioEmu contra 50 casos ATLAS e trajetórias de microssegundo da D. E.
Shaw. O que acerta: RMSF com PCC > 0,7 em 39/50; DCCM com MAE < 0,25 em 40/50;
mapas de contato preservados em 43/50.

O que erra é exatamente o caso de uso do Desenho 2:

- *"fails to predict a mutation-induced shift in conformational distribution"* —
  a divergência de KL entre ensembles **não separou** mutação *driver* de
  *passenger*, e em FGFR2/FGFR4 algumas *passenger* deram KL **maior** que
  *driver*;
- **subestima a flexibilidade global** contra MD longa;
- pior em proteínas com menos de 100 e mais de 400 resíduos;
- **viés para o estado de maior energia**, que os autores atribuem ao pré-treino
  no AlphaFold DB.

O último é o que mata se não for controlado, e é assimétrico entre os braços: **o
ancestral reconstruído não tem entrada no AFDB; as pontas extantes muitas vezes
têm.** O viés aponta na direção do resultado esperado. Um revisor que conheça
Zha 2026 derruba a figura.

### O protocolo defensável

Precedente metodológico exato: **Bhakat & Strauch 2026**, *JCIM* 66:7168–7178
(`10.1021/acs.jcim.6c01000`) — BioEmu como **semeador de MD** recupera estados
metaestáveis e deslocamentos de população por mutação melhor que semeadura por
AlphaFold2.

Então:

1. **BioEmu como triagem barata e semeador**, não como medida. Gera de uma vez:
   todas as pontas, o ancestral por ML, o `altAll`, ~100 amostras do posterior,
   o consenso do alinhamento completo e o consenso de clado restrito.
2. **OpenMM 8 é a medida que entra na figura.** Mesmo campo de força, mesma
   água, mesmo tempo, **≥3 réplicas independentes por sequência** — Knapp,
   Ospina & Deane 2018, *JCTC* 14:6127 (`10.1021/acs.jctc.8b00391`) mostra que
   diferença aparente entre sistemas emerge de simulação única por puro ruído.
   Incerteza reportada conforme LiveCoMS 1:5067 (`10.33011/livecoms.1.1.5067`).
3. **Três métricas, não uma.** ΔRMSF por resíduo com IC; divergência de
   distribuição (Jensen-Shannon / KL) sobre diedros ou coordenadas essenciais;
   entropia configuracional por MIE. Uma só não sustenta a afirmação.

**O ENCORE já está no venv.** As medidas de similaridade entre ensembles de
Lindorff-Larsen & Ferkinghoff-Borg 2009 (`10.1371/journal.pone.0004203`) estão
implementadas em Tiberti et al. 2015 (`10.1371/journal.pcbi.1004415`) como
`MDAnalysis.analysis.encore` — e o MDAnalysis 2.10.0 entrou no grupo `biofisica`
hoje. A ferramenta da métrica aceita em revisão por pares já está instalada.

Com 8 GB de VRAM o gargalo é a MD, não o BioEmu. E réplica curta em muitas
sequências é mais barata que trajetória longa em uma — o que favorece
precisamente o desenho de três braços.

## A5. Crítica metodológica ao RERconverge que precisa entrar em Métodos

Do levantamento, com DOI:

- **Partha et al. 2019**, *MBE* (`10.1093/molbev/msz107`) — os **próprios
  autores** do método mostram que ramo com média mais longa produz **maior
  variância da RER**, o que corrompe a confiança; corrigem com transformação de
  dados mais regressão ponderada. É a referência canônica do viés de comprimento
  de ramo, e ignorá-la é erro conhecido.
- **Mendes & Hahn 2016**, *Syst Biol* (`10.1093/sysbio/syw018`) — discordância
  de árvore de gene (ILS) produz **variação aparente de taxa de substituição**.
  Ataca de frente a premissa 1 do RERconverge, que é topologia única fixa para
  todos os genes. Com 250 Myr de divergência e ILS no ancestral de Archosauria,
  isto não é hipotético.
- **Treaster, Daane & Harris 2021**, *MBE* (`10.1093/molbev/msab226`) — TRACCER:
  método que compara **caminhos completos até o MRCA** de cada par, em vez de
  ramos singulares, justamente porque comparar ramos não fatora a proximidade
  filogenética entre as ocorrências do traço. Alternativa a considerar, não só a
  citar.
- **Macdonald et al. 2025**, *GBE* (`10.1093/gbe/evaf150`) — revisão PhyloG2P:
  nenhum método isolado acha todas as regiões; recomenda uso combinado.

E o tamanho de amostra do precedente, agora com origens independentes contadas:

| Estudo | Pontas | Origens independentes | Elementos testados |
|---|---|---|---|
| Chikina 2016 (vida marinha) | 59 placentários | **3** | 18.049 árvores de gene |
| Partha 2017 (subterrâneo) | 39 espécies | **4** | 18.980 genes + CNEs |
| Kowalczyk 2022 (ausência de pelo) | 62 mamíferos | binário | 19.149 genes + 343.598 CNRs |
| Gaspar 2026 (perda de voo em rails) | **7 espécies** | perdas independentes | > 11.000 alinhamentos |

**Chikina 2016 funcionou com três origens.** Isso derruba o limiar de ~5 que eu
havia fixado por conta própria e confirma que 12 origens é folgado. E Gaspar
2026 prova que 7 pontas é publicável — o que torna explícito que o obstáculo em
Archosauria **nunca foi o número de pontas, e sim a colinearidade** do caráter
com a bipartição Aves/Crocodylia.

Ponto de desenho que segue da definição do método: a RER é resíduo contra a
**média genômica** de comprimento de ramo. Testar nove genes exige, ainda assim,
**rodar milhares de árvores de gene para estimar o fundo**. Um painel de nove
genes isolado não tem denominador. O genoma inteiro entra como fundo, e os nove
genes são **lidos dentro do ranking** — não são a análise, são a leitura dela.

## A6. Os 45 PDFs já estavam no disco

`data/biblio/pdf` tinha 45 PDFs, entre eles Meredith 2009, 2011, 2013 e 2014,
Sire 2008 e Gasse 2015 — o artigo canônico deste projeto estava em disco há
meses. O índice `biblio` foi construído (4.959 chunks) e é incremental.

Consequência operacional: **antes de ir ao PubMed, consultar
`mcp__nvidia-nemo__library_search`.** Nesta sessão o `library_search` respondeu
com erro de ferramenta num dos dois levantamentos — vale conferir, porque o
acervo local é mais rápido e mais barato que qualquer busca remota.

## A7. Buscas resolvidas e as que sobraram

**Resolvidas nesta sessão:** `PhyloG2P` como termo de campo · 12 origens
independentes · data de 116 Ma com IC · painel de nove genes · o quase-caso
Abramyan 2025 · a controvérsia de hiperestabilidade com os dois lados nomeados ·
magnitude do viés de ML (−1,5 kcal/mol) · métricas de ensemble aceitas ·
veredito sobre BioEmu · posição das âncoras no jacaré e no frango.

**Sobraram, e estas ficaram com o Claude Science** (daemon de pé na porta 8765,
pendente de login no navegador). *Respondidas em 2026-10-17 por lit-scout, depois
de o Claude Science travar duas vezes — ver §A11.*

1. `ODAPH`/`C4orf26` em aves — estado não localizado.
2. `AMELX` em *Gallus*: não anotado em nenhum assembly. Onde cai o Ψ-AMEL, então?
   Sire 2008 o procurou no intervalo ARHGAP6–MID1, que na galinha é o
   cromossomo 1.
3. Mapa de recombinação do chrZ de *Gallus* com resolução suficiente para
   separar região de alta e de baixa recombinação — necessário para o confundidor
   do E1.
4. Potência de permulação em função do número de origens independentes — Redlich
   et al. 2024 (`10.1093/molbev/msae210`) traz análise de escalonamento das
   permulações categóricas com número de espécies e de categorias; conferir se
   responde.

**Nova, e é a mais importante:** o que há **imediatamente além de LPL** no chrZ
de *Gallus*, do lado oposto a PSD3? É ali que o cluster ENAM–AMBN–AMTN deveria
estar pela sintenia do jacaré. Uma chamada de API responde, e é o próximo passo
concreto do E1.

## A8. A medida que fecha a §A1 — e o próximo experimento

Feita na mesma sessão, depois do adendo. A âncora externa resolve a ambiguidade.

**Sintenia conservada nas duas espécies, mesma ordem, mesma fita (negativa):**

| | *Alligator mississippiensis* (rAllMis1, chr3) | *Gallus gallus* (GRCg7b, chrZ) |
|---|---|---|
| NRG1 | 156.331.250 – 156.595.276 | 54.642.186 – 55.109.784 |
| **ENAM** | 157.293.890 – 157.307.957 | **ausente** |
| **AMBN** | 157.326.619 – 157.341.756 | **ausente** |
| **AMTN** | 157.359.292 – 157.374.521 | **ausente** |
| LPL | 157.488.689 – 157.512.736 | 55.214.993 – 55.231.804 |
| PSD3 | 157.601.092 – 157.780.011 | 55.231.641 – 55.383.432 |

A ordem **NRG1 → [cluster de esmalte] → LPL → PSD3** é a mesma nos dois. Isso
transforma a âncora de duas em quatro genes e torna o endereço inequívoco.

**O intervalo NRG1 ↔ LPL:**

| | Tamanho | Conteúdo |
|---|---|---|
| Jacaré | **893.413 bp** | ENAM, AMBN, AMTN (três genes codificantes) |
| Frango | **105.209 bp** | **seis lncRNAs anotados, zero gene codificante** |

**O intervalo no frango tem 11,8% do tamanho do intervalo no jacaré.** Colapso
de 788.204 bp no endereço sintênico exato do cluster de esmalte.

Os seis lncRNAs do frango, naquela janela do chrZ (Ensembl GRCg7b):

```
ENSGALG00010014552   55.110.739 .. 55.112.266   (+)
ENSGALG00010014554   55.119.956 .. 55.122.615   (+)
ENSGALG00010014556   55.143.585 .. 55.145.940   (+)
ENSGALG00010014557   55.161.786 .. 55.165.166   (-)
ENSGALG00010014761   55.180.507 .. 55.195.001   (+)
ENSGALG00010014771   55.233.219 .. 55.258.832   (-)   [além de LPL, do lado de PSD3]
```

E a verificação negativa, nas duas bases: **nenhum dos nove genes do painel
(ENAM, AMBN, AMTN, MMP20, DSPP, AMELX, ODAM, ODAPH, ACP4) está anotado em
*Gallus gallus*** — nem no NCBI (nenhum assembly), nem no Ensembl (HTTP 400 para
todos os nove).

### O próximo experimento, e ele é de um dia

**Alinhar os cinco lncRNAs do intervalo NRG1↔LPL contra as sequências de ENAM,
AMBN e AMTN do jacaré.** Se algum mostrar homologia, ele *é* a relíquia — e a
relíquia é um lncRNA anotado, não um pseudogene não anotado.

Isso liga direto ao achado de Crombez et al. 2025 (§5 E1, explicação **b**), que
reporta **sobreposição de pseudogene com lncRNA** e Ka/Ks não-neutro em parte
deles. Se os lncRNAs do chrZ forem os restos do cluster, a pergunta "por que a
relíquia persiste?" ganha uma resposta mecanística candidata: ela persiste
porque **está sendo transcrita**.

Tudo o que o experimento precisa já está instalado, na camada pixi do WSL:
`blastn` (BLAST+ 2.17.0), `diamond`, `miniprot`, `mafft`, `seqkit`. Nenhum
download de genoma novo — as sequências vêm por accession.

Hipótese nula a bater: os lncRNAs não têm homologia com o cluster, e os
105 kb são sequência não relacionada que ocupou o endereço depois da deleção.
Esse resultado também é publicável, e é o controle que a afirmação exige.

---

# A9. RESULTADO — o alinhamento rodou, e a resposta é não

Executado em 2026-10-12. Scripts versionados:
`scripts/atavismo/s24_lncrna_relicto.py` (baixa), `s24_alinhar.sh` (alinha),
`s25_verificar_relicto.sh` (verifica). Saídas em `runs/atavismo/lncrna/`.
Ferramentas: BLAST+ 2.17.0, miniprot 0.18-r281, seqkit — camada pixi do WSL.

## A pergunta e a resposta

**Pergunta:** os lncRNAs que ocupam o endereço sintênico do cluster de esmalte
no chrZ do frango são os restos de ENAM, AMBN e AMTN?

**Resposta: não.** Mas as relíquias estão lá — fora dos lncRNAs, em sequência
que ninguém anotou.

## Os controles, primeiro

Nenhum hit vale sem eles.

| Controle | Resultado | Leitura |
|---|---|---|
| **Sanidade** — proteínas do jacaré contra o próprio genoma | ENAM 99,88% / AMBN 100% / AMTN 100% / LPL 100% | a janela está certa |
| **POSITIVO** — LPL do jacaré → chrZ do frango | e = **1,16e-57**, 64,9% em 154 aa, 7 HSPs; miniprot identidade **0,829** em 148.024–161.569 | o método **tem** sensibilidade a 116 Myr. Resultado negativo sobre o cluster passa a significar algo |
| **NEGATIVO** — DSPP do jacaré (chr2) → chrZ do frango | e = 8,6, 19 aa, bits 22,7 | piso de ruído em bits ≈ 21–24, e-value > 1 |

## O alvo — tblastn, proteína do jacaré contra o intervalo do frango

| Gene | Melhor e-value | Identidade | HSPs | Veredito |
|---|---|---|---|---|
| **AMBN** | **1,40e-09** | 54,8% em 31 aa | 3 | **relíquia codificante** |
| **ENAM** | **1,08e-07** | 41,2% em 68 aa | 8 | **relíquia codificante** |
| AMTN | 1,6 | 26,2% em 65 aa | — | dentro do ruído: **não detectada** |

Os HSPs de ENAM não estão espalhados: **quatro deles se agrupam em ~1,8 kb**
(posições 85.481–87.242 da janela). Hit espúrio não se agrupa.

## As três verificações

**A — número de cópias.** Cada região candidata, jogada de volta contra a
janela do próprio frango, dá **exatamente 1 hit**. Não são elementos
repetitivos. (Contraste: no primeiro alinhamento, a consulta
`ENSGALG00010014761` posição 3.469–3.580 bateu em **oito** lugares do jacaré —
assinatura de repetição, descartada.)

**B — melhor hit recíproco.** Regiões do frango traduzidas contra as cinco
proteínas do jacaré:

```
GENE   MELHOR_HIT   EVALUE      BITS   HSPs  VEREDITO
ENAM   ENAM_Amis    1.89e-07    42.7   8     RECIPROCO CONFIRMA
AMBN   AMBN_Amis    1.40e-09    44.7   3     RECIPROCO CONFIRMA
```

Para a região de ENAM, **oito HSPs apontam para ENAM**; AMBN só aparece a
e = 1,7, LPL a 7,8 e DSPP a 9,4. Para a região de AMBN, AMBN é o melhor a
e = 1,4e-09, e o próximo concorrente está a 4,2. Recíproco fecha nos dois.

**C — AMTN reprovado, e é assim que se sabe que o teste funciona.** O lncRNA
`ENSGALG00010014557` tinha dado um hit nucleotídico forte **dentro do span de
AMTN** no jacaré (72,8% em 158 bp, e = 4,41e-17, posição 157.363.280–157.363.437).
Traduzida, a mesma região não acerta nada: melhor e-value **1,2**, 8 aa, e o
acerto é em AMBN, não em AMTN. **Não é relíquia de codificante.** É elemento
não-codificante conservado ou repetição compartilhada. O teste C existia para
pegar exatamente isso, e pegou.

## Onde as relíquias estão — e onde não estão

Coordenada absoluta no chrZ, GRCg7b (janela aberta em 55.070.000):

| Elemento | Posição no chrZ |
|---|---|
| NRG1 (âncora externa) | 54.642.186 – 55.109.784 |
| lncRNA `ENSGALG00010014552` | 55.110.739 – 55.112.266 |
| lncRNA `ENSGALG00010014554` | 55.119.956 – 55.122.615 |
| lncRNA `ENSGALG00010014556` | 55.143.585 – 55.145.940 |
| **relíquia de ENAM** | **≈ 55.155.480 – 55.157.241** |
| lncRNA `ENSGALG00010014557` | 55.161.786 – 55.165.166 |
| **relíquia de AMBN** | **≈ 55.168.855 – 55.168.992** |
| lncRNA `ENSGALG00010014761` | 55.180.507 – 55.195.001 |
| LPL (âncora interna) | 55.214.993 – 55.231.804 |

**As duas relíquias caem nos vãos entre lncRNAs, não dentro deles.** A de ENAM
fica entre o 556 e o 557; a de AMBN entre o 557 e o 761. E a ordem no
cromossomo é **ENAM antes de AMBN**, que é a mesma ordem do jacaré
(ENAM 157,29 Mb → AMBN 157,33 Mb). **A colinearidade do par sobreviveu 116
milhões de anos.** Isso é um argumento independente de que os hits são reais:
ruído não preserva ordem.

## O que isto muda

**A hipótese de trabalho estava errada e o resultado é melhor.** A ideia era
que a relíquia fosse um lncRNA transcrito, com a transcrição explicando a
persistência. Não é. As relíquias estão em sequência intergênica **não anotada
em nenhuma das duas bases** — nem gene, nem pseudogene, nem lncRNA, nem no NCBI
nem no Ensembl.

Então a pergunta do E1 — *por que a relíquia persiste?* — perde a explicação
"porque está sendo transcrita" e fica com as outras duas da §5: baixa taxa de
deleção (Crombez et al. 2025) e ambiente de recombinação do chrZ
(Charlesworth 2020). As duas são sobre **ambiente genômico**, não sobre seleção.
O mapa de recombinação do chrZ deixa de ser item acessório e passa a ser o
próximo dado necessário.

**E AMTN some de verdade.** Dos três genes do cluster, dois deixaram relíquia
detectável e um não. Se a deleção foi um evento único, por que AMTN desapareceu
mais completamente? Pela ordem no jacaré, AMTN é o mais próximo de LPL — isto
é, o mais próximo da borda onde o colapso de 788 kb aconteceu. **Predição
testável: o grau de decaimento da relíquia correlaciona com a distância à borda
da deleção.** Nada disso precisa de genoma novo.

## Honestidade sobre os limites

- Identidade de 41% em 68 aa e de 55% em 31 aa é **pouca sequência**. O que
  sustenta a afirmação não é a identidade, é a conjunção: agrupamento dos HSPs,
  cópia única, recíproco confirmado, posição prevista a priori pela sintenia de
  quatro genes, ordem ENAM→AMBN preservada, e controle negativo mudo.
- Um revisor vai pedir **mais espécies**. A mesma medida em *Melopsittacus*,
  *Strigops* e *Chrysemys* — os genomas já listados na §2 do
  `PROJETO-ARCHOSAURIA.md` — diria se a relíquia é compartilhada por Neornithes
  ou é particular de *Gallus*. É o próximo passo, e é uma execução do mesmo
  script com outro accession.
- **Nenhuma destas relíquias foi anotada por ninguém.** Isso é oportunidade e é
  risco: oportunidade porque é achado novo; risco porque anotação ausente às
  vezes significa que a montagem naquele intervalo é ruim. Conferir cobertura e
  qualidade da montagem do chrZ em 55,15–55,17 Mb antes de escrever Métodos.
- `miniprot` não recuperou ENAM nem AMBN (só LPL, a 0,829 de identidade). Com
  `--outs=0.5` o limiar é alto para relíquia degradada. Rodar de novo com
  limiar mais baixo é teste adicional, não requisito.

---

# A10. O comparativo em seis espécies — a relíquia não é de *Gallus*

Executado em 2026-10-12. `scripts/atavismo/s26_relicto_comparativo.py` e
`s26_alinhar.sh`. Saídas em `runs/atavismo/lncrna/s26_*.tsv` e
`anotacao_painel.json`.

## Anotação do painel, espécie por espécie

| Espécie | Dentição | Painel de 9 anotado |
|---|---|---|
| *Alligator mississippiensis* | com dente | 7 de 9 |
| *Anolis carolinensis* | com dente | 6 de 9 |
| *Gavialis gangeticus* | com dente | 5 de 9 (scaffolds não posicionados) |
| ***Chrysemys picta*** | **sem dente** | **1 de 9 — só ENAM, no chr6** |
| *Gallus gallus* | sem dente | **0 de 9** |
| *Melopsittacus undulatus* | sem dente | **0 de 9** |
| *Strigops habroptila* | sem dente | **0 de 9** |

A tartaruga é desdentada e **mantém ENAM anotado**; nenhuma das três aves mantém
nada. Primeiro sinal de que o decaimento é graduado e mais avançado nas aves.

## Sintenia por espécie

| Espécie | Sequência | NRG1↔LPL |
|---|---|---|
| *Alligator* | chr3 | 893.413 bp |
| *Chrysemys* | chr6 | 337.095 bp |
| *Strigops* | **chrZ** | 293.404 bp |
| *Gallus* | **chrZ** | 105.218 bp |

**A sintenia cai no chrZ em duas aves independentes.** Não é peculiaridade de
*Gallus*. *Melopsittacus* ficou de fora: NRG1 não está anotado naquele assembly.
Em *Gavialis* e *Anolis* as âncoras estão em arranjo que não permite extrair
intervalo orientado, e o índice delas deve ser lido como limite inferior.

## O índice de retenção — bitscore, não e-value

O e-value escala com o tamanho do banco, e as janelas vão de 710 kb a 1,59 Mb.
Comparar e-value entre espécies seria inválido. O bitscore não depende do banco.
O índice é a soma dos bitscores dos três genes dividida pelo bitscore do
controle positivo da mesma espécie.

| Espécie | Dentição | ENAM | AMBN | AMTN | LPL (+) | DSPP (−) | **Índice** |
|---|---|---|---|---|---|---|---|
| *Alligator mississippiensis* | com | 1590 | 349 | 314 | 173 | 25,8 | **13,02** |
| *Gavialis gangeticus* | com | 1450 | 326 | 73,6 | 170 | 25,8 | **10,88** |
| *Anolis carolinensis* | com | 241 | 67,4 | 73,6 | 152 | 0 | **2,51** |
| *Chrysemys picta* | sem | 200 | 53,9 | 26,2 | 152 | 0 | **1,84** |
| *Strigops habroptila* | sem | 106 | 24,6 | 37,4 | 194 | 0 | **0,87** |
| *Gallus gallus* | sem | 42,0 | 56,6 | 23,5 | 200 | 0 | **0,61** |

**Com dente: 2,51 a 13,02. Sem dente: 0,61 a 1,84. Sem sobreposição.**

Controles em todas as seis: LPL positivo entre 152 e 200 bits; DSPP negativo
mudo ou ruído (25,8 bits nos dois crocodilianos, zero nas outras quatro).

## O controle que torna o resultado interpretável

*Anolis carolinensis* é o **mais distante** do jacaré entre as seis — Lepidosauria
está fora de Archelosauria, enquanto tartaruga e aves estão dentro. As proteínas
de consulta são do jacaré. Se o índice medisse distância filogenética, Anolis
seria o menor de todos.

Ele é **2,51** — acima da tartaruga (1,84) e muito acima das duas aves (0,87 e
0,61). **Anolis tem dentes.** O que separa os grupos não é distância; é o caráter.

## Decaimento específico por linhagem

| Gene | *Gallus* | *Strigops* |
|---|---|---|
| AMBN | **56,6** | 24,6 |
| AMTN | 23,5 | **37,4** |

Genes diferentes decaíram em linhagens diferentes de aves. A perda funcional é
única e compartilhada (Meredith 2014, ~116 Ma); o **apagamento do rastro é
independente por linhagem e por gene**. Isso é o esperado de degradação
estocástica depois da perda, e dá um fenótipo molecular **contínuo** onde antes
havia só presença e ausência — exatamente o tipo de variável que a §4 precisava
para sair do caráter binário.

## O que isto faz com o programa

1. **A relíquia é de Neornithes, não de *Gallus*.** Duas aves, duas vezes no
   chrZ, com relíquia detectável nas duas.
2. **As duas origens independentes de perda mostram retenção reduzida.** Aves e
   tartaruga caem no mesmo lado da separação, por caminhos evolutivos separados.
3. **O índice é um caráter contínuo pronto para método comparativo.** Com mais
   espécies ele entra no `phytools`/`geiger` como traço, e no RERconverge como
   fenótipo contínuo — que tem mais poder que binário e dispensa a codificação
   `foreground2Tree` que a colinearidade Aves/Crocodylia inviabilizava.
4. **O índice precisa de mais pontas antes de virar resultado.** Seis espécies
   mostram o padrão; não o testam. O teste é o Passo 5 da §A9.

## Limite honesto da normalização

Dividir pelo melhor HSP de LPL é imperfeito: o alinhamento de LPL fragmenta de
forma diferente entre espécies por causa dos íntrons, e o frango tem o **maior**
denominador do conjunto (200 bits, contra 173 do próprio jacaré). O efeito é
conservador — reduz o índice do frango, que já é o menor. A separação entre os
grupos é grande o bastante para não depender da escolha de normalização, mas
antes de publicar vale refazer com soma de todos os HSPs em vez do melhor.

---

# A11. O resto do painel e as quatro pendências (2026-10-17)

## O painel por intervalo — `s27`

O s26 só procura no intervalo NRG1..LPL (chr3 do jacaré). AMELX e MMP20 estão
no chr1, DSPP no chr2, ACP4 no chr5: rodados ali dariam zero sem significado.
`scripts/atavismo/s27_*` dá a cada gene o próprio intervalo sintênico, com três
controles: âncoras (vizinho nomeado mais próximo de cada lado, exigido nas sem
dente e no *Anolis*) como positivo; todas as proteínas do intervalo do jacaré
como **competidoras** (um trecho só conta para o alvo se nenhum vizinho alinhar
melhor ali); e **melhor hit recíproco** contra o proteoma inteiro do jacaré
(GCF_030867095.1, 48.040 proteínas). 66 de 69 trechos confirmaram; os 3
rejeitados eram parálogos de fora do intervalo.

**MMP20** — fração recuperada (soma de bits dos trechos recíprocos ÷ jacaré):

| *Alligator* | *Gavialis* | *Anolis* | *Chrysemys* | *Strigops* | *Gallus* |
|---|---|---|---|---|---|
| 100% (15) | 89% (13) | **57%** (9) | **30%** (7) | 9% (2) | 6% (1) |

Âncoras DCUN1D5/TMEM123 a 75–124 bits nas aves e na tartaruga: o intervalo está
certo. Mesma ordem do s26 num locus independente, em outro cromossomo.

**Sem conclusão:** AMELX (nem o *Anolis* sai do ruído), DSPP (controle falha em
*Gavialis* e *Anolis*; tartaruga tem 2 fragmentos recíprocos), ACP4 (sintenia
não estabelecida nas aves). ODAM e ODAPH não rodados: zero genes com esses
símbolos em Sauropsida no NCBI Gene.

## Calibração com pseudogene conhecido — `s27_calibra_amelx.sh`

O Ψ-AMEL da galinha (Sire 2008, GenBank EU340348) cai **dentro** do intervalo
MSL3..MID1 do s27, a 100% de identidade, em
**NC_052532.1:124.410.142–124.413.840** (fita −, GRCg7b) — primeira coordenada
dele numa montagem atual. O tblastn do AMELX do jacaré deixou **zero** HSPs ali;
com SEG desligado, 25,8 bits, indistinguível do ruído. Não é o filtro: é
divergência somada a decaimento. **Um zero do s27 não vale como ausência sem
controle positivo na mesma espécie.**

## As quatro pendências

| # | Resposta | Fonte | Confiança |
|---|---|---|---|
| 1. ODAPH em aves | Nada publicado sobre qualquer ave. Ortólogos em *Alligator* e *Anolis* (3 éxons; mamíferos perderam o éxon 2) achados por sintenia e RNA-seq, sem coordenada. Humano: RCHY1–THAP6–**ODAPH**–CDKL2–G3BP2 (4q21.1), fora dos clusters SCPP. Na galinha o bloco está quebrado; janela candidata chr4 ~33,6–34,2 Mb (consulta Ensembl, não dado publicado) | Liang et al. 2021, *Sci Rep*, PMID 33441959 | média |
| 2. Ψ-AMEL | Existe: inserção de 4 nt no éxon 2, stop no éxon 6. Meredith 2013 acha em mais 4 aves, com deleção de 1 pb compartilhada. Coordenada GRCg7b medida aqui (acima) | Sire et al. 2008, PMID 18775069; Meredith et al. 2013, PMID 23342979 | alta |
| 3. Recombinação do chrZ em 55,0–55,3 Mb | Sem taxa publicada nessa resolução. Melhor mapa: Elferink 2010, 1 marcador/50 kb no Z, 3,0 cM/Mb em média, só em macho, em galGal3 (janela ≈ 53,0–53,5 Mb ali). A janela fica no meio de Zq, longe de centrômero e telômero: espera-se taxa intermediária | Elferink et al. 2010, PMID 20141624 (Additional file 5, bins de ~500 kb); Groenen 2009, PMID 19088305 | média (inferida) |
| 4. Potência de permulação vs nº de origens | Não existe curva publicada. Redlich 2024 mede só tempo de execução (7 s com 2 categorias → 4.933 s com 6). A única curva contra nº de perdas é de Forward Genomics: 2–4 origens ≈ 20–40% de sensibilidade a 90% de precisão, e a idade da perda pesa mais | Redlich 2024, PMID 39404101; Prudent 2016, PMID 27222536; Saputra 2021, PMID 33739420 | alta / média |

**Correção à busca de 2026-10-16 do Claude Science:** PMID 30674270, que
apareceu como único hit de "ODAPH", é sobre
**ODAM** em 165 placentários (Springer et al. 2019, *BMC Evol Biol*). O de
C4orf26 é Springer et al. 2016, *Mol Phylogenet Evol* 95:34–45, PMID 26596502.
