# Projeto 2 — Atavismo em Archosauria: o que as aves guardaram do réptil

Plano de 2026-09-14. Irmão do [Projeto 1](PROJETO-ARCHOSAURIA.md), com o qual
compartilha espécies, grupos externos e armadilhas. Leia a §1 e a §8 antes de
tudo.

Pergunta do usuário: *analisar o que está mapeado nas aves, conferir o que elas
compartilham com outros répteis na ontogênese, e projetar in silico o que
precisaria mudar para acentuar caracteres basais como cauda e dentes.*

---

## 1. A pergunta, reformulada para ser respondível

"Reverter uma ave a um arcossauro basal" não é uma pergunta computável, pelo
mesmo motivo do Projeto 1: a trajetória histórica real não é recuperável, e
morfologia não está escrita apenas nos genes que sabemos listar. O que **é**
respondível, e é bastante:

| Pergunta | Respondível? |
|---|---|
| Quais genes do programa dentário ainda existem em ave | **Sim**, e já foi medido — §2 |
| Quais existem em crocodiliano e não em ave | **Sim**, direto |
| Se a perda foi de efetor ou de rede regulatória | **Sim**, e a resposta é nítida — §3 |
| Quais alvos moleculares plausivelmente controlam cada traço | **Sim**, com precedente publicado |
| Qual edição produziria dente/cauda num embrião real | **Não** — nenhum método prediz morfologia a partir de genótipo |
| Como era o dente do ancestral | **Não** — inferência de estado, nunca o registro |

**Reformulação adotada:** *"A perda de dentes e a redução da cauda nas aves
resultaram de perda de conteúdo gênico ou de mudança regulatória? Para cada um
dos dois traços, quais alvos moleculares teriam de mudar, em qual janela do
desenvolvimento, e qual a classe de intervenção correspondente?"*

Essa versão tem resposta, tem controle, e a §3 já responde a primeira metade.

---

## 2. O que foi medido nesta máquina

Script: `scripts/atavismo/01_inventario_genes.py`.
Saída: `runs/atavismo/inventario_genes.tsv` e `matriz_presenca.md`.
Fonte: API NCBI Datasets v2alpha, 2026-09-14, sem credencial.

Quatro módulos gênicos × dez táxons. Espécies e acessões herdadas do Projeto 1.

### Módulo A — efetores de esmalte e dentina

Genes presentes na anotação, por táxon:

| gene | Gallus | Melopsittacus | Taeniopygia | Alligator | Chrysemys | Anolis | Homo |
|---|---|---|---|---|---|---|---|
| ENAM  | – | – | – | sim | sim | sim | sim |
| AMBN  | – | – | – | sim | –   | sim | sim |
| AMELX | – | – | – | sim | –   | sim | sim |
| DSPP  | – | – | – | sim | –   | sim | sim |
| MMP20 | – | – | – | sim | –   | sim | sim |
| AMTN  | – | – | – | sim | –   | sim | sim |

**Nenhum dos seis aparece na anotação de nenhuma das três aves.** Os seis estão
íntegros e codificantes no jacaré.

> **Esta tabela mede anotação, não genoma, e foi superada.** A §4 mostra que
> `ENAM`, `AMBN` e `AMELX` **existem** em *Gallus*, como pseudogenes, num vão
> não anotado do cromossomo Z e no cromossomo 1. Leia as duas juntas: esta diz
> o que a consulta por símbolo enxerga, a §4 diz o que está no cromossomo.

Dois genes do módulo A resistem em ave — `SLC24A4` e `WDR72`. Não é
contra-exemplo: os dois têm função fora do dente, e é exatamente isso que a
§3 prevê.

### Módulo B — rede regulatória odontogênica

| clado | fração do módulo presente |
|---|---|
| Aves (3 espécies) | 25/26, 25/26, 25/26 |
| Crocodylia (*Alligator*) | 25/26 |
| Testudines | 26/26 |
| Lepidosauria | 26/26 |
| Mammalia | 26/26 |

A ave é indistinguível do crocodiliano aqui.

### Módulo C — alongamento axial e terminação da cauda

*Gallus* tem 23 dos 24, **mais** que o jacaré na mesma anotação. Todos os alvos
clássicos estão lá: `TBXT`, `GDF11`, `LIN28A/B`, `CDX1/2/4`, os quatro
parálogos `HOX13`, `WNT3A`, `ALDH1A2`, `CYP26A1`, `MSGN1`, `TBX6`, `RIPPLY2`.

### Módulo D — ciliogênese (via *talpid*)

Íntegro em todos os táxons, incluindo `C2CD3` — o gene do mutante *talpid²*,
que é o precedente experimental mais citado de germe dentário em galinha.

### Controle de convergência que valida a leitura

*Chrysemys picta* perdeu dentes **independentemente** das aves, há dezenas de
milhões de anos, num ramo separado. Ela mostra o mesmo padrão: módulo A
arrasado (3/12), módulo B intacto (26/26). Duas perdas independentes, a mesma
assinatura. Isso é o que separa sinal biológico de artefato de anotação.

---

## 3. O achado central: as aves perderam o hardware, não o software

> A ave não esqueceu **como** construir um dente. Ela perdeu o **material**.

O programa de padronização que decide onde um dente nasce está inteiro e sob
seleção — porque os mesmos genes fazem pena, escama, bico e glândula. A
pleiotropia protegeu a rede. Os genes de função única, que só serviam para
mineralizar uma coroa, não tinham o que os protegesse quando o bico queratinoso
assumiu a função, e degeneraram.

### O dado que sustenta isso, e não é narrativa

Consulta ao STRING (humano, 2026-09-14). Os dez parceiros de maior confiança de
`ENAM` são `AMBN`, `AMTN`, `AMELX`, `TUFT1`, `DSPP`, `MMP20`, `FAM20A`, `ODAM`,
`KLK4`, `FAM20C`. **Todos os dez são genes de mineralização dentária.** O
interactoma da enamelina é fechado sobre o próprio dente e não toca mais nada.

Os dez de `BMP4` são `BMPR1A`, `BMPR2`, `CHRD`, `BMPR1B`, `NOG`, `ACVR1`,
`ACVR2A`, `ACVR2B`, `BMP7`, `BMP2` — maquinaria de sinalização usada no embrião
inteiro.

A diferença não é de tamanho de rede, é de **identidade**: um módulo
auto-contido cuja única saída era a coroa do dente perde a função inteira de uma
vez e degenera junto, enquanto um nó de via geral está preso por tudo o mais que
depende dele. Isso explica por que os seis caíram **em bloco** e por que
`SLC24A4` e `WDR72`, que têm função fora do dente, ficaram.

Ressalva: o STRING devolve no máximo dez parceiros por consulta, então isto é
composição do topo do ranking, não grau de conectividade. A afirmação sustentada
é sobre com quem cada gene interage, não sobre quantos.

### As duas metades são problemas opostos

Isso não é um detalhe de classificação. É o que define o custo de cada metade
do projeto:

| | Dente | Cauda |
|---|---|---|
| Conteúdo gênico em ave | efetores mortos, alguns como ψ | tudo presente |
| Natureza da diferença | inativação de gene | regulação e tempo |
| Classe de intervenção | restaurar ORF, e/ou **knock-in** | edição *cis* / heterocronia |
| Precedente experimental | reindução parcial, sem esmalte (§6) | só *rumpless* Araucana, não validado |
| Evidência de lesão molecular | **medida** (§4) | **nenhuma publicada** (§5) |
| Risco de letalidade | menor, tecido local e tardio | alto, gene de eixo é pleiotrópico |

A inversão é o resultado mais interessante da varredura: **o traço cujos genes
sobreviveram é o mais perigoso de mexer, e o traço cujos genes morreram é o que
tem precedente de recuperação parcial.**

---

## 4. Dente — o que precisaria mudar

Três camadas, e elas costumam ser confundidas no noticiário.

**Camada 1 — induzir o germe.** Não exige gene novo, e esta é a camada com
respaldo mais forte: a rede está íntegra (§2, módulo B), demarca região
odontogênica potencial na mandíbula do pinto, e é **experimentalmente
reinduzível** — BMP exógeno restaura `Msx`, BMP+FGF gera estruturas epiteliais
`Shh`-positivas (§6, voto 3-0). O desenvolvimento está interrompido cedo, não
ausente. Alvos: o eixo `BMP4`/`SHH`/`PITX2`/`MSX1`, a via `WNT`/`LEF1`, e a via
ciliar `C2CD3`. Intervenção: *cis*-regulatória ou farmacológica local, em janela
HH definida.

Limite honesto desta camada: o que a literatura recuperou são estruturas de
identidade **ambígua entre dente e pena**, e o nó causal mais citado — perda da
expressão odontogênica de `Bmp4` — é hipótese de estudo único nunca confirmada
de forma independente.

**Camada 2 — mineralizar.**

> ⚠️ **Correção.** A primeira versão desta seção dizia que os seis genes foram
> deletados e que "não há o que corrigir". **Isso está errado, e o erro é meu.**
> A varredura da §2 mede ausência na *anotação* RefSeq, e eu a li como ausência
> no *genoma* — exatamente a armadilha de nível 1 que o agente `evo-devo` que
> escrevi nesta mesma sessão manda evitar. Pseudogene degenerado frequentemente
> não é anotado como gene, e some da consulta sem ter sumido do cromossomo.

**As relíquias foram encontradas, e as coordenadas são nossas.** Não citadas de
montagem antiga: medidas em `GRCg7b` nesta sessão, com `diamond blastx`
ultra-sensível e mascaramento desligado — proteína de matriz de esmalte é de
baixa complexidade e o filtro padrão a apaga.

Tabela completa em `runs/atavismo/reliquias_coordenadas.tsv`. Em *Gallus*:

| Relíquia | Localização (GRCg7b) | E-value | códons | STOPs internos |
|---|---|---|---|---|
| ψ-`ENAM` | `NC_052572.1` (**chrZ**) 55.155.480–55.157.294 | 4,8e-21 | 605 | **28** |
| ψ-`AMBN` | `NC_052572.1` (**chrZ**) 55.168.855–55.168.983 | 3,7e-09 | 43 | 0 (curto) |
| ψ-`AMELX` | `NC_052532.1` (**chr1**) 124.410.092–124.410.444 | 4,8e-08 | 117 | 2 |

**A lesão de `ENAM` é parada prematura múltipla no grande éxon C-terminal.** O
método não fabrica códons de parada, e isso foi controlado: jacaré e *Anolis*,
que têm dentes, dão **zero** stop nos mesmos alinhamentos; as aves dão 24 a 34.

Replicado em quatro ordens de ave, sempre no mesmo locus ancorado em `LPL`, com
a mesma ordem gênica do jacaré (`ENAM → AMBN → … → LPL`): *Meleagris*,
*Taeniopygia*, *Melopsittacus*, *Chlamydotis*. *Struthio* (Palaeognathae) tem
ψ-`AMELX` com 2 stops.

**`AMBN` não foi deletado do genoma do frango.** Era o item genuinamente em
aberto na literatura, e está resolvido: teste recíproco do fragmento de 42 aa
contra banco com a família SCPP completa de oito espécies mais iscas acerta
**só `AMBN`**, em quatro espécies independentes, e nenhum outro membro SCPP.

**Por que a varredura da §2 não os viu:** as duas relíquias do chrZ caem num vão
**não anotado** do RefSeq, entre `NRG1` e `LOC112530610`. Consulta em nível de
gene não enxerga ali. É a armadilha de nível 1, demonstrada em dado próprio.

Então a camada 2 se reclassifica:

| Gene | Estado medido em *Gallus* | Intervenção |
|---|---|---|
| `ENAM` | ψ no chrZ, 28 stops | restaurar ORF — **há o que corrigir** |
| `AMBN` | ψ no chrZ, presente | restaurar ORF; lesão ainda não localizada |
| `AMELX` | ψ no chr1, 2 stops, junto a `ARHGAP6` | restaurar ORF |
| `AMTN` | ψ achado em *Melopsittacus*, não em *Gallus* | indefinido |
| `DSPP`, `MMP20` | **não encontrados** — nível 1, não 3 | indefinido |

Restaurar 28 códons de parada não é edição de base pontual: é reescrita de éxon,
ou substituição do locus pelo ortólogo crocodiliano. A distinção prática entre
"reparo" e "knock-in" fica menor do que parece — mas **existe sequência
homóloga no lugar certo para ancorar**, e isso muda a estratégia.

Mesmo com ORF íntegra, promotor, janela e tipo celular continuam problema
separado — ORF restaurada num tecido que não a expressa não produz esmalte.

**Camada 3 — implantar e repor.** Alvéolo tecodonte e polifiodontia
(substituição contínua, que o crocodiliano tem e o mamífero não) são programas
próprios, com lâmina dentária persistente. Não vêm de brinde com as camadas 1 e 2.

O produto das três camadas não seria um dente de dinossauro. Seria um construto
quimérico ave-crocodiliano. Vale dizer isso na Discussão antes que o avaliador diga.

### O ramo, não a lista de espécies

```
                  ┌── Lepidosauria (Anolis) ······ ENAM/AMBN/AMTN íntegros, 0 stop
  Sauropsida ─────┤
                  │      ┌── Testudines (Chrysemys) ··· ψ-ENAM 19 stops, ψ-AMBN 9
                  └──────┤                                 ↑ perda INDEPENDENTE
                         │  ┌── Crocodylia (Alligator) ··· íntegros, 0 stop
                         └──┤
                            └── Aves ★ ··············· ψ-ENAM, ψ-AMBN (chrZ), ψ-AMELX (chr1)
                                Galliformes · Psittaciformes
                                Passeriformes · Otidiformes
```

A mudança é **no ramo das aves**, não na comparação ave↔jacaré. Crocodylia
retém o estado ancestral, Lepidosauria confirma que esse é o estado basal, e
Testudines dá o controle de perda independente: tartaruga e ave degradaram os
mesmos genes em separado, com o **mesmo tipo de lesão**. Convergência funcional
é muito mais parcimoniosa que artefato de anotação.

**Até onde a datação chega hoje: Neognathae**, não Aves. Quatro ordens
neognatas compartilham relíquia no mesmo locus. *Struthio* (Palaeognathae) tem
ψ-`AMELX`, mas `ENAM`/`AMBN` não foram buscados nele. Colocar a perda no
ancestral de Aves é hipótese, ainda não resultado.

### Como este trabalho reproduziu o erro de 2008, de propósito e por acidente

Vale registrar como resultado metodológico. A busca foi rodada primeiro
ancorada na sintenia **mamífera** (`SULT1E1`/`JCHAIN`, cromossomo 4) e produziu
um resultado limpo e completamente enganoso: um vão de 712 kb sem nada, com
montagem contígua e 0% de N, e zero HSP para os cinco genes — **com dois
controles positivos acertando na mesma corrida**. Um resultado assim convida a
concluir "deletado".

Foi exatamente essa inferência que levou o trabalho de 2008 ao erro que os
próprios autores corrigiram em 2010. As relíquias só apareceram ao reancorar na
sintenia **arcossauriana** medida no jacaré, onde o cluster fica a 115–181 kb de
`LPL`, e o bloco `LPL–PSD3–ABCA1–SLC44A1–FKTN` é colinear entre o jacaré e o
cromossomo Z do frango.

A lição, que vale para o Projeto 1 também: **escolher a âncora de sintenia do
clado errado produz ausência falsa com aparência de rigor.**

---

## 5. Cauda — o que precisaria mudar

Segunda campanha de pesquisa, focada só na cauda (108 agentes). A primeira
passada não cobriu este ângulo; esta cobriu.

### A resposta à pergunta que decide o projeto

**Não há evidência publicada de perda gênica, pseudogenização ou perda de
elemento não-codificante conservado no programa axial/caudal da linhagem-tronco
das aves.** Nenhuma fonte desta passada relata uma só.

Isso confirma a hipótese da §3 por caminho independente e com a força que ela
precisava: **dente e cauda são classes de problema diferentes.** O dente tem
lesão molecular documentada e agora medida; a cauda não tem nenhuma. Não existe
estratégia única que sirva para os dois.

### Contagem caudal — a escala do que mudou

| Táxon | vértebras caudais |
|---|---|
| *Alligator mississippiensis* | ~40 |
| *Anolis sagrei* | ~40 |
| *Archaeopteryx* | 20–23 |
| *Jeholornis* | 22 |
| *Sapeornis* | 6–7 livres + pigóstilo |
| *Gallus gallus* | **9** (5 livres + 4 no pigóstilo) |

*Sapeornis* é a primeira ave com cauda curta — a redução é em dois tempos, não
contínua.

### Duas correções que desmontam a proposta ingênua

**O pigóstilo não é embrionário.** As quatro vértebras que o formam estão
separadas e ossificadas aos 8 dias **pós-eclosão**; a fusão corre de distal para
proximal e só se completa por volta dos **5 meses**, desmontando e remodelando
discos intervertebrais já formados. É um **módulo separado** do que reduziu o
número de somitos caudais, e não pode ser citado como evidência sobre
eliminação de somito. Mexer no pigóstilo e mexer no número de vértebras são
dois projetos distintos, em janelas separadas por meses.

**A apoptose do broto caudal não é o gatilho.** Ela existe, mas é restrita ao
intestino caudal, ao remanescente da linha primitiva e ao mesênquima ventral — e
está **ausente** do tubo neural, da notocorda e dos somitos. Trabalho posterior
em pinto a coloca **a jusante** da subida do sinal retinoico: ela **conclui** a
terminação, não a dispara. Um desenho que ataque a apoptose ataca o epílogo.

### O elo fraco, declarado

Todo terminador molecular com evidência funcional real — o interruptor
tronco-cauda de `GDF11`, a dose `LIN28A`/*let-7* governando o tamanho do pool
progenitor caudal, a dose de `CDX1/2/4`, os `HOX13` — está estabelecido **em
camundongo**. Esta passada não achou **nenhum** dado de perda ou ganho de função
em ave para nenhum deles.

Ou seja: a §2 mostra que os genes estão todos lá em *Gallus*, mas **a
transferência camundongo→ave não está validada em nó nenhum**. Isso é o limite
central desta metade do projeto e precisa aparecer na Discussão.

A única lesão genética nativa de galinha que apareceu é o *rumpless* da raça
Araucana, ganho de função em `IRX1`/`IRX2` — mas é específica de raça, baseada
em **seis** aves sequenciadas, e **nunca validada funcionalmente**. É pista, não
alicerce.

### O que segue sem cobertura

O precedente hominoide (`TBXT`/AluY, Xia *et al.* 2024), a validade dos modelos
de linguagem genômica fora de mamífero, edição em PGC de galinha, atlas de
enhancer de frango, e a questão ética — **nenhum** produziu alegação
sobrevivente nas duas passadas. Ausência de cobertura, não resposta negativa.

### Mesmo assim: nenhum gene a acrescentar

A cauda da ave é curta porque o programa de
alongamento axial **termina cedo** e os somitos distais são eliminados; o
pigóstilo é fusão dos caudais remanescentes. O alvo é o **tempo**, não o conteúdo.

Alvancas candidatas, todas presentes em *Gallus*:

- `GDF11` — terminação do alongamento. Em camundongo, perda desloca limites e
  acrescenta vértebras.
- `LIN28A`/`let-7` — razão que governa duração do nicho progenitor.
- `CDX1/2/4` e os parálogos `HOX13` — fechamento do eixo posterior.
- `WNT3A` + `FGF8` contra `ALDH1A2`/`CYP26A1` — manutenção do nicho de
  progenitores neuromesodérmicos.
- Apoptose do broto caudal — a eliminação ativa, alvo separado dos acima.

**A armadilha que invalida a proposta ingênua:** nocaute completo de qualquer um
desses mata ou deforma o embrião muito antes do fenótipo de interesse, porque
todos agem antes e em outro tecido. A única forma defensável é **condicional no
tempo e no tecido**, declarando a janela HH. Uma proposta de nocaute simples de
`GDF11` não é um desenho ousado; é um desenho que não chega ao dia da leitura.

O efeito aqui é de **dose**, não de liga-desliga. É o oposto da camada 2 do dente.

---

## 6. O que a literatura verificada sustenta, e o que não

Primeira campanha de pesquisa adversarial de 2026-09-14: 5 ângulos, 26 fontes,
130 alegações extraídas, 25 testadas com três votos independentes cada,
**14 confirmadas, 11 derrubadas**, 10 depois da síntese. Registro completo em
`runs/atavismo/pesquisa_verificada.json`. Esta seção cobre só o dente — **a
segunda campanha, sobre a cauda, está na §5.**

### Sobreviveu com alta confiança (voto 3-0)

- **Perda única em Neornithes.** Os dentes mineralizados foram perdidos uma vez,
  no ancestral comum das aves modernas, evidenciado por mutações inativadoras
  *compartilhadas* em genes de esmalte e dentina, em quase todas as ordens.
  Note: a inferência repousa em ortólogos **presentes e lesionados**, não em
  ausência de gene — o que a torna imune ao artefato de falsa perda gênica que
  afeta microcromossomos GC-ricos em montagens de leitura curta de ave.
- **Escopo, e é onde quase todo mundo erra.** Isso vale para as aves de coroa.
  Entre as aves mesozoicas do tronco a perda de dentes foi adquirida **várias
  vezes, em separado**, e linhagens dentadas persistiram até o fim do Cretáceo.
  Dizer solto "as aves perderam os dentes uma vez só" **é excesso de afirmação**.
- **ψ-ENAM no cromossomo Z**, e a autocorreção de 2010 que desacredita o método
  de sintenia. É a base da correção da §4.
- **A rede regulatória não se perdeu.** `FGF8`, `PITX2`, `BARX1` e `PAX9` são
  expressos na mandíbula do pinto e demarcam região odontogênica potencial
  posicionalmente equivalente à lâmina dentária do camundongo. O
  desenvolvimento dentário aviário é **interrompido cedo, não ausente**.
  Isto confirma o módulo B da §2 por caminho independente.
- **A cascata latente é reinduzível.** BMP exógeno restaura expressão de `Msx`;
  BMP+FGF promove estruturas epiteliais `Shh`-positivas. Mas o que se recupera
  são estruturas de identidade ambígua entre dente e pena — **nunca dente
  mineralizado**.
- **O teto da reativação.** Como os genes estruturais estão inativados, germes
  reativados chegam no máximo a pré-dentina avançada; capa de esmalte é
  inatingível **por reativação**. Isto **não exclui** reintrodução transgênica.
  É predição a partir de conteúdo gênico, nunca demonstrada.

### Sobreviveu com confiança média (voto 2-1)

- A data de ~116 Ma é **atribuição, não fato assentado**. Cite pareada com a
  revisão para ~88–90 Ma.
- O toolkit está íntegro no grupo externo — mas a evidência verificada cobre
  **um gene**, `ENAM`, não o conjunto. Aqui a varredura da §2 acrescenta
  evidência própria: 6/6 no jacaré.
- A hipótese de Chen *et al.* 2000, de que a perda da expressão odontogênica de
  `Bmp4` é a causa específica da interrupção, é **nó causal de estudo único,
  nunca confirmado independentemente**, e há modelo concorrente posicional.

### Onze alegações derrubadas — e o que isso significa

A maioria caiu por **atribuição ou excesso de afirmação, não por ser falsa**.
Entre elas estão três alegações sobre o mutante *talpid²* (dentes de tipo
crocodiliano; suficiência de β-catenina; o modelo do centro de sinalização
lateral) e várias coordenadas de lesão. **Não cite como verificadas e não cite
como refutadas.** Se alguma virar load-bearing no texto, vá ao artigo primário.

### Advertências que viajam junto

- **Moeda de montagem.** Toda coordenada publicada vem de `WASHUC2`/`galGal3`
  (2008) ou de montagens de leitura curta de 2011–2012. Re-derive em `GRCg7b`.
- **Homologia fraca.** A relíquia de `DSPP` é alinhamento manual de ~50 pb a 54%
  de identidade, que o pipeline do GenBank não chamou. As conclusões de ausência
  de `AMBN`/`ENAM` usaram sondas de exon 2 de 54 pb.
- **Campo pequeno.** As fontes concentram-se em dois grupos de pesquisa. Há
  corroboração independente, mas a base é estreita.
- **Dois erros de citação que eu mesmo introduzi** no pedido de pesquisa e que
  contaminaram alegações: Chen *et al.* 2000 é **PNAS**, não *Development*; e o
  comparador crocodiliano de `AMEL` é ***Paleosuchus palpebrosus***, não
  *Alligator*.

### A metade da cauda ficou sem cobertura nesta passada

**Zero** alegações sobreviventes sobre `GDF11`, `CDX`, `HOX13`, `TBXT`,
`LIN28A`/*let-7*, pigóstilo ou apoptose do broto caudal. Foi o que motivou a
**segunda campanha, cujos resultados estão na §5** — e lá o ângulo foi coberto.

O saldo das duas passadas juntas, para a pergunta que decide o projeto:
**dente e cauda não têm a mesma classe de lesão molecular.** O dente tem lesão
documentada e agora medida em coordenada própria; a cauda não tem nenhuma
evidência publicada de perda gênica ou de elemento regulatório. A hipótese da
§3 deixou de ser só minha.

Seguem sem cobertura nas duas passadas: o precedente hominoide `TBXT`/AluY, a
metodologia in silico (modelos de linguagem genômica fora de mamífero, edição em
PGC, atlas de enhancer) e a questão ética. **Ausência de cobertura, não resposta
negativa.**

---

## 7. O que a stack ganhou para esta campanha

- `.claude/agents/evo-devo.md` — agente novo. A tabela de roteamento do
  `CLAUDE.md` não tinha linha para perda gênica, estado ancestral, polarização
  por grupo externo ou reversão de caráter. Codifica a escala de evidência da §7.
- `.claude/commands/atavismo.md` — comando `/atavismo`, genérico: serve para
  qualquer caráter perdido em qualquer clado, não só dente em ave.
- `scripts/atavismo/01_inventario_genes.py` — a varredura da §2, reexecutável.
- `scripts/atavismo/02_evo2_edicao.py` — prepara referência, variante e
  controles para pontuação com Evo 2. Testado.
- `runs/atavismo/` — saídas.

Endpoints NVIDIA conferidos para a fase de desenho: `evo2` e `evo2_generate`
respondem, e `enable_sampled_probs` é o caminho de log-verossimilhança para
comparar referência contra variante. `boltz2`, `rfdiffusion` e `proteinmpnn`
estão confirmados para a parte estrutural.

**Evo 2 validado ponta a ponta** nesta máquina em 2026-09-14, com sequência real
de *Gallus gallus*: a chamada retorna log-verossimilhança por base e grava o
vetor em disco. O caminho de pontuação não vai travar a fase de desenho.

**Por que o script 02 existe.** Um delta de log-verossimilhança sozinho não é
interpretável — não há escala absoluta, e o valor varia com comprimento, GC e
com o quanto o modelo conhece o clado. O script gera, junto de cada edição,
cinco controles negativos (mesmo número de substituições, sorteadas fora do
alvo, com semente fixa) e um controle positivo (lesão que sabidamente quebra).
A edição só é informativa se o delta dela sair **fora** da faixa dos controles
negativos; cair dentro significa "o modelo não distingue", e não "a edição é
segura". Ele também recusa a edição se a base de referência na coordenada não
bater — o que já pegou um erro de coordenada durante o próprio teste.

---

## 8. Armadilhas, e as quatro que já morderam nesta varredura

**Ausência de anotação não é ausência de gene.** É evidência de primeira ordem.
A escala está na Regra 2 do agente `evo-devo`. Antes de escrever "o gene não
existe", é preciso tBLASTn contra o genoma **com controle positivo**.

**O confundidor de anotação do Projeto 1 apareceu, exatamente como previsto.**
*Gavialis gangeticus* (EGAPx) e *Crocodylus porosus* (contig N50 de 32,9 kb)
deram contagens baixas em todos os quatro módulos, inclusive nos módulos que
ninguém perde. Isso é qualidade de anotação, não biologia. **Use
*Alligator mississippiensis* como referência de Crocodylia** — é a melhor
montagem anotada do clado no conjunto. Não cite os outros dois em
presença/ausência.

**Quatro artefatos, todos resolvidos — e nenhum era o que eu supus:**

- **`KLK4` é inovação de mamífero.** Zero ortólogo fora de Mammalia, incluindo
  anfíbio e peixe. **Exclua.** Cuidado extra: o Ensembl devolve um
  `ortholog_one2many` para o frango que é **parálogo** da família das
  calicreínas, não ortólogo de `KLK4`. É falso-positivo à espera.
- **`ODAM` não é inovação de mamífero — mas exclua assim mesmo, por razão
  melhor.** Tem 19 ortólogos em anfíbio e **zero em Sauropsida inteira**,
  inclusive no jacaré, que tem dentes e é polifiodonte. Ausência em ave é
  ausência em sáurópsida, não sinal de perda de dente.
- **`ACP4`, ao contrário, deve ficar.** Presente em Crocodylia e Lepidosauria,
  ausente em todas as aves e testudines. O contraste é informativo, e eu o
  tinha subestimado.
- **`FAM83H` era falha de consulta minha, nos dois sentidos.** O gene está
  íntegro em humano *e* em frango (`XP_015138735.2`, 1133 aa, chr2), com 548
  ortólogos e relação um-para-um. E as outras células vazias daquela linha têm
  explicação: em boa parte das aves o NCBI anota o ortólogo sob o símbolo
  **`SACK1H`**, não `FAM83H`. Refaça a linha com os dois símbolos.
- **`HES7`: a suspeita estava certa, a causa não.** Não é nomenclatura
  divergente — `c-hairy1` é `HES4`, e o registro `HES7` do frango foi
  **retirado pelo NCBI** por ter sido atribuído por engano a um parálogo `HES5`.
  A causa real é **buraco regional de montagem**: o bloco equivalente a 17p13.1,
  rico em GC e típico de microcromossomo, está sub-recuperado nas montagens de
  ave em conjunto. `PER1` aparece em 1 ave de 180, `AURKB` em 4 — e ninguém
  sustenta que as aves perderam Aurora B quinase. **`HES7` existe em ave**
  (anotado em *Strigops*, *Nipponia*, *Vanellus*, cinco *Falco*). Tire-o da
  matriz de perda.

A lição transferível: **resultado biologicamente implausível quase nunca é
descoberta.** Três dos quatro eram infraestrutura, não biologia.

**Limite ético e regulatório.** Tudo aqui é in silico. Qualquer passo em embrião
de vertebrado exige aprovação prévia do comitê de ética no uso de animais
(CEUA). O desenho computacional é legítimo e publicável por si; a passagem para
bancada não é decisão desta stack.

---

## 9. Itens em aberto

**Fechados nesta sessão:** `AMBN` existe (não foi deletado); coordenadas
re-derivadas em `GRCg7b`; módulo A elevado ao nível 3; cauda coberta por
pesquisa; os quatro artefatos resolvidos.

Em ordem de valor, o que sobra:

1. **Nível 4: datar a quebra.** Alinhar as quatro relíquias de ψ-`ENAM` de ave
   entre si (MUSCLE) e procurar **lesão compartilhada em coluna homóloga**. Sem
   isso só se pode dizer que cada uma está quebrada, não quando quebraram. É o
   passo mais barato com maior retorno — a evidência hoje cobre Neognathae, e
   isto a levaria ao ancestral de Aves. *Struthio* é a espécie-chave.
2. **Localizar a lesão de ψ-`AMBN`.** O fragmento recuperado tem 43 códons e
   zero stop — curto demais para carregar a lesão. Estender a janela.
3. **Resolver `DSPP` e `MMP20`.** "Não encontrados" é nível 1, não 3. A
   varredura é **ancorada em sintenia, não genômica** — o tBLASTn irrestrito do
   NCBI falhou o próprio controle positivo e foi descartado. Relíquia
   translocada para fora das janelas passaria despercebida. Caminho: montar
   banco local do genoma de *Gallus* no WSL e rodar sem âncora.
4. **Dois achados pendentes**, anotados para não virarem certeza por omissão:
   o acerto de `AMELX` em *Melopsittacus* no chrZ tem **0 stop em 156 códons** —
   não é perfil de pseudogene e não bate com `AMELX` estar no chr1 no frango;
   e não se verificou se ψ-`AMBN` de frango, peru e mandarim ocupa o mesmo éxon.
5. Complemento transgênico **mínimo** para esmalte: quais ortólogos são
   conjuntamente necessários e suficientes, e há knock-in restaurando
   mineralização em fundo edêntulo em **algum** sistema? Sem cobertura.
6. Modelos de linguagem genômica valem fora de mamífero? Sem cobertura nas duas
   passadas. Até lá o `evo2` é sinal exploratório, e o script 02 só reporta
   delta contra controle — de propósito.
7. Precedente hominoide `TBXT`/AluY, edição em PGC de galinha, atlas de enhancer
   de frango, caminho de revisão ética. Sem cobertura nas duas passadas.
8. **A ordem das mutações é inacessível, e sempre será.** Contar 28 stops diz
   que o éxon está morto há muito tempo; não diz qual veio primeiro — stops
   posteriores apagam a evidência dos anteriores. Nenhum método de sequência
   resolve isso. Também não se sabe se a pseudogenização precedeu ou seguiu a
   perda morfológica do dente, nem se houve estágio com dente sem esmalte.
9. Nada da §4 e §5 foi testado em modelo. São alvos com plausibilidade
   mecanística, não predições de fenótipo. Nenhum método prediz morfologia a
   partir de genótipo.
