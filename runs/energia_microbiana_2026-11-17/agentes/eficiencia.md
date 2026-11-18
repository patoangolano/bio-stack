# Qual microrganismo utiliza mais eficientemente a energia ambiental?

**Avaliação crítica — 2026-11-17.** Escopo: ecologia e fisiologia de organismos naturais; sem propostas de engenharia ou protocolos de cultivo. Busca dirigida, não revisão sistemática nem metanálise. Foram selecionadas **20 fontes científico-metodológicas verificadas**, distinguindo resultados experimentais, inferências de modelos e revisões. Os identificadores, níveis de acesso, critérios e números estão também em `eficiencia.json`.

## Conclusão executiva

**Não existe um vencedor universal demonstrado.** “Mais eficiente” pode significar conservar mais energia por reação, formar mais biomassa por recurso consumido, produzir mais por tempo, recuperar mais elétrons como corrente ou persistir com menor potência. São desfechos distintos, frequentemente conflitantes.

- **Luz, água e carbono inorgânico, sem doador reduzido adicional:** cianobactérias e microalgas são os candidatos diretamente pertinentes à produção primária fotossintética. Não há comparação aqui que escolha uma espécie campeã.
- **Luz e carbono orgânico:** *Rhodopseudomonas palustris* apresenta conservação de carbono elevada, mas isso não mede eficiência energética total e não a torna a fototrófica de crescimento mais rápido.
- **H₂, O₂ e CO₂ disponíveis:** *Cupriavidus necator* é um candidato pertinente à formação de biomassa no escuro; o H₂ fornece energia/redutor, o O₂ é o aceptor respiratório e o CO₂ é carbono. Os dados selecionados não estabelecem seu primeiro lugar em biomassa/J.
- **Amônia muito diluída e O₂:** *Nitrosopumilus maritimus* SCM1 destaca-se em afinidade; *Nitrospira inopinata* destaca-se pelo metabolismo comammox e por rendimento por mol de amônia frente a alguns comparadores. Afinidade e rendimento não são a mesma eficiência.
- **Anoxia com H₂/CO₂:** acetogênicos e metanogênicos são especialistas em conservação de energia sob restrições fortes; a ordem depende de condições, limiar, produto e linhagem.
- **Recuperação de elétrons como corrente a partir de orgânicos/H₂:** *Geobacter* é um candidato forte. Recuperar ~95% dos elétrons **não** significa recuperar 95% da energia.
- **Persistência com aporte mínimo:** microrganismos de sedimentos e consumidores de gases-traço são pertinentes. Sobrevivência não equivale a alta produtividade.

Versatilidade metabólica, quantidade de genes, tamanho do genoma e porcentagens de um atlas sem rastreabilidade **não estabelecem essa ordem**.

## 1. Contrato mínimo de comparação

Antes de classificar organismos, explicitar: linhagem; ambiente; recurso limitante; fonte de carbono; doador e aceptor de elétrons; fluxo de luz/energia; produto desejado; intervalo temporal; temperatura, pH e atividades; e fronteira do balanço.

| Métrica | Definição e unidade | Interpretação correta e principal armadilha |
|---|---|---|
| Energia livre da reação | ΔG = ΔG° + RT ln Q, kJ/mol de reação explicitada | Depende de estequiometria, atividades/concentrações, gases, pH e produtos; ΔG° ou ΔG°′ não é ΔG ambiental. |
| Conservação como ATP | mol ATP/mol doador ou reação; eventualmente η_ATP = n_ATP·ΔG_síntese_ATP/(−ΔG_catabolismo) | Exige acoplamento, ATP líquido e estado intracelular. Não existe ATP/doador universal por espécie. |
| Rendimento de carbono | Y_X/C = C-mol de biomassa nova/C-mol de carbono consumido | C-mol significa mol de átomos de C. CO₂ e acetato não têm o mesmo conteúdo energético. Incorporar C não é gerar energia. |
| Rendimento energético de crescimento | C-mol de biomassa/J disponibilizado, com base termodinâmica declarada | É rendimento dimensional, não percentual. Exige composição de biomassa e contabilização de todas as fontes. |
| Eficiência energética | energia ou exergia útil/entrada da mesma natureza, adimensional | Não misturar entalpia de combustão, Gibbs e energia de fótons sem uma convenção consistente. Declarar se biomassa, ATP ou produto é a saída. |
| Crescimento e produtividade | μ em h⁻¹; g de massa seca L⁻¹ d⁻¹ ou m⁻² d⁻¹ | Título em g/L não é rendimento nem produtividade; área de terreno não é área iluminada ou de eletrodo. |
| Afinidade e limiar | K_m/K_s e concentração mínima de consumo/crescimento | K_m menor não implica maior Y_X/J. Limiar de consumo não é automaticamente limiar de crescimento sustentado. |
| Manutenção | J C-mol_X⁻¹ h⁻¹, mmol ATP g_X⁻¹ h⁻¹ ou W/célula | Separar manutenção, reparo, dormência, morte/reciclagem e síntese de biomassa. Normalização por célula depende do tamanho celular. |
| Recuperação coulômbica | CE = ∫I dt/(F·n_e·n_substrato_consumido) | Quantifica carga recuperada; a energia elétrica exige tensão e tempo: E_el = ∫V_célula I dt. |
| Conversão solar da cultura | energia na produção líquida de biomassa/energia solar incidente integrada | Não é rendimento quântico de um fotossistema, nem F_v/F_m, nem eficiência do fóton já absorvido. |

Para converter massa seca em C-mol, usar `n_C = massa_seca × fração_mássica_C / 12,011`. Uma fórmula empírica genérica não substitui análise elementar; PHB/lipídios, cinzas e composição celular mudam a conversão.

### 1.1 ΔG, aceptores e concentrações

A energia vem de uma transformação acoplada, não do nome de uma molécula. Em quimiotrofia, doador e aceptor precisam ser considerados juntos; O₂, nitrato, CO₂, Fe(III) e eletrodos não são aceptores intercambiáveis. Para pares redox nas condições declaradas, ΔG = −nFΔE. Potenciais precisam estar na mesma referência; “+0,2 V versus Ag/AgCl” não é “+0,2 V versus SHE”.

Exemplos estequiométricos, **não cálculos de ΔG de um ambiente real**:

- Amônia → nitrito: NH₄⁺ + 1,5 O₂ → NO₂⁻ + 2 H⁺ + H₂O.
- Comammox: NH₄⁺ + 2 O₂ → NO₃⁻ + 2 H⁺ + H₂O.
- Metanogênese hidrogenotrófica: CO₂ + 4 H₂ → CH₄ + 2 H₂O.
- Acetogênese hidrogenotrófica: 2 CO₂ + 4 H₂ → CH₃COO⁻ + H⁺ + 2 H₂O.

As duas últimas consomem quatro H₂ por reação. Mantendo as demais atividades, baixar a atividade/fugacidade de H₂ dez vezes aumenta ΔG em **4RT ln(10) = 22,83 kJ/mol de reação a 298,15 K**. É uma derivação termodinâmica, não medição. Uma comparação baseada apenas no ΔG padrão pode, portanto, inverter sua relevância perto do limite energético. Usar atividades adimensionais, convenções coerentes para gases e solutos e não contar H⁺ duas vezes ao empregar ΔG°′ transformado.

Mais energia por mol de amônia no processo completo também pode explicar parte do maior rendimento de comammox: isso não prova melhor conservação **por joule**. O teste adequado depende da reação final, não apenas de “mol de doador”.

### 1.2 Por que ATP/doador não é uma constante universal

ATP líquido depende da cadeia respiratória, aceptores, bifurcação eletrônica, transporte reverso de elétrons, estequiometria de bombeamento de íons/ATP sintase, vazamentos, transporte, manutenção e biossíntese. O custo de sintetizar ATP depende de ATP/ADP/P_i e condições intracelulares. “−30,5 kJ/mol” não é um custo intracelular fixo aplicável a todas as células.

Um rendimento fracionário médio de ATP por reação é possível por acoplamento a ciclos iônicos e múltiplos eventos catalíticos; não se deve impor um “mínimo universal de uma molécula de ATP por reação”. Revisões de acetogênicos e metanogênicos discutem mecanismos e diferenças relevantes [S04–S05].

O exemplo quantitativo de *N. maritimus* é **anabólico**, não uma medição de ATP formado por NH₃: o ciclo HP/HB thaumarchaeal requer **4 mol de ligações anidrido de alta energia por mol de acetil-CoA**, contra **7** no cálculo do Calvin–Benson usado como comparação. Isso não inclui, como equivalentes universais, todos os custos de redutores, transporte e manutenção [S08, tabela 4].

### 1.3 Manutenção, potência e limites baixos

Uma decomposição útil é `q_S = μ/Y_X/S,max + m_S`. Mantida a manutenção, crescimento mais lento pode **reduzir** o rendimento observado por substrato, pois uma fração maior da entrada sustenta a célula em vez de formar biomassa. Adaptações que reduzem manutenção podem mudar esse resultado.

A correlação de Tijhuis et al. usa manutenção em **kJ C-mol_X⁻¹ h⁻¹**, depende fortemente de temperatura e foi ajustada a dados de crescimento de culturas [S01]. Não é um piso universal de viabilidade. Hoehler e Jørgensen descrevem ambientes com fluxos energéticos cerca de mil vezes menores que estimativas típicas de manutenção baseadas em cultura [S02].

Bradley et al. estimaram medianas de **2,23×10⁻¹⁸ W/célula** para heterótrofos aeróbios, **1,08×10⁻¹⁹** para redutores de sulfato e **1,50×10⁻²⁰** para metanogênicos em sedimentos [S03]. São resultados de um **modelo global de degradação, energia e abundância**, não calorimetria individual, mínimos universais ou um ranking experimental de espécies.

Baixa energia por reação e baixa **potência** são diferentes: `P ≈ taxa_de_reação × (−ΔG)`. Substrato muito diluído pode limitar fluxo mesmo quando sua oxidação por O₂ continua termodinamicamente favorável.

## 2. Evidência quantitativa e o que ela realmente compara

| Sistema e fonte | Resultado extraído | Alcance e ressalva decisiva |
|---|---|---|
| *R. palustris* CGA009, fotoheterotrofia com acetato [S09, tabela 2] | **1,75 ± 0,15 C-mol_X/mol acetato**, isto é, **0,875 ± 0,075 C-mol_X/C-mol acetato**; duplicação **8,4 ± 0,6 h**; médias ± DP, 8–10 réplicas biológicas | Conversão por dois C/mol de acetato calculada aqui; μ central correspondente ≈0,0825 h⁻¹. A luz é entrada adicional. O artigo usa fórmula elementar de outra linhagem, 42OL, na conversão de massa. Não é 87,5% de eficiência energética. |
| Mesmo estudo [S09] | Recuperação do balanço de C **93 ± 8%** e de elétrons **98 ± 9%**; ¹³C-MFA indica forte reassimilação de CO₂ | Fechamento do balanço não é rendimento de biomassa. Os ~94% de C destinado a biomassa no mapa de fluxos são inferência de MFA; não substituir o Y experimental de 0,875 pela porcentagem do mapa. |
| Seis culturas PNSB, incluindo *R. palustris* [S10] | Rendimento próximo de **1 g COD_biomassa/g COD_removido** em todas; *Rb. capsulatus* isolado apresentou μ de **2,3–3,8 d⁻¹**; reator inoculado com enriquecimento chegou a **1,7 g DW L⁻¹ d⁻¹** | COD mede equivalentes de oxidação, não C-mol nem joules. Luz também entra no sistema. A produtividade máxima é de comunidade em reator, **não de cultura pura de R. palustris**. Não há fundamento para proclamá-la a mais rápida. |
| *N. maritimus* SCM1 versus *Nitrosococcus oceani* [S08] | **1,3 versus 0,8 g de massa seca/mol de amônia oxidada**; duplicação **28 versus 15 h** | Comparação no mesmo estudo e condições declaradas iguais: rendimento maior com crescimento mais lento. Não converter para C-mol sem composição apropriada. É contraste local, não todos AOA versus todos AOB. |
| *N. maritimus* SCM1 [S06] | K_m de **133 nM de amônio total**; limiar de consumo **≤10 nM total** | Afinidade/limiar de oxidação, não porcentagem energética. A aproximação de **~3 nM NH₃ livre** é reportada na comparação de S07. Especiação, pH, temperatura e salinidade importam. |
| *N. inopinata* [S07] | K_m(app) **0,84 ± 0,17 µM de amônio total**, equivalente a **63 ± 10 nM NH₃**; V_max **14,8 ± 1,2 µmol N mg_proteína⁻¹ h⁻¹**; DP, n=6 | Não comparar 63 nM de NH₃ diretamente com 133 nM de amônio total. Agregados podem tornar K_m aparente uma estimativa superior. Não extrapolar a todas *Nitrospira*. |
| *N. inopinata* versus AOA não marinhos [S07] | Rendimento de proteína/mol de amônia **32,3%** maior que *N. gargensis* e **29,7%** maior que *N. viennensis*; n=4 por organismo | Oxidação completa versus parcial, temperaturas distintas e composição/medida de proteína limitam uma interpretação como eficiência por J. Não foi um confronto de rendimento com *N. maritimus*. |
| *C. necator* ATCC 17697 [S11] | Fase lag autotrófica **6–12 h**, prolongada em condições gasosas inibitórias; experimentos em duplicata | O texto reconhece a necessidade futura de medir consumo de substratos/gases durante os ensaios. Portanto, OD/lag não fornecem Y_X/H₂ ou eficiência energética fechada. Esta linhagem não deve ser automaticamente identificada com H16. |
| *C. necator* H16 [S12] | Biomassa seca total **13,1 ± 0,4 g/L**, com **79 ± 2% PHB** | Título e conteúdo de reserva/produto, não Y_X/J ou μ. Massa seca rica em PHB não é diretamente comparável a massa celular estrutural. Valores extraídos do resumo; n e definição do erro não foram auditados no suplemento. |
| *Nannochloropsis* sp. CCAP 211/78, quatro sistemas externos [S13, tabela 1] | Eficiência solar **média 1,1–2,7%**, máximas reportadas **1,5–4,2%**; produtividade média **9,7–20,5 g DW m_solo⁻² d⁻¹**, máximas **14,0–27,5** | Mesma linhagem/local, comparação de sistemas, não de espécies. A média solar mais alta foi de painel plano (2,7%; 36 dias); o máximo solar de 4,2% foi do tubular vertical (janela de 9 dias), não média anual nem máximo simultâneo de produtividade. |
| *G. sulfurreducens*, acetato → eletrodo [S16] | Recuperação média de **~95%** dos elétrons; um balanço ilustrativo do artigo dá **96,8%** | É recuperação coulômbica. Tensão, catodo, perdas e energia auxiliar determinam eficiência energética. Sistemas potenciostáticos e pilhas geradoras não são equivalentes energeticamente. |
| Gases-traço atmosféricos [S18–S19] | Revisão usa H₂ atmosférico de **~530 ppbv (0,53 ppmv)**; solos antárticos oxidaram H₂/CO e apresentaram fixação quimiossintética de CO₂ | O valor de H₂ é o contexto da revisão de 2015, não uma medição global de 2026. Evidência frequente é manutenção/persistência; S18 também sustenta produção primária comunitária. Não afirmar que todos apenas mantêm, nem que toda hidrogenotrófica cresce do ar. |

**Ausência de incerteza na extração não significa variância zero.** Não se calcularam intervalos combinados, efeitos padronizados ou médias entre estudos.

### 2.1 Eficiência solar: fotossistema não é cultura

Um rendimento quântico relaciona eventos fotoquímicos a fótons em uma etapa definida; F_v/F_m é uma métrica de fluorescência de PSII. Nenhum deles mede diretamente energia armazenada em biomassa dividida por toda a radiação solar incidente.

O balanço da cultura inclui espectro não aproveitado, reflexão, dissipação, saturação, fotoproteção, respiração diurna/noturna, alocação, autossombreamento e perdas de biomassa. Comparações com fotovoltaicos ou entre fototróficos exigem produto, período e fronteira comuns [S14].

Em S13, a eficiência usa energia de combustão assumida de **22,5 kJ/g de biomassa** e luz incidente por área ocupada, com conversão PAR/solar. Não é um balanço energético líquido completo que desconte bombeamento, aquecimento/resfriamento e colheita. A geometria e a definição de área também precisam acompanhar os números. Não transpor os valores dessa microalga para todas as cianobactérias.

### 2.2 Eletrodo como aceptor versus doador

- **S16:** o organismo oxida acetato ou H₂ e entrega elétrons ao eletrodo; o eletrodo é **aceptor** para a bactéria. A eletricidade recuperada é saída do sistema gerador, não alimento elétrico da célula.
- **S17:** elétrons fornecidos por eletrodo foram usados por *G. metallireducens* para reduzir nitrato e por *G. sulfurreducens* para reduzir fumarato. Aqui o eletrodo é **doador**. Isso não demonstra, por si, fixação autotrófica de carbono nem eficiência elétrica de crescimento.
- A energia de uma corrente é carga × diferença de potencial. CE alta com tensão útil baixa pode coexistir com recuperação energética baixa. Um potenciostato requer auditoria da energia fornecida e da contrarreação [S15–S17].

### 2.3 “Viver do ar” não é converter CO₂ em energia

H₂ e CO podem ser doadores; O₂ é aceptor; CO₂ é principalmente fonte de C nesse contexto. O estoque atmosférico pode ser continuamente reposto, mas o **fluxo local** é finito e limitado por difusão, água, superfície e atividade.

S19 descreve estratégias de persistência sob fome energética e também contribuição mixotrófica em alguns organismos. S18 mostra que comunidades antárticas podem realizar fixação quimiossintética, evitando a generalização oposta de que gases-traço só servem para manutenção.

Possuir uma hidrogenase não basta: a revisão distingue hidrogenases de alta afinidade de homólogos de *R. eutropha/C. necator* que, nos estudos examinados, não sustentavam consumo atmosférico. Isso é um exemplo direto de **homologia funcional ≠ desempenho ecológico** [S19].

## 3. Ranking condicional: ordens parciais, não um escore

“Prioridade” abaixo significa adequação ao recurso/desfecho especificado, não superioridade universal nem previsão de desempenho de campo.

| Cenário e desfecho | Prioridade condicional | O que ainda impede uma ordem completa |
|---|---|---|
| Luz + água + C inorgânico; produção líquida sem combustível reduzido externo | **Cianobactérias/microalgas oxigênicas** entre os grupos propostos | Espécie, clima, luz disponível, nutrientes, arquitetura e respiração. *R. palustris* não substitui a água como doador de elétrons fotossintético. |
| Luz + orgânicos; conservar carbono na biomassa | **PNSB, incluindo R. palustris**, são fortes candidatos | Não há vitória exclusiva: S10 mostra rendimentos similares entre PNSB. Luz deve entrar no denominador energético. |
| Mesmo contexto; maximizar taxa de formação | Não presumir *R. palustris*: **Rb. capsulatus** foi mais rápido nas culturas comparadas de S10 | Produtividade máxima de comunidade não é produtividade da espécie isolada; falta padronização energética entre todos os grupos. |
| H₂ + O₂ + CO₂ no escuro; formar biomassa | **C. necator** é candidato prioritário por elegibilidade fisiológica | Não há aqui comparação pareada e fechada em C-mol/J contra todas as alternativas; custo/origem e transferência de H₂/O₂ são essenciais. |
| NH₃ muito diluído + O₂; capacidade de captar substrato | **N. maritimus SCM1** precede *N. inopinata* nos valores de afinidade citados | Estudos/condições distintos; isso não classifica rendimento energético nem todas as linhagens AOA/comammox. |
| Amônia → nitrato; rendimento por N em comparadores estudados | **N. inopinata** supera os AOA não marinhos comparados em S07 | Reação final diferente; falta comparação por ΔG efetivo e biomassa-C. *Nitrospira* canônicas oxidantes só de nitrito não são intercambiáveis com comammox. |
| H₂/CO₂ em anoxia; persistência ou acetato/CH₄ como produto | **Acetogênicos ou metanogênicos**, conforme nicho e produto | Sem produto comum não há ordem. Entre metanogênicos, grupos com citocromos podem ter maior rendimento e maior limiar de H₂, ao contrário de grupos sem citocromos [S05]. |
| Orgânicos/H₂ e aceptor eletrodo; recuperar carga como corrente | **Geobacter**, no sistema estudado | CE não ordena energia líquida, potência, biomassa ou captação de eletricidade. |
| Potência mínima para persistir | **Especialistas de sedimentos ou gases-traço** pertinentes ao ambiente | Sem um mínimo universal ou vencedora taxonômica identificada; baixa potência não é alta produtividade. |

Se o objetivo for **alternar fontes de energia quando o ambiente muda**, *R. palustris* é uma candidata à versatilidade. Trata-se de robustez/plasticidade, não de um resultado em η, Y ou produtividade. A disponibilidade das fontes e os custos de manter a maquinaria precisam entrar na análise.

## 4. Pareto qualitativo

Dominação de Pareto exige o **mesmo cenário viável** e melhora em pelo menos um objetivo sem piora nos demais. Não foi estimada uma fronteira quantitativa: faltam vetores de recursos e medições harmonizados. O quadro identifica especializações e possíveis compromissos, não prova que todas as opções sejam não dominadas.

| Grupo | Vantagem apoiada/pertinente | Compromisso ou dimensão não estabelecida |
|---|---|---|
| *R. palustris* | Alto rendimento de carbono fotoheterotrófico e balanço redox estudado | Luz é custo adicional; alto carbono retido não demonstra máximo em C-mol/J ou velocidade. |
| *C. necator* | Autotrofia no escuro com H₂/O₂; biomassa e reservas quantificáveis | Necessita combustível e aceptor; mistura biomassa estrutural/PHB e denominadores gasosos incompletos impedem dominação demonstrada. |
| Cianobactérias/microalgas | Uso da água como doador em produção primária movida por luz | Variabilidade diária/sazonal e perdas ópticas/respiratórias; números de reação isolada não preveem cultura. |
| *N. maritimus* | Afinidade elevada e via anabólica econômica | Crescimento mais lento que AOB no contraste S08; depende de NH₃/O₂ e não vence todos os eixos. |
| *N. inopinata* | Alto rendimento por amônia e elevada afinidade frente a vários comparadores | Baixo V_max relativo e reação completa; comparação por N não resolve eficiência por energia. |
| Acetogênicos/metanogênicos | Conservação energética em anoxia; mecanismos adaptados a restrição | Limiares, rendimento e produtos variam; maior conservação por reação pode demandar maior disponibilidade de H₂. |
| *Geobacter* | Recuperação elevada de carga com eletrodo aceptor | Carga, energia elétrica útil, potência e rendimento de biomassa são objetivos diferentes. |

Um compromisso **medido localmente** é especialmente claro em S08: *N. maritimus* produz mais massa por mol de amônia, mas *N. oceani* duplica mais rápido. Não há dominação entre os dois nesses dois objetivos. Isso não transfere automaticamente uma relação de Pareto às outras espécies.

## 5. FBA: por que quatro modelos não fornecem quatro eficiências comparáveis

FBA resolve `max cᵀv`, sujeito a `Sv = 0` e limites inferiores/superiores. O resultado pertence ao conjunto **rede + meio/limites + objetivo + pressupostos**, não apenas ao nome da espécie [S20].

Conforme atualização da sessão principal, existem modelos locais `iJN678`, `iHN637`, `iAF692` e `iCN1361`, e está prevista auditoria isolada de objetivo, meio padrão, fechamento de importações e ATP. **Não executei esses modelos, não conferi suas atribuições taxonômicas e não antecipo o resultado dessa auditoria.**

Para qualquer comparação posterior seria necessário:

1. Confirmar linhagem, versão, origem, integridade e unidade de biomassa de cada modelo.
2. Fixar um **mesmo desfecho físico**: maximizar crescimento, ATP, biomassa-C/J ou produto são perguntas diferentes. Dois fluxos chamados “biomass” podem representar composições e custos diferentes.
3. Declarar e harmonizar o vetor ambiental: C, doadores, aceptores, nutrientes, luz, H₂, CO₂ e disponibilidade de O₂; não abrir todas as trocas.
4. Distinguir concentração de recurso de **limite de fluxo de importação**. `model.medium` contém limites de fluxo, não concentrações ou afinidades cinéticas.
5. Comparar por orçamento comum relevante, idealmente energia/exergia disponível e fluxos de matéria. Igualar “10 mmol de substrato” não iguala carbono, elétrons ou joules.
6. Conferir energia de luz, estequiometria fotoquímica, bombeamento de íons, P/O e transporte reverso; um fóton barato ou reação reversível indevida cria vantagem artificial.
7. Auditar manutenção não associada ao crescimento (NGAM) e ATP associado ao crescimento (GAM), nas mesmas unidades e com intervalos fisiológicos justificados.
8. Verificar balanço de massa/carga, direcionalidade, reações de demanda e ciclos geradores de energia. ATP sem entrada pode indicar artefato, mas exige primeiro conferir todas as fontes, fótons e reservatórios.
9. Não tratar remoção de loops como imposição completa de termodinâmica. ΔG e concentrações precisam de restrições adicionais; FBA convencional não os infere.
10. Avaliar soluções alternativas/FVA e sensibilidade a meio, manutenção e acoplamento. Um ótimo numérico não é intervalo de confiança nem fenótipo experimental.
11. Confrontar predições com rendimentos, fluxos e crescimento medidos em condições compatíveis; separar inferência ¹³C-MFA, FBA e medições diretas.
12. Se houver modelos viáveis e dados harmonizados, construir fronteiras condicionais por cenário. Um meio único que inviabiliza especialistas não demonstra inferioridade global.

**Exclusão explícita:** segundo a auditoria comunicada pela sessão principal, `scripts/bioenergetics_efficiency_analysis.py` não usa COBRApy nem solver; contém valores fixos e álgebra sem proveniência adequada. Portanto, os percentuais de `runs/bioenergetics_efficiency_atlas.json` **não são resultados de FBA, não são medições e não sustentam ranking**. Não foram utilizados aqui. A disponibilidade de 79 entradas no painel de genomas, ou de genes anotados nessas entradas, também não valida esses percentuais.

Fechar importações testa consistência/possibilidade de fontes indevidas, não mede manutenção real. Crescimento ou ATP previstos em um modelo também não autorizam um ranking experimental entre espécies.

## 6. Qualidade da evidência e lacunas

**Pontos fortes:** comparações dentro do mesmo estudo; distinção de denominadores; linhagens explícitas; rendimentos com balanços de C/elétrons; uso de réplicas quando reportadas; separação entre dado medido e modelo.

**Limites que afetam a conclusão:**

- Não foi encontrado um estudo que confronte todos os grupos no mesmo orçamento energético e com desfecho comum. A ausência de comparação impede declarar um campeão; não demonstra equivalência.
- Temperatura, pH, salinidade, estado fisiológico, história de crescimento, composição da biomassa e arquitetura não estão uniformizados.
- Em alguns resultados só foi acessado o resumo. A ausência de n, intervalo ou definição de erro na extração é lacuna de auditoria, não prova de que o artigo não os reporte.
- Dias consecutivos de uma instalação não são réplicas independentes de instalações; máximas de janelas curtas não são médias anuais.
- Em S09 a conversão para C-mol depende de fórmula empírica importada de outra linhagem. Em S07 rendimento é proteína/mol N; em S10 é COD; em S12 massa inclui grande reserva de PHB. Não combinar esses números como se fossem Y_X/C idênticos.
- A estimativa de potência de S03 depende de modelos de degradação, reações atribuídas e abundância celular; não identifica uma espécie com “menor manutenção medida”.
- Faltam balanços fechados e comparáveis de H₂/O₂, energia incidente/absorvida, energia de produtos e insumos auxiliares para classificar C-mol_X/J nos grupos propostos.
- Faltam ΔG ambientais calculados com atividades, produtos, aceptores e condições de cada comparação. Não se forneceram falsos ΔG in situ a partir de valores padrão.
- Não foi estabelecido um ranking independente entre espécies de cianobactérias e microalgas, entre todos os acetogênicos/metanogênicos ou entre todas as linhagens do mesmo gênero.

**Recomendação para a síntese principal:** apresentar um mapa de candidatos por recurso e objetivo e um quadro de compromissos. Manter em separado versatilidade, rendimento energético, conservação de carbono, taxa e persistência. Usar “não determinado” onde falta medida; nunca substituir por zero, escore de genes ou percentual do atlas anterior.

## 7. Fontes verificadas

Nível de acesso efetivamente usado: **T** = texto integral nas partes citadas; **R** = resumo/metadados, sem alegar revisão integral do suplemento.

1. **S01 — Tijhuis L, van Loosdrecht MCM, Heijnen JJ (1993).** *A thermodynamically based correlation for maintenance Gibbs energy requirements in aerobic and anaerobic chemotrophic growth.* Biotechnology and Bioengineering. DOI: [10.1002/bit.260420415](https://doi.org/10.1002/bit.260420415). PMID: [18613056](https://pubmed.ncbi.nlm.nih.gov/18613056/). **R**; correlação de dados de culturas, não mínimo universal.
2. **S02 — Hoehler TM, Jørgensen BB (2013).** *Microbial life under extreme energy limitation.* Nature Reviews Microbiology. DOI: [10.1038/nrmicro2939](https://doi.org/10.1038/nrmicro2939). PMID: [23321532](https://pubmed.ncbi.nlm.nih.gov/23321532/). **R**; revisão.
3. **S03 — Bradley JA et al. (2020).** *Widespread energy limitation to life in global subseafloor sediments.* Science Advances. DOI: [10.1126/sciadv.aba0697](https://doi.org/10.1126/sciadv.aba0697). PMID: [32821818](https://pubmed.ncbi.nlm.nih.gov/32821818/). PMCID: [PMC7406382](https://pmc.ncbi.nlm.nih.gov/articles/PMC7406382/). **R**; estudo primário de modelagem global.
4. **S04 — Schuchmann K, Müller V (2014).** *Autotrophy at the thermodynamic limit of life: a model for energy conservation in acetogenic bacteria.* Nature Reviews Microbiology. DOI: [10.1038/nrmicro3365](https://doi.org/10.1038/nrmicro3365). PMID: [25383604](https://pubmed.ncbi.nlm.nih.gov/25383604/). **R**; revisão mecanística.
5. **S05 — Thauer RK et al. (2008).** *Methanogenic archaea: ecologically relevant differences in energy conservation.* Nature Reviews Microbiology. DOI: [10.1038/nrmicro1931](https://doi.org/10.1038/nrmicro1931). PMID: [18587410](https://pubmed.ncbi.nlm.nih.gov/18587410/). **R**; revisão.
6. **S06 — Martens-Habbena W et al. (2009).** *Ammonia oxidation kinetics determine niche separation of nitrifying Archaea and Bacteria.* Nature. DOI: [10.1038/nature08465](https://doi.org/10.1038/nature08465). PMID: [19794413](https://pubmed.ncbi.nlm.nih.gov/19794413/). **R**; cinética experimental de SCM1.
7. **S07 — Kits KD et al. (2017).** *Kinetic analysis of a complete nitrifier reveals an oligotrophic lifestyle.* Nature. DOI: [10.1038/nature23679](https://doi.org/10.1038/nature23679). PMID: [28847001](https://pubmed.ncbi.nlm.nih.gov/28847001/). PMCID: [PMC5600814](https://pmc.ncbi.nlm.nih.gov/articles/PMC5600814/). **T**; texto, figuras e métodos de cinética/rendimento.
8. **S08 — Könneke M et al. (2014).** *Ammonia-oxidizing archaea use the most energy-efficient aerobic pathway for CO₂ fixation.* PNAS. DOI: [10.1073/pnas.1402028111](https://doi.org/10.1073/pnas.1402028111). PMID: [24843170](https://pubmed.ncbi.nlm.nih.gov/24843170/). PMCID: [PMC4050595](https://pmc.ncbi.nlm.nih.gov/articles/PMC4050595/). **T**; comparação de crescimento e tabela 4. O superlativo do título é de via aeróbia/escopo estudado, não de organismo completo.
9. **S09 — McKinlay JB, Harwood CS (2010).** *Carbon dioxide fixation as a central redox cofactor recycling mechanism in bacteria.* PNAS. DOI: [10.1073/pnas.1006175107](https://doi.org/10.1073/pnas.1006175107). PMID: [20558750](https://pubmed.ncbi.nlm.nih.gov/20558750/). PMCID: [PMC2900684](https://pmc.ncbi.nlm.nih.gov/articles/PMC2900684/). **T**; usados resultados do tipo selvagem e tabela 2, não propostas de engenharia.
10. **S10 — Alloul A et al. (2019).** *Volatile fatty acids impacting phototrophic growth kinetics of purple bacteria: Paving the way for protein production on fermented wastewater.* Water Research. DOI: [10.1016/j.watres.2018.12.025](https://doi.org/10.1016/j.watres.2018.12.025). PMID: [30665160](https://pubmed.ncbi.nlm.nih.gov/30665160/). **R**; experimento com seis culturas, inclusive comunidades.
11. **S11 — Amer A, Kim Y (2023).** *Minimizing the Lag Phase of Cupriavidus necator Growth under Autotrophic, Heterotrophic, and Mixotrophic Conditions.* Applied and Environmental Microbiology. DOI: [10.1128/aem.02007-22](https://doi.org/10.1128/aem.02007-22). PMID: [36719244](https://pubmed.ncbi.nlm.nih.gov/36719244/). PMCID: [PMC9972949](https://pmc.ncbi.nlm.nih.gov/articles/PMC9972949/). **T**; texto explicita ausência de acompanhamento completo do consumo dos substratos. Não usado como fonte de rendimento energético.
12. **S12 — Lambauer V, Kratzer R (2022).** *Lab-Scale Cultivation of Cupriavidus necator on Explosive Gas Mixtures: Carbon Dioxide Fixation into Polyhydroxybutyrate.* Bioengineering. DOI: [10.3390/bioengineering9050204](https://doi.org/10.3390/bioengineering9050204). PMID: [35621482](https://pubmed.ncbi.nlm.nih.gov/35621482/). PMCID: [PMC9138072](https://pmc.ncbi.nlm.nih.gov/articles/PMC9138072/). **R**; usados apenas título/conteúdo de biomassa natural H16; nenhum procedimento é proposto.
13. **S13 — de Vree JH et al. (2015).** *Comparison of four outdoor pilot-scale photobioreactors.* Biotechnology for Biofuels. DOI: [10.1186/s13068-015-0400-2](https://doi.org/10.1186/s13068-015-0400-2). PMID: [26689675](https://pubmed.ncbi.nlm.nih.gov/26689675/). PMCID: [PMC4683866](https://pmc.ncbi.nlm.nih.gov/articles/PMC4683866/). **T**; tabela 1 e definições da eficiência/área.
14. **S14 — Blankenship RE et al. (2011).** *Comparing photosynthetic and photovoltaic efficiencies and recognizing the potential for improvement.* Science. DOI: [10.1126/science.1200165](https://doi.org/10.1126/science.1200165). PMID: [21566184](https://pubmed.ncbi.nlm.nih.gov/21566184/). **R**; usado para fronteiras de comparação, não propostas de modificação.
15. **S15 — Logan BE et al. (2006).** *Microbial fuel cells: methodology and technology.* Environmental Science & Technology. DOI: [10.1021/es0605016](https://doi.org/10.1021/es0605016). PMID: [16999087](https://pubmed.ncbi.nlm.nih.gov/16999087/). **R**, título/autores confirmados também no Crossref; referência metodológica.
16. **S16 — Bond DR, Lovley DR (2003).** *Electricity production by Geobacter sulfurreducens attached to electrodes.* Applied and Environmental Microbiology. DOI: [10.1128/AEM.69.3.1548-1555.2003](https://doi.org/10.1128/AEM.69.3.1548-1555.2003). PMID: [12620842](https://pubmed.ncbi.nlm.nih.gov/12620842/). PMCID: [PMC150094](https://pmc.ncbi.nlm.nih.gov/articles/PMC150094/). **T**; seção “Stoichiometry and energetics”.
17. **S17 — Gregory KB, Bond DR, Lovley DR (2004).** *Graphite electrodes as electron donors for anaerobic respiration.* Environmental Microbiology. DOI: [10.1111/j.1462-2920.2004.00593.x](https://doi.org/10.1111/j.1462-2920.2004.00593.x). PMID: [15142248](https://pubmed.ncbi.nlm.nih.gov/15142248/). **R**; contraste de direção do fluxo de elétrons.
18. **S18 — Ji M et al. (2017).** *Atmospheric trace gases support primary production in Antarctic desert surface soil.* Nature. DOI: [10.1038/nature25014](https://doi.org/10.1038/nature25014). PMID: [29211716](https://pubmed.ncbi.nlm.nih.gov/29211716/). **R**; estudo comunitário, não campeão individual.
19. **S19 — Greening C et al. (2015).** *Atmospheric hydrogen scavenging: from enzymes to ecosystems.* Applied and Environmental Microbiology. DOI: [10.1128/AEM.03364-14](https://doi.org/10.1128/AEM.03364-14). PMID: [25501483](https://pubmed.ncbi.nlm.nih.gov/25501483/). PMCID: [PMC4309691](https://pmc.ncbi.nlm.nih.gov/articles/PMC4309691/). **T**; revisão, com distinção explícita de manutenção, crescimento e homólogos de hidrogenases.
20. **S20 — Orth JD, Thiele I, Palsson BØ (2010).** *What is flux balance analysis?* Nature Biotechnology. DOI: [10.1038/nbt.1614](https://doi.org/10.1038/nbt.1614). PMID: [20212490](https://pubmed.ncbi.nlm.nih.gov/20212490/). PMCID: [PMC3108565](https://pmc.ncbi.nlm.nih.gov/articles/PMC3108565/). **T**; formulação, limites de fluxo, objetivo, soluções alternativas e limitações.

## 8. Proveniência e limites de busca

- Primeiro foi chamado **ToolUniverse `find_tools`**, seguido por descoberta em `EuropePMC_search_articles`. Complementos: API REST Europe PMC com `resultType=core`, páginas HTML PMC e consulta Crossref do DOI de S15.
- As seis buscas iniciais foram recortes de 3 resultados: *Nitrospira inopinata* + oligotrophic (74 registros anunciados); *N. maritimus* + kinetics (117); *R. palustris* + energy + efficiency + yield (596); *C. necator* + hydrogen + yield (726); photosynthetic efficiency + microalgae + theoretical (416); *Geobacter* + coulombic efficiency + energy (312). **Esses totais não são estudos lidos nem estudos únicos.**
- A busca foi refinada por títulos/DOIs seminais e por termos de rendimento, com extração dirigida de tabelas. Não houve paginação exaustiva, revisão de todos os suplementos ou verificação sistemática de retratações.
- Algumas extrações XML via ToolUniverse registraram HTTP 406 no Europe PMC e 403 no OAI; houve fallback bem-sucedido ao NCBI efetch e leitura HTML PMC. A falha de uma rota não foi tratada como ausência do artigo.
- Europe PMC retornou associação ambígua do DOI de S15 a dois títulos/PMIDs; Crossref confirmou *Microbial Fuel Cells: Methodology and Technology*, autores Logan et al., páginas 5181–5192. Foi usado PMID 16999087, não o título alternativo.
- Não foi utilizado NotebookLM. Não houve execução de FBA nem uso do atlas anterior. Apenas os dois arquivos desta frente foram produzidos.

**Crédito metodológico, separado das 20 fontes de bioenergética/modelagem:** foram consultadas as skills `scientific-critical-thinking`, `cobrapy` e `paper-lookup`. Referência de Scientific Agent Skills conferida no arXiv: Kassis T, Agarwal V, He Y, Patel D, Brueckner AM (2026), *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*, [10.48550/arXiv.2609.00065](https://doi.org/10.48550/arXiv.2609.00065). Não é evidência de eficiência microbiana.
