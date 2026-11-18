# Auditoria funcional genômica da shortlist energética

**Consulta: 17/11/2026. Escopo:** representantes de módulos, não inventário exaustivo de genes, não revisão ampla de literatura e não comparação de eficiência.

## Conclusões principais

1. **Pio não é exclusivo de TIE-1.** KEGG atribui PioA/B/C tanto à CGA009 (TX73_003855/003850/003845) quanto à TIE-1 (Rpal_0817/0816/0815). Os acessos UniProt da CGA009 são A0AAE9XUC8, Q6NBT3 e Q6NBT4. Isso é evidência de anotação candidata em CGA009, não demonstração de fotoferrotrofia ou eletrotrofia nessa cepa. Os próprios KOs citam trabalhos funcionais de TIE-1: PMID 17189359 e 31690680; não os transferir entre cepas.
2. **H16 tem sistemas de H2 distintos.** HoxG/K, P31891/P31892, são subunidades da hidrogenase de membrana; HoxH/F, P22320/P22317, da hidrogenase redutora de NAD. Não contar os quatro como quatro fontes de energia.
3. **Duas CODHs não equivalentes:** Q2RUG7 em R. rubrum é anotada como família Ni-CODH, EC 1.2.7.4, com KO K00198 anaeróbico; Q0KEI3 em H16 é apenas 'aerobic-type carbon monoxide dehydrogenase homolog', predita, com domínio da família aldeído/xantina-oxidoredutase. Não demonstra CO como substrato fisiológico, centro Mo/Cu funcional nem carboxidotrofia de H16.
4. **Erro de cepa concreto:** Q6QQP7, frequentemente recuperada por 'malonyl-CoA reductase', tem STRAIN=OK-70-fl nas referências, não J-10-fl. Os candidatos nativos J-10-fl são A9WIU3/Caur_2614 e A9WEU4/Caur_0613, com K14468/K14469 no KEGG, apesar de nomes UniProt genéricos/divergentes.
5. **A evidência energética molecular mais direta desta seleção é localizada:** A9A1Y1 de SCM1 tem ECO:0000269, PMID 24843170, para a ligase 4HB-CoA ADP-forming do ciclo HP/HB. Sua reação consome ATP e forma ADP+Pi. Isto informa o custo de uma etapa, não uma eficiência global em %, rendimento de biomassa ou superioridade geral do organismo.

## Identidade das cepas e proteomas

| Alvo | Identidade registrada | Taxon ID | Proteoma | Montagem vinculada |
|---|---|---:|---|---|
| CGA009 | Rhodopseudomonas palustris (strain ATCC BAA-98 / CGA009) | 258594 | [UP000001426](https://rest.uniprot.org/proteomes/UP000001426.json) | GCA_000195775.2 |
| TIE-1 | Rhodopseudomonas palustris (strain TIE-1) | 395960 | [UP000001725](https://rest.uniprot.org/proteomes/UP000001725.json) | GCA_000020445.1 |
| H16 | Cupriavidus necator (strain ATCC 17699 / DSM 428 / KCTC 22496 / NCIMB 10442 / H16 / Stanier 337) | 381666 | [UP000008210](https://rest.uniprot.org/proteomes/UP000008210.json) | GCA_000009285.2 |
| ATCC11170 | Rhodospirillum rubrum (strain ATCC 11170 / ATH 1.1.1 / DSM 467 / LMG 4362 / NCIMB 8255 / S1) | 269796 | [UP000001929](https://rest.uniprot.org/proteomes/UP000001929.json) | GCA_000013085.1 |
| DSM2944 | Paracoccus pantotrophus; strain=ATCC 35512 / DSM 2944 / CIP 106514 / LMD 82.5 / NBRC 102493 / NCCB 82005 / GB17 | 82367 | [UP000273626](https://rest.uniprot.org/proteomes/UP000273626.json) | GCA_003633525.1 |
| PD1222 | Paracoccus denitrificans (strain Pd 1222) | 318586 | [UP000000361](https://rest.uniprot.org/proteomes/UP000000361.json) | GCA_000203895.1 |
| J-10-fl | Chloroflexus aurantiacus (strain ATCC 29366 / DSM 635 / J-10-fl) | 324602 | [UP000002008](https://rest.uniprot.org/proteomes/UP000002008.json) | GCA_000018865.1 |
| SCM1 | Nitrosopumilus maritimus (strain SCM1) | 436308 | [UP000000792](https://rest.uniprot.org/proteomes/UP000000792.json) | GCA_000018465.1 |

- DSM 2944 e GB17 são aliases no registro UP000273626; o taxon 82367 é de espécie, portanto o filtro foi **proteome:UP000273626**, não a espécie inteira. PD1222 é outra espécie/cepa, taxon 318586 ('Pd 1222'), e foi auditada separadamente. Não transferir os ensaios Sox entre elas.
- TIE-1 é taxon **395960**. O proteoma está cadastrado com 5.246 proteínas, mas a consulta funcional em UniProtKB recuperou somente cinco entradas, todas revisadas. Pio e PufL foram recuperados via KEGG Genes, não declarados ausentes. Cadastro de proteoma e disponibilidade atual no UniProtKB não são sinônimos.
- Entradas revisadas históricas de CGA009 e TIE-1 podem não ter ligação ao proteoma corrente. O JSON distingue ligação explícita de apenas identidade taxonômica exata. Tags RPA versus TX73 não foram alinhadas entre versões; entradas históricas e reanotações não foram contadas como duplicações gênicas.

## Leitura da evidência

**48 representantes**, mais **dois registros de cautela/exclusão**. Os níveis abaixo são de **anotação funcional**, não a escala de evidência de perda/pseudogenização: nenhuma perda biológica foi inferida.

- **A1:** entrada não revisada; nome importado ou inferência automática.
- **A2:** entrada revisada; sem ECO:0000269 nos campos FUNCTION/CATALYTIC ACTIVITY recuperados. Isso não nega experimentos existentes; limita esta extração.
- **A3:** ECO:0000269 nesses campos, com fonte experimental identificável, limitada à proteína/reação correspondente.
- **K1:** gene/KO do KEGG na cepa correta; o KEGG não fornece o estado Swiss-Prot/TrEMBL nem códigos ECO equivalentes.

**R/U:** revisado/não revisado. **PE:** campo 'protein existence' da entrada; PE1 é evidência em nível proteico, não comprovação de todo o comentário funcional. Nas colunas ECO, N=nome/descrição, F=função/catálise; '—' significa código não fornecido nesse campo, não ausência de função ou inexistência de experimento.

Códigos encontrados: ECO:0000313=informação importada em anotação automática; ECO:0000256=inferência por modelo/sequência em anotação automática; ECO:0000255=modelo/sequência em anotação manual; ECO:0000250=similaridade de sequência; ECO:0000305=inferência de curador; ECO:0000303=afirmação de autor; ECO:0000269=evidência experimental em anotação manual. O nome de Nmar_1500 traz **ECO:0008006, Google/ProtNLM**, não ensaio bioquímico.

## Representantes confirmados nas bases

### CGA009

| Acesso | Representante / módulo | Gene/tag na entrada | R/U; PE; nível | ECO N / F ou KO |
|---|---|---|---|---|
| [A0AAE9XUC8](https://rest.uniprot.org/uniprotkb/A0AAE9XUC8.json) | PioA-like; Transferência extracelular candidata | TX73_003855 | U; PE4; A1 | N: ECO:0000313; F: — |
| [A0AAE9XYU2](https://rest.uniprot.org/uniprotkb/A0AAE9XYU2.json) | PuhA; Fotossistema | puhA, TX73_007955 | U; PE4; A1 | N: ECO:0000313; F: — |
| [O83005](https://rest.uniprot.org/uniprotkb/O83005.json) | PufL; Fotossistema | pufL, RPA1527 | R; PE1; A2 | N: —; F: — |
| [Q6N9J0](https://rest.uniprot.org/uniprotkb/Q6N9J0.json) | CbbL / RuBisCO forma I; Fixação de carbono | cbbL, RPA1559 | R; PE3; A2 | N: ECO:0000255; F: ECO:0000255 |
| [Q6N9K0](https://rest.uniprot.org/uniprotkb/Q6N9K0.json) | BchL; Bacterioclorofila | bchL, RPA1545 | R; PE3; A2 | N: ECO:0000255; F: ECO:0000255 |
| [Q6NB67](https://rest.uniprot.org/uniprotkb/Q6NB67.json) | Uptake hydrogenase grande; Oxidação de H2 candidata | TX73_004960 | U; PE3; A1 | N: ECO:0000256; F: ECO:0000256 |
| [Q6NB68](https://rest.uniprot.org/uniprotkb/Q6NB68.json) | Uptake hydrogenase pequena; Oxidação de H2 candidata | TX73_004955 | U; PE3; A1 | N: ECO:0000256; F: ECO:0000256 |
| [Q6NBT3](https://rest.uniprot.org/uniprotkb/Q6NBT3.json) | PioB-like; Transferência extracelular candidata | TX73_003850 | U; PE4; A1 | N: ECO:0000313; F: — |
| [Q6NBT4](https://rest.uniprot.org/uniprotkb/Q6NBT4.json) | PioC-like; Transferência extracelular candidata | TX73_003845 | U; PE3; A1 | N: ECO:0000313; F: ECO:0000256 |

### TIE-1

| Acesso | Representante / módulo | Gene/tag na entrada | R/U; PE; nível | ECO N / F ou KO |
|---|---|---|---|---|
| [B3Q7C4](https://rest.uniprot.org/uniprotkb/B3Q7C4.json) | BchN; Bacterioclorofila | bchN, Rpal_1730 | R; PE3; A2 | N: ECO:0000255; F: ECO:0000255 |
| [B3Q7E1](https://rest.uniprot.org/uniprotkb/B3Q7E1.json) | CbbL / RuBisCO forma I; Fixação de carbono | cbbL, Rpal_1747 | R; PE3; A2 | N: ECO:0000255; F: ECO:0000255 |
| [rpt:Rpal_0815](https://rest.kegg.jp/get/rpt:Rpal_0815) | PioC / K28461; Transferência extracelular | Rpal_0815 | n/a; PEn/a; K1 | K28461  high potential iron sulfer protein |
| [rpt:Rpal_0816](https://rest.kegg.jp/get/rpt:Rpal_0816) | PioB / K28460; Transferência extracelular | Rpal_0816 | n/a; PEn/a; K1 | K28460  decaheme-associated outer membrane protein |
| [rpt:Rpal_0817](https://rest.kegg.jp/get/rpt:Rpal_0817) | PioA / K28459; Transferência extracelular | Rpal_0817 | n/a; PEn/a; K1 | K28459  cytochrome c-type protein PioA |
| [rpt:Rpal_1154](https://rest.kegg.jp/get/rpt:Rpal_1154) | Hidrogenase pequena / K06282; Oxidação de H2 candidata | Rpal_1154 | n/a; PEn/a; K1 | K06282  hydrogenase small subunit [EC:1.12.99.6] |
| [rpt:Rpal_1155](https://rest.kegg.jp/get/rpt:Rpal_1155) | Hidrogenase grande / K06281; Oxidação de H2 candidata | Rpal_1155 | n/a; PEn/a; K1 | K06281  hydrogenase large subunit [EC:1.12.99.6] |
| [rpt:Rpal_1716](https://rest.kegg.jp/get/rpt:Rpal_1716) | PufL / K08928; Fotossistema | Rpal_1716 | n/a; PEn/a; K1 | K08928  photosynthetic reaction center L subunit |

### H16

| Acesso | Representante / módulo | Gene/tag na entrada | R/U; PE; nível | ECO N / F ou KO |
|---|---|---|---|---|
| [P22317](https://rest.uniprot.org/uniprotkb/P22317.json) | HoxF / hidrogenase solúvel; Redução de NAD por H2 | hoxF, PHG088 | R; PE1; A2 | N: —; F: — |
| [P22320](https://rest.uniprot.org/uniprotkb/P22320.json) | HoxH / hidrogenase solúvel; Redução de NAD por H2 | hoxH, PHG091 | R; PE1; A2 | N: —; F: — |
| [P31891](https://rest.uniprot.org/uniprotkb/P31891.json) | HoxG / hidrogenase de membrana; Oxidação de H2 | hoxG, PHG002 | R; PE1; A2 | N: —; F: — |
| [P31892](https://rest.uniprot.org/uniprotkb/P31892.json) | HoxK / hidrogenase de membrana; Oxidação de H2 | hoxK, PHG001 | R; PE1; A2 | N: —; F: — |
| [Q0K1E0](https://rest.uniprot.org/uniprotkb/Q0K1E0.json) | CbbL1; Fixação de carbono | cbbL1, H16_B1395 | R; PE2; A2 | N: —; F: ECO:0000250 |
| [Q0KEI3](https://rest.uniprot.org/uniprotkb/Q0KEI3.json) | CoxL1-like; CODH aeróbica apenas homóloga | coxL1, H16_A0438 | U; PE4; A1 | N: ECO:0000313; F: — |

### ATCC11170

| Acesso | Representante / módulo | Gene/tag na entrada | R/U; PE; nível | ECO N / F ou KO |
|---|---|---|---|---|
| [Q2RQ25](https://rest.uniprot.org/uniprotkb/Q2RQ25.json) | PufL-like; Fotossistema | Rru_A2975 | U; PE1; A1 | N: ECO:0000256; F: ECO:0000256 |
| [Q2RRP5](https://rest.uniprot.org/uniprotkb/Q2RRP5.json) | CbbM / RuBisCO forma II; Fixação de carbono | cbbM, Rru_A2400 | R; PE1; A2 | N: ECO:0000255; F: ECO:0000255 |
| [Q2RUG7](https://rest.uniprot.org/uniprotkb/Q2RUG7.json) | CooS-like / Ni-CODH; CODH anaeróbica | Rru_A1427 | U; PE3; A1 | N: ECO:0000256; F: ECO:0000256 |
| [Q2RV82](https://rest.uniprot.org/uniprotkb/Q2RV82.json) | Uptake hydrogenase grande; Oxidação de H2 candidata | Rru_A1162 | U; PE3; A1 | N: ECO:0000256; F: ECO:0000256 |
| [Q2RWS1](https://rest.uniprot.org/uniprotkb/Q2RWS1.json) | BchL; Bacterioclorofila | bchL, Rru_A0620 | R; PE3; A2 | N: ECO:0000255; F: ECO:0000255 |
| [Q2RWS4](https://rest.uniprot.org/uniprotkb/Q2RWS4.json) | PuhA-like; Fotossistema | Rru_A0617 | U; PE1; A1 | N: ECO:0000313; F: — |

### DSM2944

| Acesso | Representante / módulo | Gene/tag na entrada | R/U; PE; nível | ECO N / F ou KO |
|---|---|---|---|---|
| [A0A1I5INT1](https://rest.uniprot.org/uniprotkb/A0A1I5INT1.json) | SoxY; Sox | BDE18_3405 | U; PE4; A1 | N: ECO:0000313; F: — |
| [A0A1I5INZ2](https://rest.uniprot.org/uniprotkb/A0A1I5INZ2.json) | SoxD; Sox | BDE18_3410 | U; PE4; A1 | N: ECO:0000313; F: — |
| [A0A1I5IP77](https://rest.uniprot.org/uniprotkb/A0A1I5IP77.json) | SoxA; Sox | BDE18_3407 | U; PE3; A1 | N: ECO:0000256; F: ECO:0000256 |
| [A0A1I5IPK4](https://rest.uniprot.org/uniprotkb/A0A1I5IPK4.json) | SoxC; Sox | BDE18_3409 | U; PE4; A1 | N: ECO:0000313; F: — |
| [A0A1I5J2X8](https://rest.uniprot.org/uniprotkb/A0A1I5J2X8.json) | CbbL; Fixação de carbono | cbbL, BDE18_1315 | U; PE3; A1 | N: ECO:0000256; F: ECO:0000256 |
| [A0AAE6TRU3](https://rest.uniprot.org/uniprotkb/A0AAE6TRU3.json) | Uptake hydrogenase grande; Oxidação de H2 candidata | BDE18_2448 | U; PE3; A1 | N: ECO:0000256; F: ECO:0000256 |
| [A0AAE6TS95](https://rest.uniprot.org/uniprotkb/A0AAE6TS95.json) | SoxB; Sox | BDE18_3408 | U; PE3; A1 | N: ECO:0000313; F: — |

### PD1222

| Acesso | Representante / módulo | Gene/tag na entrada | R/U; PE; nível | ECO N / F ou KO |
|---|---|---|---|---|
| [A1B2Q2](https://rest.uniprot.org/uniprotkb/A1B2Q2.json) | CbbL; Fixação de carbono | cbbL, Pden_1699 | R; PE3; A2 | N: ECO:0000255; F: ECO:0000255 |
| [A1B6N5](https://rest.uniprot.org/uniprotkb/A1B6N5.json) | Uptake hydrogenase grande; Oxidação de H2 candidata | Pden_3098 | U; PE3; A1 | N: ECO:0000256; F: ECO:0000256 |
| [A1B9M4](https://rest.uniprot.org/uniprotkb/A1B9M4.json) | SoxB; Sox | Pden_4154 | U; PE3; A1 | N: ECO:0000313; F: — |
| [A1B9T9](https://rest.uniprot.org/uniprotkb/A1B9T9.json) | NosZ; Desnitrificação, apenas uma etapa | nosZ, Pden_4219 | U; PE3; A1 | N: ECO:0000256; F: ECO:0000256 |

### J-10-fl

| Acesso | Representante / módulo | Gene/tag na entrada | R/U; PE; nível | ECO N / F ou KO |
|---|---|---|---|---|
| [A9WEU4](https://rest.uniprot.org/uniprotkb/A9WEU4.json) | Pcs segundo KEGG K14469; 3HP bi-cycle candidato | Caur_0613 | U; PE3; A1 | N: ECO:0000256; F: ECO:0000256 |
| [A9WIU3](https://rest.uniprot.org/uniprotkb/A9WIU3.json) | Mcr segundo KEGG K14468; 3HP bi-cycle candidato | Caur_2614 | U; PE1; A1 | N: ECO:0000313; F: — |
| [P11695](https://rest.uniprot.org/uniprotkb/P11695.json) | PufL; Fotossistema | pufL, Caur_1052 | R; PE1; A2 | N: —; F: — |
| [Q9F6X3](https://rest.uniprot.org/uniprotkb/Q9F6X3.json) | BchL; Bacterioclorofila | bchL, Caur_2554 | R; PE3; A2 | N: ECO:0000255; F: ECO:0000255 |

### SCM1

| Acesso | Representante / módulo | Gene/tag na entrada | R/U; PE; nível | ECO N / F ou KO |
|---|---|---|---|---|
| [A9A1Y1](https://rest.uniprot.org/uniprotkb/A9A1Y1.json) | 4HB-CoA ligase ADP-forming; HP/HB, etapa experimental | Nmar_0206 | R; PE1; A3 | N: ECO:0000269, ECO:0000303, ECO:0000305; F: ECO:0000269 |
| [A9A4U2](https://rest.uniprot.org/uniprotkb/A9A4U2.json) | AmoA-like; Oxidação de amônia candidata | Nmar_1500 | U; PE4; A1 | N: ECO:0008006; F: — |
| [A9A4U4](https://rest.uniprot.org/uniprotkb/A9A4U4.json) | AmoC-like; Oxidação de amônia candidata | Nmar_1502 | U; PE4; A1 | N: ECO:0000313; F: — |
| [A9A4U5](https://rest.uniprot.org/uniprotkb/A9A4U5.json) | AmoB-like; Oxidação de amônia candidata | Nmar_1503 | U; PE4; A1 | N: ECO:0008006; F: — |

## Auditoria por módulo: afirmação permitida e extrapolação proibida

### Luz: Puf/Puh e Bch
Há representantes de centros de reação anoxigênicos e de biossíntese de bacterioclorofila nas bactérias fototróficas selecionadas. PufL/PuhA são componentes, e BchL/BchN também: nenhum isoladamente representa fotossistema completo, quantidade de pigmento, rendimento quântico ou crescimento à luz. Não exigir a mesma arquitetura Puh de bactérias púrpuras em Chloroflexus. A associação KEGG a um módulo/pathway é mapeamento de membro, não certificado de completude. Não foi inferida fotossíntese oxigênica.

### H2: distinguir consumo, produção e regulação
Os acessos de uptake selecionados são candidatos a consumo de H2. Não confundir hidrogenases com nitrogenase produtora de H2, nem fatores Hyp de maturação com enzimas de obtenção de energia. Em CGA009, a anotação estrutural Q6NB67/Q6NB68 não resolve integridade/regulação de hupUV, promotores, expressão ou crescimento usando H2; nenhuma lesão foi mapeada nesta auditoria. Em TIE-1, os genes KEGG Rpal_1154/1155 constituem evidência anotacional, não uma medição de atividade. H16 contém os dois sistemas distintos selecionados; não se mediu acoplamento, tolerância a O2 ou fluxo.

**Cautela mesmo com Swiss-Prot:** P31891/P31892 possuem comentário genérico sobre reciclagem de H2 de nitrogenase, sem ECO nesse campo. Não usar esse comentário para afirmar diazotrofia da cepa H16. Revisado e PE1 não tornam toda frase uma medição específica da cepa.

### Fixação de carbono
CbbL forma I e CbbM forma II têm representantes, mas RuBisCO isolada não certifica o ciclo CBB completo, nem fixa seu fluxo ou prova autotrofia obrigatória. A auditoria de Chloroflexus aponta Mcr/Pcs do **3HP bi-cycle**, não exige CBB por analogia com bactérias púrpuras. O ciclo HP/HB de SCM1 também não deve ser substituído por um modelo CBB genérico. São módulos de assimilação de carbono, não fontes adicionais de energia por si mesmos.

### Sox e aceptores
Em DSM2944/GB17, SoxA/Y/B/C/D foram recuperados no proteoma correto; SoxB de PD1222 foi recuperada separadamente. São evidências candidatas para oxidação de compostos reduzidos de enxofre. Não equivalem a usar sulfato como fonte energética, e assimilação de sulfato não é Sox. Não se auditou completude de todos os pares SoxAX/SoxYZ ou sua maturação. NosZ A1B9T9 de PD1222 representa somente redução de N2O; não prova sozinho a sequência completa nitrato→N2, desnitrificação aeróbica ou taxas.

### CODH
R. rubrum Q2RUG7 tem inferência automática para família Ni-CODH, reação com ferredoxina, GO:0043885 com IEA, e KEGG K00198. H16 Q0KEI3 tem somente nome de homólogo tipo aeróbico importado do genoma; seus domínios são compartilhados com oxidoredutases de aldeído/xantina. A classificação bioquímica de CODH aeróbica Mo/Cu versus anaeróbica Ni/Fe-S não autoriza substituir uma pela outra, afirmar reversibilidade fisiológica, completar Wood–Ljungdahl, ou atribuir crescimento com CO a toda cepa com um homólogo. Não se validaram aqui motivos diagnósticos/tipo I versus tipo II, sítio catalítico, maturação ou especificidade de substrato de CoxL de H16.

### Pio, ferro e eletrotrofia
PioA é um citocromo multi-heme, PioB uma proteína de membrana externa associada e PioC uma HiPIP; **PioC não é uma oxidase de ferro isolada**. TIE-1 tem suporte de anotação KO e referências funcionais específicas nos KOs. CGA009 tem homólogos candidatos no proteoma atual. Não converter homologia Pio em fenótipo universal nem declarar Pio ausente da CGA009. Oxidação de Fe(II), transferência extracelular e entrada de elétrons de eletrodo não são três fontes independentes automaticamente demonstradas. Potencial do eletrodo, mediadores, luz, expressão e outros componentes precisariam de dados próprios; nada aqui mede eficiência coulômbica.

### Amo em SCM1
Foram retidas as três entradas genômicas Nmar_1500/1503/1502, com nomes Amo/AmoB/AmoC candidatos. Nmar_1500 tem nome inferido por ProtNLM e domínio InterPro AmoA_arc. Isso não é evidência de atividade enzimática direta; não tratar similaridade à família amônia/metano monooxigenase como prova de metanotrofia. O produto energético é uma hipótese específica do sistema de oxidação de amônia, não de versatilidade extrema com muitos substratos. A9A1Y1 é uma etapa de fixação de carbono com anotação experimental separada.

## Dois registros deliberadamente não usados para confirmar módulos

| Acesso | Problema verificável | Consequência |
|---|---|---|
| [Q6QQP7](https://rest.uniprot.org/uniprotkb/Q6QQP7.json) | Mcr de OK-70-fl, taxon de espécie 1108; referências STRAIN=OK-70-fl | Não é evidência nativa da cepa J-10-fl; usar A9WIU3 com KO K14468 |
| [Q0K5E4](https://rest.uniprot.org/uniprotkb/Q0K5E4.json) | H16: rótulo não revisado 'putative ammonia monooxygenase', gene abrB, sem demonstração do complexo AmoABC | Não promove H16 a nitrificador; não se conclui que a função rotulada seja biologicamente impossível |

## O que estes dados não podem dizer

A presença de uma entrada ou KO não demonstra ORF íntegra em todas as subculturas, via completa, disponibilidade de cofatores, expressão, fluxo, crescimento, rendimento, P/O, ATP por elétron, rendimento quântico, eficiência coulômbica ou porcentagem de conversão de energia. Não foi feita contagem de 'fontes de energia' a partir do número de genes. Não houve testes de ausência com controle positivo, genomas/leituras não foram reanalisados, e nenhum resultado zero foi tratado como ausência biológica. Não se inferiram perdas, ancestralidade, ramos, ordem de mutações ou trajetória evolutiva real. Esta é uma auditoria funcional de registros existentes, não uma reconstrução filogenética nem um desenho experimental.

## Métodos e proveniência

- Descoberta de ferramentas: ToolUniverse **find_tools**, consulta de UniProtRef_search_proteomes para R. palustris e documentação de UniProt_search; extração final pelo REST oficial público.
- Endpoints: `https://rest.uniprot.org/proteomes/{UPID}.json`, `https://rest.uniprot.org/uniprotkb/search`; KEGG `/get/{IDs}` em lotes de até dez. URLs, parâmetros, horários, release, totais e status HTTP das requisições do script estão em `anotacoes.json:request_log`.
- Consultas iniciais filtradas pela taxonomia **exata da cepa**; exceção DSM2944, que exigiu filtro por proteoma. Nomes/genéricos, tags e evidências foram verificados localmente. A busca exploratória de Mcr no taxon de espécie recuperou Q6QQP7 de outra cepa, explicitamente excluída.
- Uma busca `gene:` com prefixos BDE18/Pden devolveu zero; buscas por acesso/tag textual e descrições recuperaram SoxB/C/D. Essa falha de recuperação não foi convertida em ausência dos genes.
- O KEGG separou corretamente PioB/C (Rpal_0816/0815) de transportadores de sulfato adjacentes Rpal_0818/0819; não se usou contiguidade numérica como identidade funcional. Rpal_1715 é antena, não PufL; o representante validado é Rpal_1716.
- ECO foi extraído por campo e fonte, não inventado a partir de 'reviewed'. Foram preservados os comentários de função/catálise/cofator, versões de entrada, referências, RefSeq/EMBL/KEGG/InterPro/GO e ligação ao proteoma. GO-IEA não foi tratado como ensaio direto.
- Não foram lidos ou alterados `scripts/bioenergetics_efficiency_analysis.py` nem `runs/bioenergetics_efficiency_atlas.json`; nenhum número fixo desses arquivos foi reutilizado. Sem genomas grandes, NIM, NotebookLM ou previsões de eficiência.

**Reconciliação:** 209 hits nas oito buscas iniciais; 223 acessos candidatos únicos após buscas suplementares; 42 representantes UniProt e 6 KEGG retidos; 2 cautelas separadas. 24 requisições HTTP no script. Release UniProt registrada: 2026_04. Todos os conjuntos de busca do script couberam na primeira página, sem Link next. Isto demonstra completude das consultas registradas, **não** completude de cada via no genoma.

### Fontes-chave, com o escopo correto

- [KO K28459](https://rest.kegg.jp/get/ko:K28459), [K28460](https://rest.kegg.jp/get/ko:K28460), [K28461](https://rest.kegg.jp/get/ko:K28461): identidades Pio e referências de TIE-1; listam também genes CGA009.
- [A9A1Y1](https://rest.uniprot.org/uniprotkb/A9A1Y1.json), [PMID 24843170](https://pubmed.ncbi.nlm.nih.gov/24843170/): anotação experimental da ligase ADP-forming de SCM1.
- [Caur_2614](https://rest.kegg.jp/get/cau:Caur_2614) e [Caur_0613](https://rest.kegg.jp/get/cau:Caur_0613): Mcr/Pcs candidatos na linhagem genômica J-10-fl.
- [Q2RUG7](https://rest.uniprot.org/uniprotkb/Q2RUG7.json) e [Rru_A1427](https://rest.kegg.jp/get/rru:Rru_A1427): Ni-CODH anaeróbica candidata de ATCC11170; [Q0KEI3](https://rest.uniprot.org/uniprotkb/Q0KEI3.json): homólogo tipo aeróbico de H16, com incerteza de substrato mantida.
