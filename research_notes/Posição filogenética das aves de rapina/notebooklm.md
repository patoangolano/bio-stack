# NotebookLM: posição filogenética das aves de rapina

Braço NotebookLM do levantamento sobre a posição das aves de rapina (Cathartiformes, Accipitriformes, Falconiformes, Strigiformes) entre si e dentro de Telluraves/Neoaves, incluindo fósseis.

- **Caderno:** "Aves de rapina — filogenia, fósseis e Telluraves (bio-stack 2026-10-17)"
- **ID do caderno:** `1d1f534b-6c17-4613-ad6e-816739390010`
- **URL:** https://notebook.google.com/notebook/1d1f534b-6c17-4613-ad6e-816739390010
- **Conta ativa:** perfil `default` (limatheus458@gmail.com), único perfil salvo. Nenhuma chamada devolveu `active_profile_note` ou `user_notice`.
- **Data:** 2026-10-17
- **Pesquisa profunda:** `research_start` (mode=deep), task `0c00e563-048b-4473-aabb-249d6dc2fec5`. Achou 42 itens e eu importei 11 acadêmicos, dos quais 10 entraram com texto.
- **Consultas:** as 10 perguntas foram feitas em conversas separadas e **restritas às 41 fontes com conteúdo** (`source_ids`), para que páginas de erro não entrassem nas respostas. O Briefing Doc usou o mesmo conjunto.

## Como o texto chegou ao caderno

O NotebookLM não consegue importar páginas do PMC (`pmc.ncbi.nlm.nih.gov` devolve um reCAPTCHA). Também falham OUP, Royal Society e PNAS-via-doi (desafio Cloudflare "Just a moment..."), BioOne (bloqueio Incapsula) e link.springer.com por URL (HTTP 406). Por isso:

- **Artigos com PMCID:** baixei o JATS do **Europe PMC** (`https://www.ebi.ac.uk/europepmc/webservices/rest/<PMCID>/fullTextXML`), converti para `.txt` com `C:/Users/Manec/AppData/Local/Temp/claude/C--Users-Manec-bio-stack/81d34c2d-4fea-4683-acee-96da124ac70d/scratchpad/scripts/epmc2txt.py` (temporário; o mesmo diretório tem `qa2md.py`) (texto, legendas, sem lista de referências) e subi como arquivo. Jarvis 2014 não está no Europe PMC (HTTP 500) e veio do `efetch` do NCBI (`db=pmc&id=4405904`).
- **Springer de acesso aberto:** baixei o PDF de `https://link.springer.com/content/pdf/<doi>.pdf` com curl e subi como arquivo.
- **Arquivos de origem:** `C:\Users\Manec\AppData\Local\Temp\claude\C--Users-Manec-bio-stack\81d34c2d-4fea-4683-acee-96da124ac70d\scratchpad\epmc\` e `...\scratchpad\springer\`.

## Fontes

O caderno tem 59 entradas: **41 com conteúdo** (40 obras distintas, já que a *Phylogenetic definitions* entrou em duas versões), **1 duplicata** (Hackett 2008, segunda cópia) e **17 entradas que falharam** (páginas de captcha/erro ou conteúdo vazio). Deixei as entradas que falharam no caderno, porque `source_delete` exige confirmação do usuário. Elas podem ser apagadas pelos IDs listados abaixo.

### Fontes importadas com conteúdo (usadas nas consultas)

| # | Obra | URL / origem | Importada | Obs. |
|---|---|---|---|---|
| 1 | Hackett et al. 2008, *Science* 320:1763 — A phylogenomic study of birds reveals their evolutionary history | https://earlybird.biology.ufl.edu/wp-content/uploads/sites/557/Science-2008-Hackett-1763-8.pdf | sim (`cad68562`) | o PDF traz também outro artigo da mesma edição, sobre deslocamento altitudinal de plantas |
| 1b | idem, cópia | https://pages.stat.wisc.edu/~larget/Genetics629/Science-2008-Hackett-1763-8.pdf | sim (`dbf6865d`) | duplicata, fora das consultas |
| 2 | Suh et al. 2011, *Nat Commun* 2:443 — Mesozoic retroposons… (Eufalconimorphae/Psittacopasserae) | PMC3265382 → Europe PMC fullTextXML | sim, como arquivo (`7cfbfef1`) | a URL do PMC falhou (reCAPTCHA) |
| 3 | Suh et al. 2015, *PLoS Biol* — ILS across the neoavian radiation | PMC4540587 → Europe PMC | sim, como arquivo (`2d01758b`) | a URL do PMC falhou (reCAPTCHA) |
| 4 | Jarvis et al. 2014, *Science* 346:1320 — Whole-genome analyses… | PMC4405904 → NCBI efetch | sim, como arquivo (`c92cd31c`) | a URL do PMC falhou; Europe PMC deu HTTP 500 |
| 5 | Prum et al. 2015, *Nature* 526:569 | https://prumlab.yale.edu/sites/default/files/prum_et_al_2015.pdf | sim (`74f2801f`) | |
| 6 | Reddy et al. 2017, *Syst Biol* 66:857 | https://repository.si.edu/handle/10088/32489 | **só resumo** (`bc2e5a15`) | https://doi.org/10.1093/sysbio/syx041 falhou (Cloudflare) |
| 7 | Kimball et al. 2019, *Diversity* 11:109 — A phylogenomic supertree of birds | https://doi.org/10.3390/d11070109 | sim (`6ed04532`) | |
| 8 | Kuhl et al. 2021, *MBE* (doi 10.1093/molbev/msaa191) — 3′-UTRs | PMC7783168 → Europe PMC | sim, como arquivo (`8527a287`) | a URL do PMC falhou (reCAPTCHA) |
| 9 | Braun & Kimball 2021, *Birds* 2:1 — Data Types and the Phylogeny of Neoaves | https://www.mdpi.com/2673-6004/2/1/1 | sim (`b0e43963`) | |
| 10 | Wu et al. 2024, *PNAS* 121:e2319696121 — Genomes, fossils… Late Cretaceous | PMC10895254 → Europe PMC | sim, como arquivo (`cda64d61`) | a URL do PMC falhou (reCAPTCHA) |
| 11 | Stiller et al. 2024, *Nature* — family-level genomes (B10K) | PMC11111414 → Europe PMC | sim, como arquivo (`e9c5ae30`) | |
| 12 | Catanach, Halley & Pirro 2024, *Biol J Linn Soc* 144(2) — Enigmas no longer (Accipiter, UCE) | https://par.nsf.gov/servlets/purl/10511513 | sim (`1d4669b1`) | |
| 13 | Fuchs, Johnson & Mindell 2015, *MPE* 82:166 — Falconidae | https://doi.org/10.1016/j.ympev.2014.08.010 | **só resumo/prévia** ScienceDirect (`f82bdd8e`) | sem texto completo |
| 14 | Salter et al. 2020, *Auk* 137:ukz070 — paraphyly in Strigidae | https://par.nsf.gov/servlets/purl/10183927 | sim (`1e8f2a66`) | |
| 15 | McClure et al. 2019, *J Raptor Res* 53:419 — Defining raptors | https://www.hawkmountain.org/download/?id=5956 | **só página de citação** (`7823ecf1`) | BioOne falhou (Incapsula) |
| 16 | Ericson 2012, *J Biogeogr* 39:813 | https://doi.org/10.1111/j.1365-2699.2011.02650.x | sim, página Wiley (`254c47a3`) | o resumo do NotebookLM indica que traz resultados |
| 17 | Claramunt & Cracraft 2015, *Sci Adv* (doi 10.1126/sciadv.1501005) — A new time tree… | PMC4730849 → Europe PMC | sim, como arquivo (`fcdd26dd`) | |
| 18 | Mayr & Kitchener 2022, *JVP* 41(6):e2083515 — Masillaraptoridae = stem Falconiformes (*Danielsraptor*) | https://doi.org/10.1080/02724634.2021.2083515 | sim, texto completo T&F (`53dae551`) | |
| 19 | Mayr, Gingerich & Smith 2020, *JVP* 40:e1769116 — *Primoptynx poliotauros* | https://library.naturalsciences.be/pdfs-open-access/2020/Mayr%20et%20al%202020%20-Primoptynx%20from%20early%20Eocene%20of%20Wyoming.pdf/@@download/file/Mayr et al 2020 (Primoptynx from early Eocene of Wyoming).pdf | sim (`2b83db39`) | |
| 20 | Mayr 2021, *PalZ* 95:337 — *Tynskya waltonensis* (Messelasturidae) | PDF de https://link.springer.com/content/pdf/10.1007/s12542-020-00541-8.pdf | sim, como arquivo (`a681b573`) | https://doi.org/10.1007/s12542-020-00541-8 falhou (HTTP 406) |
| 21 | Mayr & Hurum 2020, *Sci Nat* — acipitrídeo pequeno do Oligoceno inferior (*Aviraptor*) | PDF de https://link.springer.com/content/pdf/10.1007/s00114-020-01703-z.pdf | sim, como arquivo (`32bd8904`) | https://doi.org/10.1007/s00114-020-01703-z falhou (HTTP 406) |
| 22 | Degrange et al. 2010, *PLoS ONE* (doi 10.1371/journal.pone.0011856) — *Andalgalornis* (Phorusrhacidae) | PMC2923598 → Europe PMC | sim, como arquivo (`b9c81379`) | |
| 23 | Ksepka, Stidham & Williamson 2017, *PNAS* 114(30) — *Tsidiiyazhi* | https://doi.org/10.1073/pnas.1700188114 | sim (`b84ec51f`) | |
| 24 | Livezey & Zusi 2007, *Zool J Linn Soc* — Higher-order phylogeny of modern birds… II (2954 caracteres morfológicos) | PMC2517308 → Europe PMC | sim, como arquivo (`341dae82`) | |
| 25 | Wang et al. 2014, *GBE* (doi 10.1093/gbe/evu016) — Phylogenetic position of avian nocturnal and diurnal raptors (mitogenomas) | PMC3942033 → Europe PMC | sim, como arquivo (`3cc767df`) | |
| 26 | Chen et al. 2025, *Integr Org Biol* 7:obaf029 — matriz anatômica (cintura peitoral/membro anterior) | PMC12409627 → Europe PMC | sim, como arquivo (`6c89ed1d`) | 2025 |
| 27 | Buthasane et al. 2026, *ZooKeys* 1281:105 — genoma de *Aquila nipalensis* (UCE em Aquilinae) | PMC13254548 → Europe PMC | sim, como arquivo (`21098e3a`) | 2026 |
| 28 | *Genes* 2022 (doi 10.3390/genes13071167) — Phylogenomic coalescent analyses of avian retroelements… introgression in Afroaves | PMC9324441 → Europe PMC | sim, como arquivo (`41a1c3e3`) | achado pela pesquisa profunda; a URL do PMC não foi usada |
| 29 | Mather et al. 2021 — *Archaehierax sylvestris*, acipitrídeo do Oligoceno superior da Austrália | https://newsroom.taylorandfrancisgroup.com/wp-content/uploads/2021/09/basal-raptor.pdf | sim (`8516fa6d`) | |
| 30 | AOS/NACC, proposta 2025-C — limites genéricos de *Accipiter* | https://americanornithology.org/wp-content/uploads/2025/05/2025-C.pdf | sim (`b5e56e8d`) | documento oficial |
| 31 | The Timetree of Life (2009), cap. 63 — Falconiformes (arquivo "Brown2009Chap63") | https://TimeTree.org/public/data/pdf/Brown2009Chap63.pdf | sim (`3af8052d`) | |

**Importadas da pesquisa profunda (`research_import`, só acadêmicas):**

| # | Obra | URL | Importada | Obs. |
|---|---|---|---|---|
| 32 | Psittacomimidae fam. nov., *Vertebrate Zoology* (zigodáctilos do Eoceno) | https://vertebrate-zoology.arphahub.com/article/182505/ | sim (`5102aa83`) | |
| 33 | A record of a landbird (Telluraves) from the Eocene Ikovo locality (Ucrânia) | https://www.researchgate.net/publication/308402706_A_record_of_a_landbird_Telluraves_from_the_Eocene_Ikovo_locality_East_Ukraine | sim (`6ccf22e4`) | ResearchGate |
| 34 | Sangster et al. 2022, *Avian Research* — Phylogenetic definitions for 25 higher-level clade names (HTML) | https://www.sciopen.com/article/10.1016/j.avrs.2022.100027 | sim (`1428d7cb`) | traz tabela de suporte por estudo |
| 35 | idem (PDF) | http://www.j-avianres.com/cn/article/pdf/preview/10.1016/j.avrs.2022.100027.pdf | sim (`2be08df4`) | |
| 36 | Braun, Cracraft et al. 2019 — Resolving the Avian Tree of Life from Top to Bottom (cap. em *Avian Genomics in Ecology and Evolution*, pp. 151–210; doi 10.1007/978-3-030-16477-5_6) | https://www.researchgate.net/publication/334115560_Resolving_the_Avian_Tree_of_Life_from_Top_to_Bottom_The_Promise_and_Potential_Boundaries_of_the_Phylogenomic_Era | sim (`47d8810c`) | **a página inclui trechos de artigos que citam o capítulo (ver auditoria Q9)** |
| 37 | Psittacopasseres: Psittaciformes and Passeriformes (capítulo) | https://www.researchgate.net/publication/358558725_Psittacopasseres_Psittaciformes_Parrots_and_Passeriformes_Passerines | sim (`0eef4b43`) | ResearchGate |
| 38 | Tip-dating and the Origin of Telluraves | https://www.researchgate.net/publication/328590668_Tip-dating_and_the_Origin_of_Telluraves | sim, **só resumo** (`79ba75e1`) | ResearchGate |
| 39 | Chen (A.) 2026 — Birds, Diversification of, *Encyclopedia of Evolutionary Biology*, 2ª ed., pp. 39–49 (doi 10.1016/B978-0-443-15750-9.00054-9) | https://www.researchgate.net/publication/384931118_Birds_Diversification_of | sim (`d3d5a7ff`) | ResearchGate; a página inclui também os **resumos das referências** (por exemplo, Mirarab et al. 2024 e Brocklehurst & Field 2024) |
| 40 | Mather 2021 — Taxonomy of fossil eagles and vultures (Accipitridae) from Australia (tese) | https://flex.flinders.edu.au/file/3a909d97-a072-4635-a1df-69df771bb4b0/1/Mather%202021.pdf | sim (`c24ce4bf`) | |
| 41 | Dalsätt — Fossil birds: contributions to the understanding of avian evolution (tese) | https://su.diva-portal.org/smash/get/diva2:513689/FULLTEXT01.pdf | sim (`6899a064`) | inclui corujas fósseis |

### Fontes pedidas ou tentadas que **falharam**

| Obra | URL tentada | Motivo | ID no caderno (lixo) |
|---|---|---|---|
| Suh 2011, Suh 2015, Jarvis 2014, Kuhl 2021, Wu 2024 (via PMC) | https://pmc.ncbi.nlm.nih.gov/articles/PMC3265382/, …/PMC4540587/, …/PMC4405904/, …/PMC7783168/, …/PMC10895254/ | o PMC devolve "Checking your browser - reCAPTCHA" (refeitas via Europe PMC/efetch) | `ad5ef5a1`, `35f4a2c8`, `18ad5efb`, `ae91a8cf`, `0f1c6b8d` |
| Reddy et al. 2017 (texto completo) | https://doi.org/10.1093/sysbio/syx041 | desafio Cloudflare ("Just a moment...") no OUP. Ficou só o resumo, pelo repositório SI | `9a89527e` |
| **Mindell, Fuchs & Johnson 2018**, capítulo em *Birds of Prey* (Springer) | https://doi.org/10.1007/978-3-319-73745-4_1 | HTTP 406 por URL. Pelo curl o capítulo é pago (devolve HTML, não PDF). Sem cópia aberta encontrada | `31bc9ab4` |
| Mayr 2021 (PalZ) e Mayr & Hurum 2020 (Sci Nat) por URL | https://doi.org/10.1007/s12542-020-00541-8, https://doi.org/10.1007/s00114-020-01703-z | HTTP 406 (refeitas por PDF, ver itens 20 e 21) | `e055bd75`, `d52c907f` |
| **Wink et al. 2009**, *Ardea* 97:581 — corujas, cyt b + RAG-1 | BioOne `.full` e `.pdf` (10.5253/078.097.0425) | bloqueio Incapsula / conteúdo vazio | `06eaf3b1`, `1a353216` |
| **McClure et al. 2019** (texto completo) | BioOne `.full` (10.3356/0892-1016-53.4.419) | bloqueio Incapsula. Ficou só a página de citação | `90b2bd24` |
| **Degrange et al. 2015**, *JVP* 35:e912656 — Mesembriornithinae / filogenia de Phorusrhacidae | complete.bioone.org …/10.1080/02724634.2014.912656.full | bloqueio Incapsula | `2df01110` |
| **Mayr 2011**, *J Syst Palaeontol* 9:159 — *Messelastur* (Messelasturidae + Halcyornithidae) | https://doi.org/10.1080/14772019.2010.505252 | caiu no portal genérico da T&F, sem o artigo | `a451f2f5` |
| **Ericson et al. 2006**, *Biol Lett* (doi 10.1098/rsbl.2006.0523) — Diversification of Neoaves | https://doi.org/10.1098/rsbl.2006.0523 | desafio Cloudflare. O Europe PMC deu HTTP 500 e o efetch não tem corpo de texto | `6166068e` |
| Catanach et al. 2024, preprint bioRxiv | https://www.biorxiv.org/content/10.1101/2023.07.13.548898.full.pdf | conteúdo vazio (a versão publicada entrou pelo PAR NSF) | `00ab325f` |
| Comparação de métodos filogenômicos (bioRxiv 2024, da pesquisa profunda) | https://www.biorxiv.org/content/10.1101/2024.02.28.582627v1.full.pdf | conteúdo vazio | `775187d3` |
| Tese de Oxford sobre raptorialidade e tarsometatarso (Telluraves) | https://ora.ox.ac.uk/objects/uuid:90912ef0-a777-4dbc-b0d9-7ff42d0c15b7 | não importada: HTTP 403 na verificação | — |
| Mayr 2009 (*Masillaraptor*, 2º esqueleto), Lerner & Mindell 2005, Griffiths et al. 2004 | — | não importadas: não achei texto aberto. Os pontos centrais aparecem citados em Mayr & Kitchener 2022, Catanach 2024 e Mather 2021 | — |

**Itens da pesquisa profunda descartados por não serem acadêmicos:** Wikipedia (Telluraves, Falconiformes, Messel, versões ga/pl), Reddit, blog WordPress, BirdForum, Wikiwand, All Birds Wiki, Center for Biological Diversity (Fossil Owls), Facebook e páginas de figura do ResearchGate. **Descartados por duplicata:** cópias de Wu 2024 (Edwards Lab e ResearchGate), Jarvis 2014 (ResearchGate e Griffith), Stiller 2024 (ResearchGate e PMC), Masillaraptoridae, Braun & Kimball, Chen 2025 (PMC) e Ksepka (PNAS PDF). **Descartado por estar fora do tema:** galiformes do London Clay (T&F 2024).

**Relatório da pesquisa profunda.** O relatório que o NotebookLM gerou ("Phylogenomic Systematics of Avian Apex Predators…") **não foi importado**. Ele cita Reddit, wikis e blogs. Entre os erros que ele contém está *Tynskya waltonensis* atribuída também à Green River Formation. Não usar.

---

## Perguntas e respostas (verbatim, com as citações do NotebookLM)

Os trechos citados foram truncados para caber no arquivo. Os números entre colchetes são os do NotebookLM. A coluna "fonte" vem do mapeamento `source_id` → título, e o código entre crases é o prefixo do ID da fonte no caderno.

### Q1 — Accipitrimorphae, Strigiformes e Coraciimorphae dentro de Afroaves

**Pergunta:** Which relationship do the major phylogenomic studies recover among Accipitrimorphae, Strigiformes and Coraciimorphae within Afroaves, with support per study?

conversation_id: `9dd437f7-18a8-4497-94d0-a9af8e96600a`

**Resposta do NotebookLM (verbatim):**

> Major phylogenomic studies present three principal topological arrangements among **Accipitrimorphae** (Accipitriformes: hawks, eagles, New World vultures, secretarybird, osprey), **Strigiformes** (owls), and **Coraciimorphae** (mousebirds and cavity nesters like woodpeckers, kingfishers, and rollers), alongside evidence for reticulate evolution [1-5]:
>
> ---
>
> ### 1. Accipitrimorphae Sister to (Strigiformes + Coraciimorphae) within Afroaves
> In this classical "Afroaves" topology, Accipitrimorphae forms the basal-most branch of Afroaves, leaving Strigiformes as the sister group to Coraciimorphae [2, 6]:
>
> * **Jarvis et al. (2014) (Total Evidence Nucleotide Tree / TENT):** Recovers Accipitrimorphae as sister to all other Afroaves, with Strigiformes sister to Coraciimorphae [2, 6]. Received **100% exaML bootstrap support** for Afroaves and **100%** for Coraciimorphae [2]. *(Note: In their binned MP-EST species tree analysis of TENT/intron data, an Accipitriformes + Strigiformes clade emerged instead [7]).*
> * **Kuhl et al. (2021) (3′UTR Sequences):** Supports Accipitrimorphae as sister to (Strigiformes + Coraciimorphae) within Afroaves [2, 8]. Received **100% SH-aLRT support** for Afroaves, Coraciimorphae, and Telluraves [2].
> * **Suh et al. (2015) / Retroposon Insertion (RI) Analysis:** Multi-locus retroposon insertion (MPRE) and ASTRAL_BP species trees place Strigiformes as sister to Cavitaves/Coraciimorphae, with Accipitrimorphae branching deeper [2, 9]. Monophyly of Afroaves is supported by **16 conflict-free retroposon markers** [2].
> * **Reddy et al. (2017) (Noncoding Data):** Concatenated noncoding analyses recovered Accipitrimorphae sister to (Strigiformes + Coraciimorphae), though with **<70% ML bootstrap support** for Afroaves and Coraciimorphae [2, 7].
>
> ---
>
> ### 2. Accipitrimorphae + Strigiformes as Sister Taxa (Hieraves)
> More recent large-scale genomic analyses recover a direct sister-group relationship between diurnal raptors (Accipitrimorphae) and nocturnal raptors (Strigiformes), designated as the clade **Hieraves** [4, 5, 10]:
>
> * **Stiller et al. (2024) (B10K Family-Level Genomes):** Recovers Hieraves (Accipitrimorphae + Strigiformes) as sister to the remaining Afroaves (Coraciimorphae) in their main ASTRAL species tree [4, 10].
>   * **Support & Discordance:** Quartet frequencies for Hieraves (35.0%) versus an alternative topology placing Accipitrimorphae sister to remaining Afroaves (34.6%) were statistically indistinguishable (\\(\chi^2\\) test, \\(p = 0.130\\)), contradicting expectations of incomplete lineage sorting (ILS) alone [4]. This non-ILS pattern was attributed to ancestral hybridization or taxon-sampling sensitivity (requiring \\(\ge 138\\) Passeriformes to favor Hieraves) [4, 11]. Their concatenated ML analysis placed Accipitrimorphae alone as sister to remaining Afroaves with **32% bootstrap support** [4, 12].
> * **Wu et al. (2024) (Coding Sequences / 5,756 CDS):** Also recovers the **Hieraves** clade (Accipitrimorphae + Strigiformes) [5, 10]. However, it places Hieraves as sister to **Australaves** (rather than inside Afroaves), leaving Coraciimorphae as the basal-most lineage of Telluraves [5, 13].
>
> ---
>
> ### 3. Accipitrimorphae Sister to All Core Landbirds (Afroaves Non-Monophyletic)
> * **Prum et al. (2015) (Anchored Hybrid Enrichment / 259 Loci):** Places Strigiformes as sister to Coraciimorphae (**1.0 Bayesian posterior probability**; **41% ML bootstrap** for broader nodes), but places Accipitrimorphae as sister to *all* other core landbirds (Telluraves), rendering Afroaves non-monophyletic [2, 14-16].
>
> ---
>
> ### 4. Reticulate Evolution and Ancient Hybridization Scenarios
> * **Genes / Genes2022 Retroelement Analyses:** Re-analyses of retroelement insertions (RIs) reveal asymmetric character support between minority resolutions at the Strigiformes node within Afroaves [17, 18]. Network likelihood tests (**KKSC**, \\(p = 0.046 / 0.014\\); **4-LIN**) significantly reject strictly bifurcating models in favor of an **ancient introgressive hybridization network**, inferring that **Strigiformes** (owls) originated as a hybrid lineage with genetic contributions from the stem of **Coraciimorphae (~55%)** and the stem of **Accipitrimorphae (~45%)** [18, 19].
>
> 🦉 *Would you like to explore the morphological synapomorphies associated with these alternative hypotheses, or examine the divergence timing estimates across these studies?*

**Citações (número → fonte; trecho citado, truncado):**

- [1] Resolving the Avian Tree of Life from Top to Bottom (ResearchGate, `47d8810c`) — "Relationships within core landbirds depend on data type and analytical approach. We present four candidate topologies… (a) Division into two major clades (Australaves and Afroaves…) present in the Jarvis et al. (2014) TENT and the Reddy et al. (2017) analyses of concatenated noncoding data. (b) Accipitriformes sister to all other core landbirds… Prum et al. (2015)… (c) An acciptriform + owl clade sister to all other core landbirds, found in binned MP-EST analysis of the TENT data… (d) … mousebirds sister to all other Afroaves, found in the Jarvis et al. (2014) concatenated UCE tree and the Suh et al. (2015) analysis of TE insertions."
- [2] Phylogenetic definitions for 25 higher-level clade names of birds (SciOpen, `1428d7cb`) — tabela de suporte por estudo. Linhas relevantes: Afroaves — Hackett 2008 <70% ML; Jarvis 2014 100% exaML; Burleigh 2015 50% ML; Prum 2015 not recovered; Suh 2015 16 retroposons; Reddy 2017 <70% ML; Liu 2018 not recovered; Kimball 2019 100% MRP; Kuhl 2021 100% SH-aLRT. Coraciimorphae — Hackett not recovered; Jarvis 100%; Prum 1.0 PP; Suh 2015 not recovered; Reddy <70%; Liu <70%; Kimball 72% MRP; Kuhl 100%.
- [3] idem (`1428d7cb`) — "Etymology: Derived from the word Africa… Reference phylogeny: … Figure 3 in Kuhl et al. (2021) should be regarded as the primary reference phylogeny… Composition: Includes Accipitriformes, Strigiformes, Coliiformes, Leptosomiformes, Trogoniformes, Bucerotiformes, Upupiformes, Coraciiformes, and Piciformes."
- [4] Stiller et al. 2024 (`e9c5ae30`) — "Our tree grouped Accipitriformes and Strigiformes as the sister to the remaining Afroaves… Concatenated analyses, including ours, supported Accipitriformes alone as sister to the remaining Afroaves… quartet frequencies… statistically indistinguishable (35 versus 34.6%, χ2 test, P BH adjusted = 0.130)… the non-ILS patterns could be indicative of ancestral hybridization… The main topology was obtained only when at least 138 Passeriformes were included."
- [5] Wu et al. 2024 (`cda64d61`) — "Within Telluraves, the optimum tree supports… Coraciimorphae as the most basal clade… uncovers a novel monophyletic clade, uniting Accipitriformes (hawks) and Strigiformes (owls), which we designate Hieraves. The Hieraves are sister to Australaves."
- [6] Resolving the Avian Tree (`47d8810c`) — "the falcons are placed in Australaves sister to parrots and passerines whereas the hawks, eagles, New World vultures, and owls are placed in Afroaves as the successive sister groups of a clade comprising mousebirds, and… woodpeckers and their allies (Fig. 2a). However, the Jarvis et al. (2014) TENT conflicts with the Prum et al. (2015) topology with respect to the position of… Accipitriformes."
- [7] Resolving the Avian Tree (`47d8810c`) — "acciptriform + owl clade sister to all other core landbirds, found in binned MP-EST analysis of the TENT data and unbinned MP-EST analysis of introns by Jarvis et al. (2014). Reddy et al. (2017) also found this topology in their concatenated analyses of a 104-locus coding data matrix…"
- [8] Chen et al. 2025, IOB (`6c89ed1d`) — "Within Telluraves, Jarvis et al. (2014), Kuhl et al. (2021), and Stiller et al. (2024) recovered a clade containing Accipitrimorphae, Strigiformes, and Coraciimorphae, forming a group called Afroaves…"
- [9] Genes 2022 retroelementos (`41a1c3e3`) — "The position of Strigiformes (owl) is scattered across published trees… Both RI trees place Strigiformes as sister to Cavitaves (cuckoo-roller, trogon, hornbill, bee-eater, woodpecker), with Accipitriformes… and Coliiformes… branching from successively deeper nodes within a monophyletic Afroaves."
- [10] Chen et al. 2025 (`6c89ed1d`) — "Stiller et al. (2024) and Wu et al. (2024b) instead recovered Accipitrimorphae and Strigiformes as sister groups, forming a clade called Hieraves…"
- [11] Stiller et al. 2024 (`e9c5ae30`) — "…14 nodes had quartet support below 37%… the alternative placement of Afroaves + Accipitriformes had higher quartet support when only a few species were sampled but declined with increasing taxon sampling… The main topology dominated when 138 or more passerines were sampled."
- [12] Stiller et al. 2024 (`e9c5ae30`) — "Our main tree placed Strisores… with Phaethoquornithes with moderate support (posterior probability = 0.90…), but the concatenated tree grouped them as sister to Telluraves with low support (bootstrap = 32…)."
- [13] Chen et al. 2025 (`6c89ed1d`) — "Whereas Stiller et al. (2024) placed Hieraves within Afroaves, Wu et al. (2024b) recovered it as the sister group to Australaves…"
- [14] Resolving the Avian Tree (`47d8810c`) — "TENT (and Reddy et al. 2017) places accipitriforms sister to all other Afroaves whereas Prum et al. (2015) places them sister to all other core landbirds… The Jarvis et al. (2014) analysis of first and second codon positions supported a third topology, with a clade comprising accipitriforms and owls sister to all other core landbirds."
- [15] Phylogenetic definitions (`1428d7cb`) — "Most studies support a dichotomy between Accipitriformes and all other Afroaves (Hackett 2008; Jarvis 2014; Burleigh 2015; Reddy 2017; Kimball 2019; Kuhl 2021) but one study inferred a dichotomy between Coliiformes and all other Afroaves (Suh et al. 2015…). … In two recent phylogenomic studies, Accipitriformes was sister to a large clade that included… also Australaves (Prum et al., 2015; Liu et al., 2018)."
- [16] Chen et al. 2025 (`6c89ed1d`) — "Prum et al. (2015) instead found a Strigiformes + Coraciimorphae clade as sister to Australaves, forming a group called Eutelluraves…"
- [17] Genes 2022 (`41a1c3e3`) — "…skewed RI support for minority resolutions at the controversial Strigiformes + Cavitaves node within Afroaves… an excess of support for Strigiformes + Accipitriformes relative to Cavitaves + Accipitriformes…"
- [18] Genes 2022 (`41a1c3e3`) — "The KKSC and 4-LIN approaches… inferred networks with T. alba (owl…) interpreted as a hybrid taxon with genetic contributions from L. discolor (cuckoo-roller; Cavitaves) and H. leucocephalus (eagle 2; Accipitriformes). For KKSC, a network was supported (p = 0.046)… KKSC was again significant for owl, eagle 2, and cuckoo-roller (p = 0.014)."
- [19] Genes 2022 (`41a1c3e3`) — "…the 4-LIN likelihood-based approach estimated approximately equal genetic contributions from cuckoo-roller (55%) and eagle 2 (45%)…"

---

### Q2 — Falconiformes como irmão de Psittacopasserae; posição de Cariamiformes

**Pergunta:** What evidence (data types, retroposons, support) places Falconiformes as sister to Psittacopasserae, and where do Cariamiformes fall?

conversation_id: `31eb97e3-f1de-4388-89fe-658780152b98`

**Resposta do NotebookLM (verbatim):**

> Molecular and genomic evidence robustly supports placing **Falconiformes** (falcons) as the sister group to **Psittacopasseres** (parrots + passerines), forming the clade **Eufalconimorphae** [1-3]. In turn, **Cariamiformes** (seriemas and extinct relatives like "terror birds") branch off as the **basal sister group to Eufalconimorphae**, uniting all four orders into the major landbird clade **Australaves** [1, 4-6].
>
> ---
>
> ### 1. Evidence Placing Falconiformes as Sister to Psittacopasseres (**Eufalconimorphae**)
>
> #### **Retroposon Insertions (Rare Genomic Changes)**
> * **Conflict-Free Markers:** Retroposon presence/absence insertions (such as CR1 LINEs and LTR retrotransposons) are virtually homoplasy-free markers [7, 8]. The landmark study by Suh et al. (2011) identified **7 conflict-free retroposon insertions** shared exclusively by falcons, parrots, and passerines to the exclusion of all other birds, providing statistically significant support (\\(P = 0.0005\\)) for the clade **Eufalconimorphae** [3, 9].
> * **Expanded Genomic Support:** Subsequent genome-wide retrotransposon screenings identified **up to 88 retroposon markers** supporting Eufalconimorphae [10-12]. Additionally, **3 retroposons** (and up to 16 in later analyses) specifically unite parrots and passerines as **Psittacopasseres** (\\(P = 0.0370\\)) [3, 9, 10].
>
> #### **Data Types & Phylogenomic Support**
> The sister relationship between Falconiformes and Psittacopasseres has been consistently recovered across a broad variety of molecular data types:
> * **Nuclear Introns & 19-Loci Datasets:** First recovered with moderate-to-high support in multigene nuclear studies (e.g., Hackett et al. 2008) [3, 13, 14].
> * **Ultraconserved Elements (UCEs):** Sequence capture of over 1,500 UCE loci strongly resolves Eufalconimorphae and Australaves [10, 11, 15].
> * **Anchored Hybrid Enrichment (AHE):** Hundreds of genomic loci (e.g., Prum et al. 2015) yield maximal posterior probability (1.0 PP) for Eufalconimorphae and Australaves [10, 11, 16].
> * **3′-UTRs & Whole Genomes:** Non-coding 3′-untranslated regions (Kuhl et al. 2021) and family-level whole-genome coalescent analyses (Jarvis et al. 2014, Stiller et al. 2024, Wu et al. 2024) consistently obtain maximal statistical support (100% bootstrap / SH-aLRT) for Eufalconimorphae [6, 10, 17-19].
> * *Note on Data Type Effects:* While non-coding sequences (introns, UCEs, 3′-UTRs) overwhelmingly support Eufalconimorphae, unpartitioned or uncorrected **protein-coding exons** and **mitochondrial genomes** sometimes suffer from compositional biases or long-branch attraction, occasionally misplacing falcons near other raptorial groups or mousebirds [20-22].
>
> #### **Phenotypic and Behavioral Traits**
> * **Morphology:** Classical osteological matrices (e.g., Livezey & Zusi 2007) fail to identify clear skeletal synapomorphies for Eufalconimorphae due to extreme morphological divergence between raptorial falcons and songbirds/parrots [11, 12, 23].
> * **Wing Molt Sequence:** Pyle (2013) identified a distinctive **wing molt sequence** shared between falcons (Falconiformes) and parrots (Psittaciformes), which serves as a potential behavioral/pterylo-morphological synapomorphy of Eufalconimorphae (lost or modified in passerines) [11, 12].
>
> ---
>
> ### 2. Where Cariamiformes Fall
>
> #### **Phylogenetic Position**
> * **Basal Divergence of Australaves:** Cariamiformes (comprising extant seriemas and extinct families like Phorusrhacidae, Bathornithidae, and Idiornithidae) form the **deepest branch (sister taxon)** to Eufalconimorphae within **Australaves** [4, 5, 11].
> * **Genomic Support:** The monophyly of **Australaves** (Cariamiformes + [Falconiformes + Psittacopasseres]) is supported by **2 initial retroposons** (Suh et al. 2011) and **up to 176 retroposons** in broader screenings (Suh et al. 2015), as well as near-100% support across all major phylogenomic sequence analyses [10, 24, 25].
>
> #### **Biogeographical and Evolutionary Context**
> * **Gondwanan Radiation:** Historical biogeographical reconstructions indicate that Australaves represents a major Gondwanan landbird radiation, with the **Cariamiformes and Falconiformes lineages diversifying in South America**, while the **Psittacopasseres lineage diversified in Australasia** [4, 26].
> * **Fossil Links:** Early Paleogene fossil discoveries (such as Eocene *Masillaraptoridae*) reveal a mosaic of falconid features (beak and raptorial claws) and cariamiform-like plesiomorphic characters (such as hypotarsus structure), providing morphological evidence that stems near the split between Cariamiformes and Falconiformes [20].
>
> ---
>
> Would you like to explore the fossil record of early stem-falcons and stem-parrots (such as Masillaraptoridae or Messelasturidae), or dig deeper into how retroposon incomplete lineage sorting affected early Neoavian branches?

**Citações (número → fonte; trecho citado, truncado):**

- [1] Psittacomimidae fam. nov., *Vertebrate Zoology* (`5102aa83`) — "A clade including both taxa was first found in analyses of Ericson et al. (2006) and has since been recovered in all phylogenies based on nuclear gene sequences (Hackett 2008; Jarvis 2014; Prum 2015; Kuhl 2021; Stiller 2024; Wu 2024)… termed Psittacopasserae (Suh et al. 2011), emended to… Psittacopasseres by Sangster et al. (2022). Sequence-based studies, furthermore, identify the Falconiformes as the sister taxon of the Psittacopasseres. Together with the Cariamiformes… form the clade Australaves."
- [2] Phylogenetic definitions (SciOpen, `1428d7cb`) — composição de Eufalconimorphae (Falconiformes 64 spp., Psittaciformes 376, Passeriformes 6063) e táxons extintos (Quercypsittidae, Halcyornithidae, Messelasturidae, Vastanavidae, Psittacopedidae, Zygodactylidae).
- [3] Suh et al. 2011 (`7cfbfef1`) — "Seven retroposon insertions are exclusively present in falcons, parrots and passerines, but absent in hawks, woodpeckers and other 'landbirds' (… seven REs, P=0.0005…); we therefore suggest the new name Eufalconimorphae… three retroposon insertions solely in parrots and passerines (… P=0.0370…)… Psittacopasserae."
- [4] Birds, Diversification of (ResearchGate, `d3d5a7ff`) — "Psittacopasseres, which is closely related to Cariamiformes and Falconiformes, forming the group Australaves…"
- [5] Phylogenetic definitions (`1428d7cb`) — composição de Australaves, incluindo Idiornithidae, Bathornithidae e Phorusrhacidae como grupo-tronco de Cariamiformes (Angst 2013; Degrange 2015; Mayr 2017, 2022).
- [6] Wu et al. 2024 (`cda64d61`) — "…The Hieraves are sister to Australaves (34), a known clade consisting of Passeriformes, Psittaciformes, Falconiformes, and Cariamiformes."
- [7] Genes 2022 retroelementos (`41a1c3e3`) — resumo: "Retroelement insertions (RIs) are low-homoplasy characters…"
- [8] Suh et al. 2011 (`7cfbfef1`) — resumo: "…virtually homoplasy-free retroposon insertions… ~200,000 retroposon-containing loci… 51 markers…"
- [9] Suh et al. 2011 (`7cfbfef1`) — legenda da Fig. 1 (árvore de presença/ausência de retroposons; Eufalconimorphae e Psittacopasserae).
- [10] Phylogenetic definitions (`1428d7cb`) — tabela. Australaves: Hackett 64% ML; McCormack 2013 (em branco); Jarvis 100%; Burleigh 72%; Prum 1.0 PP; Suh 2015 176 retroposons; Reddy ≥95%; Liu <70%; Kimball 100%; Kuhl 100%. Eufalconimorphae: Hackett 73%; McCormack >0.95 PP/>70%; Jarvis 100%; Burleigh 59%; Prum 1.0; Suh 2015 88 retroposons; Reddy <70%; Liu 85%; Kimball 100%; Kuhl 100%. Psittacopasserae: Hackett 77%; Jarvis 100%; Burleigh 47%; Prum 1.0; Suh 2015 16; Reddy ≥95%; Liu 92%; Kimball 100%; Kuhl 100%.
- [11] Phylogenetic definitions (`1428d7cb`) — "Suh et al. (2011) identified three retroposons that are apomorphic for this clade [Psittacopasseres]… Livezey and Zusi… revealed no morphological apomorphies. However, Pyle (2013) observed that falcons and parrots share a distinctive wing molt sequence… synapomorphy of Eufalconimorphae which is lost in Passeriformes… The name Australaves was used by Ericson (2012) for a clade that also included Cariamidae, the latter being depicted as the sister taxon of Falconidae… all subsequent phylogenomic studies have placed Cariamidae as the sister to the falcon-parrot-passerine clade."
- [12] idem, versão PDF da *Avian Research* (`2be08df4`).
- [13] Resolving the Avian Tree (`47d8810c`) — "the falcons are placed in Australaves sister to parrots and passerines…"
- [14] Hackett et al. 2008 (`cad68562`) — "One of the most unexpected findings was the sister relationship between Passeriformes and Psittaciformes (node A), with Falconidae (falcons) sister to this clade. This relationship varied slightly among analyses and gene-jackknifing…"
- [15] Resolving the Avian Tree (`47d8810c`) — legenda das quatro topologias de core landbirds (mesmo trecho de Q1 [1]).
- [16] Phylogenetic definitions (`1428d7cb`) — "…multiple genomic support for the basal dichotomy within Australaves (i.e. between Cariamidae and Eufalconimorphae) (Hackett 2008; Jarvis 2014; Burleigh 2015; Prum 2015; Suh 2015; Liu 2018; Kimball 2019; Kuhl 2021)."
- [17] Jarvis et al. 2014 (`c92cd31c`) — "Within core landbirds, we found 100% BS for a previously more weakly supported clade (Australaves) containing seriemas…, falcons…, parrots…, and Passeriformes and a sister clade (Afroaves)…"
- [18] Phylogenetic definitions (*Avian Research* PDF, `2be08df4`) — "…congruent support for the basal dichotomy within Eufalconimorphae (Hackett 2008; Suh 2011, 2015; McCormack 2013; Jarvis 2014; Prum 2015; Kimball 2019; Kuhl 2021)."
- [19] Stiller et al. 2024 (`e9c5ae30`) — "Within Neoaves we resolve four major clades… Telluraves (higher landbirds including Afroaves and Australaves)… Elementaves…"
- [20] Mayr & Kitchener 2022 (`53dae551`) — "…In most current sequence-based analyses, the Cariamiformes and Falconiformes are successive sister taxa of the clade including the Psittaciformes and Passeriformes… By contrast, unpartitioned coding exons support a sister group relationship between the Cariamiformes and the Falconiformes (Braun and Kimball 2021: fig. 4D). The mosaic of derived falconiform features and presumably plesiomorphic cariamiform characters (e.g., a block-like hypotarsus) present in masillaraptorids… lends some support to a sister group relationship between the Cariamiformes and Falconiformes… Alternatively… parallelisms… or plesiomorphic for the clade including the Cariamiformes, Falconiformes, Psittaciformes, and Passeriformes, the stem species of which would then most likely have been a long-legged bird with a raptorial way of living."
- [21] Kuhl et al. 2021 (`8527a287`) — "In the CODON, CODON12, and the AAS trees, the Falconiformes and Cariamiformes were moved away from the Passeriformes and Psittaciformes to form an assemblage of birds of prey… In the AAS tree, even the Strigiformes were enclosed in the birds of prey assemblage… In the CODON12 and AAS tree, the Coliiformes were the closest relatives of the Passeriformes and Psittaciformes."
- [22] Wang et al. 2014, GBE (`3cc767df`) — "…a relatively well-supported clade (PP = 1.0, ML = 63%) between Strigiformes (owls) and Psittaciformes (parrots)…"
- [23] Mayr & Kitchener 2022 (`53dae551`) — "Falcons were long considered to be closely related to other diurnal birds of prey. However, most current molecular analyses recover them in a clade together with seriemas…, parrots…, and passerines…"
- [24] Phylogenetic definitions (`1428d7cb`) — "[Australaves] No morphological apomorphies are known. Suh et al. (2011) identified two retroposons that are apomorphic for this clade… first documented by Ericson et al. (2006)…"
- [25] idem, versão PDF (`2be08df4`).
- [26] Ericson 2012 (`254c47a3`) — "…three reciprocally monophyletic groups of terrestrial birds have diversified in the Gondwanan land areas of Australia, South America and Africa… the only survivors today from the Australian radiation are the passerines and parrots, while the falcons and seriemas have survived from the South American radiation…"

---

### Q3 — Monofilia de Accipitrimorphae e ordem Cathartidae / Sagittariidae / Pandionidae / Accipitridae

**Pergunta:** Is Accipitrimorphae (Cathartidae + Accipitriformes) monophyletic, and what is the branching order of Cathartidae, Sagittariidae, Pandionidae and Accipitridae?

conversation_id: `21426d54-430d-4629-b961-d2ef1523eee5`

**Resposta do NotebookLM (verbatim):**

> **Yes, Accipitrimorphae (Cathartidae + Accipitriformes) is monophyletic.** Modern phylogenomic analyses, ultraconserved element (UCE) studies, and morphological assessments consistently recover **Cathartidae** (New World vultures) as the monophyletic sister group to all other **Accipitriformes** [1-5]. Morphological synapomorphies—such as caudally inflected acromial processes on the furcula—as well as extensive molecular data support this clade [1, 5].
>
> ---
>
> ### **Branching Order**
>
> The four families form a clear series of successive sister-group relationships, listed below from the **most basal/outermost split to the most deeply nested clade**:
>
> 1. **Cathartidae** (Cathartiformes / New World vultures): Diverged first as the sister group to all other Accipitriformes (~61.3 Mya) [4-6].
> 2. **Sagittariidae** (Secretarybird): Diverged next as the sister group to the combined **Pandionidae + Accipitridae** clade (~60.9 Mya) [4, 6-8].
> 3. **Pandionidae** (Osprey): Diverged next (~50.8 Mya), forming the immediate sister group to **Accipitridae** [4, 6, 7, 9].
> 4. **Accipitridae** (hawks, eagles, Old World vultures, harriers, kites): The terminal, hyper-diverse crown clade [4, 6, 7].
>
> In cladogram notation, the topology is:
> \\[\mathbf{(Cathartidae,\ (Sagittariidae,\ (Pandionidae,\ Accipitridae)))}\\] [4, 7, 8]
>
> ---
>
> 🦅 *Would you like to explore how the major subfamilies within **Accipitridae** (such as Elaninae, Perninae, Aegypiinae, and Buteoninae) branch off, or examine why **falcons (Falconidae)** were split into a completely separate order?*

**Citações (número → fonte; trecho citado, truncado):**

- [1] Chen et al. 2025, IOB (`6c89ed1d`) — "Accipitrimorphae: Seven character states were optimized as potential synapomorphies of Accipitrimorphae, of which four were consistently inferred as such across all alternative molecular topologies studied. Acromial processes on furcula strongly caudally inflected (char. 6: 0 > 1…). Median trabecula of sternum rounded at caudal termination (char. 61: 1 > 0…)."
- [2] Mather 2021, tese (Flinders, `c24ce4bf`) — "…the Cathartidae within the Accipitriformes, placing them as the most basal lineage within the group (Hackett 2008; Jarvis 2014; Prum 2015). As such, the Cathartidae can be classified as either a primitive lineage of the Accipitriformes or form a very closely related sister group as the Cathartiformes…"
- [3] Hackett et al. 2008 (`cad68562`) — "Cathartidae (New World vultures)—traditionally placed within Falconiformes—were previously allied with Ciconiidae (storks). However, our phylogeny demonstrated no affinity with storks… and instead strongly supported placement of Cathartidae within the land birds (usually with Accipitridae)."
- [4] Catanach et al. 2024 (PAR NSF, `1d4669b1`) — "The order Cathartiformes formed a clade that was sister to Accipitriformes. Within Accipitriformes, the family Sagitariidae was sister to the rest of Accipitriformes, and the family Pandionidae was sister to the rest of Accipitridae. Several large clades within Accipitridae corresponded roughly… to the subfamilies recognized by Mindell et al. (2018) and Lerner and Mindell (2005)."
- [5] Catanach et al. 2024 (`1d4669b1`) — "Our analysis, which included UCE data from all five extant genera and legacy data from all currently recognized extant species, recovered Cathartiformes as the sister group of Accipitriformes. This finding is supported by multiple morphological synapomorphies (e.g. Griffiths 1994)."
- [6] Catanach et al. 2024 (`1d4669b1`) — "The orders Cathartiformes and Accipitriformes were estimated to have diverged at 61.3 Mya. The common ancestor of extant taxa within Cathartidae was estimated at 16.2 Mya. The family Sagittariidae diverged from the rest of Accipitriformes at 60.9 Mya, and the family Pandionidae diverged from Accipitridae at 50.8 Mya. The clades corresponding to subfamilies within Accipitridae diverged between 25.1 and 33 Mya, except for Elaninae, which split from the rest of Accipitridae at 45.3 Mya…"
- [7] Mather 2021 (`c24ce4bf`) — "Pandionidae is directly sister to the Accipitridae clade, and the Sagittariidae fall sister to the clade comprising Pandionidae and Accipitridae…"
- [8] Catanach et al. 2024 (`1d4669b1`) — "The phylogenetic placement of these monotypic families was identical to previous studies. We estimated that Sagittariidae diverged from the Pandionidae + Accipitridae clade ~60.9 Mya… and that Pandionidae diverged from Accipitridae ~50.8 Mya. These estimates were similar to those inferred by Johnson et al. (2016) and Knapp et al. (2019), but considerably older than those of Mindell et al. (2018) (~45 and ~38 Mya, respectively) and Prum et al. (2015) (~40 and ~28 Mya, respectively). These discrepancies probably resulted from different calibration methods…"
- [9] Mather 2021 (`c24ce4bf`) — "Pandionidae, Sagittariidae and Cathartidae are successive sister-taxa to the Accipitridae within the Accipitriformes…"

---

### Q4 — Relações internas de Accipitridae e consequências da não-monofilia de *Accipiter*

**Pergunta:** What are the internal relationships of Accipitridae (subfamilies) and the consequences of Accipiter non-monophyly?

conversation_id: `0666af6b-455e-47fb-a3e8-d00a198091b6`

**Resposta do NotebookLM (verbatim):**

> ### Internal Relationships of Accipitridae (Subfamilies)
>
> Phylogenomic and molecular studies (such as UCE backbone trees and multi-gene matrices) have clarified the internal branching order of the family **Accipitridae**, resolving several major subfamily clades from most basal to most derived [1-4]:
>
> * **Elaninae** (kites such as *Elanus*, *Gampsonyx*, and *Chelictinia*): Forms the earliest-diverging, most basal clade within Accipitridae (splitting ~34–45 Ma) [1, 5-7].
> * **Perninae & Gypaetinae**: A combined clade containing pernine kites (*Pernis*, *Elanoides*, *Chondrohierax*, *Aviceda*) and Old World vultures/harrier-hawks (*Gypaetus*, *Neophron*, *Gypohierax*, *Polyboroides*) forms the sister group to all remaining accipitrids [6, 8, 9].
> * **Circaetinae & Aegypiinae**: These two subfamilies are reciprocally monophyletic sister groups [10]. Circaetinae includes snake eagles (*Circaetus*, *Spilornis*, *Terathopius*) and the Philippine eagle (*Pithecophaga*), while Aegypiinae comprises the obligate Old World scavenging vultures (*Gyps*, *Aegypius*, *Torgos*, *Sarcogyps*) [8-10].
> * **Aquilinae** (Booted eagles): A monophyletic subfamily characterised by feathered tarsi, encompassing genera such as *Aquila*, *Hieraaetus*, *Spizaetus*, *Nisaetus*, *Clanga*, and *Stephanoaetus* [10-12].
> * **Harpiinae**: Large forest eagles (*Harpia*, *Harpyopsis*, *Morphnus*, *Macheiramphus*) form a sister clade to the remaining core crown accipitrids [4, 8].
> * **Lophospizinae & Harpaginae**: Deep lineages sitting basal to the core hawks/buzzards; *Lophospiza* (*L. trivirgatus*, *L. griseiceps*) diverges at the base of the Buteoninae + Accipitrinae radiation [13-15], while Harpaginae includes *Harpagus* and *Microspizias* [4, 16].
> * **Accipitrinae & Buteoninae**: The most species-rich, highly derived radiation [4, 17]. Buteoninae includes buzzards, hawks, and sea-eagles (*Buteo*, *Buteogallus*, *Geranoaetus*, *Haliaeetus*, *Milvus*), while Accipitrinae contains the core true hawks and harriers [4, 8, 18].
>
> ---
>
> ### Non-Monophyly of *Accipiter* and Its Taxonomic Consequences
>
> Historically, ***Accipiter*** served as a megadiverse "catch-all" genus for ~46–50 species of forest hawks and sparrowhawks [19-21]. Dense phylogenomic sampling (including ultraconserved elements) revealed that traditional *Accipiter* is heavily **polyphyletic and paraphyletic** [20, 22-24]:
>
> 1. **Embedded & Interspersed Genera**:
>    * Morphologically and ecologically distinct genera—such as the harriers (***Circus***), Doria's goshawk (***Megatriorchis***), and the red goshawk (***Erythrotriorchis***)—are nested deeply *within* the main *Accipiter* complex [3, 24-26].
>    * Other distinct African genera (***Kaupifalco***, ***Melierax***, ***Urotriorchis***, ***Micronisus***) branch off between different *Accipiter* lineages [4, 25, 27].
>    * Certain species traditionally assigned to *Accipiter* actually belong to distant lineages (e.g., *Accipiter trivirgatus* sits near the base of all Accipitrinae + Buteoninae) [4, 13, 28].
>
> 2. **Taxonomic Revisions & Breakdown into Monophyletic Genera**:
>    To establish monophyly without lumping distinct genera like harriers (*Circus*) into *Accipiter*, taxonomic authorities (including the IOC World Bird List, Clements Checklist, WGAC, and AOS/NACC following Catanach et al. 2024 and Sangster et al. 2021) have split former *Accipiter* into **5–6 distinct genera** [20, 22, 25, 26, 29]:
>
>    * ***Lophospiza***: Re-erected for the Crested Goshawk (*L. trivirgatus*) and Sulawesi Goshawk (*L. griseiceps*), placed in the separate subfamily **Lophospizinae** [13-15].
>    * ***Microspizias***: Erected for the Tiny Hawk (*M. superciliosus*) and Semicollared Hawk (*M. collaris*), placed in **Harpaginae** sister to *Harpagus* [4, 16, 30].
>    * ***Tachyspiza***: Accommodates ~27 species of mostly Indomalayan, Australasian, and African hawks with relatively short toes (e.g., *T. badius*, *T. soloensis*, *T. novaehollandiae*, *T. fasciatus*, *T. hiogaster*) [31, 32].
>    * ***Aerospiza***: Accommodates 3 medium-sized African forest species (*A. tachiro*, *A. castanilius*, *A. toussenelii*) [33].
>    * ***Accipiter* (sensu stricto)**: Restricted to the type species (*A. nisus*, Eurasian Sparrowhawk) and 5 small bird-hunting species (*A. striatus*, *A. rufiventris*, *A. madagascariensis*, *A. ovampensis*, and tentatively *A. poliogaster*) [34-36].
>    * ***Astur***: Resurrected for ~9 larger goshawk species with stout bills and strong halluces, including the Northern/American Goshawk (*Astur gentilis* / *A. atricapillus*), Cooper's Hawk (*Astur cooperii*), *A. bicolor*, and *A. melanoleucus* [37, 38].
>
> 3. **Re-sequencing of Harriers (*Circus*)**:
>    Because *Circus* is sister to the *Astur* + *Megatriorchis* clade, global linear checklists have updated their sequence to place harriers immediately following *Astur* [25, 26, 38-40].
>
> 🦅 Would you like to explore how fossil raptors (such as *Archaehierax* or *Masillaraptoridae*) fit into this new phylogenetic framework?

**Citações (número → fonte; trecho citado, truncado):**

- [1] Mayr & Hurum 2020, *Sci. Nat.* (`32bd8904`) — "Accipitrid diurnal birds of prey include 237 extant species… Well-resolved molecular phylogenies (Lerner and Mindell 2005; Griffiths et al. 2007; Mindell et al. 2018) show that various small clades with geographically restricted distributions diverged early in the phylogeny of the Accipitridae, whereas the species-rich taxa Aquilinae…, Buteoninae…"
- [2] Mather 2021, tese (`c24ce4bf`) — "…more derived subfamilies (Aquilinae, Buteoninae and Haliaeetinae) tended to date between 17–10 Ma. Interestingly, the Aegypiinae dated between 13–7 Ma… The most recent published phylogenetic work was done by Mindell et al. (2018)…"
- [3] Catanach et al. 2024 (`1d4669b1`) — "The UCE backbone phylogeny of 134 species… was fully supported at virtually every node (bootstrap values = 100). Only five nodes had less than perfect support… Of the 20 genera for which multiple species were sampled…, all but Accipiter and Circaetus were monophyletic. The non-monophyly of Accipiter was caused primarily by the embedded placement of Circus, Erythrotriorchis, and Megatriorchis. Additionally, Dryotriorchus was embedded in the genus Circaetus…"
- [4] Catanach et al. 2024 (`1d4669b1`) — "…a clade composed of genera traditionally placed in the subfamily Harpiinae (Macheiramphus, Harpyopsis, Morphnus, and Harpia) was the sister group of the rest of Accipitridae… the monotypic [Accipiter] trivirgatus was the sister group to the rest of Accipitridae, which was divided into three major clades. The first of these was the 'Accipiter' complex (Kaupifalco, Micronisus, Melierax, Urotriorchis, Erythrotriorchis, Megatriorchis, Circus, and Accipiter)… The final two major clades were reciprocally monophyletic. One contained the genus Harpagus and two 'Accipiter' species… Microspizias…, and the other clade corresponded to the subfamily Buteoninae… Milvini (Milvus, Haliastur, and Haliaeetus) and Buteonini (…)."
- [5] Mather 2021 (`c24ce4bf`) — "…the most basal living genera, Gampsonyx and Elanus, diverged from the rest of the Accipitridae roughly 34 Ma, and the split between the Perninae+Gypaetinae clade and the rest of the accipitrid subfamilies roughly 23 Ma, in agreement with Nagy and Tökölyi (2014)…"
- [6] Catanach et al. 2024 (`1d4669b1`) — "The most basal split within Accipitridae separated the monophyletic subfamily Elaninae (Elanus, Gampsonyx, and Chelictinia)… a clade containing all genera traditionally placed in Gypaetinae… and Perninae… formed the sister group of the remainder of Accipitridae. However, the African harrier-hawk, Polyboroides typus…, was placed as sister to the rest of the Gypaetinae + Perninae clade, which rendered Gypaetinae nonmonophyletic."
- [7] Catanach et al. 2024 (`1d4669b1`) — datas (mesmo trecho de Q3 [6]): "…subfamilies within Accipitridae diverged between 25.1 and 33 Mya, except for Elaninae, which split from the rest of Accipitridae at 45.3 Mya…"
- [8] Mather 2021 (`c24ce4bf`) — lista de subfamílias e legenda da Fig. 1.1 (Lerner & Mindell 2005; Nagy & Tökölyi 2014).
- [9] Mather 2021 (`c24ce4bf`) — "Old World vultures fall strictly within Accipitridae. They are not monophyletic, being divided into two subfamilies, the Aegypiinae and Gypaetinae, the latter of which forms a clade with the pernine kites and the former with the serpent eagles…"
- [10] Catanach et al. 2024 (`1d4669b1`) — "…a clade containing all members of the subfamilies Circaetinae (Spilornis, Pithecophaga, Terathopius, Circaetus, and Dryotriorchis) and Aegypiinae (Sarcogyps, Trigonoceps, Torgos, Aegypius, Necrosyrtes, and Gyps), which were reciprocally monophyletic, formed the sister group to the rest of Accipitridae… a clade containing the subfamily Aquilinae (Stephanoaetus, Nisaetus, Spizaetus, Lophotriorchis, Polemaetus, Lophaetus, Ictinaetus, Clanga, Hieraaetus, and Aquila) formed the sister group to the rest of Accipitridae."
- [11] Mather 2021 (`c24ce4bf`) — "Their common name of 'booted eagles' refers to the fact that all species in this subfamily have feathered tarsi…"
- [12] Catanach et al. 2024 (`1d4669b1`) — "…the reciprocal monophyly of Aquila and Hieraaetus was fully supported (bootstrap value = 100)… This differed from… Haring et al. (2007), Lerner et al. (2017), and Mindell et al. (2018)… several legacy markers, particularly from the mitochondrial genome, produced topologies that conflicted with the one generated from UCE data, pointing to possible mito-nuclear discord… Our estimate for the common ancestor of modern Aquilinae (16 Mya) was similar to that by Knapp et al. (2019; ~17 Mya), but much older than the estimate by Mindell et al. (2018; ~10 Mya)."
- [13] AOS/NACC proposta 2025-C (`b5e56e8d`) — "Crested Goshawk (Accipiter trivirgatus) was recovered as sister to a large clade that included all of subfamilies Accipitrinae and Buteoninae with strong support… Catanach et al. (2024) recommended following Sangster et al. (2021) in adopting Lophospiza Kaup 1844…"
- [14] Catanach et al. 2024 (`1d4669b1`) — "…we join Sangster et al. (2021) in applying the generic name Lophospiza to [A.] trivirgatus and the Sulawesi goshawk, [A.] griseiceps… We also recommend recognizing these two species in their own subfamily (Lophospizinae)…"
- [15] Catanach et al. 2024 (`1d4669b1`) — "We recognize Lophospiza trivirgatus, which occupies a long branch that is sister to the diverse Buteoninae + Accipitrinae clade, in a monotypic subfamily… We also advocate expanding the subfamily Accipitrinae to include the genera Kaupifalco, Melierax, Urotriorchis, and Micronisus…"
- [16] Catanach et al. 2024 (`1d4669b1`) — "We recovered a sister relationship between the genera Harpagus and Microspizias… Like Oatley et al. (2015), we found the Microspizias + Harpagus clade to be sister to Buteoninae, whereas Mindell et al. (2018) found it to be sister to a clade containing Buteoninae and the rest of Accipitrinae (s.l.)… we recognize the genera Microspizias and Harpagus as members of the subfamily Harpaginae."
- [17] Mather 2021 (`c24ce4bf`) — "…A. trivirgatus, A. griseiceps and the lineage comprising A. collaris and A. superciliosus, fall outside of the Accipiter clade and are candidates for generic revision (Breman et al. 2013; Mindell et al. 2018)…"
- [18] Mather 2021 (`c24ce4bf`) — sobre *Haliaeetus*/*Icthyophaga* e "Buzzards and true hawks (Buteoninae and Accipitrinae)… reasonably closely related…"
- [19] AOS 2025-C (`b5e56e8d`) — "Revise the generic limits of Accipiter. Background: The genus Accipiter was until recently the largest genus in Accipitridae, with 46–50 species…"
- [20] Catanach et al. 2024 (`1d4669b1`) — resumo: "…ultraconserved elements… whole genome sequencing of 134 species… time-calibrated phylogeny of 237 extant or recently extinct species… includes 90% of recognized species, confirms the non-monophyly of Accipiter…"
- [21] Catanach et al. 2024 (`1d4669b1`) — "Accipiter, a catch-all genus for (usually) forest-dwelling hawks, into which ~50 species have historically been placed (Peters 1931)… Olson (1987) noted that the procoracoid foramen…"
- [22] AOS 2025-C (`b5e56e8d`) — "In 2024, WGAC voted to break Accipiter into 5 genera, following the recommendations of Catanach et al. (2024); this was followed by Clements et al. (2024) and the IOC World Bird List (Gill et al. 2025). **The AOS Checklist has not yet considered these changes to the genus Accipiter, still recognizing the expanded view of the genus.**"
- [23] AOS 2025-C (`b5e56e8d`) — "Catanach et al. (2024) combined whole genome sequence data, UCE data, as well as… Sanger… The results… largely corroborate what had been found in earlier studies, albeit with much greater support and… more complete taxon sampling…"
- [24] Catanach et al. 2024 (`1d4669b1`) — "…Circus, Megatriorchis, and Erythrotriorchis are nested within a larger 'Accipiter' clade… chewing lice (Degeeriella) obtained from the northern goshawk… and Cooper's hawk… were shown to be more closely related to lice from Circus species than from other 'Accipiter' species (Catanach and Johnson 2015)."
- [25] AOS 2025-C (`b5e56e8d`) — "…the genus Accipiter was found to form 5 distinct clades… successively sister, with other smaller genera (Circus, Megatriorchis, Erythrotriorchis, Kaupifalco, Melierax, Urotriorchis, and Micronisus) interspersed…"
- [26] Catanach et al. 2024 (`1d4669b1`) — "…we could either (i) synonymize Circus, Erythrotriorchis, and Megatriorchis with Accipiter…; or (ii) split Accipiter into four genera while retaining Circus, Kaupifalco, Erythrotriorchis, and Megatriorchis… each of these lineages diverged between 15 and 25 Mya… we opt for the second approach…"
- [27] Catanach et al. 2024 (`1d4669b1`) — "Kaupifalco monogrammicus was recovered (bootstrap value = 100) as sister to a clade containing the remainder of 'Accipiter'… This placement conflicted with Griffiths et al. (2007) and Mindell et al. (2018)…"
- [28] Catanach et al. 2024 (`1d4669b1`) — "…the subfamily Accipitrinae… was rendered polyphyletic by the embedded placement of Circus (subfamily Circinae)… our analysis placed the crested goshawk, Accipiter trivirgatus…, as sister to a large clade containing the Accipitrinae (s.l.) and Buteoninae (bootstrap value = 100)."
- [29] AOS 2025-C (`b5e56e8d`) — "…either by lumping many genera into Accipiter, including the highly distinctive harriers (Circus), or breaking Accipiter into at least 4 or 5 different genera…"
- [30] Catanach et al. 2024 (`1d4669b1`) — "…previous authors have been hesitant to revise the genus-level taxonomy… The only case in which 'Accipiter' species were reclassified in a new genus… tiny hawk… and semicollared hawk…, which are now placed by some authors in the genus Microspizias Sangster et al. (2021)."
- [31] AOS 2025-C (`b5e56e8d`) — lista de espécies de *Tachyspiza* (badius, butleri, brevipes, soloensis… rhodogaster).
- [32] Catanach et al. 2024 (`1d4669b1`) — "…[Accipiter] soloensis…, type species of the genus Tachyspiza Kaup, 1844…, and 26 other species formerly placed in Accipiter… members tend to have relatively short toes and talons… We estimated that the common ancestor of Aerospiza and Tachyspiza diverged at 17.9 Mya, whereas Mindell et al. (2018) gave an estimate of ~10 Mya."
- [33] Catanach et al. 2024 (`1d4669b1`) — "One subclade, restricted to sub-Saharan Africa, contained [Accipiter] tachiro…, type species of the genus Aerospiza Roberts, 1922…: Aerospiza tachiro, Aerospiza castanilius… and Aerospiza toussenelii…"
- [34] AOS 2025-C (`b5e56e8d`) — "These 6 species include Eurasian Sparrowhawk (Accipiter nisus)… madagascariensis, ovampensis, rufiventris, striatus, and poliogaster… The inclusion of poliogaster in this clade is provisional…"
- [35] Catanach et al. 2024 (`1d4669b1`) — "…a clade that included Accipiter nisus… and several other small, primarily bird-catching species… we tentatively place [A.] poliogaster in Accipiter…"
- [36] Catanach et al. 2024 (`1d4669b1`) — "…A. madagascariensis and A. ovampensis… within the Accipiter (s.s.) clade… supported by morphology…"
- [37] AOS 2025-C (`b5e56e8d`) — "The final clade of 'Accipiter' hawks identified by Catanach et al. (2024) was sister to Megatriorchis and Circus. The oldest available genus for this clade is Astur Lacépède 1799… This clade includes 9 species: cooperii, gundlachi, bicolor, chilensis, melanoleucus, henstii, gentilis, atricapillus, and meyerianus."
- [38] Catanach et al. 2024 (`1d4669b1`) — "Within the sister group of Accipiter (s.s.), we recovered a clade of 'Accipiter' species that was sister to Megatriorchis + Circus… we apply the generic name Astur Lacépède, 1799… the age of the split (12 Mya)… the name Cooperastur Bonaparte, 1854 is available for the clade containing Astur cooperii…"
- [39] AOS 2025-C (`b5e56e8d`) — "The harriers (Circus) should now follow American Goshawk (Astur atricapillus) in the linear sequence, following Catanach et al. (2024). Please vote on the following subproposals…"
- [40] AOS 2025-C (`b5e56e8d`) — "(d) accept the new linear placement of Circus, following Astur atricapillus. References: …"

---

### Q5 — Relações internas de Falconidae e de Strigiformes

**Pergunta:** What are the internal relationships of Falconidae and of Strigiformes?

conversation_id: `11e77f64-a157-4faf-bc0d-6843a729a1ee`

**Resposta do NotebookLM (verbatim):**

> ### Internal Relationships of Falconidae
>
> The family **Falconidae** (caracaras, forest falcons, and typical falcons) comprises 64 recognized extant species divided into three main subfamilies whose basal divergences date back to the Early–Middle Miocene [1, 2]:
>
> 1. **Herpetotherinae**: A New World endemic lineage [2] comprising two genera:
>    * **`Herpetotheres`**: Monotypic (*Herpetotheres cachinnans*, the Laughing Falcon), specialized in preying on snakes [2, 3].
>    * **`Micrastur`**: Neotropical forest falcons, which underwent radiation during the Pliocene [1, 2, 4].
> 2. **Polyborinae**: The Neotropical caracaras and allies [1, 2], featuring both diet generalists (*Phalcoboenus*, *Caracara*) and specialists (*Ibycter*, *Spiziapteryx*) [1-3].
> 3. **Falconinae**: A widespread, cosmopolitan subfamily comprising three primary genera [1, 5]:
>    * **`Microhierax`**: Five species of Oriental falconets [1, 5].
>    * **`Polihierax`**: Pygmy falcons (African and Oriental); genomic data indicate that *Polihierax* is non-monophyletic, with the Oriental *P. insignis* often placed in its own genus (*Neohierax*) [5, 6].
>    * **`Falco`**: Typical falcons (~40 species, representing ~60% of the entire family) [1, 5].
>
> * **Late Miocene Diversification**: Within *Falco*, a major burst in diversification rates occurred around **7.5 Mya** (Miocene/Pliocene transition) [1]. This rapid adaptive radiation correlates with the expansion of global C4 open grasslands/savannas, accompanied by dietary and migratory behavioral shifts [1, 6].
> * **Fossil Stem Lineages**: Fossil taxa from the Early Eocene (such as the European **Masillaraptoridae**) represent stem-group representatives of Falconiformes, establishing a long evolutionary history distinct from other diurnal birds of prey [7].
>
> ---
>
> ### Internal Relationships of Strigiformes
>
> Extant owls (**Strigiformes**) underwent a basal bifurcation during the Paleogene/Oligocene, splitting into two distinct families [8-10]:
>
> 1. **Tytonidae (Barn Owls and Bay Owls)**:
>    * Characterized by specialized osteological traits (such as a thickened distal tarsometatarsus) [11].
>    * Split during the Oligocene into two extant subfamilies: **Tytoninae** (barn owls, *Tyto*) and **Phodilinae** (bay owls, *Phodilus*) [8, 9, 12].
>
> 2. **Strigidae (Typical Owls)**:
>    * Represents the vast majority of modern owl diversity (~194 species), having undergone extensive adaptive radiation starting in the Miocene [10, 13].
>    * Recent phylogenomic analyses using ultraconserved elements (UCEs) organize Strigidae into two major sister clades [14-17]:
>      * **Clade A (Surniinae / Ieraglaucinae)**: Includes hawk-owls and pygmy-owls divided into major tribes such as **Ninoxini** (*Ninox*, *Sceloglaux*, *Uroglaux*) and **Surniini** (*Glaucidium*, *Athene*, *Aegolius*, *Surnia*, *Micrathene*, *Xenoglaux*, *Margarobyas*) [14, 16, 18].
>      * **Clade B (Striginae)**: Typical screech owls, horned owls, and wood owls, divided into tribes **Otini** (*Otus*, *Pyrroglaux*), **Asionini** (*Asio*, *Ptilopsis*, *Nesasio*, *Pseudoscops*), **Bubonini** (*Bubo*, *Ketupa*, *Scotopelia*, *Jubula*), **Megascopini** (*Megascops*, *Psiloscops*), **Pulsatrigini**, and **Strigini** (*Strix*, *Ciccaba*) [15, 16, 19].
>
> * **Morphological Convergence & Taxonomic Paraphyly**:
>   * Phylogenomic studies reveal that traditional morphology has obscured owl systematics [20-22]. Multiple traditional genera are **paraphyletic** (*Athene*, *Otus*, *Asio*, *Megascops*, *Bubo*, *Strix*) or **polyphyletic** (*Ninox*, *Glaucidium*) [20, 22, 23].
>   * This mismatch is largely caused by rapid adaptation to island environments ("island syndrome"), where colonizing species evolved distinct body sizes and leg proportions, causing them to be misclassified into separate monotypic genera (e.g., *Sceloglaux*, *Nesasio*, *Pseudoscops*, *Margarobyas*, *Uroglaux*) [21].
>   * Modern taxonomic proposals recommend subsuming several island/monotypic genera into larger clades (e.g., merging *Ciccaba* into *Strix*, *Ketupa* and *Scotopelia* into *Bubo*, *Mimizuku* and *Pyrroglaux* into *Otus*, and *Nesasio* and *Pseudoscops* into *Asio*) [23].
>
> * **Stem Fossil Owls**:
>   * Early Cenozoic stem-group owls include fossil families such as **Protostrigidae**, **Ogygoptyngidae**, **Palaeoglaucidae**, and **Sophiornithidae** [24-26].
>   * Eocene fossil owls like ***Primoptynx*** demonstrate an ancestral, accipitrid-like prey-killing mechanism (using enlarged talons to subdue large mammals) before modern owl perching and grasping tarsometatarsal modifications fully evolved [26, 27].
>
> 🦉 *Would you like me to synthesize these phylogenetic findings into a structured comparison table or draft a detailed report on the evolutionary history of Telluraves?*

**Citações (número → fonte; trecho citado, truncado):**

- [1] Fuchs et al. 2015 (ScienceDirect, só resumo/prévia, `f82bdd8e`) — "We analyzed eight independent loci for 1–7 individuals from each of the 64 currently recognized Falconidae species, together with two fossil falconid temporal calibrations… similar diversification ages in the Early to Middle Miocene for the three traditional subfamilies, Herpetotherinae, Polyborinae and Falconinae… divergences within the subfamily Falconinae began about 16 mya and divergences within… Falco, including about 60% of all Falconidae species, began about 7.5 mya… significant increase in diversification rate at the basal phylogenetic node for the genus Falco… correlates generally with expansion of C4 grasslands…"
- [2] Fuchs et al. 2015 (`f82bdd8e`) — "The Falconidae… began to diversify during the late Oligocene–early Miocene (Fuchs et al., 2011a, 2012). Sixty-four species, divided into eleven genera… Three subfamilies (Herpetotherinae, Polyborinae and Falconinae) are widely recognized and their relationships and generic limits are well established… The former two subfamilies are endemic to the New World…"
- [3] Fuchs et al. 2015 (`f82bdd8e`) — "New Falconidae classification: Subfamily Herpetotherinae — Herpetotheres…, Micrastur… Subfamily Polyborinae — Spiziapteryx…, Caracara…"
- [4] Fuchs et al. 2015 (`f82bdd8e`) — lista de referências: "J. Fuchs et al. Pliocene diversification within the South American forest-falcons (Falconidae: Micrastur) Mol. Phylogenet. Evol. (2011)…"
- [5] Fuchs et al. 2015 (`f82bdd8e`) — "The Falconinae includes three genera… Polihierax (two species…) and Microhierax (five species, Oriental), whereas the genus Falco (39 species) is cosmopolitan… whether the Oriental species of Polihierax (P. insignis) merits its own genus (Neohierax)… kestrels (17 species) are often considered as the ancestral lineage…"
- [6] Fuchs et al. 2015 (`f82bdd8e`) — destaques: "Polihierax is not monophyletic. Falcons diversified at the end of the Miocene. Falco experienced an increase in diversification rate… Species inhabiting open habitat and migratory species have higher speciation rates."
- [7] Mayr & Kitchener 2022 (`53dae551`) — introdução: "Falcons… were long considered to be closely related to other diurnal birds of prey. However, most current molecular analyses recover them in a clade together with seriemas…, parrots…, and passerines…"
- [8] Dalsätt, tese DiVA (`6899a064`) — "During the Oligocene (35.5 – 23.5 Ma) the family Tytonidae split into the two extant subfamilies Tytoninae and Phodilinae (Mourer-Chauvirè 1987)…"
- [9] Livezey & Zusi 2007 (`341dae82`) — "…the present analysis reaffirmed a basal bifurcation between barn-owls (Tytonidae) and typical owls (Strigidae), with the former including Phodilus…"
- [10] Salter et al. 2020 (`1e8f2a66`) — "Early studies of owls separated them into 2 taxonomic groups: the barn owls (Tytonidae…; 19 species…) and the typical owls (Strigidae…; 194 species…)…"
- [11] Dalsätt (`6899a064`) — "…Morphologically the new owl is most similar to the genus Prosybris (Tytonidae). A general characteristic for Tytonidae is the thick distal part…"
- [12] Livezey & Zusi 2007 (`341dae82`) — classificação: "Order Strigiformes… Family Tytonidae… Subfamily Tytoninae… Subfamily Phodilinae… Family Strigidae…"
- [13] Dalsätt (`6899a064`) — "…Strigidae, is not known from the Paleogene and began to radiate extensively only in the Miocene (23.5 – 5.2 Ma)…"
- [14] Salter et al. 2020 (`1e8f2a66`) — "These concatenated analyses resolve barn owls and typical owls as reciprocally monophyletic and suggest that typical owls comprise 2 sister clades… Clades A and B… Clade A comprises (1) a small Australasian group, tribe Ninoxini…, containing Uroglaux dimorpha…"
- [15] Salter et al. 2020 (`1e8f2a66`) — "Clade B corresponds largely to a group of taxa known as Striginae (Wink et al. 2008)… a clade of Otus relatives (tribe Otini, including Pyrroglaux podargina…) sister to the remainder… a clade of Ptilopsis granti + Asio relatives (tribe Asionini…)… 4,235 nuclear ultraconserved element loci… All nodes received 100% bootstrap support / BPP unless otherwise labeled…"
- [16] Salter et al. 2020 (`1e8f2a66`) — "…Clade A… generally corresponds to subfamilies Surniinae + Ieraglaucinae… Our coalescent-based, SVDQuartets analysis makes interpreting exact membership of Clade A more complicated, however, because we resolved Uroglaux dimorpha as a representative of a group that is sister to Clades A+B… neither concatenated nor coalescent-based analyses… recover Clade C (Margarobyas, Ninox, and Sceloglaux) as sister to Clades A+B (Wink et al. 2008, 2009) or Clade A as sister to Clades B+C (Wink 2016)…"
- [17] Salter et al. 2020 (`1e8f2a66`) — "…sampling thousands of nuclear loci from all but one of the 28 described genera… typical owls comprise 2 major sister clades that generally correspond to subfamilies Surniiniae… and Striginae… as many as 8 typical owl genera are not monophyletic…"
- [18] Salter et al. 2020 (`1e8f2a66`) — "…Uroglaux dimorpha… as sister to a clade of Sceloglaux albifacies… and Ninox boobook…, and (2) a larger clade that includes… Margarobyas lawrencii of Cuba as sister to… Surniinae… Xenoglaux loweryi belongs within Surniinae… Glaucidium brodiei is…"
- [19] Salter et al. 2020 (`1e8f2a66`) — "…sister to a large clade containing species in tribes Bubonini, Megascopini, Pulsatrigini, and Strigini… Jubula lettii as sister to tribe Bubonini… Margarobyas lawrencii within Surniinae rather than tribe Otini…"
- [20] Salter et al. 2020 (`1e8f2a66`) — resumo: "…suggest the genera Athene, Otus, Asio, Megascops, Bubo, and Strix are paraphyletic, whereas Ninox and Glaucidium are polyphyletic… both divergence and convergence in morphological characters have obscured our understanding… particularly those with insular distributions."
- [21] Salter et al. 2020 (`1e8f2a66`) — "…of the 12 monotypic genera of typical owls…, 7 of these species have insular distributions (Sceloglaux, Pyrroglaux, Mimizuku, Nesasio, Pseudoscops, Margarobyas, and Uroglaux)… at least 5 of these species (Sceloglaux, Pyrroglaux, Mimizuku, Nesasio, and Pseudoscops) are actually nested within larger genera… 'island syndrome'… larger body size, longer legs, and reduced flight morphology."
- [22] Salter et al. 2020 (`1e8f2a66`) — "…affecting the taxonomy of as many as 177 of the 194 described species of typical owls…"
- [23] Salter et al. 2020 (`1e8f2a66`) — recomendações taxonômicas: *Ninox jacquinoti* → *Athene*; *Sceloglaux* → *Ninox*; *Pyrroglaux* e *Mimizuku* → *Otus*; *Nesasio* e *Pseudoscops* → *Asio*; *Scotopelia* e *Ketupa* → *Bubo*; *Ciccaba* → *Strix*; *Megascops nudipes* → *Psiloscops*; *Glaucidium brodiei* → *Taenioptynx*.
- [24] Dalsätt (`6899a064`) — "…the fossil family Sophiornithidae (Mourer-Chauviré 1987)… Eocene (56.5 – 35.5 Ma). Other families to which fossil strigiforms from Eocene deposits have been assigned are the Protostrigidae…, Palaeglaucidae…, and the extant Tytonidae…"
- [25] Phylogenetic definitions (`1428d7cb`) — "The clade also includes the extinct Teratornithidae…, Sandcoleidae…, the extinct Ogygoptyngidae and Protostrigidae, and other groups related to owls (Mayr 2009, 2017, 2022)…"
- [26] Mayr et al. 2020, *Primoptynx* (`2b83db39`) — "Berruornis, Ogygoptynx, and the species assigned to the Protostrigidae are outside crown group Strigiformes…, but the phylogenetic affinities of early Cenozoic Strigiformes are poorly resolved and have not yet been subjected to a formal phylogenetic analysis. Mourer-Chauviré (1994) hypothesized that… the tarsometatarsal trochleae… gradually became more plantarly deflected…"
- [27] Mayr et al. 2020 (`2b83db39`) — "Primoptynx and protostrigid owls appear to have been specialized in foraging on prey items that required an accipitridlike killing strategy, such as larger-sized or more defensive mammals. Extinction of these peculiar owls could have been related to the radiation of accipitrid diurnal birds of prey, which commenced toward the late Eocene and early Oligocene (Mayr, 2009)…"

---

### Q6 — Fósseis mais antigos, idades e posição de grupo-tronco

**Pergunta:** What are the oldest fossils of each raptor lineage, their ages, and stem-group placements (Masillaraptoridae, Messelasturidae, stem owls, phorusrhacids)?

conversation_id: `f66d98be-076f-4f97-b9c3-38ac71cc9e40`

**Resposta do NotebookLM (verbatim):**

> Here is a breakdown of the oldest known fossils, ages, and evolutionary stem-group placements for each of the specified raptorial lineages:
>
> ---
>
> ### **1. Masillaraptoridae**
> * **Oldest Fossils & Taxa:** The oldest well-preserved masillaraptorid specimens include ***Danielsraptor phorusrhacoides*** and ***Masillaraptor* cf. *M. parvunguis*** from the London Clay Formation at Walton-on-the-Naze (Essex, UK) [1-3].
> * **Age:** Early Eocene (**~55–54.6 Ma**) [4, 5].
> * **Stem-Group Placement:** Stem-group representatives of **Falconiformes** (falcons) [1, 6, 7]. Despite their falconid relationship, they display remarkable morphological convergences in their deep beaks and long tarsometatarsi with the Cariamiformes (specifically the Phorusrhacidae) [1, 6].
>
> ---
>
> ### **2. Messelasturidae**
> * **Oldest Fossils & Taxa:** The genus ***Tynskya*** contains the oldest messelasturid remains, represented by ***Tynskya waltonensis*** from the London Clay (UK) and ***Tynskya eocaena*** from the Green River Formation (Wyoming, USA) [5, 8, 9].
> * **Age:** Early Eocene (**~55 Ma** for *T. waltonensis*; **~52 Ma** for *T. eocaena*) [5, 8].
> * **Stem-Group Placement:** Stem-group **Psittacopasseres** (the clade uniting parrots and passerines) or stem **Psittaciformes** / sister group to Halcyornithidae [8, 10-12]. Messelasturids represent an "owl/parrot mosaic," combining raptor-like traits (hooked beaks, long supraorbital processes, specialized claws) with semi-zygodactyl feet [8, 13-15].
>
> ---
>
> ### **3. Stem Owls (Strigiformes)**
> * **Oldest Fossils & Taxa:** The earliest confirmed owl fossils are ***Berruornis halbedeli*** from Germany and ***Ogygoptynx wetmorei*** (Ogygoptyngidae) from Colorado, USA (alongside *Berruornis orbisantiqui* from France) [16-18].
> * **Age:** Middle to Late Paleocene (**~61 Ma** for *B. halbedeli*; **~60–59 Ma** for *O. wetmorei*) [16, 17].
> * **Stem-Group Placement:** Stem-group **Strigiformes** [19, 20]. These Paleogene taxa branch outside crown group Strigiformes (which is formed by the extant Tytonidae and Strigidae) [19, 20]. Early Eocene stem owls like *Primoptynx poliotauros* and protostrigids (*Minerva*, *Eostrix*) retained unreduced supraorbital processes and raptorial talon proportions similar to modern hawks [21-24].
>
> ---
>
> ### **4. Phorusrhacidae ("Terror Birds")**
> * **Oldest Fossils & Taxa:** The earliest phorusrhacid fossils belong to a predominantly South American radiation of flightless apex terrestrial predators [25, 26].
> * **Age:** Middle Paleocene (**~60 Ma**) through to the lower Pleistocene [26].
> * **Stem-Group Placement:** Clade within **Cariamiformes**, with extant seriemas (**Cariamidae**) confirmed as their closest living relatives [25-27].
>
> ---
>
> 💡 *Would you like to explore the ecological adaptations or foot morphology trade-offs (such as semizygodactyly vs. anisodactyly) across these Paleogene raptorial groups?*

**Citações (número → fonte; trecho citado, truncado):**

- [1] Mayr & Kitchener 2022 (`53dae551`) — resumo: "The Eocene taxon Masillaraptoridae includes long-legged, raptorial birds, the phylogenetic affinities of which are poorly resolved… fossils from the London Clay of Walton-on-the-Naze… corroborate the hypothesis that masillaraptorids are stem group representatives of the Falconiformes… Danielsraptor phorusrhacoides, gen. et sp. nov.… a smaller species closely resembles Masillaraptor parvunguis… the long and deep beak of masillaraptorids is remarkably similar to that of the cariamiform Phorusrhacidae… These similarities may be plesiomorphic for a more inclusive clade containing falconiforms and cariamiforms, but it is also possible that masillaraptorids represent a direct evolutionary link between the Falconiformes and Cariamiformes."
- [2] Mayr & Kitchener 2022 (`53dae551`) — "Arguably the best-preserved fossils of Masillaraptor-like birds come from the lower Eocene London Clay of Walton-on-the-Naze… found… by the late Michael Daniels…"
- [3] Mayr & Kitchener 2022 (`53dae551`) — material referido de *Masillaraptor* cf. *M. parvunguis* (NMS.Z.2021.40.14, .15) e legenda da Fig. 6.
- [4] Mayr & Kitchener 2022 (`53dae551`) — "The strata of the London Clay at Walton-on-the-Naze are from division A2… and date back to 55–54.6 Ma (Collinson et al. 2016). The London Clay exposed at the Isle of Sheppey, which yielded the holotype of Stintonornis mitchelli, is from the younger divisions C–F (51.5–53.5 Ma…). The holotype of Antarctoboenus carlinii has an age of 52.8–49 million years…"
- [5] Mayr 2021, PalZ (`a681b573`) — "…An assignment of the ca. 55 million years old London Clay fossil to the taxon Tynsyka is also in agreement with its closer temporal proximity to the 52 myo strata of the Green River Formation (the type locality of T. eocaena…) than to the 48 myo Messel oilshale (the type locality of M. gratulator)."
- [6] Mayr & Kitchener 2022 (`53dae551`) — "Masillaraptor parvunguis was considered to be of uncertain affinities when first described (Mayr 2006), but subsequently masillaraptorids were regarded as stem group representatives of the Falconiformes (Mayr 2009, 2022). The published character evidence for this placement is still weak but the new fossils exhibit previously unknown characters… Daniels (1995) proposed… affinities to the cariamiform Phorusrhacidae. Even though this hypothesis is not confirmed herein…"
- [7] Mayr & Kitchener 2022 (`53dae551`) — "Strict consensus tree of the two most parsimonious trees (L = 64, CI = 0.67, RI = 0.73)… Falconiform affinities of masillaraptorids are also suggested by… Stintonornis mitchelli… and… Antarctoboenus carlinii…"
- [8] Mayr 2021, PalZ (`a681b573`) — resumo: "Tynskya eocaena is an early Eocene bird with a raptor-like skull and semi-zygodactyl feet… T. waltonensis… from the London Clay… Phylogenetic analyses… did not conclusively resolve the higher-level affinities of messelasturids and the closely related halcyornithids… An analysis that was constrained to a molecular scaffold, however, recovered messelasturids as the sister taxon of a clade including psittaciform and passeriform birds…"
- [9] Mayr 2021 (`a681b573`) — "The first formally described messelasturid fossils are two skulls of Messelastur gratulator from the early Eocene of Messel… originally considered to be from a hawk-like diurnal bird of prey (Accipitridae)…"
- [10] Mayr 2021 (`a681b573`) — "…a phylogenetic analysis of 77 morphological characters… In a second analysis, 20 'dummy' characters were added to constrain the topology of the extant taxa according to… molecular analyses…"
- [11] Mayr 2021 (`a681b573`) — "…The present analyses do not support closer affinities between messelasturids and strigiform birds… If messelasturids are stem group representatives of this latter clade [Psittaciformes + Passeriformes]…, some of their raptor-like features… may be plesiomorphic. At present, however, this assumption remains mere speculation, not least because close affinities between messelasturids and the clade including Psittaciformes and Passeriformes are far from being well-established."
- [12] Psittacomimidae fam. nov., *Vertebrate Zoology* (`5102aa83`) — "…Messelastur was initially assigned to the Accipitridae…, whereas Tynskya was likened to the Strigiformes… messelasturids were subsequently considered to be most closely related to the Halcyornithidae and both taxa were tentatively identified as stem group representatives of the Psittaciformes (Mayr 2011b). However… the placement of messelasturids is poorly resolved in current analyses, some of which do not even support psittacopasserine affinities…"
- [13] Mayr 2021 (`a681b573`) — "Early Cenozoic fossils often exhibit mosaic character distributions… Messelasturidae and Halcyornithidae… combine derived characteristics of owls…, diurnal birds of prey…, and parrots…"
- [14] Mayr 2021 (`a681b573`) — "…a phylogenetic analysis supported a sister group relationship between Messelasturidae and Strigiformes, with the clade including both taxa being nested amongst diurnal birds of prey (Mayr 2005)… raptor-like traits, such as a short beak and long supraorbital processes, with a 'parrot-like' tarsometatarsus morphology indicating at least semi-zygodactyl feet…"
- [15] Mayr 2021 (`a681b573`) — resultados da análise: 14 árvores mais parcimoniosas (L = 232; CI 0.36; RI 0.60); consenso pouco resolvido; Tynskya + Messelastur irmãos; Messelasturidae + Halcyornithidae (não mantido no bootstrap).
- [16] Dalsätt, tese (`6899a064`) — "…the Strigiformes has a long evolutionary history, spanning at least 60 mya… The oldest confirmed fossil strigiform is the Paleocene Ogygoptynx wetmorei (65 – 56.5 Ma) described from a tarsometatarsus collected in Colorado (Rich and Bohaska 1976)…"
- [17] Mayr et al. 2020, *Primoptynx* (`2b83db39`) — "Ogygoptynx wetmorei…, from the mid-Tiffanian (60–59 Ma) of Colorado… only known from a single tarsometatarsus. Berruornis includes two species, B. orbisantiqui…, from the upper Paleocene (late Thanetian; ∼57.5 Ma) of France…, and B. halbedeli Mayr, 2002, from the middle Paleocene (∼61 Ma; Mayr and Smith, 2019) of Germany… Berruornis was assigned to the taxon Sophiornithidae…"
- [18] Mayr et al. 2020 (`2b83db39`) — lista de referências (Rich & Bohaska 1976 "The world's oldest owl…"; 1981 "The Ogygoptyngidae…").
- [19] Mayr et al. 2020 (`2b83db39`) — "Owls (Strigiformes) have an extensive fossil record, beginning in the mid-Paleocene… most early Cenozoic taxa are represented by only a few bones…"
- [20] Mayr et al. 2020 (`2b83db39`) — "Berruornis, Ogygoptynx, and the species assigned to the Protostrigidae are outside crown group Strigiformes…, but the phylogenetic affinities of early Cenozoic Strigiformes are poorly resolved and have not yet been subjected to a formal phylogenetic analysis…"
- [21] Mayr et al. 2020 (`2b83db39`) — resumo: "…As in the protostrigid taxon Minerva from the late early to early middle Eocene…, the ungual phalanges of the hallux and the second toe of the new species are distinctly larger… Even though the new species may well be a stem group representative of protostrigid owls, current data do not allow an unambiguous phylogenetic placement… the feet of P. poliotauros correspond to those of extant hawks and allies (Accipitridae)…"
- [22] Mayr et al. 2020 (`2b83db39`) — legenda da Fig. 5 (tarsometatarsos de *Primoptynx*, *Minerva*, *Ogygoptynx*, *Tyto*, *Eostrix gulottai*).
- [23] Mayr et al. 2020 (`2b83db39`) — "Primoptynx poliotauros… and M. antiqua share a distinctive morphology of their feet, in which the ungual phalanges of the first and second toes are greatly enlarged. **The pedal phalanges are unknown for Berruornis, Ogygoptynx, and the species of Eostrix**, but in Palaeoglaux artophoron from Messel the ungual phalanges of the hallux and the second toe are not larger than those of the other toes… the new taxon is therefore clearly outside the clade including Minerva and Oligostrix…"
- [24] Mayr et al. 2020 (`2b83db39`) — "…Fossils from the London Clay (Mayr 2017) and the Bridger Formation (Fowler et al. 2018) show that the skulls of early and middle Eocene owls possessed well-developed supraorbital processes… likely to be plesiomorphic for strigiform birds and suggest that early Cenozoic owl were less nocturnal than their living relatives…"
- [25] Birds, Diversification of (ResearchGate, `d3d5a7ff`) — "…the mostly flightless **Eocene–Pleistocene** Phorusrhacidae ('terror birds'; Alvarenga and Höfling 2003; Angst and Buffetaut 2017) among Cariamiformes, both of which are known from the Americas…"
- [26] Degrange et al. 2010 (`b9c81379`) — "Phorusrhacids were a predominantly South American radiation of gruiform birds known from middle Paleocene to lower Pleistocene deposits… and most closely related to seriemas (Cariamidae) among extant birds… Andalgalornis steulleti…, from the upper Miocene–lower Pliocene (≈6 million years ago) of Argentina…"
- [27] Mayr & Kitchener 2022 (`53dae551`) — métodos da análise filogenética (matriz de Mayr 2016, 42 caracteres; Phorusrhacidae representado por *Llallawavis*, *Psilopterus*, *Patagornis*).

---

### Q7 — O ancestral comum de Telluraves (ou Australaves) era um predador raptorial?

**Pergunta:** Was the last common ancestor of Telluraves (or Australaves) a raptorial predator? Evidence for and against.

conversation_id: `be1d178e-a08b-4c5d-98d5-a2eeba072714`

**Resposta do NotebookLM (verbatim):**

> Whether the last common ancestor of **Telluraves** (core landbirds) or **Australaves** was a raptorial predator is one of the most intriguing and actively debated questions in modern avian phylogenomics and paleontology [1-4].
>
> Traditional taxonomies grouped diurnal birds of prey together, but modern phylogenomic studies have demonstrated that raptorial birds are polyphyletic, scattered across the landbird tree [1, 5, 6]. This dispersal led to two competing hypotheses: a **"raptorial grade" hypothesis** (where the common ancestor was an apex predator, and non-raptorial groups lost these traits) versus an **independent convergence hypothesis** (where raptorial lifestyles evolved multiple times independently) [2-4, 7].
>
> ---
>
> ### Evidence FOR a Raptorial Common Ancestor
>
> 1. **Phylogenomic "Raptorial Grade" Topology**:
>    * Across both main branches of Telluraves (**Afroaves** and **Australaves**), predatory or raptorial lineages occupy the deepest, earliest-diverging nodes [1, 3, 5, 8].
>    * In **Afroaves**, Accipitriformes (hawks, eagles, vultures) and Strigiformes (owls) diverge near the base [1, 5, 9].
>    * In **Australaves**, Cariamiformes (seriemas and extinct giant apex predator "terror birds") diverges first, followed by Falconiformes (falcons) as the sister group to Psittacopasseres (parrots and passerines) [3, 8, 10].
>    * Whole-genome phylogenomic analyses (e.g., Jarvis et al. 2014, Prum et al. 2015) argue that the presence of raptorial lineages at the deepest divergences of both major clades suggests the last common ancestor of Telluraves was an apex predator, followed by gradual losses of raptorial traits in non-predatory lineages (such as mousebirds, cavity nesters, parrots, and songbirds) [3, 4].
>
> 2. **Fossil "Mosaics" with Raptor-Like Features**:
>    * Stem members of non-raptorial groups in the Paleogene fossil record retain raptor-like anatomical traits [3, 11, 12].
>    * The early Eocene **Messelasturidae** (e.g., *Messelastur*, *Tynskya*), recovered in several analyses as stem group representatives of Psittacopasseres (parrots + passerines), combine a hooked "raptor-like" beak, long supraorbital processes, and curved ungual talons with semi-zygodactyl "parrot-like" feet [11-15].
>    * If Messelasturidae represent stem-psittacopassereans, their raptor-like skull and claw morphology suggest that predatory features are plesiomorphic (ancestral) for Australaves [14, 16].
>    * Similarly, early Eocene stem mousebirds and cuckoo-roller relatives in Afroaves exhibited raptor-like foot structures or specialized predatory behaviors [3].
>
> 3. **Stronger Case for Australaves Specifically**:
>    * Because both Cariamiformes and Falconiformes form successive sister lineages to parrots and songbirds, a raptorial ancestor for **Australaves** (or at least Eufalconimorphae) is considered structurally parsimonious by several authors, requiring only a single loss of raptoriality in the common ancestor of parrots and passerines [3, 10, 16].
>
> ---
>
> ### Evidence AGAINST a Raptorial Common Ancestor
>
> 1. **High Propensity for Convergent Evolution (Homoplasy)**:
>    * Predatory adaptations—such as hooked beaks, hypertrophied pedal claws, and specialized flight musculature—are classic examples of convergent evolution driven by biomechanical demands [17, 18].
>    * Opponents argue that it is morphologically and ecologically simpler for raptorial traits to evolve independently in a few predatory lineages than for complex carnivorous and hunting adaptations to be lost independently across dozens of diverse non-raptorial landbird families (e.g., woodpeckers, toucans, hornbills, mousebirds, songbirds, parrots) [2, 3].
>
> 2. **Early Non-Raptorial Fossils**:
>    * The oldest well-corroborated fossil representative of Telluraves discovered to date is ***Tsidiiyazhi*** (early Paleocene, ~62 Ma) from New Mexico, which is a stem mousebird [19, 20]. *Tsidiiyazhi* was a small, non-raptorial arboreal bird, showing that non-predatory landbird niches were already established at the very beginning of the telluravian radiation following the K–Pg extinction [19-21].
>
> 3. **Sensitivity to Phylogenomic Data Types & Analytical Artifacts**:
>    * Topologies placing raptors in a paraphyletic grade at the base of Telluraves (e.g., Accipitriformes sister to all other landbirds) often depend heavily on the specific dataset used (coding exons vs. introns, UCEs, or retrotransposon insertions) and whether concatenation or coalescent species-tree models are applied [7, 22-24].
>    * Analyses controlling for incomplete lineage sorting (ILS) or retroelement insertions frequently resolve a clean split between Afroaves and Australaves, making independent gains of raptoriality within each clade a highly supported alternative [7, 23-25].
>
> ---
>
> ### Summary
>
> * **For Australaves**: There is **moderate-to-strong support** for a raptorial ancestor, given that Cariamiformes and Falconiformes branch sequentially at the base of the clade and early fossil stem-parrots (*Messelasturidae*) retain raptorial skulls and talons [3, 10, 11, 14, 16].
> * **For Telluraves as a whole**: The hypothesis of a raptorial ancestral landbird remains **plausible but unresolved**. While phylogenomic topologies placing raptor lineages at deep nodes support an ancestral apex predator, extensive morphological convergence and non-raptorial early Cenozoic fossils keep the alternative hypothesis—multiple independent origins of raptoriality—equally prominent [2-4, 7, 19].
>
> 🦅 Would you like to explore the specific morphological adaptations of Eocene fossil groups like the Messelasturidae or Masillaraptoridae in more detail?

**Citações (número → fonte; trecho citado, truncado):**

- [1] Birds, Diversification of (ResearchGate, `d3d5a7ff`) — "(8) Telluraves… Among these are several primarily raptorial (vertebrate-eating) lineages, namely Accipitriformes…, Strigiformes…, Cariamiformes…, and Falconiformes… Phylogenomic analyses consistently find that these groups are not each other's closest relatives within…"
- [2] idem (`d3d5a7ff`) — "…Telluraves, which may suggest that the raptorial lifestyle evolved multiple times. However, the possibility that raptoriality is ancestral to Telluraves has also been proposed (Jarvis et al. 2014; Prum et al. 2015; McClure et al. 2019)…"
- [3] Jarvis et al. 2014 (`c92cd31c`) — "The non-monophyly of the birds of prey at the deepest branches of the Australaves and Afroaves radiations suggests that the common ancestor of core landbirds may have been an apex predator, followed by two losses of the raptorial trait. Seriema at the deepest branch of Australaves could be considered to belong to a raptorial taxon because they kill vertebrate prey… and are the sole living relatives of the extinct giant 'terror birds'… The deepest branches after Accipitriformes and owl among the Afroaves, the mousebirds and cuckoo-roller, have Eocene relatives with raptor-like feet…, and the cuckoo-roller specializes on chameleon prey… This suggests that losses of the predatory phenotype were gradual across successive divergences of each of the two core landbird radiations."
- [4] Prum et al. 2015 (`74f2801f`) — "Our findings also support the emerging pattern that landbirds evolved from a raptorial grade. The sister group relationships of hawks to the rest of the landbirds, of owls to the diverse coraciimorph clade, and of seriemas and falcons to the parrots and passerines indicate the persistence of a raptorial ecology among ancestral landbirds."
- [5] Resolving the Avian Tree (`47d8810c`) — "The Jarvis et al. (2014) TENT, like a number of prior analyses (Ericson 2006; Hackett 2008; Kimball 2013), splits core landbirds into two clades that Ericson (2012) named Australaves and Afroaves… That split renders the raptorial lineages para- or polyphyletic…"
- [6] Catanach et al. 2024 (`1d4669b1`) — "…Sibley (1960) proposed that 'Falconiformes [sensu lato] may be polyphyletic'… confirmed repeatedly with molecular data… falcons are more closely related to perching birds… and parrots… than to hawks… Degeeriella… the clade parasitizing falcons is more closely related to the genus Picicola…"
- [7] Braun & Kimball 2021 (`b0e43963`) — "The hypothesis that most raptorial landbirds form a paraphyletic assemblage at the base of core landbirds… can probably be rejected. Instead, the majority of recent analyses favor a topology that divides raptors (and core landbirds overall) into two major clades (Australaves and Afroaves), **each of which has two raptorial lineages at their base**. However, two alternatives… One splits raptors into four independent lineages, placing hawks, eagles, and New World vultures sister to the remaining core landbirds… The third possibility divides raptors into three major clades, a 'core' raptorial clade at the base of landbirds and separate lineages for seriemas and falcons within Australaves… it seems inappropriate to reject any of these three hypotheses at this time. Additionally, the position of mousebirds is also unclear…"
- [8] Birds, Diversification of (`d3d5a7ff`) — "Psittacopasseres, which is closely related to Cariamiformes and Falconiformes, forming the group Australaves…"
- [9] Wu et al. 2024 (`cda64d61`) — Coraciimorphae basal; Hieraves irmão de Australaves (mesmo trecho de Q1 [5]).
- [10] Psittacomimidae fam. nov. (`5102aa83`) — Psittacopasseres, Falconiformes como irmão, Australaves (mesmo trecho de Q2 [1]).
- [11] Mayr 2021, PalZ (`a681b573`) — "…Mayr (2000, 2005, 2011) emphasized the mosaic nature of the skeletal morphology of Tynskya, which combines several 'raptor-like' traits, such as a short beak and long supraorbital processes, with a 'parrot-like' tarsometatarsus morphology indicating at least semi-zygodactyl feet…"
- [12] Mather 2021 (`c24ce4bf`) — "…Horusornithidae, thought to be an extinct lineage of Accipitriformes (Mourer-Chauviré 1991), and the Messelasturidae, predators in a lineage now thought to be related to stem-group parrots (Mayr 2006a; Mayr 2011). The Eocene raptor Masillaraptor parvunguis… was later considered to be a falconid (Mayr 2009, 2017)."
- [13] Mayr 2021 (`a681b573`) — resumo (mesmo trecho de Q6 [8]).
- [14] Mayr 2021 (`a681b573`) — "…If messelasturids are stem group representatives of this latter clade…, some of their raptor-like features… may be plesiomorphic. At present, however, this assumption remains mere speculation, not least because close affinities between messelasturids and the clade including Psittaciformes and Passeriformes are far from being well-established."
- [15] Mather 2021 (`c24ce4bf`) — "…Messelastur gratulator Peters, 1994. It was described as a member of Accipitridae based on morphological characters of the skull. However, the discovery of a complete skeleton resulted in this species being placed within its own family, Messelasturidae (Mayr 2005)…"
- [16] Psittacomimidae fam. nov. (`5102aa83`) — "If the Falconiformes are the extant sister taxon of the Psittacopasseres, the stem species of the latter may have been raptorial and the similarities of the Messelasturiformes to the Falconiformes probably are plesiomorphic. However, further data on critical aspects of the skeleton are needed… **Moreover, a monophyletic Eufalconimorphae is not recovered in all sequence-based analyses (Stiller et al. 2024: extended data, fig. 4), in which case the raptor-like features of the Messelasturiformes would have to be judged differently.**"
- [17] Ksepka et al. 2017 (`b84ec51f`) — "Enforcing our preferred molecular scaffold… results in six MPTs… 31 steps longer… This highlights a potentially high level of homoplasy within some morphological characters, including many related to beak shape and pedal digit organization that contribute to the clustering of the raptorial taxa in the morphological tree…"
- [18] Catanach et al. 2024 (`1d4669b1`) — "…ornithologists have struggled to resolve relationships within Accipitridae because of high levels of apparent convergence in external morphology and skeletal structure, which might be caused by shared selective pressures that are inherent to a predatory lifestyle (Holdaway 1994, Pecsics et al. 2019)."
- [19] Birds, Diversification of (`d3d5a7ff`) — "…a putative member of the mousebird lineage, Tsidiiyazhi from the Paleocene of North America, is the oldest well-corroborated example of a telluravian yet discovered…"
- [20] Ksepka et al. 2017 (`b84ec51f`) — "All phylogenetic analyses using the full taxon sample recovered Tsidiiyazhi as a stem mousebird… the morphology-based tree conflicts with recent molecular results in uniting Passeriformes with a Piciformes + Upupiformes clade and supporting a 'raptorial clade' (Accipitriformes, Falconiformes, and Strigiformes)… mousebirds have proven highly labile…"
- [21] Tip-dating and the Origin of Telluraves (ResearchGate, só resumo, `79ba75e1`) — "…Our analyses recover a Late Cretaceous origin for crown Telluraves, with a few lineages crossing the K-Pg boundary. Following the K-Pg boundary, our results show the group underwent rapid diversification… We find very little confidence for the precise topological placement of many extinct taxa…"
- [22] Resolving the Avian Tree (`47d8810c`) — legenda das quatro topologias de core landbirds (mesmo trecho de Q1 [1]).
- [23] Braun & Kimball 2021 (`b0e43963`) — "…We recovered the 'raptors sister' topology (i.e., a clade of all non-raptorial landbirds nested within a paraphyletic assemblage of raptors) in the unpartitioned analysis of coding exons… In contrast, all other analyses of nucleotide data for the complete dataset and the coding exons placed Accipitrimorphae… sister to all other landbirds, placing owls and falcons elsewhere."
- [24] Genes 2022 retroelementos (`41a1c3e3`) — legenda da Fig. 2 (hipóteses publicadas para Telluraves; Afroaves não monofilético em várias) e discussão sobre erro de reconstrução de árvores de genes ("55–79% of internal nodes in gene trees are reconstructed inaccurately").
- [25] Stiller et al. 2024 (`e9c5ae30`) — "Within Telluraves our main tree supported the proposed split into Australaves and Afroaves…" (mesmo trecho de Q1 [4]).

---

### Q8 — Tempos de divergência: Accipitridae coroa, Falconidae, Strigiformes e espinha dorsal de Telluraves

**Pergunta:** What divergence times are estimated for crown Accipitridae, Falconidae, Strigiformes and the Telluraves backbone?

conversation_id: `0f9bbec7-a5d4-486e-b327-42231f3f3c92`

**Resposta do NotebookLM (verbatim):**

> Estimates for the divergence times of **Telluraves**, **Accipitridae**, **Falconidae**, and **Strigiformes** vary across different molecular datasets, fossil calibration strategies, and dating methodologies in your sources:
>
> ---
>
> ### 1. Telluraves Backbone (Core Landbirds)
> The timing of the initial divergence and backbone radiation of Telluraves is subject to two main evolutionary models:
>
> * **Post–K–Pg Radiation (~62–66 Ma):** Phylogenomic analyses (such as those based on whole-genome or UCE matrices) support a rapid post-extinction radiation following the Cretaceous–Paleogene (K–Pg) mass extinction (~66 Ma) [1-3]. The discovery of the Paleocene mousebird fossil *Tsidiiyazhi abini* (**62.221–62.517 Ma**) establishes a hard minimum fossil age of ~62.2–62.5 Ma for the lineages leading to Australaves, Accipitriformes, Strigiformes, and Cavitaves, narrowing the post-K–Pg window for major landbird backbone divergences to less than 4 million years [4].
> * **Late Cretaceous / Pre–K–Pg Radiation (~80–86 Ma):**
>   * Total-evidence tip-dating analyses incorporating fossil characters recover a **Late Cretaceous origin** for crown Telluraves, with a few lineages crossing the K–Pg boundary prior to rapid post-extinction diversification [5].
>   * Coalescence-based analyses by Wu et al. (2024) place the basal divergence of Telluraves at **~86 Ma** (95% CI: 80.4–91.5 Ma), supporting a continuous radiation across the K–Pg boundary [6, 7].
>   * 3′-UTR analyses under relaxed-clock models also yield basal landbird divergence dates preceding the K–Pg boundary (~66–80+ Ma) [8].
>
> ---
>
> ### 2. Crown Accipitridae (Hawks, Eagles, Harriers, and Old World Vultures)
> While the order Accipitriformes split from Cathartiformes (~61.3 Ma) and Pandionidae (~50.8 Ma) early in the Cenozoic [9], modern crown Accipitridae diverged later in the Cenozoic [9, 10]:
>
> * **Basal Split (Elaninae vs. Remaining Accipitrids):**
>   * Ultraconserved element (UCE) time-calibrated studies estimate that the basal subfamily **Elaninae** split from the rest of Accipitridae at **45.3 Mya** (Middle Eocene) [9].
>   * Other nuclear/mitochondrial clock studies estimate this basal split at **~33.7–34 Ma** near the Eocene–Oligocene transition [10, 11].
> * **Subfamily Diversification (25.1–33 Mya / ~10–24 Mya):**
>   * UCE-based estimates place the divergence of the remaining major subfamily clades (Perninae, Gypaetinae, Circaetinae, Aegypiinae, Aquilinae, Harpiinae) between **25.1 and 33 Mya** [9].
>   * Earlier multi-locus studies estimate that the Perninae + Gypaetinae clade split from other accipitrids at **~23–24 Ma**, followed by the emergence of more derived subfamilies (Aquilinae, Buteoninae, Accipitrinae, and Aegypiinae) during the Miocene (**10–17 Ma**) [11-14]. Global fossil evidence confirms that all extant subfamilies had evolved by **~20 Ma** [15].
> * **Genus & Species Radiations:** Most genus-level splits within subfamilies occurred during the Late Miocene to Pliocene (≥5 Ma), with species-level radiations continuing into the Plio-Pleistocene [9, 13].
>
> ---
>
> ### 3. Falconidae (Falcons, Caracaras, and Forest-Falcons)
> The Falconidae are phylogenetically distinct from Accipitridae and belong to Australaves [16]. Their internal crown diversification is calibrated as follows:
>
> * **Crown Falconidae Origin:** Earliest divergences among the three subfamilies (Herpetotherinae, Polyborinae, and Falconinae) began in the **Late Oligocene to Early Miocene**, with a maximum constraint for crown Falconidae estimated at **34.2 Mya** (Earliest Oligocene) [17-19].
> * **Subfamily Falconinae:** Internal divergences within Falconinae began about **16 Mya** (Middle Miocene) [18].
> * **Genus *Falco* Radiation:** Divergences within the genus *Falco* (comprising ~60% of falconid species) began **~7.5 Mya** (Late Miocene), coinciding with the expansion of open C4 grasslands and savannas near the Miocene–Pliocene transition [18].
>
> ---
>
> ### 4. Strigiformes (Owls)
> * **Stem Lineage Origin:** The Strigiformes lineage spans at least **60–65 million years** [20, 21]. The oldest confirmed fossil strigiform, *Ogygoptynx wetmorei*, dates to the Paleocene (**65–56.5 Ma**) [21], and *Tsidiiyazhi abini* establishes an early Paleocene hard minimum age constraint of **~62.2–62.5 Ma** for the owl stem lineage [4].
> * **Crown Strigidae Radiation:** While Paleogene fossil owl families (such as Sophiornithidae and Protostrigidae) were diverse, the extant family **Strigidae** (typical owls) is absent from the Paleogene fossil record and underwent its primary radiation later in the **Miocene (23.5–5.2 Ma)** [22, 23].
>
> Would you like to explore the phylogenetic placement or fossil evidence for any of these raptorial clades in more detail?

**Citações (número → fonte; trecho citado, truncado):**

- [1] Stiller et al. 2024 (`e9c5ae30`) — "…calibration densities for 34 nodes using 187 fossil occurrences… Bayesian sequential-subtree framework… branch lengths from intergenic regions… In contrast to a recent study proposing a diversification of Neoaves during the Upper Cretaceous, we found that the early divergences in Neoaves were tightly associated with the K–Pg boundary. Only two divergences occurred before the boundary: Mirandornithes diverged from the remaining Neoaves 67.4 Ma (95% CrI 66.2–68.9) and Columbaves diverged 66.5 Ma (65.2–67.9). All subsequent divergences postdate the boundary, although the 95% credible interval of the divergence time between Telluraves and Elementaves and the crown age of Elementaves spans the K–Pg boundary…"
- [2] Prum et al. 2015 (`74f2801f`) — "…The results of our divergence time analyses are congruent with the palaeontological record, supporting a major radiation of crown birds in the wake of the Cretaceous–Palaeogene (K–Pg) mass extinction."
- [3] Jarvis et al. 2014 (`c92cd31c`) — "…the Neognathae diverged into… Galloanseres and Neoaves about 88 Ma… an explosive diversification within Neoaves were dated to occur around the K-Pg boundary… the most basal divergences within Neoaves (Columbea, Passerea, and two more) occurred before the K-Pg transition (67 to 69 Ma) and all others after, with nearly all ordinal divergences completed by 50 Ma… The estimated age for the basal split of Passeriformes… was around 39 Ma…"
- [4] Ksepka et al. 2017 (`b84ec51f`) — "Under our preferred tree, Tsidiiyazhi pushes the hard minimum fossil ages of the lineages leading to four major landbird clades into the early Paleocene: Australaves, Accipitriformes, Strigiformes, and Cavitaves… Tsidiiyazhi (62.221–62.517 Ma) is older than the stem penguin Waimanu (60.5–61.6 Ma)… This narrows the post-Cretaceous mass extinction window for the divergences among most major neoavian clades to less than 4 million years…"
- [5] Tip-dating and the Origin of Telluraves (ResearchGate, resumo, `79ba75e1`) — "…Our analyses recover a Late Cretaceous origin for crown Telluraves, with a few lineages crossing the K-Pg boundary…"
- [6] Wu et al. 2024 (`cda64d61`) — "We performed dating analysis on the optimum tree using MCMCTree with CDS data and the full set of 20 fossil constraints… interordinal diversification of modern birds largely occurred in the Late Cretaceous and underwent a continuous radiation across the KPg boundary… approximately 15 avian orders that likely originated in the Paleogene, whereas the origin of the remaining 24 orders falls within the Late Cretaceous…"
- [7] Wu et al. 2024 (`cda64d61`) — "…the basal divergences of both Neoaves (~87.2 Ma, 95% CI: 81.6 to 92.8) and placental mammals (~87.3 Ma…) occurred almost concurrently… the basal divergence of Telluraves (landbirds, ~86 Ma, 95% CI: 80.4 to 91.5)… and the basal divergence of Aquaterraves (~81 Ma, 95% CI: 75.6 to 86.4)…"
- [8] Kuhl et al. 2021 (`8527a287`) — "…divergence dates for Palaeognathae and Neognathae (94 Ma) and Galloanserae and Neoaves (81 Ma) that are much earlier than those suggested previously (Prum et al. 2015). The divergence dates of the Mirandornithes and those of the Basal Landbirds and of the Aquatic & Semiaquatic Bird lineages precede the K-Pg boundary. The initial divergences within many other neoavian lineages occurred ∼10 My after the K-Pg boundary, in the Eocene… During the Oligocene epoch, a second major diversification event occurred… between 35 and 25 Ma…"
- [9] Catanach et al. 2024 (`1d4669b1`) — datas (mesmo trecho de Q3 [6]).
- [10] Mather 2021 (`c24ce4bf`) — "…the most basal living genera, Gampsonyx and Elanus, diverged from the rest of the Accipitridae roughly 34 Ma, and the split between the Perninae+Gypaetinae clade and the rest… roughly 23 Ma, in agreement with Nagy and Tökölyi (2014)…"
- [11] Mather 2021 (`c24ce4bf`) — "…Dated molecular phylogenies imply that most of the extant accipitrid subfamilies had not diverged by the late Oligocene, with only the Elaninae, which diverged at 33.7 Ma (Mindell et al. …"
- [12] Mather 2021 (`c24ce4bf`) — "Nagy and Tökölyi (2014)… estimate likely divergence dates of the accipitrid subfamilies… most of the more basal subfamilies (Elaninae, Perninae, Gypaetinae, Harpiinae, Circaetinae) originated between 30–13 Ma,"
- [13] Mather 2021 (`c24ce4bf`) — "while more derived subfamilies (Aquilinae, Buteoninae and Haliaeetinae) tended to date between 17–10 Ma. Interestingly, the Aegypiinae dated between 13–7 Ma…"
- [14] Mather 2021 (`c24ce4bf`) — "…the Perninae + Gypaetinae clade diverged at 23.8 Ma from remaining accipitrids (Mindell et al. 2018). Other lineages emerged during or after the middle Miocene (Nagy and Tökölyi 2014; Oatley 2015; Prum 2015; Mindell 2018)…"
- [15] Mather 2021 (`c24ce4bf`) — "The global fossil record, and molecular divergence dating, suggests that by the early Miocene, about 20 Ma, all extant subfamilies of the Accipitridae had evolved (Walker and Dyke 2006; Zhang et al. 2012; Nagy and Tökölyi 2014; Mindell et al. 2018)."
- [16] Mayr & Kitchener 2022 (`53dae551`) — introdução (mesmo trecho de Q2 [23]).
- [17] Mayr & Kitchener 2022 (`53dae551`) — "…the last common ancestor of crown group Falconiformes is likely to have lived in the New World, because the South American Herpetotherinae and Polyborinae are successive sister taxa of other crown group Falconiformes (Fuchs et al. 2015). Molecular data yielded a 34.2-million-year maximum time estimate for the earliest divergences within crown group Falconiformes, in the earliest Oligocene (Fuchs et al. 2015). This date is much younger than the latest (middle Eocene) occurrence of the Masillaraptoridae. The Antarctic taxon Antarctoboenus is likewise at least 15 million years older than the presumed split of crown group Falconiformes…"
- [18] Fuchs et al. 2015 (`f82bdd8e`) — resumo (mesmo trecho de Q5 [1]).
- [19] Fuchs et al. 2015 (`f82bdd8e`) — introdução (mesmo trecho de Q5 [2]).
- [20] Dalsätt, tese (`6899a064`) — "…One such group is the Strigiformes (owls) with the oldest confirmed fossil from the Paleocene. We describe a new species from the Eocene Green River Formation… closely related to… Prosybris antique and P. medius…"
- [21] Dalsätt (`6899a064`) — "…the Strigiformes has a long evolutionary history, spanning at least 60 mya… The oldest confirmed fossil strigiform is the Paleocene Ogygoptynx wetmorei (65 – 56.5 Ma)…"
- [22] Dalsätt (`6899a064`) — "…Sophiornithidae… also recorded from the Eocene (56.5 – 35.5 Ma). Other families… Protostrigidae…, Palaeglaucidae…, and the extant Tytonidae…"
- [23] Dalsätt (`6899a064`) — "…Strigidae, is not known from the Paleogene and began to radiate extensively only in the Miocene (23.5 – 5.2 Ma)…"

---

### Q9 — Nós ainda não resolvidos e a quem os autores atribuem a culpa (ILS, tipo de dado, modelo)

**Pergunta:** Which raptor-relevant nodes remain unresolved, and what do authors blame (ILS, data type, model)?

conversation_id: `d2c7c53f-3d4f-4018-91fa-51cfe8761b3d`

**Resposta do NotebookLM (verbatim):**

> ### Unresolved Raptor-Relevant Nodes
>
> Phylogenomic studies consistently demonstrate that raptors are **polyphyletic** across **Telluraves** (core landbirds), with raptorial lifestyles arising in distinct, distantly related lineages: **Accipitriformes** (hawks, eagles, New World vultures, secretarybird), **Falconiformes** (falcons and caracaras), **Strigiformes** (owls), and **Cariamiformes** (seriemas) [1-4]. Several key nodes involving these groups remain deeply contentious across modern studies:
>
> 1. **Placement of Accipitriformes relative to Telluraves**:
>    * Concatenated and species-tree analyses of total evidence or non-coding data often place Accipitriformes as the sister group to all other **Afroaves** [5-7].
>    * Conversely, whole-genome, exon, or reweighted coalescent analyses (e.g., Prum et al. 2015, wASTRAL) frequently recover Accipitriformes as the **sister to all other core landbirds (Telluraves)** [8-10].
>
> 2. **Placement of Strigiformes (Owls) and "Hieraves"**:
>    * Owls shift dramatically depending on the dataset and method:
>      * Placed as sister to **Cavitaves** (cuckoo-rollers, trogons, woodpeckers, etc.) in retroposon (RI) and 3'-UTR trees [6, 11, 12].
>      * Placed as sister to **Accipitriformes**—forming a clade termed **Hieraves**—in recent whole-genome and morphological matrix analyses [7, 13, 14].
>      * Placed sister to a broader clade including Coliiformes or Coraciimorphae in alternative nuclear sequence analyses [6, 15].
>
> 3. **Monophyly and Internal Structure of Afroaves**:
>    * Whether **Afroaves** forms a valid monophyletic group or is rendered paraphyletic by the basal branching of Accipitriformes remains unresolved [2, 6, 9, 15, 16].
>    * The relative branching order among four ancient afroavian lineages—Accipitriformes, Strigiformes, Coliiformes (mousebirds), and Cavitaves—differs across nearly every major published tree [6, 15].
>
> 4. **Position of Coliiformes (Mousebirds)**:
>    * Mousebirds frequently shift between being sister to Cavitaves or sister to all other Afroaves, directly altering the inferred relative positions of hawks and owls [6, 15, 17-19].
>
> 5. **Genus-Level Non-Monophyly (e.g., *Accipiter*)**:
>    * Dense ultraconserved element (UCE) and legacy sequencing confirms that the traditional genus ***Accipiter*** (accipitrid hawks) is non-monophyletic, as harriers (*Circus*) and the genera *Erythrotriorchis* and *Megatriorchis* are embedded within it [20-24].
>
> ---
>
> ### Causes Identified by Authors
>
> Authors attribute these recalcitrant nodes to a combination of biological processes, data biases, and analytical limitations:
>
> #### 1. Incomplete Lineage Sorting (ILS) & the "Anomaly Zone"
> * **Rapid Adaptive Radiation**: Early Neoaves and Telluraves underwent explosive speciation near the Cretaceous–Paleogene (K-Pg) boundary, resulting in extremely short internal branch lengths relative to large ancestral effective population sizes (\\(N_e\\)) [25-29].
> * **High Hemiplasy**: Up to **40%–100% of loci** across deep neoavian internodes are affected by ILS, where ancestral polymorphisms persist across multiple successive speciation events before sorting stochastically [26, 28-30].
> * **The Anomaly Zone**: In extreme cases, short internodes and long terminal branches create an "anomaly zone," where the most common gene tree topology differs from the true species tree [31-33].
>
> #### 2. Data Type Effects
> * **Signal Conflict Across Genomic Regions**: Authors report stark topological incongruences when analyzing **exons vs. introns vs. UCEs vs. 3'-UTRs vs. Rare Genomic Changes (RGCs)** like retroposon insertions (RIs) or indels [16, 34-37].
> * **Coding Sequence Biases**: Protein-coding exons suffer from **GC-content shifts**, compositional heterogeneity, rate variation, functional/morphological convergence, and alignment artifacts [37-41].
> * **Non-Coding & Recoded Signals**: Non-coding loci (introns, UCEs, intergenic regions) or RY-coded alignments often yield distinct topologies from unrecoded exon sequences [34, 36, 37, 39, 40].
>
> #### 3. Model Misspecification & Analytical Method Conflicts
> * **Concatenation vs. Multispecies Coalescent (MSC)**: Concatenation (supermatrix Maximum Likelihood) ignores ILS and can produce artificially high bootstrap support for incorrect topologies [16, 42, 43].
> * **Gene Tree Reconstruction Error**: MSC summary methods (e.g., ASTRAL, MP-EST, SVD-quartets) are highly sensitive to gene tree estimation error when individual gene alignments are short or weakly informative [16, 42, 43].
>
> #### 4. Ancient Introgressive Hybridization (Gene Flow)
> * **Non-ILS Asymmetries**: Statistical tests on retroposons (e.g., KKSC, 4-LIN) and site patterns reveal asymmetrical quartet frequencies at the node uniting owls, hawks, and Cavitaves/Afroaves [7, 44-46]. This departure from standard coalescent expectations suggests **ancient gene flow or reticulate evolution** during early afroavian diversification [7, 15, 44-47].
>
> #### 5. Taxon Sampling Effects
> * Stiller et al. (2024) discovered that the relationship between Accipitriformes, Strigiformes, and Afroaves depended heavily on the number of passerines sampled: the main "Hieraves + Afroaves" topology required sampling at least 138 passeriform taxa, whereas sparse sampling shifted Accipitriformes to sister of remaining Afroaves [7, 48].
>
> ---
>
> Would you like to explore a detailed comparison of specific raptorial clades, or look into generating a visual summary or mind map of these phylogenomic conflicts?

**Citações (número → fonte; trecho citado, truncado):**

- [1] Birds, Diversification of (ResearchGate, capítulo de enciclopédia) (`d3d5a7ff`) — "a remarkable example of convergent evolution among Phaethoquornithes ( Mayr et al., 2021). (8) Telluraves: “ Core landbirds”, the most speciose of the eight major neoavian groups and encompassing an ecologically disparate array of mainly terrestrial and arboreal birds. Among these are several primarily raptorial (vertebrate-eating) lineages, namely Accipitriformes (hawks, eagles, vultures and their close relatives), Strigiformes (owls), Cariamifo…"
- [2] Resolving the Avian Tree of Life from Top to Bottom (ResearchGate) (`47d8810c`) — "(Passeriformes, typically called passerines), parrots, and the woodpeckers and their allies (such as rollers, king ﬁ shers, bee-eaters, hornbills, hoopoes and wood hoopoes, trogons, and other lineages). The Jarvis et al. ( 2014) TENT, like a number of prior analyses (e.g., Ericson et al. 2006; Hackett et al. 2008; Kimball et al. 2013), splits core landbirds into two clades that Ericson ( 2012) named Australaves and Afroaves (Table 3). That split …"
- [3] Wang et al. 2014, mitogenomas de rapinantes (GBE) (`3cc767df`) — "Phylogenetic Position of Avian Nocturnal and Diurnal Raptors PMCID: PMC3942033 DOI: 10.1093/gbe/evu016 Fonte: Europe PMC fullTextXML Abstract We report three new avian mitochondrial genomes, two from widely separated groups of owls and a falcon relative (the Secretarybird). We then report additional progress in resolving Neoavian relationships in that the two groups of owls do come together (it is not just long-branch attraction), and the Secreta…"
- [4] Catanach et al. 2024, Accipitridae UCE (Biol J Linn Soc) (`1d4669b1`) — "Before the mid-20th century, most systematists assumed that falcons (Falconidae) and hawks and eagles (Accipitridae) were sister groups and thus classified both in the order Falconiformes (e.g. Mayr 1959, but see Starck 1959). However, based on electrophoretic profiles of egg-white proteins, Sibley (1960) pro- posed that ‘Falconiformes [sensu lato] may be polyphyletic, the Falconidae possibly being unrelated to the other diurnal birds of prey’. T…"
- [5] Resolving the Avian Tree of Life from Top to Bottom (ResearchGate) (`47d8810c`) — "the falcons are placed in Australaves sister to parrots and passerines whereas the hawks, eagles, New World vultures, and owls are placed in Afroaves as the succes- sive sister groups of a clade comprising mousebirds, and a diverse assemblage that Resolving the Avian Tree of Life from Top to Bottom: The Promise and... 173 includes woodpeckers and their allies (Fig. 2a). However, the Jarvis et al. ( 2014) TENT con ﬂ icts with the Prum et al. ( 201…"
- [6] Genes 2022, retroelementos e introgressão em Afroaves (`41a1c3e3`) — "In addition to nearly universal incongruence with previous phylogenomic work at the earliest branching events in Neoaves (Figure 4), the quartet-based species tree for RI characters also conflicts with most prior genome-scale hypotheses for Telluraves, an ecologically diverse clade of core landbirds. Among DNA-sequence based trees, just a single ML concatenation analysis of UCEs [1] was wholly congruent for the major lineages of Telluraves shown …"
- [7] Stiller et al. 2024 (Nature, B10K) (`e9c5ae30`) — "Difficult placement of owls and hawks Within Telluraves our main tree supported the proposed split into Australaves and Afroaves1,3 in contrast to other studies2,11. Our tree grouped Accipitriformes and Strigiformes as the sister to the remaining Afroaves, similar to previous coalescent-based analyses1. Concatenated analyses1,3, including ours, supported Accipitriformes alone as sister to the remaining Afroaves (Extended Data Fig. 4). This node a…"
- [8] Resolving the Avian Tree of Life from Top to Bottom (ResearchGate) (`47d8810c`) — "TENT (and Reddy et al. 2017) places accipitriforms sister to all other Afroaves whereas Prum et al. ( 2015) places them sister to all other core landbi rds (Fig. 2b). The Jarvis et al. ( 2014) analysis of ﬁ rst and second codon positions suppor ted a third topology, with a clade comprising accipitriforms and owls siste r to all other core landbirds. Although there is some con ﬂ ict between the Jarvis et al. ( 2014) exon analysis and the Prum et a…"
- [9] Resolving the Avian Tree of Life from Top to Bottom (ResearchGate) (`47d8810c`) — "... Patterns for branch 27, which unites the putative clade Afroaves (Fig. 4), were more complex because it is so unstable that a wASTRAL reanalysis of the Stiller et al. (2024) gene trees fails to recover it. Instead, it puts Accipitriformes sister to all landbirds, a pattern observed before (Braun et al. 2019 , Stiller et al. 2024. wASTRAL analysis of alternative state gene trees moved Accipitriformes to a clade including Strigiformes and sever…"
- [10] Braun & Kimball 2021, Data Types and the Phylogeny of Neoaves (`b0e43963`) — "The other part of the bird tree where these analyses provide insights is the base of the core landbirds. The hypothesis that most raptorial landbirds form a paraphyletic assemblage at the base of core landbirds (e.g., the 1st and 2nd codon position exon tree from Jarvis et al. [ 11 ] and Figure 3 B from this study) can probably be rejected. Instead, the majority of recent analyses [ 6 , 8 , 11 , 13 , 16 ] favor a topology that divides raptors (an…"
- [11] Chen et al. 2025, matriz anatômica cintura peitoral (IOB) (`6c89ed1d`) — "Strict consensus of most parsimonious tree resulting from Analysis 10 (equal weights parsimony analysis run with the topology constrained to the consensus tree proposed by Braun and Kimball 2021). Numbers on nodes represent absolute bootstrap frequencies (left of bar) and Bremer support values (right of bar). A hyphen (-) denotes values of 0. Nodes lacking support values were constrained to be monophyletic. Analysis 11 recovered one MPT of 2755 s…"
- [12] Genes 2022, retroelementos e introgressão em Afroaves (`41a1c3e3`) — "Conflicting RI evidence for deep relationships within Afroaves. (A) The ASTRAL_BP coalescent analysis of RIs supports Strigiformes (owl) as sister to Cavitaves; local branch swaps relative to the optimal hypothesis are also shown with Bayesian local PPs above internodes and perfectly congruent RI synapomorphies below internodes. (B) Positive partitioned coalescence support (PCS) scores for the optimal tree (green) and negative PCS scores for the …"
- [13] Chen et al. 2025, matriz anatômica cintura peitoral (IOB) (`6c89ed1d`) — "Stiller et al. (2024) and Wu et al. (2024b) instead recovered Accipitrimorphae and Strigiformes as sister groups, forming a clade called Hieraves. Six character states were optimized as potential synapomorphies of this clade, though only one was inferred as such under both topologies. Keel apex on sternum caudal to main body of sternum (char. 43: 1 > 0, also found in Eudromia, Galliformes, Monias, Opisthocomus, Aramus, Sarothrura, Rallus, Charadr…"
- [14] Wu et al. 2024 (PNAS) (`cda64d61`) — "Discussion Phylogenetic Relationships of Neoaves. Our optimum tree provides consistent resolution of phylogenetic relationships within Neoaves, which comprise the bulk of modern birds and whose relationships have long been phylogenetically controversial, dividing them into two basal clades, a previously recognized Telluraves (landbirds) (33) and a newly circumscribed Aquaterraves (waterbirds and relatives) (Fig. 2). Within Telluraves, the optimum…"
- [15] Genes 2022, retroelementos e introgressão em Afroaves (`41a1c3e3`) — "Informal consensus species tree for Neoaves based on three recent summaries [3,4,8]. Terminal taxa indicate genomes analyzed by Jarvis et al. [1], and higher-level taxa are shown at internodes and to the right of terminals (delimited by color). Three clades, Afroaves, Musophagotides, and Coraciimorphae + Strigiformes (owls) are not resolved in the consensus tree but were featured in the summaries of Reddy et al. [3], Houde et al. [8], and Braun e…"
- [16] Genes 2022, retroelementos e introgressão em Afroaves (`41a1c3e3`) — "Previous phylogenetic hypotheses for Telluraves based on large molecular datasets. The RI polymorphism parsimony tree [5] is congruent with just one phylogenomic analysis, ML concatenation of UCEs (A), but conflicts with other published studies (B–I), many of which place Coliiformes (mousebirds) sister to Cavitaves (B–F). Clades that conflict with the RI parsimony tree are marked by red circles at nodes. Afroaves (green background) is not monophy…"
- [17] Resolving the Avian Tree of Life from Top to Bottom (ResearchGate) (`47d8810c`) — "concatenated analyses of exon data, either placing an accipitriform-owl clade or accipitriforms alone sister to all other landbirds (Fig. 5). Thus, the position of accipitriforms and owls represents an interesting case of uncertainty relat ed to two different factors: (1) data type and (2) whether the analytical method assumes a single underlying tree or a mixture of trees (for additional details regard ing analytical methods, see below in the ne…"
- [18] Genes 2022, retroelementos e introgressão em Afroaves (`41a1c3e3`) — "4.5. Weighted ASTRAL Analyses and the Phylogenetic Placement of Mousebird ‘Total evidence’ coalescent analyses of >2000 RIs plus >6000 gene trees recorded striking changes in both the topology and clade support within Telluraves. More severe downweighting of poorly supported gene tree nodes (9B–D) yielded higher congruence with the ASTRAL_BP species tree for RIs (Figure 9A). In particular, the high local PPs for mousebird + Cavitaves (Figure 8A a…"
- [19] Suh et al. 2015 (PLoS Biol) (`2d01758b`) — "The Mousebird Conflict Is Not Caused by Hybridization It is striking that the conflicting placements of mousebirds are well-supported in the main Jarvis et al. tree [4] on the one side and our MPRE tree and the genome-scale UCE tree [4] on the other side, respectively. One explanation for this could be hybridization of two distinct and diverged ancestral species, e.g., the ancestor of mousebirds and the ancestor of the woodpecker/bee-eater/hornbi…"
- [20] AOS/NACC proposta 2025-C (`b5e56e8d`) — "“Accipiter” plus Megatriorchis and Circus. These 6 species include Eurasian Sparrowhawk (Accipiter nisus), which is the type of the genus Accipiter, so these 6 species would retain the name Accipiter. In addition to nisus, it also includes madagascariensis, ovampensis, rufiventris, striatus, and poliogaster (Catanach et al. 2024). The inclusion of poliogaster in this clade is provisional; although it received some support in Catanach et al. (2024…"
- [21] AOS/NACC proposta 2025-C (`b5e56e8d`) — "finding poliogaster instead sister to Megatriorchis, although this was not well-supported. 4. The final clade of “Accipiter” hawks identified by Catanach et al. (2024) was sister to Megatriorchis and Circus. The oldest available genus for this clade is Astur Lacépède 1799, with Astur gentilis designated as the type for the genus. This clade includes 9 species: cooperii, gundlachi, bicolor, chilensis, melanoleucus, henstii, gentilis, atricapillus,…"
- [22] Catanach et al. 2024, Accipitridae UCE (Biol J Linn Soc) (`1d4669b1`) — "A B ST R A CT Hawks, eagles, and their relatives (Accipitriformes: Accipitridae) are a diverse and charismatic clade of modern birds, with many members that are instantly recognized by the general public. However, surprisingly little is known about the relationships among genera within Accipitridae, and several studies have suggested that some genera (in particular, the megadiverse genus Accipiter) are not monophyletic. Here, we combine a large n…"
- [23] Catanach et al. 2024, Accipitridae UCE (Biol J Linn Soc) (`1d4669b1`) — "D ow nloaded from https://academ ic.oup.com /biolinnean/advance-article/doi/10.1093/biolinnean/blae028/7633854 by M edical Sciences Library of Texas A&M U niversity user on 22 M arch 2024 Accipiter and Circus, which convinced Olson (2006), but this finding was soon followed by a wave of evidence that Accipiter, as traditionally defined (Peters 1931), does not form a mono- phyletic group (e.g. Griffiths et al. 2007, Hugall and Stuart-Fox 2012, Oat…"
- [24] Catanach et al. 2024, Accipitridae UCE (Biol J Linn Soc) (`1d4669b1`) — "R E SU LTS Phylogenetic hypotheses The UCE backbone phylogeny of 134 species (Fig. 2; for detailed divergence estimates, see Supporting Information, Fig. S2) was fully supported at virtually every node (bootstrap values = 100). Only five nodes had less than perfect support, and four of those represented splits between closely related species or groups of species (two within Buteo, one within Geranoaetus, and one within Gyps). The fifth node with …"
- [25] Resolving the Avian Tree of Life from Top to Bottom (ResearchGate) (`47d8810c`) — "from the non-DC trees. Mixtures with multiple topologies arise when DC gene trees are discordant with the species tree (e.g., the dashed gene tree in Fig. 6e). There is strong direct evidence for heterotachy (e.g., note the very different terminal branch lengths for BDNF and PPP2CB in Fig. 8), although the impact of heterotachy on estimates of avian phylogeny using standard ML methods remains unclear. The evidence for discordance among gene trees…"
- [26] Jarvis et al. 2014 (Science) (`c92cd31c`) — "Indels suggest a high degree of ILS at the earliest branches of the Neoaves tree We further assessed ILS using insertions and deletions (indels) (69), because they have less homoplasy (convergence) than single nucleotides (SM9), and unlike gene trees, indels can be examined as discrete characters mapped to a reference tree without the added inference of constructing trees from them. We scored 5.7 million indels from the TENT alignment, of which 2…"
- [27] Suh et al. 2015 (PLoS Biol) (`2d01758b`) — "Introduction The rich biodiversity of many organismal groups is the result of bursts of rapid species diversifications, with extreme examples in angiosperms [1] and vertebrates [2]. Among the latter, birds are one of the most speciose groups with a total of >10,500 recognized species that are proposed to be the result of mostly recent accelerations of diversification rates [3]. Nevertheless, the deep roots of 95% of these species lie within the a…"
- [28] Suh et al. 2015 (PLoS Biol) (`2d01758b`) — "Temporal Extent of Incomplete Lineage Sorting Mapping RE markers on a dated time tree of the main Jarvis et al. analysis [4] enabled us to estimate the temporal dynamics of ILS across the very short internodes of the neoavian radiation (Fig 3A). Jarvis et al. infer the onset of Neoaves diversification at around the K-Pg boundary [4], which is in stark contrast to most mitochondrial and multilocus nuclear studies (but see refs. [17,26]) that estim…"
- [29] Suh et al. 2015 (PLoS Biol) (`2d01758b`) — "Phylogenetic Extent of Incomplete Lineage Sorting Complex gene tree–species tree incongruences of retrotransposon markers might be more accurately represented in phylogenetic networks where data conflicts are evident as reticulate relationships among taxa and alternative topologies are visible even within well-supported lineages [12,29]. The resultant neighbor-net [31] (Fig 4A) illustrates the differential distribution of ILS-derived incongruence…"
- [30] Suh et al. 2015 (PLoS Biol) (`2d01758b`) — "The remaining 745 retrotransposon markers show different degrees of gene tree–species tree incongruence. This is best explained by the persistence of ancestral polymorphisms across successive speciation events, followed by reciprocal allele fixation in each of the descendant lineages, i.e., ILS. We define the extent of ILS as corresponding to weak conflict (persistence across two speciation events; Fig 1C), moderate conflict (three events; Fig 1D…"
- [31] Jarvis et al. 2014 (Science) (`c92cd31c`) — "Discussion Our study is an example of the extraordinary amount of genomic sequence data required to produce a highly supported phylogeny spanning a rapid radiation. The conflict we observe with other data types (14,15, 24) can no longer be considered to be due to error from smaller amounts of sequence data (8,17) nor to differences in concatenation versus coalescence methods (27, 28). The absence of a single gene tree identical to the avian speci…"
- [32] Suh et al. 2015 (PLoS Biol) (`2d01758b`) — "Genomic Extent of Incomplete Lineage Sorting Retrotransposon loci affected by ILS are distributed across the avian genome, irrespective of the duration of ILS per intronic or intergenic marker (Fig 5A). This situation likely applies to the ancestral Neoaves genome, as the avian karyotype is unusually stable, including conserved synteny of the Z sex chromosome and ubiquitous presence of numerous microchromosomes [32,33]. The Z chromosome (no recom…"
- [33] Suh et al. 2015 (PLoS Biol) (`2d01758b`) — "Conclusions We conclude that Neoaves diversification is more complex than can be shown in fully bifurcating trees and exhibits a dynamic picture of ILS. The timing of the highly ILS-affected initial super-radiation coincides with the K-Pg extinction of nonavian dinosaurs and archaic birds [42], suggesting that the abrupt availability of ecological niches [4] was followed by near-simultaneous population isolations [41] via specializations and led …"
- [34] Resolving the Avian Tree of Life from Top to Bottom (ResearchGate) (`47d8810c`) — "remain poorly resolved [Fig. 1; also see Thomas ( 2015) for another compa rison of the Jarvis et al. ( 2014) and Prum et al. ( 2015) trees]. Moreover, differences across Resolving the Avian Tree of Life from Top to Bottom: The Promise and... 165 recent studies are sensitive to the data type used for analyses (i.e., the estimates of phylogeny differ depending on whether exons, introns, noncoding ultraconserved elements, conserved non-exonic elemen…"
- [35] Resolving the Avian Tree of Life from Top to Bottom (ResearchGate) (`47d8810c`) — "2014). In contrast, Prum et al. ( 2015) places Caprimulgiformes sister to all other Neoaves (Fig. 2). Both Jarvis et al. ( 2014) and Prum et al. ( 2015) dispersed the other metavian taxa across several other basal neoavian lineages. Reddy et al. ( 2017) analyzed the topological con ﬂ icts between Jarvis et al. ( 2014) and Prum et al. ( 2015) and concluded that the observed con ﬂ icts re ﬂ ect “ data-type effects,” which they de ﬁ ned as the obser…"
- [36] Resolving the Avian Tree of Life from Top to Bottom (ResearchGate) (`47d8810c`) — "TENT) place the tropicbird-eurypygiform clade as sister to waterbirds, whereas analyses of exclusively or almost exclusively noncoding data matrices [the analyses of introns and noncoding UCEs in Jarvis et al. ( 2014) and Reddy et al. ( 20 17)] place the tropicbird-eurypygiform clade sister to landbirds. However, if this is a data-type effect, it is unusual. Reddy et al. ( 2017) proposed the terminology data-type effects in order to discuss the c…"
- [37] Stiller et al. 2024 (Nature, B10K) (`e9c5ae30`) — "Conclusions Relationships along the backbone of Neoaves have long been contentious, with various analyses yielding incongruent results. At the heart of the disagreements has been a long-standing question: is it better to sample many taxa at a few loci (typically conserved regions, such as exons and UCEs) or sample many loci widely across the genome, even if available from fewer species? We can finally answer this question because our data provide…"
- [38] Resolving the Avian Tree of Life from Top to Bottom (ResearchGate) (`47d8810c`) — "“birds of prey” (owls, hawks, eagles, & New World vultures) “birds of prey” (falcons & seriemas) Fig. 7 Phylogram of the Prum et al. ( 2015) ML tree emphasizing shifts in GC-content and evolutionary rate. Taxa with the most extreme values for the GC-content of the parsimony informative sites and the “ fast sites”(sites in the 95th percentile for number of MP steps given this tree) are indicated using colored arrows (red for high values and blue f…"
- [39] Braun & Kimball 2021, Data Types and the Phylogeny of Neoaves (`b0e43963`) — "3.3. RY Coding also Affects Species Tree Analyses Prum et al. [ 12 ] also estimated avian phylogeny using several multispecies coalescent (“species tree”) methods, such as ASTRAL. Given the increased congruence between the non-coding trees and the trees estimated with RY coding, we compared the performance of ASTRAL with and without RY-coding before estimating gene trees. Both ASTRAL trees had limited support, but RY coding had a substantial impa…"
- [40] Genes 2022, retroelementos e introgressão em Afroaves (`41a1c3e3`) — "1. Introduction The early radiation of Neoaves, the clade that includes all extant birds except Palaeognathae (ratites, tinamous) and Galloanserae (land- and water-fowl), has defied convincing systematic resolution. The early diversification of the group near the Cretaceous–Paleogene boundary exemplifies the challenges that confront genome-scale analysis when the splitting of evolutionary lineages is rapid and deep in time [1,2]. Phylogenomic app…"
- [41] Genes 2022, retroelementos e introgressão em Afroaves (`41a1c3e3`) — "One response to phylogenomic inaccuracy that is driven by extensive homoplasy in DNA sequence data is to collect character data that are not impacted so severely by multiple overlapping mutations, convergence, reversal, base-compositional shifts, and lineage-specific differences in rate [47,48,49]. Nucleotide sites only have a limited number of states (i.e., G, A, T, C), and with long terminal branches are expected to accrue extensive homoplasy t…"
- [42] Resolving the Avian Tree of Life from Top to Bottom (ResearchGate) (`47d8810c`) — "species tree. The impact of ILS on estimates of the bird tree obtained using standard ML methods also remains uncertain. The theoretical and computational phylogenetics community has put a tremen- dous amount of effort into the development of MSC (“ species tree”) methods over the past decade (reviewed by Edwards 2009; Edwar ds et al. 2016; Liu et al. 2009; Warnow 2018). MSC methods are designed to infer the correct speci es tree given mixtures o…"
- [43] Genes 2022, retroelementos e introgressão em Afroaves (`41a1c3e3`) — "4.7. Estimating Large Species Trees with RIs: Prospects and Challenges Resolving ancient rapid radiations is a major challenge for the field of phylogenomics because supermatrix (concatenation) methods do not explicitly account for ILS that is expected to be rampant in such speciation scenarios [23]. Summary (shortcut) coalescence methods that are based on the MSC are statistically consistent when evolution is completely neutral, loci are not clo…"
- [44] Birds, Diversification of (ResearchGate, capítulo de enciclopédia) (`d3d5a7ff`) — "Alcheringa 47, 315 e 326. McClure, C.J.W., Schulwitz, S.E., Anderson, D.L.,, et al., 2019. De ﬁ ning raptors and birds of prey. J. Raptor Res. 53, 419 e 430. Mirarab, S., Rivas-González, I., Feng, S., et al., 2024. A region of suppressed recombination misleads neoavian phylogenomics. Proc. Natl. Acad. Sci. U. S. A. 121, e2319506121. Musser, G.M., Ksepka, D.T., Field, D.J., 2019. New material of Paleocene-Eocene Pellornis (Aves: Gruiformes) clari …"
- [45] Genes 2022, retroelementos e introgressão em Afroaves (`41a1c3e3`) — "Phylogenomic Coalescent Analyses of Avian Retroelements Infer Zero-Length Branches at the Base of Neoaves, Emergent Support for Controversial Clades, and Ancient Introgressive Hybridization in Afroaves PMCID: PMC9324441 DOI: 10.3390/genes13071167 Fonte: Europe PMC fullTextXML Abstract AbstractRetroelement insertions (RIs) are low-homoplasy characters that are ideal data for addressing deep evolutionary radiations, where gene tree reconstruction e…"
- [46] Genes 2022, retroelementos e introgressão em Afroaves (`41a1c3e3`) — "3.4. Evidence for Introgressive Hybridization Asymmetrical RI support for minority resolutions at the controversial Strigiformes + Cavitaves node within Afroaves was detected using a variety of methods (Figure 10). For the full RI dataset, related indices (Bayesian local PPs, perfectly congruent RI synapomorphies, PCS scores, quartet frequencies) measured an excess of support for Strigiformes + Accipitriformes relative to Cavitaves + Accipitrifor…"
- [47] Suh et al. 2015 (PLoS Biol) (`2d01758b`) — "Complexity of the Neoavian Radiation The probability for the occurrence of ILS depends on N e in relation to the time between consecutive speciation events [29,30], with N e correlating positively and time negatively with the expected extent of ILS, respectively. The observed complex genealogical fates of ancestral RE insertion polymorphisms during the initial super-radiation (Fig 4D) therefore suggest that the onset of neoavian diversification w…"
- [48] Stiller et al. 2024 (Nature, B10K) (`e9c5ae30`) — "Abundant discordance among gene trees Assessing the level of incongruence between gene trees (GTs) across the tree, order-level relationships ranged from showing little or no discordance to high levels of discordance (measured by the quartet score; Fig. 2a). The percentage of GT quartets matching a species-tree branch at the ordinal level ranged from 99.9 to 33.7% (close to one in three, which corresponds to a polytomy). In particular, 14 nodes h…"

---

### Q10 — Filogenias morfológicas e caracteres convergentes que as enganaram

**Pergunta:** How did morphology-based phylogenies place raptors, and which convergent characters misled them?

conversation_id: `1e3b0500-ddf3-4a97-a2f7-0c7f0de0e1ac`

**Resposta do NotebookLM (verbatim):**

> Traditional morphology-based phylogenies consistently united predatory landbirds into a single major evolutionary branch, but modern phylogenomics has revealed that this arrangement was an artifact of extensive evolutionary convergence [1-3].
>
> ---
>
> ### How Morphological Phylogenies Placed Raptors
>
> 1. **Monophyletic Diurnal Raptors (*Falconiformes*)**:
>    Classical anatomical classifications and cladistic studies based on osteology and soft tissues (e.g., Livezey & Zusi 2007) grouped all diurnal birds of prey into a single monophyletic order, **Falconiformes** [4-7]. This included falcons (*Falconidae*), hawks, eagles, and Old World vultures (*Accipitridae*), ospreys (*Pandionidae*), secretarybirds (*Sagittariidae*), and New World vultures (*Cathartidae*) [6, 7].
>
> 2. **Sister Relationship Between Diurnal Raptors and Owls (*Falconimorphae*)**:
>    Morphological analyses frequently placed nocturnal raptors (**Strigiformes** / owls) as the direct sister group to Falconiformes, creating a superordinal raptorial clade (often termed "Falconimorphae" or a unified raptorial grade) [5, 8, 9].
>
> 3. **Placement of New World Vultures**:
>    New World vultures (*Cathartidae*) were positioned morphologically either as the basal-most lineage within Falconiformes [6, 10] or occasionally allied with storks (*Ciconiiformes*) based on superficial postcranial and physiological resemblances [10, 11].
>
> ---
>
> ### The Phylogenomic Disruption
>
> Genomic-scale sequence analyses (e.g., Hackett et al. 2008, Jarvis et al. 2014, Prum et al. 2015, Stiller et al. 2024) dismantled the monophyly of raptors, demonstrating that raptorial birds are polyphyletic and distributed across two major, distantly related landbird clades (**Telluraves**) [1, 3, 12]:
>
> * **Australaves**: Contains falcons (*Falconidae*), which are closely related to parrots (*Psittaciformes*) and songbirds (*Passeriformes*), as well as seriemas (*Cariamiformes*) [1, 12-14].
> * **Afroaves**: Contains Accipitriformes (hawks, eagles, Old World vultures, ospreys, secretarybirds, and New World vultures) as well as Strigiformes (owls, which are sister to or nested near *Coraciimorphae*/*Cavitaves*) [1, 15-17].
>
> ---
>
> ### Convergent Characters That Misled Morphological Trees
>
> The clustering of raptors in morphological trees was driven by strong functional selection for a predatory or scavenging lifestyle, resulting in parallel evolution across multiple skeletal and soft-tissue systems [2, 18]:
>
> * **Beak and Skull Morphology**:
>   * **Sharply Hooked Beak (*Rostrum Maxillare*)**: Decisive, decurved, sharp beaks adapted for tearing flesh and severing prey tissues or neck vertebrae evolved independently across multiple lineages [2, 19-22].
>   * **Supraorbital Processes**: Long, prominent bony projections above the orbits evolved repeatedly in falcons, owls, accipitrids, and extinct fossil lineages (such as *Messelasturidae*) to shield the eyes or anchor jaw musculature [19, 23, 24].
>
> * **Feet, Talons, and Pedal Digit Configurations**:
>   * **Enlarged Claws (*Ungual Phalanges*) and Flexor Tubercles**: Highly hypertrophied, curved claws with prominent flexor tubercles for clutching, pinning down, and killing struggling prey [2, 16, 21, 24].
>   * **Pedal Digit Proportions**: Shortened middle phalanges on Digit IV and modified basal phalanges on Digit II evolved repeatedly to reinforce foot strength under mechanical strain [2, 21, 25].
>   * **Phalangeal Fusion**: Parallel fusion of basal pedal phalanges (e.g., in Digit II) occurred independently across distinct predatory groups [26, 27].
>
> * **Pectoral Girdle and Flight Adaptations**:
>   * Osteological features associated with powerful flapping or precise aerial maneuvering—such as the presence or placement of the *foramen nervi supracoracoideus* in the coracoid, the shape of the acromial process on the furcula, and the sternal crest shape—were repeatedly misread as shared ancestral traits rather than independent locomotor adaptations [22, 24, 28-30].
>
> * **Scavenging Adaptations**:
>   * Unfeathered heads and necks, long bill proportions, and soaring wing structures in New World vultures (*Cathartidae*) and Old World vultures (*Accipitridae*) represent one of the most striking examples of convergent evolution in response to obligate scavenging [7, 10, 31].
>
> ---
>
> 💡 *Would you like to explore how early Eocene fossil taxa (such as Masillaraptoridae or Messelasturidae) help bridge the morphological gaps between falcons, parrots, and stem landbirds?*

**Citações (número → fonte; trecho citado, truncado):**

- [1] Resolving the Avian Tree (ResearchGate, `47d8810c`) — "The Jarvis et al. (2014) TENT… splits core landbirds into two clades that Ericson (2012) named Australaves and Afroaves… That split renders the raptorial lineages para- or polyphyletic…"
- [2] Ksepka et al. 2017 (`b84ec51f`) — "Enforcing our preferred molecular scaffold… results in six MPTs… 31 steps longer… This highlights a potentially high level of homoplasy within some morphological characters, including many related to beak shape and pedal digit organization that contribute to the clustering of the raptorial taxa in the morphological tree. Applying the scaffold breaks up the 'raptorial clade'…"
- [3] Jarvis et al. 2014 (`c92cd31c`) — "…The non-monophyly of the birds of prey at the deepest branches of the Australaves and Afroaves radiations suggests that the common ancestor of core landbirds may have been an apex predator… the Columbea and Passerea clades appear to have many ecologically driven convergent traits that have led previous studies to group them into supposed monophyletic taxa…"
- [4] Livezey & Zusi 2007 (`341dae82`) — "…2954 morphological characters were defined… 97 well-resolved, most-parsimonious trees… A strict consensus tree of MPTs included only six trichotomies…"
- [5] Livezey & Zusi 2007 (`341dae82`) — "Raptors or birds of prey – comprising the diurnal Falconiformes and (principally) nocturnal Strigiformes… The sister-relationship of these raptorial orders inferred herein… and by Mayr et al. (2003) has been the subject of suspicion based on phenetic tallies… and speculations concerning convergences and raptorial specializations… **This clade is first in a sequence of four – the birds of prey, Opisthocomus, Cuculiformes, Psittaciformes and Columbiformes – that are sequential sister-groups of remaining Neornithes**… rendering the branching sequence provisional."
- [6] Livezey & Zusi 2007 (`341dae82`) — "Cathartidae was resolved as the sister-group of other Falconiformes… the Sagittariidae was sister-group of the order exclusive of the Cathartidae. The Falconiformes, exclusive of the foregoing two families, comprised a pair of sister-clades: (i) the Accipitridae, including Old World vultures…, and (ii) **a clade comprising the Pandionidae and its sister-group the Falconidae**, the latter including the caracaras."
- [7] TimeTree of Life cap. 63 (Brown 2009, `3af8052d`) — "The diurnal birds of prey constitute the Order Falconiformes, and are generally classified into five reciprocally monophyletic families: Cathartidae…, Sagittariidae…, Pandionidae…, Accipitridae…, and Falconidae… Falconiform taxa are generally characterized by morphological adaptations to predation, be it active hunting (hooked bills and strong talons) or eating carrion (long necks and unfeathered heads)…"
- [8] Livezey & Zusi 2007 (`341dae82`) — "…The current analysis strongly confirmed a sister-group relationship between the Strigiformes and the Falconiformes (Fig. 16), a union also supported by Cracraft (1988), Mayr & Clarke (2003) and Mayr et al. (2003)… Fossils that exhibit generalized raptorial characters or those of both Strigiformes and Falconiformes also have been described…"
- [9] TimeTree cap. 63 (`3af8052d`) — "The recent comprehensive morphological analysis of Neornithes (10) instead strongly supports a close relationship for Falconiformes and Strigiformes, forming the proposed Superorder 'Falconimorphae.' Published molecular genetic studies with broad taxon and character sampling fail to ally falconiform and strigiform taxa in any arrangement…"
- [10] TimeTree cap. 63 (`3af8052d`) — "…(Family Ciconiidae, Order Ciconiiformes)… neither morphological nor genetic studies aimed at discerning cathartid affinities have succeeded in recovering the Cathartidae–Ciconiidae pairing. An early mtDNA study… included erroneous sequences… and has long since been retracted. In morphological studies that recover a monophyletic Falconiformes, Cathartidae is found to be the basal-most lineage in the clade… there is also no support for a close relationship between New World (Cathartidae) and Old World (Accipitridae) vultures, and thus ecological similarities between them provide a striking example of evolutionary convergence."
- [11] Mather 2021 (`c24ce4bf`) — "…Ligon (1967) proposed that the Cathartidae… were a suborder of the Ciconiiformes (storks) based on certain morphological and physiological similarities. A relationship with the Ciconiiformes has been supported by some studies using molecular data (Sibley and Ahlquist 1990; Avise et al. 1994;…)…"
- [12] AOS 2025-C (`b5e56e8d`) — "Morphological convergence across Accipitridae and even Accipitriformes has long confused our understanding of the relationships among hawks, eagles, and other diurnal raptors… falcons and caracaras (Falconidae) are not actually closely related to the hawks and eagles…, something first suggested by Jollie (1977)…"
- [13] Mayr & Kitchener 2022 (`53dae551`) — introdução (mesmo trecho de Q2 [23]).
- [14] Phylogenetic definitions (`1428d7cb`) — apomorfias de Eufalconimorphae (mesmo trecho de Q2 [11]).
- [15] Birds, Diversification of (`d3d5a7ff`) — "…Phorusrhacidae… among Cariamiformes… Coraciimorphae, which phylogenomic analyses generally recover as close relatives to Accipitriformes and Strigiformes, forming the group Afroaves…"
- [16] Mather 2021 (`c24ce4bf`) — "…Falconidae (falcons), Strigiformes (owls), Pandionidae (ospreys) and occasionally the Cathartidae… Of these families, only the Pandionidae and Cathartidae are closely related in a clade with Accipitridae…"
- [17] Mather 2021 (`c24ce4bf`) — Cathartidae como linhagem mais basal de Accipitriformes (mesmo trecho de Q3 [2]).
- [18] Catanach et al. 2024 (`1d4669b1`) — "…ornithologists have struggled to resolve relationships within Accipitridae because of high levels of apparent convergence in external morphology and skeletal structure, which might be caused by shared selective pressures that are inherent to a predatory lifestyle (Holdaway 1994, Pecsics et al. 2019)."
- [19] Mayr 2021, PalZ (`a681b573`) — mosaico de *Tynskya*: "raptor-like" (bico curto, processos supraorbitais longos) + tarsometatarso "parrot-like".
- [20] Mayr & Hurum 2020 (`32bd8904`) — "Diurnal birds of prey exhibit diverse feeding habits… disparate dietary preferences are reflected by differences in tarsometatarsus and beak proportions… we take the combination of a small size, very long legs, and a comparatively short beak as indicative of an avivorous diet of Aviraptor longicrus…"
- [21] Mayr et al. 2020, *Primoptynx* (`2b83db39`) — "Owls and falcons (Falconidae) immobilize quarry by severing the skull or the neck vertebrae with the beak… Hawks and allies (Accipitridae), by contrast, use their feet and body weight to pin down prey… owls and falcons have pedal claws of equal size, whereas in accipitrids the claws of the first and second digits are hypertrophied… Primoptynx poliotauros… and M. antiqua… agree with species of the Accipitridae in that the ungual phalanges of the first and second toes are greatly enlarged…"
- [22] Mather et al. 2021, *Archaehierax* (`8516fa6d`) — diagnose do fóssil (hipotarso; esterno; úmero; processo procoracóide; rostro com ponta em gancho 30–40°…).
- [23] Mayr 2021 (`a681b573`) — "[supraorbital] processes evolved several times independently in taxa of the Telluraves, which diminishes their phylogenetic value. The present analyses do not support closer affinities between messelasturids and strigiform birds…"
- [24] Psittacomimidae fam. nov. (`5102aa83`) — "…Both taxa agree with the Strigiformes and Falconiformes in various features, such as the long supraorbital processes (which are present in stem group Strigiformes, but were reduced in the crown group taxa…), the asymmetric processus oticus of the quadrate (only Falconiformes), and the presence of a foramen nervi supracoracoidei… Messelasturids have a sharply hooked beak… **However, the flexor tubercles are less developed than in most raptorial birds that use their feet to subdue prey, and it is the extensor tubercle, which is prominent in messelasturids.**"
- [25] Mather 2021 (`c24ce4bf`) — "…Masillaraptor was identified by Mayr as a 'falconiform', as it possesses the shortened middle phalanges of the fourth toe and first phalanx of the…"
- [26] Mather 2021 (`c24ce4bf`) — "fused condition in the basal phalanges of pedal digit II found in some genera of Accipitridae. Olson (1982) found this condition to be present in the milvine accipitrids (Milvus and Haliastur), haliaeetine accipitrids (Haliaeetus and Icthyophaga), Ictinia, and Busarellus…"
- [27] Mather 2021 (`c24ce4bf`) — "fused phalanges condition is actually convergent, and that Ictinia and Busarellus do not form a clade with the milvine and haliaeetine accipitrids (Holdaway 1994; Lerner et al. 2008; Nagy and Tökölyi 2014)…"
- [28] Ikovo, coracoide de Telluraves (`6ccf22e4`) — "…The processus lateralis is reduced and has an obtuse angulus lateralis, like in Falconidae and Strigidae… The impressio musculi sternocoracoidei is very shallow and has no pneumatic foramina, unlike in Cathartidae…"
- [29] Chen et al. 2025 (`6c89ed1d`) — sinapomorfias potenciais de Accipitrimorphae (mesmo trecho de Q3 [1]).
- [30] Catanach et al. 2024 (`1d4669b1`) — "…Olson (1987) noted that the procoracoid foramen, which is absent in 'Accipiter' but present in nearly all other hawk genera, is also absent in some members of Harpagus ('kites') and Circus (harriers)… Holdaway 1994… concluded that Accipiter and Harpagus were closely related to the exclusion of Circus…"
- [31] Wang et al. 2014, GBE (`3cc767df`) — "Falconidae (falcons) appear to have shared a common ancestor with SPPC (PP = 0.92)… We never recovered a direct sister relationship between accipitrids and falconids… Cathartidae is possibly not a genuine raptor and its resemblance to Old World vultures is an example of convergent evolution (Tagliarini et al. 2009)…"

---


---

## Briefing Doc (Studio → report, formato "Briefing Doc")

- **Artefato:** `6c992653-30ad-458f-9f71-f42dc8362474`, gerado com as 41 fontes com conteúdo e baixado em `C:\Users\Manec\Downloads\gemini-notebook\aves-de-rapina-2026-10-17\briefing.md`.
- **Atenção:** o formato "Briefing Doc" **não aceita instrução de foco** (o `custom_prompt` só vale para "Create Your Own"). O resultado ficou **genérico sobre a evolução das aves** e é dominado pelo capítulo de enciclopédia de Chen et al. 2026: origem nos dinossauros e aves mesozoicas. Sobre rapinantes, quase só traz uma linha de tabela. Fiz uma segunda tentativa de "Briefing Doc" com prompt de foco (`5ba59a20-78b8-4afd-93f4-22c5240bc973`, `...\briefing-focado.md`), que também saiu genérica e **não vai transcrita aqui**. Por isso gerei também um relatório "Create Your Own" com o formato de briefing e foco em rapinantes, transcrito depois deste.

Texto verbatim do Briefing Doc:

### Comprehensive Briefing Document: Avian Diversification, Phylogenomics, and Evolutionary History

#### Executive Summary

This briefing document synthesizes current scientific knowledge regarding the evolutionary history, phylogenetic relationships, and diversification timeline of birds (**Aves**). Originating within maniraptoran theropod dinosaurs during the Mesozoic Era, birds represent the sole surviving lineage of dinosaurs. Over more than 150 million years, avian evolution transitioned from early, non-volant feathered theropods to specialized Mesozoic flyers, surviving the end-Cretaceous (**K–Pg**) mass extinction (~66 million years ago) and radiating into the ~11,000 extant species observed today.

Recent integration of paleontology, comparative anatomy, and genome-scale sequencing (**phylogenomics**) has clarified long-standing controversies surrounding avian origins and higher-level systematics. While early molecular studies suggested deep Cretaceous origins for modern bird orders, recent tip-dating and phylogenomic supertree syntheses align closer to the fossil record: crown birds (**Neornithes**) originated during the mid-to-late Cretaceous (~100 Ma), but the major ordinal diversification of modern birds (**Neoaves**) was catalyzed by the K–Pg extinction event. This rapid radiation was facilitated by ecological opportunities opened up by the collapse of Mesozoic ecosystems and global forest habitats.

Despite major phylogenomic breakthroughs, recalcitrant nodes in the avian tree of life persist due to rapid speciation intervals, incomplete lineage sorting, convergent evolution, and complex genomic phenomena such as ancient chromosomal rearrangements with suppressed recombination. Resolving these relationships requires combining high-throughput sequence capture, comprehensive taxon sampling, fossil calibrations, and computationally efficient supertree methods.

---

#### Detailed Analysis of Key Themes

##### 1. Dinosaurian Origin and Avian Ancestry

###### Theropod Ancestry and Miniaturization
Rigorous anatomical and phylogenetic analyses demonstrate that birds are living theropod dinosaurs. Specifically, birds belong to **Avialae**, a clade nested within **Paraves** (which also includes dromaeosaurids such as *Deinonychus* and *Velociraptor*, and troodontids) inside the broader group **Maniraptora**. 

The evolution of birds was driven by sustained evolutionary miniaturization (body size reduction) within meat-eating theropod dinosaurs. While body mass decreased, paravian dinosaurs maintained forelimb proportions, providing early avialans with relatively long forelimbs that facilitated the evolution of powered flight.

```
Dinosauria
 └── Theropoda
      └── Maniraptora
           └── Paraves
                ├── Deinonychosauria (Dromaeosauridae & Troodontidae)
                └── Avialae (Stem & Crown Birds)
```

###### Pre-Avian Assembly of "Bird" Traits
Fossil discoveries—particularly exceptionally preserved specimens from China discovered since the late 1990s—demonstrate that iconic "bird" features evolved long before powered flight in non-avian theropod ancestors:
* **Feathers:** Filamentous and pinnate feather structures initially evolved for non-flight functions (such as thermoregulation, display, or brooding) in non-avian theropods. Modern asymmetrical flight feather anatomy emerged later within volant lineages.
* **Skeletal Modifications:** The wishbone (**furcula**) and the specialized **semilunate carpal** (wrist bone) originated in non-avian theropods.
* **Physiological and Behavioral Traits:** Respiratory air sacs, avian-like alimentary and reproductive structures (including partially buried eggs and contact incubation precursors), and stereotypical avian sleeping postures (documented in troodontids like *Mei* and *Sinornithoides*, and the alvarezsaurid *Jaculinykus*) were already present in non-avian maniraptorans.

###### Aerodynamic Experimentation
Flight potential evolved across multiple paravian lineages during the Jurassic and Cretaceous periods. For instance, microraptorian dromaeosaurids possessed elongate aerodynamic feathers on both forelimbs and hindlimbs ("four-winged" configurations), utilizing leg feathers as supplemental lifting surfaces. Debate remains regarding whether avian flight originated through ground-up takeoff or tree-down gliding, and whether flight evolved multiple times or was secondarily lost in certain non-avian maniraptorans.

---

##### 2. Diversity and Evolution of Mesozoic (Pre-Modern) Birds

The vernacular term "bird" generally corresponds to **Avialae**—the lineage encompassing modern birds and all extinct taxa more closely related to modern birds than to dromaeosaurids and troodontids.

```
Mesozoic Avialan Lineages
 ├── Non-Volant / Weakly Volant Early Avialans (Anchiornithidae, Archaeopteryx, Jeholornis)
 └── Short-Tailed Avialans (Pygostylia)
      ├── Primitive Short-Tailed Forms (Sapeornis, Confuciusornithiformes)
      └── Ornithothoraces
           ├── Enantiornithes ("Opposite Birds" - Extinct at K-Pg)
           └── Euornithes
                ├── Basal Euornithines (Gansus, Yanornithidae)
                └── Ornithurae
                     ├── Ichthyornithes & Hesperornithes (Toothed Stem Groups)
                     └── Neornithes (Crown Birds)
```

###### Early Avialans
* **Anchiornithidae:** Middle-to-Late Jurassic forms (e.g., *Anchiornis*) representing the oldest known avialans; they retained non-avialan paravian traits and lacked powered flight capabilities.
* **Long-Tailed Avialans:** Late Jurassic *Archaeopteryx* and Early Cretaceous *Jeholornis* retained long, bony tails and possessed weak flight capabilities.
* **Early Short-Tailed Forms:** Shortened bony tails ending in a pygostyle appeared in the Early Cretaceous, represented by *Sapeornis*, Jinguofortisidae, and *Confuciusornithiformes* (known from thousands of Chinese fossils).

###### Ornithothoraces and the Enantiornithes vs. Euornithes Split
Further structural adaptations for powered flight arose in **Ornithothoraces**, characterized by an enlarged, keeled sternum for flight muscle attachment and upward-facing shoulder sockets for higher wing elevation:
1. **Enantiornithes ("Opposite Birds"):** Highly diverse and abundant in Cretaceous terrestrial ecosystems with over 100 described species. They retained teeth and clawed wings and occupied diverse ecological niches, but became completely extinct at the end of the Cretaceous.
2. **Euornithes:** The lineage that led to modern birds. They evolved an aerodynamic fan of tail feathers controlled by tail musculature, along with toothless jaw regions that formed an early avian beak.
   * **Cretaceous Euornithines:** Included terrestrial and semi-aquatic forms such as *Archaeorhynchus*, *Gansus*, *Iteravis*, and secondarily flightless species like *Patagopteryx*.
   * **Ornithurae:** Subgroup containing modern birds (**Neornithes**) and their closest extinct relatives, such as the toothed, tern-like *Ichthyornis* and the flightless, foot-propelled diving *Hesperornithes*.

---

##### 3. K–Pg Extinction Dynamics and the "Rocks vs. Clocks" Timeline Controversy

###### Survival Across the K–Pg Boundary
All toothed birds vanished during the end-Cretaceous mass extinction ~66 million years ago. Survival of early modern bird lineages (**Neornithes**) across the K–Pg boundary was facilitated by specific ecological and physiological advantages over other Mesozoic avialans:
* **Dietary Flexibility:** Ability to consume seeds and hard foods during global ecosystem disruption.
* **Habitats and Forest Collapse:** Global forest destruction favored ground-dwelling or water-associated birds over arboreal enantiornithes.
* **Body Size and Reproduction:** Smaller body sizes (an avian "Lilliput Effect") and distinct reproductive strategies (such as fully subaerial contact incubation).

###### Harmonizing Fossils and Molecular Clocks
Historical discordance between fossil dates ("Rocks") and molecular clock estimates ("Clocks") has narrowed significantly:
* **Strict Cretaceous Radiation Hypothesis (Rejected):** Older molecular estimates proposing deep Early Cretaceous origins (130 Ma) or widespread Cretaceous diversification of modern orders are unsupported by tip-dating combined-evidence analyses and Bayes factor model tests.
* **Emerging Scientific Consensus:** 
  * Crown **Neornithes** originated in the mid-to-late Cretaceous (~100 Ma), marking the primary split between **Palaeognathae** and **Neognathae**.
  * Very few modern bird lineages are documented in Late Cretaceous fossils (e.g., *Vegavis* from Antarctica, *Asteriornis* and *Janavis* from Europe).
  * The explosive radiation of major **Neoavian** orders occurred rapidly in the immediate aftermath of the K–Pg mass extinction (~66 Ma), driven by cleared ecological space, followed by further interfamilial diversification around the Eocene–Oligocene transition (~34 Ma).

---

##### 4. High-Level Phylogeny and Systematics of Modern Birds (Neornithes)

Extant birds are divided into two main superorders: **Palaeognathae** and **Neognathae**.

```
Neornithes (Crown Birds)
 ├── Palaeognathae (Tinamous & Ratites)
 └── Neognathae
      ├── Galloanserae (Landfowl & Waterfowl)
      └── Neoaves (8 Major Clades)
```

###### Summary of Major Avian Lineages and Representative Taxa

| Major Lineage | Sub-Clades / Orders Included | Key Characteristics & Representative Extinct/Extant Taxa |
| :--- | :--- | :--- |
| **Palaeognathae** | Struthioniformes, Rheiformes, Casuariiformes, Apterygiformes, Aepyornithiformes, Dinornithiformes, Tinamiformes | Flightless "ratites" are **polyphyletic**; flightlessness evolved independently multiple times. Volant tinamous are nested within ratites. Early volant forms include *Lithornithidae*. |
| **Galloanserae** | Galliformes (Landfowl), Anseriformes (Waterfowl) | Basal neognaths. Oldest known fossils include Late Cretaceous *Asteriornis* and *Vegavis*, and Paleocene *Conflicto*. Includes giant extinct flightless forms (*Gastornithidae*, *Dromornithidae*, *Sylviornithidae*). |
| **Neoaves: Mirandornithes** | Phoenicopteriformes (Flamingos), Podicipediformes (Grebes) | Highly specialized aquatic birds; ancestral morphology illustrated by Eocene *Juncitarsus*. Phylogenomic studies place them sister to all other Neoaves. |
| **Neoaves: Columbaves** | Columbimorphae (Pigeons, Sandgrouse, Mesites), Otidimorphae (Bustards, Cuckoos, Turacos) | Poor early fossil record; Eocene *Foro panarium* represents an early stem turaco. |
| **Neoaves: Strisores** | Caprimulgiformes, Steatornithiformes, Nyctibiiformes, Podargiformes, Aegotheliformes, Apodiformes | Primarily nocturnal/crepuscular insectivores and diurnal swifts/hummingbirds. Eocene fossils (*Eocypselus*, *Scaniacypselus*) document wing and bill evolution. |
| **Neoaves: Opisthocomiformes** | Opisthocomidae | Single extant species: Hoatzin (*Opisthocomus hoazin*) of South America (folivorous, juvenile wing claws). Miocene African fossils suggest non-South American origin. |
| **Neoaves: Gruiformes** | Gruoidea (Cranes, Trumpeters), Ralloidea (Rails, Finfoots, Adzebills) | Eocene *Pellornis*. Secondary flightlessness evolved dozens of times independently (e.g., extinct New Zealand adzebills, *Aptornithidae*). |
| **Neoaves: Charadriiformes** | Charadrii (Plovers), Scolopaci (Sandpipers), Lari (Gulls, Auks) | Shorebirds, pelagic, and arid-adapted birds. Paleogene record sparse; Eocene *Scandiavis* is an early representative. |
| **Neoaves: Phaethoquornithes** | Phaethontimorphae (Tropicbirds, Sunbittern, Kagu), Aequornithes (Loons, Penguins, Tubenoses, Storks, Pelicans, Herons) | Core waterbirds and tropicbirds. Includes Paleocene *Waimanu*, giant extinct penguins (*Kumimanu fordycei*, ~148–160 kg), giant bony-toothed *Pelagornithidae* (>6m wingspan), and convergent flightless diving *Plotopteridae*. |
| **Neoaves: Telluraves** | Afroaves (Accipitriformes, Strigiformes, Coraciimorphae) & Australaves (Cariamiformes, Falconiformes, Psittacopasseres) | Core landbirds. Includes multiple raptorial lineages (*Teratornithidae*, *Phorusrhacidae* "terror birds"), Paleocene *Tsidiiyazhi* (oldest telluravian), and **Passeriformes** (~60% of all bird species, originating in Australasia from zygodactyl ancestors). |

---

##### 5. Methodological Frameworks in Avian Phylogenomics

###### Supertree vs. Supermatrix Methodologies
Constructing a complete species-level avian tree requires processing massive genomic datasets. Two primary analytical strategies are utilized:

```
Phylogenomic Approaches
 ├── Supermatrix Methods
 │    ├── Pros: Direct alignment of sequence characters; detects hidden support.
 │    └── Cons: Extremely sparse matrices, "terraces in treespace," high compute burden (>400 CPU years).
 └── Supertree Methods (MRP / MRL)
      ├── Pros: Highly efficient (<1 week compute), integrates diverse datasets & backbones.
      └── Cons: Decoupled from raw sequence sites; sensitive to missing taxon overlap and backbone conflict.
```

1. **Supermatrix Approaches:** Concatenate raw sequence alignment character data across sampled loci. 
   * *Limitations:* Sparse data matrices with uneven locus overlap create "terraces in phylogenetic treespace" where thousands of distinct tree topologies exhibit identical optimization scores. Furthermore, whole-genome supermatrices impose massive computational demands (e.g., >400 CPU years for whole-genome datasets).
2. **Supertree Approaches (MRP and MRL):** Combine previously estimated subtrees (source trees) into a single overarching topology using algorithms like **Matrix Representation with Parsimony (MRP)** or **Matrix Representation with Likelihood (MRL)**.
   * *Advantages:* Supertrees do not require overlapping loci among studies and are computationally efficient (requiring less than a week on standard desktop computers for 707+ taxa). They allow integration of target capture data, whole-genome alignments, and taxon-rich megaphylogenies (e.g., *Jetz*, *BigBird*, and *Brown* backbones).
   * *Branch Length & Time Calibration:* Because supertrees lack intrinsic sequence substitution rates, mitochondrial gene alignments (such as *CYB* and *ND2*) are extracted to calculate branch lengths. Autocorrelated penalized likelihood (via programs like `treePL`) calibrated with vetted fossil dates can then generate time-calibrated trees.

###### Analytical Challenges and Genomic Discordance
* **Hard Polytomies and Discordance:** Neoaves underwent an extremely rapid radiation near the K–Pg boundary. Short inter-node intervals caused extensive **incomplete lineage sorting (ILS)** and ancient hybridization, leaving recalcitrant nodes across the tree.
* **Genomic Anomaly Zones:** A 21-megabase region mapped to chicken chromosome 4 exhibits suppressed recombination that persisted for millions of years. This single locus carries a strong, discordance-free signal that differs from the true species tree, misleading phylogenomic analyses that lack sufficient taxon sampling.
* **Taxonomic Flux:** Reconciling species names across historical datasets represents a major logistical bottleneck. Ongoing taxonomic splits emphasize the necessity of linking all genomic sequences directly to museum vouchers.

---

#### Important Quotes with Context

##### 1. On Birds as Living Dinosaurs and Paravian Miniaturization
> **Quote:** *"Current analyses of comparative anatomy that have included the known varieties of feathered and non-feathered dinosaurs, modern and ancient birds, and other ancient reptile lineages, consistently show that birds are descendants of a particular lineage of dinosaurs, the maniraptorans... Thus, based on current scientific consensus, birds are living dinosaurs."*
* **Context:** Derived from *Chen et al. (2026)*, this passage encapsulates the scientific consensus regarding avian origins. It highlights that defining "bird" traits (such as feathers, furculae, and wrist adaptations) originated deep within theropod dinosaurs, and that modern birds represent a surviving lineage of miniaturized maniraptoran dinosaurs.

##### 2. On Post-K–Pg Extinction Diversification Dynamics
> **Quote:** *"We discovered sharp increases in effective population size, substitution rates and relative brain size following the Cretaceous–Palaeogene extinction event, supporting the hypothesis that emerging ecological opportunities catalysed the diversification of modern birds."*
* **Context:** From *Stiller et al. (2024)* (referenced in ResearchGate citations), this finding provides genomic proof for the ecological release hypothesis. Following the end-Cretaceous extinction of non-avian dinosaurs and archaic toothed birds, surviving crown bird lineages expanded into vacant ecological niches, driving rapid morphological, neural, and genetic evolution.

##### 3. On Resolving "Rocks vs. Clocks" Divergence Dates
> **Quote:** *"Our results allow us to reject scenarios of crown bird origins deep in the Early Cretaceous, as well as an origin of crown birds within the last 10 million years of the Cretaceous. The scenario best supported by our analyses is one where Neornithes originated between the Early and Late Cretaceous (~ 100 MYA), while numerous divergences within major neoavian clades either spans or postdates the K–Pg."*
* **Context:** 来自 *Brocklehurst & Field (2024)*, this study uses tip-dating and Bayes Factors to reconcile molecular clock estimates with morphological fossil datasets. It refutes early molecular estimates claiming ancient Mesozoic ordinal diversification, placing crown origins at ~100 Ma and attributing major ordinal expansions to the post-K–Pg period.

##### 4. On Computational Efficiency of Supertree Construction
> **Quote:** *"Supertrees get around many of the limitations of supermatrices; they do not require overlapping loci and they can be very computationally efficient... The analyses for this study were conducted on a desktop computer... and they used less than one week of total compute time, as opposed to the >400 years used to analyze the Jarvis et al. data."*
* **Context:** From *Kimball et al. (2019)*, this quote underscores the practicality of supertree methods (MRP/MRL) in synthesizing large-scale phylogenomic studies. While whole-genome supermatrices require massive computing clusters, supertrees quickly weave disparate target-capture and genomic trees together using megaphylogeny backbones.

##### 5. On Misleading Genomic Signals from Ancient Chromosomal Events
> **Quote:** *"We found that a 21-Mb region in avian genomes, mapped to chicken chromosome 4, shows an extremely strong and discordance-free signal for a history different from that of the inferred species tree... We show that the presence of this region has misled previous phylogenomic efforts with lower taxon sampling..."*
* **Context:** From *Mirarab et al. (2024)*, this revelation explains why different high-impact phylogenomic studies obtained contradictory Neoavian tree topologies. A suppressed recombination event 65 million years ago locked a large genomic block into a non-species-tree topology, demonstrating that locus sampling alone without broad taxon sampling can bias phylogenomic inference.

---

#### Actionable Insights

##### For Systematic and Evolutionary Biologists
1. **Prioritize Locus Quality and Analytical Models over Raw Matrix Size:** To resolve deep, rapid radiations within Neoaves, studies must account for incomplete lineage sorting, rate heterogeneity, and genomic regions with suppressed recombination (such as the 21-Mb region on chromosome 4). Combining coalescent-based species tree methods with structural genomic partitions is essential.
2. **Standardize Fossil Calibrations:** Divergence time estimations should strictly adopt vetted fossil calibrations with explicitly justified minimum and maximum bounds (such as *Asteriornis*, *Vegavis*, *Tsidiiyazhi*, *Waimanu*, and *Lithoptila*) following established paleontological protocols.
3. **Expand Targeted Sequencing in Data-Deficient Groups:** Large-scale initiatives like **OpenWings** (utilizing ultraconserved elements and mitochondrial sequences from museum specimens) offer a rapid pathway to complete a species-level tree of life, complementing whole-genome efforts (**B10K**).

##### For Computational Biologists and Software Developers
1. **Deploy Supertree Frameworks for Large Meta-Analyses:** Given the exponential growth of reduced-representation (UCE) and whole-genome datasets, supertree construction algorithms (e.g., MRP, MRL, ASTRAL) should be leveraged to rapidly update global phylogenies without re-analyzing massive, sparse supermatrices.
2. **Automate Taxonomic Reconciliation Pipelines:** Develop automated tools to track sample provenance, mapping historical species names, splits, and lumps directly to persistent museum voucher identifiers (e.g., frozen tissues or study skins) to eliminate taxonomic reconciliation errors during meta-analyses.

##### For Conservationists and Policy Makers
1. **Incorporate Deep Phylogenetic Diversity into Conservation Priorities:** Evolutionary lineages representing isolated, ancient Neoavian branches with low species richness (such as the Hoatzin [*Opisthocomus hoazin*], Kagu, Sunbittern, and Seriemas) possess unique evolutionary history and require targeted conservation funding.
2. **Mitigate Accelerated Modern Extinction Rates:** Projections indicate substantial loss of avian biodiversity over the next 50 years due to climate change and habitat destruction. Conservation strategies must prioritize families at highest risk of extinction, leveraging detailed phylogenetic maps to preserve global evolutionary potential.
---

## Relatório focado em rapinantes (Studio → report, formato "Create Your Own" com estrutura de briefing)

- **Artefato:** `bee72a25-62ac-44b9-8810-4662e05b96b4`, mesmas 41 fontes, baixado em `C:\Users\Manec\Downloads\gemini-notebook\aves-de-rapina-2026-10-17\briefing-rapinantes.md`.
- **Prompt:** briefing restrito a aves de rapina: posição de Cathartiformes, Accipitriformes, Falconiformes e Strigiformes em Telluraves, suporte por estudo, relações internas, fósseis do Paleógeno, datações, hipótese do ancestral raptorial e nós não resolvidos. Pedi atribuição de cada afirmação a uma fonte nomeada e aviso explícito quando o texto viesse de artigo citante reproduzido numa página do ResearchGate.
- **Atenção:** o gerador se apoiou quase só em Chen et al. 2026, Kimball et al. 2019 e nos trechos das páginas do ResearchGate. Ele chega a afirmar que Masillaraptoridae, Messelasturidae e *Archaehierax* "do not appear in the primary ground-truth source text", o que é **falso** para o caderno (ver auditoria). Use com cautela.

Texto verbatim:

### Technical Briefing Document: Phylogenetic Position, Fossil History, and Evolutionary Dynamics of Birds of Prey (Telluraves)

---

##### 1. Executive Summary

The evolutionary history and phylogenetic relationships of birds of prey within the core landbird clade (**Telluraves**) represent one of the most dynamic and heavily scrutinized topics in modern avian systematics. Comparative genomic studies and fossil evidence confirm that modern birds (*Neornithes*) are living maniraptoran theropod dinosaurs that survived the Cretaceous-Paleogene (K-Pg) mass extinction event (~66 Ma). Within Neoaves, Telluraves encompasses the vast majority of arboreal and terrestrial species. However, high-density phylogenomic analyses consistently establish that classical birds of prey do not form a monophyletic group. Instead, raptorial lineages are polyphyletically distributed across two major landbird super-clades: **Afroaves** (comprising **Accipitriformes** [hawks, eagles, Old World vultures, ospreys, secretarybirds] and **Strigiformes** [owls]) and **Australaves** (comprising **Falconiformes** [falcons, caracaras] and **Cariamiformes** [seriemas and extinct terror birds]).

This profound structural dispersion has catalyzed a primary macroevolutionary debate in neoavian biology: whether the suite of morphological and ecological adaptations defining birds of prey—including hooked beaks, sharp talons, forward-facing eyes, and carnivorous habits—evolved independently through convergent evolution across distinct lineages, or whether predatory, vertebrate-eating habits represent the primitive, ancestral condition of all core landbirds (the **raptorial-ancestor hypothesis**). Phylogenomic frameworks recovering highly speciose, non-raptorial groups (such as Coraciimorphae within Afroaves and Psittacopasseres within Australaves) intercalated between raptorial orders indicate that if raptoriality were the ancestral landbird condition, secondary losses of predatory traits must have occurred repeatedly across songbirds, parrots, woodpeckers, and mousebirds.

Resolving these deep node topologies has historically been hindered by the explosive radiation of landbirds near or immediately following the K-Pg boundary. Early Paleocene fossil discoveries, such as the stem mousebird *Tsidiiyazhi abini* (~62.2 Ma) and the stem owl *Ogygoptynx wetmorei* (~58.5 Ma), provide mandatory calibration points that establish tight minimum divergence bounds for Telluravian diversification. Concurrently, advancements in genomic sequencing—utilizing ultraconserved elements (UCEs), whole-genome alignments, and multispecies coalescent frameworks—have exposed critical biological and analytical hurdles, including extensive incomplete lineage sorting (ILS), variable substitution rates, hemiplasy, and structural genomic anomalies like non-recombining chromosomal regions that distort concatenated phylogenetic trees.

***

###### Key Takeaways
*   **Polyphyletic Raptor Distribution:** "Birds of prey" do not constitute a natural evolutionary group; raptorial lineages are divided across **Afroaves** (*Accipitriformes* and *Strigiformes*) and **Australaves** (*Falconiformes* and *Cariamiformes*).
*   **The Raptorial-Ancestor Debate:** The basal placement of raptorial orders in both Afroaves and Australaves supports the hypothesis of a predatory Telluravian ancestor, though this model necessitates widespread secondary losses across non-raptorial clades like songbirds, parrots, and woodpeckers.
*   **Rapid Post-K-Pg Radiation:** Core landbird diversification was concentrated in a tight temporal window following the end-Cretaceous extinction (~66 Ma), anchored by early Paleocene calibration fossils (~62.2 Ma).
*   **Genomic Discordance & Hemiplasy:** Short internodes during rapid post-K-Pg speciation produced severe incomplete lineage sorting and hemiplasy, while localized genomic phenomena (such as a 21-Mb non-recombining block on chicken chromosome 4) systematically mislead concatenated tree estimation.

---

##### 2. Higher-Level Systematics: Raptor Lineages within Telluraves

The clade **Telluraves** ("core landbirds") encompasses the most speciose radiation within Neoaves (*Chen, van Tuinen, & Field, 2026*; *Kimball et al., 2019*). Modern phylogenomic consensus divides Telluraves into two principal superordinal subclades: **Afroaves** and **Australaves**.

###### Afroaves Systematics
Afroaves unites two major raptorial orders with the species-rich, non-raptorial clade **Coraciimorphae** (comprising mousebirds [Coliiformes], cuckoo-rollers [Leptosomiformes], trogons [Trogoniformes], hornbills and hoopoes [Bucerotiformes], rollers and kingfishers [Coraciiformes], and woodpeckers and toucans [Piciformes]):
1.  **Accipitriformes:** Comprises hawks, eagles, and Old World vultures (**Accipitridae**), ospreys (**Pandionidae**), and the terrestrial secretarybird (**Sagittariidae**). Modern phylogenomic analyses also nest or closely align New World vultures (**Cathartidae** / **Cathartiformes**) within or sister to this order, alongside the extinct giant soaring teratorns (**Teratornithidae**) (*Chatterjee et al., 2007*, cited in *Chen et al., 2026*).
2.  **Strigiformes:** True owls (**Strigidae**) and barn owls (**Tytonidae**).
3.  *Taxonomic Note on "Hieraves":* In secondary systematic literature, the sister-group relationship between Accipitriformes and Strigiformes within Afroaves is frequently designated under the node label **Hieraves**. While the sister relationship between Accipitriformes and Strigiformes (and their collective sister relationship to Coraciimorphae) is explicitly supported by the primary text of *Chen et al. (2026)*, the formal clade name "Hieraves" represents external secondary terminology not explicitly utilized in *Chen et al. (2026)*.

###### Australaves Systematics
Australaves links two raptorial orders with **Psittacopasseres** (parrots [Psittaciformes] and perching birds/songbirds [Passeriformes]):
1.  **Falconiformes:** Comprising falcons and caracaras (**Falconidae**). Phylogenomic studies demonstrate that falcons are not closely related to hawks or eagles, but are instead the sister group to Psittacopasseres, together forming the superordinal clade **Eufalconimorphae**.
2.  **Cariamiformes:** Extant seriemas (**Cariamidae**) and extinct flightless apex predators, including the terror birds (**Phorusrhacidae**). Cariamiformes represents an early-diverging lineage within Australaves, positioned sister to Eufalconimorphae.

| Order / Taxon | Clade Membership | Sister Group | Underlying Dataset Type |
| :--- | :--- | :--- | :--- |
| **Accipitriformes** | Afroaves | Strigiformes (or Coraciimorphae) | Phylogenomic UCE, Whole Genome, Supertree |
| **Strigiformes** | Afroaves | Accipitriformes (or Coraciimorphae) | Phylogenomic UCE, Whole Genome, Supertree |
| **Cathartiformes** | Afroaves (Accipitriformes *s.l.*) | Accipitridae / Teratornithidae | Phylogenomic UCE, Morphological, Legacy DNA |
| **Falconiformes** | Australaves / Eufalconimorphae | Psittacopasseres | Phylogenomic UCE, Whole Genome |
| **Cariamiformes** | Australaves | Eufalconimorphae (*Falconiformes* + *Psittacopasseres*) | Phylogenomic UCE, Whole Genome, Morphological |
| **Coraciimorphae** *(Non-Raptorial)* | Afroaves | Accipitriformes + Strigiformes | Phylogenomic UCE, Whole Genome |
| **Psittacopasseres** *(Non-Raptorial)* | Australaves / Eufalconimorphae | Falconiformes | Phylogenomic UCE, Whole Genome |

---

##### 3. Internal Systematics of Primary Raptor Families

Modern phylogenomic investigations detail the internal branch topologies within primary extant raptor families and their positions relative to allied landbird groups (*Chen et al., 2026*; *Kimball et al., 2019*).

*   **Order Accipitriformes**
    *   **Family Sagittariidae**
        *   Taxon: Secretarybird (*Sagittarius serpentarius*).
        *   Systematic Position: Represents the earliest-diverging extant lineage within total-clade Accipitriformes.
    *   **Family Pandionidae**
        *   Taxon: Osprey (*Pandion haliaetus*).
        *   Systematic Position: Forms the well-supported sister group to core Accipitridae.
    *   **Family Accipitridae**
        *   Taxa: Diurnal hawks, eagles, and Old World vultures.
        *   Systematic Position: Core monophyletic family within Accipitriformes; Old World scavenging vultures are nested directly within accipitrid diversity rather than forming an independent order.
    *   **Clade / Order Cathartiformes**
        *   Taxa: New World vultures (**Cathartidae**) and extinct teratorns (**Teratornithidae**).
        *   Systematic Position: Closely affiliated with or nested within Accipitriformes *s.l.*; teratorns represent giant Cenozoic soaring predators/scavengers sharing close affinities with cathartids (*Chatterjee et al., 2007*).
*   **Order Falconiformes**
    *   **Family Falconidae**
        *   Taxa: Falcons and caracaras (subfamilies Polyborinae and Falconinae).
        *   Systematic Position: The sole extant family of Falconiformes.
    *   **Clade Eufalconimorphae Structure**
        *   Systematic Position: Falconidae forms the sister group to **Psittacopasseres**:
            *   **Order Psittaciformes:** Parrots, possessing zygodactyl feet.
            *   **Order Passeriformes:** Perching birds and songbirds, accounting for ~60% of all living avian species. Developmental and fossil evidence indicates passeriforms evolved from ancestors with zygodactyl feet retained by parrots (*Botelho et al., 2014*; *Mayr, 2015*, cited in *Chen et al., 2026*).
*   **Order Strigiformes**
    *   **Family Tytonidae**
        *   Taxa: Barn owls.
        *   Systematic Position: Represents an early divergence within crown Strigiformes.
    *   **Family Strigidae**
        *   Taxa: True owls / typical owls.
        *   Systematic Position: Encompasses the majority of extant nocturnal raptor diversity.
    *   **Afroavian Structural Linkage**
        *   Systematic Position: Strigiformes is recovered as the sister group to Accipitriformes within Afroaves, positioned relative to the non-raptorial Coraciimorphae.

---

##### 4. Paleogene Raptor Fossil Record and Stem Lineages

The fossil record of Paleogene landbirds provides direct anatomical evidence establishing the temporal boundaries, morphological transitions, and body-size dynamics of early raptorial lineages (*Chen et al., 2026*; *Kimball et al., 2019*).

```
+-----------------------------------------------------------------------------------+
|                            PALEOGENE LANDBIRD TIMELINE                            |
|                                                                                   |
|  Paleocene (~66-56 Ma)           Eocene (~56-34 Ma)          Oligocene-Pleistocene  |
|  --------------------           -------------------          ---------------------|
|  • Tsidiiyazhi (~62.2 Ma)        • Eozygodactylus (~51.8 Ma)  • Giant Teratorns    |
|    Stem Mousebird (Telluraves)     Stem Passeriform             (Teratornithidae) |
|  • Ogygoptynx (~58.5 Ma)         • Phorusrhacidae             • Phorusrhacidae     |
|    Oldest Afroaves / Stem Owl      "Terror Birds" Apex         Apex Predators     |
|    (Berruornis, Primoptynx)        Predators (Cariamiformes)   in Americas        |
+-----------------------------------------------------------------------------------+
```

###### Stem Owls (*Strigiformes*)
The oldest documented representative of Afroaves in the fossil record is the Paleocene stem owl ***Ogygoptynx wetmorei*** (~58.5 Ma; cited in *Kimball et al., 2019*, Appendix A.18). Additional early Paleogene stem owls, including *Berruornis* and *Primoptynx*, demonstrate that specialized nocturnal raptorial adaptations were already established in the early Cenozoic.

###### Giant Extinct Predators of Australaves and Afroaves
*   **Phorusrhacidae ("Terror Birds"):** Predominantly flightless, giant apex predators within Cariamiformes known from the Cenozoic of the Americas (*Alvarenga & Höfling, 2003*; *Angst & Buffetaut, 2017*, cited in *Chen et al., 2026*). They demonstrate extreme ecological and body-size diversification in early Australaves.
*   **Teratornithidae:** Giant soaring predatory/scavenging birds from the Oligocene to Pleistocene of the Americas (*Chatterjee et al., 2007*, cited in *Chen et al., 2026*), affiliated with Cathartiformes / Accipitriformes.

###### Clarification on Unmentioned Fossil Taxa
*   *Masillaraptoridae* & *Messelasturidae*: While frequently discussed in external literature as Eocene stem landbirds or stem raptors, these two families do not appear in the primary ground-truth source text (*Chen et al., 2026*; *Kimball et al., 2019*).
*   *Archaehierax*: The Australian Eocene/Oligocene accipitrid *Archaehierax* is likewise unmentioned in the provided primary source text.

###### Bounding Calibration Anchors for Telluravian Divergence
To establish the temporal framework of raptor evolution, systematists rely on stratigraphically constrained non-raptorial landbird fossils that anchor key bounding nodes across Telluraves:
*   ***Tsidiiyazhi abini*** (Paleocene, ~62.2 Ma; *Ksepka et al., 2017*, cited in *Chen et al., 2026* and *Kimball et al., 2019*): Positioned as a stem mousebird (*Coliiformes* / Coraciimorphae), *Tsidiiyazhi* is the oldest well-corroborated telluravian fossil yet discovered. It establishes a hard minimum divergence age of 62.221 Ma for the Coliiformes–Cavitaves split, proving that the diversification of Afroaves and Telluraves was underway immediately following the K-Pg extinction.
*   ***Eozygodactylus americanus*** (Early Eocene, ~51.81 Ma; *Kimball et al., 2019*): Positioned as a stem passeriform / zygodactylid, *Eozygodactylus* anchors the minimum divergence age for Psittacopasseres (Psittaciformes–Passeriformes split), providing a mandatory calibration bound for Eufalconimorphae and Australaves.

| Geological Epoch / Age | Fossil Taxon / Family | Taxonomic Affinity | Evolutionary Significance |
| :--- | :--- | :--- | :--- |
| **Early Paleocene** (~62.2 Ma) | *Tsidiiyazhi abini* (NMMNH P-54128) | Stem Coliiformes (Afroaves / Coraciimorphae) | Oldest well-corroborated Telluravian; anchors minimum landbird divergence age (*Ksepka et al., 2017*). |
| **Late Paleocene** (~58.5 Ma) | *Ogygoptynx wetmorei* | Stem Strigiformes (Afroaves) | Oldest documented representative of Afroaves and stem owl lineage (*Kimball et al., 2019*). |
| **Early Eocene** (~51.8 Ma) | *Eozygodactylus americanus* (USNM 299821) | Stem Passeriformes / Zygodactylidae | Anchors minimum divergence bound for Psittacopasseres and Australaves (*Kimball et al., 2019*). |
| **Eocene – Pleistocene** | Phorusrhacidae ("Terror Birds") | Cariamiformes (Australaves) | Giant flightless apex predators demonstrating extreme body-size diversification (*Alvarenga & Höfling, 2003*). |
| **Oligocene – Pleistocene** | Teratornithidae | Accipitriformes / Cathartiformes | Giant flying predatory/scavenging birds exhibiting upper limits of body mass in volant raptors (*Chatterjee et al., 2007*). |

---

##### 5. Divergence Timing and Evolutionary Timeline

Establishing the temporal framework for core landbird diversification relies on integrating molecular clock models with stratigraphically constrained fossil calibrations (*Kimball et al., 2019*; *Brocklehurst & Field, 2024*).

```
  Cretaceous                     K-Pg Mass Extinction (~66 Ma)     Paleocene / Eocene
  =======================================|==============================================>
  Crown Neornithes Origin (~100 Ma)      |  Explosive Telluraves Radiation (~66-60 Ma)
  [Brocklehurst & Field, 2024]           |  • Tsidiiyazhi (62.2 Ma)
                                         |  • Ogygoptynx (58.5 Ma)
                                         |  • Eozygodactylus (51.8 Ma)
```

###### The Post-K-Pg Radiation Model
Combined-evidence tip-dating analyses utilizing Bayes Factors indicate that crown birds (*Neornithes*) originated between the Early and Late Cretaceous (~100 Ma) (*Brocklehurst & Field, 2024*). However, both fossil and phylogenomic data demonstrate that major inter-ordinal divergences within Neoaves and Telluraves spanned or immediately postdated the Cretaceous-Paleogene (K-Pg) mass extinction event (~66 Ma) (*Brocklehurst & Field, 2024*; *Stiller et al., 2024*; *Kimball et al., 2019*). Ecological vacuums created by the end-Cretaceous extinction catalyzed rapid diversification across landbird lineages within a brief window of a few million years (*Stiller et al., 2024*).

###### Computational Mechanics of Alignment Data Effects in Timetrees
Divergence time estimation in landbirds exhibits extreme sensitivity to data matrix completeness and rate-smoothing algorithms (*Kimball et al., 2019*). In penalized likelihood analyses (e.g., `treePL`), sparse data alignments containing missing sequence data introduce artificial rate heterogeneity across matrix partitions. Smoothing algorithms penalize these localized rate variations by smoothing branch lengths globally, which artificially pushes deep node divergences deep into the Cretaceous (>70 Ma for crown Neoaves) (*Kimball et al., 2019*). Conversely, when matrices are filtered for high sequence completeness (such as mitochondrial alignments restricted to <25% missing data), `treePL` rate smoothing yields timetrees with deep Neoavian and Telluravian divergences clustered tightly around the K-Pg boundary (~66 Ma) (*Kimball et al., 2019*).

###### Key Numerical Divergence Estimates
*   **Crown Neornithes Origin:** ~100 Ma (Mid-to-Late Cretaceous) (*Brocklehurst & Field, 2024*).
*   **Telluraves Radiation:** ~66.0 – 62.0 Ma (Early Paleocene) (*Ksepka et al., 2017*; *Stiller et al., 2024*).
*   **Stem Coliiformes / Cavitaves Split (*Tsidiiyazhi*):** Hard minimum of **62.221 Ma**; soft maximum of **66.5 Ma** (*Kimball et al., 2019*).
*   **Afroaves / Stem Owl Divergence (*Ogygoptynx*):** Hard minimum of **58.5 Ma** (*Kimball et al., 2019*).
*   **Psittacopasseres Divergence (*Eozygodactylus*):** Hard minimum of **51.81 Ma**; soft maximum of **66.5 Ma** (*Kimball et al., 2019*).

---

##### 6. The Raptorial-Ancestor Hypothesis: Arguments and Evidence

The polyphyletic distribution of raptorial orders across Afroaves (*Accipitriformes*, *Strigiformes*) and Australaves (*Cariamiformes*, *Falconiformes*) has generated opposing hypotheses regarding trait evolution across core landbirds (*Jarvis et al., 2014*; *Prum et al., 2015*; *McClure et al., 2019*; *Chen et al., 2026*).

| Arguments & Evidence **FOR** Ancestral Raptoriality | Arguments & Evidence **AGAINST** Ancestral Raptoriality |
| :--- | :--- |
| **Basal Clade Topology:** Early-diverging branches in both Afroaves (*Accipitriformes*, *Strigiformes*) and Australaves (*Cariamiformes*, *Falconiformes*) are raptorial or predatory (*Jarvis et al., 2014*; *Prum et al., 2015*; *McClure et al., 2019*). | **Intercalated Non-Raptorial Clades:** Phylogenomic trees consistently recover highly speciose, non-raptorial lineages (*Coraciimorphae* in Afroaves; *Psittacopasseres* in Australaves) intercalated between raptorial orders (*Chen et al., 2026*). |
| **Complex Anatomical Suite:** Hooked bills, curved sharp talons, forward-facing eyes, and carnivorous feeding mechanics across disparately placed orders suggest a single ancestral origin rather than four independent acquisitions. | **Requirement for Multiple Losses:** An ancestral raptorial state necessitates multiple independent secondary losses of raptorial adaptations across songbirds, parrots, woodpeckers, and mousebirds (*Chen et al., 2026*). |
| **Macroevolutionary Hemiplasy Mechanism:** Under rapid post-K-Pg speciation, ancestral polymorphic trait variants (e.g., raptorial bill/claw loci) can be sorted differential across lineages (*Somogyi et al., 2026*), creating a scenario where non-raptorial groups retain alternative alleles while raptorial groups inherit ancestral predatory traits without true convergent homoplasy. | **Functional Morphological Convergence:** Mechanical constraints on vertebrate-eating habits exert strong selection for hooked beaks and raptorial talons, driving homoplasy across distinct lineages independently adapted to carnivory (*Somogyi et al., 2026*). |

---

##### 7. Unresolved Nodes, Methodological Conflicts, and Underlying Causes

Despite genome-scale datasets, deep nodes within Neoaves and Telluraves remain recalcitrant. Modern systematists attribute this persistent phylogenetic discordance to distinct biological phenomena and analytical artifacts (*Stiller et al., 2024*; *Mirarab et al., 2024*; *Somogyi et al., 2026*; *Kimball et al., 2019*).

###### Biological Causes: Explosive Radiations and Incomplete Lineage Sorting
The explosive diversification of Neoaves near the K-Pg boundary produced extremely compressed deep internodes. During rapid speciation, ancestral genetic polymorphism fails to fix prior to subsequent lineage splitting, resulting in widespread **incomplete lineage sorting (ILS)** (*Stiller et al., 2024*). Under severe ILS, individual gene trees frequently conflict with one another and with the true species tree topology.

###### Hemiplasy and Morphological Character Evolution
Apparent morphological convergence among landbirds is frequently driven by **hemiplasy**—the phenomenon where a polymorphic genetic trait variant is sorted across rapid speciation events, producing a morphological character distribution that conflicts with the species tree (*Somogyi, Lee, Chen, & Field, 2026*). In Telluraves, hemiplasy provides a precise molecular mechanism explaining why predatory traits (e.g., hooked bills, curved talons) may appear as convergent homoplasies or ancestral retained traits across Afroaves and Australaves when they actually reflect polymorphic ancestral variants sorted across rapid post-K-Pg speciation windows.

###### Genomic Anomalies: The Chromosome 4 Non-Recombining Block
Phylogenomic trees can be systematically hijacked by structural genomic anomalies. A **21-Mb region on chicken chromosome 4** exhibits complete suppression of recombination dating back ~65 Ma (*Mirarab et al., 2024*). Because recombination was blocked across millions of base pairs, this multi-megabase block inherited a single, strongly supported evolutionary history that conflicts with the true species tree.

```
   True Species Tree History             Chromosome 4 Block (21 Mb)
    (Multispecies Coalescent)              (Suppressed Recombination)
           |                                     |
    +------+------+                       +------+------+
    |             |                       |             |
Lineage A     Lineage B               Lineage A     Lineage C  <-- (Distorts Concatenated ML)
```

###### Concatenated ML vs. Multispecies Coalescent (MSC) Frameworks
The 21-Mb chromosomal anomaly exposes a critical methodological divide:
*   **Concatenated Maximum Likelihood (ML):** Concatenating alignment partitions pools all genomic loci into a single matrix. Highly informative, non-recombining structural blocks (like the 21-Mb chromosome 4 region) contribute thousands of concordant site patterns, overwhelming true species-tree signal in concatenated ML analyses under low taxon sampling (*Mirarab et al., 2024*).
*   **Multispecies Coalescent (MSC) Frameworks:** Summary species-tree approaches (e.g., ASTRAL) treat individual non-coding partitions (UCEs, introns) as independent gene trees. MSC models account for ILS and isolate structural anomalies to single locus inputs, preventing localized genomic distortions from hijacking the species tree topology (*Stiller et al., 2024*; *Mirarab et al., 2024*).

---

##### 8. Key Supporting Quotes and Named Source Attributions

The following table compiles verbatim quotes and precise claims from the literature, adhering strictly to mandatory attribution rules for citing papers listed on ResearchGate source pages.

| Topic / Issue | Verbatim Quote / Precise Claim | Exact Source & Attribution Citation |
| :--- | :--- | :--- |
| **Polyphyly of Raptors** | *"Phylogenomic analyses consistently find that these groups are not each other's closest relatives within Telluraves, which may suggest that the raptorial lifestyle evolved multiple times."* | *Chen, van Tuinen, & Field (2026)*, p. 8 |
| **Raptorial Ancestor Hypothesis** | *"...the possibility that raptoriality is ancestral to Telluraves has also been proposed (Jarvis et al., 2014; Prum et al., 2015; McClure et al., 2019)."* | *Chen, van Tuinen, & Field (2026)*, p. 8 |
| **Post-K-Pg Radiation Dynamics** | *"The tree confirms that Neoaves experienced rapid radiation at or near the Cretaceous–Palaeogene boundary... We discovered sharp increases in effective population size, substitution rates and relative brain size following the Cretaceous–Palaeogene extinction event..."* | **This claim originates from a citing paper quoted on the ResearchGate source page:** *Stiller et al. (2024)* |
| **Chromosome 4 Genomic Anomaly** | *"We found that a 21-Mb region in avian genomes, mapped to chicken chromosome 4, shows an extremely strong and discordance-free signal for a history different from that of the inferred species tree... We show that the presence of this region has misled previous phylogenomic efforts..."* | **This claim originates from a citing paper quoted on the ResearchGate source page:** *Mirarab et al. (2024)* |
| **Morphological Hemiplasy** | *"Hemiplasy helps explain high rates of apparent morphological convergence in neoavian birds..."* | **This claim originates from a citing paper quoted on the ResearchGate source page:** *Somogyi, Lee, Chen, & Field (2026)* |
| **Timing of Crown Neornithes** | *"The scenario best supported by our analyses is one where Neornithes originated between the Early and Late Cretaceous (~ 100 MYA), while numerous divergences within major neoavian clades either spans or postdates the K–Pg."* | **This claim originates from a citing paper quoted on the ResearchGate source page:** *Brocklehurst & Field (2024)* |
| **Oldest Telluravian Fossil** | *"...a putative member of the mousebird lineage, Tsidiiyazhi from the Paleocene of North America, is the oldest well-corroborated example of a telluravian yet discovered (Ksepka et al., 2017)."* | *Chen, van Tuinen, & Field (2026)*, p. 8; *Kimball et al. (2019)*, App. A.17 |
| **Sensitivity to Missing Data in Timetrees** | *"treePL analysis of the reduced taxon set tree (which had limited missing data) resulted in a timetree that corroborated several recent studies... showing that Neoaves underwent an explosive radiation close to the K-Pg boundary... use of the more taxon-complete sequence matrix that included more missing data resulted in an older estimated radiation..."* | *Kimball et al. (2019)*, p. 6 & 8 |

---

##### 9. Open Questions and Priority Research Directives

To resolve remaining uncertainties in core landbird systematics and raptor evolutionary dynamics, strategic future research should prioritize the following directives:

1.  **Deployment of Coalescent-Based Genomic Filtering:** Execute whole-genome multi-species coalescent (MSC) analyses across all landbird families while explicitly masking localized structural anomalies—specifically the 21-Mb non-recombining region on chromosome 4—to stabilize the deep backbone topologies of Eufalconimorphae and Afroaves.
2.  **Quantitative Separation of Hemiplasy from Homoplasy:** Integrate high-density morphological datasets with locus-specific gene trees under site-heterogeneous models to quantitatively differentiate predatory anatomical traits produced by true convergent homoplasy from those inherited via hemiplasy during post-K-Pg radiation windows.
3.  **Targeted Paleocene Stratigraphic Field Sampling:** Intensify paleontological collection efforts in early Paleocene terrestrial deposits (~66–60 Ma) globally to uncover missing stem fossils of Falconiformes, Accipitriformes, and Strigiformes, bridging the gap between molecular clock divergence estimates and the physical fossil record.
4.  **Comparative Functional Genomics of Predatory Suites:** Perform comparative genomic and developmental analyses of the gene networks governing bill keratinization, digital tendon mechanics, and retinal foveal development across falcons, owls, hawks, and seriemas to test the raptorial-ancestor hypothesis at the molecular developmental level.
---

## Auditoria

**Método.** Para cada resposta, comparei cada afirmação com o **trecho que o próprio NotebookLM devolveu** para a citação (`references[].cited_text`). Quando foi preciso, li também o texto indexado da fonte (`source_get_content`) ou os arquivos `.txt` que eu mesmo subi. Em cada item marquei:

- (a) afirmações cuja fonte citada **não confere** com o trecho, ou que saíram de algo que não é a obra (lista de referências, trechos de artigos citantes);
- (b) afirmações que parecem **inferência do NotebookLM além das fontes**.

O que não foi marcado confere com os trechos. Nenhuma citação foi inventada. Quando não consegui conferir algo, isso está dito.

### Problemas que atravessam todas as respostas

1. **As páginas do ResearchGate não são só a obra.** A página de "Resolving the Avian Tree of Life from Top to Bottom" (2019, `47d8810c`) e a de Chen et al. 2026 "Birds, Diversification of" (`d3d5a7ff`) incluem a seção "Citations" (trechos de **artigos que citam** a obra) e, na de Chen, os **resumos das referências**. O NotebookLM atribui esse texto à obra da página. Caso documentado: a reanálise wASTRAL dos dados de Stiller 2024, atribuída ao capítulo de 2019 (Q9 [9]). Os resumos de Mirarab et al. 2024 (região de 21 Mb), Brocklehurst & Field 2024 e Somogyi et al. 2026 (hemiplasia) aparecem nos relatórios do Studio como se fossem fontes do caderno, mas nenhum deles foi importado.
2. **Fontes só com resumo:** Fuchs et al. 2015 (Falconidae; só a prévia do ScienceDirect), Reddy et al. 2017, "Tip-dating and the Origin of Telluraves" e McClure et al. 2019 (só a página de citação). As respostas sobre Falconidae (Q5, Q8) dependem só do resumo de Fuchs.
3. **Lacunas importantes no caderno:** Mindell et al. 2018 (capítulo pago), Wink et al. 2009, Degrange et al. 2015, Mayr 2011 (*Messelastur*) e Ericson et al. 2006 não entraram (ver tabela de falhas). A tese de Oxford que testa a raptorialidade ancestral também não. Isso empobrece Q5 (corujas: só Salter 2020), Q6 (forarracídeos: só Degrange 2010) e Q7 (contra-evidência).
4. **Respostas com mais inferência própria:** Q7 (ancestral raptorial) e Q10 (caracteres convergentes). As que mais aderem às fontes: Q3, Q4 (com dois erros de fato) e Q8.

### Por pergunta

**Q1**

- *Erro de atribuição.* "Concatenated ML analysis placed Accipitrimorphae alone as sister to remaining Afroaves with **32% bootstrap** [4, 12]": o 32 do trecho [12] é o bootstrap de **Strisores** como irmão de Telluraves, não de Accipitriformes. Stiller et al. 2024 diz só que a análise concatenada pôs Accipitriformes como irmão do resto de Afroaves; o valor de suporte desse nó não aparece no trecho.
- *Rótulo errado.* "Wu et al. (2024) (Coding Sequences / 5,756 CDS)": a árvore ótima de Wu et al. sai de uma análise de sensibilidade sobre quatro classes de DNA (25 460 locos), não só de CDS. "5,756 CDS" é uma das classes. Chamar a topologia de "CDS" é inferência do NotebookLM.
- *Sem fonte.* "Prum… **41% ML bootstrap** for broader nodes": nenhum trecho citado traz 41%. Na tabela [2], Prum 2015 tem Coraciimorphae = 1.0 PP e Afroaves = "Not recovered".
- *Conflação.* "Suh et al. (2015)… place Strigiformes as sister to **Cavitaves/Coraciimorphae**": as árvores de retroposons põem as corujas como irmãs de **Cavitaves**, com Coliiformes mais profundo, e então **Coraciimorphae não é recuperado** (tabela [2], Suh 2015: "Not recovered"). "MPRE" e "conflict-free" não aparecem nos trechos; a tabela diz só "16 retroposons" para Afroaves.
- *Exagero leve.* Retroelementos de 2022 [18, 19]: as contribuições de 55% e 45% vêm de táxons terminais (*Leptosomus discolor*, de Cavitaves, e *Haliaeetus leucocephalus*) num teste com quatro espécies. "Stem of Coraciimorphae / stem of Accipitrimorphae" é extrapolação do NotebookLM, e o próprio trabalho chama o segundo teste 4-LIN de ambíguo.
- *Parcial.* Kuhl 2021 com "Accipitrimorphae sister to (Strigiformes + Coraciimorphae)" está citado como [2, 8], mas quem sustenta essa ordem interna é o trecho [15], da mesma fonte de [2]. As fontes conferem; só o número da citação está deslocado.

**Q2**

- *Contradito pela própria fonte citada.* "First recovered … (e.g., Hackett et al. 2008)": o trecho [1] diz que o clado papagaio+passeriforme "was first found in analyses of Ericson et al. (2006)", e o [24] atribui Australaves também a Ericson et al. 2006.
- *Sem fonte.* "UCEs … strongly resolves Eufalconimorphae **and Australaves**": na tabela [10], McCormack 2013 (UCEs) tem Eufalconimorphae >0.95 PP/>70%, mas Australaves está em branco.
- *Sem fonte.* "maximal statistical support (100% …) for Eufalconimorphae" em Stiller 2024 e Wu 2024: nenhum trecho traz valor de suporte desses dois estudos. Jarvis 2014 (TENT) é análise concatenada, não "coalescent". Só Jarvis 100% e Kuhl 100% aparecem na tabela.
- *Falso segundo a tabela citada.* "near-100% support across all major phylogenomic sequence analyses" para Australaves: na mesma tabela [10], Hackett tem 64%, Burleigh 72% e Liu 2018 <70%.
- *Inferência.* "CR1 LINEs **and LTR retrotransposons**": os trechos não falam em LTR. Em "falcons near … mousebirds" [20–22], quem cai junto dos colíiformes nas árvores CODON12/AAS de Kuhl são passeriformes+papagaios [21], não falcões.
- *Omissão relevante, sem erro.* O trecho [20] registra a hipótese alternativa: éxons codificantes não particionados (Braun & Kimball 2021) e o mosaico dos masillaraptorídeos sustentam **Cariamiformes + Falconiformes** como grupos-irmãos. O [11] lembra que Ericson 2012 desenhou Cariamidae como irmão de Falconidae. A resposta apresenta Cariamiformes como basal sem ressalva.

**Q3**

- *Correspondência boa.* As quatro afirmações de topologia batem com os trechos de Catanach 2024 [4–8] e da tese de Mather 2021 [7, 9].
- *Omissão que muda a leitura.* As idades (61,3 / 60,9 / 50,8 Ma) são de Catanach 2024 e aparecem como se fossem consenso. O trecho [8], que a própria resposta cita, diz que Mindell et al. 2018 estimam ~45 e ~38 Ma e Prum et al. 2015 ~40 e ~28 Ma para os mesmos nós, e atribui a diferença à calibração (prior no nó Cathartiformes–Accipitriformes).
- *Generalização.* "Morphological assessments consistently recover…": a fonte morfológica citada [1] (Chen 2025) otimiza sinapomorfias **sobre topologias moleculares fixas** (molecular scaffold). Não é uma análise morfológica independente que recupere o clado. O [5] cita Griffiths 1994, que não está entre as fontes.

**Q4**

- *Contradito pela fonte citada.* "taxonomic authorities (including … **AOS/NACC**) have split former *Accipiter*": o trecho [22] diz o contrário, "The AOS Checklist has not yet considered these changes… still recognizing the expanded view of the genus". O documento 2025-C é uma **proposta** submetida a voto [39, 40], e o resultado não está nas fontes. Quem adotou a divisão foram WGAC (2024), Clements (2024) e IOC (2025).
- *Topologia invertida.* "Because *Circus* is sister to the *Astur* + *Megatriorchis* clade": os trechos [37, 38] dizem que **Astur é irmão de (Megatriorchis + Circus)**. Circus é irmão de Megatriorchis.
- *Posição errada de Harpaginae.* "Lophospizinae & Harpaginae: deep lineages sitting basal to the core hawks/buzzards": pelo [16], em Catanach 2024 Harpaginae (*Harpagus* + *Microspizias*) é **irmão de Buteoninae**. Só em Mindell et al. 2018 aparece como irmão de Buteoninae + Accipitrinae. Basal ao conjunto é só *Lophospiza* [15, 28].
- *Omissão.* *Polyboroides* fica fora de Gypaetinae + Perninae e torna Gypaetinae não monofilético [6]. A resposta o lista dentro de Gypaetinae sem ressalva. *Circaetus* também não é monofilético (*Dryotriorchis* aninhado) [3].
- *Inferência leve.* "~27 species of mostly Indomalayan, Australasian, and African hawks" (*Tachyspiza*): a distribuição não aparece nos trechos [31, 32].

**Q5**

- *Fonte fraca.* Para Falconidae o NotebookLM só tem o **resumo e a prévia** de Fuchs et al. 2015 no ScienceDirect, não o texto completo. A resposta lista as três subfamílias, mas não dá a ordem de ramificação entre elas, porque os trechos não a trazem.
- *Citação de lista de referências.* "*Micrastur*… underwent radiation during the Pliocene [1, 2, 4]": o [4] é o título de um artigo de 2011 na lista de referências de Fuchs 2015. A afirmação é plausível, mas não tem trecho de resultado que a sustente.
- *Sem fonte.* "Extant owls underwent a basal bifurcation **during the Paleogene/Oligocene**": os trechos datam no Oligoceno só a separação **Tytoninae/Phodilinae** dentro de Tytonidae [8]. A idade da separação Tytonidae/Strigidae não aparece em nenhuma fonte citada.
- *Generalização.* Composição das tribos de Clade A ("Surniini: *Glaucidium, Athene, Aegolius, Surnia, Micrathene, Xenoglaux, Margarobyas*"): os trechos [14–19] não listam esses membros. A legenda da figura de Salter traz **Aegolini** como tribo à parte. A resposta também omite um ponto que o [16] destaca: na análise coalescente (SVDquartets) *Uroglaux* sai como irmão de A+B, o que torna instável a composição de Clade A.
- *Inferência leve.* "before modern owl perching and grasping tarsometatarsal modifications fully evolved": é a hipótese de Mourer-Chauviré 1994 citada por Mayr et al. 2020 [26]. O mesmo trecho diz que as afinidades dos estrigiformes do Cenozoico inicial "have not yet been subjected to a formal phylogenetic analysis".

**Q6**

- *Contradito pela fonte citada.* "Early Eocene stem owls like *Primoptynx* and protostrigids (*Minerva*, ***Eostrix***) retained… raptorial talon proportions": o trecho [23] diz "The pedal phalanges are unknown for… the species of *Eostrix*". O mesmo trecho mostra que *Palaeoglaux* **não** tem as garras aumentadas. As garras grandes nos dedos I e II estão documentadas só em *Primoptynx* e *Minerva*, e o artigo trata isso como derivado, não como ancestral. Os processos supraorbitais [24] vêm de crânios do London Clay e da Bridger Formation, não de *Primoptynx*. *Minerva* é do Eoceno inferior tardio ao médio inicial, não do Eoceno inferior.
- *Conflito entre fontes não sinalizado.* Idade de Phorusrhacidae: Degrange et al. 2010 [26] diz "middle Paleocene to lower Pleistocene", e a fonte [25] diz "Eocene–Pleistocene". O "~60 Ma" é conversão do NotebookLM, e nenhuma fonte nomeia o fóssil forarracídeo mais antigo.
- *Interpretação além da fonte.* "remarkable morphological **convergences**" entre masillaraptorídeos e forarracídeos: Mayr & Kitchener [1, 6] deixam três hipóteses abertas (plesiomorfia, elo evolutivo direto, paralelismo) e não falam em convergência. Os próprios autores dizem que "the published character evidence for this placement is still weak".
- *Incerteza atenuada.* Messelasturidae como "stem-group Psittacopasseres or stem Psittaciformes": o [11] chama isso de "mere speculation… far from being well-established". A análise sem restrição não resolve a posição, e o [12] diz que algumas análises "do not even support psittacopasserine affinities".
- *Omissão.* A pergunta pedia o fóssil mais antigo de **cada** linhagem. A resposta não traz Accipitridae, Cathartiformes nem Falconidae coroa, embora o caderno tenha material: o pequeno acipitrídeo do Oligoceno inferior da Polônia (Mayr & Hurum 2020), *Archaehierax* do Oligoceno superior (Mather et al. 2021) e *Stintonornis* / *Antarctoboenus* [4, 7] como possíveis falconiformes de 53,5–49 Ma.

**Q7** (a resposta com mais inferência própria)

- *Lógica invertida em relação às fontes.* "Analyses… frequently resolve a clean split between Afroaves and Australaves, making independent gains of raptoriality within each clade a highly supported alternative [7, 23–25]": o trecho [7] diz que, na divisão Afroaves/Australaves, **cada clado tem duas linhagens raptoriais na base**. É exatamente nessa topologia que Jarvis et al. 2014 [3] apoiam o ancestral predador. Nenhuma fonte diz que a divisão favorece origens independentes, e "highly supported" é do NotebookLM. O [7] conclui que "it seems inappropriate to reject any of these three hypotheses at this time".
- *Argumento sem fonte.* "Opponents argue that it is morphologically and ecologically simpler for raptorial traits to evolve independently… than… to be lost independently across dozens of… families [2, 3]": o [2] diz apenas "may suggest that the raptorial lifestyle evolved multiple times", e o [3] (Jarvis) defende o contrário. O argumento de parcimônia é construção do NotebookLM.
- *Inferência sobre Tsidiiyazhi.* "showing that non-predatory landbird niches were already established at the very beginning of the telluravian radiation": nem Ksepka 2017 [20] nem a fonte [19] usam o fóssil como evidência contra o ancestral raptorial. Ksepka 2017 diz que o foco é a posição de grupo-tronco de Coliiformes. A fonte [21] (tip-dating) põe a origem de Telluraves coroa no **Cretáceo Superior**, então ~62 Ma não é "the very beginning". O próprio [3] diz que os colíiformes do Eoceno têm "raptor-like feet".
- *Generalização.* "considered structurally parsimonious by several authors, requiring only a single loss": o [16] formula a hipótese com condicional ("If… may have been raptorial") e uma ressalva que a resposta **omite**: "a monophyletic Eufalconimorphae is not recovered in all sequence-based analyses (Stiller et al. 2024: extended data, fig. 4)". O [14] chama a hipótese de "mere speculation".
- *Sem fonte.* "specialized flight musculature" como exemplo de convergência [17, 18]: os trechos falam de forma do bico, organização dos dígitos e morfologia externa e esquelética, não de musculatura. "Whole-genome phylogenomic analyses (e.g., … Prum et al. 2015)": Prum usa enriquecimento por hibridização ancorada (259 locos), não genoma inteiro.
- *Juízo do NotebookLM.* "**moderate-to-strong support** for a raptorial ancestor [de Australaves]" e "equally prominent": nenhuma fonte gradua o suporte. São avaliações do modelo.
- *Lacuna do caderno.* O teste ecomorfológico de raptorialidade por tarsometatarso (tese de Oxford, ora.ox.ac.uk, uuid:90912ef0…) concluiu por convergência, mas **não foi importado** (HTTP 403). A herbivoria do cariamiforme *Strigogyps*, citada como contra-argumento em fontes secundárias, também não está em nenhuma fonte do caderno.

**Q8**

- *Método errado.* "**Coalescence-based** analyses by Wu et al. (2024) place the basal divergence of Telluraves at ~86 Ma": o trecho [6] diz que a datação foi feita com **MCMCTree sobre dados CDS** e 20 calibrações fósseis, sobre a árvore ótima. Não é uma datação coalescente.
- *Conflação e número sem fonte.* "3′-UTR analyses… yield basal landbird divergence dates preceding the K–Pg boundary (**~66–80+ Ma**) [8]": Kuhl 2021 fala de "Basal Landbirds" como um grupo informal do próprio artigo, que não é a espinha dorsal de Telluraves (Telluraves corresponde ao que Kuhl chama de Higher Landbirds). O intervalo 66–80+ Ma não aparece no trecho.
- *Generalização.* "Post–K–Pg Radiation… Phylogenomic analyses (such as those based on whole-genome or UCE matrices)": Stiller 2024 [1] usou regiões intergênicas e 187 ocorrências fósseis. O trecho não dá idade para a coroa de Telluraves, só diz que o IC 95% de Telluraves–Elementaves atravessa o limite K–Pg. A faixa "~62–66 Ma" é síntese do NotebookLM.
- *Lacuna na resposta.* A idade da coroa de **Strigiformes** (separação Tytonidae–Strigidae) não aparece. A resposta dá só a idade-tronco (fósseis do Paleoceno, ≥60 Ma) e a radiação de Strigidae no Mioceno. A idade de Accipitrinae em "10–17 Ma" também não está no [13], que cita Aquilinae, Buteoninae e Haliaeetinae.
- *Fontes do caderno não usadas.* Claramunt & Cracraft 2015 (árvore temporal), Kimball 2019 (supertree datada) e o capítulo do TimeTree of Life (Brown 2009, cap. 63) têm datações relevantes e não foram consultados nesta resposta.
- *Conferem com os trechos.* Accipitridae: 45,3 Ma (Catanach) contra 33,7–34 Ma (Mindell 2018; Nagy & Tökölyi). Falconidae: máximo de 34,2 Ma, Falconinae ~16 Ma, *Falco* ~7,5 Ma. Telluraves em Wu: ~86 Ma (80,4–91,5).

**Q9**

- *Contaminação de fonte (vale para todo o caderno).* O trecho [9] ("a wASTRAL reanalysis of the Stiller et al. (2024) gene trees fails to recover [Afroaves]. Instead, it puts Accipitriformes sister to all landbirds") está atribuído a "Resolving the Avian Tree of Life from Top to Bottom" (Braun, Cracraft & Houde, **2019**), que não poderia citar um trabalho de 2024. Conferi o texto indexado (`source_get_content`): a página do ResearchGate inclui a seção "Citations (66)", com trechos de **artigos que citam o capítulo**. A afirmação vem de um artigo posterior não identificado. **Toda citação ao `47d8810c` (e, por precaução, às outras páginas do ResearchGate: `d3d5a7ff`, `0eef4b43`, `79ba75e1`) pode ter esse problema.**
- *Citação de bibliografia.* "Ancient introgressive hybridization… [7, 44–46]": o [44] é só a lista de referências do capítulo "Birds, Diversification of" (McClure 2019; Mirarab et al. 2024, "A region of suppressed recombination misleads neoavian phylogenomics"). Não é evidência. Mirarab et al. 2024 não está no caderno como fonte própria. O resumo dele (a região de 21 Mb no cromossomo 4 da galinha, com recombinação suprimida) aparece apenas embutido na lista de referências da página ResearchGate de Chen 2026. É esse trabalho que propõe a recombinação suprimida como causa de conflito, e a resposta Q9 não o usa.
- *Sem fonte / errado.* Corujas como irmãs de Cavitaves "in retroposon (RI) **and 3′-UTR** trees [6, 11, 12]": o [11] é uma árvore de Chen 2025 **restrita** ao consenso de Braun & Kimball 2021, e nenhum trecho mostra Kuhl 2021 (3′-UTR) com essa topologia. Hieraves "in recent whole-genome **and morphological matrix** analyses [7, 13, 14]": Chen 2025 [13] não recupera Hieraves, só otimiza caracteres sobre as topologias de Stiller e Wu. "whole-genome, exon… (e.g., Prum et al. 2015…)": Prum não é genoma inteiro.
- *Exagero leve.* "Anomaly zone": Jarvis 2014 [31] diz "It is possible that some of the deep branches… are in the anomaly zone, although the gene tree support is not high enough to confidently test this idea". A resposta apresenta isso como causa. "40%–100% of loci": no [28] o número se refere a **marcadores de retroposons** nos ramos cujo intervalo cruza o limite K–Pg, não a "loci" em geral.
- *Conferem.* ILS e ramos curtos (Suh 2015 [27–30, 47]), efeito do tipo de dado (Braun & Kimball [39]; Reddy 2017), erro de reconstrução de árvores de genes ([16, 43]), assimetria KKSC/4-LIN ([46]) e efeito da amostragem de passeriformes ≥ 138 (Stiller [7, 48]).

**Q10**

- *Omissões centrais da própria fonte citada.* Livezey & Zusi 2007 [6] recuperaram **Pandionidae + Falconidae** como grupos-irmãos dentro de Falconiformes, um resultado morfológico que contraria todas as filogenômicas e que a resposta não menciona. O [5] diz que o clado aves de rapina (Falconiformes + Strigiformes) foi o **primeiro de uma série de grupos-irmãos sucessivos do resto de Neornithes**, uma posição quase basal que a resposta também omite.
- *Convergência mal atribuída.* "Phalangeal fusion… occurred independently **across distinct predatory groups** [26, 27]": os trechos tratam da fusão do dígito II **dentro de Accipitridae** (milvíneos e haliaeetíneos contra *Ictinia*/*Busarellus*), não entre grupos raptoriais distintos.
- *Contradito pela fonte.* "Enlarged claws… with **prominent flexor tubercles** [2, 16, 21, 24]": o [24] diz que nos messelasturídeos os tubérculos flexores são **menos** desenvolvidos que na maioria dos raptoriais e que proeminente é o extensor. O [21] mostra que corujas e falcões têm garras de tamanho igual e só os acipitrídeos têm as dos dígitos I e II hipertrofiadas. Isso é diferença entre linhagens, não convergência que tenha unido os grupos.
- *Inferência.* "Shortened middle phalanges on Digit IV… evolved repeatedly [25]": o [25] apresenta isso como caráter **falconiforme** usado para classificar *Masillaraptor*. Mayr & Kitchener (Q2 [20]) dizem que Cariamiformes e Falconiformes compartilham o caráter, o que pode ser sinapomorfia de clado e não convergência. As funções ("to shield the eyes or anchor jaw musculature", "to reinforce foot strength") não aparecem nos trechos.
- *Sem fonte.* "Pectoral girdle… features associated with powerful flapping or precise aerial maneuvering… repeatedly misread as shared ancestral traits [22, 24, 28–30]": os trechos descrevem caracteres (forame do supracoracoide, processo acromial, forame procoracoide) e sua homoplasia, mas nenhum diz que foram interpretados erradamente como ancestrais ou que são adaptações de voo.
- *Conferem.* "Falconimorphae" [9]; Cathartidae basal em Falconiformes morfológico [6, 10]; a ligação Cathartidae–cegonhas não foi recuperada e o estudo de mtDNA que a sustentava foi retratado [10]; convergência entre urubus e abutres do Velho Mundo [10, 31]; homoplasia no bico e nos dígitos dos pés na árvore morfológica de Ksepka [2]. Fonte relevante pouco usada: Chen et al. 2025 atribui o conflito morfologia × genômica a "functional convergence and rapid cladogenesis".

### Relatórios do Studio

**Briefing Doc (`6c992653`)**

- *Fora do tema pedido.* A maior parte trata de origem nos dinossauros, aves mesozoicas, ratitas e aves aquáticas, a partir de Chen et al. 2026. Sobre a posição das aves de rapina há só uma linha de tabela (Telluraves = Afroaves + Australaves).
- *Citações conferidas.* A frase de "Chen et al. (2026)" ("birds are living dinosaurs") está no texto indexado do capítulo, cujos autores na página do ResearchGate são Albert Chen, Marcel van Tuinen et al. A de Stiller et al. 2024 ("sharp increases in effective population size…") está no texto completo de Stiller.
- *Fonte que não é do caderno.* As citações de **Brocklehurst & Field (2024)** e **Mirarab et al. (2024)** são verbatim, mas vêm dos **resumos de referências embutidos na página ResearchGate** de Chen 2026. Nenhum dos dois artigos foi importado. Há ainda um resíduo de outra língua no texto ("来自 *Brocklehurst & Field (2024)*").
- *Não conferido.* A citação de Kimball et al. 2019 (">400 years… Jarvis et al. data") é plausível, mas não abri o trecho.

**Relatório focado (`bee72a25`)**

- *Falso para o caderno.* "*Masillaraptoridae* & *Messelasturidae*… do not appear in the primary ground-truth source text" e "*Archaehierax*… unmentioned": o caderno tem Mayr & Kitchener 2022, Mayr 2021, o artigo de Psittacomimidae, Mather et al. 2021 e a tese de Mather. O gerador usou só um subconjunto das fontes (sobretudo Chen 2026 e Kimball 2019). Os fósseis estão tratados nas respostas Q6 e Q7.
- *Atribuição errada.* Diz que a citação de Stiller et al. 2024 "originates from a citing paper quoted on the ResearchGate source page", mas Stiller 2024 está no caderno com texto completo (`e9c5ae30`). Somogyi et al. 2026, Mirarab et al. 2024 e Brocklehurst & Field 2024 vêm de fato só das páginas do ResearchGate, como o próprio relatório avisa.
- *Erros de fato em relação a outras fontes do caderno:*
  - "Sagittariidae… earliest-diverging extant lineage within total-clade Accipitriformes" e a linha "Cathartiformes… sister group: Accipitridae / Teratornithidae" contrariam Catanach 2024 e Mather 2021: Cathartiformes é irmão de Accipitriformes, depois vêm Sagittariidae, depois Pandionidae + Accipitridae.
  - Falconidae "subfamilies Polyborinae and Falconinae" omite Herpetotherinae (Fuchs 2015).
  - "*Ogygoptynx*… oldest documented representative of Afroaves" (58,5 Ma, número de calibração de Kimball 2019): Mayr et al. 2020 dá 60–59 Ma para *Ogygoptynx* e ~61 Ma para *Berruornis halbedeli*.
  - "*Archaehierax*: The Australian Eocene/Oligocene accipitrid": é do Oligoceno superior (26–24 Ma) segundo Mather et al. 2021.
- *Inferência sem fonte.* "Complex Anatomical Suite… suggest a single ancestral origin" e o parágrafo sobre hemiplasia de "raptorial bill/claw loci": Somogyi et al. 2026 só existe no caderno como título e resumo na página do ResearchGate, e a aplicação a bico e garras de rapinantes é extrapolação. Os "Open Questions" são sugestões do modelo.
- *Confere.* As citações de Chen 2026 (p. 8) sobre a polifilia dos rapinantes e a proposta de ancestral raptorial batem com o trecho [1]/[2] de Q7. *Tsidiiyazhi* (62,221 Ma) e a sensibilidade das datações do treePL a dados faltantes (Kimball 2019) são plausíveis, mas não conferi o texto de Kimball.
