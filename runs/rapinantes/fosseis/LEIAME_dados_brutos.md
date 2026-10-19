# Matrizes morfológicas e árvore de referência — braço paleontológico (rapinantes)

Este LEIAME descreve `data/raw/rapinantes/morfologia/`. Ele foi escrito aqui
porque o hook `guard_data.py` bloqueia escrita de arquivo novo por ferramenta
dentro de `data/raw/`; para tê-lo junto dos dados, copie-o à mão para
`data/raw/rapinantes/morfologia/LEIAME.md`.

Baixado em 2026-10-17 (hora de Brasília; 2026-10-18 UTC). Arquivos brutos:
não editar nem sobrescrever. Checksums em `log/sha256_brutos.txt` (nesta
pasta). Todos os DOIs abaixo foram conferidos na Crossref ou na DataCite em
2026-10-17; nenhum foi escrito de memória.

| Pasta / arquivo | O que é | Origem (URL exata do download) | Artigo (DOI conferido) | Licença |
|---|---|---|---|---|
| `lloyd/Ksepka_etal_2019a.nex` | 62 táxons × 146 caracteres; Telluraves atuais (Falco, Buteo, Sagittarius, Cathartes, Tyto, Bubo, Cariama, papagaios, passeriformes, Coraciimorphae) + 28 fósseis paleógenos (halcionitídeos, zigodactilídeos, *Eofringillirostrum*...) | http://www.graemetlloyd.com/nexus/Ksepka_etal_2019a.nex | Ksepka, Grande & Mayr 2019, Curr Biol 29:657 — 10.1016/j.cub.2018.12.040 | coleção de G. T. Lloyd (recodificação da matriz publicada) |
| `lloyd/Ksepka_etal_2017a.nex` | 48 × 111 (*Tsidiiyazhi*); baixado, não analisado (contida na de 2019) | http://www.graemetlloyd.com/nexus/Ksepka_etal_2017a.nex | 10.1073/pnas.1700188114 | idem |
| `lloyd/Mayr_2011c.nex` | 33 terminais de nível família × 100; Messelasturidae, Halcyornithidae, Sandcoleidae + todas as famílias de rapinantes | http://www.graemetlloyd.com/nexus/Mayr_2011c.nex | Mayr 2011, J Syst Palaeontol 9:159 — 10.1080/14772019.2010.505252 | idem |
| `lloyd/Mayr_2005a.nex` | 23 × 110; *Messelastur*, *Tynskya*, Pseudasturidae, corujas, falconiformes s.l., Cariamidae | http://www.graemetlloyd.com/nexus/Mayr_2005a.nex | Mayr 2005, J Vertebr Paleontol 25:635 — 10.1671/0272-4634(2005)025[0635:tpoapp]2.0.co;2 | idem |
| `lloyd/Mayr_2016a.nex` | 13 × 40; *Bathornis*, *Paracrax*, Idiornithidae, fororracídeos | http://www.graemetlloyd.com/nexus/Mayr_2016a.nex | Mayr 2016, J Paleontol 90:357 — 10.1017/jpa.2016.45 | idem |
| `lloyd/Mayr_et_Clarke_2003a.nex` | 46 × 148; Neornithes em nível família | http://www.graemetlloyd.com/nexus/Mayr_et_Clarke_2003a.nex | Mayr & Clarke 2003, Cladistics 19:527 — 10.1016/j.cladistics.2003.10.003 | idem |
| `lloyd/Livezey_et_Zusi_2007a.nex` | 188 × 2954 (usados os 154 Neornithes) | http://www.graemetlloyd.com/nexus/Livezey_et_Zusi_2007a.nex | Livezey & Zusi 2007, Zool J Linn Soc 149:1 — 10.1111/j.1096-3642.2006.00293.x | idem |
| `lloyd/Mayr_etal_2003a.nex` | 29 × 89 (“higher land birds”); baixado, não analisado | http://www.graemetlloyd.com/nexus/Mayr_etal_2003a.nex | 10.1046/j.1439-0469.2003.00230.x | idem |
| `lloyd/matrbird.html`, `matrdino.html` | páginas-índice da coleção no dia do download | http://www.graemetlloyd.com/matrbird.html ; http://www.graemetlloyd.com/matrdino.html | — | — |
| `morphobank_P510/` | a mesma matriz de Livezey & Zusi, com lista de caracteres (conferido: 188 × 2954, 100% das células iguais à versão de Lloyd) | https://www.morphobank.org/services/public/projects/510/matrices/1228/download?format=nexus&notes= e https://www.morphobank.org/services/public/projects/510/matrices/1228/download/characters?notes= | projeto 10.7934/P510 | termos de uso do MorphoBank |
| `morphobank_P3896/` | Degrange et al. 2021: 72 × 207 (lista de caracteres em português); Cathartidae (7), Accipitridae (24), Falconidae (8), *Pandion*, *Sagittarius*, cegonhas, garças, íbis + fósseis *Pampagyps imperator* e *Dryornis pampeanus* | https://www.morphobank.org/services/public/projects/3896/matrices/27138/download?format=nexus&notes= e .../download/characters?notes= | Degrange et al. 2021, Pap Palaeontol — 10.1002/spp2.1361; projeto 10.7934/P3896 | CC0 (projeto marcado publish_cc0) |
| `mayr_kitchener_2022/` | 10 × 42; *Masillaraptor*, *Danielsraptor*, Falconidae, Cariamiformes extintos; descrição dos caracteres (.docx) | https://ndownloader.figshare.com/files/36204378 (nex), https://ndownloader.figshare.com/files/36204381 (docx); registro 10.6084/m9.figshare.20267628.v2 | Mayr & Kitchener 2022, J Vertebr Paleontol 41:e2083515 — 10.1080/02724634.2021.2083515 | CC BY 4.0 |
| `mayr_2026_vz/` | 42 × 105 (tabela .docx; File S1 = caracteres, S2 = matriz); *Messelastur*, *Tynskya*, halcionitídeos, Falconidae, Strigiformes, Psittacopasserae | https://vertebrate-zoology.arphahub.com/article/182505/download/suppl/31/ | Mayr 2026, Vertebr Zool 76:121 — 10.3897/vz.76.e182505 | ODbL |
| `chen_etal_2025/` | 77 × 203, cintura escapular e membro anterior, codificados em boa parte em TC (Bjarnason & Benson 2021; Benito et al. 2022); 75 atuais incl. *Pandion*, *Elanus*, *Coragyps*, *Micrastur*, *Tyto*, *Ninox* | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12409627/supplementaryFiles | Chen et al. 2025, Integr Org Biol 7 — 10.1093/iob/obaf029 | artigo em acesso aberto (PMC) |
| `degrange_etal_2015/` | 19 × 84 (interno a Phorusrhacidae); baixado, não analisado (só Cariamidae como atual) | https://ndownloader.figshare.com/files/1963424 ; registro 10.6084/m9.figshare.1344833.v2 | 10.1080/02724634.2014.912656 | CC BY 4.0 |
| `stiller_2024/63K.tre` | árvore de espécies principal (ASTRAL, 63 430 regiões intergênicas, 363 espécies): topologia do esqueleto molecular | https://erda.ku.dk/archives/341f72708302f1d0c461ad616e783b86/B10K/data_upload/03_species_trees/63K.tre (o README do repositório aponta o mesmo arquivo em https://sid.erda.dk/share_redirect/ENhZODU9YE/03_species_trees/63K.tre) | Stiller et al. 2024, Nature 629 — 10.1038/s41586-024-07323-1; dados 10.17894/ucph.85624f66-c8e5-4b89-8e8a-fe984ca89e4a | ver repositório |
| `stiller_2024/63K_concatenation.tre`, `README_B10K.txt` | árvore concatenada e README do repositório | mesma pasta do ERDA | idem | idem |

## Procurado e não obtido

- **Corujas-tronco** (*Ogygoptynx*, *Berruornis*, Protostrigidae, *Primoptynx*,
  *Sophiornis*, *Ypresiglaux*): nenhuma matriz publicada que as inclua (coleção
  de Lloyd; listagem completa da API pública do MorphoBank; figshare; Europe
  PMC; web).
- ***Danielsavis*** não é rapinante: é Danielsavidae (Galloanseres; MorphoBank
  P5084). O masilarraptorídeo de Walton-on-the-Naze é ***Danielsraptor
  phorusrhacoides*** (`mayr_kitchener_2022/`).
- ***Masillaraptor buchheimi*** (Li, Musser, Aase & Clarke 2025, Contrib Mus
  Paleontol Univ Michigan — 10.7302/27584): PDF no Deep Blue atrás de
  verificação anti-robô; não baixado.
- **Listas de caracteres de Ksepka et al. 2019/2017 e de Mayr 2005/2011, Mayr &
  Clarke 2003**: suplementos atrás de paywall ou anti-robô; as matrizes vieram de
  Lloyd sem descrição de caráter.
- ***Horusornis*, fósseis de Pandionidae**: nenhuma matriz encontrada.
