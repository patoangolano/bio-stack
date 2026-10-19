# Posição filogenética das aves de rapina

Meta-análise e reanálise computacional de onde ficam as aves de rapina na árvore das aves, entre si e em relação às
demais Neoaves: Cathartidae, Accipitriformes (Sagittariidae, Pandionidae, Accipitridae), Strigiformes e Falconiformes.
Quatro frentes de evidência, pesadas pela independência entre si: filogenômica nova, morfologia de fósseis, estado
ancestral do hábito rapinante e reanálise de árvores e quartetos publicados.

O texto principal está em [`reports/Posição filogenética das aves de rapina.md`](reports/). A versão 2 é de 19 de
outubro de 2026. Ela conclui que as linhagens rapinantes estão todas em Telluraves, mas não formam um clado, e
sustenta Australaves (Cariamiformes + Falconiformes + Psittacopasserae) e Accipitrimorphae (Cathartidae +
Accipitriformes). Também mostra que os pés rapinantes surgiram três vezes de forma independente.

A posição das corujas (irmãs de Accipitrimorphae ou de Coraciimorphae) continua sem solução.

## Estado em 19 de outubro de 2026

| Frente | Onde | Estado |
|---|---|---|
| Filogenômica: 300 locos de cópia única, 57 táxons | `runs/filogenia/20261018-*`, `runs/rapinantes/filogenomica/`, `scripts/rapinantes/filogenomica/` | pronta: árvores de genes e wASTRAL nas três execuções; supermatriz, gCF/sCFL e testes AU na execução MAFFT |
| Fósseis e morfologia | `runs/rapinantes/fosseis/`, `scripts/rapinantes/fosseis/` | pronta |
| Estado ancestral | `runs/rapinantes/ancestral/` | pronta |
| Reanálise de B10K e Open Tree | `runs/rapinantes/reanalise/`, `scripts/rapinantes/reanalise/` | pronta |

As três execuções do pipeline:

- `20261018-084205-rapinantes_DNA12`: nucleotídeos, 1ª e 2ª posições do códon. Completa.
- `20261018-112050-rapinantes_AA_famsa_clipkit`: aminoácidos, alinhados com FAMSA e limpos com ClipKIT. Árvores de
  gene e wASTRAL, sem supermatriz.
- `20261018-020210-rapinantes_AA_mafft_clipkit`: aminoácidos, alinhados com MAFFT e limpos com ClipKIT. Completa:
  árvores de genes, wASTRAL, supermatriz particionada (IQ-TREE 3.1.4, 48 grupos de modelos, UFBoot com `--bnni`,
  SH-aLRT, aBayes), gCF e sCFL. A supermatriz foi refeita com os modelos fixados depois de uma interrupção. O
  [`RETOMADA_MANUAL.md`](runs/filogenia/20261018-020210-rapinantes_AA_mafft_clipkit/RETOMADA_MANUAL.md) conta o que
  aconteceu e qual ressalva vai para os Métodos.

### Supermatriz contra árvores de espécies

Os testes de topologia estão em `runs/rapinantes/filogenomica/topologia_AA/`, e os nós-chave em `nos_chave_AA.tsv`.
A supermatriz recupera Australaves (UFBoot 96) e Eufalconimorphae (95), e põe Accipitrimorphae na base de Afroaves,
com as corujas irmãs de Coraciimorphae (UFBoot 90).

O teste AU sobre a supermatriz dá estes resultados:

- rejeita *Cariama* + Falconiformes (p = 0,019);
- rejeita a árvore do wASTRAL inteira (p = 0,0007);
- rejeita Falconiformes s.l., corujas + bacurau, urubus + cegonha e as alternativas para *Sagittarius* (p < 0,001);
- não rejeita Hieraves (p = 0,13).

Na base de Afroaves, os sítios se dividem em 43 % para corujas + Coraciimorphae, 32 % para Hieraves e 24 % para
Accipitrimorphae + Coraciimorphae. O puxão para *Cariama* + Falconiformes aparece nas árvores de genes (41–45 de 300,
contra 6–10 para Eufalconimorphae), e não nos sítios (29,6 % contra 28,6 % para a outra alternativa).

### Árvores de espécies (wASTRAL)

Probabilidade posterior local dos clados, nas árvores enraizadas em *Gallus*. Saída de
`scripts/rapinantes/filogenomica/f07a_clados_wastral.py` e `f07b_irmaos_wastral.py`.

| Clado | FAMSA | MAFFT | DNA12 |
|---|---|---|---|
| Accipitriformes, Accipitrimorphae, Strigiformes, Falconiformes (cada um) | 1,00 | 1,00 | 1,00 |
| Strigiformes + Accipitrimorphae | 0,69 | 0,74 | não (*Leptosomus* + Strigiformes, 0,61) |
| Strigiformes + Accipitrimorphae + *Leptosomus* | 0,59 | 0,56 | 0,90 |
| Falconiformes + *Cariama* | 0,77 | 0,73 | 0,61 |
| Todos os rapinantes + *Cariama* + *Leptosomus* | 0,83 | 0,92 | 0,66 |
| Australaves / Eufalconimorphae | não | não | não |
| Psittacopasserae | 0,91 | 1,00 | 1,00 |

Nas três árvores de espécies, Falconiformes + *Cariama* sai como grupo-irmão de Strigiformes + Accipitrimorphae +
*Leptosomus*. O ramo que define esse agrupamento mede 0,04 a 0,09 unidades de coalescência. A supermatriz e o teste
AU rejeitam esse arranjo (seção anterior). É o artefato de dado codificante que o relatório descreve para os éxons do
B10K e para Kuhl et al. 2021.

A posição de *Apus* + *Calypte* muda com o alinhador. No FAMSA e no DNA12 eles caem como irmãos de Psittacopasserae
(PP 0,59 e 0,67). No MAFFT ficam na base de Neoaves (PP 1,00).

## Estrutura

| Pasta | Conteúdo |
|---|---|
| `reports/` | relatório |
| `research_notes/` | notas de pesquisa por frente, com as fontes |
| `docs/metodos_rapinantes.md` | Métodos |
| `scripts/rapinantes/` | código de cada frente; `filogenomica/retomada_020210/` refaz a supermatriz MAFFT |
| `runs/filogenia/` | as três execuções do pipeline: alinhamentos, árvores, concordância, resumo, proveniência (RO-Crate, BCO) |
| `runs/rapinantes/` | saídas das frentes de fósseis, estado ancestral e reanálise; na filogenômica, as etapas feitas antes do pipeline, os nós-chave e os testes de topologia |

## Como reproduzir

O pipeline de filogenia (Nextflow + IQ-TREE 3 + wASTRAL) não está aqui. Ele fica em
[`bio-stack/pipelines/filogenia`](https://github.com/patoangolano/bio-stack/tree/main/pipelines/filogenia), e as
correções usadas nestas execuções entraram no commit `1ea0923`. Cada execução guarda a linha de comando em
`executar.sh`, os parâmetros em `pipeline_info/parametros.json` e as versões em `pipeline_info/versoes.yml`.

## O que não está aqui

- **Figuras por árvore de gene.** São 600 PNG + SVG por execução. Ficam só as figuras da árvore de espécies e da
  concatenada. As demais se redesenham a partir de `04_arvores_genes/` com `filogenia.py desenhar`.
- **Dados de terceiros baixados.** AVONET, EltonTraits, ERDA/B10K, Kuhl et al. 2021 e Prum et al. 2015 não são
  redistribuídos. O [`MANIFESTO_downloads.tsv`](runs/rapinantes/ancestral/downloads/MANIFESTO_downloads.tsv) dá a URL de
  origem e o md5 de cada arquivo.
- **Sequências de entrada, supermatrizes e checkpoints do IQ-TREE.** Todos se reconstroem a partir dos alinhamentos.

## Origem deste repositório

Este repositório é um espelho da parte "rapinantes" do [bio-stack](https://github.com/patoangolano/bio-stack), que é a
fonte. O espelho é gerado por `scripts/rapinantes/exportar_repositorio.sh`, e cada commit aqui cita o commit do
bio-stack de onde veio.
