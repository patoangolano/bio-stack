# Auditoria computacional de modelos — interpretação

## O que foi executado

COBRApy 0.32.1, solver GLPK, Python 3.11 em ambiente isolado de `uv`.
Não foi necessário alterar o lock ou a venv de análises da bio-stack.
Os quatro arquivos-modelo já estavam disponíveis localmente; a origem
externa e a equivalência exata à versão publicada não foram certificadas.
Os hashes, versões, objetivos e limites foram preservados no JSON.

| Arquivo | Identificador lido | Reações | Metabólitos | Genes | Objetivo nativo |
|---|---|---:|---:|---:|---|
| iCN1361.xml | CnecH16 | 1292 | 1265 | 1362 | R_Biomass |
| iJN678.json | iJN678 | 863 | 795 | 622 | BIOMASS_Ec_SynHetero |
| iHN637.json | iHN637 | 785 | 698 | 637 | BIOMASS_Cl_DSM_WT_46p666M1 |
| iAF692.json | iAF692 | 690 | 628 | 692 | BIOMASS_Mb_30 |

Todas as otimizações nativas foram `optimal`. Isso significa que o problema
matemático tem solução sob os limites de cada arquivo; não significa que as
espécies tiveram taxas de crescimento medidas nesta pesquisa.

## Por que não há ranking FBA

1. Fontes de energia e limites de importação diferem entre arquivos.
2. As equações de biomassa e as exigências de manutenção diferem.
3. iJN678 chegou com objetivo heterotrófico: seu resultado padrão não é uma
   comparação fotossintética com os demais organismos.
4. FBA estequiométrica não determina afinidade, cinética, abundância enzimática,
   custo do fornecimento de gases, iluminação ou eficiência de um reator.
5. Comparar os valores máximos desses objetivos seria comparar problemas
   diferentes, não estabelecer a espécie biologicamente mais eficiente.

## Teste limitado de energia sem entrada

Foram fechadas as importações de reações de fronteira monometabólito,
incluindo possíveis fontes de demanda/sumidouro reversíveis. Com as exigências
nativas de manutenção, iCN1361, iHN637 e iAF692 ficaram inviáveis; iJN678
teve objetivo zero. Após reduzir o piso de ATPM a zero e maximizar sua
hidrólise, iCN1361, iHN637 e iAF692 retornaram zero. iJN678 não tinha uma
ATPM inequívoca reconhecida pelo roteiro, portanto não foi testado nesse item.

Fluxo zero nesses testes não certifica ausência de todos os ciclos indevidos
sob outros meios. Nenhuma capacidade nova ou eficiência percentual foi inferida.

## Literatura verificada para os modelos

- Pearcy et al. (2022), iCN1361:
  https://doi.org/10.1371/journal.pcbi.1010106 — PMID 35604933, PMC9166356.
- Nogales et al. (2012), Synechocystis:
  https://doi.org/10.1073/pnas.1117907109 — PMID 22308420, PMC3289291.
- Nagarajan et al. (2013), iHN637:
  https://doi.org/10.1186/1475-2859-12-118 — PMID 24274140, PMC4222884.
- Feist et al. (2006), iAF692:
  https://doi.org/10.1038/msb4100046 — PMID 16738551, PMC1681478.
- Tec-Campos et al. (2023), iDT1294, outra linhagem de R. palustris:
  https://doi.org/10.1371/journal.pcbi.1011371 — PMID 37556472, PMC10441798.
  Este quinto modelo não foi executado; a referência demonstra a importância
  de uma reconstrução especificamente validada para Bis A53.

Metadados e resumos recuperados via ToolUniverse/Europe PMC. Os quatro
primeiros trabalhos e o quinto também foram adicionados ao NotebookLM
para aprofundamento. Não afirmar equivalência binária entre artigos e
arquivos locais apenas por compartilharem o nome de modelo.
