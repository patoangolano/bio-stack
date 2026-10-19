# Retomada manual da execução 20261018-020210 (AA, MAFFT + ClipKIT)

Esta nota registra por que a supermatriz desta execução não saiu do processo `IQTREE_CONCAT` do Nextflow, e como ela foi
refeita. Os passos seguintes (concordância e resumo) também rodaram fora do Nextflow.

## O que aconteceu

| Quando | Evento |
|---|---|
| 2026-10-18 08:38 | `IQTREE_CONCAT` começa: `iqtree -m MFP+MERGE -B 1000 --bnni --alrt 1000 --abayes --seed 12345 -o Gallus_gallus -st AA -s supermatriz.fa -spp particoes.nex -pre concatenado -T 8` (IQ-TREE 3.1.4, imagem `quay.io/biocontainers/iqtree:3.1.4--h8471819_0`). Os dados são 57 táxons, 215 366 colunas e 300 partições. |
| 2026-10-18 14:25 | O ModelFinder termina. O agrupamento relaxado (RClusterF, 10 %, até 3000 pares) para em 6 passos, com **48 grupos**: BIC 4467486.907, LnL -2230458.528, df 535. Grava `concatenado.best_scheme.nex`. Tempo de CPU do ModelFinder: 43 h 19 min. |
| 2026-10-18 17:16 | A busca de árvore termina na iteração 100, com o melhor LogL das candidatas em -2230208.324. Começa o refinamento dos 1000 UFBoot (`--bnni`). |
| 2026-10-18 21:09 | **SIGHUP.** A sessão do Claude Code que tinha lançado o Nextflow foi fechada; o `wsl.exe` morre e o sinal desce até o IQ-TREE. O checkpoint `concatenado.ckp.gz` fica com `refined_samples: 377`. |
| 2026-10-19 03:37 | Retomada a partir do checkpoint, com o mesmo comando e no mesmo diretório de trabalho (`work/81/e3fdfb…`), por tarefa agendada que sobrevive ao fim da sessão. Antes, os arquivos foram copiados para `*.ate_sighup_20261018`. |
| 2026-10-19 05:57 | A retomada é interrompida (detalhes abaixo). |
| 2026-10-19 05:58 | Começa a nova execução com os 48 modelos fixos (`manual/020210/fixo/run`), com `-T 8`. |
| 2026-10-19 09:41 | Teste numa cópia descartável (`fixo/teste_t18`): retomar o checkpoint da execução fixa com `-T 18` restaura a árvore inicial, os 1000 UFBoot, os parâmetros do modelo (LogL -2230241.611) e as árvores candidatas. A busca continua na iteração 74, sem reotimizar parâmetros. No teste com `-p`, descrito abaixo, o modelo mudava; aqui muda só o número de threads. |
| 2026-10-19 09:42 | O FAMSA 112050 terminou e liberou os núcleos. A execução fixa é parada depois da iteração 76 (cópia do checkpoint em `fixo/ckp_antes_T18/`) e relançada com `-T 18`. Ela retoma na iteração 76. |

## Por que a retomada foi abandonada

Ao retomar com `-m MFP+MERGE`, o IQ-TREE 3.1.4 entrou de novo no ModelFinder. O cache `concatenado.model.gz` guarda
cada partição isolada apenas com os 2 modelos da fase de agrupamento (LG+I+G4 e LG+F+I+G4); o resultado completo dessas
partições não está nele. O IQ-TREE passou então a testar de novo cada partição isolada com o conjunto completo (CCP110 subiu para 271
entradas no cache, HEATR5A para 208), a cerca de 7 min por partição com 8 threads. Em 2 h 20 min avançou 27 partições de
cerca de 290; projeção: ~32 h só de ModelFinder, antes de voltar à árvore. Havia ainda o risco de esse ModelFinder
chegar a um esquema diferente daquele sob o qual a árvore e os bootstraps do checkpoint foram calculados.

Testou-se também, numa cópia descartável, retomar o checkpoint com os modelos fixos (`-p`). O IQ-TREE aceita a linha de
comando diferente (só avisa) e restaura os 1000 UFBoot e as 20 árvores candidatas, mas reinicia a busca
(`Number of iterations: 0`, `OPTIMIZING CANDIDATE TREE SET`). Não economiza tempo, e misturaria bootstraps de duas
buscas.

## O que roda agora

Execução nova e limpa, com o mesmo alinhamento, o mesmo esquema e os mesmos modelos que o ModelFinder escolheu na
execução interrompida, a mesma semente e as mesmas opções. Só a busca de modelo fica de fora:

```
iqtree -p modelos_fixos.nex -B 1000 --bnni --alrt 1000 --abayes --seed 12345 -o Gallus_gallus -st AA \
  -s supermatriz.fa -pre concatenado -T 8
```

- `modelos_fixos.nex` foi reconstruído por `modelos_48.py` a partir do próprio `concatenado.model.gz`. Para cada um dos
  48 grupos do `best_scheme.nex`, o script escolhe o modelo de menor BIC, o mesmo critério do ModelFinder. A soma
  reproduz exatamente o original (LnL -2230458.528, df 535); sem essa conferência o script não grava nada. As matrizes
  escolhidas: Q.BIRD em 27 grupos, Q.MAMMAL em 7, JTT em 4, VT em 1 e Q.PLANT em 1, quase sempre com taxas livres (+R).
  O grupo AP2M1/GNB4/KLHL15/YWHAE ficou com JTT puro, sem heterogeneidade de taxa. Os 8 grupos de uma partição só ficaram
  com LG+I+G4.
- `-p` é o modelo de partição proporcional nas arestas, o mesmo do `-spp` original.
- Até a iteração 76 a busca rodou com `-T 8`; dali em diante, com `-T 18` (`roda_fixo.sh` recebe o número de threads como
  argumento). O número de threads muda a ordem das somas em ponto flutuante. A busca não muda, mas repetir a execução
  inteira com um único `-T` pode dar comprimentos de ramo diferentes na 8.ª casa decimal e SH-aLRT diferente em até
  ~1 ponto. A semente é a mesma nas duas fases.
- Ao terminar, `pos_concat_020210.sh` publica `05_concatenado/` e roda a concordância (gCF e sCFL, como no módulo
  `IQTREE_CONCORDANCIA`) em `06_concordancia/`. Depois roda o `RESUMO` (`resumo.tsv`, `resumo.json`, `metodos.md`) com
  as mesmas imagens e argumentos do pipeline. A única diferença está no sCFL: ele roda com `-T 8`, e no módulo roda com
  `-T 1`. Os valores podem diferir de uma execução com `-T 1` na mesma escala citada acima. O gCF só conta árvores e não
  muda.

## Onde está cada coisa

| O quê | Onde |
|---|---|
| Scripts | `scripts/rapinantes/filogenomica/retomada_020210/` (`modelos_48.py`, `roda_fixo.sh`, `pos_concat_020210.sh`, `lanca_fixo.ps1`) |
| Execução nova | WSL `~/nf-work/filogenia/manual/020210/fixo/run/`; log em `fixo/supermatriz_fixa.log` |
| Checkpoint antes da troca para `-T 18` | WSL `~/nf-work/filogenia/manual/020210/fixo/ckp_antes_T18/concatenado.ckp.gz` (09:41:28) |
| Execução interrompida | WSL `~/nf-work/filogenia/work/81/e3fdfb681f0da353c2d65ff89f587d/`, com backups `*.ate_sighup_20261018` |
| Log da retomada abandonada | WSL `~/nf-work/filogenia/retomada_concat_020210.log` |
| Tarefa agendada | `\bio-stack\supermatriz_fixa_020210` (remover ao terminar) |

## Consequência para os Métodos

A árvore concatenada desta execução vem de uma análise de máxima verossimilhança com o modelo particionado selecionado
pelo ModelFinder com agrupamento (MFP+MERGE), reaplicado como modelo fixo numa segunda corrida idêntica nas demais opções.
O `metodos.md` gerado pelo `RESUMO` descreve o comando do pipeline (MFP+MERGE); esta nota é a ressalva que o acompanha.
