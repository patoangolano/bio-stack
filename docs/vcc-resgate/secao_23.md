
## 23. A v5, e uma predição que foi falsificável e sobreviveu

`submissao_v5.h5ad` carrega `delta_calibrado_kd4.npz`. Aprovada nos dois portões:

```
contrato   28 de 28 itens, varredura de 1.982.711.741 entradas em 117,4 s
numérico   D +0,9982  E +0,9986  F +1,0010   viés −7,2e-04 a −8,1e-04
```

360.000 células, 1.982.711.741 não nulos (5.507,5/célula, **41,74%** do teto de
4,75e9), 3.881,0 MiB, 334,6 s, pico de 1.561,5 MiB. Maior soma por célula 66.829
contra teto de 1.000.000.

### A correlação numérica subiu de 0,43 para 0,87, e isso era uma aposta

A seção 20 registrou que a correlação de **0,43** no validador numérico da v4 não
era infidelidade do arquivo — a inclinação de 0,998 já provava fidelidade em
expectativa. A explicação dada foi que o rms de amostragem por grupo de 400
células (0,038) excede o |delta| típico (0,0018), então o resíduo domina a
correlação sem deslocar a inclinação.

Essa explicação faz uma predição quantitativa e falsificável. Se a correlação é
`s/√(s²+n²)` com `n` fixo pela amostragem, amplificar `s` por 4 tem de levar a
correlação a um valor calculável **antes** de medir:

| contexto | `s/n` implicado na v4 | corr. prevista a 4× | **corr. medida na v5** | amplificação efetiva |
|---|---:|---:|---:|---:|
| D | 0,4763 | 0,8854 | **0,8701** | 3,71× |
| E | 0,4584 | 0,8779 | **0,8588** | 3,66× |
| F | 0,4729 | 0,8841 | **0,8650** | 3,65× |

A correlação foi para onde a explicação dizia. O resto da diferença — 3,68× em
vez de 4× — tem a causa certa e no sentido certo: `k_desvio` escala **só o
desvio específico**, e a resposta média continua em ×1. A razão 3,68/4 = 0,92
implica que a parcela da média responde por ~8% da variância por alvo dentro do
universo de DE, o que é coerente com a decomposição medida na §21 (`|média|` 2,35
contra `|desvio|` mediano 5,46).

Isso não prova que o modelo está certo. Prova que o diagnóstico da correlação
baixa estava certo, e fecha a última ressalva que o validador numérico carregava.

### Nada submetido

Cinco arquivos em disco, 2 submissões por dia, placar oculto até 05/11. A
decisão de enviar é do Мухамед, e `runs/vcc/submissao/LEIA-ME.md` diz qual
`--delta` casa com qual arquivo — passar o par errado reprova um arquivo bom.
