
## 25. A v6, e por que ela não pode ser comparada à v5 no instrumento que temos

`submissao_v6.h5ad` carrega `delta_calibrado_kd4_gate.npz` — o mesmo `k_desvio = 4`
da v5 sobre o agrupamento de média corrigido (gate separado, §24).

```
360.000 células, 1.982.874.616 não nulos (5.508,0/célula, 41,74% do teto)
3.881,4 MiB, 366,2 s, pico de 1.561 MiB
contagem máxima por célula 72.187 contra teto de 1.000.000
```

### Quase idêntica à v5, e esse é o dado

| contexto | rms(v5) | rms(v6) | rms(v6−v5) | correlação |
|---|---:|---:|---:|---:|
| D | 0,06745 | 0,06747 | 0,00612 | **+0,9959** |
| E | 0,06542 | 0,06543 | 0,00565 | **+0,9963** |
| F | 0,06597 | 0,06598 | 0,00567 | **+0,9963** |

A correção do gate muda a submissão em **9% em rms**. Por construção: o desvio
específico é idêntico (verificado a 6 partes por milhão, §24) e só a parcela da
média se move — e a média é a parcela de magnitude menor por entrada, mesmo
valendo para os 300 alvos contra 157 do desvio.

### A comparação que NÃO é possível, e dizer isso é o ponto

O `avaliador_proxy` tem verdade Orion/HCT116 transplantada, que **contém** o
desvio específico e **não contém** a resposta média de D. A §19 provou isso em
número: a correção do viés de seleção mudou a magnitude da média de 0,063 para
0,019 e moveu `avg_score` em +0,0028, dentro do ruído.

A v6 difere da v5 **exclusivamente** na parcela da média. Então pontuar a v6 no
proxy mede, na melhor hipótese, nada — e na pior, confunde: um ±0,005 de ruído
leria como veredito. O que o proxy pode dar aqui é **teste nulo**, não
comparação: se a v6 cair materialmente, algo vazou para o desvio. O teste nulo
já foi feito de forma mais direta e mais barata em espaço de delta (§24, resíduo
de 2,27e-07 contra rms 0,0354).

**A evidência a favor da v6 é `pool_media_sem_mesmo_repo.json`, não o proxy.**
Tirar o gate ganha em 4 de 4 configurações de corte, com ganho mediano de
+0,0983 contra +0,0524 do esquema antigo. É medida em dado real, de transferência
entre contextos que existem, sem atalho de mesmo-estudo. Não é `avg_score`, e não
há conversão entre as duas escalas.

### A cauda subiu, e a razão confirma que a correção funcionou

`escala_de_gene_maxima` vai de 199,17 (v5) para **230,16** (v6), e segue sendo a
mesma entrada: `SIN3B` → `DDIT3`. O perfil é igual (47 e 55 entradas acima de
10× contra 46 e 54; p99,99 em 3,52 e 3,61 contra 3,50 e 3,59).

A diferença em log2 é +0,209, cerca de nove vezes o rms da própria média (0,023).
Isso não é ruído, é a correção operando: `DDIT3` é CHOP, e o que entrou no
agrupamento foi justamente o Replogle — rastreio de genes essenciais, com
assinatura de estresse forte e |média| de 0,18 a 0,27 contra 0,033 do
`arc_vcc_25`. Com 4,4% de peso e uma média 5 a 8 vezes maior, o Replogle move um
gene de estresse exatamente nessa ordem de magnitude. O mecanismo e o sinal
batem.

### Qual arquivo enviar, e por que a resposta não é automática

A v6 domina a v5 **se** o resultado de `pool_media` transferir para D/E/F. As
duas melhorias são ortogonais por construção — `k` escala só o desvio, o gate
toca só a média — então não competem, somam.

O que pesa contra: o ganho de `k=4` foi medido no instrumento **certo para a
pergunta** (+0,372 de `avg_score`, 133× o ruído), enquanto o ganho do gate foi
medido num instrumento **analógico** (transferência entre contextos doadores),
em unidade que não converte, sobre uma diferença de 9% em rms.

Nada aqui é submetido. As duas estão em disco, as duas passam os dois portões, e
`runs/vcc/submissao/LEIA-ME.md` diz qual `--delta` casa com qual arquivo.
