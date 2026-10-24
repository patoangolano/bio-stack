
## 26. O validador numérico não distinguia v5 de v6, e isso era medível

A §20 criou o `validar_numerico.py` porque os 28 itens de contrato aprovavam o
`stub.h5ad` de delta zero. A pergunta que ele responde é "este arquivo carrega o
delta do modelo atual", e o mecanismo é regressão: recompõe o pseudobulk emitido
e regride contra o delta pretendido, exigindo inclinação em [0,85; 1,15].

Com cinco versões em disco esse teste tinha sido exercitado só contra erros
grandes: v1 dá 0,809, o `stub` dá −0,001. A v5 e a v6 diferem **apenas** na
parcela da média, com correlação +0,996 — então valia perguntar se ele ainda
discrimina.

### Previsto antes de medir, e medido depois

A inclinação é `(pedido · recuperado) / (pedido · pedido)`, e o ruído de
amostragem não a desloca. Com pedido = delta da v5 e recuperado ≈ v6, isso é
`(v5 · v6)/(v5 · v5)`, calculável do disco:

```
previsto    D +0,9962   E +0,9965   F +0,9964
medido      D +0,9943   E +0,9951   F +0,9975
```

**O validador aprovava a v6 contra o delta da v5.** Viés de −0,0004 a −0,0011
contra um limite de 0,005 — também passava.

### Por que apertar o limite é o conserto errado

Separar 0,996 de 1,000 exigiria uma janela de ±0,004. O rms do resíduo de
amostragem por grupo de 400 células é 0,038 contra |delta| típico de 0,002, e é
esse ruído que a janela de [0,85; 1,15] existe para tolerar. Qualquer janela
estreita o bastante para separar v5 de v6 reprovaria o par correto na próxima
emissão com outra semente.

A pergunta não é de magnitude, é de **identidade** — e identidade se checa com
impressão digital, não com regressão.

### O conserto

`emissor_400.py` passa a gravar `delta_sha256` no `medicao_<nome>.json`, ao lado
do caminho que já gravava. `validar_numerico.py` passa a ler esse arquivo e
comparar antes de qualquer aritmética: sha256 quando existe, nome do arquivo
quando não. Ausência do arquivo de medição ou do campo não reprova — é a situação
das emissões anteriores a esta sessão — mas imprime que a cobertura então é só a
da regressão, em vez de deixar implícito.

Medido depois do conserto, no mesmo par cruzado:

```
[FALHA] o emissor usou delta_calibrado_kd4_gate.npz e o --delta passado e
        delta_calibrado_kd4.npz; sao arquivos DIFERENTES
REPROVADO   (saída 1)
```

A regressão continua rodando e continua reportando, porque ela pega uma coisa que
a impressão digital não pega: emissão **infiel** do delta certo — arredondamento
errado, normalização errada, gate acidental. Os dois testes cobrem falhas
diferentes, e é por isso que o conserto adiciona em vez de substituir.

### O padrão, porque é o terceiro da mesma família nesta sessão

Três vezes um portão aprovou o que não devia, e nas três a causa foi a mesma
forma: **o teste media uma propriedade correlacionada com a que importava, em vez
da que importava.**

| portão | aprovava | porque media |
|---|---|---|
| 28 itens de contrato | `stub.h5ad`, delta zero | estrutura, não conteúdo |
| checagem biológica da v3 | 87% de artefato de cobertura | se os genes fazem sentido, não a fração de sinais |
| `validar_numerico` | v6 contra o delta da v5 | magnitude por regressão, não identidade |

A pergunta que gera o portão útil não é "o resultado parece certo" — é "o que eu
veria se estivesse errado, e esse teste veria?".
