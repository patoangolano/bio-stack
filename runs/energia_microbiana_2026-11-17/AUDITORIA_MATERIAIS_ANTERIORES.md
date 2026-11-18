# Auditoria de materiais preexistentes

Os arquivos anteriores foram mantidos intactos. Esta nota limita como podem ser
utilizados como evidência; não atribui autoria nem presume a origem deles.

## 1. Percentuais fixos não são medições

`scripts/bioenergetics_efficiency_analysis.py` cria registros com valores
constantes para eficiência termodinâmica, rendimento faradaico, ATP por doador
e velocidade de crescimento, sem vincular cada valor a um estudo/condição.
Esses números não podem sustentar uma comparação quantitativa entre espécies.

## 2. O nome da função não demonstra uso de COBRApy

`run_cobrapy_bioenergetic_benchmark()` contém valores de regimes previamente
fixados e expressões algébricas simplificadas. O arquivo não importa COBRApy
nem carrega modelos metabólicos; a função não executa otimização de uma rede
genômica. `scipy.optimize.linprog` é importada, mas não chamada.

Além disso, a comparação fotoautotrófica adiciona luz sem contabilizar seu
custo energético no mesmo denominador usado para os demais cenários.
Retenção de elétrons em acetato/biomassa e conservação de energia em ATP
também são grandezas distintas.

## 3. Linhagens e montagens precisam ser separadas

- A entrada `CGA009 / TIE-1` funde linhagens distintas de R. palustris.
- NCBI confirmou CGA009 em `GCF_000195775.2`, enquanto o arquivo antigo usa
  a versão anterior `.1`.
- Para C. necator H16, `GCF_000009285.1` foi recuperado. Não se deve formar
  `GCF_000009285.2` trocando automaticamente o prefixo de um acesso GenBank.
- Há montagens históricas `suppressed`; isso não significa que a espécie
  nunca tenha sido sequenciada.
- Número atual de CDS é dependente da versão da anotação, não uma
  constante biológica nem prova de maior eficiência.

## 4. Substituição feita nesta pesquisa

A auditoria real encontra-se em `modelos/auditoria_cobrapy.json`, gerada por
`scripts/energia_microbiana/auditar_modelos.py` com COBRApy 0.32.1 e GLPK.
Ela registra SHA-256 dos quatro modelos locais preexistentes, os objetivos e
os limites de importação. Os arquivos-modelo não foram alterados.

As quatro otimizações sob limites nativos foram viáveis. Com importações
fechadas, três modelos ficaram inviáveis por exigirem manutenção mínima;
isso não foi convertido em crescimento zero. Em três testes separados com
a manutenção mínima relaxada para zero, a maximização da hidrólise de ATP
retornou zero. No modelo restante não foi identificada uma reação ATPM
inequívoca e esse teste não foi realizado.

Esses resultados são uma verificação limitada de consistência, não uma
validação termodinâmica completa nem um ranking de eficiência.
