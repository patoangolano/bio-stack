# Energia ambiental em microrganismos — protocolo de escopo

Data de início: 2026-11-17. Idioma da síntese: português.

## Questão

Quais microrganismos naturais possuem estratégias diversas de obtenção de
energia ambiental, quais têm genoma público, como funcionam seus mecanismos,
e sob quais condições algum é mais apropriado/eficiente?

## Contrato da revisão

- Revisão de escopo aprofundada e exploratória, não revisão sistemática
  exaustiva nem metanálise. Sem promessa de cobrir toda a diversidade microbiana.
- Separar variedade **dentro de uma linhagem** da variedade entre espécies.
- Cobrir Bacteria, Archaea e microrganismos eucarióticos.
- Evidência principal: trabalhos primários, genomas NCBI e anotações UniProt.
  Revisões para contexto e rastreamento de referências.
- Excluir engenharia, otimização de virulência, protocolos de cultivo e
  organismos de risco elevado. A inclusão ecológica não certifica biossegurança.
- Preservar distinções: fóton / doador de elétrons / aceptor de elétrons /
  fonte de carbono; energia de manutenção / crescimento / produto útil.
- Não confundir CO2 fixado ou N2 assimilado com geração de energia.
- Não comparar percentuais de eficiência com denominadores diferentes.
- Não inferir eficiência de contagem de genes, tamanho de genoma, fotossistema
  isolado, estrutura predita ou anotação automática.
- Fontes consultadas até a data real de execução, sem limite inferior de ano.
  Incluir trabalhos seminais e atualizações verificáveis.

## Frentes independentes

1. Literatura sobre generalistas e versatilidade: agente lit-scout.
2. Mapa amplo de estratégias e especialistas: agente general.
3. Critérios de eficiência e limitações: agente general.
4. Auditoria de anotações funcionais por linhagem: agente evo-devo.
5. Inventário NCBI, NotebookLM e integração: sessão principal.

## Busca e genomas

As buscas de literatura, seus resultados e limitações ficam nos arquivos dos
agentes. Consulta de descoberta inicial ToolUniverse/Europe PMC:
`("Rhodopseudomonas palustris" AND genome AND versatile)`,
total de 339 registros, somente três retornados (não são 339 estudos lidos).

Inventário: lista versionada `scripts/energia_microbiana/candidatos.tsv`;
NCBI Datasets CLI via WSL Ubuntu-26.04, até dez montagens por consulta.
As consultas guardam nome solicitado, filtro de linhagem, data, versão,
retorno bruto, total anunciado e SHA-256. Escolha de representante apenas
entre os resultados retornados; não declara a melhor montagem de todo o NCBI.
Não há download de genomas inteiros nem reanotação de novo.

## NotebookLM

Notebook criado por MCP: `4cb14b1e-7be6-439c-8ec6-f24686551300`.
Pesquisa web profunda:
`ChA1MDQxOTM0NjIyM2I4OWM0EAgaAjBjKgN1c2M`.
URL: https://notebooklm.google.com/notebook/4cb14b1e-7be6-439c-8ec6-f24686551300

O NotebookLM é ferramenta auxiliar de descoberta e crítica. Suas respostas
não substituem os estudos originais; suas referências precisam de conferência.
Enviar apenas literatura pública e sínteses sem informações privadas.

## Materiais anteriores

Arquivos anteriores em `scripts/microbios_energia/`,
`runs/microbios_energia/`, `scripts/bioenergetics_efficiency_analysis.py`
e `runs/bioenergetics_efficiency_atlas.json` serão preservados.
Valores numéricos fixos sem rastreabilidade desse atlas não são evidência
empírica nem resultado de FBA. O relatório novo não os adotará como medições.

## Eficiência: comparações permitidas

- Rendimento de biomassa por mol de substrato e por energia fornecida.
- Produtividade por tempo/volume/área, com unidade e condições.
- Eficiência energética com entrada e saída explicitadas.
- Afinidade e limiar de uso de substrato, separados de rendimento.
- Sobrevivência/manutenção versus crescimento balanceado.
- Recursos e aceptores realmente presentes; custo de fornecer luz/H2/O2.
- Não criar vencedor universal nem um escore numérico sem validação.
