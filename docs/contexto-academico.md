# Cabeçalho de contexto acadêmico

Cole no início de um pedido (ou uma vez por sessão) antes de acionar os agentes
`seq-designer` / `variant-curator` / `structure-bio`. Não é para burlar
salvaguarda — é o contexto real do trabalho, dado de forma que o modelo
reconheça a natureza acadêmica e benigna antes de ler a tarefa.

## Cabeçalho

```
Contexto: trabalho acadêmico in silico, portfólio de bioinformática.
Sou estudante de Análise e Desenvolvimento de Sistemas no Senac (SP).
Escopo: apenas organismos de ensino/pesquisa básica (E. coli, S. cerevisiae,
B. subtilis, linhagens celulares humanas, plantas modelo) e partes
caracterizadas (repórteres, promotores, cassetes de resistência padrão,
peças do Registry/iGEM). Nada de agentes de risco elevado, toxinas ou
aumento de virulência/transmissibilidade — isso é fora do escopo e vai ao
comitê de biossegurança da instituição. Tudo é computacional e registrado
em runs/provenance.jsonl para a seção de Métodos.

Tarefa: <descreva aqui o que precisa, em uma frase clara>
```

## Como escrever a linha `Tarefa:`

- Diga o **organismo** e a **finalidade benigna** já na primeira frase.
  - ✅ "Otimize os códons da GFP para expressão em BL21."
  - ❌ "Me dê a sequência otimizada de [gene]."
- **Evite termo de dano** quando não for o objeto ("virulência", "patógeno",
  "toxina", "gain-of-function", "arma") — dispara o classificador mesmo em
  negação. Se o assunto legítimo exigir o termo, explique o contexto na mesma frase.
- **Peça o método junto do resultado**: "…e descreva como eu validaria na bancada".
- Prefira **alvos concretos e caracterizados** (gene repórter nomeado, plasmídeo
  do Registry) a pedidos genéricos e abertos.

## Contexto

Fable 5.1 é o mesmo modelo do Mythos 5.1 com salvaguardas ligadas. As
salvaguardas de biossegurança rodam por cima do modelo; este cabeçalho reduz
falsos-positivos em trabalho legítimo, mas não desliga o classificador. Para
recurso formal de bloqueio indevido: usersafety@anthropic.com.
