# Como escrever um pedido de bioinformática com escopo explícito

Um pedido que declara organismo, finalidade e critério de validação produz
resposta melhor — e é o mesmo hábito que um caderno de laboratório exige. Este
cabeçalho existe para tornar isso automático, não para ser burocracia.

Cole no início de um pedido (ou uma vez por sessão) antes de acionar os agentes
`seq-designer` / `variant-curator` / `structure-bio`.

## Cabeçalho

```
Contexto: trabalho acadêmico in silico, portfólio de bioinformática.
Estudante de Análise e Desenvolvimento de Sistemas no Senac (SP).
Escopo: organismos de ensino/pesquisa básica (E. coli, S. cerevisiae,
B. subtilis, linhagens celulares humanas, plantas modelo) e partes
caracterizadas (repórteres, promotores, cassetes de resistência padrão,
peças do Registry/iGEM). Tudo é computacional e registrado em
runs/provenance.jsonl para a seção de Métodos.

Tarefa: <descreva aqui o que precisa, em uma frase clara>
```

## Como escrever a linha `Tarefa:`

- Diga o **organismo** e a **finalidade** já na primeira frase.
  - ✅ "Otimize os códons da GFP para expressão em BL21."
  - ❌ "Me dê a sequência otimizada de [gene]."
- **Peça o método junto do resultado**: "…e descreva como eu validaria na
  bancada". Um resultado sem critério de validação não entra no relatório.
- Prefira **alvos concretos e caracterizados** (gene repórter nomeado,
  plasmídeo do Registry) a pedidos genéricos e abertos. Pedido vago devolve
  resposta vaga.
- Declare o **critério de sucesso** quando houver: qual número, qual limiar,
  qual controle decide se deu certo.

## Fora de escopo

Trabalho com agentes de risco elevado, toxinas, ou qualquer coisa que aumente
virulência ou transmissibilidade **não é objeto deste projeto** e tem via
institucional própria: o comitê de biossegurança da instituição. Isso não é
uma formalidade contornável — é o procedimento correto, e este repositório
não substitui nem antecipa essa análise.
