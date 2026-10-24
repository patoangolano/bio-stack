<!-- Auditoria de skills e conectores contra os gargalos medidos do VCC 2026. -->

# Ferramentas auditadas contra os gargalos medidos

Dados de 2026-10-23. A auditoria foi feita contra **gargalo medido**, não
contra novidade: cada item abaixo diz qual número ele move.

---

## Tier 1 — resolve um gargalo medido

### 1. Colab por CLI — `google-colab-cli` 0.7.5

**Estado: instalado, roteiro escrito, bloqueado em autenticação interativa.**

Instalado como `uv tool` (`colab` em `~/.local/bin`). É o **sexto** uv tool da
máquina, todos no cpython-3.13 — relevante porque `winget install uv` zera
`AppData/Roaming/uv` e já aconteceu antes.

**O gargalo que ele move.** O `avaliador_proxy` cabe **48 perturbações** nos 15
GiB da VM do WSL (medido: pico 10,855 GiB para 31.200 células, 0,3303 GiB por
mil células). O painel oficial de um contexto projeta **46,3 GiB** — não cabe
nem nos 31,4 GiB do host. Isso é a diferença entre pontuar 16% do painel e
pontuar o painel.

**Duas decisões de configuração, e a segunda não era óbvia:**

1. **`--high-mem` sem acelerador.** O gargalo é RAM, não GPU.
2. **GPU quebraria a comparação.** `resolved_device` é um dos **oito campos de
   identidade** que `score --real-bundle` compara, e divergência é fatal sem
   perdão. O bundle foi construído em CPU e os oito batem
   (`todos_os_oito_peers_batem: true`). Subir numa GPU mudaria esse campo.
   Então CPU de RAM alta é mais barato **e** mais correto que A100.

**Terceira decisão: regenerar em vez de subir.** `submissao_v3.h5ad` tem 3,9
GiB. O delta que a gerou tem **27 MB**, e o `emissor_400` é streaming. Então
sobe-se o delta mais os controles (~672 MB no total) e emite-se lá.

**O que falta, e só você pode fazer.** A CLI pede autorização de forma
interativa — imprime URL, espera código colado. Rode **uma vez**:

```
colab usage
```

Depois disso `scripts/vcc/colab_avaliador.py --preparar` e `--rodar` são
headless. O script detecta a falta de autorização e diz isso em vez de falhar
de forma obscura.

### 2. Levantamento de literatura — PubMed, bioRxiv, Consensus, alphaXiv

**Estado: em curso.** Saída em `docs/VCC-2026-LITERATURA.md`.

**O gargalo que move.** Na validação de 2025 o 1º lugar fez 0,472 e o 2º ao 10º
ficaram entre 0,32 e 0,35. Nosso modelo rende ~4,5% de melhora de MSE sobre a
resposta média. Saber **o que o vencedor de fato fez** vale mais que qualquer
ajuste incremental, e decide se vale investir nos itens do Tier 2.

Perguntas encomendadas: método das equipes de topo; pontuação do baseline
oficial; se modelos profundos batem baseline linear em contexto novo *e*
perturbação nova (os dois lados, com números); como as métricas de DE são
ganhas; e o que prediz transferibilidade, já que medimos que proximidade de
linhagem **não** prediz.

---

## Tier 2 — melhoria plausível, com ressalva registrada

### 3. Modelos de perturbação — `envs/perturb-scvi/`

**Estado: em construção.** `pertpy` 1.4.0, `scvi-tools` 1.5.1, `scgen` 2.1.0,
torch 2.11.0+cu128, Python 3.12.

A verificação do `montar.sh` **não** usa `cuda.is_available()`, que mente numa
sm_120 com wheel errado: roda matmul fp16 na GPU mais um passo de autograd. Um
wheel cu121 instala sem reclamar nesta máquina e falha no primeiro kernel.

### 4. CPA — `envs/perturb-cpa/`

**Estado: configurado, NÃO construído.** E o motivo é um conflito duro que vale
registrar:

```
cpa-tools 0.8.8   exige python <3.11,>=3.9
scvi-tools 1.5.1  exige python >=3.12
```

**Mutuamente exclusivos.** Não há versão de Python que satisfaça os dois.
Precedente no projeto: há três OpenSim nesta máquina em ambientes separados
exatamente por isso. Não unifique.

Não construí porque o custo é ~3 GB e o retorno é o mais incerto dos dois: CPA
foi desenhado para dose-resposta e combinação em contexto **já visto**, que não
é o nosso regime.

### 5. Geneformer e scGPT como features dos alvos

**Estado: caminho documentado, não executado.**

`geneformer` **não está no PyPI** — vive em `ctheodoris/Geneformer` no Hugging
Face, acessível pelo conector já autenticado. `scgpt` 0.2.4 está no PyPI.

**Por que não executei, e é a ressalva mais importante desta página.** Nosso
teste de features (`scripts/vcc/alvos_sem_evidencia.py`) usou `esm2`,
`depmap_gene_effect` e `string_space`, que cobrem **300 de 300** alvos, e deu
**R² negativo em todos os seis valores de α e nos três contextos** — a predição
por ridge é pior que prever a carga média.

R² negativo em toda a varredura não é feature ruim, é **rótulo ruim**: em 148
dos 157 casos o rótulo de treino vem do Orion, com mediana de **dois** genes
significativos por perturbação. Trocar a feature não conserta o rótulo.

O arnês de validação cruzada já existe e aceita qualquer matriz de feature, então
testar Geneformer custa pouco — mas **só vale depois** de melhorar o rótulo, não
antes.

### 6. Revisão adversarial do código

**Estado: em curso.** Seis scripts numéricos sob revisão.

**O gargalo que move.** Eu mesmo achei três bugs reais hoje: centragem
degenerada quando o contexto doador tem uma só perturbação; gate duplicado no
emissor zerando a resposta média; e `glob` misturando `.npz` de um pipeline
anterior (31 direções onde existem 24). Cada bug aqui custa uma submissão, e há
duas por dia com placar oculto.

---

## Tier 4 — considerado e recusado, com o motivo

### Synthesize Bio (GEM-1) — dado de expressão gerado por IA

**Estado: conector autenticado e disponível. Recusado por ora.**

A ideia era tentadora: gerar perfis de knockdown para os **143 alvos sem
evidência** e usá-los como doador. Duas objeções, e a segunda é decisiva.

1. **Circularidade.** O VCC é pontuado contra dado **real**. Treinar numa
   predição de outro modelo é prever o modelo, não a célula.
2. **O prior já foi medido e é ruim.** Seria uma quarta fonte de conhecimento,
   e as três que cobrem 300/300 deram R² negativo. Adicionar uma fonte com viés
   desconhecido a um conjunto que já não carrega sinal não é aposta favorável.

**O teste que mudaria a decisão** está pronto: a mesma validação cruzada de
`alvos_sem_evidencia.py`, com a matriz de feature trocada. Se o levantamento de
literatura apontar que fontes de conhecimento funcionam no nosso regime, roda-se
em minutos.

### NVIDIA NIM — `mcp__nvidia-bio__*`

**Recusado: domínio errado.** São estrutura de proteína, docking e geração de
molécula. A tarefa é expressão. Nenhum dos 11 endpoints vivos toca resposta
transcricional a perturbação.

### BioRender e `dataviz`

**Recusado para a nota, guardado para o artigo.** Figura não pontua no VCC. Vale
quando a seção de Métodos for escrita — e aí `runs/provenance.jsonl` mais
`docs/VCC-2026-EXECUCAO.md` são a matéria-prima.

### `biocontext`, `OpenTargets`, `gget`, `tooluniverse`

**Recusados por redundância.** São lookup de gene e de alvo. O que o projeto
precisava deles (`gget_enrichr` nos 300 alvos) já foi feito, e está em
`runs/vcc/caracterizacao/`.

---

## O que de fato muda o número, em uma linha

O avaliador no Colab. Todo o resto é incremental sobre um modelo cujo limite
medido não é a arquitetura — é que 148 dos 157 alvos com evidência têm essa
evidência num ensaio com mediana de dois genes significativos por perturbação.
