# Estilos de citação (CSL)

| Arquivo | Origem | Commit | sha256 | Licença |
|---|---|---|---|---|
| `associacao-brasileira-de-normas-tecnicas.csl` | [citation-style-language/styles](https://github.com/citation-style-language/styles/blob/c7de5bedd4caf9c2ce193e4ccc536f8b5899707b/associacao-brasileira-de-normas-tecnicas.csl) | `c7de5bed` (2025-05-15) | `4f3f2ead…e2163eaae` | CC BY-SA 3.0 (aviso no próprio arquivo) |

Fixado por commit, não pelo `master`: um estilo que muda depois do texto
pronto muda a lista de referências sem ninguém editar nada.

Uso com Quarto ou pandoc — **nenhum dos dois está instalado ainda** (ver a
recomendação de Quarto no relatório de pesquisa), e o texto precisa citar por
chave (`[@chave]`) com um `.bib` ao lado. Hoje o `metodos.md` do pipeline e o
de `/metodos` trazem as referências por extenso, com DOI; o `.bib` sai desses
DOIs por content negotiation do Crossref (`habanero.cn.content_negotiation(ids=doi,
format="bibtex")`, grupo `pesquisa` do pyproject):

```bash
pandoc metodos.md --citeproc --bibliography refs.bib --csl docs/csl/associacao-brasileira-de-normas-tecnicas.csl -o metodos.docx
```

```yaml
# _quarto.yml
csl: docs/csl/associacao-brasileira-de-normas-tecnicas.csl
```

Para conferir uma lista de DOIs antes de citar, o grupo `pesquisa` do
`pyproject.toml` traz `habanero` (Crossref) e `pyalex` (OpenAlex).
