# Transporte Colab -> bio-stack (medido em 2026-10-10)

O MCP do Colab roda numa sessao de navegador NA MAQUINA do usuario. Logo
`google.colab.files.download()` dispara um download real do navegador e o
arquivo cai em `C:\Users\Manec\Downloads`.

Verificado ponta a ponta com um PDB sintetico de 450 residuos:

    Colab   /content/teste_transporte.pdb   248.854 bytes
            sha256 06312a10ef4a7aec30b0ca143fc063e587d08496a4882e976fbcd105e412f36c
    Local   Downloads/teste_transporte.pdb  248.854 bytes
            sha256 06312a10ef4a7aec30b0ca143fc063e587d08496a4882e976fbcd105e412f36c

Integro. Nenhum byte da estrutura passou pelo contexto do agente, o que
atende a regra do CLAUDE.md ("coordenadas atomicas nunca entram no contexto";
as ferramentas devolvem so o caminho e as metricas).

## Alternativas descartadas e por que

- **base64 inline na saida da celula.** Funciona (450 residuos = 28.619 bytes
  em gzip, 38.160 chars em base64) mas a estrutura inteira entraria no
  contexto do agente a cada dobramento. Medido, rejeitado por politica.
- **Google Drive montado.** Nao ha Drive sincronizado nesta maquina, nem
  rclone nem gcloud instalados. Exigiria setup novo para ganhar nada que o
  download ja da.

## Contrato do lado do Colab

A celula grava o arquivo, imprime caminho, bytes e sha256, e chama
`files.download()`. O lado local confere o sha256 antes de mover de
`Downloads` para `runs/colab/`. Sha diferente = download truncado, e o
arquivo e descartado em vez de entrar na analise.
