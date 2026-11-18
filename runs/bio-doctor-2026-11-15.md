# Bio-doctor — 15/11/2026

**Resultado: base computacional operacional; integração MCP no OpenCode ausente e R fora do PATH.**

Nenhuma instalação, restauração de ambiente ou correção foi aplicada. Os comandos de correção abaixo são recomendações, não comandos executados.

A distribuição realmente instalada é `Ubuntu-26.04`, e a GPU é uma **NVIDIA GeForce RTX 5060**, não a RTX 5070 Laptop citada na documentação. A sonda NVIDIA foi registrada em `2026-11-16T01:13:08Z` (15/11 no horário local).

## Relatório por impacto

| Componente | Estado | O que fazer |
|---|---|---|
| MCP `nvidia-bio` no OpenCode | **BLOQUEADO no cliente.** `opencode mcp list`: `No MCP servers configured`. O servidor funciona fora do cliente: handshake MCP e catálogo com 14 ferramentas passaram. | Cadastrar, após autorização: `opencode mcp add nvidia-bio -- uv run --script C:/Users/user/bio-stack/servers/nvidia_bio_mcp/server.py`. A chave existente é herdada pelo processo; não colocá-la na linha de comando. |
| MCP `gget` no OpenCode | **BLOQUEADO no cliente**, pela mesma ausência de cadastro. Servidor isolado passou no handshake MCP e expôs 19 ferramentas. | Cadastrar, após autorização: `opencode mcp add gget -- uvx --from gget-mcp gget-mcp server`. Depois conferir ambos com `opencode mcp list`. |
| R / Rscript no PowerShell | **AJUSTE DE PATH.** R 4.6.2 instalado em `C:\Program Files\R\R-4.6.2`; executáveis não são resolvidos pelo PATH. Não é necessário reinstalar R. | Para corrigir apenas a sessão PowerShell: `$env:Path += ';C:\Program Files\R\R-4.6.2\bin'`. Conferir com `Rscript --version` e `R.exe --version`; `R` pode ser um alias do PowerShell. |
| NIM DiffDock / MolMIM | **INCONCLUSIVO.** Ambos devolveram HTTP 500 ao POST com `{}`. Isso não demonstra falha com entradas válidas. | Validar com um payload válido quando houver uma tarefa real para esses modelos. Nenhuma correção local está indicada por esta sonda; testes com payload real não foram disparados. |
| NIM ESMFold / AlphaFold2 / AlphaFold2-Multimer | **APOSENTADOS — HTTP 410**, conforme já marcado no registro do projeto. Não são slugs novos nem falha da chave. | Não selecionar esses endpoints hospedados. Consultar `nim_list_models(task="folding")` após o cadastro do MCP; o registro oferece Boltz2 e OpenFold2/3. Não instalar esses modelos na GPU de 8 GB. |
| Referências de hardware e WSL | **DOCUMENTAÇÃO DESATUALIZADA.** `Ubuntu-24.04` não está instalada. `BIO_STACK_WSL_DISTRO=Ubuntu-26.04` já está correto; os wrappers e o teste de filogenia usam essa distribuição. | Não reinstalar WSL nem reaplicar bootstrap por esse motivo. Usar o prefixo medido, por exemplo: `wsl -d Ubuntu-26.04 -- bash -lc "samtools --version"`. Revisar referências documentais em uma alteração separada. |
| GPU / CUDA / PyTorch | **OK.** RTX 5060, 8.151 MiB de VRAM total, 1.261 MiB ocupados na sonda; compute capability 12.0; driver 621.24. Torch `2.11.0+cu128`, runtime CUDA 12.8 e `sm_120` presentes. | Nenhuma ação. `cuda.is_available()` passou e um produto de matrizes na GPU produziu o resultado esperado (`4096.0`). |
| Windows: ferramentas de base | **OK.** Git 2.56.1; Node 26.13.1; npm 11.20.0; uv/uvx 0.13.7; Python global 3.15.2. | Nenhuma instalação. Para análises, usar o Python 3.11 do projeto, não o Python global 3.15. |
| Python do projeto | **OK.** `.venv\Scripts\python.exe`: Python 3.11.18. `uv pip check --python .venv\Scripts\python.exe`: 282 pacotes compatíveis. | Nenhuma ação. |
| Pacotes Python centrais | **OK.** Biopython 1.88; Scanpy 1.11.5; PyDESeq2 0.5.4; RDKit 2026.3.6; NumPy 2.4.6; SciPy 1.17.1; statsmodels 0.15.0. Imports de `Bio`, `scanpy`, `pydeseq2`, `rdkit`, `anndata`, `scipy` e `statsmodels` passaram. | Nenhuma ação. |
| `gget` / `mcp` ausentes no venv principal | **NÃO É LACUNA OBRIGATÓRIA.** Os servidores usam ambientes isolados. Em ambiente já existente: gget 0.30.8 e gget-mcp 0.1.3 importaram; os dois servidores passaram no protocolo MCP. | Não instalar essas dependências no venv principal apenas para eliminar um aviso de presença. O ajuste necessário é o cadastro no OpenCode, acima. |
| Pacotes adicionais do venv | **METADADOS PRESENTES.** MDAnalysis 2.10.0; MuJoCo 3.15.0; PyVista 0.49.1; Open3D 0.20.0; gmsh 4.15.2; morphosource 1.2.0. | Nenhuma ação por ausência. Imports, aplicações gráficas e cálculos desses pacotes não integraram este teste. |
| Chaves NVIDIA / NCBI | **DEFINIDAS** no ambiente do processo e no escopo User. Valores e prefixos não foram publicados. A sonda NVIDIA não retornou 401/403. | Nenhuma ação de configuração. A presença da chave NCBI foi verificada, mas sua autenticação remota não foi testada. Ausência no escopo Machine não é defeito. |
| WSL / Pixi | **OK.** WSL2 `Ubuntu-26.04`, systemd ativo; Pixi 0.81.0; 54 ambientes instalados. O manifesto instalado é idêntico a `bootstrap/pixi-global.toml`; os 467 comandos expostos foram resolvidos no PATH de login. | Nenhuma ação. Resolução no PATH não equivale a executar todos os 467 binários; os executáveis centrais também foram exercitados separadamente. |
| HTS / CLIs centrais | **OK.** samtools, bcftools e HTSlib 1.24; samtools/bcftools vêm de `~/.pixi/bin`. IQ-TREE3, MAFFT, HyPhy, codeml, BLAST, DIAMOND, miniprot, Datasets e nf-test foram encontrados. | Nenhuma ação. Manter `bash -lc` / `bash -l` nas chamadas WSL. |
| Nextflow / nf-core / Java | **OK.** Nextflow 26.04.6; nf-core 4.1.0; OpenJDK 21.0.13.0.3. | Nenhuma ação. |
| Docker | **OK.** Docker Engine 29.9.3 respondeu dentro do WSL, com 12 CPUs e aproximadamente 15,6 GiB de RAM disponíveis ao daemon. O teste Nextflow/Docker concluiu. | Nenhuma ação. Docker no PATH do Windows não é requisito. |
| Filogenia: lint | **OK.** `nextflow lint main.nf nextflow.config conf modules/local`: 15 arquivos sem erros. | Nenhuma ação. |
| Filogenia: teste ponta a ponta | **OK.** Execução de 12 s, 16 verificações aprovadas: 5 locos de entrada, 4 inferidos, descarte esperado de `loco_E`; concatenada e árvore de espécies com RF `0/18`; gCF 97,22 e sCF 64,13; Métodos sem pendências. | Nenhuma ação. O teste realmente executou; não foi um caso de `PULADO` com saída 0. |
| Biblioteca R / renv | **OK.** Os 280 pacotes do lock foram encontrados nas bibliotecas do projeto e do R, sem divergências de versão. Imports de renv, DESeq2, edgeR, limma, ape, phytools e ggtree passaram. | Nenhuma restauração de pacotes. O runtime é R 4.6.2, enquanto `r/renv.lock` registra R 4.6.1; registrar essa diferença nos Métodos. Não regravar o lock automaticamente. |
| Disco / memória | **OK.** C: com 203,93 GiB livres; WSL com aproximadamente 15 GiB de RAM e 4 GiB de swap. Docker reportou 5 imagens, aproximadamente 1,03 GB. | Planejar downloads pelo espaço físico de C:, não pelos 936 GB virtuais livres reportados dentro do VHD do WSL. |

## Sonda NVIDIA: resultado por rota

`nim_config` foi chamado antes de `nim_probe_endpoints`, usando as funções existentes de `servers/nvidia_bio_mcp/server.py`. Modo: `hosted`; base: `https://health.api.nvidia.com/v1`; chave presente; timeout configurado: 900 s.

| Componente | Estado | O que fazer |
|---|---|---|
| Boltz2 | HTTP 422 — rota alcançada, corpo vazio rejeitado | Usar o payload próprio do modelo. Inferência com entrada válida não foi testada. |
| OpenFold2 | HTTP 422 — rota alcançada | Usar o payload próprio, incluindo MSA/template quando requerido. |
| OpenFold3 | HTTP 422 — rota alcançada | Usar o payload próprio do modelo. |
| MSA-Search | HTTP 422 — rota alcançada | Usar o payload próprio do modelo. |
| RFdiffusion | HTTP 422 — rota alcançada | Usar o payload próprio do modelo. |
| ProteinMPNN | HTTP 422 — rota alcançada | Usar o payload próprio do modelo. |
| Evo2 forward | HTTP 422 — rota alcançada | Usar o payload próprio do modelo. |
| Evo2 generate | HTTP 422 — rota alcançada | Usar o payload próprio do modelo. |
| GenMol | HTTP 200 — aceitou o corpo vazio | Nenhum ajuste de rota indicado. |
| DiffDock | HTTP 500 — inconclusivo com `{}` | Testar com entrada válida em uma tarefa real. |
| MolMIM | HTTP 500 — inconclusivo com `{}` | Testar com entrada válida em uma tarefa real. |
| ESMFold | HTTP 410 — aposentado | Usar um modelo de dobramento ativo. |
| AlphaFold2 | HTTP 410 — aposentado | Usar um modelo de dobramento ativo. |
| AlphaFold2-Multimer | HTTP 410 — aposentado | Usar Boltz2 ou OpenFold3 para complexos. |

**Totais:** 9 rotas respondentes, 2 inconclusivas, 3 aposentadas. Nenhum HTTP 404, 401 ou 403. Responder a uma sonda de rota não comprova uma inferência completa nem a disponibilidade de todo o catálogo NVIDIA.

## Evidências e limites

- Sonda NIM, com chave/prefixo omitidos: `C:\Users\user\AppData\Local\Temp\opencode\bio_doctor_nim_20261115.json`.
- Handshakes e catálogos MCP diretos: `C:\Users\user\AppData\Local\Temp\opencode\bio_doctor_mcp_20261115.json`. NVIDIA: 14 ferramentas; gget: 19 ferramentas; protocolo negociado `2025-11-25`.
- Resultado da filogenia: `C:\Users\user\AppData\Local\Temp\opencode\filogenia_check_47051gec\out`.
- Teste de filogenia disparado com `uv run --offline --no-sync --script tests/check_filogenia.py`, usando a distribuição configurada `Ubuntu-26.04`. O uv informou que `--no-sync` não se aplica a scripts com metadados inline; o teste terminou sem mensagem de instalação de dependências.
- As sondas MCP diretas usaram interpretadores isolados já instalados; não configuraram servidores no OpenCode e não repetiram a sonda remota.
- A configuração existente `.mcp.json` não estava registrada como servidores do OpenCode. Os comandos propostos usam a interface documentada em <https://opencode.ai/docs/mcp-servers>.
- Arquivos científicos brutos, configurações de ambiente, lockfiles e alterações de trabalho preexistentes não foram editados por este diagnóstico. Este relatório é o único arquivo acrescentado ao repositório.
