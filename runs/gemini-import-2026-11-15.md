# Importação Gemini → OpenCode

Validação final: 2026-11-15, 22:57, UTC−03:00.

## Resultado

Configuração portátil de `C:/Users/user/.gemini` importada globalmente para
OpenCode 2.3.0. A distribuição adotada é **Ubuntu-26.04**, já instalada.

| Componente importado | Resultado verificado |
| --- | --- |
| MCPs | 14 cadastrados: 13 conectados; fal aguarda autenticação |
| Agentes | 11 disponíveis pela API |
| Skills | 21 IDs importados disponíveis pela API |
| Comandos | 21 comandos disponíveis pela API |
| Plugin bio-stack | `gemini-bio-stack` ativo |
| Plugins de autenticação existentes | 3 preservados e ativos |

Configuração principal:
`C:/Users/user/.config/opencode/opencode.jsonc`.

## Arquivos e backup

- Backup anterior à migração:
  `C:/Users/user/.config/opencode/backups/gemini-import-20261116T013751Z/`.
  Contém os originais `opencode.json` e `opencode.jsonc`.
- Manifesto:
  `C:/Users/user/.config/opencode/imports/gemini/manifest.json`.
- Skills e regras:
  `C:/Users/user/.config/opencode/imports/gemini/skills/` e `rules/`.
- Agentes convertidos:
  `C:/Users/user/.config/opencode/agents/`.
- Adaptador nativo de hooks:
  `C:/Users/user/.config/opencode/imports/gemini/bio-stack-plugin/`.
- Ponte de lint:
  `C:/Users/user/.config/opencode/imports/gemini/lint_current.py`.
- Inventário preparado:
  `C:/Users/user/.config/opencode/imports/gemini/config.fragment.json`.
  É um fragmento, não um substituto completo da configuração ativa.

O arquivo ativo `opencode.json` permaneceu idêntico ao backup, confirmado por
hash. Os hashes dos 11 agentes de origem também permanecem iguais aos
registrados no manifesto. A origem Gemini não foi editada.

## MCPs

Conectados na verificação final:

`biocontext`, `biomcp`, `blender`, `chimerax`, `claude-science`,
`gemini-notebook-mcp`, `gget`, `nexus`, `nvidia-bio`, `nvidia-nemo`,
`slicer`, `tooluniverse` e `totalsegmentator`.

O estado conectado confirma a conexão MCP; não significa que todos os
aplicativos gráficos ou serviços externos tenham sido exercitados.

`nvidia-bio.nim_config` foi chamado pela conexão real do OpenCode e confirmou:
modo hospedado, chave presente, saída em `runs/nim/` e proveniência em
`runs/provenance.jsonl`, sem exibir credenciais. Não foram executadas novas
predições ou análises biológicas nesta migração.

### Pendência: fal

A configuração de origem usava `FAL_KEY`, ausente no ambiente. Foi adotado
o endpoint OAuth `https://mcp.fal.ai/mcp-relay`, cujo desafio de autenticação
foi confirmado. O estado final é `needs authentication`.

**No OpenCode, abra `/mcps`, selecione `fal` e conclua o login no navegador.**
Nenhum login foi feito em nome do usuário.

## Compatibilidade e segurança

- Os 13 MCPs locais recebem `BIO_STACK_WSL_DISTRO=Ubuntu-26.04`;
  `claude-science` também recebe `CLAUDE_SCIENCE_WSL_DISTRO`.
- Agentes, regras e skills copiados foram normalizados para a distro atual.
  Nextflow e demais CLIs Linux usam WSL com shell de login.
- Os oito agentes bio e `mod-orchestrator` usam `mode: all`; os dois
  mod-workers usam `mode: subagent`.
- `swarm` foi substituído pela delegação nativa do OpenCode, com profundidade
  máxima configurada em 2. O runtime Gemini não foi carregado como plugin.
- O rótulo de modelo `Orchestrator (Argon)` foi apenas registrado no manifesto:
  não houve troca do modelo atual nem tentativa de inventar um ID equivalente.
- Workspaces confiáveis foram registrados, sem convertê-los em permissões
  irrestritas.
- Providers e os plugins existentes `m365-copilot`,
  `cortexkit.antigravity-auth` e `cortexkit-openai-auth` foram preservados.
- A credencial de API MCP importada é referenciada por arquivo com ACL
  restritiva, não exibida nem escrita inline na configuração.
- Conversas, histórico, caches, bancos e tokens de login Gemini não foram
  importados.
- O toolkit universal-modder foi registrado como referência para manter acesso
  a seus playbooks e recursos locais.

## Hooks e testes

O plugin `gemini-bio-stack` adapta os hooks existentes para a API nativa V2:
contexto, proteção de dados brutos, proveniência, feedback de lint e marcadores
de fim de turno. Esses hooks ficam restritos ao diretório de bio-stack e seus
subdiretórios, embora o plugin esteja registrado globalmente.

Verificações concluídas:

1. Ruff: ambos os scripts de migração e lint passaram.
2. Basedpyright: 0 erros, 0 avisos, 0 notas.
3. `node --check` do plugin: passou.
4. `test.mjs`: registros V2, escopo de projeto, normalização da distro e
   bloqueios de edição/patch em dados brutos passaram.
5. `lint_current.py` com payload `Edit` para
   `pipelines/filogenia/main.nf`: saída 0, sem diagnóstico. Usa
   `wsl.exe -d Ubuntu-26.04 -- bash -l` e não executa o pipeline.
6. APIs `/api/agent`, `/api/skill` e `/api/command`: todos os IDs importados
   encontrados.
7. API `/api/plugin`: adaptador e três plugins anteriores ativos.
8. `opencode mcp list`: 13 conectados e fal aguardando autenticação.
9. Proveniência real e retorno automático de lint foram observados durante a
   migração.

As alterações preexistentes no repositório foram preservadas. Não houve
commit, push, instalação de distro ou instalação de Docker Desktop.

## Referências

- <https://opencode.ai/docs/config>
- <https://opencode.ai/docs/mcp-servers>
- <https://opencode.ai/docs/build/plugins/server>
