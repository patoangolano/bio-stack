-- Expressao dos 300 alvos de CRISPRi em cada contexto do VCC 2026.
--
-- CRISPRi reprime transcricao. Alvo que o contexto nao transcreve nao tem o
-- que derrubar: a perturbacao tende ao controle, e o modelo nao deve inventar
-- efeito ali. Entrada: saida de scripts/vcc/caracterizar_contextos.py.
--
-- Rodar no WSL (duckdb vem do pixi, exige shell de login):
--   wsl -d Ubuntu-24.04 -- bash -lc "duckdb < /mnt/c/Users/Manec/bio-stack/scripts/vcc/alvos_por_contexto.sql"

CREATE TABLE alvos AS
SELECT *
FROM read_csv('/mnt/c/Users/Manec/bio-stack/runs/vcc/caracterizacao/alvos_por_contexto.tsv',
              delim = '\t', header = true);

.print '== 1. Alvos por faixa de deteccao (fracao de celulas com >= 1 UMI) =='
SELECT contexto,
       count(*) FILTER (WHERE fracao_expressa >= 0.5)                            AS ubiquo_ge50,
       count(*) FILTER (WHERE fracao_expressa >= 0.1  AND fracao_expressa < 0.5) AS parcial_10a50,
       count(*) FILTER (WHERE fracao_expressa >= 0.02 AND fracao_expressa < 0.1) AS raro_2a10,
       count(*) FILTER (WHERE fracao_expressa < 0.02)                            AS ausente_lt2
FROM alvos
GROUP BY contexto
ORDER BY contexto;

.print '== 2. Alvos praticamente ausentes nos tres contextos (< 2% das celulas) =='
SELECT alvo, round(max(fracao_expressa), 4) AS fracao_max
FROM alvos
GROUP BY alvo
HAVING max(fracao_expressa) < 0.02
ORDER BY fracao_max;

.print '== 3. Alvos cuja deteccao muda >= 50 pontos entre contextos =='
SELECT count(*) AS alvos_dependentes_de_contexto
FROM (SELECT alvo FROM alvos GROUP BY alvo
      HAVING max(fracao_expressa) - min(fracao_expressa) >= 0.5);

SELECT alvo,
       round(max(fracao_expressa) FILTER (WHERE contexto = 'D'), 3) AS D,
       round(max(fracao_expressa) FILTER (WHERE contexto = 'E'), 3) AS E,
       round(max(fracao_expressa) FILTER (WHERE contexto = 'F'), 3) AS F,
       round(max(fracao_expressa) - min(fracao_expressa), 3)        AS amplitude
FROM alvos
GROUP BY alvo
HAVING max(fracao_expressa) - min(fracao_expressa) >= 0.5
ORDER BY amplitude DESC
LIMIT 20;
