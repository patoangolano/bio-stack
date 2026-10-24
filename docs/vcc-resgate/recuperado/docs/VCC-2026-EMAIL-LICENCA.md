<!-- Rascunho de e-mail sobre licenca dos checkpoints PIE. NAO ENVIADO. -->

# Rascunho — pergunta de licença à Arc Institute

> **Estado: NÃO ENVIADO.** Precisa da autorização do Мухамед. Enviar é ação
> externa e irreversível: a mensagem identifica a equipe e a modalidade de
> inscrição a quem organiza a competição.

**Para:** help@virtualcellchallenge.org
**Cópia:** legal@arcinstitute.org
**Assunto:** Virtual Cell Challenge 2026 — eligibility of PIE checkpoints under the non-commercial clause

---

Dear Virtual Cell Challenge team,

I am competing in the 2026 Virtual Cell Challenge final phase as team
**MATHEUS LIMA**, registered as an individual student entrant at SENAC, a
Brazilian educational institution. I have no commercial affiliation and no
commercial funding for this work.

I would like to confirm whether I may use the pretrained PIE checkpoints
(`arcinstitute/PIE_replogle_xdataset` and `arcinstitute/PIE_replogle_wdataset`)
for my submission. Two points in the licensing are ambiguous to me, and I would
rather ask than assume.

**1. Does participation in the Challenge itself count as a Non-Commercial
Purpose?**

The checkpoint cards permit use "if you are a non-commercial entrant", but that
term is not defined in the model card or, as far as I can find, in the
Challenge rules. Clause 1.4 of the license defines Non-Commercial Purpose and
its second sentence excludes activities "engaged in or sponsored by a
Commercial Entity".

The Challenge awards cash prizes and lists commercial sponsors (NVIDIA, 10x
Genomics, Ultima Genomics). I cannot tell whether that makes my participation
an activity "sponsored by a Commercial Entity" under 1.4, even though I am an
individual student entrant with no commercial interest in the result. A
clarification on how 1.4 applies to entrants would settle it.

**2. Does the FAQ exemption for STATE and the Arc Virtual Cell Atlas extend to
the PIE checkpoints?**

The FAQ grants an exemption that reads, to me, as covering the **code**, which
carries CC BY-NC-SA 4.0 — the same license as the PIE repositories. It does not
mention model weights. Since the checkpoints add the separate
"non-commercial entrant" condition, I am unsure whether the same exemption was
intended to reach them.

**3. Generated outputs.**

If the answer to (1) or (2) permits use, I understand from clause 2.1 that
*generating* output with the model requires a Non-Commercial Purpose, while
distributing output falls under 2.1(c)(ii). I would like to confirm that
submitting generated predictions to the Challenge is within scope.

I am proceeding in the meantime with a model trained from scratch on the
publicly released preprocessed data, so my submission does not depend on the
answer. I am asking because I would rather not use the checkpoints at all than
use them outside their terms.

Thank you for your time, and for releasing the PIE assets publicly — the
preprocessed datasets have been useful regardless of the checkpoint question.

Best regards,
Matheus Lima
Team MATHEUS LIMA — Virtual Cell Challenge 2026

---

## Notas para o Мухамед, não para o e-mail

- **Por que vale perguntar mesmo seguindo sem os checkpoints.** A pergunta
  protege a inscrição: se um dia a Arc pedir verificação — a descrição pública
  diz que eles "reserve the right to request more information for verification"
  —, existir um e-mail nosso perguntando antes vale mais que uma interpretação
  nossa depois.
- **Por que o e-mail diz que não dependemos da resposta.** Porque é verdade:
  as etapas 1 a 4 já rodaram sobre dado preprocessado público, e a etapa 7 está
  fora do caminho crítico por limite de RAM, não por licença. Dizer isso tira
  pressão da resposta e aumenta a chance de vir completa.
- **O que o e-mail não faz.** Não pede permissão especial nem exceção, não
  negocia, e não revela nada de estratégia de modelagem. Só pede leitura de
  duas cláusulas.
- **Risco de enviar:** nenhum identificado além de identificar a equipe, o que
  a inscrição já fez.
- **Risco de não enviar:** usar os checkpoints sem esclarecimento, ou deixar de
  usá-los por excesso de cautela. Hoje não usamos, então o risco imediato é
  zero — a pergunta é sobre manter a porta aberta para a etapa 7.
- **Checar a caixa de entrada da inscrição.** Esclarecimentos formais chegam
  por e-mail a todos os participantes, e as regras prevalecem sobre o FAQ.
