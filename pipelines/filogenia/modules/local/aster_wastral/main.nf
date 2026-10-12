process ASTER_WASTRAL {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/biocontainers/aster:1.25--hc52dbad_1'

    input:
    tuple val(meta), path(arvores_genes, stageAs: 'genes/*')
    val raiz
    // [n: campos no rotulo, usar: campo (1-based) que vira peso, modo: preset do wASTRAL]
    val rotulo

    output:
    tuple val(meta), path("${prefix}.tre"), emit: arvore
    tuple val(meta), path("${prefix}.log"), emit: log
    tuple val(meta), path("${prefix}.entrada.tre"), emit: entrada
    tuple val("${task.process}"), val('aster'), eval("wastral -h 2>&1 | sed -n 's/^Version: v//p' | head -n1"), topic: versions, emit: versions_aster

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // O IQ-TREE escreve um rotulo "SH-aLRT/aBayes/UFBoot" (so os pedidos), e o
    // wASTRAL le em silencio so o PRIMEIRO numero — ponderaria pelo SH-aLRT
    // sem avisar (medido no v1.25). Fica so o campo escolhido.
    def n = rotulo.n as Integer
    def campo = '([0-9.eE+-]+)'
    def padrao = ([campo] * n).join('\\/')
    def reescreve = n > 1 ? "sed -E 's/\\)${padrao}/)\\${rotulo.usar}/g'" : 'cat'
    """
    cat ${arvores_genes.join(' ')} | ${reescreve} > ${prefix}.entrada.tre

    # --root so entra se o taxon existir em alguma arvore de genes: wASTRAL
    # aborta com raiz desconhecida, e enraizar e so apresentacao.
    RAIZ_ARG=""
    if [ -n "${raiz}" ] && grep -qE "[(,]${raiz}[:,)]" ${prefix}.entrada.tre; then
        RAIZ_ARG="--root ${raiz}"
    fi

    wastral \\
        ${rotulo.modo} \\
        -t ${task.cpus} \\
        \$RAIZ_ARG \\
        ${args} \\
        -o ${prefix}.tre \\
        ${prefix}.entrada.tre \\
        2> ${prefix}.log
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '((a,b),(c,d));' > ${prefix}.tre
    touch ${prefix}.log ${prefix}.entrada.tre
    """
}
