process IQTREE_CONCORDANCIA {
    tag "${meta.id}"
    label 'process_medium'

    // Mesma imagem do modulo IQTREE (patch local para 3.1.4): concordancia
    // calculada com outra versao do IQ-TREE nao seria comparavel com a arvore.
    container 'quay.io/biocontainers/iqtree:3.1.4--h8471819_0'

    input:
    tuple val(meta), path(arvore_ref), path(arvores_genes, stageAs: 'genes/*'), path(supermatriz), path(modelo_particoes)

    output:
    tuple val(meta), path("${prefix}.gcf.cf.tree"), path("${prefix}.scfl.cf.tree"), emit: arvores
    path "${prefix}.{gcf,scfl}.cf.stat", emit: stats
    tuple val(meta), path("${prefix}.{gcf,scfl}.log"), emit: log
    tuple val("${task.process}"), val('iqtree'), eval("iqtree -version 2>&1 | grep -oE '[0-9]+[.][0-9]+[.][0-9]+' | head -n1"), topic: versions, emit: versions_iqtree

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def quartetos = task.ext.scfl ?: 100
    """
    # Ordem fixa das arvores de genes: o gCF nao depende dela, mas o arquivo
    # publicado fica identico entre execucoes.
    cat ${arvores_genes.join(' ')} > arvores_genes.tre

    # gCF: fracao das arvores de genes que contem cada ramo da arvore de referencia.
    iqtree \\
        -t ${arvore_ref} \\
        --gcf arvores_genes.tre \\
        --prefix ${prefix}.gcf \\
        -T 1 \\
        ${args}

    # sCFL (Mo et al. 2023): concordancia de sitios por verossimilhanca, com o
    # modelo por particao que a analise concatenada ja escolheu (-p best_model.nex),
    # em vez de rodar o ModelFinder de novo.
    iqtree \\
        -te ${arvore_ref} \\
        -s ${supermatriz} \\
        -p ${modelo_particoes} \\
        --scfl ${quartetos} \\
        --prefix ${prefix}.scfl \\
        -T ${task.cpus} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gcf.cf.tree ${prefix}.scfl.cf.tree ${prefix}.gcf.cf.stat ${prefix}.scfl.cf.stat
    touch ${prefix}.gcf.log ${prefix}.scfl.log
    """
}
