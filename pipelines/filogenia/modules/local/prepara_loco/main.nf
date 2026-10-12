process PREPARA_LOCO {
    tag "${meta.id}"
    label 'process_single'

    // Python 3.12 da BioContainers, fixado por digest: o script so usa a
    // biblioteca padrao, entao o container nao precisa de mais nada.
    container 'quay.io/biocontainers/python:3.12.12@sha256:4cc84261c4b7a77a23cfbaa5b5316416a507b42ac8cde330b002631666d99ca2'

    input:
    tuple val(meta), path(fasta, stageAs: 'entrada/*')
    path util
    val min_taxa
    val outgroup
    val desalinhar

    output:
    tuple val(meta), path("${meta.id}.entrada.fa"), env('N_TAXA'), env('TIPO'), env('STATUS'), env('OUTGROUP_PRESENTES'), emit: loco
    path "${meta.id}.entrada.json", emit: relatorio
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // O script entra como arquivo de entrada (nos tres modulos que o usam), e nao
    // pelo nome no PATH do bin/: (1) o git no Windows grava arquivo novo sem bit
    // de execucao, e num clone Linux o nome solto falharia com exit 126; (2)
    // como entrada, o conteudo do script entra no hash da tarefa, e o -resume
    // nao reaproveita resultado de uma versao antiga dele.
    """
    python3 ${util} prepara \\
        --loco ${meta.id} \\
        --fasta ${fasta} \\
        --prefixo ${meta.id}.entrada \\
        --min-taxa ${min_taxa} \\
        --outgroup '${outgroup}' \\
        ${desalinhar ? '--desalinhar' : ''} \\
        ${args}

    # O script devolve as decisoes num arquivo KEY=VALUE; exportar daqui faz
    # cada uma virar saida `env` sem o Groovy precisar abrir JSON.
    set -a
    . ./${meta.id}.entrada.env
    set +a
    """

    stub:
    """
    printf '>a\\nACGT\\n>b\\nACGA\\n>c\\nACTT\\n>d\\nAGTT\\n' > ${meta.id}.entrada.fa
    echo '{}' > ${meta.id}.entrada.json
    N_TAXA=4
    TIPO=DNA
    STATUS=ok
    OUTGROUP_PRESENTES=
    export N_TAXA TIPO STATUS OUTGROUP_PRESENTES
    """
}
