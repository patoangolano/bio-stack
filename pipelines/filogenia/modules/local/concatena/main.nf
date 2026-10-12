process CONCATENA {
    tag "${locos.size()} locos"
    label 'process_single'

    container 'quay.io/biocontainers/python:3.12.12@sha256:4cc84261c4b7a77a23cfbaa5b5316416a507b42ac8cde330b002631666d99ca2'

    input:
    tuple val(locos), path(alinhamentos, stageAs: 'alinhamentos/*')
    val outgroup
    path util

    output:
    tuple env('OUTGROUP_PRESENTES'), path('supermatriz.fa'), path('particoes.nex'), emit: supermatriz
    path 'ocupacao.tsv', emit: ocupacao
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // A ordem de `locos` e a de `alinhamentos` vem da mesma lista ordenada no
    // workflow; e ela que amarra cada charset do NEXUS ao arquivo certo.
    """
    python3 ${util} concatena \\
        --locos ${locos.join(' ')} \\
        --alinhamentos ${alinhamentos.join(' ')} \\
        --outgroup '${outgroup}'

    set -a
    . ./concatena.env
    set +a
    """

    stub:
    """
    printf '>a\\nACGTACGT\\n>b\\nACGAACGA\\n>c\\nACTTACTT\\n>d\\nAGTTAGTT\\n' > supermatriz.fa
    printf '#nexus\\nbegin sets;\\n  charset l1 = 1-4;\\n  charset l2 = 5-8;\\nend;\\n' > particoes.nex
    printf 'taxon\\tl1\\tl2\\tn_locos\\n' > ocupacao.tsv
    OUTGROUP_PRESENTES=
    export OUTGROUP_PRESENTES
    """
}
