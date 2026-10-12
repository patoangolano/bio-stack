process RESUMO {
    label 'process_single'

    container 'quay.io/biocontainers/python:3.12.12@sha256:4cc84261c4b7a77a23cfbaa5b5316416a507b42ac8cde330b002631666d99ca2'

    // Entradas opcionais chegam como lista vazia. Cada grupo vai para um
    // diretorio proprio; o script olha o que existe em cada um, entao etapa
    // que nao rodou (um loco so, sem concatenacao) simplesmente nao aparece.
    input:
    path parametros
    path versoes
    path preparos, stageAs: 'preparos/*'
    path alinhados, stageAs: 'alinhados/*'
    path limpos, stageAs: 'limpos/*'
    path relatorios_genes, stageAs: 'genes/*'
    path arvores_genes, stageAs: 'genes/*'
    path concatenado, stageAs: 'concat/*'
    path cf_stats, stageAs: 'cf/*'
    path especies, stageAs: 'especies/*'
    path util

    output:
    path 'resumo.tsv', emit: tsv
    path 'resumo.json', emit: json
    path 'metodos.md', emit: metodos

    when:
    task.ext.when == null || task.ext.when

    script:
    def cf_args = cf_stats ? "--cf ${cf_stats}" : ''
    """
    python3 ${util} resumo \\
        --params ${parametros} \\
        --versoes ${versoes} \\
        --versao-nextflow ${workflow.nextflow.version} \\
        --preparos ${preparos} \\
        --alinhados ${alinhados} \\
        --limpos ${limpos} \\
        --relatorios-genes ${relatorios_genes} \\
        --dir-concat concat \\
        --dir-especies especies \\
        ${cf_args} \\
        --saida-prefixo resumo
    """

    stub:
    """
    touch resumo.tsv metodos.md
    echo '{}' > resumo.json
    """
}
