#!/usr/bin/env nextflow
/*
 * bio-stack/filogenia
 *
 * Locos nao alinhados -> alinhamento (MAFFT ou FAMSA) -> limpeza (ClipKIT ou
 * trimAl) -> arvore de ML por loco (IQ-TREE 3, ModelFinder + UFBoot2 +
 * SH-aLRT) -> supermatriz particionada (MFP+MERGE) -> fatores de concordancia
 * (gCF e sCFL) -> arvore de especies (wASTRAL) -> resumo e paragrafo de Metodos.
 *
 * Nao existe pipeline nf-core publicado para esse caminho (genephylomodeler
 * esta em desenvolvimento, multiplesequencealign para no alinhamento), entao
 * o fluxo e local, mas cada ferramenta vem de modulo nf-core sempre que ha um.
 */

include { PREPARA_LOCO                 } from './modules/local/prepara_loco/main'
include { MAFFT_ALIGN                  } from './modules/nf-core/mafft/align/main'
include { FAMSA_ALIGN                  } from './modules/nf-core/famsa/align/main'
include { CLIPKIT                      } from './modules/nf-core/clipkit/main'
include { TRIMAL                       } from './modules/nf-core/trimal/main'
include { IQTREE as IQTREE_GENES       } from './modules/nf-core/iqtree/main'
include { IQTREE as IQTREE_CONCAT      } from './modules/nf-core/iqtree/main'
include { CONCATENA                    } from './modules/local/concatena/main'
include { IQTREE_CONCORDANCIA          } from './modules/local/iqtree_concordancia/main'
include { ASTER_WASTRAL                } from './modules/local/aster_wastral/main'
include { RESUMO                       } from './modules/local/resumo/main'

// Parametros TIPADOS. Sem o tipo, o Nextflow 26.04 entrega o valor da linha de
// comando como String: `--concatenar false` virava "false", que e verdadeiro
// em Groovy, e `--ufboot 0` virava "0" — nenhum dos dois desligava nada
// (medido em 2026-10-11). Com o tipo declarado, "false" vira Boolean e "0"
// vira Integer. Texto opcional usa '' como "nao informado".
params {
    input: String = ''
    outdir: String = 'resultados'

    // alinhamento e limpeza
    alinhador: String = 'mafft'          // mafft | famsa | nenhum
    limpeza: String = 'clipkit'          // clipkit | trimal | nenhuma
    clipkit_modo: String = 'smart-gap'

    // IQ-TREE
    modelo: String = 'MFP'               // MFP = ModelFinder Plus; na supermatriz vira MFP+MERGE (so MFP com --amostragem_concat)
    ufboot: Integer = 1000
    bnni: Boolean = true
    alrt: Integer = 1000
    abayes: Boolean = true               // aBayes: e o suporte que o wASTRAL usa como peso
    amostragem_concat: String = ''       // '' | GENE | GENESITE — ver docs/FILOGENIA.md
    seed: Integer = 12345
    outgroup: String = ''                // um ou mais taxons, separados por virgula
    tipo_seq: String = ''                // DNA | AA ; vazio = detecta
    nome_taxon: String = 'primeiro'      // primeiro | completo | colchete (padrao NCBI)
    min_taxa: Integer = 4
    cpus_concatenado: Integer = 1        // >1 acelera, mas perde reproducao bit a bit

    // etapas multi-loco
    concatenar: Boolean = true
    concordancia: Boolean = true
    scfl: Integer = 100
    arvore_especies: Boolean = true
    min_locos_especies: Integer = 3
}

// Nome do loco a partir do arquivo: so letra, digito e "_". O nome vira prefixo
// de arquivo, charset do NEXUS e chave do resumo; um ponto quebraria os tres, e
// "-" viraria "_" no charset (que nao aceita hifen) e dois locos poderiam
// colidir la sem colidir aqui.
def saneiaLoco(nome) {
    return nome.replaceAll(/[^A-Za-z0-9_]/, '_')
}

// Extensao de FASTA, com ou sem .gz, sem distinguir maiuscula: pasta vinda do
// Windows traz loco.FASTA, e o glob do Java (que distingue) o ignorava calado.
def ehFasta(nome) {
    return nome ==~ /(?i).+\.(fa|fasta|fas|fna|faa|aln|fsa)(\.gz)?/
}

def idLoco(arquivo) {
    return saneiaLoco(arquivo.name.replaceAll(/(?i)\.(fa|fasta|fas|fna|faa|aln|fsa)(\.gz)?$/, ''))
}

// Saida opcional de processo como lista de arquivos; lista vazia quando o
// processo nao rodou, para o RESUMO nao ficar esperando um canal que nunca emite.
def soArquivo(canal) {
    return canal.map { _meta, f -> f }.collect().ifEmpty([])
}

def validaParametros() {
    if (!params.input) {
        error("Informe --input: diretorio com um FASTA por loco, um FASTA unico, um padrao glob, ou CSV com colunas loco,fasta")
    }
    if (!(params.alinhador in ['mafft', 'famsa', 'nenhum'])) {
        error("--alinhador deve ser mafft, famsa ou nenhum (recebi '${params.alinhador}')")
    }
    if (!(params.limpeza in ['clipkit', 'trimal', 'nenhuma'])) {
        error("--limpeza deve ser clipkit, trimal ou nenhuma (recebi '${params.limpeza}')")
    }
    if ((params.min_taxa as Integer) < 4) {
        error("--min_taxa precisa ser >= 4: com menos de 4 sequencias nao ha ramo interno e o UFBoot nao faz sentido (o IQ-TREE so avisa)")
    }
    if (params.tipo_seq && !(params.tipo_seq in ['DNA', 'AA'])) {
        error("--tipo_seq deve ser DNA ou AA (ou ficar vazio para detectar)")
    }
    if (!(params.nome_taxon in ['primeiro', 'completo', 'colchete'])) {
        error("--nome_taxon deve ser primeiro, completo ou colchete")
    }
    if (params.cpus_concatenado < 1) {
        error("--cpus_concatenado deve ser >= 1")
    }
    // SITE nao e opcao: e o padrao quando --sampling nao e passado. Passar
    // --sampling SITE explicitamente muda a verossimilhanca das particoes e,
    // com --bnni, aborta o IQ-TREE 3.1.4 ("Wrong -bsam") — medido.
    if (params.amostragem_concat && !(params.amostragem_concat in ['GENE', 'GENESITE'])) {
        error("--amostragem_concat deve ser GENE ou GENESITE (ou ficar vazio)")
    }
    if (params.ufboot && params.ufboot < 1000) {
        error("--ufboot precisa ser 0 (desliga) ou >= 1000: o UFBoot recusa menos replicas")
    }
}

// Rotulo interno que o IQ-TREE escreve: "SH-aLRT/aBayes/UFBoot", so os campos
// pedidos, nessa ordem. O wASTRAL pondera por um campo so; aBayes e o que o
// artigo do wASTRAL mostrou render mais (preset -B, escala 1/3 a 1). Sem
// aBayes, UFBoot ou SH-aLRT na escala 0-100 (preset -S). Sem suporte, so
// comprimento de ramo (--mode 3).
def rotuloSuporte() {
    def ordem = []
    if (params.alrt) {
        ordem.add('alrt')
    }
    if (params.abayes) {
        ordem.add('abayes')
    }
    if (params.ufboot) {
        ordem.add('ufboot')
    }
    def escolha = params.abayes ? 'abayes' : (params.ufboot ? 'ufboot' : (params.alrt ? 'alrt' : null))
    def modo = escolha == 'abayes' ? '-B' : (escolha ? '-S' : '--mode 3')
    return [n: ordem.size(), usar: escolha ? ordem.indexOf(escolha) + 1 : 1, modo: modo]
}

def canalEntrada(entrada) {
    if (entrada.endsWith('.csv')) {
        // Caminho relativo no CSV e relativo ao CSV, nao ao diretorio de
        // lancamento (que, pelo invólucro, e ~/nf-work/filogenia/<id>).
        def base = file(entrada).parent
        return channel.fromPath(entrada, checkIfExists: true)
            .splitCsv(header: true, strip: true)
            .map { linha ->
                if (!linha.loco || !linha.fasta) {
                    error("O CSV de entrada precisa das colunas 'loco' e 'fasta'")
                }
                // file() nao expande o til (medido: file('~/x') vira
                // <lancamento>/~/x), entao o HOME entra aqui.
                def fasta = linha.fasta.startsWith('~/') ? System.getenv('HOME') + linha.fasta.substring(1) : linha.fasta
                def caminho = fasta.startsWith('/') ? file(fasta) : base.resolve(fasta)
                if (!caminho.exists()) {
                    error("CSV: arquivo do loco '${linha.loco}' nao existe: ${caminho} (use caminho do WSL, /mnt/c/...)")
                }
                [[id: saneiaLoco(linha.loco)], caminho]
            }
    }
    def pasta = file(entrada)
    if (!entrada.contains('*') && pasta.isDirectory()) {
        // Lista a pasta (sem descer em subpasta) e filtra pela extensao sem
        // distinguir maiuscula. Arquivo visivel com extensao que nao e de FASTA e
        // avisado; oculto (.*) e subpasta ficam de fora calados — sao o lixo de
        // sempre (.DS_Store, ._loco.fa do macOS, .gitkeep).
        def todos = pasta.listFiles().findAll { f -> f.isFile() && !f.name.startsWith('.') }
        def aceitos = todos.findAll { f -> ehFasta(f.name) }.sort { f -> f.name }
        def ignorados = (todos - aceitos).collect { f -> f.name }.sort()
        if (ignorados) {
            log.warn("Fora da entrada (extensao que nao e de FASTA): ${ignorados.join(', ')}")
        }
        return channel.fromList(aceitos.collect { f -> [[id: idLoco(f)], f] })
    }
    return channel.fromPath(entrada, checkIfExists: true).map { f -> [[id: idLoco(f)], f] }
}

workflow {
    main:
    validaParametros()

    def outgroup = params.outgroup ?: ''
    def util = file("${projectDir}/bin/filogenia_util.py", checkIfExists: true)
    def ordemPorId = { a, b -> a[0].id <=> b[0].id }

    //
    // 1. Entrada: saneia nomes, valida alfabeto, decide quem segue
    //
    // Dois arquivos que saneiam para o mesmo id (gene.1.fa e gene_1.fa, ou duas
    // linhas do CSV) so quebrariam la no CONCATENA, horas depois, com o
    // publishDir ja sobrescrito. Falha aqui, antes de qualquer processo.
    def ch_entrada = canalEntrada(params.input)
        .toList()
        .flatMap { lista ->
            def ids = lista.collect { par -> par[0].id }
            def repetidos = ids.findAll { i -> ids.count(i) > 1 }.unique()
            if (repetidos) {
                def nomes = lista.findAll { par -> par[0].id in repetidos }.collect { par -> par[1].name }
                error("Nomes de loco repetidos depois do saneamento: ${repetidos} (arquivos: ${nomes})")
            }
            if (!lista) {
                error("Nenhum FASTA encontrado em --input ${params.input}")
            }
            return lista
        }
    PREPARA_LOCO(ch_entrada, util, params.min_taxa, outgroup, params.alinhador != 'nenhum')

    // O PREPARA_LOCO altera dados (tira '*' interno, descarta registro vazio), e o
    // stderr de tarefa bem-sucedida o Nextflow nao mostra: o aviso vem do JSON.
    PREPARA_LOCO.out.relatorio.subscribe { f ->
        def rel = new groovy.json.JsonSlurper().parse(f.toFile())
        rel.avisos.each { aviso -> log.warn("${rel.loco}: ${aviso}") }
    }

    def ch_loco = PREPARA_LOCO.out.loco
    ch_loco
        .filter { _meta, _fa, _n, _tipo, status, _og -> status != 'ok' }
        .subscribe { meta, _fa, n, _tipo, status, _og ->
            log.warn("Loco ${meta.id} fica fora da inferencia: ${status} (${n} taxa; minimo ${params.min_taxa})")
        }
    // Outgroup pedido e ausente de todos os locos nao quebra nada — so deixa as
    // arvores sem raiz. Por isso o aviso.
    if (outgroup) {
        ch_loco
            .map { _meta, _fa, _n, _tipo, _status, og -> og }
            .collect()
            .subscribe { ogs ->
                if (ogs.every { og -> !og }) {
                    log.warn("Outgroup '${outgroup}' nao aparece em nenhum loco: arvores sem raiz definida, figuras no ponto medio")
                }
            }
    }
    def ch_ok = ch_loco
        .filter { _meta, _fa, _n, _tipo, status, _og -> status == 'ok' }
        .map { meta, fa, n, tipo, _status, og ->
            [meta + [n_seqs: n as Integer, tipo: tipo, outgroup_presentes: og], fa]
        }
        // Sem nenhum loco valido o resto rodaria vazio e o metodos.md descreveria
        // analises que nao aconteceram.
        // O canal tambem fecha vazio quando o PREPARA_LOCO falha (taxon repetido,
        // caractere invalido) — a mensagem nao pode atribuir isso ao min_taxa.
        .ifEmpty {
            error("Nenhum loco chegou a inferencia. Se o PREPARA_LOCO falhou, a causa e o erro acima; senao, todos os locos ficaram com menos de ${params.min_taxa} taxons (ver ${params.outdir}/01_entrada/*.entrada.json)")
        }

    //
    // 2. Alinhamento. Os tres caminhos sao chamados sempre; o filtro decide
    //    qual recebe dados, e processo sem entrada termina sem rodar.
    //
    MAFFT_ALIGN(
        ch_ok.filter { _m, _f -> params.alinhador == 'mafft' },
        [[:], []], [[:], []], [[:], []], [[:], []], [[:], []],
        false
    )
    FAMSA_ALIGN(ch_ok.filter { _m, _f -> params.alinhador == 'famsa' }, [[:], []], false)
    def ch_alinhados = MAFFT_ALIGN.out.fas
        .mix(FAMSA_ALIGN.out.alignment)
        .mix(ch_ok.filter { _m, _f -> params.alinhador == 'nenhum' })

    //
    // 3. Limpeza do alinhamento
    //
    CLIPKIT(ch_alinhados.filter { _m, _f -> params.limpeza == 'clipkit' }, 'fa', [])
    TRIMAL(ch_alinhados.filter { _m, _f -> params.limpeza == 'trimal' }, 'fa')
    def ch_limpos = CLIPKIT.out.clipkit
        .mix(TRIMAL.out.trimal)
        .mix(ch_alinhados.filter { _m, _f -> params.limpeza == 'nenhuma' })

    //
    // 4. Uma arvore de ML por loco
    //
    IQTREE_GENES(
        ch_limpos.map { meta, aln -> [meta, aln, []] },
        [], [], [], [], [], [], [], [], [], [], [], []
    )

    def ch_genes_lista = IQTREE_GENES.out.phylogeny
        .toSortedList(ordemPorId)
        .map { lista -> lista.collect { par -> par[1] } }

    //
    // 5. Supermatriz particionada, com fusao de particoes pelo ModelFinder
    //
    def ch_para_concatenar = ch_limpos
        .toSortedList(ordemPorId)
        .filter { lista -> params.concatenar && lista.size() >= 2 }
        .map { lista -> [lista.collect { par -> par[0].id }, lista.collect { par -> par[1] }] }
    CONCATENA(ch_para_concatenar, outgroup, util)

    IQTREE_CONCAT(
        CONCATENA.out.supermatriz.map { og, fa, _part -> [[id: 'concatenado', outgroup_presentes: og], fa, []] },
        [], [], [], [],
        CONCATENA.out.supermatriz.map { _og, _fa, part -> part },
        [], [], [], [], [], [], []
    )

    //
    // 6. Fatores de concordancia sobre a arvore concatenada
    //
    // O sCFL reaproveita o modelo por particao da analise concatenada
    // (best_model.nex, que o IQ-TREE 3.1.4 grava em toda analise particionada,
    // inclusive com modelo fixo). Na falta dele tenta o best_scheme.nex e, por
    // fim, as particoes originais, sem modelo — e ai o IQ-TREE seleciona o modelo
    // de novo. So o -stub chega as particoes: ele nao cria nenhum dos dois .nex.
    def ch_concordancia = IQTREE_CONCAT.out.phylogeny
        .join(IQTREE_CONCAT.out.nex)
        .combine(ch_genes_lista.map { arvores -> [arvores] })
        .combine(CONCATENA.out.supermatriz.map { _og, fa, part -> [fa, part] })
        .filter { _meta, _arv, _nexs, arvores, _fa, _part -> params.concordancia && arvores.size() >= 2 }
        .map { meta, arv, nexs, arvores, fa, part ->
            def lista = nexs instanceof List ? nexs : [nexs]
            def modelo = lista.find { f -> f.name.endsWith('.best_model.nex') }
                ?: lista.find { f -> f.name.endsWith('.best_scheme.nex') }
                ?: part
            [meta, arv, arvores, fa, modelo]
        }
    IQTREE_CONCORDANCIA(ch_concordancia)

    //
    // 7. Arvore de especies a partir das arvores de genes
    //
    ASTER_WASTRAL(
        ch_genes_lista
            .filter { arvores -> params.arvore_especies && arvores.size() >= (params.min_locos_especies as Integer) }
            .map { arvores -> [[id: 'especies'], arvores] },
        // Mesmo saneamento que o PREPARA_LOCO aplica aos nomes: "Homo sapiens"
        // na linha de comando vira "Homo_sapiens" nas arvores.
        (outgroup.tokenize(',').find { t -> t.trim() }?.trim() ?: '')
            .replaceAll(/[^A-Za-z0-9_.\-]+/, '_').replaceAll(/_+/, '_').replaceAll(/^_|_$/, ''),
        rotuloSuporte()
    )

    //
    // 8. Versoes, parametros e resumo
    //
    def ch_versoes = channel.topic('versions')
        .unique()
        .map { processo, ferramenta, versao -> [processo.tokenize(':').last(), "  ${ferramenta}: ${versao}"] }
        .groupTuple()
        .map { processo, linhas -> "${processo}:\n" + linhas.unique().sort().join('\n') }
        .collectFile(name: 'versoes.yml', newLine: true, sort: true, storeDir: "${params.outdir}/pipeline_info")

    def chaves = [
        'input', 'alinhador', 'limpeza', 'clipkit_modo', 'modelo', 'ufboot', 'bnni', 'alrt', 'abayes',
        'amostragem_concat', 'seed', 'outgroup', 'min_taxa', 'tipo_seq', 'nome_taxon',
        'cpus_concatenado', 'concatenar', 'concordancia', 'scfl', 'arvore_especies',
        'min_locos_especies',
    ]
    def ch_parametros = channel.of(groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params.subMap(chaves))))
        .collectFile(name: 'parametros.json', storeDir: "${params.outdir}/pipeline_info")

    RESUMO(
        ch_parametros,
        ch_versoes,
        PREPARA_LOCO.out.relatorio.collect().ifEmpty([]),
        soArquivo(ch_alinhados),
        soArquivo(ch_limpos),
        soArquivo(IQTREE_GENES.out.report),
        soArquivo(IQTREE_GENES.out.phylogeny),
        IQTREE_CONCAT.out.report
            .mix(IQTREE_CONCAT.out.phylogeny, IQTREE_CONCAT.out.nex, IQTREE_CONCAT.out.log)
            .map { _meta, f -> f }
            .mix(CONCATENA.out.supermatriz.map { _og, _fa, part -> part }, CONCATENA.out.ocupacao)
            .flatten()
            .collect()
            .ifEmpty([]),
        IQTREE_CONCORDANCIA.out.stats.collect().ifEmpty([]),
        soArquivo(ASTER_WASTRAL.out.arvore.mix(ASTER_WASTRAL.out.log)),
        util
    )

    onComplete:
    if (workflow.success) {
        log.info("Filogenia concluida. Resumo: ${params.outdir}/resumo.tsv | Metodos: ${params.outdir}/metodos.md")
    } else {
        log.error("Filogenia falhou. Veja .nextflow.log e o work/ do processo que caiu.")
    }
}
