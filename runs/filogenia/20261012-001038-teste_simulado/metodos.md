# Métodos — filogenia

As sequências de cada loco foram alinhadas com o MAFFT v7.520 (L-INS-i até 200 sequências e estratégia automática acima disso; refinamento iterativo em uma thread, para resultado determinístico). Colunas com excesso de gaps foram removidas com o ClipKIT v2.11.4 (modo smart-gap, limiar de gaps determinado dinamicamente). Para cada um dos 4 locos, a árvore de máxima verossimilhança foi inferida no IQ-TREE v3.1.4, com o modelo de substituição escolhido pelo ModelFinder (-m MFP) segundo o BIC; o suporte dos ramos foi estimado por SH-aLRT (1000 réplicas), aBayes e UFBoot2 (1000 réplicas, com otimização --bnni). Os alinhamentos limpos foram concatenados numa supermatriz de 12 táxons e 4055 sítios, particionada por loco (4 partições; 93.8% das células táxon × loco preenchidas), e analisada com modelo de partição de comprimentos de ramo proporcionalmente ligados e fusão de partições pelo ModelFinder (MFP+MERGE, algoritmo rclusterf); o esquema final teve 1 partição, com os mesmos critérios de suporte. Os fatores de concordância de genes (gCF) e de sítios por verossimilhança (sCFL, 100 quartetos por ramo) foram calculados no IQ-TREE sobre a árvore concatenada. A árvore de espécies foi estimada com o wASTRAL (ASTER v1.25.3.8) a partir das árvores de genes, ponderando os quartetos pelo suporte aBayes e pelo comprimento de ramo (modo híbrido). O fluxo foi executado no Nextflow v26.04.6 com módulos nf-core em containers BioContainers, semente aleatória 12345 e IQ-TREE em uma thread por loco (supermatriz: 1).

## Referências

- Wong TKF, Ly-Trong N, Ren H, et al. (2026) IQ-TREE 3: phylogenomic inference software using complex evolutionary models. Molecular Biology and Evolution 43(5):msag117. doi:10.1093/molbev/msag117
- Di Tommaso P, Chatzou M, Floden EW, et al. (2017) Nextflow enables reproducible computational workflows. Nature Biotechnology 35:316-319. doi:10.1038/nbt.3820
- Ewels PA, Peltzer A, Fillinger S, et al. (2020) The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology 38:276-278. doi:10.1038/s41587-020-0439-x
- Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution 30:772-780. doi:10.1093/molbev/mst010
- Steenwyk JL, Buida TJ, Li Y, Shen XX, Rokas A (2020) ClipKIT: a multiple sequence alignment trimming software for accurate phylogenomic inference. PLoS Biology 18:e3001007. doi:10.1371/journal.pbio.3001007
- Guindon S, Dufayard JF, Lefort V, et al. (2010) New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. Systematic Biology 59:307-321. doi:10.1093/sysbio/syq010
- Anisimova M, Gil M, Dufayard JF, Dessimoz C, Gascuel O (2011) Survey of branch support methods demonstrates accuracy, power, and robustness of fast likelihood-based approximation schemes. Systematic Biology 60:685-699. doi:10.1093/sysbio/syr041
- Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS (2018) UFBoot2: improving the ultrafast bootstrap approximation. Molecular Biology and Evolution 35:518-522. doi:10.1093/molbev/msx281
- Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS (2017) ModelFinder: fast model selection for accurate phylogenetic estimates. Nature Methods 14:587-589. doi:10.1038/nmeth.4285
- Lanfear R, Calcott B, Kainer D, Mayer C, Stamatakis A (2014) Selecting optimal partitioning schemes for phylogenomic datasets. BMC Evolutionary Biology 14:82. doi:10.1186/1471-2148-14-82
- Lanfear R, Frandsen PB, Wright AM, Senfeld T, Calcott B (2017) PartitionFinder 2: new methods for selecting partitioned models of evolution for molecular and morphological phylogenetic analyses. Molecular Biology and Evolution 34:772-773. doi:10.1093/molbev/msw260
- Chernomor O, von Haeseler A, Minh BQ (2016) Terrace aware data structure for phylogenomic inference from supermatrices. Systematic Biology 65:997-1008. doi:10.1093/sysbio/syw037
- Minh BQ, Hahn MW, Lanfear R (2020) New methods to calculate concordance factors for phylogenomic datasets. Molecular Biology and Evolution 37:2727-2733. doi:10.1093/molbev/msaa106
- Mo YK, Lanfear R, Hahn MW, Minh BQ (2023) Updated site concordance factors minimize effects of homoplasy and taxon sampling. Bioinformatics 39(1):btac741. doi:10.1093/bioinformatics/btac741
- Zhang C, Mirarab S (2022) Weighting by gene tree uncertainty improves accuracy of quartet-based species trees. Molecular Biology and Evolution 39:msac215. doi:10.1093/molbev/msac215
- Zhang C, Nielsen R, Mirarab S (2025) ASTER: a package for large-scale phylogenomic reconstructions. Molecular Biology and Evolution 42(8):msaf172. doi:10.1093/molbev/msaf172

## Pendências

- Nenhuma.
