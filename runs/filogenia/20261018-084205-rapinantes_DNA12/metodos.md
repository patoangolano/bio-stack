# Métodos — filogenia

Os alinhamentos fornecidos foram usados sem realinhamento. Para cada um dos 300 locos, a árvore de máxima verossimilhança foi inferida no IQ-TREE v3.1.4, com o modelo de substituição escolhido pelo ModelFinder (-m MFP) segundo o BIC; o suporte dos ramos foi estimado por SH-aLRT (1000 réplicas), aBayes e UFBoot2 (1000 réplicas, com otimização --bnni). A árvore de espécies foi estimada com o wASTRAL (ASTER v1.25.3.8) a partir das árvores de genes, ponderando os quartetos pelo suporte aBayes e pelo comprimento de ramo (modo híbrido). O fluxo foi executado no Nextflow v26.04.6 com módulos nf-core em containers BioContainers, semente aleatória 12345 e IQ-TREE em uma thread por loco.

## Referências

- Wong TKF, Ly-Trong N, Ren H, et al. (2026) IQ-TREE 3: phylogenomic inference software using complex evolutionary models. Molecular Biology and Evolution 43(5):msag117. doi:10.1093/molbev/msag117
- Di Tommaso P, Chatzou M, Floden EW, et al. (2017) Nextflow enables reproducible computational workflows. Nature Biotechnology 35:316-319. doi:10.1038/nbt.3820
- Ewels PA, Peltzer A, Fillinger S, et al. (2020) The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology 38:276-278. doi:10.1038/s41587-020-0439-x
- Guindon S, Dufayard JF, Lefort V, et al. (2010) New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0. Systematic Biology 59:307-321. doi:10.1093/sysbio/syq010
- Anisimova M, Gil M, Dufayard JF, Dessimoz C, Gascuel O (2011) Survey of branch support methods demonstrates accuracy, power, and robustness of fast likelihood-based approximation schemes. Systematic Biology 60:685-699. doi:10.1093/sysbio/syr041
- Hoang DT, Chernomor O, von Haeseler A, Minh BQ, Vinh LS (2018) UFBoot2: improving the ultrafast bootstrap approximation. Molecular Biology and Evolution 35:518-522. doi:10.1093/molbev/msx281
- Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS (2017) ModelFinder: fast model selection for accurate phylogenetic estimates. Nature Methods 14:587-589. doi:10.1038/nmeth.4285
- Zhang C, Mirarab S (2022) Weighting by gene tree uncertainty improves accuracy of quartet-based species trees. Molecular Biology and Evolution 39:msac215. doi:10.1093/molbev/msac215
- Zhang C, Nielsen R, Mirarab S (2025) ASTER: a package for large-scale phylogenomic reconstructions. Molecular Biology and Evolution 42(8):msaf172. doi:10.1093/molbev/msaf172

## Pendências

- Nenhuma.
