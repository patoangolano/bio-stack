# Métodos — filogenia

As sequências de cada loco foram alinhadas com o FAMSA v2.4.1-45c9b2b. Colunas com excesso de gaps foram removidas com o ClipKIT v2.11.4 (modo smart-gap, limiar de gaps determinado dinamicamente). Para cada um dos 300 locos, a árvore de máxima verossimilhança foi inferida no IQ-TREE v3.1.4, com o modelo de substituição escolhido pelo ModelFinder (-m MFP) segundo o BIC; o suporte dos ramos foi estimado por aBayes. A árvore de espécies foi estimada com o wASTRAL (ASTER v1.25.3.8) a partir das árvores de genes, ponderando os quartetos pelo suporte aBayes e pelo comprimento de ramo (modo híbrido). O fluxo foi executado no Nextflow v26.04.6 com módulos nf-core em containers BioContainers, semente aleatória 12345 e IQ-TREE em uma thread por loco.

## Referências

- Wong TKF, Ly-Trong N, Ren H, et al. (2026) IQ-TREE 3: phylogenomic inference software using complex evolutionary models. Molecular Biology and Evolution 43(5):msag117. doi:10.1093/molbev/msag117
- Di Tommaso P, Chatzou M, Floden EW, et al. (2017) Nextflow enables reproducible computational workflows. Nature Biotechnology 35:316-319. doi:10.1038/nbt.3820
- Ewels PA, Peltzer A, Fillinger S, et al. (2020) The nf-core framework for community-curated bioinformatics pipelines. Nature Biotechnology 38:276-278. doi:10.1038/s41587-020-0439-x
- Deorowicz S, Debudaj-Grabysz A, Gudys A (2016) FAMSA: fast and accurate multiple sequence alignment of huge protein families. Scientific Reports 6:33964. doi:10.1038/srep33964
- Steenwyk JL, Buida TJ, Li Y, Shen XX, Rokas A (2020) ClipKIT: a multiple sequence alignment trimming software for accurate phylogenomic inference. PLoS Biology 18:e3001007. doi:10.1371/journal.pbio.3001007
- Anisimova M, Gil M, Dufayard JF, Dessimoz C, Gascuel O (2011) Survey of branch support methods demonstrates accuracy, power, and robustness of fast likelihood-based approximation schemes. Systematic Biology 60:685-699. doi:10.1093/sysbio/syr041
- Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS (2017) ModelFinder: fast model selection for accurate phylogenetic estimates. Nature Methods 14:587-589. doi:10.1038/nmeth.4285
- Zhang C, Mirarab S (2022) Weighting by gene tree uncertainty improves accuracy of quartet-based species trees. Molecular Biology and Evolution 39:msac215. doi:10.1093/molbev/msac215
- Zhang C, Nielsen R, Mirarab S (2025) ASTER: a package for large-scale phylogenomic reconstructions. Molecular Biology and Evolution 42(8):msaf172. doi:10.1093/molbev/msaf172
- Tan G, Muffato M, Ledergerber C, et al. (2015) Current methods for automated filtering of multiple sequence alignments frequently worsen single-gene phylogenetic inference. Systematic Biology 64:778-791. doi:10.1093/sysbio/syv033

## Pendências

- Locos com mais de 20% dos sítios removidos na limpeza: ACVR2B, AGGF1, ATXN3, BMP2, C6orf62, CD74, CHRNA2, COA7, CTSC, DPP7, GKAP1, GPBP1, GUCD1, HDAC8, INKA2, LRWD1, LUZP2, NKIRAS1, NUBP1, OTX2, PACRGL, PLCXD1, RAB32, RAD51B, SLCO5A1, SNRNP48, SNX10, SNX5, SUMF1, SUSD4, TANGO2, TIMMDC1, TSTD1, TTF1. Limpeza agressiva tende a piorar árvores de gene único (Tan et al. 2015); confira se a conclusão sobrevive a --limpeza nenhuma.
