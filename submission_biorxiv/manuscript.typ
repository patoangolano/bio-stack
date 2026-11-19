// =============================================================================
// bioRxiv Preprint — Pure Typst 0.15 Source (Post-Adversarial Review Final)
// Continuous line numbering + Page numbering + Path-Outlined Vector SVGs + 46/46 BibTeX
// =============================================================================

#set document(
  title: "Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles",
  author: "Matheus Lima",
)

#set page(
  paper: "us-letter",
  margin: (left: 1.0in, right: 0.85in, top: 0.85in, bottom: 0.85in),
  header: context {
    if counter(page).get().first() > 1 {
      set text(font: "Lato", size: 8pt, fill: luma(110))
      grid(
        columns: (1fr, auto),
        [Lima (2026) · Enamel-Gene Relics, IDP Homology Bias, and Pseudogene Dating in Sauropsids],
        [bioRxiv Preprint],
      )
      v(-6pt)
      line(length: 100%, stroke: 0.4pt + luma(200))
    }
  },
  footer: context {
    set align(center)
    set text(font: "Lato", size: 8.5pt, fill: luma(90))
    [Page #counter(page).display("1 of 1", both: true)]
  },
)

#set text(font: "Libertinus Serif", size: 10.2pt, lang: "en", fill: luma(20))
#set par(justify: true, leading: 0.68em, spacing: 0.95em)

// Continuous line numbering across the manuscript (bioRxiv review standard)
#set par.line(
  numbering: n => text(font: "Lato", size: 7.5pt, fill: luma(140))[#n],
  number-margin: left,
  number-clearance: 11pt,
)

#show heading.where(level: 1): it => {
  v(1.1em, weak: true)
  set par.line(numbering: none)
  set text(font: "Lato", size: 12.0pt, weight: "bold", fill: rgb("#112233"))
  block(breakable: false)[
    #it.body
    #v(-5pt)
    #line(length: 100%, stroke: 0.6pt + rgb("#B8C4D0"))
  ]
  v(0.45em, weak: true)
}

#show heading.where(level: 2): it => {
  v(0.85em, weak: true)
  set par.line(numbering: none)
  set text(font: "Lato", size: 10.3pt, weight: "bold", fill: rgb("#1C3550"))
  block(breakable: false, it.body)
  v(0.3em, weak: true)
}

#show figure.where(kind: table): set figure.caption(position: top)

#show figure.caption: it => {
  set par.line(numbering: none)
  set text(font: "Lato", size: 8.5pt, fill: luma(35))
  set par(justify: true, leading: 0.55em)
  block(width: 100%, inset: (top: 3pt, bottom: 4pt))[
    *#it.supplement #context it.counter.display(it.numbering).* #it.body
  ]
}

#show table: set text(font: "Lato", size: 7.9pt)
#show table: set par(justify: false, leading: 0.45em)
#show table: set par.line(numbering: none)

#let mono(s) = text(font: "DejaVu Sans Mono", size: 7.9pt)[#s]
#let tmono(s) = text(font: "DejaVu Sans Mono", size: 7.2pt)[#s]
#let pcite(k) = cite(k, form: "prose")

// =============================================================================
// TITLE BLOCK (Line numbering disabled on front matter)
// =============================================================================
#[
  #set par.line(numbering: none)
  #set par(justify: false)

  #text(font: "Lato", size: 15.2pt, weight: "bold", fill: rgb("#0E2439"))[
    Unannotated Enamel-Gene Relics on the Avian Z Chromosome, Disordered-Protein Homology Bias, and Phylogenetic Limits of Inactivating-Mutation Dating in Birds and Turtles
  ]

  #v(0.7em)

  #text(font: "Lato", size: 10.8pt, weight: "semibold")[
    Matheus Lima#super[1,\*]
  ]

  #v(0.25em)

  #text(font: "Lato", size: 8.9pt, fill: luma(60))[
    #super[1] Independent Researcher, Brazil \
    \* Correspondence: #mono("limatheus458@gmail.com") · ORCID: #link("https://orcid.org/0009-0000-3557-8337")[#mono("0009-0000-3557-8337")] \
    Target bioRxiv Subject Area: *Evolutionary Biology* (secondary: *Genomics*)
  ]

  #v(0.7em)
  #line(length: 100%, stroke: 0.8pt + rgb("#1C3550"))
  #v(0.5em)
]

// =============================================================================
// ABSTRACT (Strictly <= 250 words)
// =============================================================================
#heading(level: 1, numbering: none)[Abstract]

Convergent tooth loss in birds and turtles inactivated enamel matrix protein genes (_ENAM_, _AMBN_, _AMTN_), leaving sequence relics whose detection and dating are confounded by annotation blind spots and structural disorder. We audited 18 sauropsid RefSeq reference assemblies alongside _Gallus gallus_ `GRCg7b` Ensembl Release 116 across the conserved _NRG1–LPL_ enamel cluster. Both NCBI RefSeq Release 106 and Ensembl Release 116 omit all enamel genes in `GRCg7b`, yet synteny-anchored searches recover single-copy, whole-proteome reciprocal-confirmed relics of _ΨENAM_ Exon 8 (`NC_052572.1:55,155,480–55,157,464`, `−`), _ΨENAM_ Exon 6 (`55,160,529–55,160,677`, `−`, bearing a 1-nt deletion and in-frame `TAA` stop), and _ΨAMBN_ Exon 3 (`55,168,855–55,168,992`, `−`) with zero annotation overlap, plus a reciprocal-confirmed _ΨAMTN_ relic in kakapo (`NC_044302.2:51,216,319–51,216,444`, `−`). Conversely, automated RefSeq pipelines misannotate pseudogenized enamel loci as `PROTEIN_CODING` in four edentulous genomes. Equidistant _Anolis_ queries (\~281 Ma) with 20 composition-matched permutations reveal that intact osteogenic intrinsically disordered proteins (_IBSP_, _SPP1_) lose 26–52% of excess `tBLASTn` bitscore in birds—overlapping decayed avian _ENAM_/_AMBN_—whereas globular _LPL_ scores 112% of crocodylian levels, demonstrating that bitscore attenuation alone cannot diagnose pseudogenization. Splice-aware `miniprot` alignment detects >150-fold coding-lesion enrichment in edentulous enamel loci (16.98 vs. 0.11 lesions/1,000 aa in toothed controls; 0.00 across 18/18 _LPL_ controls, 133/133 canonical `GT-AG` introns) and disrupted `GT→AT` splice donors. Dual-aligner (`MAFFT L-INS-i` + `MUSCLE 5`) consensus recovers 23 shared lesions dating Palaeognathae, Novaeratitae, Galloanserae, Psittacopasserae, and Durocryptodira, but 0 crown Neornithes synapomorphies, with 6/16 (37.5%) shared avian _ENAM_ premature stops representing homoplasies.

#[
  #set par.line(numbering: none)
  #v(0.2em)
  #text(font: "Lato", size: 8.6pt, fill: luma(50))[
    *Keywords:* comparative genomics; paleogenomics; edentulism; enamel matrix proteins; intrinsically disordered proteins; pseudogene dating
  ]
]

// =============================================================================
// 1. INTRODUCTION
// =============================================================================
#heading(level: 1, numbering: none)[1. Introduction]

Mineralized teeth have been lost independently across multiple vertebrate lineages, including the stem lineages of crown birds (Aves, \~116 Ma) and crown turtles (Testudines, \~157–210 Ma) as well as several placental mammal clades @davitbeal2009 @meredith2013 @meredith2014 @brocklehurst2021. Because tooth enamel mineralization requires a specialized subset of secretory calcium-binding phosphoprotein (_SCPP_) genes—chiefly enamelin (_ENAM_), ameloblastin (_AMBN_), amelotin (_AMTN_), and amelogenin (_AMELX_)—together with maturation factors such as odontogenesis-associated phosphoprotein (_ODAPH_, formerly _C4orf26_) and matrix metalloproteinase 20 (_MMP20_), relaxed purifying selection following edentulism leaves decaying molecular fossils (pseudogenes) embedded in conserved genomic neighborhoods @kawasaki2003 @kawasaki2009 @meredith2009 @springer2016 @liang2021 @kawasaki2023 @chun2025.

Foundational paleogenomic studies established the core architecture of sauropsid enamel-gene decay. #pcite(<sire2008>) identified a degraded amelogenin pseudogene (_ΨAMEL_, `EU340348`) on chicken chromosome 1 while inferring complete deletion of other enamel loci in early chicken assemblies. Shortly thereafter, #pcite(<alhashimi2010>) demonstrated that in sauropsids the _ENAM–AMBN–AMTN_ enamel subcluster resides in a translocated syntenic neighborhood flanked by _NRG1_ and _LPL_—distinct from the ancestral _SPARCL1–SIBLING_ cluster retained on mammalian chromosome 4—and reported the first chicken _ENAM_ pseudogene fragment in that region. Subsequent comparative analyses documented widespread inactivating mutations in _ENAM_, _AMBN_, _AMTN_, _AMELX_, _MMP20_, and _DSPP_ across turtles @meredith2013 and 48 avian genomes @meredith2014, inferring a single loss of mineralized teeth in the common ancestor of crown birds and documenting parallel regression across toothless mammals @springer2019 @mu2021 @randall2022 @emerling2023 @randall2024 @emerling2026. Recent developmental work has also linked ciliary and signaling regulators such as _C2CD3_ to avian tooth suppression alongside structural enamel-gene decay @kantaputra2026.

Despite this empirical foundation, three pervasive methodological vulnerabilities affect modern paleogenomic surveys of tooth-gene loss—and comparative analyses of trait loss more broadly. First, automated gene-annotation pipelines (NCBI RefSeq and Ensembl) exhibit bidirectional errors on degenerating loci: they routinely omit genuine pseudogene fragments even in modern chromosome-level reference genomes while occasionally projecting erroneous `PROTEIN_CODING` models onto frame-disrupted pseudogenes @weisman2022 @kirilenko2023. Second, enamel matrix proteins are proline- and glutamine-rich intrinsically disordered proteins (IDPs) lacking globular catalytic domains @kawasaki2006 @kawasaki2021 @chun2025. In translated sequence searches (`tBLASTn`), compositional bias, low-complexity filtering (`SEG`), and phylogenetic asymmetry systematically depress alignment bitscores even for intact, fast-evolving IDPs @wootton1993 @wootton1996 @schaffer2001 @yu2003 @yu2005 @altschul2005 @gertz2006 @weisman2020 @chow2024. Treating depressed translated-alignment bitscores as proxies for pseudogenization risks conflating evolutionary disorder with gene inactivation. Third, dating ancient trait losses via shared premature termination codons (PTCs) is vulnerable to mutational convergence (homoplasy) at hypermutable triplets (e.g., `CpG` dinucleotides and `TGG` tryptophan codons) and alignment-dependent gap placement @meredith2014 @emerling2026.

To separate empirical replication from original methodological and genomic findings, we designed a deterministic, four-level audit ("Evidence Ladder") across 18 high-contiguity sauropsid RefSeq reference genomes (6 toothed controls and 12 edentulous birds and turtles, anchored against human `GRCh38.p14`). While replicating and extending prior paleogenomic discoveries of avian and testudine enamel pseudogenes @alhashimi2010 @meredith2013 @meredith2014, this study makes four original contributions:
1. A dual-database annotation audit (NCBI RefSeq Release 106 and Ensembl Release 116) proving that _all_ enamel loci remain completely unannotated in the modern chromosome-level reference _Gallus gallus_ `GRCg7b` while establishing exact, single-copy, whole-proteome reciprocal-confirmed `GRCg7b` coordinates of _ΨENAM_ Exon 8, _ΨENAM_ Exon 6 (bearing a direct frameshift and premature stop), and _ΨAMBN_ Exon 3 on chromosome Z, alongside a reciprocal-confirmed _ΨAMTN_ relic on kakapo (_Strigops habroptila_) chromosome Z.
2. A quantitative 18-species benchmark of compositional and `SEG`-masking bias in translated homology searches, contrasting functional osteogenic IDPs (_IBSP_, _SPP1_, _SPARCL1_) against a globular flanking control (_LPL_) using equidistant _Anolis_ queries (\~281 Ma) and 20 composition-preserving permutations.
3. A splice-aware `miniprot` coding-lesion and `GT-AG` splice-site audit across all 18 _NRG1–LPL_ intervals showing that replacing a truncated legacy GenBank enamelin sequence (`QOX06214.1`, 849 aa) with full-length RefSeq _ENAM_ (`XP_019346377.1`, 1,078 aa) unlocks previously invisible pseudogene exons and canonical `GT→AT` splice-donor disruptions while exposing a 58.4-kb distal low-complexity chaining artifact in _Chrysemys picta_.
4. A dual-aligner (`MAFFT L-INS-i` + `MUSCLE 5`) consensus analysis of shared inactivating mutations and exon splice borders across avian and testudine _ENAM_ and _AMBN_, quantifying a 37.5% (6/16) homoplasy rate among shared avian _ENAM_ PTCs and demonstrating why multi-nucleotide frameshift indels and canonical splice-site disruptions outperform isolated PTCs for deep-time pseudogene dating.

// =============================================================================
// 2. RESULTS
// =============================================================================
#heading(level: 1, numbering: none)[2. Results]

#heading(level: 2, numbering: none)[2.1. Level 1 Annotation Audit Reveals Both Unannotated Relics and Automated Coding False Positives]

Across the 18 sauropsid RefSeq reference assemblies audited via NCBI Datasets v2 (all current chromosome-, contig-, or scaffold-level reference genomes; Table S1), the nine-gene canonical odontogenic panel (_ENAM_, _AMBN_, _AMTN_, _AMELX_, _MMP20_, _DSPP_, _ODAM_, _ODAPH_, _ACP4_) is completely represented in the mammalian outgroup _Homo sapiens_ `GRCh38.p14` (9/9 genes annotated as `PROTEIN_CODING`), whereas toothed sauropsids carry 4–7 officially annotated tooth genes per genome (7 in _Alligator mississippiensis_, 6 in _Crocodylus porosus_ and _Anolis carolinensis_, 5 in _Python bivittatus_, and 4 in _Gavialis gangeticus_ and _Podarcis muralis_). Missing official symbols in toothed species reflect nomenclature and pipeline gaps rather than biological absence: _Alligator_ _ODAPH_ is annotated as an uncharacterized locus (`LOC106740521`, `XP_019349100.2`, `NC_081825.1:164,672,404–164,676,187`), and functional _Anolis_ and _Podarcis_ _ODAPH_ orthologs are misannotated as non-coding RNAs (Table S1, Table S6).

Conversely, automated RefSeq annotation generated false-positive `PROTEIN_CODING` gene models on demonstrably pseudogenized loci in 4 of the 12 edentulous genomes (Table S1): _Chrysemys picta_ (_ENAM_, chr6), _Chelonia mydas_ (_AMBN_, chr5), _Struthio camelus_ (_AMELX_, chr1), and _Apteryx rowi_ (_AMELX_, unplaced scaffold). Level 3 alignment below confirms multiple frameshifts, premature stops, and a `GT→AT` splice-donor disruption directly inside the RefSeq-annotated _C. picta_ _ENAM_ and _C. mydas_ _AMBN_ coding models.

In _Gallus gallus_ `GRCg7b` (`GCF_016699485.2`), zero of the nine odontogenic genes are annotated genome-wide in either NCBI RefSeq Annotation Release 106 (released 2021-12-23) or Ensembl Release 116 (`Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.116.gtf.gz`, SHA-256 `407e3e27...`). Within the conserved 105,218-bp intergenic interval between _NRG1_ (`NC_052572.1:54,642,186–55,109,784`, `−`) and _LPL_ (`55,215,002–55,232,438`, `−`), RefSeq Release 106 annotates a single uncharacterized lncRNA (`LOC112530610`, `55,180,611–55,194,886`, `+`) and Ensembl Release 116 annotates five intergenic lncRNAs (`ENSGALG00010014552` through `ENSGALG00010014761`), none of which overlap the enamel-gene relics recovered below (Figure 1B; @tab1).

#heading(level: 2, numbering: none)[2.2. Level 2 Syntenic Recovery and Whole-Proteome Reciprocal Validation Across 18 Sauropsids]

As first noted for chicken _ENAM_ @alhashimi2010, the sauropsid _ENAM–AMBN–AMTN_ enamel cluster is conserved in strict synteny between _NRG1_ and _LPL_ (adjacent to _PSD3_) across all four extant reptile and bird clades (Figure 1A; Table S2a). Flanking _NRG1_ and _LPL_ anchors were recovered on the same chromosome or scaffold in 16 of 18 species (with single-anchor ±800-kb windows used for _Dromaius novaehollandiae_ and _Tinamus guttatus_ where one flanking symbol is unannotated in RefSeq), spanning 402–1,097 kb in toothed Crocodylia and Squamata, 160–337 kb in Testudines, and 105–330 kb on chromosome Z (or homologous scaffolds) in Aves.

#figure(
  image("figures/fig1_synteny_map.png", width: 100%),
  caption: [Conserved _NRG1–LPL_ enamel-cluster synteny, dual-database annotation audit of _Gallus gallus_ `GRCg7b` chromosome Z, and _ODAPH_ syntenic gap collapse. (*A*) Normalized intergenic positions (0% = _NRG1_ boundary, 100% = _LPL_ boundary; inverted scaffolds oriented _NRG1_ → _LPL_) of intact enamel genes (solid boxes) and whole-proteome reciprocal-confirmed pseudogene relics (_Ψ_, hatched boxes) across six representative sauropsids. (*B*) High-resolution map of the 105.2-kb _NRG1–LPL_ intergenic interval on `GRCg7b` chrZ (`NC_052572.1:55.098–55.236 Mb`), contrasting NCBI RefSeq Release 106 and Ensembl Release 116 features against three completely unannotated relics recovered in this study (_ΨENAM_ Exon 8, _ΨENAM_ Exon 6, and _ΨAMBN_ Exon 3; vertical shading confirms zero overlap). (*C*) Intergenic gap length between flanking anchors (_CDKL2_ and _THAP6_/_RCHY1_/_PARM1_/_SCLT1_) at the _ODAPH_ locus across 17 measurable sauropsids; all eight birds and _Pelodiscus sinensis_ exhibit gaps shorter than the 3.78-kb _Alligator_ _ODAPH_ gene span (dashed line), demonstrating physical deletion.],
) <fig1>

To establish orthology rigorously without small-database ascertainment bias, every candidate `tBLASTn` hit (bitscore ≥ 28 or _E_ ≤ 0.01) in all 18 _NRG1–LPL_ intervals was tested by reciprocal `BLASTx` against the entire _Alligator mississippiensis_ RefSeq proteome (48,040 proteins, `GCF_030867095.1`) both at the clustered locus-span level and at the individual top-HSP exon level (±30 bp; Table S2b). Evaluating both windows exposed a subtle analytical pitfall: when multi-exon pseudogene spans encompass >5 kb of degraded intronic sequence (_Struthio_ and _Apteryx_ _AMBN_, span 5.3 kb), intronic transposable or uncharacterized hits outscore the short \~150-bp coding exon in full-span `BLASTx`, whereas top-HSP exon reciprocal `BLASTx` unequivocally confirms _Alligator_ ameloblastin as the proteome-wide best hit (_Struthio_: `XP_019346384.1`, 51.2 bits, _E_ = 9.5×10#super[−9] vs. runner-up 28.9 bits; _Apteryx_: `XP_019346385.1`, 50.8 bits, _E_ = 8.9×10#super[−9] vs. runner-up 32.7 bits).

In _Gallus gallus_ `GRCg7b` chromosome Z (`NC_052572.1`), three distinct unannotated relics preserve the exact ancestral transcriptional orientation (`−` strand) and physical order (`ENAM` Exon 8 → `ENAM` Exon 6 → `AMBN` Exon 3; @tab1, Figure 1B):
- *_ΨENAM_ Exon 8* (`NC_052572.1:55,155,480–55,157,464`, `−`, 1,975 bp aligned to kakapo Exon 8 at _E_ = 7.3×10#super[−87], 317 bits; core `tBLASTn` vs. full RefSeq _Alligator_ `ENAM` `55,155,480–55,157,125`, 40.8 bits, _E_ = 1.2×10#super[−4]): single-copy in the _NRG1–LPL_ interval (1 hit in `dc-megablast`), with whole-proteome reciprocal `BLASTx` returning _Alligator_ enamelin (`XP_019346375.1`, 41.6 bits, _E_ = 0.004 full span; `XP_019346377.1`, 41.6 bits, _E_ = 2.7×10#super[−5] top HSP) and zero non-enamelin hits at _E_ ≤ 10 across 48,040 proteins in the full span.
- *_ΨENAM_ Exon 6* (`NC_052572.1:55,160,529–55,160,677`, `−`, 149 bp, situated 3.06 kb upstream of Exon 8 and 1.11 kb downstream of the non-overlapping Ensembl lncRNA `ENSGALG00010014557`): recovered at 78.0% nucleotide identity (_E_ = 2.7×10#super[−27], 115 bits) against kiwi _ENAM_ Exon 6 and carrying a direct 1-nt frameshift deletion plus an in-frame `TAA` premature stop codon (`cag|atttacaac-TAA-gaggcttct...`) while retaining canonical `AG...GT` splice borders (Table S5b).
- *_ΨAMBN_ Exon 3* (`NC_052572.1:55,168,855–55,168,992`, `−`, 138 bp): single-copy (1 `dc-megablast` hit), forward `tBLASTn` 56.6 bits (_E_ = 4.0×10#super[−10], 50.0% identity over 46 aa), whole-proteome reciprocal `BLASTx` confirmed to _Alligator_ ameloblastin (`XP_019346384.1`, 45.8 bits, _E_ = 5.1×10#super[−7] top HSP vs. non-ortholog runner-up 26.6 bits).

Beyond chicken, whole-proteome reciprocal validation confirmed _ΨENAM_ in 7 of 8 birds (all except _Tinamus guttatus_, whose weak 29.6-bit hit reciprocally matched kinase `DCLK3`) and in _Chrysemys picta_ and _Chelonia mydas_ (plus a marginal 30.0-bit, _E_ = 0.25 fragment in _Pelodiscus sinensis_, whereas _Gopherus evgoodei_'s 29.3-bit hit failed reciprocity). Crucially, we confirmed an unannotated *_ΨAMTN_ relic on kakapo chromosome Z* (_Strigops habroptila_ `bStrHab1.2.pri`, `NC_044302.2:51,216,319–51,216,444`, `−`, 126 bp, `tBLASTn` 37.4 bits, _E_ = 2.8×10#super[−4]; reciprocal best hit across 48,040 proteins = _Alligator_ amelotin `XP_014451634.1`, 38.1 bits, _E_ = 1.9×10#super[−4] vs. runner-up 26.9 bits) positioned in exact ancestral syntenic order downstream of _ΨENAM_ (@tab1; Figure 1A). A lower-confidence _ΨAMTN_ top-HSP candidate also occurs in exact syntenic order in mallard duck (_Anas platyrhynchos_, `NC_092621.1:53,220,153–53,220,311`, `+`, `tBLASTn` 34.3 bits, _E_ = 3.0×10#super[−3]; top-HSP reciprocal `XP_014451634.1`, 35.0 bits, _E_ = 5.0×10#super[−3], though its 7.6-kb multi-HSP span is masked by an intronic hit to `LOC132245913`, 134 bits; Table S2b).

#figure(
  table(
    columns: (1.10fr, 1.78fr, 0.72fr, 1.08fr, 1.38fr, 1.34fr),
    inset: 4.2pt,
    align: (left, left, center, right, left, left),
    stroke: (x, y) => if y == 0 { (top: 1pt, bottom: 0.7pt) } else if y == 5 { (bottom: 1pt) } else { none },
    [*Species & Relic*], [*`GRCg7b` / RefSeq Locus*], [*RefSeq / Ens.*], [*Forward Search*], [*Whole-Proteome Recip. Top Hit (48,040 prot.)*], [*Direct Lesion / Status*],
    [_G. gallus_\ _ΨENAM_ Ex8], [#tmono("NC_052572.1")\ #tmono("55,155,480–55,157,464 (−)")], [0 / 0], [40.8 b\ (_E_=1.2e−4)\*], [#tmono("XP_019346377.1")\ (41.6 b, _E_=2.7e−5)], [Galloanserae fs (c.619) + PTC (c.558)],
    [_G. gallus_\ _ΨENAM_ Ex6], [#tmono("NC_052572.1")\ #tmono("55,160,529–55,160,677 (−)")], [0 / 0], [115 b\ (_E_=2.7e−27)†], [Exon 6 MSA vs.\ _Amis_ + 4 birds], [1-nt del + in-frame `TAA` PTC],
    [_G. gallus_\ _ΨAMBN_ Ex3], [#tmono("NC_052572.1")\ #tmono("55,168,855–55,168,992 (−)")], [0 / 0], [56.6 b\ (_E_=4.0e−10)], [#tmono("XP_019346384.1")\ (45.8 b, _E_=5.1e−7)], [Single-copy (1 hit); unannotated],
    [_S. habroptila_\ _ΨAMTN_], [#tmono("NC_044302.2")\ #tmono("51,216,319–51,216,444 (−)")], [0 / N/A], [37.4 b\ (_E_=2.8e−4)], [#tmono("XP_014451634.1")\ (38.1 b, _E_=1.9e−4)], [Preserves _ENAM→AMTN_ order],
    [_G. gallus_\ _ΨAMEL_ (cal.)], [#tmono("NC_052532.1")\ #tmono("124,410,142–124,413,840 (−)")], [0 / 0], [100% BLASTn\ (`EU340348`)], [0 HSPs vs. _Amis_\ `AMELX` (`-seg yes`)], [Empirical `tBLASTn` blind spot],
  ),
  caption: [Summary of key unannotated avian enamel-gene relics audited in _Gallus gallus_ `GRCg7b` (`GCF_016699485.2`, RefSeq Release 106 & Ensembl Release 116) and _Strigops habroptila_ (`GCF_004027225.2`), alongside empirical blind-spot calibration on chicken _ΨAMEL_ (`EU340348`). \*Forward `tBLASTn` core spans `55,155,480–55,157,125` (40.8 bits); full Exon 8 alignment to _S. habroptila_ spans `55,155,480–55,157,464` (1,975 bp, 317 bits, _E_ = 7.3×10#super[−87]). †`BLASTn` against _A. rowi_ _ENAM_ Exon 6. Full 18-species reciprocal results in Table S2b.],
) <tab1>

Four auxiliary controls reinforce and delimit these Level 2 findings (Table S6):
1. *_ODAPH_ physical deletion vs. _ODAM_ uninformative negative:* Across 17 measurable sauropsids, toothed controls retain 16.8–47.8 kb between _CDKL2_ and the adjacent downstream anchor (_THAP6_/_PARM1_), housing an intact 3,784-bp _ODAPH_ gene in _Alligator_ (Figure 1C). All eight birds exhibit drastic syntenic gap collapse to 0–2,447 bp (`0 bp` in _G. gallus_, `146 bp` in _A. platyrhynchos_, `859–2,447 bp` in the remaining six birds), and _Pelodiscus sinensis_ collapses to `132 bp`—physically shorter than the _Alligator_ _ODAPH_ locus and providing Level 2 evidence of deletion rather than sequence divergence. By contrast, _ODAM_ yielded no significant hit even in toothed _Alligator_, _Anolis_, or _Python_ (best `tBLASTn` _E_ ≥ 0.21; mammalian/amphibian sensitivity test on _Xenopus_ yielded _E_ = 0.042–1.6), classifying _ODAM_ absence as an uninformative Level 1 negative.
2. *Empirical `tBLASTn` blind-spot calibration (_ΨAMEL_):* The chicken _ΨAMEL_ pseudogene sequenced by #pcite(<sire2008>) (`EU340348`) maps at 100.0% nucleotide identity across four exon blocks (513 bp) to `GRCg7b` chr1 (`NC_052532.1:124,410,142–124,413,840`, `−`) between intact flanking anchors _MSL3_ (121 bits) and _MID1_ (455 bits). Yet querying this exact interval with _Alligator_ `AMELX` yields 0 HSPs under `-seg yes` and a noise-floor 25.8 bits under `-seg no`, proving that absence of `tBLASTn` signal cannot rule out a retained pseudogene relic.
3. *Independent whole-genome alignment coverage (`multiz100way`):* Across 29 human `hg38` coding exons of _ENAM_, _AMBN_, and _AMTN_, UCSC `multiz100way` independently aligns 56.4–63.0% of _ENAM_ CDS in _C. mydas_ and _C. picta_, 41.0% in _Amazona vittata_, and 64.4% of _AMBN_ CDS in two falcons (_Falco cherrug_, _F. peregrinus_); however, functional _Alligator_ (`allMis1`) yields 54 raw frameshifts against `hg38`, restricting `multiz100way` strictly to coverage-only support (Table S6).
4. *Chromosome Z recombination context:* Anchoring `GRCg7b` chicken _ΨENAM_ and _ΨAMBN_ to `galGal3` via exact 32-mer matching (4,380/4,380 and 1,893/2,650 concordant 32-mers) places both relics inside recombination bin 1784 of the #pcite(<elferink2010>) linkage map (male recombination rate 3.94 cM/Mb; _NRG1–LPL_ window weighted average 2.27 cM/Mb, 54th percentile of 122 chrZ bins vs. chrZ mean 3.00 and autosomal mean 2.97 cM/Mb; genome-wide baseline from #pcite(<groenen2009>)), rejecting a recombination cold-spot explanation for relic persistence.

#heading(level: 2, numbering: none)[2.3. Compositional, SEG-Masking, and Phylogenetic Bias in Translated Searches of Enamel IDPs]

Sequence profiling with `segmasker` (`-window 12 -locut 2.2 -hicut 2.5`) confirms stark compositional divergence between globular flanking controls and mineralized-tissue IDPs (Figure 2A; Table S3a). Globular _LPL_ exhibits low Pro+Gln content (6.9–8.2%), high Shannon entropy (4.19–4.22 bits/aa), and 0.0% (_Anolis_) to 3.1% (_Alligator_) `SEG`-masked residues. By contrast, sauropsid enamel proteins (_ENAM_, _AMBN_, _AMTN_) contain 18.0–25.9% Pro+Gln and 46.6–57.3% disorder-promoting SCPP residues (P, Q, S, D, E, G, N), with 5.0–13.9% `SEG`-masked sequence. Functional osteogenic SCPP/SIBLING IDPs (_IBSP_, _SPP1_, _SPARCL1_)—intact across all 18 amniotes—exhibit even stronger low-complexity enrichment, reaching 25.5% (_SPP1_) and 52.8% (_IBSP_) `SEG`-masked residues in _Anolis_ (14.4% and 46.0% in _Alligator_), while the cross-chromosome negative control _DSPP_ is 47.4% `SEG`-masked.

Evaluating null-corrected retention indices ($sum "excess"("ENAM", "AMBN", "AMTN") \/ "excess"("LPL")$, using 20 composition-matched random permutations per query, seed `20261017`) across all 18 species reveals two compounding artifacts in translated bitscore metrics (Figure 2B–C; @tab2, Table S3b–c):

First, *phylogenetic query asymmetry confounds non-equidistant scores.* When queried with full-length _Alligator_ proteins (`-seg yes`), crocodylians score 12.60–15.39, toothed squamates (diverged \~281 Ma per TimeTree 5, #pcite(<kumar2022>)) drop to 1.49–2.59, and edentulous _Chelonia mydas_ (2.17) and _Chrysemys picta_ (1.55) directly overlap or exceed toothed _Python bivittatus_ (1.67) and _Podarcis muralis_ (1.49; Figure 2C1). Switching to an equidistant _Anolis carolinensis_ query—diverged \~281 Ma identically from Crocodylia, Testudines, and Aves—eliminates this phylogenetic overlap cleanly: toothed Crocodylia score 2.25–2.78, edentulous Testudines drop to 0.00–0.56, and Aves drop to 0.00–0.14 (Figure 2C2).

Second, and more fundamentally, *functional intrinsically disordered proteins suffer steep bitscore attenuation in Aves even under purifying selection.* Under the equidistant _Anolis_ query (@tab2; Figure 2B):
- Globular *_LPL_* (`0.0%` SEG) is completely invariant to `-seg yes` vs. `-seg no` and scores _higher_ in Aves (median 136.1 bits, range 132.3–143.4, *111.7%* of Crocodylia 121.8 bits), ruling out generalized avian substitution-rate acceleration or locus-wide alignment breakdown across the _NRG1–LPL_ interval.
- Functional bone *_IBSP_* (`52.8%` SEG) drops in Aves to median 40.9 bits (*73.8%* of Crocodylia 55.4) under `-seg yes`, and when low-complexity masking is disabled (`-seg no`), Crocodylia jumps to 100.3 bits while Aves rises only to 48.0 bits (*47.9%* of Crocodylia). Using closer archosaurian _Alligator_ queries (\~240 Ma), functional avian _IBSP_ drops to just 8.1–27.7 bits excess (median 17.1, *15.0%* of Crocodylia 114.1 bits; Table S3c).
- Functional bone *_SPP1_* (`25.5%` SEG) drops in Aves to median 19.0 bits (*54.9%* of Crocodylia 34.6 bits; range *6.8–25.3 bits* under `-seg yes`), directly overlapping the excess bitscore range of decayed avian _ENAM_ (0.0–14.4 bits) and _AMBN_ (0.0–15.3 bits; Figure 2B).

The cross-chromosome negative control _DSPP_ (located outside _NRG1–LPL_) confirmed near-zero median excess (0.0 bits across all four clades under both `-seg yes` and `-seg no`, max observed 29.6 bits), while exposing an edge-case `SEG` artifact where heavy masking collapsed _DSPP_ shuffle baselines to 0.0 in one turtle (Table S3c). Together, these controls prove quantitatively that depressed `tBLASTn` bitscores in Pro/Gln- or Asp/Glu-rich IDPs cannot serve as standalone evidence of pseudogenization.

#figure(
  image("figures/fig2_compositional_bias.png", width: 100%),
  caption: [Quantitative dissection of compositional, `SEG`-masking, and phylogenetic bias in translated homology searches across 18 sauropsids. (*A*) Pro+Gln fraction versus `SEG`-masked low-complexity fraction (`segmasker`) across 14 representative queries (marker area proportional to protein length). (*B*) Equidistant _Anolis carolinensis_ (\~281 Ma) `tBLASTn` excess bitscores above the 20-shuffle 95th-percentile null (`-seg yes`, filled symbols; `-seg no`, open symbols; horizontal black bars = clade medians) across Crocodylia (_n_ = 3), Testudines (_n_ = 4), and Aves (_n_ = 8). Globular _LPL_ scores 111.7% of crocodylian levels in birds with zero `SEG` effect, whereas functional osteogenic IDPs (_IBSP_, _SPP1_) lose 26–52% of excess bitscore in birds—with intact avian _SPP1_ (6.8–25.3 bits) directly overlapping decayed avian _ENAM_/_AMBN_ (shaded pink band). (*C1–C2*) Null-corrected retention indices across all 18 species under (*C1*) non-equidistant _Alligator_ full RefSeq queries (where edentulous _C. mydas_ 2.17 overlaps toothed squamates 1.49–1.67) versus (*C2*) equidistant _Anolis_ queries (clean separation: Crocodylia 2.25–2.78 vs. Testudines 0.00–0.56 and Aves 0.00–0.14).],
) <fig2>

#figure(
  table(
    columns: (1.75fr, 0.72fr, 1.18fr, 1.18fr, 1.18fr, 0.85fr),
    inset: 4.3pt,
    align: (left, left, right, right, right, right),
    stroke: (x, y) => if y == 0 { (top: 1pt, bottom: 0.7pt) } else if y == 10 { (bottom: 1pt) } else { none },
    [*Equidistant Query (_Anolis_, \~281 Ma)*], [*`SEG` Mode*], [*Crocodylia (_n_=3)*], [*Testudines (_n_=4)*], [*Aves (_n_=8)*], [*Aves / Croc (%)*],
    [_LPL_ (globular, 0.0% SEG)], [`-seg yes`], [121.8 (119.3–122.3)], [122.7 (118.6–128.9)], [136.1 (132.3–143.4)], [111.7%],
    [_LPL_ (globular, 0.0% SEG)], [`-seg no`], [121.8 (119.3–122.3)], [122.5 (118.6–128.9)], [136.1 (132.3–143.4)], [111.7%],
    [_IBSP_ (bone IDP, 52.8% SEG)], [`-seg yes`], [55.4 (54.3–55.8)], [46.3 (45.5–48.6)], [40.9 (27.7–46.1)], [73.8%],
    [_IBSP_ (bone IDP, 52.8% SEG)], [`-seg no`], [100.3 (75.8–107.4)], [84.5 (72.6–91.6)], [48.0 (38.1–64.7)], [47.9%],
    [_SPP1_ (bone IDP, 25.5% SEG)], [`-seg yes`], [34.6 (31.9–35.9)], [39.1 (30.6–45.9)], [19.0 (6.8–25.3)], [54.9%],
    [_SPP1_ (bone IDP, 25.5% SEG)], [`-seg no`], [35.8 (35.3–38.0)], [48.2 (36.8–55.1)], [20.8 (13.3–27.4)], [58.1%],
    [_ENAM_ (enamel IDP, 7.5% SEG)], [`-seg yes`], [246.2 (231.4–253.2)], [27.1 (0.0–60.0)], [1.5 (0.0–14.4)], [0.6%],
    [_ENAM_ (enamel IDP, 7.5% SEG)], [`-seg no`], [263.9 (250.4–265.2)], [30.0 (0.0–60.0)], [2.1 (0.0–13.6)], [0.8%],
    [_AMBN_ (enamel IDP, 5.0% SEG)], [`-seg yes`], [42.3 (35.6–42.6)], [5.6 (0.0–13.7)], [1.1 (0.0–15.3)], [2.6%],
    [_AMBN_ (enamel IDP, 5.0% SEG)], [`-seg no`], [41.9 (35.6–42.2)], [5.6 (0.0–12.3)], [0.2 (0.0–16.1)], [0.5%],
  ),
  caption: [Clade medians (and min–max ranges) of `tBLASTn` excess bitscores above the 20-shuffle 95th-percentile null under equidistant _Anolis carolinensis_ queries (`-seg yes` vs. `-seg no`). Values are read directly from Table S3c. Functional avian _SPP1_ (6.8–25.3 bits) overlaps decayed avian _ENAM_ (0.0–14.4 bits) and _AMBN_ (0.0–15.3 bits).],
) <tab2>

#heading(level: 2, numbering: none)[2.4. Level 3 Splice-Aware Coding Lesions and Splice-Site Disruptions]

To separate structural pseudogenization from fast-evolving functional disorder, we aligned full-length _Alligator_ RefSeq proteins (`ENAM` `XP_019346377.1`, 1,078 aa; `AMBN` `XP_019346386.1`, 390 aa; `AMTN` `XP_014451634.1`, 328 aa; flanking positive control `LPL` `XP_019346320.1`, 477 aa) across all 18 _NRG1–LPL_ intervals with `miniprot` (v0.18-r281; Figure 3A, Table S4).

The flanking positive control *_LPL_* aligned cleanly across *18/18 species* (8,246 aligned aa total: 2,695 toothed + 5,551 edentulous) with *0 frameshifts, 0 premature stop codons (0.00 lesions/1,000 aa), and 133/133 canonical `GT-AG` introns* (0 non-canonical splice junctions), ruling out systematic `miniprot` alignment artifacts or assembly-quality confounding across the exact same genomic intervals. Similarly, toothed enamel controls (_n_ = 15 alignments spanning 9,394 aligned aa across Crocodylia and Squamata, with _Gavialis_ _AMTN_ partially truncated by its 14.8-kb scaffold contig N50) exhibited only 1 isolated `miniprot` frameshift call in _Anolis_ _AMBN_ (`fs@171` at 45.6% aa identity across \~281 Ma divergence, despite an intact native RefSeq ORF `NP_001315156.1` of 377 aa), establishing an empirical cross-clade alignment false-positive floor of *0.11 lesions/1,000 aa* (0 premature stop codons; 92/93 canonical `GT-AG` introns).

In stark contrast, edentulous enamel pseudogenes exhibited *41 frameshifts and 29 premature stop codons* (70 coding lesions across 4,122 pooled aligned aa = *16.98 lesions/1,000 aa*, or *17.49 lesions/1,000 aa* across the 7 substantial lesioned alignments spanning 4,002 aa—a *>150-fold enrichment* above the empirical toothed floor; Figure 3A):
- _Chrysemys picta_ _ENAM_: 12 frameshifts + 4 PTCs (22.44/1,000 aa across 713 aa, 66.1% cov) — directly invalidating its RefSeq `PROTEIN_CODING` annotation.
- _Strigops habroptila_ _ENAM_: 8 frameshifts + 9 PTCs (19.63/1,000 aa across 866 aa, 80.3% cov).
- _Apteryx rowi_ _ENAM_: *7 frameshifts + 8 PTCs (19.18/1,000 aa across 782 aa, 72.5% cov) plus a disrupted `GT→AT` splice donor* spanning a 1,100-bp intron (Intron 6). Notably, using the legacy partial GenBank entry `QOX06214.1` (849 aa, missing 212 N-terminal residues) completely failed to align _A. rowi_ _ENAM_ (0% coverage; Table S4), demonstrating how query completeness governs pseudogene recovery.
- _Chelonia mydas_ _ENAM_: 8 frameshifts + 3 PTCs (19.16/1,000 aa across 574 aa, 53.2% cov).
- _Dromaius novaehollandiae_ _AMBN_: 1 frameshift + 3 PTCs (13.94/1,000 aa across 287 aa, 73.6% cov).
- _Chrysemys picta_ _AMBN_: 3 frameshifts + 1 PTC (10.26/1,000 aa across 390 aa, 100% cov).
- _Chelonia mydas_ _AMBN_: *2 frameshifts + 1 PTC (7.69/1,000 aa across 390 aa, 100% cov) plus a disrupted `GT→AT` splice donor* spanning a 1,552-bp intron (Intron 8) — directly invalidating its RefSeq `PROTEIN_CODING` annotation.

Because `miniprot` depends on seeded 4-mer amino-acid chains against a crocodylian query, highly degraded avian _ENAM_ loci falling below \~35% protein identity with dense frameshifts (_G. gallus_, _A. platyrhynchos_, _D. novaehollandiae_, _S. camelus_) were not spanned by `miniprot` even though intra-avian nucleotide `BLASTn` (63–89% nt identity) readily recovered 1.7–2.5 kb of orthologous Exon 8 and Exon 6 in all seven birds for Level 4 analysis below. Audit of `miniprot` intron lengths also uncovered an instructive quality-control failure mode: full-length `ENAM` initially chained an implausible 58.4-kb intron in _C. picta_ to an upstream 382-bp low-complexity hit (`NC_088796.1:129,180,576–129,180,957`), which whole-proteome reciprocal `BLASTx` proved was *protein SON* (_E_ = 0.010), not enamelin. Imposing a biological 25,000-bp maximum-intron filter segmented this distal chaining artifact cleanly without altering any genuine enamelin exon (Table S4).

#heading(level: 2, numbering: none)[2.5. Level 4 Dual-Aligner Shared Inactivating Mutations, Splice Disruptions, and Pervasive PTC Homoplasy]

To test which inactivating lesions reflect shared ancestry versus convergent decay, we extracted verified coding exons of _Alligator_ _ENAM_ (Exons 6 and 8; complete 8-exon CDS = 3,237 bp, 100% exact translation to `XP_019346377.1`) and _AMBN_ (Exons 3, 4, and 9; complete 9-exon CDS = 1,173 bp, 100% exact translation to `XP_019346386.1`), recruited homologous avian and testudine fragments, and aligned each exon independently with *MAFFT `L-INS-i`* (v7.525) and *MUSCLE 5* (v5.3) alongside all three toothed crocodylians (_Alligator_, _Crocodylus_, _Gavialis_). All three crocodylian controls yielded *0 coding lesions in both aligners across every dataset* (Table S5a).

Requiring exact reference coordinate, mutation-class, and carrier-set agreement across both aligners is intentionally conservative—filtering out alignment-sensitive indels whose gap boundaries shift by 1–2 nt between MAFFT and MUSCLE (e.g., yielding 17 dual-aligner consensus sites from 31 MAFFT and 33 MUSCLE candidates in avian _ENAM_ Exon 8). Across all datasets, this strict consensus yielded *23 robust shared coding mutations* (17 phylogenetically congruent synapomorphies and 6 incongruent homoplasies) plus *2 shared Exon 6 splice-border disruptions* (Figure 3B–C; @tab3, Table S5a–b):

1. *Palaeognathae (_Apteryx_, _Dromaius_, _Struthio_):* Anchored by a *16-nt frameshift deletion in _AMBN_ Exon 4* (codon 109 of `XP_019346386.1`, verified outside homopolymer runs), *3 congruent PTCs in _ENAM_ Exon 8* (`Q370*` `CAA→TGA`, `F601*` `TTT→TAA`, `W751*` `TGG→TGA`), and a *shared 3' splice-acceptor disruption at _ENAM_ Exon 6* (`AG→AC` in _Struthio_ and _Apteryx_, further mutated to `GT` in _Dromaius_, while _Alligator_, _Gallus_, and _Anas_ preserve canonical `AG`).
2. *Novaeratitae (_Apteryx_ + _Dromaius_):* Anchored by *5 congruent PTCs in _ENAM_ Exon 8* (`G270*`, `S401*`, `S559*`, `Q667*`, `H844*`), *2 congruent PTCs in _AMBN_ Exon 3* (`S56*` `TCA→TAA`, `Q62*` `CAA→TGA`), and a *shared canonical 5' splice-donor disruption at _ENAM_ Exon 6 (`GT→AT`)* confirming the `miniprot` lesion independently across both ratite assemblies.
3. *Galloanserae (_Gallus_ + _Anas_):* Adding _Gallus gallus_ Exon 8 (1,975 bp) directly into the 7-bird alignment recovered *1 congruent frameshift indel* (2-nt insertion relative to crocodylians at codon 619, equivalent to a 1-nt deletion within an ancestral avian 3-nt insertion) plus *1 congruent PTC* (`G558*` `GGA→TGA`), dating chicken _ENAM_ inactivation minimally to crown Galloanserae.
4. *Psittacopasserae (_Strigops_ + _Taeniopygia_):* Anchored by *1 congruent PTC* in _ENAM_ Exon 8 (`W1003*` `TGG→TGA`).
5. *Crown Durocryptodira (_Chelonia_ + _Chrysemys_):* Anchored by *1 congruent 1-nt frameshift insertion in _AMBN_ Exon 4* (codon 123) and *2 congruent PTCs in _ENAM_ Exon 8* (`Q753*` `CAG→TAG`, `A859*` `GCA→TGA`), dating testudine inactivation minimally to crown Durocryptodira (\~105–157 Ma).
6. *Crown Neornithes — Zero Congruent Lesions and 37.5% PTC Homoplasy:* Across all 7 aligned birds, *zero coding lesions or splice disruptions were shared congruently across crown Neornithes* (Figure 3B). Strikingly, *6 of the 16 dual-aligner shared avian _ENAM_ PTCs (37.5%)* are phylogenetically incongruent homoplasies spanning disjoint clades (Figure 3C; Table S5a). Their sequence context explains this convergence mechanistically: 1 occurs at a hypermutable `CpG` dinucleotide (codon 311 `CCG→TGA` in _Dromaius_ + _Struthio_ while _Apteryx_ retains `CGA`), 2 occur at single-step `TGG` tryptophan hotspots (codons 562 and 615), and *3 of the 6 homoplastic sites carry distinct stop triplets across lineages*—proving independent mutational origins directly from triplet identity. Most tellingly, at _ENAM_ codon 684 (ancestral `CAA`), all three palaeognaths share a `CAA→TAG` transition (a genuine Palaeognathae synapomorphy), whereas _Gallus gallus_ independently mutated the exact same codon to `TAA` while _Anas_ and _Strigops_ retain sense `CAA` (Figure 3C).

#figure(
  image("figures/fig3_lesions_and_phylogeny.png", width: 100%),
  caption: [Level 3 coding-lesion accumulation and Level 4 dual-aligner phylogenetic mapping of shared inactivating mutations and PTC homoplasy. (*A*) Level 3 `miniprot` coding-lesion rates (frameshifts in blue, PTCs in vermilion, per 1,000 aligned aa) contrasting 18/18 flanking _LPL_ positive controls (0.00/1,000 aa, 133/133 canonical `GT-AG` introns) and 15 toothed enamel alignments (0.11/1,000 aa empirical cross-clade floor) against pooled and individual edentulous pseudogenes (16.98/1,000 aa pooled); magenta diamonds mark disrupted `GT→AT` splice donors in _A. rowi_ _ENAM_ and _C. mydas_ _AMBN_. (*B*) Cladogram of sampled sauropsids showing all 17 dual-aligner (`MAFFT L-INS-i` + `MUSCLE 5`) phylogenetically congruent coding synapomorphies and 2 _ENAM_ Exon 6 splice-border disruptions mapped to crown nodes (0 lesions in all 3 crocodylian outgroups; 0 congruent lesions at crown Neornithes). (*C*) Per-species codon matrix across the six phylogenetically incongruent avian _ENAM_ Exon 8 PTCs (37.5% of shared PTCs), displaying exact stop triplets (`TGA`, `TAG`, `TAA`), intact sense codons, and gaps; three sites (`c.439`, `c.615`, `c.684`) carry two distinct stop codons across clades, providing direct triplet-level proof of convergent origin.],
) <fig3>

#figure(
  table(
    columns: (1.42fr, 1.05fr, 0.72fr, 0.50fr, 2.75fr),
    inset: 4.3pt,
    align: (left, left, right, right, left),
    stroke: (x, y) => if y == 0 { (top: 1pt, bottom: 0.7pt) } else if y == 7 { (bottom: 1pt) } else { none },
    [*Phylogenetic Node / Category*], [*Carriers in Sample*], [*Frameshifts*], [*PTCs*], [*Key Diagnostic Lesions & Splice Disruptions*],
    [Palaeognathae (congruent)], [_Arow, Dnov, Scam_], [1], [3], [_AMBN_ Ex4 16-nt del (c.109); _ENAM_ Ex8 Q370\*, F601\*, W751\*; Ex6 `AG→AC/GT`],
    [Novaeratitae (congruent)], [_Arow, Dnov_], [0], [7], [_ENAM_ Ex8 5 PTCs; _AMBN_ Ex3 S56\*, Q62\*; _ENAM_ Ex6 `GT→AT` donor],
    [Galloanserae (congruent)], [_Ggal, Apla_], [1], [1], [_ENAM_ Ex8 2-nt/1-nt fs (c.619); G558\* (`GGA→TGA`)],
    [Psittacopasserae (congruent)], [_Shab, Tgut_], [0], [1], [_ENAM_ Ex8 W1003\* (`TGG→TGA`)],
    [Durocryptodira (congruent)], [_Cmyd, Cpic_], [1], [2], [_AMBN_ Ex4 1-nt ins (c.123); _ENAM_ Ex8 Q753\*, A859\*],
    [Crown Neornithes (congruent)], [All 7 aligned birds], [0], [0], [No single shared lesion spans Palaeognathae + Neognathae in this sample],
    [Avian _ENAM_ Homoplasies], [Disjoint clades], [0], [6], [37.5% (6/16) of shared PTCs: 1 CpG (`c.311`), 2 `TGG` (`c.562, 615`), 3 dual-stop (`c.439, 615, 684`)],
  ),
  caption: [Summary of Level 4 dual-aligner (`MAFFT L-INS-i` + `MUSCLE 5`) consensus shared inactivating mutations and Exon 6 splice-border disruptions across avian and testudine _ENAM_ and _AMBN_ (0 lesions across all 3 crocodylian outgroups in every dataset; full per-site details in Table S5a–b).],
) <tab3>

// =============================================================================
// 3. DISCUSSION
// =============================================================================
#heading(level: 1, numbering: none)[3. Discussion]

#heading(level: 2, numbering: none)[3.1. Calibrating Replication Against Original Contributions]

Our findings build directly upon the classic paleogenomic demonstration of independent tooth loss in birds and turtles @sire2008 @alhashimi2010 @meredith2013 @meredith2014. Specifically, prior work first identified the translocated _NRG1–LPL_ location of sauropsid _ENAM–AMBN–AMTN_, recovered avian _ENAM_/_AMBN_/_AMTN_ fragments in earlier draft assemblies, and documented physical loss of _ODAPH_ in edentulous lineages @alhashimi2010 @meredith2013 @meredith2014 @springer2016. Against that established baseline, our systematic 18-genome audit delivers four distinct empirical and methodological advances: (1) dual-database proof (`RefSeq 106` and `Ensembl 116`) that _all_ enamel relics remain completely unannotated in the modern chromosome-level chicken reference `GRCg7b`, accompanied by exact single-copy, whole-proteome reciprocal-confirmed coordinates—including the lesion-bearing _ΨENAM_ Exon 6 (`NC_052572.1:55,160,529–55,160,677`, `−`) and a reciprocal-confirmed _ΨAMTN_ relic on kakapo chromosome Z (`NC_044302.2:51,216,319–51,216,444`, `−`); (2) quantitative 18-species proof that functional osteogenic IDPs (_IBSP_, _SPP1_) undergo severe translated bitscore collapse in birds overlapping decayed enamel pseudogenes; (3) demonstration that replacing truncated GenBank queries with full-length RefSeq orthologs recovers previously invisible pseudogene exons and splice disruptions (e.g., _Apteryx rowi_ _ENAM_) while requiring intron-length filters to reject low-complexity chaining artifacts (`protein SON` in _C. picta_); and (4) dual-aligner quantification showing that 37.5% (6/16) of shared avian _ENAM_ premature stop codons are homoplastic.

#heading(level: 2, numbering: none)[3.2. Bidirectional Blind Spots of Automated Genome Annotation]

Automated annotation pipelines remain the entry point for most large-scale comparative genomic scans, yet our Level 1 audit illustrates why annotation tables alone are unreliable for inferring gene presence or loss @weisman2022 @kirilenko2023. In _Gallus gallus_ `GRCg7b`, both NCBI RefSeq Release 106 and Ensembl Release 116 report zero enamel genes despite three single-copy, syntenic, whole-proteome reciprocal-confirmed exon relics residing in a normal-recombination region of chromosome Z (54th percentile of chrZ recombination bins) @elferink2010. Equally concerning is the converse error: four edentulous RefSeq assemblies (_C. picta_, _C. mydas_, _S. camelus_, _A. rowi_) annotate frame-disrupted or splice-mutated enamel loci as `PROTEIN_CODING`, while functional toothed-reptile _ODAPH_ loci are labeled uncharacterized (`LOC106740521`) or non-coding RNA. Without synteny-anchored sequence extraction and coding-frame verification, database-level presence/absence matrices systematically inject both false negatives and false positives into genotype–phenotype association frameworks @prudent2016 @kowalczyk2019 @saputra2021 @redlich2024.

#heading(level: 2, numbering: none)[3.3. Why Intrinsically Disordered Proteins Break Naive Homology Metrics]

Enamel matrix proteins evolved within the SCPP family as flexible calcium- and mineral-binding scaffolds rich in Pro, Gln, Ser, and acidic residues @kawasaki2006 @kawasaki2021 @kawasaki2023 @chun2025. Standard `tBLASTn` scoring matrices and `SEG` low-complexity masking were engineered primarily for globular domains @wootton1993 @schaffer2001 @yu2003 @gertz2006. Our equidistant _Anolis_ benchmark (\~281 Ma to Crocodylia, Testudines, and Aves) isolates the resulting failure mode cleanly: intact globular _LPL_ shows zero `SEG` sensitivity and scores 111.7% of crocodylian bitscore in birds, ruling out generalized avian rate acceleration across this locus. Yet intact, functionally indispensable osteogenic IDPs (_IBSP_, _SPP1_) lose 26–52% of their excess bitscore in birds—dropping to 6.8–25.3 bits for avian _SPP1_, completely overlapping decayed avian _ENAM_ (0.0–14.4 bits) and _AMBN_ (0.0–15.3 bits).

This empirical overlap echoes theoretical predictions of homology detection failure in fast-evolving disordered regions @weisman2020 @chow2024 and yields a concrete methodological prescription for paleogenomics: *bitscore attenuation in IDPs should never be interpreted as pseudogenization without corroborating Level 3 coding lesions or Level 2 physical gap collapse.* Our empirical blind-spot calibration on chicken _ΨAMEL_ (`EU340348`) underscores the symmetric caveat: a confirmed 513-bp pseudogene yields zero `tBLASTn` HSPs under default `-seg yes`, confirming that absence of translated sequence similarity does not prove physical deletion. Furthermore, our top-HSP vs. full-span reciprocal `BLASTx` comparison highlights that whole-proteome reciprocity should be evaluated at the exon level whenever multi-kilobase pseudogene introns can harbor transposable-element hits.

#heading(level: 2, numbering: none)[3.4. Homoplastic Premature Stops Versus Indel and Splice-Site Synapomorphies]

Dating ancient phenotypic losses from shared inactivating mutations assumes that identical lesions across descendant clades reflect a single ancestral mutational event rather than convergence @meredith2013 @meredith2014 @emerling2026. Our dual-aligner consensus analysis demonstrates why premature stop codons carry a high homoplasy burden after >100 Myr of neutral evolution: over one-third (6/16 = 37.5%) of dual-aligner shared PTCs in avian _ENAM_ Exon 8 are phylogenetically incongruent, driven by `CpG` deamination (`c.311`), single-step `TGG` tryptophan transitions (`c.562`, `c.615`), and multi-stop hits (`c.439`, `c.615`, `c.684`) where distinct stop triplets (`TAG` vs. `TAA` or `TGA`) directly betray independent origins. Notably, #pcite(<meredith2014>) detected a shared crown-avian frameshift deletion in _ENAM_ across a much broader 48-bird sampling and inferred a single avian tooth loss \~116 Ma consistent with Mesozoic fossil transitions @brocklehurst2021; in our 8-bird RefSeq panel, nucleotide divergence and partial exon deletion prevented recovering any single crown-Neornithes synapomorphy, whereas lower-homoplasy multi-nucleotide frameshift indels (Palaeognathae _AMBN_ Exon 4 16-nt deletion; Galloanserae _ENAM_ Exon 8 codon 619 frameshift) and canonical splice-junction disruptions (Palaeognathae _ENAM_ Exon 6 `AG→AC/GT`; Novaeratitae `GT→AT`) provided crisp, aligner-invariant minimum dates for crown subclades and crown Durocryptodira.

#heading(level: 2, numbering: none)[3.5. Limitations]

Five limitations delimit the scope of our conclusions and motivate targeted follow-up work:
1. *Phylogenetic sample size and inferential scope:* Although our panel spans 18 genomes across all four major sauropsid clades, birds and turtles represent only _n_ = 2 independent evolutionary origins of edentulism. Because trait-association p-values across _n_ = 2 evolutionary replicates are severely underpowered and confounded by shared ancestry @saputra2021, all clade comparisons here are intentionally descriptive and anchored by within-genome controls (_LPL_, _IBSP_, _SPP1_, crocodylian outgroups).
2. *Assembly contiguity heterogeneity:* While 14 of 18 reference genomes are chromosome-level or complete T2T assemblies (including `GRCg7b`, `IASCAAS_PekinDuck_T2T`, and `ASM5242630v1`), four assemblies remain scaffold- or contig-level (_Gavialis_ contig N50 14.8 kb; _Crocodylus_ 34.1 kb; _Apteryx rowi_ `aptRow1` 129.4 kb; _Tinamus_ 3.85 Mb; Table S1). Isolated base-call or homopolymer indels cannot be excluded in lower-contiguity assemblies; however, our key Level 4 synapomorphies are independently validated across high-contiguity assemblies (_D. novaehollandiae_ N50 24.2 Mb + _S. camelus_ 28.8 Mb + _A. rowi_; _G. gallus_ 18.8 Mb + T2T _A. platyrhynchos_ 80.2 Mb; _C. mydas_ 39.4 Mb + _C. picta_ 13.4 Mb).
3. *Testudine basal sampling:* All four sampled turtles belong to Cryptodira (three Durocryptodira + one Trionychia, _Pelodiscus sinensis_). Inclusion of chromosome-level Pleurodira assemblies is required to test whether _ENAM_/_AMBN_ frameshifts predate crown Testudines around \~157–210 Ma @shaffer2013 @meredith2013.
4. *Assembly version alignment in auxiliary _ODAPH_ panel:* Primary _NRG1–LPL_ and _IBSP–SPP1_ analyses for _Pelodiscus sinensis_ use its current T2T reference genome (`GCF_052426305.1`, released 2025-09-03, annotated 2026-10-19), whereas auxiliary _ODAPH_ gap coordinates for _P. sinensis_ (`NC_134715.1`, 132-bp gap) originate from the prior RefSeq representative (`GCF_049634645.1`), explicitly noted in Table S6.
5. *Selection-intensity and transcriptional profiling:* Branch-specific $d_N\/d_S$ relaxation modeling (e.g., `HyPhy RELAX`) and embryonic oral RNA-seq verification of residual _ΨENAM_/_ΨAMBN_ transcription were outside the scope of this genomic audit and represent logical next steps.

// =============================================================================
// 4. MATERIALS AND METHODS
// =============================================================================
#heading(level: 1, numbering: none)[4. Materials and Methods]

#heading(level: 2, numbering: none)[4.1. Taxon Sampling and Dual-Database Annotation Audit]

Eighteen sauropsid NCBI RefSeq reference genomes (3 Crocodylia, 3 Squamata, 4 Testudines, 8 Aves) plus _Homo sapiens_ `GRCh38.p14` (`GCF_000001405.40`) were audited programmatically via the NCBI Datasets v2 REST API and `datasets` CLI (`v18.38.0`, audited `2026-11-19`; complete accessions, assembly N50s, and annotation release dates in Table S1). For _Gallus gallus_ `GRCg7b` (`GCF_016699485.2`), we downloaded and locally parsed both NCBI RefSeq Annotation Release 106 (`Ggal_genomic.gff`) and Ensembl Release 116 (`Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.116.gtf.gz`, SHA-256 `407e3e270da9668eecc205b685493d157abbd8fa80bc698a65bfd3ed3a5f8009`) to check genome-wide odontogenic symbols and coordinate overlaps within `NC_052572.1:55,109,000–55,235,000`.

#heading(level: 2, numbering: none)[4.2. Syntenic Interval Extraction and Query Composition Controls]

For each of the 18 sauropsids, two genomic windows were retrieved via NCBI E-utilities (`efetch`): (i) the _NRG1–LPL_ enamel-cluster interval (±60 kb padding around _NRG1_ and _LPL_, or ±800 kb around the single annotated flanking gene in _D. novaehollandiae_ and _T. guttatus_; Table S2a), and (ii) the _IBSP–SPP1–SPARCL1_ osteogenic IDP control interval (±100 kb padding). Sixteen protein queries—full-length RefSeq _Alligator_ `ENAM` (`XP_019346377.1`, 1,078 aa), legacy partial GenBank `ENAM` (`QOX06214.1`, 849 aa), `AMBN` (`XP_019346386.1`, 390 aa), `AMTN` (`XP_014451634.1`, 328 aa; `KYO35092.1`, 375 aa), globular positive control `LPL` (`XP_019346320.1`, 477 aa), cross-chromosome negative control `DSPP` (`XP_059578269.1`, 783 aa), _Alligator_ `IBSP`/`SPP1`/`SPARCL1`, and equidistant _Anolis carolinensis_ `ENAM` (`XP_062825318.1`), `AMBN` (`NP_001315156.1`), `AMTN` (`NP_001315157.1`), `LPL` (`XP_062827593.1`), `IBSP` (`XP_003221343.2`), and `SPP1` (`XP_003221344.2`)—were profiled for amino-acid composition, Shannon entropy, and low-complexity fraction using NCBI `segmasker` `v1.0.0` (`-window 12 -locut 2.2 -hicut 2.5`; Table S3a).

For every query, 20 composition-preserving random permutations were generated deterministically in Python 3.11 using `numpy.random.default_rng(20261017)` (never global `np.random.seed()`). Translated searches were executed with NCBI `BLAST+` `tblastn` `v2.17.0+` @camacho2009 across all 18 enamel intervals and all 18 IDP intervals under both `-seg yes` and `-seg no` (`-evalue 10 -max_target_seqs 100`). For each species–query–mode combination, excess bitscore was defined as $max(0, "obs" - p_(95,"null"))$, where $p_(95,"null")$ is the 95th percentile of maximum bitscores across the 20 matched permutations. Null-corrected retention indices were computed as $("exc"_"ENAM" + "exc"_"AMBN" + "exc"_"AMTN") \/ "exc"_"LPL"$ (Table S3b–c).

#heading(level: 2, numbering: none)[4.3. Whole-Proteome Reciprocal BLASTx and Single-Copy Verification]

The full _Alligator mississippiensis_ RefSeq proteome (`GCF_030867095.1`, 48,040 protein sequences) was formatted with `makeblastdb` (`v2.17.0+`). Every candidate enamel hit (bitscore ≥ 28 or _E_ ≤ 0.01) in `tblastn -seg yes` was extracted in two ways: (i) clustered same-strand HSPs within 15 kb (±100 bp padding) and (ii) the single highest-scoring HSP exon (±30 bp padding) to eliminate intronic confounding. Both sets were queried against the 48,040-protein database with `blastx` `v2.17.0+` (`-evalue 10 -seg yes -max_target_seqs 25`; Table S2b). Single-copy specificity of _Gallus gallus_ _ΨENAM_ and _ΨAMBN_ within the chicken _NRG1–LPL_ interval was tested via `blastn -task dc-megablast` (`-evalue 1e-3 -dust yes`).

#heading(level: 2, numbering: none)[4.4. Level 3 Splice-Aware Alignment and Level 4 Dual-Aligner Shared Mutations]

Splice-aware protein-to-genome alignments were executed with `miniprot` `v0.18-r281` @li2023 (`-j 2 -k 4 -M 0 -n 1 -l 4 --outc=0.01 --outs=0.5`) across all 18 _NRG1–LPL_ intervals. Frameshifts, in-frame premature stop codons, and splice donor/acceptor dinucleotides (`GT-AG`, `GC-AG`, or non-canonical) were parsed directly from PAF `cs`/`cg` strings. Alignments containing introns exceeding 25,000 bp were segmented to the primary orthologous block to filter distal low-complexity chaining (specifically _C. picta_ _ENAM_'s 58.4-kb jump to `protein SON`; Table S4).

For Level 4 shared-mutation dating, coding exons of _Alligator_ _ENAM_ (Exons 6 and 8) and _AMBN_ (Exons 3, 4, and 9) were verified against RefSeq GFF3 CDS features (100% exact protein translation) and used to recruit orthologous fragments from _NRG1–LPL_ intervals via `blastn` (`-task dc-megablast` / `-task blastn`). Each exon dataset (including all three toothed crocodylians as 0-lesion positive controls) was aligned independently with *MAFFT `v7.525`* (`--localpair --maxiterate 1000`) @katoh2013 and *MUSCLE `v5.3`* (`-align`) @edgar2022. Frameshift indels (non-multiples of 3, excluding terminal boundary overhangs within 3 nt of exon ends) and in-frame PTCs (requiring all three non-gap bases aligned) were retained only when identically called in both aligners (same reference coordinate, lesion class, and carrier set). Phylogenetic congruence was scored against the established sauropsid topology (Durocryptodira, Palaeognathae, Novaeratitae, Galloanserae, Psittacopasserae, Neornithes; Table S5a–b). Auxiliary _ODAPH_/_ODAM_ gap measurements, UCSC `multiz100way` coverage @blanchette2004, `EU340348` _ΨAMEL_ calibration, and `galGal3` 32-mer recombination liftOver @elferink2010 are consolidated with SHA-256 checksums in Table S6.

#heading(level: 2, numbering: none)[4.5. Provenance and Reproducible Execution]

All bioinformatics binaries were executed inside WSL2 (`Ubuntu-26.04`, Pixi-managed environments) orchestrated from Python 3.11 (`uv` environment), with every script run, tool version (`tool_versions.txt`), parameter set, and SHA-256 file hash logged automatically to `runs/enamel_relics/provenance_run.jsonl` and `supplementary/MANIFEST_SHA256.tsv`. Figures were rendered directly from verified TSV outputs using `matplotlib` with the Okabe-Ito colorblind-safe palette and redundant shape/hatch encodings.

// =============================================================================
// 5. TRANSPARENCY & DECLARATIONS (bioRxiv Mandatory Disclosures)
// =============================================================================
#heading(level: 1, numbering: none)[Data and Code Availability]

All analysis scripts (`scripts/enamel_relics/01_audit_and_fetch.py` through `06_make_figures.py`), pinned software versions, provenance logs (`provenance_run.jsonl`), verified BibTeX references, vector figures, and consolidated supplementary tables (`Table_S1` through `Table_S6` with `MANIFEST_SHA256.tsv`) are openly available at `https://github.com/patoangolano/bio-stack` (archival release `enamel-relics-v1.0.0` permanently archived on Zenodo at #link("https://doi.org/10.5281/zenodo.24358067")[`https://doi.org/10.5281/zenodo.24358067`]). All genome assemblies and protein queries analyzed are publicly archived in NCBI RefSeq/GenBank and Ensembl Release 116 under the exact accessions listed in Table S1 and Table S3a.

#heading(level: 1, numbering: none)[Funding]

This research received no specific grant or external funding from any funding agency in the public, commercial, or not-for-profit sectors.

#heading(level: 1, numbering: none)[Competing Interests]

The author declares no competing financial or non-financial interests.

#heading(level: 1, numbering: none)[Declaration of Generative AI and AI-Assisted Technologies]

During the preparation of this work, AI coding and research assistants (Anthropic Claude / Google Gemini agentic workflows) were used to assist with bioinformatics script orchestration, automated verification of database accessions and primary bibliographic metadata, figure layout refinement, and English language editing. All research questions, analytical design choices, biological interpretations, manual inspection of alignments and control outputs, and final scientific verification were performed by the human author, who takes full responsibility for the integrity and accuracy of the published content.

// =============================================================================
// 6. REFERENCES (Auto-generated from verified BibTeX)
// =============================================================================
#heading(level: 1, numbering: none)[References]

#[
  #set par.line(numbering: none)
  #set text(size: 8.7pt)
  #set par(leading: 0.55em, spacing: 0.58em)
  #bibliography("references.bib", title: none, style: "apa")
]
