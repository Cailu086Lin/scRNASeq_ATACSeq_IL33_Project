# scRNASeq_ATACSeq_IL33_Project

Analysis provided in the manuscript titled "MrgprA3 neurons drive cutaneous immunity against helminths through selective control of myeloid-derived IL-33"

This repository is meant to provide the source code for the analysis of scRNASeq and ATACSeq data provided in the above manuscript by JM Inclan-Rico, CM Napuri, C Lin, LY Hung, AA. Ferguson, Q Wu, CF Pastore, A Stephenson, UM Femoe, F Musaigwa, HL Rossi, DR Reed, T Macháček, P Horák, I Abdus-Saboor, W Luo, and DR. Herbert.


Keywords: MrgprA3, macrophages, helminths, IL-33, type 17 responses, keratinization.

Abstract
Skin employs interdependent cellular networks to facilitate barrier integrity and host immunity through ill-defined mechanisms. Herein, we demonstrate that the human pathogen Schistosoma mansoni inhibits pruritus evoked by itch-sensing neurons bearing the Mas-related G protein-coupled receptor A3 (MrgprA3). MrgprA3 neurons control IL-17+ γδ T cell expansion, epidermal hyperplasia and host resistance against S. mansoni through shaping cytokine expression in cutaneous antigen-presenting cells (APC). Activated MrgprA3 neurons downregulate interleukin 33 (IL-33), but up-regulate TNFα in macrophages and cDC2 partially through the neuropeptide calcitonin gene-related peptide (CGRP). Strikingly, stimulation with MrgprA3-derived soluble mediators or myeloid-intrinsic deletion of IL-33 basally increases chromatin accessibility at inflammatory cytokine loci, promoting IL-17/23-dependent changes to the epidermis and resistance to helminth infection. This work uncovers a previously unknown mechanism of intercellular cross-talk wherein “itch” neuron activation reshapes myeloid cytokine expression patterns to alter skin composition for cutaneous immunity against invasive pathogens.


> sessionInfo()
R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: America/New_York
tzcode source: system (glibc)

attached base packages:
[1] grid      stats4    stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
 [1] Seurat_5.1.0                SeuratObject_5.0.2          sp_2.1-4                    DESeq2_1.44.0              
 [5] tidyr_1.3.1                 GenomicAlignments_1.40.0    Rsamtools_2.20.0            Biostrings_2.72.1          
 [9] XVector_0.44.0              Rsubread_2.18.0             DOSE_3.30.5                 enrichplot_1.24.4          
[13] lattice_0.22-6              data.table_1.16.0           fgsea_1.30.0                eulerr_7.0.2               
[17] ggvenn_0.1.10               limma_3.60.4                cinaR_0.2.3                 dplyr_1.1.4                
[21] rtracklayer_1.64.0          ChIPQC_1.40.0               BiocParallel_1.38.0         DiffBind_3.14.0            
[25] SummarizedExperiment_1.34.0 Biobase_2.64.0              MatrixGenerics_1.16.0       matrixStats_1.3.0          
[29] GenomicRanges_1.56.1        GenomeInfoDb_1.40.1         IRanges_2.38.1              S4Vectors_0.42.1           
[33] BiocGenerics_0.50.0         ggplot2_3.5.1              

loaded via a namespace (and not attached):
  [1] spatstat.sparse_3.1-0                     fs_1.6.4                                 
  [3] bitops_1.0-8                              httr_1.4.7                               
  [5] RColorBrewer_1.1-3                        TxDb.Hsapiens.UCSC.hg18.knownGene_3.2.2  
  [7] numDeriv_2016.8-1.1                       sctransform_0.4.1                        
  [9] tools_4.4.1                               utf8_1.2.4                               
 [11] R6_2.5.1                                  uwot_0.2.2                               
 [13] lazyeval_0.2.2                            apeglm_1.26.1                            
 [15] withr_3.0.1                               gridExtra_2.3                            
 [17] Nozzle.R1_1.1-1.1                         progressr_0.14.0                         
 [19] cli_3.6.3                                 spatstat.explore_3.3-2                   
 [21] fastDummies_1.7.4                         scatterpie_0.2.4                         
 [23] spatstat.data_3.1-2                       SQUAREM_2021.1                           
 [25] mvtnorm_1.3-1                             pbapply_1.7-2                            
 [27] ggridges_0.5.6                            mixsqp_0.3-54                            
 [29] yulab.utils_0.1.7                         R.utils_2.12.3                           
 [31] parallelly_1.38.0                         BSgenome_1.72.0                          
 [33] invgamma_1.1                              bbmle_1.0.25.1                           
 [35] rstudioapi_0.16.0                         RSQLite_2.3.7                            
 [37] generics_0.1.3                            gridGraphics_0.5-1                       
 [39] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2   BiocIO_1.14.0                            
 [41] hwriter_1.3.2.1                           spatstat.random_3.3-1                    
 [43] ica_1.0-3                                 gtools_3.9.5                             
 [45] GO.db_3.19.1                              Matrix_1.7-0                             
 [47] interp_1.1-6                              fansi_1.0.6                              
 [49] abind_1.4-5                               R.methodsS3_1.8.2                        
 [51] lifecycle_1.0.4                           yaml_2.3.10                              
 [53] edgeR_4.2.1                               Rtsne_0.17                               
 [55] gplots_3.1.3.1                            qvalue_2.36.0                            
 [57] SparseArray_1.4.8                         blob_1.2.4                               
 [59] promises_1.3.0                            crayon_1.5.3                             
 [61] bdsmatrix_1.3-7                           pwalign_1.0.0                            
 [63] miniUI_0.1.1.1                            cowplot_1.1.3                            
 [65] GenomicFeatures_1.56.0                    KEGGREST_1.44.1                          
 [67] pillar_1.9.0                              rjson_0.2.22                             
 [69] systemPipeR_2.10.0                        future.apply_1.11.2                      
 [71] codetools_0.2-20                          leiden_0.4.3.1                           
 [73] fastmatch_1.1-4                           glue_1.7.0                               
 [75] ShortRead_1.62.0                          spatstat.univar_3.0-1                    
 [77] ggfun_0.1.6                               GreyListChIP_1.36.0                      
 [79] vctrs_0.6.5                               png_0.1-8                                
 [81] treeio_1.28.0                             spam_2.10-0                              
 [83] gtable_0.3.5                              amap_0.8-19                              
 [85] emdbook_1.3.13                            cachem_1.1.0                             
 [87] chipseq_1.54.0                            mime_0.12                                
 [89] TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2   S4Arrays_1.4.1                           
 [91] tidygraph_1.3.1                           coda_0.19-4.1                            
 [93] survival_3.7-0                            statmod_1.5.0                            
 [95] fitdistrplus_1.2-1                        ROCR_1.0-11                              
 [97] nlme_3.1-165                              ggtree_3.12.0                            
 [99] bit64_4.0.5                               RcppAnnoy_0.0.22                         
[101] irlba_2.3.5.1                             KernSmooth_2.23-24                       
[103] colorspace_2.1-1                          DBI_1.2.3                                
[105] tidyselect_1.2.1                          TxDb.Celegans.UCSC.ce6.ensGene_3.2.2     
[107] bit_4.0.5                                 compiler_4.4.1                           
[109] curl_5.2.2                                httr2_1.0.3                              
[111] plotly_4.10.4                             TxDb.Mmusculus.UCSC.mm10.knownGene_3.10.0
[113] DelayedArray_0.30.1                       shadowtext_0.1.4                         
[115] scales_1.3.0                              caTools_1.18.2                           
[117] lmtest_0.9-40                             rappdirs_0.3.3                           
[119] goftest_1.2-3                             stringr_1.5.1                            
[121] digest_0.6.37                             spatstat.utils_3.1-0                     
[123] htmltools_0.5.8.1                         pkgconfig_2.0.3                          
[125] jpeg_0.1-10                               fastmap_1.2.0                            
[127] rlang_1.1.4                               htmlwidgets_1.6.4                        
[129] UCSC.utils_1.0.0                          shiny_1.9.1                              
[131] farver_2.1.2                              zoo_1.8-12                               
[133] jsonlite_1.8.8                            GOSemSim_2.30.2                          
[135] R.oo_1.26.0                               RCurl_1.98-1.16                          
[137] magrittr_2.0.3                            GenomeInfoDbData_1.2.12                  
[139] ggplotify_0.1.2                           dotCall64_1.1-1                          
[141] patchwork_1.2.0                           munsell_0.5.1                            
[143] Rcpp_1.0.13                               reticulate_1.39.0                        
[145] ape_5.8                                   viridis_0.6.5                            
[147] stringi_1.8.4                             ggraph_2.2.1                             
[149] zlibbioc_1.50.0                           MASS_7.3-61                              
[151] plyr_1.8.9                                parallel_4.4.1                           
[153] listenv_0.9.1                             ggrepel_0.9.5                            
[155] deldir_2.0-4                              graphlayouts_1.1.1                       
[157] splines_4.4.1                             tensor_1.5                               
[159] locfit_1.5-9.10                           igraph_2.0.3                             
[161] spatstat.geom_3.3-2                       RcppHNSW_0.6.0                           
[163] reshape2_1.4.4                            XML_3.99-0.17                            
[165] latticeExtra_0.6-30                       renv_1.0.7                               
[167] BiocManager_1.30.25                       httpuv_1.6.15                            
[169] tweenr_2.0.3                              RANN_2.6.2                               
[171] purrr_1.0.2                               polyclip_1.10-7                          
[173] scattermore_1.2                           future_1.34.0                            
[175] ashr_2.2-63                               ggforce_0.4.2                            
[177] xtable_1.8-4                              TxDb.Rnorvegicus.UCSC.rn4.ensGene_3.2.2  
[179] restfulr_0.0.15                           RSpectra_0.16-2                          
[181] tidytree_0.4.6                            later_1.3.2                              
[183] viridisLite_0.4.2                         truncnorm_1.0-9                          
[185] tibble_3.2.1                              aplot_0.2.3                              
[187] memoise_2.0.1                             AnnotationDbi_1.66.0                     
[189] cluster_2.1.6                             globals_0.16.3                           
[191] TxDb.Dmelanogaster.UCSC.dm3.ensGene_3.2.2
