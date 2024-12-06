This repository contains the data and code to reproduce 'Urbanisation and urban green space shape the distribution of Kestrel species at the global scale'. 
Scripts to undertake the analysis are located in the folders 'R/Habitat Association modelling', 'R/Landscape Analysis' and 'R/Occupancy modelling'.
Due to file size constraints, some of the files for the Habitat Association modelling, and all the files for the Landscape analysis are not in this repository. 
The code last worked under the following session information:
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=English_Australia.utf8  LC_CTYPE=English_Australia.utf8    LC_MONETARY=English_Australia.utf8
[4] LC_NUMERIC=C                       LC_TIME=English_Australia.utf8    

time zone: Australia/Adelaide
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] dggridR_3.1.0          sdmTMB_0.6.0           slider_0.3.2           auk_0.7.0             
 [5] tibble_3.2.1           spOccupancy_0.7.6      gt_0.11.1              gtsummary_2.0.3       
 [9] patchwork_1.3.0        performance_0.12.4     landscapemetrics_2.1.4 purrr_1.0.2           
[13] arrow_17.0.0.1         tidyr_1.3.1            terra_1.7-83           sf_1.0-19             
[17] scam_1.2-17            mgcv_1.8-42            nlme_3.1-162           stringr_1.5.1         
[21] readr_2.1.5            ranger_0.16.0          mccf1_1.1              lubridate_1.9.3       
[25] gridExtra_2.3          ggplot2_3.5.1          ebirdst_3.2022.3       dplyr_1.1.4           

loaded via a namespace (and not attached):
 [1] DBI_1.2.3          httr2_1.0.6        s2_1.1.7           rlang_1.1.4        magrittr_2.0.3    
 [6] e1071_1.7-16       compiler_4.3.1     vctrs_0.6.5        wk_0.9.4           pkgconfig_2.0.3   
[11] fastmap_1.2.0      utf8_1.2.4         tzdb_0.4.0         nloptr_2.1.1       bit_4.5.0         
[16] xfun_0.49          gert_2.1.4         jsonlite_1.8.9     collapse_2.0.17    parallel_4.3.1    
[21] R6_2.5.1           stringi_1.8.4      boot_1.3-28.1      Rcpp_1.0.13-1      assertthat_0.2.1  
[26] iterators_1.0.14   knitr_1.48         warp_0.2.1         usethis_3.0.0      gitcreds_0.1.2    
[31] Matrix_1.6-5       splines_4.3.1      timechange_0.3.0   tidyselect_1.2.1   rstudioapi_0.17.1 
[36] abind_1.4-8        doParallel_1.0.17  codetools_0.2-19   curl_6.0.0         lattice_0.21-8    
[41] withr_3.0.2        askpass_1.2.1      coda_0.19-4.1      units_0.8-5        proxy_0.4-27      
[46] xml2_1.3.6         pillar_1.9.0       spAbundance_0.2.1  KernSmooth_2.23-21 foreach_1.5.2     
[51] insight_0.20.5     generics_0.1.3     rprojroot_2.0.4    credentials_2.0.2  hms_1.1.3         
[56] munsell_0.5.1      scales_1.3.0       minqa_1.2.8        class_7.3-22       glue_1.8.0        
[61] tools_4.3.1        sys_3.4.3          lme4_1.1-35.5      RANN_2.6.2         fs_1.6.5          
[66] grid_4.3.1         gh_1.4.1           colorspace_2.1-1   cli_3.6.3          rappdirs_0.3.3    
[71] fansi_1.0.6        gtable_0.3.6       digest_0.6.37      classInt_0.4-10    farver_2.1.2      
[76] htmltools_0.5.8.1  lifecycle_1.0.4    openssl_2.2.2      bit64_4.5.2        MASS_7.3-60   
