
library(tidyverse)
library(purrr)
library(Barnard)

files <- c( 
     "raw_data/030923_DNA_Final1.csv",
     "raw_data/031623_DNA_Final1.csv",
     "raw_data/032323_DNA_Final1.csv",
     "raw_data/041323_DNA_Final1.csv"
     )
siRNA <- as_tibble(read.csv("raw_data/HatchScreen2.csv")) %>%
    reframe(siRNA, Metadata_Well)

## Make sure the Apply_Filters function is updated based on 
## the adjustments made in the Setting_Filters template
## to reflect the intensity signals of the experiment

Apply_Filters <- function() {
          df_siRNA %>%
          ungroup() %>%
          group_by(Experiment) %>%
          filter(
          Nuc_RFP <= 0.05,
         Cyto_RFP <= (med_Cyto_RFP + (4*sd_Cyto_RFP)),
         Cyto_RFP >= (0.0001),
         Nuc_DAPI <= (med_Nuc_DAPI + 3*sd_Nuc_DAPI),
         Nuc_DAPI >= (med_Nuc_DAPI - 2*sd_Nuc_DAPI),
         #Nuc_GFP <= (med_Nuc_GFP + 2*sd_Nuc_GFP),
         Nuc_GFP >= (0.0005)
     ) %>%
     group_by(Experiment, siRNA) %>%
     filter(
         Nuc_ConvexArea <= (med_Nuc_ConvexArea + (3*sd_Nuc_ConvexArea)),
         Nuc_ConvexArea >= (1000) ,
         Nuc_Area <= (med_Nuc_Area + (2*sd_Nuc_Area)),
         Nuc_Area >= (med_Nuc_Area - (2*sd_Nuc_Area)),
         Cyto_Area <= (med_Cyto_Area + (3*sd_Cyto_Area)),
         Cyto_Area >= (500)
     )
}

process_file <- function(files) {
     right_join(siRNA,
    as_tibble(read.csv(files, skip = 1)) %>%
    filter(!grepl('s31|s32|s33|s34|s35|s36', FileName_DNA)) %>%
    reframe("Experiment" = files, 
    Metadata_Well, 
    ImageNumber, 
    ObjectNumber, 
    "Cyto_Area" = AreaShape_Area.1,
    "Cyto_DAPI" = Intensity_MeanIntensity_CorrDAPI.1,
    "Cyto_GFP" = Intensity_MeanIntensity_CorrGFP.1,
    "Cyto_RFP" = Intensity_MeanIntensity_CorrRFP.1,
    "Nuc_Area" = AreaShape_Area,
    "Nuc_Orientation" = AreaShape_Orientation,
    "Nuc_ConvexArea" = AreaShape_ConvexArea,
    "Nuc_Eccentricity" = AreaShape_Eccentricity,
    "Nuc_Solidity" = AreaShape_Solidity,
    "Nuc_DAPI" = Intensity_MeanIntensity_CorrDAPI,
    "Nuc_GFP" = Intensity_MeanIntensity_CorrGFP,
    "Nuc_RFP" = Intensity_MeanIntensity_CorrRFP) %>%
    group_by(Experiment, Metadata_Well) %>%
    na.omit() %>%
    mutate(
       "med_Cyto_Area" = median(Cyto_Area),
        "med_Cyto_DAPI" = median(Cyto_DAPI),
        "med_Cyto_GFP" = median(Cyto_GFP),
        "med_Cyto_RFP" = median(Cyto_RFP),
        "med_Nuc_Area" = median(Nuc_Area),
          "med_Nuc_Orientation" = median(Nuc_Orientation),
          "med_Nuc_ConvexArea" = median(Nuc_ConvexArea),
          "med_Nuc_Eccentricity" = median(Nuc_Eccentricity),
        "med_Nuc_Solidity" = median(Nuc_Solidity),
        "med_Nuc_DAPI" = median(Nuc_DAPI),
        "med_Nuc_GFP" = median(Nuc_GFP),
        "med_Nuc_RFP" = median(Nuc_RFP),
        "sd_Cyto_Area" = sd(Cyto_Area),
        "sd_Cyto_DAPI" = sd(Cyto_DAPI),
        "sd_Cyto_GFP" = sd(Cyto_GFP),
        "sd_Cyto_RFP" = sd(Cyto_RFP),
        "sd_Nuc_Area" = sd(Nuc_Area),
          "sd_Nuc_Orientation" = sd(Nuc_Orientation),
          "sd_Nuc_ConvexArea" = sd(Nuc_ConvexArea),
          "sd_Nuc_Eccentricity" = sd(Nuc_Eccentricity),
        "sd_Nuc_Solidity" = sd(Nuc_Solidity),
        "sd_Nuc_DAPI" = sd(Nuc_DAPI),
        "sd_Nuc_GFP" = sd(Nuc_GFP),
        "sd_Nuc_RFP" = sd(Nuc_RFP)
    ) ) 
}

df_raw <- purrr::map_dfr(files, process_file)

df_siRNA <- left_join(df_raw %>%     
     group_by(Experiment, ImageNumber, Metadata_Well) %>%
     count() %>%
     filter(n >= 1,
          n <= 50) %>%
     mutate('stage_filter' = '1_50'),
     df_raw) %>%
    filter(siRNA != 'NA') %>%
    group_by(Experiment, siRNA) 

df_filtered <- Apply_Filters() %>%
     ungroup() 

df_final <- left_join(
     df_filtered %>%
          group_by(Experiment, siRNA) %>%
          mutate("Ratio_RFP" = Nuc_RFP/Cyto_RFP,
               "Ratio_GFP" = Cyto_GFP/Nuc_GFP) %>%
               ungroup(siRNA),

     df_filtered %>%
          group_by(Experiment, siRNA) %>%
          mutate("Ratio_RFP" = Nuc_RFP/Cyto_RFP,
               "Ratio_GFP" = Cyto_GFP/Nuc_GFP) %>%
          filter(grepl('CTRL', siRNA)) %>%
          ungroup() %>%
          group_by(Experiment) %>%
          mutate('Rthreshold' = (median(Ratio_RFP) + sd(Ratio_RFP)),
              'Gthreshold' = (median(Ratio_GFP) + sd(Ratio_GFP))) %>%
          select(Rthreshold, Gthreshold) %>%
          unique()
          ) %>%
     ungroup() %>%
     group_by(Experiment, siRNA) %>%
     mutate('median' = median(Ratio_RFP), 'mean' = mean(Ratio_RFP)) %>%
     mutate('ruptured' = 
          case_when(
               Ratio_RFP >= Rthreshold ~ 'y',
               Ratio_RFP < Rthreshold ~ 'n'
          ), 
          'rupturing' = 
               case_when(
               Ratio_GFP >= Gthreshold ~ 'y',
               Ratio_GFP < Gthreshold ~ 'n'
          ))

rfp_frequency_n <- df_final %>%
     ungroup() %>%
     group_by(Experiment, siRNA, ruptured) %>%
     count() %>%
     ungroup(ruptured) %>%
     reframe(Experiment, siRNA, n, 'total' = sum(n), ruptured) %>%
     filter(ruptured != 'n') %>%
     reframe(Experiment, siRNA, 'ruptured' = n, total) %>%
     mutate('rupture_frequency' = ruptured/total) %>%
     mutate('rupture_neg' = (total - ruptured))

stats_prep <- function() {
     df_final %>%
     ungroup() %>%
     group_by(siRNA, ruptured) %>%
     count() %>%
     ungroup(ruptured) %>%
     reframe(siRNA, n, 'total' = sum(n), ruptured) %>%
     filter(ruptured != 'n') %>%
     reframe(siRNA, 'ruptured' = n, total) %>%
     mutate('rupture_frequency' = ruptured/total) %>%
     mutate('rupture_neg' = (total - ruptured))
}

CTRLS <- stats_prep() %>%
     filter(grepl('CTRL', siRNA)) %>%
     reframe('NEG' = rupture_neg, 'POS' = ruptured)

rfp_frequency_pool <- stats_prep() %>%
     mutate('NEG' = CTRLS$NEG, 'POS' = CTRLS$POS) %>%
     ungroup() %>%
     rowwise() %>%
     mutate('p' = barnard.test(NEG, rupture_neg, POS, ruptured)$p.value[2]) 

df_final_RFPfrequency <- left_join(df_final, rfp_frequency_pool %>% select(siRNA, rupture_frequency))

df_final_RFPfrequency <- left_join(
          
          df_final %>%
          filter(ruptured == 'y'),

          rfp_frequency_pool %>% 
          select(siRNA, rupture_frequency)) 

stats_prep <- function() {
     df_final_RFPfrequency %>%
     ungroup() %>%
     group_by(siRNA, rupturing) %>%
     count() %>%
     ungroup(rupturing) %>%
     reframe(siRNA, rupturing, n, 'total' = sum(n)) %>%
     filter(rupturing != 'n') %>%
     reframe(siRNA, 'rupturing' = n, total) %>%
     mutate('rupturing_frequency' = rupturing/total) %>%
     mutate('neg' = (total-rupturing))
}

CTRLS <- stats_prep() %>%
     filter(grepl('CTRL', siRNA)) %>%
     reframe('NEG' = neg, 'POS' = rupturing)

gfp_frequency_pool <- stats_prep() %>%
     mutate('NEG' = CTRLS$NEG, 'POS' = CTRLS$POS) %>%
     ungroup() %>%
     rowwise() %>%
     mutate('p' = barnard.test(NEG, neg, POS, rupturing)$p.value[2]) 

gfp_frequency_n <- df_final_RFPfrequency %>%
     ungroup() %>%
     group_by(Experiment, siRNA, rupturing, rupture_frequency) %>%
     count() %>%
     ungroup(rupturing) %>%
     reframe(Experiment, siRNA, rupturing, n, 'total' = sum(n)) %>%
     filter(rupturing != 'n') %>%
     reframe(Experiment, siRNA, 'rupturing' = n, total) %>%
     mutate('rupturing_frequency' = rupturing/total) %>%
     mutate('neg' = (total-rupturing))

df_final %>% write.csv('Exported_Data/Ratio_Final.csv')

rfp_frequency_pool %>% write.csv('Exported_Data/Frequency_RFP_Pool.csv')

rfp_frequency_n %>% write.csv('Exported_Data/Frequency_RFP_N.csv')

gfp_frequency_pool %>% write.csv('Exported_Data/Frequency_GFP_Pool.csv')

gfp_frequency_n %>% write.csv('Exported_Data/Frequency_GFP_N.csv')
