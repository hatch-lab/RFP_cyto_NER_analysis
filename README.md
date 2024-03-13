---
title: Code/Data Map
author: AG
date: 041423
output: 
  html_document:
    keep_md: TRUE
---
# RFP_cyto_NER_analysis
Code used to analyze screen images in Gunn et al., Sci Reports, 2024

## raw data

The raw data directory contains the raw CellProfiler .csv output for each data. 

Row 1 headers indicate whether the information is relevant to the individual object, or the entire image. The information pulled for processing these data largely came from the Row 2 headers, as described below.

__Useful column name descriptions__

The following columns are self descriptive, but the structure (in order) below may help identify the measurements that are nucleus-specific vs cytoplasm-specific.

* Metadata_Well
* FileName_* 
* Nucleus (DNA_Final1 from CellProfiler)
     * AreaShape_Area
     * AreaShape_ConvexArea
     * AreaShape_Eccentricity
     * AreaShape_Solidity
* Nucleus (DNA_Object_Measure from CellProfiler)
     * Intensity_MeanIntensity_CorrDAPI
     * Intensity_MeanIntensity_CorrGFP
     * Intensity_MeanIntensity_CorrRFP
* Cytoplasm (Cyto_Object_Image from CellProfiler)
     * AreaShape_Area
* Cytoplasm (Cyto_Object_Measure from CellProfiler)
     * Intensity_MeanIntensity_CorrDAPI
     * Intensity_MeanIntensity_CorrGFP
     * Intensity_MeanIntensity_CorrRFP

## code

For each code chunk described below, the corresponding data output has been stored in the processed_data folder of this directory. For these scripts to run smoothly, they should be run in order.

__CellProfiler Pipeline__

CellProfiler v4.2.4

The CellProfiler pipeline is included in the code folder of this directory as `Screen.cppipe`. 

Within the Names and Types module, the single image locations should be adjusted to reflect the filepath to a blank image with the same microscope settings to be used for background subtraction and flatfield correction.

There are two sections within which minor adjustments may need to be made based on the relative intensity of the images for each iteration.

* Identify Primary Objects
    Size of object
    Threshold smoothing scale
    Threshold correction factor

* Identify Secondary Objects
    
    Threshold smoothing scale
    Threshold correction factor

__Setting Filters Template__

R v4.2.1

The purpose of this template is not to produce data output, but provide an opportunity to visualize the size, shape, and fluorescence intensity information from select wells of the CellProfiler output. The filters can be adjusted as described in the document based on experimental variation.

After using this tool, it is important to copy the filter settings over to Apply_Filters function in Export_Data.r 

__Export Data__

R v4.2.1

This script cleans, filters, processes into RFP and GFP compartment ratios, and makes rupture designations followed by calculating statistics on the pooled rupture frequency designations using Barnard's test. At the end of the script, the processed data is saved to an Exported Data folder in the working directory for future use.

Input:

raw_data, filters from Setting_Filters_Template.rmd

Output:

`Exported_Data/Frequency_GFP_N.csv` Frequency of cyto:nuc GFP over threshold by each experiment

`Exported_Data/Frequency_RFP_N.csv` Frequency of nuc:cyto RFP over threshold by each experiment

`Exported_Data/Frequency_GFP_Pool.csv` Frequency of cyto:nuc GFP over threshold, pooled with stats

`Exported_Data/Frequency_RFP_Pool.csv` Frequency of nuc:cyto RFP over threshold, pooled with stats

`Exported_Data/Ratio_Final.csv` Includes the above information along with morphological measurements and the  RFP/GFP ratios by object

__Statistics__

R v4.2.1

Although the statistics by Barnard's test have been calculated in `Export_Data.r`, the initial results do not include adjustments for multiple testing. In this script, both False Discovery Rate and Bonferroni adjustments for multiple correction are calculated and compared.

Input:

Exported_Data

Output:

`statistics.csv`

__Preliminary Plots__

R v4.2.1

This script takes the processed information and creates preliminary visualizations for comparison.

Input:

Exported_Data, `statistics.csv` (in processed_data)

Output:

visualizations with the option of saving by adjusting the ggsave commands after each plot.

__Morphology stats__

R v4.2.1

Although the morphological information for each nucleus was originally pulled with the objective of quality control, this script allows for analysis of the medians/statistics for those data.

Input: 

Exported_Data

Output:

`morphology_stats.csv`

__manual analysis__

R v4.2.1

As a method of quality control, it is important to take a sample of the total exported data and manually assess the rupture and morphology calls. This script creates an empty spreadsheet that can be used to compare CellProfiler object numbers directly to the images of those objects in the raw data output. Following manual assessment of each object, re-save the .csv as `extremes_manual_completed.csv` to re-import into the script, which will complete a quality assessment using the precision, recall, and False Discovery Rate.

Input:
Exported_data, `extremes_manual.csv` (in processed_data)

Output:
`extremes_manual.csv`, `extremes_manual_completed.csv`

__Enrichment__

R v4.2.1

This pilot screen is composed of equal parts biased genes selected by the authors ("targeted"), and genes selected by a random number generator ("untargeted"). Because the screen gene set is, accordingly, biased toward the categories chosen by the authors, hypergeometric enrichment analysis was completed _within_ the gene set rather based on comparison to the entire database. 

For this analysis, the Broad Institute C5 database was joined to the screen gene set prior to analysis. All enrichment is carried out internally, looking at the genes with statistically significant rupture frequencies relative to the screen.

Input:

Exported_Data, `statistics.csv` (in processed_data)

Output:

`enrichment_results/enrich.C5_RFP_bonf.csv`, `enrichment_results/enrichment_bonf.csv`, `results_targeted_RFP_bonf.csv`
