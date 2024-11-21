# Detection of inflammation via pan-leukocyte biomarkers
This repository contains datasets (raw and processed) and R scripts used in the manuscript: "High-throughput, pan-leukocyte biomarkers for the detection of inflammation in human breastmilk and stool" by Dunnet, M.J., Morison, I.M., Bond, D.M., Hore, T.A.

## Datasets
Information file for datasets
-	Dunnet_SequencingFileNames.txt contains a list of all the raw and processed sequencing files available in this repository. The file name in the “Raw data file” column pairs with the file name in the “Processed data file” column.

These are the raw and processed sequencing files for various amplicon sequencing analyses.
-	Folder: Raw_Files
-	Folder: Bismark_Methylation_Extraction_Data_Files

## Scripts
For each script to run correctly, the directory input and output locations must be edited within the script itself

### Methylation_Heatmaps.R
This script produces heatmaps based on individual sequencing reads. The script takes text files produced by the bismark methylation extractor as input and converts it into a binary array for plotting. Sequencing reads run from left to right across the rows, and columns represent individual CpG positions. Red indicates a methylated CpG, and blue indicates an unmethylated CpG. 

### Generate_ROC_curve.R
This script takes text files produced by the bismark methylation extractor as input. It is deisgned to find the optimum number of methylated CpGs to call a leukocyte vs non-leukocyte; therefore, input should be two amplicon datasets comprised of  only leukocyte DNA and only non-leukocyte DNA. It will calculate a binomial regression of cell types and produce a ROC curve based on methylation thresholds for leukocyte calling. The output is plots and a tab-deliminated text file

### HOXA3_46_Leukocyte_Caller.R
This script is designed to work only with HOXA3 amplicons descirbed in the manuscript. It takes text files produced by the bismark methylation extractor as input. It will determine the number of methylated CpG sites from the forward read and classify it as leukocyte derived, non-leukocyte dervied, or unknown. Thresholds for each classifcation are based upon binomial regression from the `Generate_ROC_curve.R` script. The output is a tab-separated text file.

### MAP4K1_Leukocyte_Caller.R
This script is designed to work only with MAP4K1 amplicons descirbed in the manuscript. It takes text files produced by the bismark methylation extractor as input. It will determine the number of methylated CpG sites from the forward read and classify it as leukocyte derived, non-leukocyte dervied, or unknown. Thresholds for each classifcation are based upon binomial regression from the `Generate_ROC_curve.R` script. The output is a tab-separated text file.
