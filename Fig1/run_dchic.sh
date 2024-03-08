#!/bin/bash
RPATH=$1
INPUTFILE=$2
OUTFOLDER=$3
### FOR DETAILS PLEASE VISIT https://github.com/ay-lab/dcHiC ###

echo "
###                  STEP 1                  ###
### Step 1: dcHiC will create the raw PC files for each chromosome.
### This step requires the raw Hi-C matrix files similar to HiC-Pro 
### output described in https://nservant.github.io/HiC-Pro/RESULTS.html#intra-and-inter-chromosomal-contact-maps
### For size limits we not providing the raw matrix files. 
### We are providing the *.pc.txt files from step 1 of dcHiC generated as follows (commented, otherwise will raise an error) -

$RPATH/Rscript /data/kun/Softwares/dcHiC/demo/dcHiC_demo/scripts/dchicf.r --file $INPUTFILE --pcatype cis --dirovwt T --cthread 25 --pthread 30
### STEP 1 ENDS

"
### Decomment this part to run the step 1
$RPATH/Rscript /data/kun/Softwares/dcHiC/demo/dcHiC_demo/scripts/dchicf.r --file $INPUTFILE --pcatype cis --dirovwt T --cthread 1 --pthread 1

echo "
###                  STEP 2                  ###
### Step 2: dcHiC will select the best pc out of PC1 and PC2 for each chromosome 
### by comparing each one against GC content and gene density through correlation.
### Run the step 2 of dcHiC using the following command - 

$RPATH/Rscript /data/kun/Softwares/dcHiC/demo/dcHiC_demo/scripts/dchicf.r --file $INPUTFILE --pcatype select --dirovwt T --genome mm10
### STEP 2 ENDS

"
### Decomment this part to run the step 2
$RPATH/Rscript /data/kun/Softwares/dcHiC/demo/dcHiC_demo/scripts/dchicf.r --file $INPUTFILE --pcatype select --dirovwt T --genome hg19

echo "
###                  STEP 3                  ###
### Step 3: dcHiC will use the selected PC and qunatile normalize them.
### This will be followed by mahalanobis distance calculation and outlier detection.
### If replicates are available, dcHiC will apply IHW to adjust the significance. 
### This will create a directory "DifferentialResult" and few more subfolder. 
### The differential compartment result files will be under ./DifferentialResult/ES_NPC_CN_100Kb/fdr_result .
### differential.intra_sample_group.pcOri.bedGraph shows the significance score for all the Hi-C bins.
### differential.intra_sample_group.Filtered.pcOri.bedGraph shows the filtered (padj < 0.1) the Hi-C bins.
### The differential.intra_sample_combined.*.bedGraph file shows the compartment scores and significance 
### of each replicates, differential.intra_sample_group.*.bedGraph files are the subset of these.
### *.pcOri.* are files with original pc values i.e. before quantile normalized values.
### *.pcQnm.* are files with quantile normalized values. 
### Note: All the differential compartment analysis is performed using *.pcQnm.* files. 
### *.pcOri.* files are generated for vizualization purpose.
### Run the step 3 of dcHiC using the following command -  

$RPATH/Rscript /data/kun/Softwares/dcHiC/demo/dcHiC_demo/scripts/dchicf.r --file $INPUTFILE --pcatype analyze --dirovwt T --diffdir NT_PT_RT_20000
### STEP 3 ENDS

"
### Decomment this part to run the step 3
$RPATH/Rscript /data/kun/Softwares/dcHiC/demo/dcHiC_demo/scripts/dchicf.r --file $INPUTFILE --pcatype analyze --dirovwt T --diffdir $OUTFOLDER


echo "
###                  STEP 4                  ###
### Step 4: dcHiC will perform significant loop calling using FitHiC2 (https://github.com/ay-lab/fithic).
### This step requires the raw hic files and also time consuming, thus will be skipped.
### But at the end of this step, dcHiC will create a folder named as fithic_run under ./DifferentialResult/ES_NPC_CN_100Kb.
### We are proving the output of this step under ./DifferentialResult/ES_NPC_CN_100Kb/fithic_run/FithicResult.txt as an example. 
### Users can run any of their favorite loop caller and create the FithicResult.txt file. 
### Each row if this file shows the interacting partners (chr1_start1_chr2_start2) and each column represents if that interaction is 
### significant (1) or not (0) in that sample. We are skipping this step, otherwise it will raise an error.

$RPATH/Rscript /data/kun/Softwares/dcHiC/demo/dcHiC_demo/scripts/dchicf.r --file $2 --pcatype fithic --dirovwt T --diffdir ES_NPC_CN_100Kb --maxd 10e6 --minc 0 --fithicpath \"/<path to fithic>/fithic.py\" --pythonpath \"/<path to python 3 or later>/bin/python3.X\"
### STEP 4 ENDS

"
### Decomment this part to run the step 4
#$RPATH/Rscript ../scripts/dchicf.r --file input.ES_NPC_CN.txt --pcatype fithic --dirovwt T --diffdir ES_NPC_CN_100Kb --maxd 10e6 --minc 0 --fithicpath "/<path to fithic>/fithic.py" --pythonpath "/<path to python 3 or later>/bin/python3.X"

echo "
###                  STEP 5                  ###
### Step 5: dcHiC will perform differential interaction calling using the 
### ./DifferentialResult/ES_NPC_CN_100Kb/fithic_run/FithicResult.txt file
### This will create a file differential.intra_compartmentLoops.bedpe under 
### ./DifferentialResult/ES_NPC_CN_100Kb/fdr_result directory which shows the
### differential interactions involving diffential interactions. The step 
### requires the Hi-C matrix file. We are skipping this step, otherwise it will raise an error.

$RPATH/Rscript ../scripts/dchicf.r --file input.ES_NPC_CN.txt --pcatype dloop --dirovwt T --diffdir ES_NPC_CN_100Kb --maxd 5e6 --minc 10
### STEP 5 ENDS

"
### Decomment this part to run the step 5
#$RPATH/Rscript ../scripts/dchicf.r --file input.ES_NPC_CN.txt --pcatype dloop --dirovwt T --diffdir ES_NPC_CN_100Kb --maxd 5e6 --minc 10


echo "
###                  STEP 6                  ###
### Step 6: dcHiC will generate the standalone IGV web page.
### The step reqires "create_datauri" code provided within scripts directory.
### Please create a soft-link to your ln -s ./scripts/create_datauri ~/.local/bin/
### This will create the html page under DifferentialResult/ES_NPC_CN_100Kb/viz/vizIGV_intra/intra_igv_pcOri.html .
### Run the step 6 of dcHiC using the following command -

$RPATH/Rscript /data/kun/Softwares/dcHiC/demo/dcHiC_demo/scripts/dchicf.r --file $2 --pcatype viz --dirovwt T --diffdir NT_PT_RT_20000 --genome mm10 --pcgroup pcOri
### STEP 6 ENDS

"
### Decomment this part to run the step 6
$RPATH/Rscript /data/kun/Softwares/dcHiC/demo/dcHiC_demo/scripts/dchicf.r --file $INPUTFILE --pcatype viz --dirovwt T --diffdir $OUTFOLDER --genome hg19 --pcgroup pcOri


echo "
###                  STEP 7                  ###
### Step 7: dcHiC will perform enrichment with the genes overlapping with the 
### differential A-compartments in each sample. This step will generate a 
### geneEnrichment directory under DifferentialResult/ES_NPC_CN_100Kb .
### Under geneEnrichment directory there will sub-directories named as *_geneEnrichment.
### We recomment using the *_geneList.anchor.txt (Entrez IDs) file to perform the gene-enrichment using
### https://toppgene.cchmc.org/enrichment.jsp .

$RPATH/Rscript ../scripts/dchicf.r --file $2 --pcatype enrich --genome hg19 --diffdir NT_PT_RT_20000 --exclA F --pcscore T --region anchor --pcgroup pcOri
### STEP 7 ENDS

"
### Decomment this part to run the step 7
$RPATH/Rscript /data/kun/Softwares/dcHiC/demo/dcHiC_demo/scripts/dchicf.r --file $INPUTFILE --pcatype enrich --genome hg19 --diffdir $OUTFOLDER --exclA F --pcscore T --region anchor --pcgroup pcOri
