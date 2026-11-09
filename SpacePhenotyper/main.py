import SpacePhenotyper
import pandas as pd

geneFilterTCGABRCA = pd.read_csv("../Data/Real/BulkExpr/BRCA-GeneInfo.csv", index_col=0, header=0)
geneFilterTCGACOAD = pd.read_csv("../Data/Real/BulkExpr/COAD-GeneInfo.csv", index_col=0, header=0)
percentGenes = 0.10
outputRoot = "../Output/20261104/"

Synthetic_BulkExpr = "../Data/Synthetic/syntheticBulkExpr.tsv"
Synthetic_Phenotype = "../Data/Synthetic/syntheticBulkPhenotype.tsv"
Synthetic_SRTExpr = "../Data/Synthetic/syntheticSRTExpr.tsv"
Synthetic_coords = "../Data/Synthetic/syntheticSRTCoordinates.tsv"

TCGABRCA_BulkExpr = "../Data/Real/BulkExpr/TCGA-BRCA-normalized.csv.gz"
TCGABRCA_TumorPurity = "../Data/Real/BulkPheno/TumorPurity.csv"
TCGABRCA_Proliferative = "../Data/Real/BulkPheno/Proliferative.csv"
TCGABRCA_LymphocyteActivation = "../Data/Real/BulkPheno/LymphocyteActivationGSVA.csv"
TCGABRCA_Clinical = "../Data/Real/BulkPheno/BRCA-Clinical20231212.csv"

TransNeo_BulkExpr = "../Data/Real/BulkExpr/RCB-Bulk-ExpData-FPKM_filtered.csv.gz"
TransNeo_RCB = "../Data/Real/BulkPheno/ScoreRCB.csv"

TCGACOAD_BulkExpr = "../Data/Real/BulkExpr/TCGA-COAD-normalized.csv.gz"
TCGACOAD_TumorPurity = "../Data/Real/BulkPheno/COAD-TumorPurity.csv"
TCGACOAD_Proliferative = "../Data/Real/BulkPheno/COAD-Proliferative-GSVA.csv"
TCGACOAD_LymphocyteActivation = "../Data/Real/BulkPheno/COAD-LymphocyteActivation-GSVA.csv"
TCGACOAD_Clinical = "../Data/Real/BulkPheno/COAD-Clinical20231212.csv"

BRCA110_SRTExpr = "../Data/Real/SRTExpr/BRCA-1.1.0-SRT-smooth2.csv.gz"
BRCA110_coords = "../Data/Real/SRTExpr/BRCA-1.1.0-coords.csv"
BRCA130_SRTExpr = "../Data/Real/SRTExpr/BRCA-1.3.0-SRT-smooth2.csv.gz"
BRCA130_coords = "../Data/Real/SRTExpr/BRCA-1.3.0-coords.csv"
BRCA200_SRTExpr = "../Data/Real/SRTExpr/BRCA-2.0.0-SRT-smooth2.csv.gz"
BRCA200_coords = "../Data/Real/SRTExpr/BRCA-2.0.0-coords.csv"
COAD300_SRTExpr = "../Data/Real/SRTExpr/COAD-3.0.0-SRT-smooth2.csv.gz"
COAD300_coords = "../Data/Real/SRTExpr/COAD-3.0.0-coords.csv.gz"



### Synthetic data
inputData = [
	[Synthetic_BulkExpr, Synthetic_Phenotype, Synthetic_SRTExpr, Synthetic_coords],
]
SpacePhenotyper.SpacePhenotyper_IOWrapper(inputData=inputData, outputRoot=outputRoot, geneFilter=None, CVNumFolds=5, CVNumShuffles=100, percentile=percentGenes, maxNumPCs=30)

### TCGA BRCA continuous phenotypes
inputData = [
	[TCGABRCA_BulkExpr, TCGABRCA_TumorPurity, BRCA110_SRTExpr, BRCA110_coords],
	[TCGABRCA_BulkExpr, TCGABRCA_TumorPurity, BRCA130_SRTExpr, BRCA130_coords],
	[TCGABRCA_BulkExpr, TCGABRCA_TumorPurity, BRCA200_SRTExpr, BRCA200_coords],
	[TCGABRCA_BulkExpr, TCGABRCA_Proliferative, BRCA110_SRTExpr, BRCA110_coords],
	[TCGABRCA_BulkExpr, TCGABRCA_Proliferative, BRCA130_SRTExpr, BRCA130_coords],
	[TCGABRCA_BulkExpr, TCGABRCA_Proliferative, BRCA200_SRTExpr, BRCA200_coords],
	[TCGABRCA_BulkExpr, TCGABRCA_LymphocyteActivation, BRCA110_SRTExpr, BRCA110_coords],
	[TCGABRCA_BulkExpr, TCGABRCA_LymphocyteActivation, BRCA130_SRTExpr, BRCA130_coords],
	[TCGABRCA_BulkExpr, TCGABRCA_LymphocyteActivation, BRCA200_SRTExpr, BRCA200_coords],
]
SpacePhenotyper.SpacePhenotyper_IOWrapper(inputData=inputData, outputRoot=outputRoot, geneFilter=geneFilterTCGABRCA, CVNumFolds=5, CVNumShuffles=100, percentile=percentGenes, maxNumPCs=30)

### TCGA BRCA hazard
inputData = [
	[TCGABRCA_BulkExpr, TCGABRCA_Clinical, BRCA110_SRTExpr, BRCA110_coords],
	[TCGABRCA_BulkExpr, TCGABRCA_Clinical, BRCA130_SRTExpr, BRCA130_coords],
	[TCGABRCA_BulkExpr, TCGABRCA_Clinical, BRCA200_SRTExpr, BRCA200_coords],
]
SpacePhenotyper.SpacePhenotyper_IOWrapper(inputData=inputData, outputRoot=outputRoot, geneFilter=geneFilterTCGABRCA, CVNumFolds=5, CVNumShuffles=100, percentile=percentGenes, maxNumPCs=30, phenotype="hazard")

### TransNeo RCB
inputData = [
	[TransNeo_BulkExpr, TransNeo_RCB, BRCA110_SRTExpr, BRCA110_coords],
	[TransNeo_BulkExpr, TransNeo_RCB, BRCA130_SRTExpr, BRCA130_coords],
	[TransNeo_BulkExpr, TransNeo_RCB, BRCA200_SRTExpr, BRCA200_coords],
]
SpacePhenotyper.SpacePhenotyper_IOWrapper(inputData=inputData, outputRoot=outputRoot, geneFilter=None, CVNumFolds=5, CVNumShuffles=100, percentile=percentGenes, maxNumPCs=30)

### TCGA COAD continuous phenotypes
inputData = [
	[TCGACOAD_BulkExpr, TCGACOAD_TumorPurity, COAD300_SRTExpr, COAD300_coords],
	[TCGACOAD_BulkExpr, TCGACOAD_Proliferative, COAD300_SRTExpr, COAD300_coords],
	[TCGACOAD_BulkExpr, TCGACOAD_LymphocyteActivation, COAD300_SRTExpr, COAD300_coords],
]
SpacePhenotyper.SpacePhenotyper_IOWrapper(inputData=inputData, outputRoot=outputRoot, geneFilter=geneFilterTCGACOAD, CVNumFolds=5, CVNumShuffles=100, percentile=percentGenes, maxNumPCs=30)

### TCGA COAD hazard
inputData = [
	[TCGACOAD_BulkExpr, TCGACOAD_Clinical, COAD300_SRTExpr, COAD300_coords],
]
SpacePhenotyper.SpacePhenotyper_IOWrapper(inputData=inputData, outputRoot=outputRoot, geneFilter=geneFilterTCGACOAD, CVNumFolds=5, CVNumShuffles=100, percentile=percentGenes, maxNumPCs=30, phenotype="hazard")

print('Done!')
