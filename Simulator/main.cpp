#include "DataMatrix.h"
#include "RecordSeparatedFiles.h"
#include "Simulator.h"
#include "Preprocessing.h"

int main() {
	DataMatrix SC("../Data/Real/SingleCell/krasnow_hlca_facs_counts_3cellTypesV2_TCGAgenes.tsv", '\t', FLAG_DataMatrix_HEADER_YES);
	SC.printSize("SC");
	__cpplib__io__rsf__::vec2D_Str barcodeTable = __cpplib__io__rsf__::readSeparatedRecordsFromPath("../krasnow_hlca_facs_barcodeCellTypes.tsv", '\t');
	map_Str_Str barcodeCellTypes;
	for (const i32 i : irange(1, (i32)barcodeTable.size())) {
		barcodeCellTypes[barcodeTable[i][0]] = barcodeTable[i][1];
	}
	const auto [cellTypeNames, typeIdxOfCells, nameToTypeIdx] = getColumnCellTypes(SC, barcodeCellTypes);

	// {
	// 	vec_Str geneList = __cpplib__io__rsf__::readLinesFromPath("../TCGA-LUAD_genes.tsv");
	// 	SC.keepOnlyFollowingRows(geneList);
	// 	SC.printSize("SC");

	// 	// keepOnlyFollowingCellTypes(SC, barcodeCellTypes, {"Ciliated", "B", "Airway Smooth Muscle"});
	// 	SC.printSize("SC");

	// 	SC.printToFileAsTSV("../krasnow_hlca_facs_counts_3cellTypesV3_TCGAgenes.tsv", "gene");
	// }
	
	filterOutUnderexpressedGenes(SC, typeIdxOfCells);
	SC.printSize("SC");
	filterOutMitochondrialGenes(SC);
	SC.printSize("SC");
	vec2D_i32 markerGenes = selectMarkersSharedGenes(SC, cellTypeNames, typeIdxOfCells, /*maxNumMarkersPerCellType=*/10);

	auto [Bulk, BulkPhi, BulkCellCounts, SRT, SRTPhi, SRTCellCounts, SRTCategory] = generateSyntheticData(
		SC,
		cellTypeNames,
		typeIdxOfCells,
		nameToTypeIdx.find("Ciliated")->second,
		nameToTypeIdx.find("B")->second,
		/*numSpotsRows*/ 50,
		/*numSpotsCols*/ 50,
		/*RNGSeed*/ 29,
		/*minCellsPerSpot*/ 4,
		/*maxCellsPerSpot*/ 10,
		/*numBulkSamples*/ 200,
		/*avgCellsOfEachType*/ 15000.0,
		markerGenes[ nameToTypeIdx.find("Ciliated")->second ],
		markerGenes[ nameToTypeIdx.find("B")->second ],
		/*immuneAbsCoef*/ 8.0,
		/*numThreads*/ 64
	);
	
	// DataMatrix BulkExprPhenotypeCorr = geneCorrelationWithPhenotype(
	// 	Bulk,
	// 	BulkPhi,
	// 	/*RNGseed*/ 1337,
	// 	/*numThreads*/ 64,
	// 	/*numPermutations*/ 100000
	// );
	
	Bulk.transpose();
	SRT.transpose();

	Str folder = "../Data/Synthetic/";
	system(("mkdir " + folder).c_str());
	Bulk.printToFileAsTSV(folder + "syntheticBulkExpr.tsv", "gene");
	BulkPhi.printToFileAsTSV(folder + "syntheticBulkPhenotype.tsv", "sample");
	// BulkExprPhenotypeCorr.printToFileAsTSV(folder + "syntheticBulkPhenotypeGeneCorrelation.tsv", "gene");
	BulkCellCounts.printToFileAsTSV(folder + "syntheticBulkCellCounts.tsv", "sample");
	SRT.printToFileAsTSV(folder + "syntheticSRTExpr.tsv", "gene");
	SRTPhi.printToFileAsTSV(folder + "syntheticSRTPhenotype.tsv", "sample");
	SRTCellCounts.printToFileAsTSV(folder + "syntheticSRTCellCounts.tsv", "sample");
	SRTCategory.printToFileAsTSV(folder + "syntheticSpotCategory.tsv", "sample");




	return 0;
}