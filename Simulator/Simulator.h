#ifndef __H__EIGENSTUFF__SIMULATOR__
#define __H__EIGENSTUFF__SIMULATOR__



#include "DataMatrix.h"
#include <cassert>
#include <random>
#include <tuple>
#include <thread>
#include <mutex>

template <class T>
T square(const T& x) { return x*x; }



f64 distSquared(const i32 cR, const i32 cC, const i32 i, const i32 j) {
	return square(cR - i) + square(cC - j);
}



i32 sumAllElements(const vec_i32& M) {
	i32 r = 0;
	for (const i32 val : M) r += val;
	return r;
}



void singleThreadComputePhenotypeInBulk(
	const i32 startingJobIdx,
	const i32 numBulkSamples,
	const i32 numThreads,
	const vec_i32& numCellsABulk,
	const vec_i32& numCellsBBulk,
	const vec2D_i32& sampledCellIndices,
	DataMatrix& BulkPhi,
	const vec_i32& geneIndicesA,
	const vec_i32& geneIndicesB,
	const DataMatrix& scExp,
	const f64& immuneAbsCoeff,
	const map_i32_f64& geneExpAvg,
	std::mutex& mtxCount,
	i32& printCounter
) {
	f64 avgA = 0;
	f64 avgB = 0;
	for (const i32 gIdx : geneIndicesA) {
		avgA += geneExpAvg.find(gIdx)->second;
	}
	// avgA /= (f64)geneIndicesA.size();
	for (const i32 gIdx : geneIndicesB) {
		avgB += geneExpAvg.find(gIdx)->second;
	}
	// avgB /= (f64)geneIndicesB.size();
	
	for (i32 sIdx = startingJobIdx; sIdx < numBulkSamples; sIdx += numThreads) {
		{
			std::lock_guard<std::mutex> lock(mtxCount);
			fprintf(stderr, "\r%d / %d", printCounter+1, numBulkSamples);
			++printCounter;
		}
		
		
		i32 arrIdx = 0;
		i32 expA = 0;
		i32 expB = 0;
		for (i32 _ = numCellsABulk[sIdx]; _--; ) {
			const i32 cellIdx = sampledCellIndices[sIdx][arrIdx++];
			for (const i32 gIdx : geneIndicesA) {
				if (scExp[gIdx, cellIdx] >= 1.5*geneExpAvg.find(gIdx)->second) {
					expA += scExp[gIdx, cellIdx];
				}
			}
		}
		for (i32 _ = numCellsBBulk[sIdx]; _--; ) {
			const i32 cellIdx = sampledCellIndices[sIdx][arrIdx++];
			for (const i32 gIdx : geneIndicesB) {
				if (scExp[gIdx, cellIdx] >= geneExpAvg.find(gIdx)->second) {
					expB += scExp[gIdx, cellIdx];
				}
			}
		}
		// if (expA < numCellsABulk[sIdx] * avgA/4.0) expA = 0;
		// if (expB < numCellsBBulk[sIdx] * avgB/4.0) expB = 0;
		
		BulkPhi[sIdx, 0] = expA - immuneAbsCoeff*expB;
		// if (BulkPhi[sIdx, 0] > 0) BulkPhi[sIdx, 0] = log10(BulkPhi[sIdx, 0]);
		// else if (BulkPhi[sIdx, 0] < 0) BulkPhi[sIdx, 0] = -log10(-BulkPhi[sIdx, 0]);
	}
}



void singleThreadComputeBulkPerGene(
	const i32 startingGeneIdx,
	const i32 numThreads,
	const vec_i32& numCellsABulk,
	const vec_i32& numCellsBBulk,
	const vec_i32& numCellsOtherBulk,
	const vec2D_i32& sampledCellIndices,
	DataMatrix& Bulk,
	const DataMatrix& scExp,
	std::mutex& mtxCount,
	i32& printCounter
) {
	const i32 numGenes = Bulk.numRows;
	const i32 numSamples = Bulk.numCols;
	for (i32 gIdx = startingGeneIdx; gIdx < numGenes; gIdx += numThreads) {
		{
			std::lock_guard<std::mutex> lock(mtxCount);
			fprintf(stderr, "\r%d / %d", printCounter+1, numGenes);
			++printCounter;
		}

		for (const i32 sIdx : irange(0, numSamples)) {
			i32 arrIdx = 0;
			
			for (i32 _ = numCellsABulk[sIdx]; _--; ) {
				// const i32 cellIdx = sampledCellIndices[sIdx].at(arrIdx++);
				const i32 cellIdx = sampledCellIndices[sIdx][arrIdx++];
				Bulk[gIdx, sIdx] += scExp[gIdx, cellIdx];
			}
			
			for (i32 _ = numCellsBBulk[sIdx]; _--; ) {
				// const i32 cellIdx = sampledCellIndices[sIdx].at(arrIdx++);
				const i32 cellIdx = sampledCellIndices[sIdx][arrIdx++];
				Bulk[gIdx, sIdx] += scExp[gIdx, cellIdx];
			}
			
			for (i32 _ = numCellsOtherBulk[sIdx]; _--; ) {
				// const i32 cellIdx = sampledCellIndices[sIdx].at(arrIdx++);
				const i32 cellIdx = sampledCellIndices[sIdx][arrIdx++];
				Bulk[gIdx, sIdx] += scExp[gIdx, cellIdx];
			}
		}
	}
}







std::tuple<DataMatrix, DataMatrix, DataMatrix, DataMatrix, DataMatrix, DataMatrix, DataMatrix> generateSyntheticData(
	const DataMatrix& scExp,
	const vec_Str& cellTypeNames,
	const vec_i32& typeIdxOfCells,
	const i32 typeAIdx,
	const i32 typeBIdx,
	const i32 numSpotRows,
	const i32 numSpotCols,
	const u32 RNGSeed,
	const i32 minCellsPerSpot,
	const i32 maxCellsPerSpot,
	const i32 numBulkSamples,
	const f64& avgCellsOfEachType,
	const vec_i32& geneIndicesA,
	const vec_i32& geneIndicesB,
	const f64& immuneAbsCoeff,
	const i32 numThreads
) {
	const vec_Str& geneNames = scExp.rowNames;
	const i32 numGenes = scExp.numRows;
	assert(numGenes > 0);
	const i32 numCells = scExp.numCols;
	assert(numCells > 0);
	assert((i32)typeIdxOfCells.size() == numCells);
	assert(numSpotRows > 0);
	assert(numSpotCols > 0);
	const i32 numSpots = numSpotRows * numSpotCols;
	assert(minCellsPerSpot > 0);
	assert(maxCellsPerSpot >= minCellsPerSpot);
	const set_i32 setCellTypes(typeIdxOfCells.begin(), typeIdxOfCells.end());
	assert(setCellTypes.find(typeAIdx) != setCellTypes.end());
	assert(setCellTypes.find(typeBIdx) != setCellTypes.end());
	const set_Str setGeneNames(geneNames.begin(), geneNames.end());
	assert(setGeneNames.size() == geneNames.size());

	vec_Str sampleNames;
	vec_Str spotNames;
	for (const i32 i : irange(0, numBulkSamples)) sampleNames.push_back("BS" + std::to_string(i));
	for (const i32 i : irange(0, numSpots)) spotNames.push_back("SS" + std::to_string(i));

	// Construct output data matrices
	DataMatrix Bulk(geneNames, sampleNames);
	DataMatrix BulkPhi(sampleNames, {"PhenotypeQuantity"});
	DataMatrix BulkCellCounts(sampleNames, cellTypeNames);
	DataMatrix SRT(geneNames, spotNames);
	DataMatrix SRTPhi(spotNames, {"PhenotypeQuantity"});
	DataMatrix SRTCellCounts(spotNames, cellTypeNames);
	DataMatrix SRTCategory(spotNames, {"SpotCategory"});

	vec_i32 typeACells;
	vec_i32 typeBCells;
	vec_i32 otherCells;
	const i32 cIdxRow = 0;
	const i32 cIdxCol = 0;
	const f64 distCore = numSpotRows/2.0;
	const f64 distInnerBelt = distCore + numSpotRows/4.0;
	const f64 distOuterBelt = distInnerBelt + numSpotRows/4.0;
	// vec_i32 spotType(numSpots);
	vec_i32 numCellsASRT(numSpots);
	vec_i32 numCellsBSRT(numSpots);
	vec_i32 numCellsOtherSRT(numSpots);
	vec_i32 numCellsABulk(numBulkSamples);
	vec_i32 numCellsBBulk(numBulkSamples);
	vec_i32 numCellsOtherBulk(numBulkSamples);
	std::mt19937 RNG(RNGSeed);
	std::uniform_int_distribution<i32> distNumCells(minCellsPerSpot, maxCellsPerSpot);


	// Store cell indices per category
	i32 otherTypeIdx = -1;
	for (const i32 cIdx : irange(0, numCells)) {
		if (typeIdxOfCells[cIdx] == typeAIdx) typeACells.push_back(cIdx);
		else if (typeIdxOfCells[cIdx] == typeBIdx) typeBCells.push_back(cIdx);
		else {
			otherCells.push_back(cIdx);
			otherTypeIdx = typeIdxOfCells[cIdx];
		}
	}
	std::uniform_int_distribution<i32> distTypeA(0, (i32)typeACells.size()-1);
	std::uniform_int_distribution<i32> distTypeB(0, (i32)typeBCells.size()-1);
	std::uniform_int_distribution<i32> distTypeOther(0, (i32)otherCells.size()-1);
	
	
	// Get average count of each marker in its cell type
	map_i32_f64 geneExpAvg;
	for (const i32 gIdx : geneIndicesA) {
		f64 avgCnt = 0;
		for (const i32 cellIdx : typeACells) {
			avgCnt += scExp[gIdx, cellIdx];
		}
		avgCnt /= (f64)typeACells.size();
		geneExpAvg[gIdx] = avgCnt;
	}
	for (const i32 gIdx : geneIndicesB) {
		f64 avgCnt = 0;
		for (const i32 cellIdx : typeBCells) {
			avgCnt += scExp[gIdx, cellIdx];
		}
		avgCnt /= (f64)typeBCells.size();
		geneExpAvg[gIdx] = avgCnt;
	}


	// Determine spatial cell type category distribution
	for (const i32 rIdx : irange(0, numSpotRows)) {
		for (const i32 cIdx : irange(0, numSpotCols)) {
			const i32 sIdx = rIdx * numSpotCols + cIdx;
			const i32 totalNumCells = distNumCells(RNG);
			if (distSquared(cIdxRow, cIdxCol, rIdx, cIdx) <= square(distCore)) {
				// inside the core
				SRTCategory[sIdx, 0] = 3;
				numCellsASRT[sIdx] = totalNumCells * 0.8;
				numCellsBSRT[sIdx] = 0;
			}
			else if (distSquared(cIdxRow, cIdxCol, rIdx, cIdx) <= square(distInnerBelt)) {
				// inside the inner belt
				SRTCategory[sIdx, 0] = 2;
				numCellsASRT[sIdx] = totalNumCells * 0.4;
				numCellsBSRT[sIdx] = totalNumCells * 0.4;
			}
			else if (distSquared(cIdxRow, cIdxCol, rIdx, cIdx) <= square(distOuterBelt)) {
				// inside the outer belt
				SRTCategory[sIdx, 0] = 1;
				numCellsASRT[sIdx] = 0;
				numCellsBSRT[sIdx] = totalNumCells * 0.8;
			}
			else {
				// neutral cells
				SRTCategory[sIdx, 0] = 0;
				numCellsASRT[sIdx] = 0;
				numCellsBSRT[sIdx] = 0;
			}
			numCellsOtherSRT[sIdx] = totalNumCells -numCellsASRT[sIdx] -numCellsBSRT[sIdx];
			SRTCellCounts[sIdx, typeAIdx] = numCellsASRT[sIdx];
			SRTCellCounts[sIdx, typeBIdx] = numCellsBSRT[sIdx];
			SRTCellCounts[sIdx, otherTypeIdx] = numCellsOtherSRT[sIdx];
		}
	}


	const f64 cellCountMean = avgCellsOfEachType;
	const f64 exponentialDistLambda = 1.0L/cellCountMean;
	std::exponential_distribution<f64> distCellQuantity(exponentialDistLambda);
	// Determine bulk cell type category distribution
	for (const i32 sIdx : irange(0, numBulkSamples)) {
		numCellsABulk[sIdx] = distCellQuantity(RNG);
		numCellsBBulk[sIdx] = distCellQuantity(RNG);
		numCellsOtherBulk[sIdx] = distCellQuantity(RNG);
		// fprintf(stderr, "%d %d %d\n", numCellsABulk[sIdx], numCellsBBulk[sIdx], numCellsOtherBulk[sIdx]);
		BulkCellCounts[sIdx, typeAIdx] = numCellsABulk[sIdx];
		BulkCellCounts[sIdx, typeBIdx] = numCellsBBulk[sIdx];
		BulkCellCounts[sIdx, otherTypeIdx] = numCellsOtherBulk[sIdx];
	}


	fprintf(stderr, "Sampling gene expression for spots...\n");
	// Sample expression for spots
	f64 avgA = 0;
	f64 avgB = 0;
	for (const i32 gIdx : geneIndicesA) {
		avgA += geneExpAvg.find(gIdx)->second;
	}
	// avgA /= (f64)geneIndicesA.size();
	for (const i32 gIdx : geneIndicesB) {
		avgB += geneExpAvg.find(gIdx)->second;
	}
	fprintf(stderr, "%d\n", __LINE__);
	// avgB /= (f64)geneIndicesB.size();
	for (const i32 sIdx : irange(0, numSpots)) {
		assert(numCellsASRT[sIdx] >= 0);
		assert(numCellsBSRT[sIdx] >= 0);
		assert(numCellsOtherSRT[sIdx] >= 0);
		
		i32 expA = 0;
		i32 expB = 0;
		for (i32 _ = numCellsASRT[sIdx]; _--; ) {
			// const i32 cellIdx = typeACells.at( distTypeA(RNG) );
			const i32 cellIdx = typeACells[ distTypeA(RNG) ];
			for (const i32 gIdx : irange(0, numGenes)) {
				SRT[gIdx, sIdx] += scExp[gIdx, cellIdx];
			}
			for (const i32 gIdx : geneIndicesA) {
				if (scExp[gIdx, cellIdx] >= 1.5*geneExpAvg.find(gIdx)->second) {
					expA += scExp[gIdx, cellIdx];
				}
			}
		}
		for (i32 _ = numCellsBSRT[sIdx]; _--; ) {
			// const i32 cellIdx = typeBCells.at( distTypeB(RNG) );
			const i32 cellIdx = typeBCells[ distTypeB(RNG) ];
			for (const i32 gIdx : irange(0, numGenes)) {
				SRT[gIdx, sIdx] += scExp[gIdx, cellIdx];
			}
			for (const i32 gIdx : geneIndicesB) {
				if (scExp[gIdx, cellIdx] >= geneExpAvg.find(gIdx)->second) {
					expB += scExp[gIdx, cellIdx];
				}
			}
		}
		for (i32 _ = numCellsOtherSRT[sIdx]; _--; ) {
			// const i32 cellIdx = otherCells.at( distTypeOther(RNG) );
			const i32 cellIdx = otherCells[ distTypeOther(RNG) ];
			for (const i32 gIdx : irange(0, numGenes)) {
				SRT[gIdx, sIdx] += scExp[gIdx, cellIdx];
			}
		}
		// if (expA < numCellsASRT[sIdx] * avgA/4.0) expA = 0;
		// if (expB < numCellsBSRT[sIdx] * avgB/4.0) expB = 0;
		SRTPhi[sIdx, 0] = expA - immuneAbsCoeff*expB;
		// if (SRTPhi[sIdx, 0] > 0) SRTPhi[sIdx, 0] = log10(SRTPhi[sIdx, 0]);
		// else if (SRTPhi[sIdx, 0] < 0) SRTPhi[sIdx, 0] = -log10(-SRTPhi[sIdx, 0]);
	}


	fprintf(stderr, "Sampling cell indices for bulk samples...\n");
	vec2D_i32 sampledCellIndices(numBulkSamples);
	for (const i32 sIdx : irange(0, numBulkSamples)) {
		assert(numCellsABulk[sIdx] >= 0);
		assert(numCellsBBulk[sIdx] >= 0);
		assert(numCellsOtherBulk[sIdx] >= 0);
		const i32 numCellsTotal = numCellsABulk[sIdx] + numCellsBBulk[sIdx] + numCellsOtherBulk[sIdx];
		sampledCellIndices[sIdx].reserve(numCellsTotal);
		
		for (i32 _ = numCellsABulk[sIdx]; _--; ) {
			// const i32 cellIdx = typeACells.at( distTypeA(RNG) );
			const i32 cellIdx = typeACells[ distTypeA(RNG) ];
			sampledCellIndices[sIdx].push_back(cellIdx);
		}
		for (i32 _ = numCellsBBulk[sIdx]; _--; ) {
			// const i32 cellIdx = typeBCells.at( distTypeB(RNG) );
			const i32 cellIdx = typeBCells[ distTypeB(RNG) ];
			sampledCellIndices[sIdx].push_back(cellIdx);
		}
		for (i32 _ = numCellsOtherBulk[sIdx]; _--; ) {
			// const i32 cellIdx = otherCells.at( distTypeOther(RNG) );
			const i32 cellIdx = otherCells[ distTypeOther(RNG) ];
			sampledCellIndices[sIdx].push_back(cellIdx);
		}
	}
	
	
	fprintf(stderr, "Calculating gene expression for bulk samples...\n");
	assert(numThreads >= 1);
	std::vector< std::thread > threads(numThreads);
	i32 mutexedCounter = 0;
	std::mutex mtxCount;
	for (const i32 tIdx : irange(0, numThreads)) {
		threads[tIdx] = std::thread(
			singleThreadComputeBulkPerGene,
			/*startingGeneIdx*/ tIdx,
			/*numThreads*/ numThreads,
			std::ref(numCellsABulk),
			std::ref(numCellsBBulk),
			std::ref(numCellsOtherBulk),
			std::ref(sampledCellIndices),
			std::ref(Bulk),
			std::ref(scExp),
			std::ref(mtxCount),
			std::ref(mutexedCounter)
		);
	}
	for (const i32 tIdx : irange(0, numThreads)) {
		threads[tIdx].join();
	}
	fprintf(stderr, "\n");
	
	
	fprintf(stderr, "Calculating phenotype quantities for bulk samples...\n");
	mutexedCounter = 0;
	for (const i32 tIdx : irange(0, numThreads)) {
		threads[tIdx] = std::thread(
			singleThreadComputePhenotypeInBulk,
			/*startingGeneIdx*/ tIdx,
			numBulkSamples,
			/*numThreads*/ numThreads,
			std::ref(numCellsABulk),
			std::ref(numCellsBBulk),
			std::ref(sampledCellIndices),
			std::ref(BulkPhi),
			std::ref(geneIndicesA),
			std::ref(geneIndicesB),
			std::ref(scExp),
			std::ref(immuneAbsCoeff),
			std::ref(geneExpAvg),
			std::ref(mtxCount),
			std::ref(mutexedCounter)
		);
	}
	for (const i32 tIdx : irange(0, numThreads)) {
		threads[tIdx].join();
	}
	fprintf(stderr, "\n");
	
	
	// for (const i32 sIdx : irange(0, SRTPhi.numRows)) {
	// 	if (SRTPhi[sIdx, 0] > 0) SRTPhi[sIdx, 0] = 1;
	// 	else if (SRTPhi[sIdx, 0] < 0) SRTPhi[sIdx, 0] = -1;
	// }
	// for (const i32 sIdx : irange(0, BulkPhi.numRows)) {
	// 	if (BulkPhi[sIdx, 0] > 0) BulkPhi[sIdx, 0] = 1;
	// 	else if (BulkPhi[sIdx, 0] < 0) BulkPhi[sIdx, 0] = -1;
	// }

	return {Bulk, BulkPhi, BulkCellCounts, SRT, SRTPhi, SRTCellCounts, SRTCategory};
}







// std::tuple<DataMatrix, DataMatrix, DataMatrix, DataMatrix, DataMatrix, DataMatrix> generateSyntheticDataDumb(
// 	const DataMatrix& scExp,
// 	const vec_Str& cellTypeNames,
// 	const vec_i32& typeIdxOfCells,
// 	const i32 typeAIdx,
// 	const i32 typeBIdx,
// 	const i32 numSpotRows,
// 	const i32 numSpotCols,
// 	const u32 RNGSeed,
// 	const i32 minCellsPerSpot,
// 	const i32 maxCellsPerSpot,
// 	const i32 numBulkSamples,
// 	const f64& avgCellsOfEachType,
// 	const vec_i32& geneIndicesA,
// 	const vec_i32& geneIndicesB,
// 	const i32 numThreads
// ) {
// 	const f64 immuneAbsCoeff = 1;
	
// 	const vec_Str& geneNames = scExp.rowNames;
// 	const i32 numGenes = scExp.numRows;
// 	assert(numGenes > 0);
// 	const i32 numCells = scExp.numCols;
// 	assert(numCells > 0);
// 	assert((i32)typeIdxOfCells.size() == numCells);
// 	assert(numSpotRows > 0);
// 	assert(numSpotCols > 0);
// 	const i32 numSpots = numSpotRows * numSpotCols;
// 	assert(minCellsPerSpot > 0);
// 	assert(maxCellsPerSpot >= minCellsPerSpot);
// 	const set_i32 setCellTypes(typeIdxOfCells.begin(), typeIdxOfCells.end());
// 	assert(setCellTypes.find(typeAIdx) != setCellTypes.end());
// 	assert(setCellTypes.find(typeBIdx) != setCellTypes.end());
// 	const set_Str setGeneNames(geneNames.begin(), geneNames.end());
// 	assert(setGeneNames.size() == geneNames.size());

// 	vec_Str sampleNames;
// 	vec_Str spotNames;
// 	for (const i32 i : irange(0, numBulkSamples)) sampleNames.push_back("BS" + std::to_string(i));
// 	for (const i32 i : irange(0, numSpots)) spotNames.push_back("SS" + std::to_string(i));

// 	// Construct output data matrices
// 	DataMatrix Bulk(geneNames, sampleNames);
// 	DataMatrix BulkPhi(sampleNames, {"PhenotypeQuantity"});
// 	DataMatrix BulkCellCounts(sampleNames, cellTypeNames);
// 	DataMatrix SRT(geneNames, spotNames);
// 	DataMatrix SRTPhi(spotNames, {"PhenotypeQuantity"});
// 	DataMatrix SRTCategory(spotNames, {"SpotCategory"});

// 	vec_i32 typeACells;
// 	vec_i32 typeBCells;
// 	vec_i32 otherCells;
// 	const i32 cIdxRow = 0;
// 	const i32 cIdxCol = 0;
// 	const f64 distCore = numSpotRows/2.0;
// 	const f64 distInnerBelt = distCore + numSpotRows/4.0;
// 	const f64 distOuterBelt = distInnerBelt + numSpotRows/4.0;
// 	// vec_i32 spotType(numSpots);
// 	vec_i32 numCellsASRT(numSpots);
// 	vec_i32 numCellsBSRT(numSpots);
// 	vec_i32 numCellsOtherSRT(numSpots);
// 	vec_i32 numCellsABulk(numBulkSamples);
// 	vec_i32 numCellsBBulk(numBulkSamples);
// 	vec_i32 numCellsOtherBulk(numBulkSamples);
// 	std::mt19937 RNG(RNGSeed);
// 	std::uniform_int_distribution<i32> distNumCells(minCellsPerSpot, maxCellsPerSpot);


// 	// Store cell indices per category
// 	i32 otherTypeIdx = -1;
// 	for (const i32 cIdx : irange(0, numCells)) {
// 		if (typeIdxOfCells[cIdx] == typeAIdx) typeACells.push_back(cIdx);
// 		else if (typeIdxOfCells[cIdx] == typeBIdx) typeBCells.push_back(cIdx);
// 		else {
// 			otherCells.push_back(cIdx);
// 			otherTypeIdx = typeIdxOfCells[cIdx];
// 		}
// 	}
// 	vec_f64 avgA(numGenes, 0);
// 	vec_f64 avgB(numGenes, 0);
// 	vec_f64 avgC(numGenes, 0);
// 	for (const i32 gIdx : irange(0, numGenes)) {
// 		for (const i32 cellIdx : typeACells) avgA[gIdx] += scExp[gIdx, cellIdx];
// 		for (const i32 cellIdx : typeBCells) avgB[gIdx] += scExp[gIdx, cellIdx];
// 		for (const i32 cellIdx : otherCells) avgC[gIdx] += scExp[gIdx, cellIdx];
// 		avgA[gIdx] /= (i32)typeACells.size();
// 		avgB[gIdx] /= (i32)typeBCells.size();
// 		avgC[gIdx] /= (i32)otherCells.size();
// 	}


// 	// Determine spatial cell type category distribution
// 	for (const i32 rIdx : irange(0, numSpotRows)) {
// 		for (const i32 cIdx : irange(0, numSpotCols)) {
// 			const i32 sIdx = rIdx * numSpotCols + cIdx;
// 			const i32 totalNumCells = distNumCells(RNG);
// 			if (distSquared(cIdxRow, cIdxCol, rIdx, cIdx) <= square(distCore)) {
// 				// inside the core
// 				SRTCategory[sIdx, 0] = 3;
// 				numCellsASRT[sIdx] = totalNumCells * 0.8;
// 				numCellsBSRT[sIdx] = 0;
// 			}
// 			else if (distSquared(cIdxRow, cIdxCol, rIdx, cIdx) <= square(distInnerBelt)) {
// 				// inside the inner belt
// 				SRTCategory[sIdx, 0] = 2;
// 				numCellsASRT[sIdx] = totalNumCells * 0.4;
// 				numCellsBSRT[sIdx] = totalNumCells * 0.4;
// 			}
// 			else if (distSquared(cIdxRow, cIdxCol, rIdx, cIdx) <= square(distOuterBelt)) {
// 				// inside the outer belt
// 				SRTCategory[sIdx, 0] = 1;
// 				numCellsASRT[sIdx] = 0;
// 				numCellsBSRT[sIdx] = totalNumCells * 0.8;
// 			}
// 			else {
// 				// neutral cells
// 				SRTCategory[sIdx, 0] = 0;
// 				numCellsASRT[sIdx] = 0;
// 				numCellsBSRT[sIdx] = 0;
// 			}
// 			numCellsOtherSRT[sIdx] = totalNumCells -numCellsASRT[sIdx] -numCellsBSRT[sIdx];
// 		}
// 	}


// 	const f64 cellCountMean = avgCellsOfEachType;
// 	const f64 exponentialDistLambda = 1.0L/cellCountMean;
// 	std::exponential_distribution<f64> distCellQuantity(exponentialDistLambda);
// 	// Determine bulk cell type category distribution
// 	for (const i32 sIdx : irange(0, numBulkSamples)) {
// 		numCellsABulk[sIdx] = distCellQuantity(RNG);
// 		numCellsBBulk[sIdx] = distCellQuantity(RNG);
// 		numCellsOtherBulk[sIdx] = distCellQuantity(RNG);
// 		// fprintf(stderr, "%d %d %d\n", numCellsABulk[sIdx], numCellsBBulk[sIdx], numCellsOtherBulk[sIdx]);
// 		BulkCellCounts[sIdx, typeAIdx] = numCellsABulk[sIdx];
// 		BulkCellCounts[sIdx, typeBIdx] = numCellsBBulk[sIdx];
// 		BulkCellCounts[sIdx, otherTypeIdx] = numCellsOtherBulk[sIdx];
// 	}


// 	fprintf(stderr, "Sampling gene expression for spots...\n");
// 	// Sample expression for spots
// 	for (const i32 sIdx : irange(0, numSpots)) {
// 		assert(numCellsASRT[sIdx] >= 0);
// 		assert(numCellsBSRT[sIdx] >= 0);
// 		assert(numCellsOtherSRT[sIdx] >= 0);

// 		for (const i32 gIdx : irange(0, numGenes)) {
// 			SRT[gIdx, sIdx] += avgA[gIdx] * numCellsASRT[sIdx];
// 			SRT[gIdx, sIdx] += avgB[gIdx] * numCellsBSRT[sIdx];
// 			SRT[gIdx, sIdx] += avgC[gIdx] * numCellsOtherSRT[sIdx];
// 		}
// 		for (const i32 gIdx : geneIndicesA) {
// 			SRTPhi[sIdx, 0] += avgA[gIdx] * numCellsASRT[sIdx];
// 		}

// 		// if (SRTPhi[sIdx, 0] > 0) SRTPhi[sIdx, 0] = log10(SRTPhi[sIdx, 0]);
// 		// else if (SRTPhi[sIdx, 0] < 0) SRTPhi[sIdx, 0] = -log10(fabs(SRTPhi[sIdx, 0]));
// 	}
	
	
// 	// Sample expression for bulk samples
// 	for (const i32 sIdx : irange(0, numBulkSamples)) {
// 		assert(numCellsASRT[sIdx] >= 0);
// 		assert(numCellsBSRT[sIdx] >= 0);
// 		assert(numCellsOtherSRT[sIdx] >= 0);

// 		for (const i32 gIdx : irange(0, numGenes)) {
// 			Bulk[gIdx, sIdx] += avgA[gIdx] * numCellsABulk[sIdx];
// 			Bulk[gIdx, sIdx] += avgB[gIdx] * numCellsBBulk[sIdx];
// 			Bulk[gIdx, sIdx] += avgC[gIdx] * numCellsOtherBulk[sIdx];
// 		}
// 		for (const i32 gIdx : geneIndicesA) {
// 			BulkPhi[sIdx, 0] += avgA[gIdx] * numCellsABulk[sIdx];
// 		}

// 		// if (SRTPhi[sIdx, 0] > 0) SRTPhi[sIdx, 0] = log10(SRTPhi[sIdx, 0]);
// 		// else if (SRTPhi[sIdx, 0] < 0) SRTPhi[sIdx, 0] = -log10(fabs(SRTPhi[sIdx, 0]));
// 	}


	
	
	
// 	// for (const i32 sIdx : irange(0, SRTPhi.numRows)) {
// 	// 	if (SRTPhi[sIdx, 0] > 0) SRTPhi[sIdx, 0] = 1;
// 	// 	else if (SRTPhi[sIdx, 0] < 0) SRTPhi[sIdx, 0] = -1;
// 	// }
// 	// for (const i32 sIdx : irange(0, BulkPhi.numRows)) {
// 	// 	if (BulkPhi[sIdx, 0] > 0) BulkPhi[sIdx, 0] = 1;
// 	// 	else if (BulkPhi[sIdx, 0] < 0) BulkPhi[sIdx, 0] = -1;
// 	// }

// 	return {Bulk, BulkPhi, BulkCellCounts, SRT, SRTPhi, SRTCategory};
// }



#endif