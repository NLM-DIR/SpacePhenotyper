#ifndef __H__EIGENSTUFF__PREPROCESSING__
#define __H__EIGENSTUFF__PREPROCESSING__



#include "DataMatrix.h"
#include <algorithm>
#include <cmath>
#include <tuple>
#include <random>
#include <thread>
#include <mutex>



std::tuple<vec_Str, vec_i32, map_Str_i32> getColumnCellTypes(
	const DataMatrix& M,
	const map_Str_Str& barcodeCellTypes
) {
	vec_Str cellTypeNames;
	vec_i32 typeIdxOfCells;
	map_Str_i32 nameToTypeIdx;
	
	for (const Str& barcode : M.colNames) {
		const auto it = barcodeCellTypes.find(barcode);
		assert(it != barcodeCellTypes.end());
		const Str& typeName = it->second;
		if (!nameToTypeIdx.count(typeName)) {
			const i32 idx = (i32)nameToTypeIdx.size();
			nameToTypeIdx[typeName] = idx;
			cellTypeNames.push_back(typeName);
		}
		typeIdxOfCells.push_back( nameToTypeIdx[typeName] );
	}
	
	fprintf(stderr, "%d cell types:\n", (i32)nameToTypeIdx.size());
	for (const auto& [name, idx] : nameToTypeIdx) {
		fprintf(stderr, "%d\t%s\n", idx, name.c_str());
	}
	
	return {cellTypeNames, typeIdxOfCells, nameToTypeIdx};
}



void keepOnlyFollowingCellTypes(
	DataMatrix& M,
	const map_Str_Str& barcodeCellTypes,
	const vec_Str& cellTypesToKeep
) {
	auto [cellTypeNames, typeIdxOfCells, nameToTypeIdx] = getColumnCellTypes(M, barcodeCellTypes);
	set_i32 setTypeIndicesToKeep;
	for (const Str& name : cellTypesToKeep) {
		assert(nameToTypeIdx.count(name));
		setTypeIndicesToKeep.insert(nameToTypeIdx[name]);
	}
	vec_i32 barcodeIndicesToKeep;
	for (const i32 cIdx : irange(0, M.numCols)) {
		if (setTypeIndicesToKeep.count(typeIdxOfCells[cIdx])) {
			barcodeIndicesToKeep.push_back(cIdx);
		}
	}
	M.keepOnlyColumnsWithIndices(barcodeIndicesToKeep);
	{
		vec_i32 geneIndicesToKeep;
		for (const i32 gIdx : irange(0, M.numRows)) {
			bool keepIt = false;
			for (const i32 cIdx : irange(0, M.numCols)) {
				if (M[gIdx, cIdx] > 0) {
					keepIt = true;
					break;
				}
			}
			if (keepIt) {
				geneIndicesToKeep.push_back(gIdx);
			}
		}
		M.keepOnlyRowsWithIndices(geneIndicesToKeep);
	}
	const i32 numTypes = (i32)typeIdxOfCells.size();
	vec_i32 typeCellCount(numTypes, 0);
	for (const i32 cIdx : irange(0, M.numCols)) {
		// ++typeCellCount.at( typeIdxOfCells.at(cIdx) );
		++typeCellCount[ typeIdxOfCells[cIdx] ];
	}
	const set_Str setCellTypesToKeep(cellTypesToKeep.begin(), cellTypesToKeep.end());
	for (const Str& cellType : setCellTypesToKeep) {
		// const i32 typeIdx = nameToTypeIdx.at(cellType);
		const i32 typeIdx = nameToTypeIdx[cellType];
		// fprintf(stderr, "%s\t%d\n", cellType.c_str(), typeCellCount.at(typeIdx));
		fprintf(stderr, "%s\t%d\n", cellType.c_str(), typeCellCount[typeIdx]);
	}
}



// std::map<Str, vec_Str> findExclusiveGenes(
// 	const DataMatrix& M,
// 	const vec_Str& colCellTypes
// ) {
// 	struct Entry {
// 		Str gene;
// 		f64 medianCount;
// 		f64 avgCount;

// 		bool operator<(const Entry& other) const {
// 			if (medianCount != other.medianCount) return (medianCount > other.medianCount);
// 			else if (avgCount != other.avgCount) return (avgCount > other.avgCount);
// 			else return gene < other.gene;
// 		}
// 	};

// 	std::map<Str, i32> cellCount;
// 	for (const i32 cIdx : irange(0, M.numCols)) {
// 		cellCount[ colCellTypes[cIdx] ]++;
// 	}
// 	fprintf(stderr, "%d different cell types.\n", (i32)cellCount.size());

// 	std::map< Str, std::set<Entry> > exclusiveGenes;
// 	vec_f64 row(M.numCols, 0);
// 	for (const i32 gIdx : irange(0, M.numRows)) {
// 		set_Str typesExpressed;
// 		for (const i32 cIdx : irange(0, M.numCols)) {
// 			if (M[gIdx, cIdx] > 0) {
// 				typesExpressed.insert(colCellTypes[cIdx]);
// 			}
// 		}
// 		if (typesExpressed.size() == 1) {
// 			const Str cellType = *typesExpressed.begin();

// 			f64 avgCount = 0;
// 			row.clear();
// 			for (const i32 cIdx : irange(0, M.numCols)) {
// 				if (colCellTypes[cIdx] == cellType) {
// 					row.push_back(M[gIdx, cIdx]);
// 					avgCount += M[gIdx, cIdx];
// 				}
// 			}

// 			avgCount /= cellCount[cellType];
// 			std::sort(row.begin(), row.end());
// 			f64 medianCount = row[(i32)cellCount[cellType]/2];
// 			exclusiveGenes[cellType].insert({M.rowNames[gIdx], medianCount, avgCount});
// 		}
// 	}

// 	fprintf(stderr, "%d\n", (i32)exclusiveGenes.size());

// 	std::map<Str, vec_Str> selectedGenes;
// 	for (const auto& [cellType, entryList] : exclusiveGenes) {
// 		fprintf(stderr, "Cell type '%s'\n", cellType.c_str());
// 		for (const auto& [gene, medianCount, avgCount] : entryList) {
// 			if (medianCount < 1) break;
// 			selectedGenes[cellType].push_back(gene);
// 			fprintf(stderr, "%s\t%lf\t%lf\n", gene.c_str(), medianCount, avgCount);
// 		}
// 		fprintf(stderr, "\n");
// 	}

// 	return selectedGenes;
// }



void filterOutUnderexpressedGenes(
	DataMatrix& M,
	const vec_i32& typeIdxOfCells
) {
	const i32 numTypes = *std::max_element(typeIdxOfCells.begin(), typeIdxOfCells.end()) + 1;
	vec2D_i32 typeCellIndices(numTypes);
	for (const i32 cIdx : irange(0, M.numCols)) {
		// assert(typeIdxOfCells.at(cIdx) < numTypes);
		assert(typeIdxOfCells[cIdx] < numTypes);
		// typeCellIndices.at( typeIdxOfCells[cIdx] ).push_back(cIdx);
		typeCellIndices[ typeIdxOfCells[cIdx] ].push_back(cIdx);
	}
	
	vec_f64 row(M.numCols, 0);
	vec_i32 keepGeneIndices;
	for (const i32 gIdx : irange(0, M.numRows)) {
		vec_i32 typeMedians(numTypes);
		for (const i32 tIdx : irange(0, numTypes)) {
			i32 n = 0;
			for (const i32 cIdx : typeCellIndices[tIdx]) {
				row[n++] = M[gIdx, cIdx];
			}
			std::sort(row.begin(), row.begin() + n);
			typeMedians[tIdx] = row[n/2];
		}
		bool keepIt = false;
		for (const i32 val : typeMedians) {
			keepIt |= (val > 0);
		}
		if (keepIt) keepGeneIndices.push_back(gIdx);
	}
	fprintf(stderr, "%d out of %d genes have median count > 0 in at least one cell type. The rest are removed.\n", (i32)keepGeneIndices.size(), M.numRows);
	M.keepOnlyRowsWithIndices(keepGeneIndices);
}



void filterOutMitochondrialGenes(
	DataMatrix& M
) {
	vec_i32 keepGeneIndices;
	for (const i32 gIdx : irange(0, M.numRows)) {
		if (!((M.rowNames[gIdx].size() > 3) and (M.rowNames[gIdx][0] == 'M') and (M.rowNames[gIdx][1] == 'T') and (M.rowNames[gIdx][2] == '-'))) keepGeneIndices.push_back(gIdx);
	}
	fprintf(stderr, "%d out of %d genes are mitochondrial, and are removed.\n", M.numRows - (i32)keepGeneIndices.size(), M.numRows);
	M.keepOnlyRowsWithIndices(keepGeneIndices);
}



vec_f64 getVarCoef(const DataMatrix& M) {
	vec_f64 std(M.numRows, 0);
	vec_f64 avg(M.numRows, 0);
	for (const i32 gIdx : irange(0, M.numRows)) {
		for (const i32 cIdx : irange(0, M.numCols)) {
			avg[gIdx] += M[gIdx, cIdx];
		}
		avg[gIdx] /= M.numCols;
		for (const i32 cIdx : irange(0, M.numCols)) {
			const f64 x = avg[gIdx] -  M[gIdx, cIdx];
			std[gIdx] += x*x;
		}
		std[gIdx] /= M.numCols;
	}
	for (const i32 gIdx : irange(0, M.numRows)) {
		assert(avg[gIdx] > 0);
		std[gIdx] /= avg[gIdx];
	}
	return std;
}



vec2D_i32 selectMarkersExclusiveGenes(
	const DataMatrix& SC,
	const vec_Str& cellTypeNames,
	const vec_i32& typeIdxOfCells
) {
	struct Entry {
		i32 gIdx;
		f64 varCoef;
		
		bool operator<(const Entry& other) const {
			if (varCoef != other.varCoef) return varCoef > other.varCoef;
			else return gIdx < other.gIdx;
		}
	};
	const vec_f64 varCoefs = getVarCoef(SC);
	const i32 numTypes = (i32)cellTypeNames.size();
	std::vector< std::set<Entry> > exclusiveGenes(numTypes);
	for (const i32 gIdx : irange(0, SC.numRows)) {
		set_i32 typesExpressed;
		for (const i32 cIdx : irange(0, SC.numCols)) {
			if (SC[gIdx, cIdx] > 0) typesExpressed.insert(typeIdxOfCells[cIdx]);
		}
		if (typesExpressed.size() == 1) {
			const i32 typeIdx = *typesExpressed.begin();
			// exclusiveGenes.at(typeIdx).insert({gIdx, varCoefs[gIdx]});
			exclusiveGenes[typeIdx].insert({gIdx, varCoefs[gIdx]});
		}
	}
	
	for (const i32 typeIdx : irange(0, numTypes)) {
		fprintf(stderr, "%s\n", cellTypeNames[typeIdx].c_str());
		// for (const auto& [gIdx, varCoef] : exclusiveGenes.at(typeIdx)) {
		for (const auto& [gIdx, varCoef] : exclusiveGenes[typeIdx]) {
			fprintf(stderr, "%s\t%le\n", SC.rowNames[gIdx].c_str(), varCoef);
		}
		fprintf(stderr, "\n");
	}
	vec2D_i32 markers(numTypes);
	for (const i32 typeIdx : irange(0, numTypes)) {
		// for (const auto& [gIdx, varCoef] : exclusiveGenes.at(typeIdx)) {
		for (const auto& [gIdx, varCoef] : exclusiveGenes[typeIdx]) {
			// markers.at(typeIdx).push_back(gIdx);
			markers[typeIdx].push_back(gIdx);
		}
	}
	return markers;
}



vec2D_i32 selectMarkersSharedGenes(
	const DataMatrix& SC,
	const vec_Str& cellTypeNames,
	const vec_i32& typeIdxOfCells,
	const i32 maxNumMarkersPerCellType
) {
	struct Entry {
		i32 gIdx;
		f64 varCoef;
		
		bool operator<(const Entry& other) const {
			if (varCoef != other.varCoef) return varCoef > other.varCoef;
			else return gIdx < other.gIdx;
		}
	};
	const vec_f64 varCoefs = getVarCoef(SC);
	const i32 numTypes = (i32)cellTypeNames.size();
	std::set<Entry> sharedGenes;
	for (const i32 gIdx : irange(0, SC.numRows)) {
		set_i32 typesExpressed;
		for (const i32 cIdx : irange(0, SC.numCols)) {
			if (SC[gIdx, cIdx] > 0) typesExpressed.insert(typeIdxOfCells[cIdx]);
		}
		if (typesExpressed.size() == 3) {
			sharedGenes.insert({gIdx, varCoefs[gIdx]});
		}
	}
	vec2D_i32 markers(numTypes);
	const i32 numMarkers = std::min(maxNumMarkersPerCellType, ((i32)sharedGenes.size())/3);
	std::set<Entry>::iterator it = sharedGenes.begin();
	for (const i32 typeIdx : irange(0, numTypes)) {
		for (i32 _ = numMarkers; _--; ) {
			markers[typeIdx].push_back(it->gIdx);
			++it;
		}
	}
	for (const i32 typeIdx : irange(0, numTypes)) {
		fprintf(stderr, "%s\n", cellTypeNames[typeIdx].c_str());
		// for (const auto& [gIdx, varCoef] : exclusiveGenes.at(typeIdx)) {
		for (const i32 gIdx : markers[typeIdx]) {
			fprintf(stderr, "%s\t%le\n", SC.rowNames[gIdx].c_str(), varCoefs[gIdx]);
		}
		fprintf(stderr, "\n");
	}
	return markers;
}



f64 pearsonR(const vec_f64& A, const vec_f64& B) {
	f64 avgA = 0;
	f64 avgB = 0;
	f64 normA = 0;
	f64 normB = 0;
	const i32 n = (i32)A.size();
	for (const i32 i : irange(0, n)) {
		avgA += A[i];
		avgB += B[i];
	}
	avgA /= n;
	avgB /= n;
	f64 r = 0;
	for (const i32 i : irange(0, n)) {
		normA += (A[i] - avgA)*(A[i] - avgA);
		normB += (B[i] - avgB)*(B[i] - avgB);
		r += (A[i] - avgA)*(B[i] - avgB);
	}
	return r / sqrt(normA * normB);
}



void singleThreadPermutationTest(
	const i32 startingGeneIdx,
	const i32 numThreads,
	std::mutex& mtxCount,
	i32& printCounter,
	std::mt19937& RNG,
	const DataMatrix& Bulk,
	const vec_f64& phi,
	DataMatrix& R,
	const i32 numPermutations
) {
	const i32 numGenes = Bulk.numRows;
	const i32 numSamples = Bulk.numCols;
	vec_f64 Bg(numSamples);
	for (i32 gIdx = startingGeneIdx; gIdx < numGenes; gIdx += numThreads) {
		{
			std::lock_guard<std::mutex> lock(mtxCount);
			fprintf(stderr, "\r%d / %d", printCounter+1, numGenes);
			++printCounter;
		}
		for (const i32 sIdx : irange(0, numSamples)) Bg[sIdx] = Bulk[gIdx, sIdx];
		R[gIdx, 0] = pearsonR(phi, Bg);
		i32 numAtLeastAbs = 0;
		for ([[maybe_unused]] const i32 pIdx : irange(0, numPermutations)) {
			std::shuffle(Bg.begin(), Bg.end(), RNG);
			const f64 r = pearsonR(phi, Bg);
			if ((r * R[gIdx, 0] > 0) and (fabs(r) >= fabs(R[gIdx, 0]))) ++numAtLeastAbs;
		}
		R[gIdx, 1] = f64(numAtLeastAbs) / f64(numPermutations);
	}
}



DataMatrix geneCorrelationWithPhenotype(
	const DataMatrix& Bulk,
	const DataMatrix& BulkPhi,
	const u32 RNGSeed,
	const i32 numThreads,
	const i32 numPermutations
) {
	DataMatrix R(Bulk.rowNames, {"pearsonR", "permPVal"});
	vec_f64 phi(Bulk.numCols);
	vec_f64 Bg(Bulk.numCols);
	for (const i32 sIdx : irange(0, Bulk.numCols)) phi[sIdx] = BulkPhi[sIdx, 0];
	
	fprintf(stderr, "Calculating p-values for gene-expression-to-phenotype pearson correlation...\n");
	
	std::vector< std::thread > threads(numThreads);
	i32 mutexedCounter = 0;
	std::mutex mtxCount;
	std::mt19937 TheRNG(RNGSeed);
	std::vector<std::mt19937> RNGs;
	for ([[maybe_unused]] const i32 tIdx : irange(0, numThreads)) {
		RNGs.emplace_back(TheRNG());
	}
	for (const i32 tIdx : irange(0, numThreads)) {
		threads[tIdx] = std::thread(
			singleThreadPermutationTest,
			/*startingGeneIdx*/ tIdx,
			numThreads,
			std::ref(mtxCount),
			std::ref(mutexedCounter),
			std::ref(RNGs[tIdx]),
			std::ref(Bulk),
			std::ref(phi),
			std::ref(R),
			numPermutations
		);
	}
	for (const i32 tIdx : irange(0, numThreads)) {
		threads[tIdx].join();
	}
	fprintf(stderr, "\n");
	
	return R;
}



#endif