#include "DataMatrix.h"
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <cstring>
#include <cmath>



void DataMatrix::allocateDataMemory() {
	data = new f64 [numRows*numCols];
	dataRows = new f64* [numRows];
	for (const i32 i : irange(0, numRows)) {
		dataRows[i] = data + i*numCols;
	}
	std::fill(data, data + numRows*numCols, f64(0));
}



void DataMatrix::deallocateDataMemory() {
	delete[] dataRows;
	delete[] data;
}



DataMatrix::DataMatrix(
	const DataMatrix& other,
	const DataMatrix_CopyFlag copyVals
) {
	numRows = other.numRows;
	numCols = other.numCols;
	rowNames = other.rowNames;
	colNames = other.colNames;
	rowNameToIdx = other.rowNameToIdx;
	colNameToIdx = other.colNameToIdx;
	allocateDataMemory();
	if (copyVals == FLAG_DataMatrix_COPYVALS_YES) {
		std::memcpy(data, other.data, sizeof(data[0])*numRows*numCols);
	}
}



DataMatrix::DataMatrix(
	const vec_Str& _rowNames,
	const vec_Str& _colNames
) {
	rowNames = _rowNames;
	colNames = _colNames;
	numRows = (i32)rowNames.size();
	numCols = (i32)colNames.size();
	for (const i32 i : irange(0, numRows)) {
		const Str& rowName = rowNames[i];
		rowNameToIdx[rowName] = i;
	}
	for (const i32 j : irange(0, numCols)) {
		const Str& colName = colNames[j];
		colNameToIdx[colName] = j;
	}
	allocateDataMemory();
}



DataMatrix::DataMatrix(
	const Str& path,
	const char sep,
	const DataMatrix_HeaderFlag& headerFlag
) : 
	data(nullptr),
	dataRows(nullptr)
{
	if (headerFlag == FLAG_DataMatrix_HEADER_YES) readFromPath(path, sep);
	else readFromPathNoHeader(path, sep);
}



f64& DataMatrix::operator() (const i32 rowIdx, const i32 colIdx) {
	if ((rowIdx < 0) or (rowIdx >= numRows)) {
		fprintf(stderr, "\n\nERROR in DataMatrix::operator() : Row index %d is out of bounds [0, %d]. Exiting.\n", rowIdx, numRows);
		exit(1);
	}
	else if ((colIdx < 0) or (colIdx >= numCols)) {
		fprintf(stderr, "\n\nERROR in DataMatrix::operator() : Column index %d is out of bounds [0, %d]. Exiting.\n", colIdx, numCols);
		exit(1);
	}
	// else return dataRows[rowIdx][colIdx];
	else return operator[](rowIdx, colIdx);
}



const f64& DataMatrix::operator() (const i32 rowIdx, const i32 colIdx) const {
	if ((rowIdx < 0) or (rowIdx >= numRows)) {
		fprintf(stderr, "\n\nERROR in DataMatrix::operator() : Row index %d is out of bounds [0, %d]. Exiting.\n", rowIdx, numRows);
		exit(1);
	}
	else if ((colIdx < 0) or (colIdx >= numCols)) {
		fprintf(stderr, "\n\nERROR in DataMatrix::operator() : Column index %d is out of bounds [0, %d]. Exiting.\n", colIdx, numCols);
		exit(1);
	}
	// else return dataRows[rowIdx][colIdx];
	else return operator[](rowIdx, colIdx);
}



i32 DataMatrix::getRowIdx(const Str& rowName) const {
	auto it = rowNameToIdx.find(rowName);
	if (it == rowNameToIdx.end()) {
		fprintf(stderr, "\n\nERROR in DataMatrix::getRowIdx : No row found with the name '%s'. Exiting.\n", rowName.c_str());
		exit(1);
	}
	return it->second;
}



i32 DataMatrix::getColIdx(const Str& colName) const {
	auto it = colNameToIdx.find(colName);
	if (it == colNameToIdx.end()) {
		fprintf(stderr, "\n\nERROR in DataMatrix::getColIdx : No column found with the name '%s'. Exiting.\n", colName.c_str());
		exit(1);
	}
	return it->second;
}



f64 DataMatrix::getVal(const Str& rowName, const Str& colName) const {
	const i32 rowIdx = getRowIdx(rowName);
	const i32 colIdx = getColIdx(colName);
	// return dataRows[rowIdx][colIdx];
	return operator[](rowIdx, colIdx);
}



void DataMatrix::printSize(const Str& matrixName) const {
	fprintf(stderr, "%s (%d x %d)\n", matrixName.c_str(), numRows, numCols);
}



void DataMatrix::validateNonNegativeValues(const Str& nameForPrint) const {
	for (const i32 i : irange(0, numRows)) {
		for (const i32 j : irange(0, numCols)) {
			// if (dataRows[i][j] < 0) {
			if (operator[](i, j) < 0) {
				fprintf(stderr, "\n\nERROR in DataMatrix::checkNegativeValues : '%s'[%d][%d] = %.14lf < 0. Exiting.\n", nameForPrint.c_str(), i, j, operator[](i, j));
				exit(1);
			}
		}
	}
}



void DataMatrix::setVal(const Str& rowName, const Str& colName, const f64& x) {
	const i32 rowIdx = getRowIdx(rowName);
	const i32 colIdx = getColIdx(colName);
	// dataRows[rowIdx][colIdx] = x;
	operator[](rowIdx, colIdx) = x;
}



void DataMatrix::resetToZero() {
	std::fill(data, data + numRows*numCols, f64(0));
}



void DataMatrix::transpose() {
	f64* newData = new f64 [numRows*numCols];
	f64** newDataRows = new f64* [numCols];
	for (const i32 j : irange(0, numCols)) {
		newDataRows[j] = newData + j*numRows;
	}

	for (const i32 i : irange(0, numRows)) {
		for (const i32 j : irange(0, numCols)) {
			// newDataRows[j][i] = dataRows[i][j];
			newData[j*numRows + i] = data[i*numCols + j];
			// std::swap(data[j*numRows + i], data[i*numCols + j]);
		}
	}
	delete[] data;
	delete[] dataRows;

	data = newData;
	dataRows = newDataRows;
	std::swap(numCols, numRows);
	std::swap(colNames, rowNames);
	std::swap(colNameToIdx, rowNameToIdx);
}



void DataMatrix::lexicographicallyOrderRows() {
	DataMatrix origCopy(*this, FLAG_DataMatrix_COPYVALS_YES);
	std::sort(rowNames.begin(), rowNames.end());
	vec_i32 myDataIsAt(numRows);
	vec_i32 dataBelongsTo(numRows);
	for (const i32 i : irange(0, numRows)) {
		const Str& rowName = rowNames[i];
		const i32 oldRowIdx = rowNameToIdx[rowName];
		myDataIsAt[i] = oldRowIdx;
		dataBelongsTo[oldRowIdx] = i;
	}
	for (const i32 i : irange(0, numRows)) {
		if (i != dataBelongsTo[i]) {
			for (const i32 j : irange(0, numCols)) {
				std::swap(data[i*numCols + j], data[myDataIsAt[i]*numCols + j]);
			}
			const i32 myDataWasAt = myDataIsAt[i];
			const i32 dataBelongedTo = dataBelongsTo[i];
			myDataIsAt[dataBelongedTo] = myDataWasAt;
			dataBelongsTo[myDataWasAt] = dataBelongedTo;
			myDataIsAt[i] = i;
			dataBelongsTo[i] = i;
		}
	}
	for (const i32 i : irange(0, numRows)) {
		const Str& rowName = rowNames[i];
		rowNameToIdx[rowName] = i;
	}
	// Validate correctness
	for (const i32 i : irange(0, numRows)) {
		const Str& rowName = this->rowNames[i];
		for (const i32 j : irange(0, numCols)) {
			const Str& colName = this->colNames[j];
			assert(this->getVal(rowName, colName) == origCopy.getVal(rowName, colName));
		}
	}
}



void DataMatrix::lexicographicallyOrderColumns() {
	DataMatrix origCopy(*this, FLAG_DataMatrix_COPYVALS_YES);
	f64* newData = new f64 [numRows*numCols];
	f64** newDataRows = new f64* [numRows];
	for (const i32 i : irange(0, numRows)) {
		newDataRows[i] = newData + i*numCols;
	}
	std::sort(colNames.begin(), colNames.end());
	for (const i32 j : irange(0, numCols)) {
		const Str& colName = colNames[j];
		const i32 oldColIdx = colNameToIdx[colName];
		for (const i32 i : irange(0, numRows)) {
			// newDataRows[i][j] = dataRows[i][oldColIdx];
			newData[i*numCols + j] = data[i*numCols + oldColIdx];
		}
	}
	for (const i32 j : irange(0, numCols)) {
		const Str& colName = colNames[j];
		colNameToIdx[colName] = j;
	}
	delete[] data;
	delete[] dataRows;
	data = newData;
	dataRows = newDataRows;
	// Validate correctness
	for (const i32 i : irange(0, numRows)) {
		const Str& rowName = this->rowNames[i];
		for (const i32 j : irange(0, numCols)) {
			const Str& colName = this->colNames[j];
			assert(this->getVal(rowName, colName) == origCopy.getVal(rowName, colName));
		}
	}
}



void DataMatrix::lexicographicallyOrderRowsAndColumns() {
	this->lexicographicallyOrderRows();
	this->lexicographicallyOrderColumns();
}



void DataMatrix::keepOnlyRowsWithIndices(const vec_i32& indices) {
	set_i32 keepOnly(indices.begin(), indices.end());
	i32 newSize = 0;
	for (const i32 i : irange(0, numRows)) {
		if (keepOnly.count(i)) {
			if (i != newSize) {
				rowNames[newSize] = rowNames[i];
				memcpy(data + newSize*numCols, data + i*numCols, sizeof(data[0])*numCols);
			}
			newSize++;
		}
	}
	numRows = newSize;
	rowNames.resize(newSize);
	rowNameToIdx.clear();
	for (const i32 i : irange(0, numRows)) {
		rowNameToIdx[ rowNames[i] ] = i;
	}
	// Internal check
	for (const i32 i : irange(0, numRows)) {
		const Str& rowName = rowNames[i];
		assert(this->getRowIdx(rowName) == i);
		for (const i32 j : irange(0, numCols)) {
			const Str& colName = colNames[j];
			assert(this->getColIdx(colName) == j);
			// assert(this->getVal(rowName, colName) == dataRows[i][j]);
			assert(this->getVal(rowName, colName) == operator[](i, j));
		}
	}
}



void DataMatrix::keepOnlyFollowingRows(const vec_Str& keepRows) {
	vec_bool keepIt(numRows, false);
	for (const Str& row : keepRows) {
		if (rowNameToIdx.count(row)) {
			keepIt[ rowNameToIdx.find(row)->second ] = true;
		}
	}

	i32 newSize = 0;
	for (const i32 i : irange(0, numRows)) {
		if (keepIt[i]) {
			if (i != newSize) {
				rowNames[newSize] = rowNames[i];
				memcpy(data + newSize*numCols, data + i*numCols, sizeof(data[0])*numCols);
			}
			newSize++;
		}
	}
	numRows = newSize;
	rowNames.resize(newSize);
	rowNameToIdx.clear();
	for (const i32 i : irange(0, numRows)) {
		rowNameToIdx[ rowNames[i] ] = i;
	}
	// Internal check
	for (const i32 i : irange(0, numRows)) {
		const Str& rowName = rowNames[i];
		assert(this->getRowIdx(rowName) == i);
		for (const i32 j : irange(0, numCols)) {
			const Str& colName = colNames[j];
			assert(this->getColIdx(colName) == j);
			// assert(this->getVal(rowName, colName) == dataRows[i][j]);
			assert(this->getVal(rowName, colName) == operator[](i, j));
		}
	}
}



void DataMatrix::keepOnlyColumnsWithIndices(const vec_i32& indices) {
	vec_bool keepIt(numCols, false);
	// for (const i32 idx : indices) keepIt.at(idx) = true;
	for (const i32 idx : indices) keepIt[idx] = true;
	
	i32 newSize = 0;
	for (const i32 i : irange(0, numRows)) {
		newSize = 0;
		for (const i32 j : irange(0, numCols)) {
			if (keepIt[j]) {
				operator[](i, newSize) = operator[](i, j);
				++newSize;
			}
		}
	}
	// Move now-shrunk row blocks so that they are contiguous in memory
	for (const i32 i : irange(0, numRows)) {
		memcpy(data + i*newSize, data + i*numCols, sizeof(data[0])*newSize);
		dataRows[i] = data + i*newSize;
	}
	// Fix colNames and colNameToIdx
	{
		newSize = 0;
		for (const i32 j : irange(0, numCols)) {
			if (keepIt[j]) {
				colNames[newSize++] = colNames[j];
			}
		}
		numCols = newSize;
		colNames.resize(newSize);
		colNameToIdx.clear();
		for (const i32 j : irange(0, numCols)) {
			colNameToIdx[ colNames[j] ] = j;
		}
	}
	// Internal check
	for (const i32 i : irange(0, numRows)) {
		const Str& rowName = rowNames[i];
		assert(this->getRowIdx(rowName) == i);
		for (const i32 j : irange(0, numCols)) {
			const Str& colName = colNames[j];
			assert(this->getColIdx(colName) == j);
			// assert(this->getVal(rowName, colName) == dataRows[i][j]);
			assert(this->getVal(rowName, colName) == operator[](i, j));
		}
	}
}



void DataMatrix::keepOnlyFollowingColumns(const vec_Str& keepCols) {
	vec_bool keepIt(numCols, false);
	for (const Str& col : keepCols) {
		if (colNameToIdx.count(col)) {
			keepIt[ colNameToIdx.find(col)->second ] = true;
		}
	}

	i32 newSize = 0;
	for (const i32 j : irange(0, numCols)) {
		if (keepIt[j]) {
			if (j != newSize) {
				// copy the column and the column name
				for (const i32 i : irange(0, numRows)) {
					// dataRows[i][newSize] = dataRows[i][j];
					operator[](i, newSize) = operator[](i, j);
				}
				colNames[newSize] = colNames[j];
			}
			newSize++;
		}
	}
	for (const i32 i : irange(0, numRows)) {
		memcpy(data + i*newSize, data + i*numCols, sizeof(data[0])*newSize);
		dataRows[i] = data + i*newSize;
	}
	numCols = newSize;
	colNames.resize(newSize);
	colNameToIdx.clear();
	for (const i32 j : irange(0, numCols)) {
		colNameToIdx[ colNames[j] ] = j;
	}
	// Internal check
	for (const i32 i : irange(0, numRows)) {
		const Str& rowName = rowNames[i];
		assert(this->getRowIdx(rowName) == i);
		for (const i32 j : irange(0, numCols)) {
			const Str& colName = colNames[j];
			assert(this->getColIdx(colName) == j);
			// assert(this->getVal(rowName, colName) == dataRows[i][j]);
			assert(this->getVal(rowName, colName) == operator[](i, j));
		}
	}
}



void DataMatrix::makeColumnsFractions() {
	validateNonNegativeValues("");
	for (const i32 j : irange(0, numCols)) {
		f64 colSum = 0;
		for (const i32 i : irange(0, numRows)) {
			// colSum += dataRows[i][j];
			colSum += operator[](i, j);
		}
		for (const i32 i : irange(0, numRows)) {
			// dataRows[i][j] /= colSum;
			operator[](i, j) /= colSum;
		}
	}
}



void DataMatrix::normalizeColumns(const f64& scale) {
	makeColumnsFractions();
	vec_i32 keepRowIndices;
	for (const i32 i : irange(0, numRows)) {
		bool keepIt = false;
		for (const i32 j : irange(0, numCols)) {
			// dataRows[i][j] *= scale;
			operator[](i, j) *= scale;
			// if (dataRows[i][j] < 1) dataRows[i][j] = 0;
			if (operator[](i, j) < 1) operator[](i, j) = 0;
			else keepIt = true;
		}
		if (keepIt == true) keepRowIndices.push_back(i);
	}
	keepOnlyRowsWithIndices(keepRowIndices);
}



void DataMatrix::addMatrix(const DataMatrix& other) {
	if (this->numRows != other.numRows) {
		fprintf(stderr, "\n\nError in DataMatrix::addMatrix : The two matrices do not have the same number of rows. Exiting.\n");
		exit(1);
	}
	if (this->numCols != other.numCols) {
		fprintf(stderr, "\n\nError in DataMatrix::addMatrix : The two matrices do not have the same number of columns. Exiting.\n");
		exit(1);
	}
	if (this->rowNames != other.rowNames) {
		fprintf(stderr, "\n\nError in DataMatrix::addMatrix : The two matrices do not have the same row names (the exact order DOES matter). Exiting.\n");
		exit(1);
	}
	if (this->colNames != other.colNames) {
		fprintf(stderr, "\n\nError in DataMatrix::addMatrix : The two matrices do not have the same column names (the exact order DOES matter). Exiting.\n");
		exit(1);
	}
	for (const i32 i : irange(0, numRows)) {
		for (const i32 j : irange(0, numCols)) {
			// dataRows[i][j] += other.dataRows[i][j];
			operator[](i, j) += other[i, j];
		}
	}
}



void DataMatrix::readFromPath(const Str& path, const char sep, const bool noName) {
	// fprintf(stderr, "%d\n", __LINE__);
	deallocateDataMemory();
	// fprintf(stderr, "%d\n", __LINE__);
	numCols = 0;
	numRows = 0;
	i32 maxRecordLen = 0;
	char cc;
	FILE * fin = fopen(path.c_str(), "r");
	if (fin == nullptr) {
		fprintf(stderr, "\n\nERROR in DataMatrix::readFromPath : Cannot open [%s] for reading. Exiting.\n", path.c_str());
		exit(1);
	}
	i32 recordCharLen = 0;
	while (((cc = getc(fin)) != EOF) and (cc != '\n')) {
		if (cc == sep) {
			numCols++;
			if (recordCharLen > maxRecordLen) maxRecordLen = recordCharLen;
			recordCharLen = 0;
		}
		else {
			recordCharLen++;
		}
	}
	if (recordCharLen > maxRecordLen) maxRecordLen = recordCharLen;	// We need this line because the last row ended with '\n' instead of sep
	if (cc == EOF) {
		fprintf(stderr, "\n\nERROR in DataMatrix::readFromPath : The data matrix in [%s] has no data rows. Exiting application.\n", path.c_str());
		exit(1);
	}
	if (noName == true) numCols++;

	recordCharLen = 0;
	i32 currentCols = 0;
	bool doneWithTheName = false;
	while ((cc = getc(fin)) != EOF) {
		if (cc == sep) {
			currentCols++;
			if (!doneWithTheName) {
				if (recordCharLen > maxRecordLen) maxRecordLen = recordCharLen;
				recordCharLen = 0;
				doneWithTheName = true;
			}
		}
		else if (cc == '\n') {
			if ((currentCols > 0) or (recordCharLen > 0)) {
				if ((recordCharLen > 0) or (currentCols != numCols)) {
					fprintf(stderr, "\n\nERROR in DataMatrix::readFromPath : Row %d (0-indexed) of the file [%s] has %d columns, vs %d in the header. Exiting application.\n", numRows + 1, path.c_str(), currentCols, numCols);
					exit(1);
				}
				currentCols = 0;
				numRows++;
				doneWithTheName = false;
			}
		}
		else {
			if (!doneWithTheName) recordCharLen++;
		}
	}
	if (currentCols > 0) {	// In case the last row ended with EOF instead of '\n'
		if ((recordCharLen > 0) or (currentCols != numCols)) {
			fprintf(stderr, "\n\nERROR in DataMatrix::readFromPath : Row %d (0-indexed) of the file [%s] has %d columns, vs %d in the header. Exiting application.\n", numRows + 1, path.c_str(), currentCols, numCols);
			exit(1);
		}
		numRows++;
	}
	rowNames.resize(numRows);
	colNames.resize(numCols);
	allocateDataMemory();

	rewind(fin);
	char* const temp = new char[maxRecordLen + 1];
	// Skip the first word in the header (table name)
	if (noName == false) while (((cc = getc(fin))) != EOF and (cc != sep));
	// Read the column names from the rest of the header
	for (const i32 j : irange(0, numCols)) {
		i32 len = 0;
		while (((cc = getc(fin)) != EOF) and (cc != sep) and (cc != '\r') and (cc != '\n')) temp[len++] = cc;
		temp[len] = 0;
		if (colNameToIdx.count(temp)) {
			fprintf(stderr, "\n\nERROR in DataMatrix::readFromPath : Column name '%s' appears multiple times in the header row of [%s]. Exiting application.\n", temp, path.c_str());
			exit(1);
		}
		colNameToIdx[temp] = j;
		colNames[j] = temp;
	}

	i32 len;
	for (const i32 i : irange(0, numRows)) {
		// Read the row name
		len = 0;
		while (((cc = getc(fin)) != EOF) and (cc != sep) and (cc != '\r') and (cc != '\n')) temp[len++] = cc;
		temp[len] = 0;
		if (rowNameToIdx.count(temp)) {
			fprintf(stderr, "\n\nERROR in DataMatrix::readFromPath : Row name '%s' appears multiple times in the first column of [%s]. Exiting application.\n", temp, path.c_str());
			exit(1);
		}
		rowNameToIdx[temp] = i;
		rowNames[i] = temp;
		// Read the row data
		for (const i32 j : irange(0, numCols)) {
			// Read the current record into the buffer
			len = 0;
			while (((cc = getc(fin)) != EOF) and (cc != sep) and (cc != '\r') and (cc != '\n')) temp[len++] = cc;
			temp[len] = 0;
			// dataRows[i][j] = std::strtod(temp, nullptr);
			operator[](i, j) = std::strtod(temp, nullptr);
		}
		// Skip to the new line if there's trailing whitespace
		if (cc != EOF and cc != '\n') {
			while (((cc = getc(fin)) != EOF) and (cc != '\n'));
		}
	}

	delete[] temp;
	fclose(fin);
}



void DataMatrix::readFromPathNoHeader(const Str& path, const char sep) {
	deallocateDataMemory();
	i32 maxRecordLen = 0;
	i32 numSeps = 0;
	char cc;
	FILE * fin = fopen(path.c_str(), "r");
	if (fin == nullptr) {
		fprintf(stderr, "\n\nERROR in DataMatrix::readFromPathNoHeader : Cannot open [%s] for reading. Exiting.\n", path.c_str());
		exit(1);
	}
	i32 recordCharLen = 0;
	while (((cc = getc(fin)) != EOF) and (cc != '\n')) {
		if (cc == sep) {
			numSeps++;
			if (recordCharLen > maxRecordLen) maxRecordLen = recordCharLen;
			recordCharLen = 0;
		}
		else {
			recordCharLen++;
		}
	}
	numCols = numSeps + 1;

	numRows = 1;
	i32 currentSeps = 0;
	while ((cc = getc(fin)) != EOF) {
		if (cc == sep) {
			currentSeps++;
			if (recordCharLen > maxRecordLen) maxRecordLen = recordCharLen;
			recordCharLen = 0;
		}
		else if (cc == '\n') {
			if (currentSeps > 0) {
				if (currentSeps != numSeps) {
					fprintf(stderr, "\n\nERROR in DataMatrix::readFromPathNoHeader : Row %d (0-indexed) of the file [%s] has %d columns, vs %d in the first line. Exiting application.\n", numRows + 1, path.c_str(), currentSeps + 1, numSeps + 1);
					exit(1);
				}
				currentSeps = 0;
				numRows++;
			}
			if (recordCharLen > maxRecordLen) maxRecordLen = recordCharLen;
			recordCharLen = 0;
		}
		else {
			recordCharLen++;
		}
	}
	if (currentSeps > 0) {
		if (currentSeps != numSeps) {
			fprintf(stderr, "\n\nERROR in DataMatrix::readFromPathNoHeader : Row %d (0-indexed) of the file [%s] has %d columns, vs %d in the first line. Exiting application.\n", numRows + 1, path.c_str(), currentSeps + 1, numSeps + 1);
			exit(1);
		}
		numRows++;
		if (recordCharLen > maxRecordLen) maxRecordLen = recordCharLen;
	}
	rowNames.resize(numRows);
	colNames.resize(numCols);
	for (const i32 i : irange(0, numRows)) {
		rowNames[i] = "Row" + std::to_string(i);
		rowNameToIdx[ rowNames[i] ] = i;
	}
	for (const i32 j : irange(0, numCols)) {
		colNames[j] = "Column" + std::to_string(j);
		colNameToIdx[ colNames[j] ] = j;
	}
	allocateDataMemory();

	rewind(fin);
	char * temp = new char[maxRecordLen + 1];
	i32 len;
	for (const i32 i : irange(0, numRows)) {
		// Read the row data
		for (const i32 j : irange(0, numCols)) {
			len = 0;
			while (((cc = getc(fin)) != EOF) and (cc != sep) and (cc != '\r') and (cc != '\n')) temp[len++] = cc;
			temp[len] = 0;
			// dataRows[i][j] = std::strtod(temp, nullptr);
			operator[](i, j) = std::strtod(temp, nullptr);
		}
		// Skip to the new line if there's trailing whitespace
		if (cc != EOF and cc != '\n') {
			while (((cc = getc(fin)) != EOF) and (cc != '\n'));
		}
	}

	delete[] temp;
	fclose(fin);
}



void DataMatrix::printToFilestreamAsTSV(
	FILE * const stream,
	const Str& firstHeaderCol,
	const DataMatrix* const dmR /*= nullptr*/
) const {
	const f64 EPSILON = 1e-6;
	fprintf(stream, "%s", firstHeaderCol.c_str());
	for (const i32 j : irange(0, numCols)) {
		fprintf(stream, "\t%s", colNames[j].c_str());
	}
	fprintf(stream, "\n");
	for (const i32 i : irange(0, numRows)) {
		fprintf(stream, "%s", rowNames[i].c_str());
		for (const i32 j : irange(0, numCols)) {
			// if ((dmR != nullptr) and (dmR->dataRows[i][j] == 0)) {
			if ((dmR != nullptr) and (dmR->operator[](i, j) == 0)) {
				fprintf(stream, "\tNaN");
			}
			else {
				const bool isInteger = (fabs(round(operator[](i, j)) - operator[](i, j)) < EPSILON);

				if (isInteger) fprintf(stream, "\t%.0lf", operator[](i, j));
				else if (fabs(operator[](i, j)) >= 10) fprintf(stream, "\t%.1lf", operator[](i, j));
				else fprintf(stream, "\t%.6lf", operator[](i, j));
			}
		}
		fprintf(stream, "\n");
	}
}



void DataMatrix::printToFileAsTSV(
	const Str& filename,
	const Str& firstHeaderCol,
	const DataMatrix* const dmR /*= nullptr*/
) const {
	if (FILE * fout = fopen(filename.c_str(), "w"); fout) {
		printToFilestreamAsTSV(fout, firstHeaderCol, dmR);
		fclose(fout);
		fprintf(stderr, "Data matrix written to [%s].\n", filename.c_str());
	}
	else {
		fprintf(stderr, "\n\nERROR in DataMatrix::printToFile : Cannot open [%s] for writing. Exiting.\n", filename.c_str());
	}
}



void intersectMatrixRows(DataMatrix& A, DataMatrix& B) {
	vec_Str rowsIntersected;
	for (const Str& row : A.rowNames) {
		if (B.rowNameToIdx.count(row)) {
			rowsIntersected.push_back(row);
		}
	}
	A.keepOnlyFollowingRows(rowsIntersected);
	B.keepOnlyFollowingRows(rowsIntersected);
}



void intersectMatrixColumns(DataMatrix& A, DataMatrix& B) {
	vec_Str colsIntersected;
	for (const Str& col : A.colNames) {
		if (B.colNameToIdx.count(col)) {
			colsIntersected.push_back(col);
		}
	}
	A.keepOnlyFollowingColumns(colsIntersected);
	B.keepOnlyFollowingColumns(colsIntersected);
}