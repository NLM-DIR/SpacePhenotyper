#ifndef __H__LIB__DATAMATRIX__
#define __H__LIB__DATAMATRIX__

#include "CommonTypes.h"



enum DataMatrix_CopyFlag : bool {
	FLAG_DataMatrix_COPYVALS_YES	= true,
	FLAG_DataMatrix_COPYVALS_NO		= false
};



enum DataMatrix_HeaderFlag : bool {
	FLAG_DataMatrix_HEADER_YES		= true,
	FLAG_DataMatrix_HEADER_NO		= false
};



struct DataMatrix {
	i32 numRows;
	i32 numCols;
	vec_Str rowNames;
	vec_Str colNames;
	map_Str_i32 rowNameToIdx;
	map_Str_i32 colNameToIdx;
	f64* data;
	f64** dataRows;

	void allocateDataMemory();
	void deallocateDataMemory();

	DataMatrix() : numRows(0), numCols(0), data(nullptr), dataRows(nullptr) {}
	DataMatrix(const DataMatrix& other, const DataMatrix_CopyFlag copyVals);
	DataMatrix(const vec_Str& _rowNames, const vec_Str& _colNames);
	DataMatrix(const Str& path, const char sep, const DataMatrix_HeaderFlag& headerFlag);

	f64* operator[] (const i32 rowIdx) { return dataRows[rowIdx]; }
	const f64* operator[] (const i32 rowIdx) const { return dataRows[rowIdx]; }

	f64& operator[] (const i32 rowIdx, const i32 colIdx) { return data[rowIdx*numCols + colIdx]; }
	const f64& operator[] (const i32 rowIdx, const i32 colIdx) const { return data[rowIdx*numCols + colIdx]; }
	
	f64& operator() (const i32 rowIdx, const i32 colIdx);
	const f64& operator() (const i32 rowIdx, const i32 colIdx) const;

	i32 getRowIdx(const Str& rowName) const;
	i32 getColIdx(const Str& colName) const;
	f64 getVal(const Str& rowName, const Str& colName) const;
	void printSize(const Str& matrixName) const;
	void validateNonNegativeValues(const Str& nameForPrint) const;

	void setVal(const Str& rowName, const Str& colName, const f64& x);
	void resetToZero();
	void transpose();
	void lexicographicallyOrderRows();
	void lexicographicallyOrderColumns();
	void lexicographicallyOrderRowsAndColumns();
	void keepOnlyRowsWithIndices(const vec_i32& indices);
	void keepOnlyFollowingRows(const vec_Str& keepRows);
	void keepOnlyColumnsWithIndices(const vec_i32& indices);
	void keepOnlyFollowingColumns(const vec_Str& keepCols);
	void makeColumnsFractions();
	void normalizeColumns(const f64& scale);

	void addMatrix(const DataMatrix& other);

	void readFromPath(const Str& path, const char sep, const bool noName = false);
	void readFromPathNoHeader(const Str& path, const char sep);

	void printToFilestreamAsTSV(FILE * const stream, const Str& firstHeaderCol, const DataMatrix* const B = nullptr) const;
	void printToFileAsTSV(const Str& filename, const Str& firstHeaderCol, const DataMatrix* const B = nullptr) const;
};



void intersectMatrixRows(DataMatrix& A, DataMatrix& B);
void intersectMatrixColumns(DataMatrix& A, DataMatrix& B);

#endif