/*
	Part of Ermin's lightweight CPP Library collection

	IO_RecordSeparatedFiles: Provides functions to read record-separated files and store the results into a 2D string vector.

	Version: 1.00
*/

#ifndef __header__cpplib__io__rsf__
#define __header__cpplib__io__rsf__

#include <cstdio>
#include <tuple>
#include <string>
#include <vector>
#include <ranges>

namespace __cpplib__io__rsf__ {



typedef int32_t										i32;
typedef std::string									Str;
typedef std::vector<Str>							vec_Str;
typedef std::vector<vec_Str>						vec2D_Str;
typedef std::ranges::iota_view<i32, i32>			irange;

/*	RETURN VALUES
	For a record-separated file, this function returns a triplet of i32 values representing, in order:
		The number of rows in the file
		The number of columns in the file
		The maximum length of any single record in the file
*/

/*	ERRORS AND EXCEPTIONS
	If the file does not have the same number of columns in every row (except the last row in the case that it is just a single empty line), an error message will be printed and the program will end.
*/

std::tuple<i32, i32, i32> getRSFInfo(
	FILE* const fin,
	const char recordSeparatorChar
) {
	if (recordSeparatorChar == '\n') {
		fprintf(stderr, "\n\nERROR in getRSFInfo: The record separator character cannot be the newline character. Exiting.\n");
		exit(1);
	}
	char cc;
	i32 maxRecordLen = 0;
	i32 currentRecordLen = 0;
	i32 currentNumSeps = 0;

	// Read the first line and identify the number of columns
	rewind(fin);
	while (true) {
		cc = getc(fin);
		if ((cc == '\n') or (cc == EOF) or (cc == recordSeparatorChar)) {
			// Reached the end of the current record
			if (currentRecordLen > maxRecordLen) maxRecordLen = currentRecordLen;
			currentRecordLen = 0;
			if (cc == recordSeparatorChar) currentNumSeps++;
			else break;	// Either reached EOF or '\n'
		}
		else {
			// Still within the record
			currentRecordLen++;
		}
	}
	const i32 numCols = currentNumSeps + 1;

	// Now read the remaining rows
	i32 numRows = 1;
	i32 lineLen = 0;
	currentNumSeps = 0;
	bool colNumMismatch = false;
	while (true) {
		cc = getc(fin);
		if ((cc == '\n') or (cc == EOF) or (cc == recordSeparatorChar)) {
			// Reached the end of the current record
			if (currentRecordLen > maxRecordLen) maxRecordLen = currentRecordLen;
			currentRecordLen = 0;
			if (cc == recordSeparatorChar) {
				currentNumSeps++;
				lineLen++;
			}
			else {
				if ((currentNumSeps != numCols - 1) and ((cc == '\n') or (lineLen > 0))) {
					colNumMismatch = true;
					break;
				}
				if ((cc == '\n') or (lineLen > 0)) {
					lineLen = 0;
					currentNumSeps = 0;
					numRows++;
				}
				if (cc == EOF) break;
			}
		}
		else {
			// Still within the record
			currentRecordLen++;
			lineLen++;
		}
	}
	if (colNumMismatch) {
		fprintf(stderr,
			"\n\nERROR in getRSFInfo: Row %d (0-indexed) of the record-separated file has %d columns, vs %d in the first line. Exiting.\n",
			numRows,
			currentNumSeps + 1,
			numCols
		);
		exit(1);
	}
	return {numRows, numCols, maxRecordLen};
}





std::tuple<i32, i32> getLineInfo(FILE* const fin) {
	char cc;
	i32 maxLineLen = 0;
	i32 lineLen = 0;
	i32 numRows = 0;
	rewind(fin);
	while (true) {
		cc = getc(fin);
		if ((cc == '\n') or (cc == EOF)) {
			// Reached the end of the current line
			if (lineLen > maxLineLen) maxLineLen = lineLen;
			if ((cc == '\n') or (lineLen > 0)) {
				lineLen = 0;
				numRows++;
			}
			lineLen = 0;
			if (cc == EOF) break;
		}
		else {
			// Still within the line
			lineLen++;
		}
	}
	return {numRows, maxLineLen};
}





vec2D_Str readSeparatedRecordsFromPath(
	const Str& filePath,
	const char recordSeparatorChar
) {
	char cc = 0;
	FILE* fin = fopen(filePath.c_str(), "r");
	if (fin == nullptr) {
		fprintf(stderr, "\n\nERROR in readSeparatedRecordsFromPath: Cannot open [%s] for reading. Exiting.\n", filePath.c_str());
		exit(1);
	}
	const auto rsfInfo = getRSFInfo(fin, recordSeparatorChar);
	const i32 numRows = std::get<0>(rsfInfo);
	const i32 numCols = std::get<1>(rsfInfo);
	const i32 maxRecordLen = std::get<2>(rsfInfo);
	const char* outString = "\n\nERROR in readSeparatedRecordsFromPath: File [%s] has 0 %s. Exiting.\n";
	if (numRows == 0) {
		fprintf(stderr, outString, filePath.c_str(), "rows");
		exit(1);
	}
	if (numCols == 0) {
		fprintf(stderr, outString, filePath.c_str(), "columns");
		exit(1);
	}

	// Allocate memory
	vec2D_Str table(numRows, vec_Str(numCols));
	char* temp = new char[maxRecordLen + 1];
	i32 len = 0;
	rewind(fin);
	for (const i32 i : irange(0, numRows)) {
		// Read the row data
		for (const i32 j : irange(0, numCols)) {
			len = 0;
			while (((cc = getc(fin)) != EOF) and (cc != recordSeparatorChar) and (cc != '\n')) temp[len++] = cc;
			if ((cc == '\n') and (len > 0) and (temp[len-1] == '\r')) len--;
			temp[len] = 0;
			table[i][j] = temp;
		}
		// Skip to the new line
		while ((cc != EOF) and (cc != '\n')) cc = getc(fin);
	}
	delete[] temp;

	fclose(fin);
	fprintf(stderr, "%d rows and %d columns read from [%s]\n", numRows, numCols, filePath.c_str());
	return table;
}





vec_Str readLinesFromPath(const Str& filePath) {
	char cc = 0;
	FILE* fin = fopen(filePath.c_str(), "r");
	if (fin == nullptr) {
		fprintf(stderr, "\n\nERROR in readLinesFromPath: Cannot open [%s] for reading. Exiting.\n", filePath.c_str());
		exit(1);
	}
	const auto lineInfo = getLineInfo(fin);
	const i32 numRows = std::get<0>(lineInfo);
	const i32 maxLineLen = std::get<1>(lineInfo);
	if (numRows == 0) {
		fprintf(stderr, "\n\nERROR in readLinesFromPath: File [%s] has 0 lines. Exiting.\n", filePath.c_str());
		exit(1);
	}

	// Allocate memory
	vec_Str records(numRows);
	char* temp = new char[maxLineLen + 1];
	i32 len;
	rewind(fin);
	for (const i32 i : irange(0, numRows)) {
		len = 0;
		while (((cc = getc(fin)) != EOF) and (cc != '\n')) temp[len++] = cc;
		if ((len > 0) and (temp[len-1] == '\r')) len--;
		temp[len] = 0;
		records[i] = temp;
	}
	delete[] temp;

	fclose(fin);
	fprintf(stderr, "%d rows read from [%s]\n", numRows, filePath.c_str());
	return records;
}



}

#endif