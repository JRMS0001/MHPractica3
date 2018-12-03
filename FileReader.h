#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

class FileReader {
public:
	FileReader(std::string path);
	~FileReader();

	void displayMatrices();

	int getMatrixSize();
	int** getFlowMatrix();
	int** getDistanceMatrix();

private:
	void readMatrixFromDataFile(std::string path);
	size_t split(const std::string &txt, std::vector<std::string> &strs, char ch);
	
	int matrixSize;
	int** flowMatrix;
	int** distanceMatrix;
};

