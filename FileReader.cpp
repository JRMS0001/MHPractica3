#include "FileReader.h"


FileReader::FileReader(std::string path) {
	readMatrixFromDataFile(path);
}

FileReader::~FileReader()
{
}


int FileReader::getMatrixSize() {
	return matrixSize;
}

int** FileReader::getFlowMatrix() {
	return flowMatrix;
}

int** FileReader::getDistanceMatrix() {
	return distanceMatrix;
}

void FileReader::readMatrixFromDataFile(std::string path) {
	matrixSize = 0;
	flowMatrix = nullptr;
	distanceMatrix = nullptr;

	std::string line;
	std::ifstream myfile(path);

	if (myfile.is_open())
	{

		// Getting the matrix size
		getline(myfile, line);
		matrixSize = std::stoi(line);

		// Declaring the matrixes
		flowMatrix = new int*[matrixSize]; //col
		for (int i = 0; i < matrixSize; ++i)
			flowMatrix[i] = new int[matrixSize]; //row

		distanceMatrix = new int*[matrixSize]; //col
		for (int i = 0; i < matrixSize; ++i)
			distanceMatrix[i] = new int[matrixSize]; //row

		// Filling the matrixes in
		int lineNumber = 2;
		while (getline(myfile, line))
		{
			// Flow matrix
			if (3 <= lineNumber && lineNumber < 3 + matrixSize) {
				int i = 0;
				do {
					std::vector<std::string> v;
					v.clear();
					split(line, v, ' ');

					// Clearing blank data (due to double spaces)
					std::vector<std::string> valuesFlow;
					valuesFlow.clear();
					for (int j = 0; j < v.size(); j++) {
						if (!(v[j] == "")) {
							valuesFlow.push_back(v[j]);
						}
					}
					for (int j = 0; j < matrixSize; j++) {
						if (valuesFlow.size() != 0)
							flowMatrix[i][j] = std::stoi(valuesFlow[j]);
					}

					// New line
					getline(myfile, line);
					lineNumber++;
					i++;
				} while (i < matrixSize);
			}

			// Distance matrix
			if (4 + matrixSize <= lineNumber && lineNumber < 4 + matrixSize + matrixSize) {
				int i = 0;
				do {
					std::vector<std::string> v;
					v.clear();
					split(line, v, ' ');

					// Clearing blank data (due to double spaces)
					std::vector<std::string> valuesDistances;
					valuesDistances.clear();
					for (int j = 0; j < v.size(); j++) {
						if (!(v[j] == "")) {
							valuesDistances.push_back(v[j]);
						}
					}
					for (int j = 0; j < matrixSize; j++) {
						if (valuesDistances.size() != 0)
							distanceMatrix[i][j] = std::stoi(valuesDistances[j]);
					}

					// New line
					getline(myfile, line);
					lineNumber++;
					i++;
				} while (i < matrixSize);
			}

			lineNumber++;

		}
		myfile.close();

	}
	else {
		std::cout << "Unable to open file";
	}
}


size_t FileReader::split(const std::string &txt, std::vector<std::string> &strs, char ch) {
	size_t pos = txt.find(ch);
	size_t initialPos = 0;
	strs.clear();

	// Decompose statement
	while (pos != std::string::npos) {
		strs.push_back(txt.substr(initialPos, pos - initialPos));
		initialPos = pos + 1;

		pos = txt.find(ch, initialPos);
	}

	// Add the last one
	strs.push_back(txt.substr(initialPos, std::min(pos, txt.size()) - initialPos + 1));
	return strs.size();
}


void FileReader::displayMatrices() {
	if (matrixSize != 0) {
		for (int i = 0; i < matrixSize; i++) {
			for (int j = 0; j < matrixSize; j++) {
				std::cout << " " << flowMatrix[i][j];
			}
			std::cout << std::endl;
		}

		std::cout << std::endl;
		for (int i = 0; i < matrixSize; i++) {
			for (int j = 0; j < matrixSize; j++) {
				std::cout << " " << distanceMatrix[i][j];
			}
			std::cout << std::endl;
		}
	}
}
