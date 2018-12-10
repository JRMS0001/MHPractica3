#pragma once

#include <fstream>
#include <string>
#include <vector>
#include "Structures.h"

#define INPUTS_FILE_PATH "InputsFile.txt"

class InputsFileReader
{
public:
	InputsFileReader();
	~InputsFileReader();

	ALGORITHM algorithm;
	CROSSOVER crossover;
	std::vector<std::string> inputs;
	std::vector<int> seeds;
	int populationSize;
	std::vector<int> iterationsBL;
	float probabilityCrossover;
	float probabilityMutation;
	int eliteNumber;

private:
	void readInputsFromDataFile(std::string path);
	size_t split(const std::string &txt, std::vector<std::string> &strs, char ch);
};
