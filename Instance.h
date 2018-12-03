#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

#include "Structures.h"
#include "FileReader.h"
#include "InputsFileReader.h"

#define PROB_CROSSOVER_STATIONARY 1

class Instance {

private:
	int matrixSize;
	FileReader* matricesFileReader;
	InputsFileReader* inputsFileReader;

	void displaySolution(int* solution);
	void displayPopulationCosts(std::vector<Element> population, std::ofstream &outfile);
	bool checkDLB(int * DLB);
	int checkMove(int * sol, int i, int j);
	bool checkAlreadySearched(std::vector<int> alreadySearched, int number);
	
public:
	Instance(std::string path);
	~Instance();

	int * AGE(InputsFileReader *ifl, int * cost ,std::ofstream &outfile );
	int * AGG(InputsFileReader *ifl, int * cost ,std::ofstream &outfile );

	int * bestFirst(int * cost,int* initialSolution, int intensity,std::ofstream &outfile);

	int evaluateSolution(int* solution);
	static bool compareElements(Element i, Element j);
	bool areElementsEquals(Element firstElement, Element secondElement);

	// Interval : [ intervalBegining ; intervalEnding [
	int* PMXCrossover(int* father, int* mother, int intervalBegining, int intervalEnding);
	int* OXCrossover(int* father, int* mother, int intervalBegining, int intervalEnding);
};
