#include "Instance.h"


/* CONSTRUCTOR AND DESTRUCTOR */

Instance::Instance(std::string path) {
	inputsFileReader = new InputsFileReader();
	matricesFileReader = new FileReader(path);
	matrixSize = matricesFileReader->getMatrixSize();
}

Instance::~Instance() {
}






/* STATIONARY ALGORITHM */

int* Instance::AGE(CROSSOVER crossoverType, int * cost , std::ofstream &outfile ){
	int** flowMatrix = matricesFileReader->getFlowMatrix();
	int** distanceMatrix = matricesFileReader->getDistanceMatrix();

	std::vector<Element> population;

	// Generate the random population
	int* unitAndLocationAssociation = new int[matrixSize];
	for (int i=0; i< inputsFileReader->populationSize; i++){
		*cost=0;
		/*Random solution*/
		for(int j=0; j<matrixSize; j++){
			unitAndLocationAssociation[j]= j+1;
		}
		std::random_shuffle(&unitAndLocationAssociation[0],&unitAndLocationAssociation[matrixSize]); // matrixSize - 1 ?

		Element element;
		element.solution = unitAndLocationAssociation;
		element.cost = evaluateSolution(unitAndLocationAssociation);

		population.push_back(element);
	}
	std::sort(population.begin(),population.end(), &compareElements);
	outfile << "Best cost of the initial population: " << population.at(0).cost << std::endl;


	int it=1;
	while(it < 50000){
	outfile << "Iteration n°" << it << std::endl;

		/* SELECTION */

		//First binary tournament
		int r = rand() % inputsFileReader->populationSize;
		int s = rand() % inputsFileReader->populationSize;
		Element firstFather;
		Element secondFather;
		if(population.at(r).cost < population.at(s).cost ){
			firstFather = population.at(r);
		}
		else{
			firstFather = population.at(s);
		}

		//Second binary tournament
		r = rand() % inputsFileReader->populationSize;
		s = rand() % inputsFileReader->populationSize;
		if(population.at(r).cost < population.at(s).cost ){
			secondFather = population.at(r);
		}
		else{
			secondFather = population.at(s);
		}



		/* CROSSOVER */

		Element firstSon;
		Element secondSon;

		double pCrossover = rand() / (double)RAND_MAX; // between 0 and 1
		if (pCrossover <= PROB_CROSSOVER_STATIONARY) {

			// Initialization
			int* solutionSon1;
			int* solutionSon2;
			int intervalBegining;
			int intervalEnd;

			// Choosing interval beginning and end
			int intervalSize = 0;
			while (intervalSize <= 1) {
				int intervalPoint1 = (rand() % (matrixSize - 4)) + 2; // Between 2 and matrixSize-2
				int intervalPoint2 = (rand() % (matrixSize - 4)) + 2; // Between 2 and matrixSize-2

				if (intervalPoint1 < intervalPoint2) {
					intervalBegining = intervalPoint1;
					intervalEnd = intervalPoint2;
				}
				else {
					intervalBegining = intervalPoint2;
					intervalEnd = intervalPoint1;
				}
				intervalSize = intervalEnd - intervalBegining;
			}

			// Crossover
			if (crossoverType == PMX) {
				solutionSon1 = PMXCrossover(firstFather.solution, secondFather.solution, intervalBegining, intervalEnd);
				solutionSon2 = PMXCrossover(secondFather.solution, firstFather.solution, intervalBegining, intervalEnd);
			}
			else if (crossoverType == OX) {
				solutionSon1 = OXCrossover(firstFather.solution, secondFather.solution, intervalBegining, intervalEnd);
				solutionSon2 = OXCrossover(secondFather.solution, firstFather.solution, intervalBegining, intervalEnd);
			}
			else {
				solutionSon1 = nullptr;
				solutionSon2 = nullptr;
			}

			firstSon.solution = solutionSon1;
			firstSon.cost = evaluateSolution(solutionSon1);
			it++;
			secondSon.solution = solutionSon2;
			secondSon.cost = evaluateSolution(solutionSon2);
			it++;
		}
		else {
			firstSon = firstFather;
			secondSon = secondFather;
		}
				



		/* MUTATION */
		for (int gen = 0; gen < matrixSize; gen++) {
			double pMut = rand() / (double)RAND_MAX; // between 0 and 1
			if (pMut < inputsFileReader->probabilityMutation * (double)matrixSize) {
				int random = rand() % matrixSize;
				//Swapping elements in the solution
				int swap = firstSon.solution[gen];
				firstSon.solution[gen] = firstSon.solution[random];
				firstSon.solution[random] = swap;
				//Factorization
				firstSon.cost = evaluateSolution(firstSon.solution);
				it++;
			}
		}

		for (int gen = 0; gen < matrixSize; gen++) {
			double pMut = rand() / (double)RAND_MAX; // between 0 and 1
			if (pMut < inputsFileReader->probabilityMutation * (double)matrixSize) {
				int random = rand() % matrixSize;
				//Swapping elements in the solution
				int swap = secondSon.solution[gen];
				secondSon.solution[gen] = secondSon.solution[random];
				secondSon.solution[random] = swap;
				//Factorization
				secondSon.cost = evaluateSolution(secondSon.solution);
				it++;
			}
		}



		/* REPLACEMENT */
		std::sort(population.begin(),population.end(), &compareElements);

		std::vector<Element> replacementElements;
		// 2 worst elements
		replacementElements.push_back(population.at(inputsFileReader->populationSize - 2));
		replacementElements.push_back(population.at(inputsFileReader->populationSize - 1));
		// 2 new elements
		replacementElements.push_back(firstSon);
		replacementElements.push_back(secondSon);
		// Replacement
		std::sort(replacementElements.begin(), replacementElements.end(), &compareElements);
		population.at(inputsFileReader->populationSize - 2) = replacementElements.at(0);
		population.at(inputsFileReader->populationSize - 1) = replacementElements.at(1);


		// Display current population costs
		std::sort(population.begin(), population.end(), &compareElements);

		displayPopulationCosts(population, outfile);
		outfile << "Best solution's cost : " << population.at(0).cost << std::endl;
		outfile << std::endl;
	}

	/* RETURN */
	std::sort(population.begin(),population.end(), &compareElements);
	*cost= population.at(0).cost;
	return population.at(0).solution;
}






/* GENERATIONAL ALGORITHM */

int* Instance::AGG(CROSSOVER crossoverType, int * cost , std::ofstream &outfile ) {
	int** flowMatrix = matricesFileReader->getFlowMatrix();
	int** distanceMatrix = matricesFileReader->getDistanceMatrix();

	std::vector<Element> population;

	// Generate the random population
	int* unitAndLocationAssociation = new int[matrixSize];
	for (int i = 0; i < inputsFileReader->populationSize; i++) {
		*cost = 0;
		/*Random solution*/
		for (int j = 0; j < matrixSize; j++) {
			unitAndLocationAssociation[j] = j + 1;
		}
		std::random_shuffle(&unitAndLocationAssociation[0], &unitAndLocationAssociation[matrixSize]); // matrixSize - 1 ?

		Element element;
		element.solution = unitAndLocationAssociation;
		element.cost = evaluateSolution(unitAndLocationAssociation);

		population.push_back(element);
	}
	std::sort(population.begin(),population.end(), &compareElements);
	outfile << "Best cost of the initial population: " << population.at(0).cost << std::endl;


	int it = 1;
	while (it < 50000) {
		outfile << "Iteration n°" << it << std::endl;

		/* SELECTION */
		std::sort(population.begin(), population.end(), &compareElements);
		Element elite = population.at(0);
		//elite.cost = evaluateSolution(elite.solution);

		std::vector<Element> selectedPopulation;
		for (int i = 0; i < inputsFileReader->populationSize; i++) {
			int r = rand() % inputsFileReader->populationSize;
			int s = rand() % inputsFileReader->populationSize;
			Element element;
			if (population.at(r).cost < population.at(s).cost) {
				element = population.at(r);
			}
			else {
				element = population.at(s);
			}
			//element.cost = evaluateSolution(element.solution);
			selectedPopulation.push_back(element);
		}
		population = selectedPopulation;



		/* CROSSOVER */

		std::vector<Element> crossoveredPopulation;

		for (int i = 0; i < inputsFileReader->populationSize; i+=2) {

			Element firstSon;
			Element secondSon;

			double pCrossover = rand() / (double)RAND_MAX; // between 0 and 1
			if (pCrossover <= inputsFileReader->probabilityCrossover) {

				// Initialization
				int* solutionSon1;
				int* solutionSon2;
				int intervalBegining;
				int intervalEnd;

				// Choosing interval beginning and end
				int intervalSize = 0;
				while (intervalSize <= 1) {
					int intervalPoint1 = (rand() % (matrixSize - 4)) + 2; // Between 2 and matrixSize-2
					int intervalPoint2 = (rand() % (matrixSize - 4)) + 2; // Between 2 and matrixSize-2

					if (intervalPoint1 < intervalPoint2) {
						intervalBegining = intervalPoint1;
						intervalEnd = intervalPoint2;
					}
					else {
						intervalBegining = intervalPoint2;
						intervalEnd = intervalPoint1;
					}
					intervalSize = intervalEnd - intervalBegining;
				}

				// Crossover
				if (crossoverType == PMX) {
					solutionSon1 = PMXCrossover(population.at(i).solution, population.at(i + 1).solution, intervalBegining, intervalEnd);
					solutionSon2 = PMXCrossover(population.at(i + 1).solution, population.at(i).solution, intervalBegining, intervalEnd);
				}
				else if (crossoverType == OX) {
					solutionSon1 = OXCrossover(population.at(i).solution, population.at(i + 1).solution, intervalBegining, intervalEnd);
					solutionSon2 = OXCrossover(population.at(i + 1).solution, population.at(i).solution, intervalBegining, intervalEnd);
				}
				else {
					solutionSon1 = nullptr;
					solutionSon2 = nullptr;
				}
				
				firstSon.solution = solutionSon1;
				firstSon.cost = evaluateSolution(solutionSon1);
				it++;
				secondSon.solution = solutionSon2;
				secondSon.cost = evaluateSolution(solutionSon2);
				it++;
			}
			else {
				firstSon = population.at(i);
				secondSon = population.at(i + 1);
			}

			crossoveredPopulation.push_back(firstSon);
			crossoveredPopulation.push_back(secondSon);

		}
		population = crossoveredPopulation;



		/* MUTATION */

		for (int i = 0; i < inputsFileReader->populationSize; i++) {
			for (int gen = 0; gen < matrixSize; gen++) {
				double pMut = rand() / (double)RAND_MAX; // between 0 and 1
				if (pMut < inputsFileReader->probabilityMutation * (double)matrixSize) {
					int random = rand() % matrixSize;
					//Swapping elements in the solution
					int swap = population.at(i).solution[gen];
					population.at(i).solution[gen] = population.at(i).solution[random];
					population.at(i).solution[random] = swap;
					//Factorization
					population.at(i).cost = evaluateSolution(population.at(i).solution);
					it++;
				}
			}
		}



		/* REPLACEMENT */

		std::sort(population.begin(), population.end(), &compareElements);

		if (elite.cost < population.at(0).cost) {
			population.at(inputsFileReader->populationSize - 1) = elite;
		}
		else if (elite.cost > population.at(inputsFileReader->populationSize - 1).cost) {
			population.at(inputsFileReader->populationSize - 1) = elite;
		}
		else {
			bool containsElite = false;
			bool costEqualsElite = false;
			for (int i = 0; i < inputsFileReader->populationSize; i++) {
				if (population.at(i).cost == elite.cost) {
					costEqualsElite = true;
					if (areElementsEquals(population.at(i), elite)) {
						containsElite = true;
						break;
					}
				}
				else {
					if (costEqualsElite == true) {
						break;
					}
				}
			}
			if (!containsElite) {
				population.at(inputsFileReader->populationSize - 1) = elite;
			}
		}

		// Display current population costs
		std::sort(population.begin(), population.end(), &compareElements);

		displayPopulationCosts(population, outfile);
		outfile << "Best solution's cost : " << population.at(0).cost << std::endl;
		outfile << std::endl;

		it++;
	}



	/* RETURN */
	std::sort(population.begin(), population.end(), &compareElements);
	*cost = population.at(0).cost;
	return population.at(0).solution;
}






/* USEFUL FUNCTIONS */

int Instance::evaluateSolution(int* solution) {
	//Calculation of the cost of the generated solution
	int** flowMatrix = matricesFileReader->getFlowMatrix();
	int** distanceMatrix = matricesFileReader->getDistanceMatrix();

	int cost = 0;
	for (int i = 0; i < matrixSize; i++) {
		for (int j = 0; j < matrixSize; j++) {
			if (i != j) {
				cost += flowMatrix[i][j] * distanceMatrix[solution[i] - 1][solution[j] - 1];
			}	
		}
	}
	return cost;
}

bool Instance::compareElements(Element i, Element j) {
	return i.cost < j.cost;
}

bool Instance::areElementsEquals(Element firstElement, Element secondElement) {
	bool areEquals = true;
	for (int i = 0; i < matrixSize; i++) {
		if (!(firstElement.solution[i] == secondElement.solution[i])) {
			areEquals = false;
			break;
		}
	}
	return areEquals;
}

void Instance::displaySolution(int* solution) {
	for (int i = 0; i < matrixSize; i++) {
		std::cout << solution[i] << " ";
	}
	std::cout << std::endl;
}

void Instance::displayPopulationCosts(std::vector<Element> population, std::ofstream &outfile) {
	for (Element solution : population) {
		outfile << "Cost " << solution.cost << std::endl;
	}
	outfile << std::endl;
}






/* OX AND PMX CROSSOVERS */

int* Instance::OXCrossover(int* father, int* mother, int intervalBegining, int intervalEnd) {

	int intervalSize = intervalEnd - intervalBegining;
	if (intervalBegining >= 2 && intervalEnd <= matrixSize - 2) {

		/* VARIABLES */
		int *son = new int[matrixSize];
		std::vector<int> fatherIntervalValues;
		std::vector<int> motherSortedValues;
		bool* motherMask = new bool[matrixSize];
		// Initializing mother mask
		for (int i = 0; i < matrixSize; i++) {
			motherMask[i] = 0;
		}

		/* VALUES FROM FATHER TO SON */
		for (int i = intervalBegining; i < intervalEnd; i++) {
			son[i] = father[i];
			fatherIntervalValues.push_back(father[i]);
		}

		/* VALUES FROM MOTHER TO SON */

		// Fill in mother mask
		for (int value : fatherIntervalValues) {
			for (int i = 0; i < matrixSize; i++) {
				if (value == mother[i]) {
					motherMask[i] = 1;
					break;
				}
			}
		}

		// Taking (and sorting) values from mother
		for (int i = intervalEnd; i < matrixSize; i++) {
			if (motherMask[i] == 0) {
				motherSortedValues.push_back(mother[i]);
			}
		}
		for (int i = 0; i < intervalEnd; i++) {
			if (motherMask[i] == 0) {
				motherSortedValues.push_back(mother[i]);
			}
		}

		// Fill in son
		int iterationCounter = 0;
		for (int i = intervalEnd; i < matrixSize; i++) {
			son[i] = motherSortedValues[iterationCounter];
			iterationCounter++;
		}
		for (int i = 0; i < intervalBegining; i++) {
			son[i] = motherSortedValues[iterationCounter];
			iterationCounter++;
		}

		return son;

	}
	else {
		std::cout << "WARNING! The crossover interval is not correct." << std::endl;
		std::cout << "Interval beginning: " << intervalBegining << " and interval end: " << intervalEnd << std::endl;
		std::cout << "Interval size: " << intervalSize << " and matrix size: " << matrixSize << std::endl;
		return NULL;
	}
}



int* Instance::PMXCrossover(int* father, int* mother, int intervalBegining, int intervalEnd){

	int intervalSize = intervalEnd - intervalBegining;
	if (intervalBegining >= 2 && intervalEnd <= matrixSize - 2) {

		/* VARIABLES */
		int *son = new int[matrixSize];
		bool* motherMask = new bool[matrixSize];
		// Initializing mother mask
		for (int i = 0; i < matrixSize; i++) {
			motherMask[i] = 0;
		}

		/* VALUES COPIED AT THE SAME PLACE FROM FATHER TO SON */
		for (int i = intervalBegining; i < intervalEnd; i++) {
			son[i] = father[i];
		}

		/* VALUES COPIED FROM MOTHER TO SON */

		// Outside the father interval (on the left)
		for (int indexFather = 0; indexFather < intervalBegining; indexFather++) {
			// Outside the mother interval (on the left)
			for (int indexMother = 0; indexMother < intervalBegining; indexMother++) {
				if (father[indexFather] == mother[indexMother]) {
					son[indexMother] = mother[indexMother];
					break;
				}
			}
			// Outside the mother interval (on the right)
			for (int indexMother = intervalEnd; indexMother < matrixSize; indexMother++) {
				if (father[indexFather] == mother[indexMother]) {
					son[indexMother] = mother[indexMother];
					break;
				}
			}
			// Inside the mother interval
			for (int indexMother = intervalBegining; indexMother < intervalEnd; indexMother++) {
				if (father[indexFather] == mother[indexMother]) {
					int indexSon = indexMother;
					while (intervalBegining <= indexSon && indexSon < intervalEnd) {
						for (int i = 0; i < matrixSize; i++) {
							if (father[indexSon] == mother[i]) {
								indexSon = i;
								break;
							}
						}
					}
					son[indexSon] = mother[indexMother];
					break;
				}
			}
		}
		// Outside the father interval (on the right)
		for (int indexFather = intervalEnd; indexFather < matrixSize; indexFather++) {
			// Outside the mother interval (on the left)
			for (int indexMother = 0; indexMother < intervalBegining; indexMother++) {
				if (father[indexFather] == mother[indexMother]) {
					son[indexMother] = mother[indexMother];
					break;
				}
			}
			// Outside the mother interval (on the right)
			for (int indexMother = intervalEnd; indexMother < matrixSize; indexMother++) {
				if (father[indexFather] == mother[indexMother]) {
					son[indexMother] = mother[indexMother];
					break;
				}
			}
			// Inside the mother interval
			for (int indexMother = intervalBegining; indexMother < intervalEnd; indexMother++) {
				if (father[indexFather] == mother[indexMother]) {
					int indexSon = indexMother;
					while (intervalBegining <= indexSon && indexSon < intervalEnd) {
						for (int i = 0; i < matrixSize; i++) {
							if (father[indexSon] == mother[i]) {
								indexSon = i;
								break;
							}
						}
					}
					son[indexSon] = mother[indexMother];
					break;
				}
			}
		}

		return son;

	}
	else {
		std::cout << "WARNING! The crossover interval is not correct." << std::endl;
		std::cout << "Interval beginning: " << intervalBegining << " and interval end: " << intervalEnd << std::endl;
		std::cout << "Interval size: " << intervalSize << " and matrix size: " << matrixSize << std::endl;
		return NULL;
	}
}
