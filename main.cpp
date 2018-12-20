#include <iostream>
#include <sstream>
#include <time.h>
#include <string>

#include "Instance.h"
#include "FileReader.h"
#include "InputsFileReader.h"
#include "Structures.h"

int main(int argc,const char * argv[]) {


	// Run the algorithm
	InputsFileReader *ifl= new InputsFileReader();

	for (std::string input : ifl->inputs) {

		std::cout << "For the file: " << input << "." << std::endl;

		std::string path = "";
		path.append("datos/").append(input);

		Instance* instance = new Instance(path);
		for(int seed : ifl->seeds){
			std::cout << "For the seed: " << seed << "." << std::endl;
			std::srand(seed);

			// Stationary
			if (ifl->algorithm == AGE) {
				std::cout << "Executing Stationary algorithm." << std::endl;
				int stationaryCost = 0;
				std::ofstream stationaryOutfile;
				if(ifl->crossover == OX)
					stationaryOutfile.open("logs/" + input + "_AGEOX_" + std::to_string(seed) + ".log");
				else
					stationaryOutfile.open("logs/" + input + "_AGEPMX_" + std::to_string(seed) + ".log");
				const clock_t stationary_begin_time = clock();
				instance->AGE(ifl->crossover, &stationaryCost, stationaryOutfile);
				std::cout << "Stationary execution time: " << float(clock() - stationary_begin_time) / CLOCKS_PER_SEC << std::endl;
				std::cout << "Stationary cost: " << stationaryCost << std::endl;
				stationaryOutfile.close();
			}

			// Generational
			else if (ifl->algorithm == AGG) {
				std::cout << "Executing Generational algorithm." << std::endl;
				int generationalCost = 0;
				std::ofstream generationalOutfile;
				if(ifl->crossover == OX && ifl->techniqueBL == ALL)
					generationalOutfile.open("logs/" + input + "_AGGOX_ALL_" + std::to_string(seed) + ".log");
				else if(ifl->crossover == OX && ifl->techniqueBL == TENPERCENT)
					generationalOutfile.open("logs/" + input + "_AGGOX_TENPERCENT_" + std::to_string(seed) + ".log");
				else if(ifl->crossover == OX && ifl->techniqueBL == ELITE)
					generationalOutfile.open("logs/" + input + "_AGGOX_ELITE_" + std::to_string(seed) + ".log");
				else if(ifl->crossover == PMX && ifl->techniqueBL == ALL)
					generationalOutfile.open("logs/" + input + "_AGGPMX_ALL_" + std::to_string(seed) + ".log");
				else if(ifl->crossover == PMX && ifl->techniqueBL == TENPERCENT)
					generationalOutfile.open("logs/" + input + "_AGGPMX_TENPERCENT_" + std::to_string(seed) + ".log");
				else if(ifl->crossover == PMX && ifl->techniqueBL == ELITE)
					generationalOutfile.open("logs/" + input + "_AGGPMX_ELITE_" + std::to_string(seed) + ".log");
				else{
					std::cout << "WARNING! Crossover or technique not supported, check the user manual." << std::endl;
				}

				const clock_t generational_begin_time = clock();
				instance->AGG(ifl->crossover, &generationalCost, generationalOutfile);
				std::cout << "Generational execution time: " << float(clock() - generational_begin_time) / CLOCKS_PER_SEC << std::endl;
				std::cout << "Generational cost: " << generationalCost << std::endl;
				generationalOutfile.close();
			}
		}
	}

	std::cout << "Execution over." << std::endl;

	//system("pause");
}
