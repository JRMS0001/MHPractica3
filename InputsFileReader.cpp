#include "InputsFileReader.h"
#include <algorithm>
#include <iostream>


InputsFileReader::InputsFileReader()
{
	readInputsFromDataFile(INPUTS_FILE_PATH);
}

InputsFileReader::~InputsFileReader()
{
}

void InputsFileReader::readInputsFromDataFile(std::string path) {

	std::string line;
	std::ifstream myfile(path);
	std::vector<std::string> res;

	if (myfile.is_open()){
		getline(myfile, line);
		this->split(line,res,'=');
		if (res.at(1) == "AGE") {
			algorithm = AGE;
		}
		else if (res.at(1) == "AGG") {
			algorithm = AGG;
		}
		else {
			std::cout << "WARNING: The algorithm is not supported, only AGE and AGG are";
		}

		getline(myfile, line);
		this->split(line, res, '=');
		if (res.at(1) == "PMX") {
			crossover = PMX;
		}
		else if (res.at(1) == "OX") {
			crossover = OX;
		}
		else {
			std::cout << "WARNING: The crossover is not supported, only PMX and OX are";
		}

		getline(myfile, line);
		this->split(line, res, '=');
		std::string input = res.at(1);
		this->split(input, res, ' ');
		for (int i = 0; i < res.size(); i++) {
			inputs.push_back(res.at(i));
		}

		getline(myfile, line);
		this->split(line,res,'=');
		std::string input_seeds = res.at(1);
		this->split(input_seeds, res, ' ');
		for (int i = 0; i < res.size(); i++) {
			seeds.push_back(stoi(res.at(i)));
		}

		getline(myfile, line);
		this->split(line,res,'=');
		populationSize = stoi(res.at(1));

		getline(myfile, line);
		this->split(line, res, '=');
		if (res.at(1) == "ELITE") {
			techniqueBL = ELITE;
		}
		else if (res.at(1) == "TENPERCENT") {
			techniqueBL = TENPERCENT;
		}
		else if (res.at(1) == "ALL") {
			techniqueBL = ALL;
		}
		else {
			std::cout << "WARNING: The BL techinque is not supported, only ELITE, TENPERCENT and ALL are";
		}

		getline(myfile, line);
		this->split(line, res, '=');
		iterationsBL = stoi(res.at(1));

		/*
		getline(myfile, line);
		this->split(line, res, '=');
		std::string input_iterationsBL = res.at(1);
		this->split(input_iterationsBL, res, ' ');
		for (int i = 0; i < res.size(); i++) {
			iterationsBL.push_back(stoi(res.at(i)));
		}
		*/

		getline(myfile, line);
		this->split(line,res,'=');
		probabilityCrossover = stof(res.at(1));

		getline(myfile, line);
		this->split(line,res,'=');
		probabilityMutation = stof(res.at(1));

		getline(myfile, line);
		this->split(line,res,'=');
		eliteNumber = stoi(res.at(1));
	}

}

size_t InputsFileReader::split(const std::string &txt, std::vector<std::string> &strs, char ch) {
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


