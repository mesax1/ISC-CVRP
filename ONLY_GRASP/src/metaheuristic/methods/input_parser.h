#pragma once
#include "../model/instance.h"
#include "../model/solution.h"
#include <string>
#include <iostream>
#include <fstream>

struct SampleInfo{
    std::string neighborhood;
	std::string neighborhood_1;
    std::string neighborhood_2;
    std::string neighborhood_3;
    std::string neighborhood_4;
    std::string neighborhood_5;
	int sample_number;
    std::shared_ptr <Solution> tsp_solution;
    std::shared_ptr <Solution> initial_solution;
    float final_cost;
    float final_cost_1;
    float final_cost_2;
    float final_cost_3;
    float final_cost_4;
    float final_cost_5;
    std::vector<float> final_costs;
    std::vector<int> neighborhoods;

};


std::vector<SampleInfo*> read_dataset(Instance& instance, std::string& dataset_file);
std::vector<SampleInfo*> read_select_vnd_dataset(Instance& instance, std::string& dataset_file, int number_of_vnds);
std::vector<SampleInfo*> read_validation_dataset(Instance& instance, std::string& dataset_file);
std::vector<int> read_classification(std::string& dataset_file);
