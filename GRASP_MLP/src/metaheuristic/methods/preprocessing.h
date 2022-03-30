#pragma once
#include "../model/instance.h"
#include "input_parser.h"
#include <vector>
#include <string>

struct DTOP_values{
    std::vector<std::vector<float>> features;
    std::vector<std::string> features_names;
};

void print_features(std::vector<std::shared_ptr <SampleInfo>> dataset_train, Instance& instance, std::string& method, std::string& file_parameters, std::string& direc_output, int number_of_vnds);
void print_validation_features(std::vector<std::shared_ptr <SampleInfo>> dataset_train, Instance& instance, std::string& method, std::string& file_parameters, std::string& direc_output);
DTOP_values print_vnd_features(std::vector<std::shared_ptr <SampleInfo>> dataset_train, Instance& instance, std::string& method, std::string& file_parameters, std::string& direc_output, int number_of_vnds);
DTOP_values create_solution_features(std::vector<std::shared_ptr <SampleInfo>> dataset_train, Instance& instance, std::string& method);

