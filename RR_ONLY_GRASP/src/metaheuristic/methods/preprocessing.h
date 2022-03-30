#pragma once
#include "../model/instance.h"
#include "input_parser.h"
#include <vector>
#include <string>

void print_features(std::vector<SampleInfo*>& dataset_train, Instance& instance, std::string& method, std::string& file_parameters, std::string& direc_output, int number_of_vnds);
void print_validation_features(std::vector<SampleInfo*>& dataset_train, Instance& instance, std::string& method, std::string& file_parameters, std::string& direc_output);
void print_vnd_features(std::vector<SampleInfo*>& dataset_train, Instance& instance, std::string& method, std::string& file_parameters, std::string& direc_output, int number_of_vnds);

