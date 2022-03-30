#pragma once

#include <string>
#include <vector>

struct InstanceValues{
    int n_cust;
    std::vector<std::vector<int>> cust_locations;
    int Q;
    std::vector<std::vector<float>> dist_matrix;
	std::vector<std::vector<unsigned int>> sorted_dist_matrix;
    //std::vector<std::vector<std::vector<float>>> insertion_matrix;
    //std::vector<std::vector<std::vector<float>>> remotion_matrix;
    std::vector<float> dist_warehouse;
    std::vector<int> demands;
    int max_demand;
    float mean_demands;
};

InstanceValues read_input_cvrp(const std::string& filename);