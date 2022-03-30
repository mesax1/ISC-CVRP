#pragma once
#include <vector>
#include <string>

class Instance {
    public:
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


    int depot = 0;
    int customers_num = 0;
    int customers_begin = 0;
    int customers_end = 0;

    Instance();
    Instance(const std::string& vrp_file);
    ~Instance();

    float get_cost(int i, int j) const;
    int get_demand(int i) const ;
};