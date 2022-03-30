#include "read_vrp.h"
#include "instance.h"
#include <vector>

using namespace std;

Instance::Instance(){}

Instance::Instance(const string& vrp_file){
    InstanceValues values = read_input_cvrp(vrp_file);
    this->n_cust = values.n_cust;
    this->cust_locations = values.cust_locations;
    this->Q = values.Q;
    this->dist_matrix = values.dist_matrix;
	this->sorted_dist_matrix = values.sorted_dist_matrix;
    //this->insertion_matrix = values.insertion_matrix;
    //this->remotion_matrix = values.remotion_matrix;
    this->dist_warehouse = values.dist_warehouse;
    this->demands = values.demands;
    this->max_demand = values.max_demand;
    this->mean_demands = values.mean_demands;

    this->depot = 0;
    this->customers_num = this->demands.size() - 1;
    this->customers_begin = 1;
    this->customers_end = this->demands.size();
}

Instance::~Instance(){}

float Instance::get_cost(int i, int j) const {

    return this->dist_matrix[i][j];

};

int Instance::get_demand(int i) const {
    return this->demands[i];
};
