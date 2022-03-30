/*
    Created on Sun Aug 2 2020
    @author: 
*/

#include "./read_vrp.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <string>
#include <cctype>
#include <cmath>
#include <boost/algorithm/string.hpp>

using namespace std;

template <typename T>
vector<unsigned int> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<unsigned int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

float compute_dist(float xi, float xj, float yi, float yj){
    //float exact_dist = sqrt(pow(xi - xj, 2) + pow(yi - yj, 2));
    float round_dist = round(sqrt((xi - xj) * (xi - xj) + (yi - yj) * (yi - yj)));
    //float round_dist = round(sqrt(pow(xi - xj, 2) + pow(yi - yj, 2)));
    return round_dist;
}

// Computes the distance matrix
vector<vector<float>> compute_distance_matrix(vector<int>& customers_x, vector<int>& customers_y) {
    int nb_customers = customers_x.size();
    vector<vector<float>> distance_matrix(nb_customers, vector<float>(nb_customers, NULL));

    for(int i = 0; i < nb_customers; i++){
        distance_matrix[i][i] = 0;
        for(int j = 0; j < nb_customers; j++){
            float dist = compute_dist(customers_x[i], customers_x[j], customers_y[i], customers_y[j]);
            distance_matrix[i][j] = dist;
            distance_matrix[j][i] = dist;
        }
    }

    return distance_matrix;
}

// Computes the insertion matrix
vector<vector<vector<float>>> compute_insertion_matrix(vector<vector<float>>& dist_matrix){
    int nb_customers = dist_matrix.size();
    vector<vector<vector<float>>> insertion_matrix(nb_customers, vector<vector<float>>(nb_customers, vector<float>(nb_customers)));
    for (int i = 0; i < nb_customers; i++) 
        { 
            for (int j = 0; j < nb_customers; j++) 
            { 
                for (int k = 0; k < nb_customers; k++) 
                { 
                    insertion_matrix[i][j][k] = -dist_matrix[i][j] + dist_matrix[i][k] + dist_matrix[k][j]; 
                } 
            } 
        }
    return insertion_matrix; 
}

// Computes the remotion matrix
vector<vector<vector<float>>> compute_remotion_matrix(vector<vector<float>>& dist_matrix){
    int nb_customers = dist_matrix.size();
    vector<vector<vector<float>>> remotion_matrix(nb_customers, vector<vector<float>>(nb_customers, vector<float>(nb_customers)));
    for (int i = 0; i < nb_customers; i++) 
        { 
            for (int j = 0; j < nb_customers; j++) 
            { 
                for (int k = 0; k < nb_customers; k++) 
                { 
                    remotion_matrix[i][j][k] = -dist_matrix[i][k] - dist_matrix[k][j] + dist_matrix[i][j]; 
                } 
            } 
        }
    return remotion_matrix; 
}


//Sort distance matrix for Giant tour
vector<vector<unsigned int>> compute_sorted_distance_matrix(vector<vector<float>>& dist_matrix){
	vector<vector<unsigned int>> sorted_distance_matrix(dist_matrix.size());
	for(int i = 0; i < dist_matrix.size(); i++){
		sorted_distance_matrix[i] = sort_indexes(dist_matrix[i]);
	}
	
	return sorted_distance_matrix;
}

// Computes the distances to warehouse
vector<float> compute_distance_warehouses(int depot_x, int depot_y, vector<int>& customers_x, vector<int>& customers_y){
    int nb_customers = customers_x.size();
    vector<float> distance_warehouses(nb_customers, NULL);

    for(int i = 0; i < nb_customers; i++){
        float dist = compute_dist(depot_x, customers_x[i], depot_y, customers_y[i]);
        distance_warehouses[i] = dist;
    }

    return distance_warehouses;
}


// This function uses boost library downloaded from: https://www.boost.org/users/history/version_1_73_0.html
vector<string> read_elem(const string& filename) {
    string line;
    vector<string> elements;
    vector<string> result;
    ifstream f{filename};
	
	if(!f){
		cout << "File not found" << endl;
	}else{
		while(!f.eof()){
			getline(f, line);
			boost::split(result, line, [](unsigned char x){return std::isspace(x);});
			for(int i = 0; i < result.size(); i++){
				if(result[i].size() != 0){
					elements.push_back(result[i]);
				}
			}
		}
	}
	
	f.close();

    return elements;
}

// The input files follow the "Augerat" format.
InstanceValues read_input_cvrp(const string& filename){
	cout << filename << endl;
    vector<string> elements = read_elem(filename);
	vector<string>::iterator ptr = elements.begin();

    int nb_customers = 0;
    int nb_nodes = 0;
    int truck_capacity;
    string token;

    while(1){
        token = *(next(ptr, 1));
		advance(ptr, 1);
        if(token == "DIMENSION"){
            //next(ptr, 1); // Removes the ":"
			advance(ptr, 1);
            nb_nodes = stoi(*next(ptr, 1));
			advance(ptr, 1);
            nb_customers = nb_nodes;
        }else if(token == "CAPACITY"){
            //next(ptr, 1); // Removes the ":"
			advance(ptr, 1);
            truck_capacity = stoi(*next(ptr, 1));
			advance(ptr, 1);
        }else if(token == "EDGE_WEIGHT_TYPE"){
            //next(ptr, 1); // Removes the ":"
			advance(ptr, 1);
            token = *next(ptr, 1);
			advance(ptr, 1);
            if(token != "EUC_2D"){
                cout << "Edge Weight Type " << token << " is not supported (only EUC_2D)" << endl;
                //exit(1);
            }
        }else if(token == "NODE_COORD_SECTION"){
            break;
        }
    }

    vector<int> customers_x(nb_customers, NULL);
    vector<int> customers_y(nb_customers, NULL);
    int depot_x = 0;
    int depot_y = 0;
    int node_id;
	
    for(int n = 0; n < nb_nodes; n++){
        node_id = stoi(*next(ptr, 1));
		advance(ptr, 1);
        
        if(node_id != (n + 1)){
            cout << "Unexpected index" << endl;
            exit(1);
        }

        if(node_id == 1){
            depot_x = stoi(*next(ptr, 1));
			advance(ptr, 1);
            depot_y = stoi(*next(ptr, 1));
			advance(ptr, 1);
            customers_x[node_id - 1] = depot_x;
            customers_y[node_id - 1] = depot_y; 
        }else{
            customers_x[node_id - 1] = stoi(*next(ptr, 1));
			advance(ptr, 1);
            customers_y[node_id - 1] = stoi(*next(ptr, 1));
			advance(ptr, 1); 
        }
    }

    // Compute distance matrix
    vector<vector<float>> distance_matrix = compute_distance_matrix(customers_x, customers_y);
	vector<vector<unsigned int>> sorted_dist_matrix = compute_sorted_distance_matrix(distance_matrix);
    //vector<vector<vector<float>>> insertion_matrix = compute_insertion_matrix(distance_matrix);
    //vector<vector<vector<float>>> remotion_matrix = compute_remotion_matrix(distance_matrix);
    vector<float> distance_warehouses = compute_distance_warehouses(depot_x, depot_y, customers_x, customers_y);
	
    token = *next(ptr, 1);
	advance(ptr, 1); 
	
    if(token != "DEMAND_SECTION"){
        cout << "Expected token DEMAND_SECTION" << endl;
        exit(1);
    }

    vector<int> demands(nb_customers, NULL);
    for(int n = 0; n < nb_nodes; n++){
        node_id = stoi(*next(ptr, 1));
		advance(ptr, 1);

        if(node_id != (n + 1)){
            cout << "Unexpected index" << endl;
            exit(1);
        }

        if(node_id == 1){
            if(stoi(*next(ptr, 1)) != 0){
                cout << "Demand for depot should be 0" << endl;
                exit(1);
            }
			advance(ptr, 1);
            demands[node_id - 1] = 0;
        }else{
            demands[node_id - 1] = stoi(*next(ptr, 1));
			advance(ptr, 1);
        }
    }

    token = *next(ptr, 1);
	advance(ptr, 1);
    if(token != "DEPOT_SECTION"){
        cout << "Expected token DEPOT_SECTION" << endl;
        exit(1);
    }

    int warehouse_id = stoi(*next(ptr, 1));
	advance(ptr, 1);
    if(warehouse_id != 1){
        cout << "Warehouse id is supposed to be 1" << endl;
        exit(1);
    }

    int end_of_depot_section = stoi(*next(ptr, 1));
	advance(ptr, 1);
    if(end_of_depot_section != -1){
        cout << "Expecting only one warehouse, more than one found" << endl;
        exit(1);
    }
    
    int max_demand = *max_element(demands.begin(), demands.end());
    float sum_demands = std::accumulate(demands.begin(), demands.end(), 0.0);
    float mean_demands = sum_demands / demands.size();

    vector<vector<int>> cust_locations = {customers_x, customers_y};
    
    return {nb_customers, cust_locations, truck_capacity, distance_matrix, sorted_dist_matrix,  distance_warehouses, demands, max_demand, mean_demands};
}

