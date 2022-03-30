#include "../model/solution.h"
#include "preprocessing.h"
#include <iostream>
#include <fstream>
#include <boost/geometry.hpp>
#include <boost/geometry/arithmetic/cross_product.hpp>
#include <boost/geometry/arithmetic/arithmetic.hpp>
#include <numeric>
#include <cmath>
#include <map>
#include <algorithm> 

using namespace std;
namespace bg = boost::geometry;

struct DDOP_values{
	vector<int> solution_intersections;
	vector<vector<float>> longest_customers_distances;
	vector<vector<float>> smallest_customers_distances;
	vector<vector<float>> depot_customers_distances;
	vector<float> all_average_gravity;
	vector<vector<float>> all_routes_gravity_x;
	vector<vector<float>> all_routes_gravity_y;
	vector<vector<vector<float>>> all_central_line_distances;
	vector<vector<vector<float>>> all_routes_radians;
	vector<vector<vector<float>>> all_gravity_radians;
	vector<vector<vector<float>>> all_route_depot_distances;
	vector<vector<int>> all_routes_customers;
	vector<vector<int>> all_new_parameters;
	vector<string> new_parameters_names;
	vector<vector<float>> all_routes_distances;
	vector<vector<float>> depot_customers_distances_divided;
	vector<vector<float>> longest_interior_distances;
	vector<vector<float>> smallest_interior_distances;
	vector<vector<float>> all_first_last_cust_demands;
	vector<vector<float>> all_furthest_client_demand;
	vector<vector<float>> all_degree_of_neighborhood;
	vector<vector<float>> all_routes_degree_of_neighborhood;
	vector<vector<float>> all_routes_two_don;
    vector<float> all_routes_average_tsp_area ;
    vector<float> all_routes_largest_tsp_area ;
    vector<float> all_routes_smallest_tsp_area ;

};

/*
struct DTOP_values{
	vector<vector<float>> features;
	vector<string> features_names;
};
 */

struct STOP_values{
	vector<float> features;
	vector<string> features_names;
};

template<typename Iter_T>
long double vectorNorm(Iter_T first, Iter_T last) {
  return sqrt(inner_product(first, last, first, 0.0L));
}
  
// Given three colinear points p, q, r, the function checks if 
// point q lies on line segment 'pr' 
bool onSegment(Point p, Point q, Point r) { 
    if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) && 
        q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y)) 
       return true; 
  
    return false; 
} 
  
// To find orientation of ordered triplet (p, q, r). 
// The function returns following values 
// 0 --> p, q and r are colinear 
// 1 --> Clockwise 
// 2 --> Counterclockwise 
int orientation(Point p, Point q, Point r) { 
    // See https://www.geeksforgeeks.org/orientation-3-ordered-points/ 
    // for details of below formula. 
    int val = (q.y - p.y) * (r.x - q.x) - 
              (q.x - p.x) * (r.y - q.y); 
  
    if (val == 0) return 0;  // colinear 
  
    return (val > 0)? 1: 2; // clock or counterclock wise 
} 
  
// The main function that returns true if line segment 'p1q1' 
// and 'p2q2' intersect. 
bool doIntersect(Point p1, Point q1, Point p2, Point q2) { 
    // Find the four orientations needed for general and 
    // special cases 
    int o1 = orientation(p1, q1, p2); 
    int o2 = orientation(p1, q1, q2); 
    int o3 = orientation(p2, q2, p1); 
    int o4 = orientation(p2, q2, q1); 
  
    // General case 
    if (o1 != o2 && o3 != o4) 
        return true; 
  
    // Special Cases 
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1 
    if (o1 == 0 && onSegment(p1, p2, q1)) return true; 
  
    // p1, q1 and q2 are colinear and q2 lies on segment p1q1 
    if (o2 == 0 && onSegment(p1, q2, q1)) return true; 
  
    // p2, q2 and p1 are colinear and p1 lies on segment p2q2 
    if (o3 == 0 && onSegment(p2, p1, q2)) return true; 
  
     // p2, q2 and q1 are colinear and q1 lies on segment p2q2 
    if (o4 == 0 && onSegment(p2, q1, q2)) return true; 
  
    return false; // Doesn't fall in any of the above cases 
}

float sum(vector<float>& array_to_sum){
	float sumatory = 0;
	for(int i = 0; i < array_to_sum.size(); i++){
		sumatory += array_to_sum[i];
	}
	
	return sumatory;
}

double sum(vector<double>& array_to_sum){
	double sumatory = 0;
	for(int i = 0; i < array_to_sum.size(); i++){
		sumatory += array_to_sum[i];
	}
	
	return sumatory;
}

int sum_int(vector<int>& array_to_sum){
	int sumatory = 0;
	for(int i = 0; i < array_to_sum.size(); i++){
		sumatory += array_to_sum[i];
	}
	
	return sumatory;
}

float polygonArea(vector<float> X, vector<float> Y, int n){ 
    // Initialze area 
    float area = 0.0; 
  
    // Calculate value of shoelace formula 
    int j = n - 1; 
    for (int i = 0; i < n; i++) 
    { 
        area += (X[j] + X[i]) * (Y[j] - Y[i]); 
        j = i;  // j is previous vertex to i 
    } 
  
    // Return absolute value 
    return abs(area / 2.0); 
} 

float desv_est(vector<float>& v_desv){
	float sum = std::accumulate(v_desv.begin(), v_desv.end(), 0.0); 
	float mean = sum/v_desv.size(); 

	std::vector<float> diff(v_desv.size()); 
	std::transform(v_desv.begin(), v_desv.end(), diff.begin(), [mean](float x) { return x - mean; }); 
	float sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0); 
	return std::sqrt(sq_sum/v_desv.size());
}

DTOP_values dataset_determine_client_pairs_parameters(vector<shared_ptr <SampleInfo>> dataset_train, Instance& instance){
	
	vector<vector<float>> all_new_parameters = vector<vector<float>>();
	vector<string> new_parameters_names = vector<string>();
	
	std::shared_ptr <Solution> solution;
	map<string, float> new_parameters;
	
	for(int i = 0; i < dataset_train.size(); i++){
		solution = dataset_train[i]->initial_solution;
		int n_vehicles = solution->n_vehicles;
		
		new_parameters = {};
		string key;
		for(int x = 1; x < instance.n_cust; x++){
			for(int y = x; y < instance.n_cust; y++){
				if(x != y){
					key = to_string(x) + "-" + to_string(y);
					new_parameters[key] = 0.0;
				}
			}
		}
		
		std::shared_ptr <Vehicle> vehicle;
		for(int n = 0; n < n_vehicles; n++){
			vehicle = solution->vehicles[n];
			if(vehicle->customers.size() > 2){
				for(const auto& customer : vehicle->customers){
					for(const auto& customer_2 : vehicle->customers){
						if(customer->id != customer_2->id && customer_2->id != 0 && customer->id != 0){
							key = to_string(customer->id) + "-" + to_string(customer_2->id);
							if(new_parameters.count(key) == 1){
								new_parameters[key] = 1.0;
							}
						}
					}
				}
			}
		}
	
		vector<float> values;
		
		for(std::map<string, float>::iterator it = new_parameters.begin(); it != new_parameters.end(); ++it) {
			values.push_back(it->second);
		}
		all_new_parameters.push_back(values);
	}
	
	
	for(std::map<string, float>::iterator it = new_parameters.begin(); it != new_parameters.end(); ++it) {
		new_parameters_names.push_back(it->first);
	}
	
	vector<vector<float>> route_parameters;
	vector<float> parameters_values;
	for(int i = 0; i < dataset_train.size(); i++){
		parameters_values = {};
		parameters_values.insert(parameters_values.end(), all_new_parameters[i].begin(), all_new_parameters[i].end());
		route_parameters.push_back(parameters_values);
	}
	
	cout << "Pasa de definir route_parameters" << endl;
					
	return {route_parameters, new_parameters_names};
}

DDOP_values dataset_determine_vrp_parameters(vector<shared_ptr <SampleInfo>> dataset_train, Instance& instance){
	DDOP_values resulting_data;
	
	// Declaration of lists
	vector<int> solution_intersections = vector<int>();
	vector<vector<float>> longest_customers_distances = vector<vector<float>>();
	vector<vector<float>> smallest_customers_distances = vector<vector<float>>();
	vector<vector<float>> longest_interior_distances = vector<vector<float>>();
	vector<vector<float>> smallest_interior_distances = vector<vector<float>>();
	vector<vector<float>> depot_customers_distances = vector<vector<float>>();
	vector<float> all_average_gravity = vector<float>();
	
	vector<vector<float>> all_routes_gravity_x = vector<vector<float>>();
	vector<vector<float>> all_routes_gravity_y = vector<vector<float>>();
	vector<vector<vector<float>>> all_central_line_distances = vector<vector<vector<float>>>();
	vector<vector<vector<float>>> all_routes_radians = vector<vector<vector<float>>>();
	vector<vector<vector<float>>> all_gravity_radians = vector<vector<vector<float>>>();
	vector<vector<vector<float>>> all_route_depot_distances = vector<vector<vector<float>>>();
	vector<vector<int>> all_routes_customers = vector<vector<int>>();
	vector<vector<int>> all_new_parameters = vector<vector<int>>();
	vector<string> new_parameters_names = vector<string>();
	vector<vector<float>> all_routes_distances = vector<vector<float>>();
	vector<vector<float>> depot_customers_distances_divided = vector<vector<float>>();
	vector<vector<float>> all_first_last_cust_demands = vector<vector<float>>();
	vector<vector<float>> all_furthest_client_demand = vector<vector<float>>();
	vector<vector<float>> all_degree_of_neighborhood = vector<vector<float>>(); // Revisar
	vector<vector<float>> all_routes_degree_of_neighborhood = vector<vector<float>>();
	vector<vector<float>> all_routes_two_don = vector<vector<float>>();
    vector<float> all_routes_average_tsp_area = vector<float>();
    vector<float> all_routes_largest_tsp_area = vector<float>();
    vector<float> all_routes_smallest_tsp_area = vector<float>();
	
	std::shared_ptr <Solution> solution;
	vector<vector<float>> dist_matrix;
	vector<vector<unsigned int>> D;
	vector<vector<int>> position_matrix;
	vector<float> degree_of_neighborhood;
	map<string, int> new_parameters;
	
	for(int i = 0; i < dataset_train.size(); i++){
		solution = dataset_train[i]->initial_solution;
		dist_matrix = instance.dist_matrix;
		//position_matrix is used to get the degree of neighborhood
		D = instance.sorted_dist_matrix;
		position_matrix = vector<vector<int>>(D.size(), vector<int>(D.size(), 0));
		degree_of_neighborhood = vector<float>(D.size(), 0); //revisar
		
		for(int u = 0; u < position_matrix.size(); u++){
			for(int v = 0; v < position_matrix.size(); v++){
				unsigned int position = D[u][v];
				position_matrix[u][position] = v;
			}
		}
		
		vector<vector<pair<std::shared_ptr <Customer>, std::shared_ptr <Customer>>>> vehicles_edge_list = vector<vector<pair<std::shared_ptr <Customer>, std::shared_ptr <Customer>>>>();
		int j = 0;
		int n_vehicles = solution->vehicles.size();
		
		vector<float> longest_dist_per_route(n_vehicles);
		vector<float> smallest_dist_per_route(n_vehicles);
		vector<float> longest_interior_dist_per_route(n_vehicles);
		vector<float> smallest_interior_dist_per_route(n_vehicles);
		
		vector<float> average_depot_distance(n_vehicles);
		vector<float> average_depot_distance_divided(n_vehicles);
		vector<float> routes_gravity_x(n_vehicles);
		vector<float> routes_gravity_y(n_vehicles);
		
		vector<vector<float>> routes_central_line_distance(n_vehicles);
		vector<vector<float>> each_route_radians(n_vehicles);
		vector<vector<float>> each_gravity_radians(n_vehicles);
		vector<vector<float>> each_route_depot_distance(n_vehicles);


		
		vector<int> route_customers(n_vehicles);
		new_parameters = {};
		
		for(int x = 1; x < instance.n_cust; x++){
			for(int y = x; y < instance.n_cust; y++){
				if(x != y){
					string key = x + "-" + y;
					new_parameters[key] = 0;
				}
			}
		}
		
		vector<float> each_route_distance(n_vehicles);
		vector<float> average_first_last_cust_demand(n_vehicles);
		vector<float> demand_furthest(n_vehicles);
		vector<float> each_route_don(n_vehicles);
		vector<float> each_route_two_don(n_vehicles);
		
		std::shared_ptr <Vehicle> vehicle;
		std::shared_ptr <Customer> curr_customer;
		std::shared_ptr <Customer> prev_customer;

        vector<float> each_route_tsp_area = vector<float>();
		
		for(int n = 0; n < n_vehicles; n++){
			vehicle = solution->vehicles[n];
			if(vehicle->customers.size() > 2){
				vehicles_edge_list.push_back(vector<pair<std::shared_ptr <Customer>, std::shared_ptr <Customer>>>());
				
				curr_customer = nullptr;
				prev_customer = nullptr;
				vector<float> edge_features = vector<float>();
				vector<float> first_and_last_cust = vector<float>();
				vector<float> interior_edge_features = vector<float>();
				vector<float> depot_features = vector<float>();
				
				float furthest_distance = 0;
				float furthest_client_demand = 0;
				
				float gravity_center_x = 0;
				float gravity_center_y = 0;
				float vehicle_cost = 0;
                vector<float> vector_x = vector<float>();
                vector<float> vector_y = vector<float>();
                float tsp_area = 0;
				
				vector<int> vehicle_customers_list = vector<int>();
				float curr_distance;
				
				for(const auto& customer : vehicle->customers){
					curr_customer = customer;
					curr_distance = dist_matrix[0][curr_customer->id];
					vehicle_customers_list.push_back(curr_customer->id);
					
					if(prev_customer != nullptr){
						gravity_center_x += curr_customer->x;
						gravity_center_y += curr_customer->y;
						vehicles_edge_list[j].push_back(make_pair(prev_customer, curr_customer));
						vehicle_cost += dist_matrix[prev_customer->id][curr_customer->id];

                        vector_x.push_back(curr_customer->x);
                        vector_y.push_back(curr_customer->y);
						
						if(prev_customer->id == 0){
							float distance = dist_matrix[prev_customer->id][curr_customer->id];
							depot_features.push_back(distance);
							edge_features.push_back(distance);
							first_and_last_cust.push_back(curr_customer->demand);
						}
						
						if(curr_customer->id == 0){
							float distance = dist_matrix[prev_customer->id][curr_customer->id];
							depot_features.push_back(distance);
							edge_features.push_back(distance);
							first_and_last_cust.push_back(prev_customer->demand);
						}
						
						if(curr_distance > furthest_distance){
							furthest_distance = curr_distance;
							furthest_client_demand = curr_customer->demand;
						}
					}
					
					if(prev_customer != nullptr && prev_customer->id != 0 && curr_customer->id != 0){
						float distance = dist_matrix[prev_customer->id][curr_customer->id];
						interior_edge_features.push_back(distance);
						edge_features.push_back(distance);
					}
					
					prev_customer = curr_customer;
				}
				
				j++;
				average_depot_distance_divided[n] = sum(depot_features)/vehicle_cost;
				average_depot_distance[n] = sum(depot_features);
				average_first_last_cust_demand[n] = sum(first_and_last_cust)/2;
				demand_furthest[n] = furthest_client_demand;
				
				vector<float> clients_don(vehicle_customers_list.size() - 2);
				for(int i = 1; i < vehicle_customers_list.size() - 1; i++){
					clients_don[i-1] = (position_matrix[vehicle_customers_list[i-1]][vehicle_customers_list[i]] + position_matrix[vehicle_customers_list[i]][vehicle_customers_list[i+1]])/2;
				}
				
				each_route_two_don[n] = sum(clients_don)/clients_don.size();
				
				if(interior_edge_features.size() > 0){
					longest_dist_per_route[n] = (*max_element(edge_features.begin(), edge_features.end()))/vehicle_cost;
					smallest_dist_per_route[n] = (*min_element(edge_features.begin(), edge_features.end()))/vehicle_cost;
					longest_interior_dist_per_route[n] = (*max_element(interior_edge_features.begin(), interior_edge_features.end()))/vehicle_cost;
					smallest_interior_dist_per_route[n] = (*min_element(interior_edge_features.begin(), interior_edge_features.end()))/vehicle_cost;
				}
				
				int route_cardinality = vehicle->customers.size() - 1;
				
				gravity_center_x = gravity_center_x/route_cardinality;
				gravity_center_y = gravity_center_y/route_cardinality;
				routes_gravity_x[n] = gravity_center_x;
				routes_gravity_y[n] = gravity_center_y;
				each_route_distance[n] = vehicle_cost;
				
				vector<float> central_line_distance = vector<float>();
				vector<float> route_radians = vector<float>();
				vector<float> gravity_radians = vector<float>();
				vector<float> each_depot_distance = vector<float>();
				vector<float> route_degree_of_neighborhood = vector<float>();
				
				int depot_id = 0;
				Point depot_pos = vehicle->customers[0]->pos;
				Point grav_center_pos = {gravity_center_x, gravity_center_y};
				
				for(const auto& customer : vehicle->customers){
					
					bg::model::point<double, 2, bg::cs::cartesian> vector_1;
					bg::model::point<double, 2, bg::cs::cartesian> vector_grav;
					float grav_dot_product;
					float angle_grav;
					
					if(customer->id != 0){
						bg::model::point<double, 2, bg::cs::cartesian> depot_pos_p(depot_pos.x, depot_pos.y);
						bg::model::point<double, 2, bg::cs::cartesian> grav_center_pos_p(grav_center_pos.x, grav_center_pos.y);
						bg::model::point<double, 2, bg::cs::cartesian> customer_pos_p(customer->pos.x, customer->pos.y);
						//cout << customer_pos_p.get<0>() << " " << customer_pos_p.get<1>() << endl;
						//cout << depot_pos_p.get<0>() << " " << depot_pos_p.get<1>() << endl;
						bg::model::point<double, 2, bg::cs::cartesian> grav_minus_depot = grav_center_pos_p;
						bg::model::point<double, 2, bg::cs::cartesian> customer_minus_depot = customer_pos_p;
						
						bg::subtract_point(grav_minus_depot, depot_pos_p);
						bg::subtract_point(customer_minus_depot, depot_pos_p);
						
						bg::model::point<double, 2, bg::cs::cartesian> cross_product_result;
						//cout << grav_minus_depot.get<0>() << " " << grav_minus_depot.get<1>() << endl;
						//cout << customer_minus_depot.get<0>() << " " << customer_minus_depot.get<1>() << endl;
						cross_product_result = bg::cross_product(grav_minus_depot, customer_minus_depot);
						
						vector<double> grav_minus_depot_vector = {grav_minus_depot.get<0>(), grav_minus_depot.get<1>()};
						vector<double> customer_minus_depot_vector = {customer_minus_depot.get<0>(), customer_minus_depot.get<1>()};
						
						//cout << cross_product_result.get<0>() << endl;
						central_line_distance.push_back(cross_product_result.get<0>()/vectorNorm(grav_minus_depot_vector.begin(), grav_minus_depot_vector.end()));
						
						vector_1 = customer_minus_depot;
						bg::divide_value(vector_1, vectorNorm(customer_minus_depot_vector.begin(), customer_minus_depot_vector.end()));
						
						vector_grav = grav_minus_depot;
						bg::divide_value(vector_grav, vectorNorm(grav_minus_depot_vector.begin(), grav_minus_depot_vector.end()));
						
						grav_dot_product = bg::dot_product(vector_1, vector_grav);
						angle_grav = acos(grav_dot_product);
						gravity_radians.push_back(angle_grav);
						
						float distance = dist_matrix[depot_id][customer->id];
						each_depot_distance.push_back(distance);
					}
					
					vector<int> curr_neigh_degree = vector<int>();
					for(const auto& customer_2 : vehicle->customers){
						if(customer->id != customer_2->id && customer->id != 0){
							curr_neigh_degree.push_back(position_matrix[customer->id][customer_2->id]);
						}
						
						if(customer->id != customer_2->id && customer_2->id != 0 && customer->id != 0){
							bg::model::point<double, 2, bg::cs::cartesian> depot_pos_p(depot_pos.x, depot_pos.y);
							bg::model::point<double, 2, bg::cs::cartesian> customer_2_minus_depot(customer_2->pos.x, customer_2->pos.y);
							bg::subtract_point(customer_2_minus_depot, depot_pos_p);
							
							vector<double> customer_2_minus_depot_vector = {customer_2_minus_depot.get<0>(), customer_2_minus_depot.get<1>()};
							
							bg::model::point<double, 2, bg::cs::cartesian> vector_2 = customer_2_minus_depot;
							bg::divide_value(vector_2, vectorNorm(customer_2_minus_depot_vector.begin(), customer_2_minus_depot_vector.end()));
							
							float dot_product = bg::dot_product(vector_1, vector_2);
							float angle = acos(dot_product);
							route_radians.push_back(angle);
							string key = customer->id + "-" + customer_2->id;
							if(new_parameters.count(key) == 1){
								new_parameters[key] = 1;
							}
						}
					}
					
					degree_of_neighborhood[customer->id] = sum_int(curr_neigh_degree)/(vehicle->customers.size() - 1);
					route_degree_of_neighborhood.push_back(degree_of_neighborhood[customer->id]);
				}
				
				each_route_radians[n] = route_radians;
				routes_central_line_distance[n] = central_line_distance;
				each_gravity_radians[n] = gravity_radians;
				each_route_depot_distance[n] = each_depot_distance;
				route_customers[n] = vehicle->customers.size() - 2;
				each_route_don[n] = sum(route_degree_of_neighborhood)/route_degree_of_neighborhood.size();
                tsp_area = polygonArea(vector_x, vector_y, vector_x.size());

                if (tsp_area > 0) {
                    each_route_tsp_area.push_back(tsp_area);
                }
			}
		}
		
		vector<float> gravity_distances = vector<float>();
		for(int i = 0; i < routes_gravity_x.size(); i++){
			for(int j = 0; j < routes_gravity_x.size(); j++){
				if(i != j){
					bg::model::point<double, 2, bg::cs::cartesian> a(routes_gravity_x[i], routes_gravity_y[i]);
					bg::model::point<double, 2, bg::cs::cartesian> b(routes_gravity_x[j], routes_gravity_y[j]);
					bg::subtract_point(a, b);
					
					vector<double> a_vector = {a.get<0>(), a.get<1>()};
					gravity_distances.push_back(vectorNorm(a_vector.begin(), a_vector.end())); 
				}
			}
		}

		all_routes_average_tsp_area.push_back(sum(each_route_tsp_area)/solution->vehicles.size());
        all_routes_largest_tsp_area.push_back(*max_element(each_route_tsp_area.begin(), each_route_tsp_area.end()));
        all_routes_smallest_tsp_area.push_back(*min_element(each_route_tsp_area.begin(), each_route_tsp_area.end()));



		all_average_gravity.push_back(sum(gravity_distances)/solution->vehicles.size());
		
		all_routes_customers.push_back(route_customers);
		all_route_depot_distances.push_back(each_route_depot_distance);
		all_gravity_radians.push_back(each_gravity_radians);
		all_routes_radians.push_back(each_route_radians);
		all_central_line_distances.push_back(routes_central_line_distance);
		all_routes_gravity_x.push_back(routes_gravity_x);
		all_routes_gravity_y.push_back(routes_gravity_y);
		longest_customers_distances.push_back(longest_dist_per_route);
		smallest_customers_distances.push_back(smallest_dist_per_route);
		longest_interior_distances.push_back(longest_interior_dist_per_route);
		smallest_interior_distances.push_back(smallest_interior_dist_per_route);
		depot_customers_distances.push_back(average_depot_distance);
		depot_customers_distances_divided.push_back(average_depot_distance_divided);
		int number_intersections = 0;
		all_routes_distances.push_back(each_route_distance);
		all_first_last_cust_demands.push_back(average_first_last_cust_demand);
		all_furthest_client_demand.push_back(demand_furthest);
		all_degree_of_neighborhood.push_back(degree_of_neighborhood);
		all_routes_degree_of_neighborhood.push_back(each_route_don);
		all_routes_two_don.push_back(each_route_two_don);
		
		for(int i = 0; i < vehicles_edge_list.size() - 1; i++){
			for(const auto& edges : vehicles_edge_list[i]){
				for(int j = 0; j < vehicles_edge_list.size(); j++){
					if(j != i){
						for(const auto& other_edges : vehicles_edge_list[j]){
							if(doIntersect(edges.first->pos, edges.second->pos, other_edges.first->pos, other_edges.second->pos)){
								number_intersections++;
							}
						}
					}
				}
			}
		}
		
		solution_intersections.push_back(number_intersections/2);
		vector<int> values;
		
		for(std::map<string,int>::iterator it = new_parameters.begin(); it != new_parameters.end(); ++it) {
			values.push_back(it->second);
		}
		all_new_parameters.push_back(values);
	}
	
	
	for(std::map<string,int>::iterator it = new_parameters.begin(); it != new_parameters.end(); ++it) {
		new_parameters_names.push_back(it->first);
	}
	
	resulting_data = {solution_intersections, longest_customers_distances, smallest_customers_distances, depot_customers_distances, 
	all_average_gravity, all_routes_gravity_x, all_routes_gravity_y, all_central_line_distances, all_routes_radians, all_gravity_radians, 
	all_route_depot_distances, all_routes_customers, all_new_parameters, new_parameters_names, all_routes_distances, depot_customers_distances_divided, 
	longest_interior_distances, smallest_interior_distances, all_first_last_cust_demands, all_furthest_client_demand, all_degree_of_neighborhood, 
	all_routes_degree_of_neighborhood, all_routes_two_don, all_routes_average_tsp_area, all_routes_largest_tsp_area, all_routes_smallest_tsp_area};

	//cout << "Done with dataset_determine_vrp_parameters" << endl;
	return resulting_data;
}

DTOP_values dataset_determine_tsp_parameters(vector<shared_ptr <SampleInfo>> dataset_train, Instance& instance, int N=5){
	
	// Declaration of lists
	vector<vector<float>> all_routes_longest_customers_distances = vector<vector<float>>();
	vector<vector<float>> all_routes_smallest_customers_distances = vector<vector<float>>();
	vector<vector<float>> all_routes_largest_don = vector<vector<float>>();
	vector<vector<float>> all_routes_smallest_don = vector<vector<float>>();
	vector<vector<float>> depot_customers_distances = vector<vector<float>>();
	vector<float> all_routes_gravity_x = vector<float>();
	vector<float> all_routes_gravity_y = vector<float>();

	vector<float> all_route_depot_distances = vector<float>();

	vector<float> all_average_gravity = vector<float>();

	vector<vector<float>> all_routes_distances = vector<vector<float>>();
	vector<vector<float>> all_first_last_cust_demands = vector<vector<float>>();
	vector<vector<float>> all_routes_two_don = vector<vector<float>>();
	
	vector<float> all_routes_tsp_area = vector<float>();
	
	std::shared_ptr <Solution> solution;
	vector<vector<float>> dist_matrix;
	vector<vector<unsigned int>> D;
	vector<vector<int>> position_matrix;
	
	//cout << "Define los vectores" << endl;
	
	for(int i = 0; i < dataset_train.size(); i++){
		solution = dataset_train[i]->tsp_solution;
		dist_matrix = instance.dist_matrix;
		//position_matrix is used to get the degree of neighborhood
		D = instance.sorted_dist_matrix;
		position_matrix = vector<vector<int>>(D.size(), vector<int>(D.size(), 0));
		
		for(int u = 0; u < position_matrix.size(); u++){
			for(int v = 0; v < position_matrix.size(); v++){
				unsigned int position = D[u][v];
				position_matrix[u][position] = v;
			}
		}
		
		vector<vector<pair<std::shared_ptr <Customer>, std::shared_ptr <Customer>>>> vehicles_edge_list = vector<vector<pair<std::shared_ptr <Customer>, std::shared_ptr <Customer>>>>();
		int j = 0;
		int n_vehicles = solution->n_vehicles;
		
		vector<float> longest_dist_per_route(N);
		vector<float> smallest_dist_per_route(N);		
		vector<float> average_depot_distance(n_vehicles);
		
		vector<float> vector_x = vector<float>();
		vector<float> vector_y = vector<float>();
		
		vector<float> each_route_distance(n_vehicles);
		vector<float> average_first_last_cust_demand(n_vehicles);
		vector<float> each_route_two_don(n_vehicles);
		
		std::shared_ptr <Vehicle> vehicle;
		std::shared_ptr <Customer> curr_customer;
		std::shared_ptr <Customer> prev_customer;
		float tsp_area = 0;
		vector<float> largest_don(N);
		vector<float> smallest_don(N);
		vector<float> customer_to_grav;
		float depot_to_grav = 0;
		
		vehicle = solution->vehicles[0];
		vehicles_edge_list.push_back(vector<pair<std::shared_ptr <Customer>, std::shared_ptr <Customer>>>());
		
		curr_customer = nullptr;
		prev_customer = nullptr;
		vector<float> edge_features = vector<float>();
		vector<float> first_and_last_cust = vector<float>();
		vector<float> depot_features = vector<float>();
		
		float gravity_center_x = 0;
		float gravity_center_y = 0;
		float vehicle_cost = 0;
		
		vector<int> vehicle_customers_list = vector<int>();
		int customers_visited = 0;
		
		for(const auto& customer : vehicle->customers){
			curr_customer = customer;
			vehicle_customers_list.push_back(curr_customer->id);
			
			if(prev_customer != nullptr){
				vehicles_edge_list[j].push_back(make_pair(prev_customer, curr_customer));
				vehicle_cost += dist_matrix[prev_customer->id][curr_customer->id];
				
				vector_x.push_back(curr_customer->x);
				vector_y.push_back(curr_customer->y);
				
				// "determine gravity center according to half of the customers already routed"
				customers_visited++;
				
				if(customers_visited == (instance.n_cust/2)){
					gravity_center_x = curr_customer->x;
					gravity_center_y = curr_customer->y;
				}
				
				if(prev_customer->id == 0){
					float distance = dist_matrix[prev_customer->id][curr_customer->id];
					depot_features.push_back(distance);
					edge_features.push_back(distance);
					first_and_last_cust.push_back(curr_customer->demand);
				}
				
				if(curr_customer->id == 0){
					float distance = dist_matrix[prev_customer->id][curr_customer->id];
					depot_features.push_back(distance);
					edge_features.push_back(distance);
					first_and_last_cust.push_back(prev_customer->demand);
				}
			}
			
			if(prev_customer != nullptr && prev_customer->id != 0 && curr_customer->id != 0){
				float distance = dist_matrix[prev_customer->id][curr_customer->id];
				edge_features.push_back(distance);
			}
			
			prev_customer = curr_customer;
		}
		
		gravity_center_x = abs(gravity_center_x);
		gravity_center_y = abs(gravity_center_y);
		
		j++;
		average_depot_distance[0] = sum(depot_features);
		average_first_last_cust_demand[0] = sum(first_and_last_cust)/2;
		
		vector<float> clients_don(vehicle_customers_list.size() - 2);
		for(int i = 1; i < vehicle_customers_list.size() - 1; i++){
			clients_don[i-1] = (position_matrix[vehicle_customers_list[i-1]][vehicle_customers_list[i]] + position_matrix[vehicle_customers_list[i]][vehicle_customers_list[i+1]])/2;
		}
		
		sort(clients_don.begin(), clients_don.end());
		
		int count = 0;
		for(int i = clients_don.size()-N; i < clients_don.size(); i++){
			largest_don[count] = clients_don[i];
			count++;
		}
		
		for(int i = 0; i < N; i++){
			smallest_don[i] = clients_don[i];
		}
		
		each_route_two_don[0] = sum(clients_don)/clients_don.size();
		
		sort(edge_features.begin(), edge_features.end());
		count = 0;
		for(int i = edge_features.size()-N; i < edge_features.size(); i++){
			longest_dist_per_route[count] = edge_features[i];
			count++;
		}
		
		for(int i = 0; i < N; i++){
			smallest_dist_per_route[i] = edge_features[i];
		}
		
		int route_cardinality = vehicle->customers.size() - 1;
		
		gravity_center_x = gravity_center_x/route_cardinality;
		gravity_center_y = gravity_center_y/route_cardinality;
		each_route_distance[0] = vehicle_cost;
		
		customer_to_grav = {};
		
		for(const auto& customer : vehicle->customers){
			bg::model::point<float, 2, bg::cs::cartesian> grav_center_pos_p(gravity_center_x, gravity_center_y);
			if(customer->id == 0){
				bg::model::point<float, 2, bg::cs::cartesian> depot_pos_p(customer->x, customer->y);
				bg::model::point<float, 2, bg::cs::cartesian> grav_center_minus_depot = grav_center_pos_p;
				bg::subtract_point(grav_center_minus_depot, depot_pos_p);
				vector<float> grav_minus_depot_vector = {grav_center_minus_depot.get<0>(), grav_center_minus_depot.get<1>()};
				depot_to_grav = vectorNorm(grav_minus_depot_vector.begin(), grav_minus_depot_vector.end());
			}
			if(customer->id != 0){
				bg::model::point<float, 2, bg::cs::cartesian> customer_pos_p(customer->x, customer->y);
				bg::subtract_point(customer_pos_p, grav_center_pos_p);
				vector<float> cust_minus_grav_vector = {customer_pos_p.get<0>(), customer_pos_p.get<1>()};
				customer_to_grav.push_back(vectorNorm(cust_minus_grav_vector.begin(), cust_minus_grav_vector.end()));
			}
		}
		
		tsp_area = polygonArea(vector_x, vector_y, vector_x.size());
		
		
		
		all_routes_tsp_area.push_back(tsp_area);
		all_routes_largest_don.push_back(largest_don);
		all_routes_smallest_don.push_back(smallest_don);
		all_average_gravity.push_back(sum(customer_to_grav));
		
		all_route_depot_distances.push_back(depot_to_grav);
		
		all_routes_gravity_x.push_back(gravity_center_x);
		all_routes_gravity_y.push_back(gravity_center_y);
		all_routes_longest_customers_distances.push_back(longest_dist_per_route);
		all_routes_smallest_customers_distances.push_back(smallest_dist_per_route);
		depot_customers_distances.push_back(average_depot_distance);
		all_routes_distances.push_back(each_route_distance);
		all_first_last_cust_demands.push_back(average_first_last_cust_demand);
		
		all_routes_two_don.push_back(each_route_two_don);
	}
	//cout << "Sale del ciclo for" << endl;
	
	vector<vector<float>> route_parameters;
	vector<float> parameters_values;
	for(int i = 0; i < dataset_train.size(); i++){
		parameters_values = {};
		parameters_values.insert(parameters_values.end(), all_routes_longest_customers_distances[i].begin(), all_routes_longest_customers_distances[i].end());
		//cout << "size: " << all_routes_longest_customers_distances[i].size() << " longest_customers " << parameters_values.size() << endl;
		parameters_values.insert(parameters_values.end(), all_routes_smallest_customers_distances[i].begin(), all_routes_smallest_customers_distances[i].end());
		//cout << "size: " << all_routes_smallest_customers_distances[i].size() << " smallest customers " << parameters_values.size() << endl;
		parameters_values.insert(parameters_values.end(), depot_customers_distances[i].begin(), depot_customers_distances[i].end());
		//cout << "depot_customer_dsitances " << parameters_values.size() << endl;
		parameters_values.push_back(all_average_gravity[i]);
		//cout << "average gravity" << parameters_values.size() << endl;
		parameters_values.push_back(all_routes_gravity_x[i]);
		//cout << "gravity_x " << parameters_values.size() << endl;
		parameters_values.push_back(all_routes_gravity_y[i]);
		//cout << "gravity_y " << parameters_values.size() << endl;
		parameters_values.push_back(all_route_depot_distances[i]);
		//cout << "route_depot_distance " << parameters_values.size() << endl;
		parameters_values.insert(parameters_values.end(), all_routes_distances[i].begin(), all_routes_distances[i].end());
		//cout << "route_distances " << parameters_values.size() << endl;
		parameters_values.insert(parameters_values.end(), all_first_last_cust_demands[i].begin(), all_first_last_cust_demands[i].end());
		//cout << "first and last " << parameters_values.size() << endl;
		parameters_values.insert(parameters_values.end(), all_routes_two_don[i].begin(), all_routes_two_don[i].end());
		//cout << "two_don " << parameters_values.size() << endl;
		parameters_values.insert(parameters_values.end(), all_routes_largest_don[i].begin(), all_routes_largest_don[i].end());
		//cout << "largest_don " << parameters_values.size() << endl;
		parameters_values.push_back(all_routes_tsp_area[i]);
		route_parameters.push_back(parameters_values);
	}
	
	//cout << "Pasa de definir route_parameters" << endl;
	
	vector<string> parameters_and_names = {"lcd_1", "lcd_2", "lcd_3", "lcd_4", "lcd_5",
											"scd_1", "scd_2", "scd_3", "scd_4", "scd_5",
											"depot_custom_dist", "average_gravity", "gravity_x", "gravity_y", 
											"depot_gravity_dist", "initial_distance", "first_last_demand", "sum_DON",
											"ldon_1", "ldon_2", "ldon_3", "ldon_4", "ldon_5", "tsp_area"};
					
	return {route_parameters, parameters_and_names};
}

DTOP_values dataset_to_vrp_parameters(vector<shared_ptr <SampleInfo>> dataset_train, Instance& instance){
	DDOP_values all_parameters = dataset_determine_vrp_parameters(dataset_train, instance);
	//cout << "Obtiene all_parameters" << endl;
	vector<float> s1_solution_intersections(all_parameters.solution_intersections.size());
	vector<float> s2_longest_distance_route(all_parameters.solution_intersections.size());
	vector<float> s3_average_depot_distance_route(all_parameters.solution_intersections.size());
	vector<float> s4_average_gravity_distance = all_parameters.all_average_gravity;
	vector<vector<vector<float>>> s5 = all_parameters.all_central_line_distances;
	
	/*for(int i = 0; i < all_parameters.all_central_line_distances.size(); i++){
		for(int j = 0; j < all_parameters.all_central_line_distances[i].size(); j++){
			for(int k = 0; k < all_parameters.all_central_line_distances[i][j].size(); k++){
				cout << all_parameters.all_central_line_distances[i][j][k] << " ";
			}
			cout << endl;
		}
	}*/
	vector<float> s5_average_width_route(all_parameters.solution_intersections.size());
	vector<vector<vector<float>>> s6 = all_parameters.all_routes_radians;
	vector<float> s6_average_span_route(all_parameters.solution_intersections.size());
	vector<float> s7_average_width_compactness(all_parameters.solution_intersections.size());
	vector<vector<vector<float>>> s8 = all_parameters.all_gravity_radians;
	vector<float> s8_average_span_compactness(all_parameters.solution_intersections.size());
	vector<vector<vector<float>>> s9 = all_parameters.all_route_depot_distances;
	vector<float> s9_average_depth_route(all_parameters.solution_intersections.size());
	vector<vector<int>> s10 = all_parameters.all_routes_customers;
	vector<float> s10_std_dev_customers_route(all_parameters.solution_intersections.size());
	
	vector<float> routes_min(all_parameters.solution_intersections.size());
	vector<float> routes_max(all_parameters.solution_intersections.size());
	vector<float> routes_mean(all_parameters.solution_intersections.size());
	vector<float> routes_std(all_parameters.solution_intersections.size());
	
	vector<float> initial_distances(all_parameters.solution_intersections.size());
	
	for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
		s1_solution_intersections[i] = all_parameters.solution_intersections[i]/instance.n_cust;
		s2_longest_distance_route[i] = sum(all_parameters.longest_customers_distances[i])/all_parameters.longest_customers_distances[i].size();
		s3_average_depot_distance_route[i] = sum(all_parameters.depot_customers_distances[i])/(2*all_parameters.depot_customers_distances[i].size());
		routes_mean[i] = sum(all_parameters.all_routes_distances[i])/all_parameters.all_routes_distances[i].size();
		routes_min[i] = *min_element(all_parameters.all_routes_distances[i].begin(), all_parameters.all_routes_distances[i].end());
		routes_max[i] = *max_element(all_parameters.all_routes_distances[i].begin(), all_parameters.all_routes_distances[i].end());
	}
	
	try{
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			vector<float> v_desv = all_parameters.all_routes_distances[i];
			routes_std[i] = desv_est(v_desv);
		}
	}catch(...){
		routes_std = vector<float>(all_parameters.solution_intersections.size());
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			routes_std[i] = 0;
		}
	}
	
	vector<float> s3_mod(all_parameters.solution_intersections.size());
	vector<float> s2_longest_mod(all_parameters.solution_intersections.size());
	vector<float> s2_shortest_mod(all_parameters.solution_intersections.size());
	vector<float> first_last_cust_demands(all_parameters.solution_intersections.size());
	vector<float> furthest_customers_demand(all_parameters.solution_intersections.size());
	
	for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
		s3_mod[i] = sum(all_parameters.depot_customers_distances_divided[i])/(2*all_parameters.depot_customers_distances_divided[i].size());
		s2_longest_mod[i] = sum(all_parameters.longest_interior_distances[i])/all_parameters.longest_interior_distances[i].size();
		s2_shortest_mod[i] = sum(all_parameters.smallest_interior_distances[i])/all_parameters.smallest_interior_distances[i].size();
		first_last_cust_demands[i] = sum(all_parameters.all_first_last_cust_demands[i])/all_parameters.all_first_last_cust_demands[i].size();
		furthest_customers_demand[i] = sum(all_parameters.all_furthest_client_demand[i])/all_parameters.all_furthest_client_demand[i].size();
	}
	
	vector<float> furthest_customers_demand_std(all_parameters.solution_intersections.size());
	try{
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			vector<float> v_desv = all_parameters.all_furthest_client_demand[i];
			furthest_customers_demand_std[i] = desv_est(v_desv);
		}
	}catch(...){
		furthest_customers_demand_std = vector<float>(all_parameters.solution_intersections.size());
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			furthest_customers_demand_std[i] = 0;
		}
	}
	
	vector<vector<float>> degree_of_neighborhood = all_parameters.all_degree_of_neighborhood;
	
	vector<string> degree_names(degree_of_neighborhood[0].size()); // duda acá
	for(int i = 0; i < degree_of_neighborhood[0].size(); i++){
		degree_names[i] = "DON_" + i;
	} // Duda aca
	
	vector<float> routes_degree_of_neighborhood(all_parameters.solution_intersections.size());
	vector<float> routes_degree_of_neighborhood_std(all_parameters.solution_intersections.size());
	vector<float> routes_two_degree_of_neighborhood(all_parameters.solution_intersections.size());
	vector<float> routes_two_degree_of_neighborhood_std(all_parameters.solution_intersections.size());


	for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
		routes_degree_of_neighborhood[i] = sum(all_parameters.all_routes_degree_of_neighborhood[i])/all_parameters.all_routes_degree_of_neighborhood[i].size();
		routes_two_degree_of_neighborhood[i] = sum(all_parameters.all_routes_two_don[i])/all_parameters.all_routes_two_don[i].size();

	}

	try{
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			vector<float> v_desv = all_parameters.all_routes_degree_of_neighborhood[i];
			routes_degree_of_neighborhood_std[i] = desv_est(v_desv);
		}
	}catch(...){
		routes_degree_of_neighborhood_std = vector<float>(all_parameters.solution_intersections.size());
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			routes_degree_of_neighborhood_std[i] = 0;
		}
	}
	
	try{
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			vector<float> v_desv = all_parameters.all_routes_two_don[i];
			routes_two_degree_of_neighborhood_std[i] = desv_est(v_desv);
		}
	}catch(...){
		routes_two_degree_of_neighborhood_std = vector<float>(all_parameters.solution_intersections.size());
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			routes_two_degree_of_neighborhood_std[i] = 0;
		}
	}
	
	// Average Smallest distance
	vector<float> s2_smallest_distance_route(all_parameters.solution_intersections.size());
	// Min largest distance
	vector<float> s21(all_parameters.solution_intersections.size());
	// Max largest distance
	vector<float> s22(all_parameters.solution_intersections.size());
	// Min smallest distance
	vector<float> s23(all_parameters.solution_intersections.size());
	// Max smallest distance
	vector<float> s24(all_parameters.solution_intersections.size());
	// Standard deviation largest distance
	vector<float> std_dev_largest_distance(all_parameters.solution_intersections.size());
	// Standard deviation smallest distance
	vector<float> std_dev_smallest_distance(all_parameters.solution_intersections.size());
	// Min distance to depot
	vector<float> s31(all_parameters.solution_intersections.size());
	// Max distance to depot
	vector<float> s32(all_parameters.solution_intersections.size());
	
	for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
		s2_smallest_distance_route[i] = sum(all_parameters.smallest_customers_distances[i])/all_parameters.smallest_customers_distances[i].size();
		s21[i] = *min_element(all_parameters.longest_customers_distances[i].begin(), all_parameters.longest_customers_distances[i].end());
		s22[i] = *max_element(all_parameters.longest_customers_distances[i].begin(), all_parameters.longest_customers_distances[i].end());
		s23[i] = *min_element(all_parameters.smallest_customers_distances[i].begin(), all_parameters.smallest_customers_distances[i].end());
		s24[i] = *max_element(all_parameters.smallest_customers_distances[i].begin(), all_parameters.smallest_customers_distances[i].end());
		s31[i] = (*min_element(all_parameters.depot_customers_distances[i].begin(), all_parameters.depot_customers_distances[i].end()))/2;
		s32[i] = (*max_element(all_parameters.depot_customers_distances[i].begin(), all_parameters.depot_customers_distances[i].end()))/2;
	}
	
	try{
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			vector<float> v_desv = all_parameters.longest_customers_distances[i];
			std_dev_largest_distance[i] = desv_est(v_desv);
		}
	}catch(...){
		std_dev_largest_distance = vector<float>(all_parameters.solution_intersections.size());
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			std_dev_largest_distance[i] = 0;
		}
	}
	
	try{
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			vector<float> v_desv = all_parameters.smallest_customers_distances[i];
			std_dev_smallest_distance[i] = desv_est(v_desv);
		}
	}catch(...){
		std_dev_smallest_distance = vector<float>(all_parameters.solution_intersections.size());
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			std_dev_smallest_distance[i] = 0;
		}
	}
	
	// Std width per route
	vector<float> s5_std(all_parameters.solution_intersections.size());
	// Min width per route
	vector<float> s51(all_parameters.solution_intersections.size());
	// Max width per route
	vector<float> s52(all_parameters.solution_intersections.size());
	// Std span per route
	vector<float> s6_std(all_parameters.solution_intersections.size());
	// Min span per route
	vector<float> s61(all_parameters.solution_intersections.size());
	// Max span per route
	vector<float> s62(all_parameters.solution_intersections.size());
	// Std width compactness per route
	vector<float> s7_std(all_parameters.solution_intersections.size());
	// Min width compactness per route
	vector<float> s71(all_parameters.solution_intersections.size());
	// Max width compactness per route
	vector<float> s72(all_parameters.solution_intersections.size());
	// Std radian compactness per route
	vector<float> s8_std(all_parameters.solution_intersections.size());
	// Min radian compactness per route
	vector<float> s81(all_parameters.solution_intersections.size());
	// Max radian compactness per route
	vector<float> s82(all_parameters.solution_intersections.size());
	// Std depth per route
	vector<float> s9_std(all_parameters.solution_intersections.size());
	// Min depth per route
	vector<float> s91(all_parameters.solution_intersections.size());
	// Max depth per route
	vector<float> s92(all_parameters.solution_intersections.size());
	
	vector<vector<int>> client_pairs = all_parameters.all_new_parameters;
	vector<string> client_pairs_names = all_parameters.new_parameters_names;
	//cout << "Llega hasta aca" << endl;
	for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
		vector<float> average_width_route;
		vector<float> average_span_route(s6[i].size());
		vector<float> average_width_compactness;
		vector<float> average_span_compactness;
		vector<float> average_depth_route;
		vector<double> std_dev_route;
		
		for(int j = 0; j < s5[i].size(); j++){
			average_width_route.push_back(*max_element(s5[i][j].begin(), s5[i][j].end()) - *min_element(s5[i][j].begin(), s5[i][j].end()));
			vector<float> width_compactness;
			for(int z = 0; z < s5[i][j].size(); z++){
				width_compactness.push_back(abs(s5[i][j][z]));
			}
			average_width_compactness.push_back(sum(width_compactness));
		}
		
		s5_average_width_route[i] = sum(average_width_route)/s5[i].size();
		s51[i] = *min_element(average_width_route.begin(), average_width_route.end());
		s52[i] = *max_element(average_width_route.begin(), average_width_route.end());
		try{
			s5_std[i] = desv_est(average_width_route);
		}catch(...){
			s5_std[i] = 0;
		}
		
		s7_average_width_compactness[i] = sum(average_width_compactness)/instance.n_cust;
		//cout << sum(average_width_compactness) << endl;
		s71[i] = *min_element(average_width_compactness.begin(), average_width_compactness.end());
		s72[i] = *max_element(average_width_compactness.begin(), average_width_compactness.end());
		try{
			s7_std[i] = desv_est(average_width_compactness);
		}catch(...){
			s7_std[i] = 0;
		}
		
		for(int w = 0; w < s6[i].size(); w++){
			try{
				if(s6[i][w].size() > 0){
					average_span_route[w] = *max_element(s6[i][w].begin(), s6[i][w].end());
				}
			}catch(...){
				average_span_route[w] = 0.0;
			}
		}
		
		s6_average_span_route[i] = sum(average_span_route)/s6[i].size();
		s61[i] = *min_element(average_span_route.begin(), average_span_route.end());
		s62[i] = *max_element(average_span_route.begin(), average_span_route.end());
		try{
			s6_std[i] = desv_est(average_span_route);
		}catch(...){
			s6_std[i] = 0;
		}
		
		for(int k = 0; k < s8[i].size(); k++){
			average_span_compactness.push_back(sum(s8[i][k]));
		}
		
		s8_average_span_compactness[i] = sum(average_span_compactness)/instance.n_cust;
		s81[i] = *min_element(average_span_compactness.begin(), average_span_compactness.end());
		s82[i] = *max_element(average_span_compactness.begin(), average_span_compactness.end());
		try{
			s8_std[i] = desv_est(average_span_compactness);
		}catch(...){
			s8_std[i] = 0;
		}
		
		for(int q = 0; q < s9[i].size(); q++){
			try{
				average_depth_route.push_back(*max_element(s9[i][q].begin(), s9[i][q].end()));
			}catch(...){
				average_depth_route.push_back(0.0);
			}
		}
		
		s9_average_depth_route[i] = sum(average_depth_route)/s9[i].size();
		s91[i] = *min_element(average_depth_route.begin(), average_depth_route.end());
		s92[i] = *max_element(average_depth_route.begin(), average_depth_route.end());
		try{
			s9_std[i] = desv_est(average_depth_route);
		}catch(...){
			s9_std[i] = 0;
		}
		float clients_number = instance.n_cust;
		for(int h = 0; h < s10[i].size(); h++){
			std_dev_route.push_back(pow(s10[i][h] - (clients_number/s10[i].size()), 2));
		}
		s10_std_dev_customers_route[i] = sqrt(sum(std_dev_route)/s10[i].size());
		initial_distances[i] = dataset_train[i]->initial_solution->cost;
	}
	
	vector<vector<float>> features;
	vector<string> features_names;

    vector<float> routes_average_tsp_area(all_parameters.all_routes_average_tsp_area.size());
    vector<float> routes_largest_tsp_area(all_parameters.all_routes_average_tsp_area.size());
    vector<float> routes_smallest_tsp_area(all_parameters.all_routes_average_tsp_area.size());

    /*
    cout << "Arrived here 2" << endl;
    for(int i = 0; i < dataset_train.size(); i++){
        routes_average_tsp_area[i] = all_parameters.all_routes_average_tsp_area[i];
        routes_largest_tsp_area[i] = all_parameters.all_routes_largest_tsp_area[i];
        routes_smallest_tsp_area[i] = all_parameters.all_routes_smallest_tsp_area[i];
    }
    cout << "Arrived here 3" << endl;
     */

    //cout << " all_parameters.all_routes_average_tsp_area.size() " << all_parameters.all_routes_average_tsp_area.size() << endl;
    //cout << " all_parameters.all_routes_largest_tsp_area.size() " << all_parameters.all_routes_largest_tsp_area.size() << endl;
    //cout << " all_parameters.all_routes_smallest_tsp_area.size() " << all_parameters.all_routes_smallest_tsp_area.size() << endl;
	
	for(int i = 0; i < dataset_train.size(); i++){
		features.push_back({s1_solution_intersections[i], s2_longest_distance_route[i], s3_average_depot_distance_route[i], s4_average_gravity_distance[i], 
							s5_average_width_route[i], s6_average_span_route[i], s7_average_width_compactness[i], s8_average_span_compactness[i], 
							s9_average_depth_route[i], s10_std_dev_customers_route[i], s21[i], s22[i], s23[i], s24[i],
							s31[i], s32[i], s51[i], s52[i], s61[i], s62[i], s71[i], s72[i], 
							s81[i], s82[i], s91[i], s92[i], 
							std_dev_largest_distance[i], std_dev_smallest_distance[i],
							s5_std[i], s6_std[i], s7_std[i], s8_std[i], s9_std[i], 
							routes_min[i], routes_max[i], routes_mean[i], routes_std[i], 
							s3_mod[i], s2_longest_mod[i], s2_shortest_mod[i], first_last_cust_demands[i],
							furthest_customers_demand[i], furthest_customers_demand_std[i], routes_degree_of_neighborhood[i], 
							routes_degree_of_neighborhood_std[i], routes_two_degree_of_neighborhood[i], routes_two_degree_of_neighborhood_std[i], initial_distances[i],
                          all_parameters.all_routes_average_tsp_area[i], all_parameters.all_routes_largest_tsp_area[i], all_parameters.all_routes_smallest_tsp_area[i]});
	}

	
	features_names = {"S1", "S2", "S3", "S4", "S5", "S6", "S7",
					 "S8", "S9", "S10", "S21", "S22", "S23", "S24", 
					 "S31", "S32", "S51",
					 "S52", "S61", "S62", "S71", "S72", "S81", "S82", "S91", "S92",
					 "S2_L_STD", "S2_S_STD", "S5_STD", "S6_STD", 
					 "S7_STD", "S8_STD", "S9_STD", 
					 "RL_MIN", "RL_MAX", "RL_MEAN", "RL_STD",
					 "S3_*", "S2_L_MOD", "S2_S_MOD", "F_L_CD", 
					 "FC_DEM_MEAN", "FC_DEM_STD", "R_DON", "R_DON_STD",
					 "R_2DON", "R_2DON_STD", "INITIAL_DISTANCE", "AVERAGE_ROUTE_AREA", "LARGEST_ROUTE_AREA", "SMALLEST_ROUTE_AREA"};
					 
	DTOP_values parameters_and_names = {features, features_names};
	return parameters_and_names;
}

DTOP_values dataset_determine_tsp_vrp_parameters(vector<shared_ptr <SampleInfo>>dataset_train, Instance& instance, int N=5){
	DTOP_values tsp_values = dataset_determine_tsp_parameters(dataset_train, instance, N);
	DTOP_values vrp_values = dataset_to_vrp_parameters(dataset_train, instance);
	
	for(int i = 0; i < tsp_values.features.size(); i++){
		tsp_values.features[i].insert(tsp_values.features[i].end(), vrp_values.features[i].begin(), vrp_values.features[i].end());
	}
	
	tsp_values.features_names.insert(tsp_values.features_names.end(), vrp_values.features_names.begin(), vrp_values.features_names.end());
	
	return tsp_values;
}

DTOP_values dataset_determine_tsp_client_pairs_parameters(vector<shared_ptr <SampleInfo>> dataset_train, Instance& instance, int N=5){
	DTOP_values tsp_values = dataset_determine_tsp_parameters(dataset_train, instance, N);
	DTOP_values client_pairs_values = dataset_determine_client_pairs_parameters(dataset_train, instance);
	
	for(int i = 0; i < tsp_values.features.size(); i++){
		tsp_values.features[i].insert(tsp_values.features[i].end(), client_pairs_values.features[i].begin(), client_pairs_values.features[i].end());
	}
	
	tsp_values.features_names.insert(tsp_values.features_names.end(), client_pairs_values.features_names.begin(), client_pairs_values.features_names.end());
	
	return tsp_values;
}

DTOP_values dataset_determine_vrp_client_pairs_parameters(vector<shared_ptr <SampleInfo>> dataset_train, Instance& instance){
	DTOP_values vrp_values = dataset_to_vrp_parameters(dataset_train, instance);
	DTOP_values client_pairs_values = dataset_determine_client_pairs_parameters(dataset_train, instance);
	
	for(int i = 0; i < vrp_values.features.size(); i++){
		vrp_values.features[i].insert(vrp_values.features[i].end(), client_pairs_values.features[i].begin(), client_pairs_values.features[i].end());
	}
	
	vrp_values.features_names.insert(vrp_values.features_names.end(), client_pairs_values.features_names.begin(), client_pairs_values.features_names.end());
	
	return vrp_values;
}

DTOP_values dataset_determine_tsp_vrp_client_pairs_parameters(vector<shared_ptr <SampleInfo>> dataset_train, Instance& instance, int N=5){
	DTOP_values tsp_values = dataset_determine_tsp_parameters(dataset_train, instance, N);
	DTOP_values vrp_values = dataset_to_vrp_parameters(dataset_train, instance);
	DTOP_values client_pairs_values = dataset_determine_client_pairs_parameters(dataset_train, instance);
	
	for(int i = 0; i < tsp_values.features.size(); i++){
		tsp_values.features[i].insert(tsp_values.features[i].end(), vrp_values.features[i].begin(), vrp_values.features[i].end());
		tsp_values.features[i].insert(tsp_values.features[i].end(), client_pairs_values.features[i].begin(), client_pairs_values.features[i].end());
	}
	
	tsp_values.features_names.insert(tsp_values.features_names.end(), vrp_values.features_names.begin(), vrp_values.features_names.end());
	tsp_values.features_names.insert(tsp_values.features_names.end(), client_pairs_values.features_names.begin(), client_pairs_values.features_names.end());
	
	return tsp_values;
}

DDOP_values solution_determine_other_parameters(std::shared_ptr <Solution> initial_solution, Instance& instance){
	DDOP_values resulting_data;
	
	// Declaration of lists
	vector<int> solution_intersections = vector<int>();
	vector<vector<float>> longest_customers_distances = vector<vector<float>>();
	vector<vector<float>> smallest_customers_distances = vector<vector<float>>();
	vector<vector<float>> longest_interior_distances = vector<vector<float>>();
	vector<vector<float>> smallest_interior_distances = vector<vector<float>>();
	vector<vector<float>> depot_customers_distances = vector<vector<float>>();
	vector<float> all_average_gravity = vector<float>();
	
	vector<vector<float>> all_routes_gravity_x = vector<vector<float>>();
	vector<vector<float>> all_routes_gravity_y = vector<vector<float>>();
	vector<vector<vector<float>>> all_central_line_distances = vector<vector<vector<float>>>();
	vector<vector<vector<float>>> all_routes_radians = vector<vector<vector<float>>>();
	vector<vector<vector<float>>> all_gravity_radians = vector<vector<vector<float>>>();
	vector<vector<vector<float>>> all_route_depot_distances = vector<vector<vector<float>>>();
	vector<vector<int>> all_routes_customers = vector<vector<int>>();
	vector<vector<int>> all_new_parameters = vector<vector<int>>();
	vector<string> new_parameters_names = vector<string>();
	vector<vector<float>> all_routes_distances = vector<vector<float>>();
	vector<vector<float>> depot_customers_distances_divided = vector<vector<float>>();
	vector<vector<float>> all_first_last_cust_demands = vector<vector<float>>();
	vector<vector<float>> all_furthest_client_demand = vector<vector<float>>();
	vector<vector<float>> all_degree_of_neighborhood = vector<vector<float>>(); // Revisar
	vector<vector<float>> all_routes_degree_of_neighborhood = vector<vector<float>>();
	vector<vector<float>> all_routes_two_don = vector<vector<float>>();
	
	std::shared_ptr <Solution> solution = initial_solution;
	vector<vector<float>> dist_matrix = instance.dist_matrix;
	vector<vector<unsigned int>> D = instance.sorted_dist_matrix;
	//position_matrix is used to get the degree of neighborhood
	vector<vector<int>> position_matrix = vector<vector<int>>(D.size(), vector<int>(D.size(), 0));
	vector<float> degree_of_neighborhood = vector<float>(D.size(), 0);
	
	for(int u = 0; u < position_matrix.size(); u++){
		for(int v = 0; v < position_matrix.size(); v++){
			unsigned int position = D[u][v];
			position_matrix[u][position] = v;
		}
	}
	
	vector<vector<pair<std::shared_ptr <Customer>, std::shared_ptr <Customer>>>> vehicles_edge_list = vector<vector<pair<std::shared_ptr <Customer>, std::shared_ptr <Customer>>>>();
	int j = 0;
	int n_vehicles = solution->n_vehicles;
	
	vector<float> longest_dist_per_route(n_vehicles);
	vector<float> smallest_dist_per_route(n_vehicles);
	vector<float> longest_interior_dist_per_route(n_vehicles);
	vector<float> smallest_interior_dist_per_route(n_vehicles);
	
	vector<float> average_depot_distance(n_vehicles);
	vector<float> average_depot_distance_divided(n_vehicles);
	vector<float> routes_gravity_x(n_vehicles);
	vector<float> routes_gravity_y(n_vehicles);
	
	vector<vector<float>> routes_central_line_distance(n_vehicles);
	vector<vector<float>> each_route_radians(n_vehicles);
	vector<vector<float>> each_gravity_radians(n_vehicles);
	vector<vector<float>> each_route_depot_distance(n_vehicles);
	
	vector<int> route_customers(n_vehicles);
	map<string, int> new_parameters;
	
	for(int x = 1; x < instance.n_cust; x++){
		for(int y = x; y < instance.n_cust; y++){
			if(x != y){
				string key = x + "-" + y;
				new_parameters[key] = 0;
			}
		}
	}
	
	vector<float> each_route_distance(n_vehicles);
	vector<float> average_first_last_cust_demand(n_vehicles);
	vector<float> demand_furthest(n_vehicles);
	vector<float> each_route_don(n_vehicles);
	vector<float> each_route_two_don(n_vehicles);
	
	std::shared_ptr <Vehicle> vehicle;
	std::shared_ptr <Customer> curr_customer;
	std::shared_ptr <Customer> prev_customer;
	
	for(int n = 0; n < n_vehicles; n++){
		vehicle = solution->vehicles[n];
		if(vehicle->customers.size() > 2){
			vehicles_edge_list.push_back(vector<pair<std::shared_ptr <Customer>, std::shared_ptr <Customer>>>());
			
			curr_customer = nullptr;
			prev_customer = nullptr;
			vector<float> edge_features = vector<float>();
			vector<float> first_and_last_cust = vector<float>();
			vector<float> interior_edge_features = vector<float>();
			vector<float> depot_features = vector<float>();
			
			float furthest_distance = 0;
			float furthest_client_demand = 0;
			
			float gravity_center_x = 0;
			float gravity_center_y = 0;
			float vehicle_cost = 0;
			
			vector<int> vehicle_customers_list = vector<int>();
			float curr_distance;
			
			for(const auto& customer : vehicle->customers){
				curr_customer = customer;
				curr_distance = dist_matrix[0][curr_customer->id];
				vehicle_customers_list.push_back(curr_customer->id);
				
				if(prev_customer != nullptr){
					gravity_center_x += curr_customer->x;
					gravity_center_y += curr_customer->y;
					vehicles_edge_list[j].push_back(make_pair(prev_customer, curr_customer));
					vehicle_cost += dist_matrix[prev_customer->id][curr_customer->id];
					
					if(prev_customer->id == 0){
						float distance = dist_matrix[prev_customer->id][curr_customer->id];
						depot_features.push_back(distance);
						edge_features.push_back(distance);
						first_and_last_cust.push_back(curr_customer->demand);
					}
					
					if(curr_customer->id == 0){
						float distance = dist_matrix[prev_customer->id][curr_customer->id];
						depot_features.push_back(distance);
						edge_features.push_back(distance);
						first_and_last_cust.push_back(prev_customer->demand);
					}
					
					if(curr_distance > furthest_distance){
						furthest_distance = curr_distance;
						furthest_client_demand = curr_customer->demand;
					}
				}
				
				if(prev_customer != nullptr && prev_customer->id != 0 && curr_customer->id != 0){
					float distance = dist_matrix[prev_customer->id][curr_customer->id];
					interior_edge_features.push_back(distance);
					edge_features.push_back(distance);
				}
				
				prev_customer = curr_customer;
			}
			
			j++;
			average_depot_distance_divided[n] = sum(depot_features)/vehicle_cost;
			average_depot_distance[n] = sum(depot_features);
			average_first_last_cust_demand[n] = sum(first_and_last_cust)/2;
			demand_furthest[n] = furthest_client_demand;
			
			vector<float> clients_don(vehicle_customers_list.size() - 2);
			for(int i = 1; i < vehicle_customers_list.size() - 1; i++){
				clients_don[i-1] = (position_matrix[vehicle_customers_list[i-1]][vehicle_customers_list[i]] + position_matrix[vehicle_customers_list[i]][vehicle_customers_list[i+1]])/2;
			}
			
			each_route_two_don[n] = sum(clients_don)/clients_don.size();
			
			if(interior_edge_features.size() > 0){
				longest_dist_per_route[n] = (*max_element(edge_features.begin(), edge_features.end()))/vehicle_cost;
				smallest_dist_per_route[n] = (*min_element(edge_features.begin(), edge_features.end()))/vehicle_cost;
				longest_interior_dist_per_route[n] = (*max_element(interior_edge_features.begin(), interior_edge_features.end()))/vehicle_cost;
				smallest_interior_dist_per_route[n] = (*min_element(interior_edge_features.begin(), interior_edge_features.end()))/vehicle_cost;
			}
			
			int route_cardinality = vehicle->customers.size() - 1;
			
			gravity_center_x = gravity_center_x/route_cardinality;
			gravity_center_y = gravity_center_y/route_cardinality;
			routes_gravity_x[n] = gravity_center_x;
			routes_gravity_y[n] = gravity_center_y;
			each_route_distance[n] = vehicle_cost;
			
			vector<float> central_line_distance = vector<float>();
			vector<float> route_radians = vector<float>();
			vector<float> gravity_radians = vector<float>();
			vector<float> each_depot_distance = vector<float>();
			vector<float> route_degree_of_neighborhood = vector<float>();
			
			for(const auto& customer : vehicle->customers){
				int depot_id;
				Point depot_pos;
				Point grav_center_pos;
				
				
				if(customer->id == 0){
					depot_pos = customer->pos;
					depot_id = customer->id;
					grav_center_pos = {gravity_center_x, gravity_center_y};
				}
				
				bg::model::point<double, 2, bg::cs::cartesian> vector_1;
				bg::model::point<double, 2, bg::cs::cartesian> vector_grav;
				float grav_dot_product;
				float angle_grav;
				
				if(customer->id != 0){
					bg::model::point<double, 2, bg::cs::cartesian> depot_pos_p(depot_pos.x, depot_pos.y);
					bg::model::point<double, 2, bg::cs::cartesian> grav_center_pos_p(grav_center_pos.x, grav_center_pos.y);
					bg::model::point<double, 2, bg::cs::cartesian> customer_pos_p(customer->pos.x, customer->pos.y);
					
					bg::model::point<double, 2, bg::cs::cartesian> grav_minus_depot = grav_center_pos_p;
					bg::model::point<double, 2, bg::cs::cartesian> customer_minus_depot = customer_pos_p;
					
					bg::subtract_point(grav_minus_depot, depot_pos_p);
					bg::subtract_point(customer_minus_depot, depot_pos_p);
					
					bg::model::point<double, 2, bg::cs::cartesian> cross_product_result;
					cross_product_result = bg::cross_product(grav_minus_depot, customer_minus_depot);
					
					vector<double> grav_minus_depot_vector = {grav_minus_depot.get<0>(), grav_minus_depot.get<1>()};
					vector<double> customer_minus_depot_vector = {customer_minus_depot.get<0>(), customer_minus_depot.get<1>()};
					
					central_line_distance.push_back(cross_product_result.get<0>()/vectorNorm(grav_minus_depot_vector.begin(), grav_minus_depot_vector.end()));
					
					vector_1 = customer_minus_depot;
					bg::divide_value(vector_1, vectorNorm(customer_minus_depot_vector.begin(), customer_minus_depot_vector.end()));
					
					vector_grav = grav_minus_depot;
					bg::divide_value(vector_grav, vectorNorm(grav_minus_depot_vector.begin(), grav_minus_depot_vector.end()));
					
					grav_dot_product = bg::dot_product(vector_1, vector_grav);
					angle_grav = acos(grav_dot_product);
					gravity_radians.push_back(angle_grav);
					
					float distance = dist_matrix[depot_id][customer->id];
					each_depot_distance.push_back(distance);
				}
				
				vector<int> curr_neigh_degree = vector<int>();
				for(const auto& customer_2 : vehicle->customers){
					if(customer->id != customer_2->id && customer->id != 0){
						curr_neigh_degree.push_back(position_matrix[customer->id][customer_2->id]);
					}
					
					if(customer->id != customer_2->id && customer_2->id != 0 && customer->id != 0){
						bg::model::point<double, 2, bg::cs::cartesian> depot_pos_p(depot_pos.x, depot_pos.y);
						bg::model::point<double, 2, bg::cs::cartesian> customer_2_minus_depot(customer_2->pos.x, customer_2->pos.y);
						bg::subtract_point(customer_2_minus_depot, depot_pos_p);
						
						vector<double> customer_2_minus_depot_vector = {customer_2_minus_depot.get<0>(), customer_2_minus_depot.get<1>()};
						
						bg::model::point<double, 2, bg::cs::cartesian> vector_2 = customer_2_minus_depot;
						bg::divide_value(vector_2, vectorNorm(customer_2_minus_depot_vector.begin(), customer_2_minus_depot_vector.end()));
						
						float dot_product = bg::dot_product(vector_1, vector_2);
						float angle = acos(dot_product);
						route_radians.push_back(angle);
						string key = customer->id + "-" + customer_2->id;
						if(new_parameters.count(key) == 1){
							new_parameters[key] = 1;
						}
					}
				}
				
				degree_of_neighborhood[customer->id] = sum_int(curr_neigh_degree)/(vehicle->customers.size() - 1);
				route_degree_of_neighborhood.push_back(degree_of_neighborhood[customer->id]);
				//Sum of neighborhood degree, divided by number of customers in route - 1 because depot visited twice on each route
			}
			
			each_route_radians[n] = route_radians;
			routes_central_line_distance[n] = central_line_distance;
			each_gravity_radians[n] = gravity_radians;
			each_route_depot_distance[n] = each_depot_distance;
			route_customers[n] = vehicle->customers.size() - 2;
			each_route_don[n] = sum(route_degree_of_neighborhood)/route_degree_of_neighborhood.size();
		}
	}
	
	vector<float> gravity_distances = vector<float>();
	for(int i = 0; i < routes_gravity_x.size(); i++){
		for(int j = 0; j < routes_gravity_x.size(); j++){
			if(i != j){
				bg::model::point<double, 2, bg::cs::cartesian> a(routes_gravity_x[i], routes_gravity_y[i]);
				bg::model::point<double, 2, bg::cs::cartesian> b(routes_gravity_x[j], routes_gravity_y[j]);
				bg::subtract_point(a, b);
				
				vector<double> a_vector = {a.get<0>(), a.get<1>()};
				gravity_distances.push_back(vectorNorm(a_vector.begin(), a_vector.end())); 
			}
		}
	}
	
	all_average_gravity.push_back(sum(gravity_distances)/solution->vehicles.size());
	
	all_routes_customers.push_back(route_customers);
	all_route_depot_distances.push_back(each_route_depot_distance);
	all_gravity_radians.push_back(each_gravity_radians);
	all_routes_radians.push_back(each_route_radians);
	all_central_line_distances.push_back(routes_central_line_distance);
	all_routes_gravity_x.push_back(routes_gravity_x);
	all_routes_gravity_y.push_back(routes_gravity_y);
	longest_customers_distances.push_back(longest_dist_per_route);
	smallest_customers_distances.push_back(smallest_dist_per_route);
	longest_interior_distances.push_back(longest_interior_dist_per_route);
	smallest_interior_distances.push_back(smallest_interior_dist_per_route);
	depot_customers_distances.push_back(average_depot_distance);
	depot_customers_distances_divided.push_back(average_depot_distance_divided);
	int number_intersections = 0;
	all_routes_distances.push_back(each_route_distance);
	all_first_last_cust_demands.push_back(average_first_last_cust_demand);
	all_furthest_client_demand.push_back(demand_furthest);
	all_degree_of_neighborhood.push_back(degree_of_neighborhood);
	all_routes_degree_of_neighborhood.push_back(each_route_don);
	all_routes_two_don.push_back(each_route_two_don);
	
	for(int i = 0; i < vehicles_edge_list.size() - 1; i++){
		for(const auto& edges : vehicles_edge_list[i]){
			for(int j = 0; j < vehicles_edge_list.size(); j++){
				if(j != i){
					for(const auto& other_edges : vehicles_edge_list[j]){
						if(doIntersect(edges.first->pos, edges.second->pos, other_edges.first->pos, other_edges.second->pos)){
							number_intersections++;
						}
					}
				}
			}
		}
	}
	
	solution_intersections.push_back(number_intersections/2);
	vector<int> values;
	
	for(std::map<string,int>::iterator it = new_parameters.begin(); it != new_parameters.end(); ++it) {
		new_parameters_names.push_back(it->first);
		values.push_back(it->second);
	}
	all_new_parameters.push_back(values);
	resulting_data = {solution_intersections, longest_customers_distances, smallest_customers_distances, depot_customers_distances, 
	all_average_gravity, all_routes_gravity_x, all_routes_gravity_y, all_central_line_distances, all_routes_radians, all_gravity_radians, 
	all_route_depot_distances, all_routes_customers, all_new_parameters, new_parameters_names, all_routes_distances, depot_customers_distances_divided, 
	longest_interior_distances, smallest_interior_distances, all_first_last_cust_demands, all_furthest_client_demand, all_degree_of_neighborhood, 
	all_routes_degree_of_neighborhood, all_routes_two_don};
	
	return resulting_data;
}

STOP_values solution_to_other_parameters(std::shared_ptr <Solution> initial_solution, Instance& instance){
	
	DDOP_values all_parameters = solution_determine_other_parameters(initial_solution, instance);
	
	vector<float> s1_solution_intersections(all_parameters.solution_intersections.size());
	vector<float> s2_longest_distance_route(all_parameters.solution_intersections.size());
	vector<float> s3_average_depot_distance_route(all_parameters.solution_intersections.size());
	vector<float> s4_average_gravity_distance = all_parameters.all_average_gravity;
	vector<vector<vector<float>>> s5 = all_parameters.all_central_line_distances;
	vector<float> s5_average_width_route(all_parameters.solution_intersections.size());
	vector<vector<vector<float>>> s6 = all_parameters.all_routes_radians;
	vector<float> s6_average_span_route(all_parameters.solution_intersections.size());
	vector<float> s7_average_width_compactness(all_parameters.solution_intersections.size());
	vector<vector<vector<float>>> s8 = all_parameters.all_gravity_radians;
	vector<float> s8_average_span_compactness(all_parameters.solution_intersections.size());
	vector<vector<vector<float>>> s9 = all_parameters.all_route_depot_distances;
	vector<float> s9_average_depth_route(all_parameters.solution_intersections.size());
	vector<vector<int>> s10 = all_parameters.all_routes_customers;
	vector<float> s10_std_dev_customers_route(all_parameters.solution_intersections.size());
	
	vector<float> routes_min(all_parameters.solution_intersections.size());
	vector<float> routes_max(all_parameters.solution_intersections.size());
	vector<float> routes_mean(all_parameters.solution_intersections.size());
	vector<float> routes_std(all_parameters.solution_intersections.size());
	
	for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
		s1_solution_intersections[i] = all_parameters.solution_intersections[i]/instance.n_cust;
		s2_longest_distance_route[i] = sum(all_parameters.longest_customers_distances[i])/all_parameters.longest_customers_distances[i].size();
		s3_average_depot_distance_route[i] = sum(all_parameters.depot_customers_distances[i])/(2*all_parameters.depot_customers_distances[i].size());
		routes_mean[i] = sum(all_parameters.all_routes_distances[i])/all_parameters.all_routes_distances[i].size();
		routes_min[i] = *min_element(all_parameters.all_routes_distances[i].begin(), all_parameters.all_routes_distances[i].end());
		routes_max[i] = *max_element(all_parameters.all_routes_distances[i].begin(), all_parameters.all_routes_distances[i].end());
	}
	
	try{
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			vector<float> v_desv = all_parameters.all_routes_distances[i];
			routes_std[i] = desv_est(v_desv);
		}
	}catch(...){
		routes_std = vector<float>(all_parameters.solution_intersections.size());
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			routes_std[i] = 0;
		}
	}
	
	vector<float> s3_mod(all_parameters.solution_intersections.size());
	vector<float> s2_longest_mod(all_parameters.solution_intersections.size());
	vector<float> s2_shortest_mod(all_parameters.solution_intersections.size());
	vector<float> first_last_cust_demands(all_parameters.solution_intersections.size());
	vector<float> furthest_customers_demand(all_parameters.solution_intersections.size());
	
	for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
		s3_mod[i] = sum(all_parameters.depot_customers_distances_divided[i])/(2*all_parameters.depot_customers_distances_divided[i].size());
		s2_longest_mod[i] = sum(all_parameters.longest_interior_distances[i])/all_parameters.longest_interior_distances[i].size();
		s2_shortest_mod[i] = sum(all_parameters.smallest_interior_distances[i])/all_parameters.smallest_interior_distances[i].size();
		first_last_cust_demands[i] = sum(all_parameters.all_first_last_cust_demands[i])/all_parameters.all_first_last_cust_demands[i].size();
		furthest_customers_demand[i] = sum(all_parameters.all_furthest_client_demand[i])/all_parameters.all_furthest_client_demand[i].size();
	}
	
	vector<float> furthest_customers_demand_std(all_parameters.solution_intersections.size());
	try{
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			vector<float> v_desv = all_parameters.all_furthest_client_demand[i];
			furthest_customers_demand_std[i] = desv_est(v_desv);
		}
	}catch(...){
		furthest_customers_demand_std = vector<float>(all_parameters.solution_intersections.size());
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			furthest_customers_demand_std[i] = 0;
		}
	}
	
	vector<vector<float>> degree_of_neighborhood = all_parameters.all_degree_of_neighborhood;
	
	vector<string> degree_names(degree_of_neighborhood[0].size()); // duda acá
	for(int i = 0; i < degree_of_neighborhood[0].size(); i++){
		degree_names[i] = "DON_" + i;
	} // Duda aca
	
	vector<float> routes_degree_of_neighborhood(all_parameters.solution_intersections.size());
	vector<float> routes_degree_of_neighborhood_std(all_parameters.solution_intersections.size());
	vector<float> routes_two_degree_of_neighborhood(all_parameters.solution_intersections.size());
	vector<float> routes_two_degree_of_neighborhood_std(all_parameters.solution_intersections.size());
	
	for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
		routes_degree_of_neighborhood[i] = sum(all_parameters.all_routes_degree_of_neighborhood[i])/all_parameters.all_routes_degree_of_neighborhood[i].size();
		routes_two_degree_of_neighborhood[i] = sum(all_parameters.all_routes_two_don[i])/all_parameters.all_routes_two_don[i].size();
	}
	
	try{
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			vector<float> v_desv = all_parameters.all_routes_degree_of_neighborhood[i];
			routes_degree_of_neighborhood_std[i] = desv_est(v_desv);
		}
	}catch(...){
		routes_degree_of_neighborhood_std = vector<float>(all_parameters.solution_intersections.size());
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			routes_degree_of_neighborhood_std[i] = 0;
		}
	}
	
	try{
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			vector<float> v_desv = all_parameters.all_routes_two_don[i];
			routes_two_degree_of_neighborhood_std[i] = desv_est(v_desv);
		}
	}catch(...){
		routes_two_degree_of_neighborhood_std = vector<float>(all_parameters.solution_intersections.size());
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			routes_two_degree_of_neighborhood_std[i] = 0;
		}
	}
	
	// Average Smallest distance
	vector<float> s2_smallest_distance_route(all_parameters.solution_intersections.size());
	// Min largest distance
	vector<float> s21(all_parameters.solution_intersections.size());
	// Max largest distance
	vector<float> s22(all_parameters.solution_intersections.size());
	// Min smallest distance
	vector<float> s23(all_parameters.solution_intersections.size());
	// Max smallest distance
	vector<float> s24(all_parameters.solution_intersections.size());
	// Standard deviation largest distance
	vector<float> std_dev_largest_distance(all_parameters.solution_intersections.size());
	// Standard deviation smallest distance
	vector<float> std_dev_smallest_distance(all_parameters.solution_intersections.size());
	// Min distance to depot
	vector<float> s31(all_parameters.solution_intersections.size());
	// Max distance to depot
	vector<float> s32(all_parameters.solution_intersections.size());
	
	for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
		s2_smallest_distance_route[i] = sum(all_parameters.smallest_customers_distances[i])/all_parameters.smallest_customers_distances[i].size();
		s21[i] = *min_element(all_parameters.longest_customers_distances[i].begin(), all_parameters.longest_customers_distances[i].end());
		s22[i] = *max_element(all_parameters.longest_customers_distances[i].begin(), all_parameters.longest_customers_distances[i].end());
		s23[i] = *min_element(all_parameters.smallest_customers_distances[i].begin(), all_parameters.smallest_customers_distances[i].end());
		s24[i] = *max_element(all_parameters.smallest_customers_distances[i].begin(), all_parameters.smallest_customers_distances[i].end());
		s31[i] = (*min_element(all_parameters.depot_customers_distances[i].begin(), all_parameters.depot_customers_distances[i].end()))/2;
		s32[i] = (*max_element(all_parameters.depot_customers_distances[i].begin(), all_parameters.depot_customers_distances[i].end()))/2;
	}
	
	try{
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			vector<float> v_desv = all_parameters.longest_customers_distances[i];
			std_dev_largest_distance[i] = desv_est(v_desv);
		}
	}catch(...){
		std_dev_largest_distance = vector<float>(all_parameters.solution_intersections.size());
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			std_dev_largest_distance[i] = 0;
		}
	}
	
	try{
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			vector<float> v_desv = all_parameters.smallest_customers_distances[i];
			std_dev_smallest_distance[i] = desv_est(v_desv);
		}
	}catch(...){
		std_dev_smallest_distance = vector<float>(all_parameters.solution_intersections.size());
		for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
			std_dev_smallest_distance[i] = 0;
		}
	}
	
	// Std width per route
	vector<float> s5_std(all_parameters.solution_intersections.size());
	// Min width per route
	vector<float> s51(all_parameters.solution_intersections.size());
	// Max width per route
	vector<float> s52(all_parameters.solution_intersections.size());
	// Std span per route
	vector<float> s6_std(all_parameters.solution_intersections.size());
	// Min span per route
	vector<float> s61(all_parameters.solution_intersections.size());
	// Max span per route
	vector<float> s62(all_parameters.solution_intersections.size());
	// Std width compactness per route
	vector<float> s7_std(all_parameters.solution_intersections.size());
	// Min width compactness per route
	vector<float> s71(all_parameters.solution_intersections.size());
	// Max width compactness per route
	vector<float> s72(all_parameters.solution_intersections.size());
	// Std radian compactness per route
	vector<float> s8_std(all_parameters.solution_intersections.size());
	// Min radian compactness per route
	vector<float> s81(all_parameters.solution_intersections.size());
	// Max radian compactness per route
	vector<float> s82(all_parameters.solution_intersections.size());
	// Std depth per route
	vector<float> s9_std(all_parameters.solution_intersections.size());
	// Min depth per route
	vector<float> s91(all_parameters.solution_intersections.size());
	// Max depth per route
	vector<float> s92(all_parameters.solution_intersections.size());
	
	vector<vector<int>> client_pairs = all_parameters.all_new_parameters;
	vector<string> client_pairs_names = all_parameters.new_parameters_names;
	
	for(int i = 0; i < all_parameters.solution_intersections.size(); i++){
		vector<float> average_width_route;
		vector<float> average_span_route;
		vector<float> average_width_compactness;
		vector<float> average_span_compactness;
		vector<float> average_depth_route;
		vector<float> std_dev_route;
		
		for(int j = 0; j < s5[i].size(); j++){
			average_width_route.push_back(*max_element(s5[i][j].begin(), s5[i][j].end()) - *min_element(s5[i][j].begin(), s5[i][j].end()));
			vector<float> width_compactness;
			for(int z = 0; z < s5[i][j].size(); z++){
				width_compactness.push_back(abs(s5[i][j][z]));
			}
			average_width_compactness.push_back(sum(width_compactness));
		}
		
		s5_average_width_route[i] = sum(average_width_route)/s5[i].size();
		s51[i] = *min_element(average_width_route.begin(), average_width_route.end());
		s52[i] = *max_element(average_width_route.begin(), average_width_route.end());
		try{
			s5_std[i] = desv_est(average_width_route);
		}catch(...){
			s5_std[i] = 0;
		}
		
		s7_average_width_compactness[i] = sum(average_width_compactness)/instance.n_cust;
		s71[i] = *min_element(average_width_compactness.begin(), average_width_compactness.end());
		s72[i] = *max_element(average_width_compactness.begin(), average_width_compactness.end());
		try{
			s7_std[i] = desv_est(average_width_compactness);
		}catch(...){
			s7_std[i] = 0;
		}
		
		for(int w = 0; w < s6[i].size(); w++){
			try{
				average_span_route.push_back(*max_element(s6[i][w].begin(), s6[i][w].end()));
			}catch(...){
				average_span_route.push_back(0.0);
			}
		}
		
		s6_average_span_route[i] = sum(average_span_route)/s6[i].size();
		s61[i] = *min_element(average_span_route.begin(), average_span_route.end());
		s62[i] = *max_element(average_span_route.begin(), average_span_route.end());
		try{
			s6_std[i] = desv_est(average_span_route);
		}catch(...){
			s6_std[i] = 0;
		}
		
		for(int k = 0; k < s8[i].size(); k++){
			average_span_compactness.push_back(sum(s8[i][k]));
		}
		
		s8_average_span_compactness[i] = sum(average_span_compactness)/instance.n_cust;
		s81[i] = *min_element(average_span_compactness.begin(), average_span_compactness.end());
		s82[i] = *max_element(average_span_compactness.begin(), average_span_compactness.end());
		try{
			s8_std[i] = desv_est(average_span_compactness);
		}catch(...){
			s8_std[i] = 0;
		}
		
		for(int q = 0; q < s9[i].size(); q++){
			try{
				average_depth_route.push_back(*max_element(s9[i][q].begin(), s9[i][q].end()));
			}catch(...){
				average_depth_route.push_back(0.0);
			}
		}
		
		s9_average_depth_route[i] = sum(average_depth_route)/s9[i].size();
		s91[i] = *min_element(average_depth_route.begin(), average_depth_route.end());
		s92[i] = *max_element(average_depth_route.begin(), average_depth_route.end());
		try{
			s9_std[i] = desv_est(average_depth_route);
		}catch(...){
			s9_std[i] = 0;
		}
		
		for(int h = 0; h < s10[i].size(); h++){
			std_dev_route.push_back(pow(s10[i][h] - (instance.n_cust/s10[i].size()), 2));
		}
		s10_std_dev_customers_route[i] = sqrt(sum(std_dev_route)/s10[i].size());
	}
	
	vector<float> features;
	vector<string> features_names;
	
	int i = all_parameters.solution_intersections.size() - 1;
	
	features = {s1_solution_intersections[i], s2_longest_distance_route[i], s3_average_depot_distance_route[i], s4_average_gravity_distance[i], 
						s5_average_width_route[i], s6_average_span_route[i], s7_average_width_compactness[i], s8_average_span_compactness[i], 
						s9_average_depth_route[i], s10_std_dev_customers_route[i], s21[i], s22[i], s23[i], s24[i],
						s31[i], s32[i], s51[i], s52[i], s61[i], s62[i], s71[i], s72[i], 
						s81[i], s82[i], s91[i], s92[i], 
						std_dev_largest_distance[i], std_dev_smallest_distance[i],
						s5_std[i], s6_std[i], s7_std[i], s8_std[i], s9_std[i], 
						routes_min[i], routes_max[i], routes_mean[i], routes_std[i], 
						s3_mod[i], s2_longest_mod[i], s2_shortest_mod[i], first_last_cust_demands[i],
						furthest_customers_demand[i], furthest_customers_demand_std[i], routes_degree_of_neighborhood[i], 
						routes_degree_of_neighborhood_std[i], routes_two_degree_of_neighborhood[i], routes_two_degree_of_neighborhood_std[i]};

	
	features_names = {"S1", "S2", "S3", "S4", "S5", "S6", "S7",
					 "S8", "S9", "S10", "S21", "S22", "S23", "S24", 
					 "S31", "S32", "S51",
					 "S52", "S61", "S62", "S71", "S72", "S81", "S82", "S91", "S92",
					 "S2_L_STD", "S2_S_STD", "S5_STD", "S6_STD", 
					 "S7_STD", "S8_STD", "S9_STD", 
					 "RL_MIN", "RL_MAX", "RL_MEAN", "RL_STD",
					 "S3_*", "S2_L_MOD", "S2_S_MOD", "F_L_CD", 
					 "FC_DEM_MEAN", "FC_DEM_STD", "R_DON", "R_DON_STD",
					 "R_2DON", "R_2DON_STD"};
					 
	STOP_values parameters_and_names = {features, features_names};
	return parameters_and_names;
}

void print_features(vector<shared_ptr <SampleInfo>> dataset_train, Instance& instance, string& method, string& file_parameters, string& direc_output){
	DTOP_values values;
	string filename;
	if(method == "tsp"){
		cout << "Entra a TSP" << endl;
		values = dataset_determine_tsp_parameters(dataset_train, instance);
		filename = direc_output + "train_tsp_features_" + file_parameters; 
	}else if(method == "vrp"){
		cout << "Entra a VRP" << endl;
		values = dataset_to_vrp_parameters(dataset_train, instance);
		filename = direc_output + "train_vrp_features_" + file_parameters; 
	}else if(method == "client-pairs"){
		cout << "Entra a Client" << endl;
		values = dataset_determine_client_pairs_parameters(dataset_train, instance);
		filename = direc_output + "train_client-pairs_features_" + file_parameters; 
	}else if(method == "tsp-vrp"){
		values = dataset_determine_tsp_vrp_parameters(dataset_train, instance);
		filename = direc_output + "train_tsp-vrp_features_" + file_parameters; 
	}else if(method == "tsp-client-pairs"){
		values = dataset_determine_tsp_client_pairs_parameters(dataset_train, instance);
		filename = direc_output + "train_tsp-client-pairs_features_" + file_parameters; 
	}else if(method == "vrp-client-pairs"){
		values = dataset_determine_vrp_client_pairs_parameters(dataset_train, instance);
		filename = direc_output + "train_vrp-client-pairs_features_" + file_parameters; 
	}else if(method == "tsp-vrp-client-pairs"){
		values = dataset_determine_tsp_vrp_client_pairs_parameters(dataset_train, instance);
		filename = direc_output + "train_tsp-vrp-client-pairs_features_" + file_parameters; 
	}else{
		cout << "Invalid method" << endl;
		exit(1);
	}
	
	ofstream file(filename);
	for(int i = 0; i < values.features_names.size() - 1; i++){
		file << values.features_names[i] << ",";
	}
	
	file << values.features_names[values.features_names.size() - 1] << ",";
	file << "final_cost" << endl;  //Add final cost ONLY FOR TRAINING SET
	
	for(int i = 0; i < dataset_train.size(); i++){
		file << i << ",";
		for(int j = 0; j < values.features[i].size() - 1; j++){
			file << values.features[i][j] << ",";
		}
		file << values.features[i][values.features[i].size() - 1] << ",";
		file << to_string(dataset_train[i]->final_cost) << endl; //Add final cost ONLY FOR TRAINING SET
	}
	
	file.close();
	cout << "Termina de crear archivo" << endl;
}

DTOP_values print_vnd_features(vector<shared_ptr <SampleInfo>> dataset_train, Instance& instance, string& method, string& file_parameters, string& direc_output, int number_of_vnds){
    DTOP_values values;
    string filename;
    if(method == "tsp"){
        cout << "Entra a TSP" << endl;
        values = dataset_determine_tsp_parameters(dataset_train, instance);
        filename = direc_output + "train_tsp_features_" + file_parameters;
    }else if(method == "vrp"){
        cout << "Entra a VRP" << endl;
        values = dataset_to_vrp_parameters(dataset_train, instance);
        filename = direc_output + "train_vrp_features_" + file_parameters;
    }else if(method == "client-pairs"){
        cout << "Entra a Client" << endl;
        values = dataset_determine_client_pairs_parameters(dataset_train, instance);
        filename = direc_output + "train_client-pairs_features_" + file_parameters;
    }else if(method == "tsp-vrp"){
        cout << "Entra a TSP-VRP" << endl;
        values = dataset_determine_tsp_vrp_parameters(dataset_train, instance);
        filename = direc_output + "train_tsp-vrp_features_" + file_parameters;
    }else if(method == "tsp-client-pairs"){
        values = dataset_determine_tsp_client_pairs_parameters(dataset_train, instance);
        filename = direc_output + "train_tsp-client-pairs_features_" + file_parameters;
    }else if(method == "vrp-client-pairs"){
        values = dataset_determine_vrp_client_pairs_parameters(dataset_train, instance);
        filename = direc_output + "train_vrp-client-pairs_features_" + file_parameters;
    }else if(method == "tsp-vrp-client-pairs"){
        values = dataset_determine_tsp_vrp_client_pairs_parameters(dataset_train, instance);
        filename = direc_output + "train_tsp-vrp-client-pairs_features_" + file_parameters;
    }else{
        cout << "Invalid method" << endl;
        exit(1);
    }

    ofstream file(filename);
    for(int i = 0; i < values.features_names.size() - 1; i++){
        file << values.features_names[i] << ",";
    }

    file << values.features_names[values.features_names.size() - 1] << ",";
    for (int w=0; w < number_of_vnds; w++){
        //Add final cost ONLY FOR TRAINING SET
        file << "VND_" << w << "," << "final_cost_" << w;
        if (w != number_of_vnds-1){
            file << ",";
        }
    }
    file << endl;
    //file << "VND_1,final_cost_1,VND_2,final_cost_2,VND_3,final_cost_3,VND_4,final_cost_4,VND_5,final_cost_5" << endl;

    for(int i = 0; i < dataset_train.size(); i++){
        file << i << ",";
        for(int j = 0; j < values.features[i].size() - 1; j++){
            file << values.features[i][j] << ",";
        }
        file << values.features[i][values.features[i].size() - 1] << ",";
        for (int w=0; w < number_of_vnds; w++){
            //Add final cost ONLY FOR TRAINING SET
            file << w << "," <<dataset_train[i]->final_costs[w];
            if (w != number_of_vnds-1){
                file << ",";
            }
        }
        file << endl;
        /*
        file << dataset_train[i]->neighborhood_1 << "," <<dataset_train[i]->final_cost_1<<"," <<
        dataset_train[i]->neighborhood_2 << "," <<dataset_train[i]->final_cost_2<<"," <<
        dataset_train[i]->neighborhood_3 << "," <<dataset_train[i]->final_cost_3<<"," <<
        dataset_train[i]->neighborhood_4 << "," <<dataset_train[i]->final_cost_4<<"," <<
        dataset_train[i]->neighborhood_5 << "," <<dataset_train[i]->final_cost_5<<endl;
         */
    }

    file.close();
    cout << "Termina de crear archivo de train features" << endl;
    return values;
}

DTOP_values create_solution_features(vector<shared_ptr <SampleInfo>> dataset_train, Instance& instance, string& method){
    DTOP_values values;
    if(method == "tsp"){
        values = dataset_determine_tsp_parameters(dataset_train, instance);
    }else if(method == "vrp"){
        values = dataset_to_vrp_parameters(dataset_train, instance);
    }else if(method == "client-pairs"){
        values = dataset_determine_client_pairs_parameters(dataset_train, instance);
    }else if(method == "tsp-vrp"){
        values = dataset_determine_tsp_vrp_parameters(dataset_train, instance);
    }else if(method == "tsp-client-pairs"){
        values = dataset_determine_tsp_client_pairs_parameters(dataset_train, instance);
    }else if(method == "vrp-client-pairs"){
        values = dataset_determine_vrp_client_pairs_parameters(dataset_train, instance);
    }else if(method == "tsp-vrp-client-pairs"){
        values = dataset_determine_tsp_vrp_client_pairs_parameters(dataset_train, instance);
    }else{
        cout << "Invalid method" << endl;
        exit(1);
    }
    return values;
}

void print_validation_features(vector<shared_ptr <SampleInfo>> dataset_train, Instance& instance, string& method, string& file_parameters, string& direc_output){
	DTOP_values values;
	string filename;
	
	if(method == "tsp"){
		cout << "Entra a TSP" << endl;
		values = dataset_determine_tsp_parameters(dataset_train, instance);
		filename = direc_output + "validation_tsp_features_" + file_parameters; 
	}else if(method == "vrp"){
		cout << "Entra a VRP" << endl;
		values = dataset_to_vrp_parameters(dataset_train, instance);
		filename = direc_output + "validation_vrp_features_" + file_parameters; 
	}else if(method == "client-pairs"){
		cout << "Entra a Client" << endl;
		values = dataset_determine_client_pairs_parameters(dataset_train, instance);
		filename = direc_output + "validation_client-pairs_features_" + file_parameters; 
	}else if(method == "tsp-vrp"){
        cout << "Entra a TSP-VRP" << endl;
		values = dataset_determine_tsp_vrp_parameters(dataset_train, instance);
		filename = direc_output + "validation_tsp-vrp_features_" + file_parameters; 
	}else if(method == "tsp-client-pairs"){
		values = dataset_determine_tsp_client_pairs_parameters(dataset_train, instance);
		filename = direc_output + "validation_tsp-client-pairs_features_" + file_parameters; 
	}else if(method == "vrp-client-pairs"){
		values = dataset_determine_vrp_client_pairs_parameters(dataset_train, instance);
		filename = direc_output + "validation_vrp-client-pairs_features_" + file_parameters; 
	}else if(method == "tsp-vrp-client-pairs"){
		values = dataset_determine_tsp_vrp_client_pairs_parameters(dataset_train, instance);
		filename = direc_output + "validation_tsp-vrp-client-pairs_features_" + file_parameters; 
	}else{
		cout << "Invalid method" << endl;
		exit(1);
	}
	
	ofstream file(filename);
	for(int i = 0; i < values.features_names.size() - 1; i++){
		file << values.features_names[i] << ",";
	}
	
	file << values.features_names[values.features_names.size() - 1] << endl;
	
	
	for(int i = 0; i < dataset_train.size(); i++){
		file << i << ",";
		for(int j = 0; j < values.features[0].size() - 1; j++){
			file << values.features[i][j] << ",";
		}
		file << values.features[i][values.features[0].size() - 1] << endl;
	}
	
	file.close();
    cout << "Termina de crear archivo de validation features" << endl;
}





