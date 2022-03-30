#pragma once
#include <vector>
#include <string> 
#include <sstream>
#include "./instance.h"
#include "../methods/CircleSector.h"
#include <memory>

struct Point{
	float x;
	float y;
};

class Customer {
    public:
    int id;
    int demand;
	int vehicle_route;
    bool isRouted;
    float x;
    float y;
    Point pos;
    // new attributes
    float dist_to_origin;
    float dist_to_prev;
    float dist_to_forw;
    int cumu_demand;
    float cumu_dist;
    int polar_angle;

    Customer* customer_address;

    Customer();
	Customer(std::shared_ptr <Customer> other_customer);
    Customer(int id, int demand, std::vector<std::vector<int>> cust_locations);
    ~Customer();
};

class Vehicle{
    public:
    int Q;
    std::vector<std::shared_ptr <Customer>> customers;
    std::vector<int> cummulative_demands;
    int load;
    int currentCustomer;
    float cost;
    float polar_angle_barycenter;
    CircleSector sector;

    Vehicle();
	Vehicle(std::shared_ptr <Vehicle> other_vehicle);
    Vehicle(int Q);
    ~Vehicle();

    Vehicle* vehicle_address;

    void add_customer(std::shared_ptr <Customer> customer);
    void insert_customer(std::shared_ptr <Customer> customer, int position);
    void remove_customer(int position);
    void temporary_remove_customer(int position);
    bool check_constraint(int newLoad);
    bool current_load_below_Q();
    void compute_cummulative_demands();
    int get_demand_between(int left, int right);
    float compute_cost(const Instance& instance);
    void update_load();
};

class Solution {
    public:
    //std::vector<Vehicle*> vehicles;
    std::vector<std::shared_ptr <Vehicle>> vehicles;
    //std::vector<Customer*> customers;
    std::vector<std::shared_ptr <Customer>> customers;
    float cost;
    int n_vehicles;
	
	std::shared_ptr <Solution> create_giant_tour(const Instance& instance);


    static const int dummy_vertex;
    bool contains_vertex(int route, int vertex) const;
    int get_route_index(int customer) const;
    int get_next_vertex(int route, int vertex) const;
    int get_prev_vertex(int route, int vertex) const;
    int get_prev_vertex(int customer) const;
    int get_next_vertex(int customer) const;
    bool is_customer_in_solution(int customer) const;
    void cloneSolution(std::shared_ptr<Solution> otherSolution);
    void append_route(const int route, const int route_to_append, const Instance& instance);





    Solution();
    Solution(const Instance& instance, int n_vehicles);
    Solution(std::shared_ptr<Solution> otherSolution);
    ~Solution();


    Solution* solution_address;

    bool is_feasible();
    float compute_cost(const Instance& instance);
    float update_cost(const Instance& instance);
    void delete_empty_vehicles();
    void compute_cummulative_demands();
    //const std::string& str() const;
    std::ostringstream local_solutions() const;

};