#include "../methods/CircleSector.h"
#include "./solution.h"
#include <boost/format.hpp>
#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <memory>

#define unlikely(condition) __builtin_expect(static_cast<bool>(condition), 0)

using namespace std;

Vehicle::Vehicle()
{
}

Vehicle::Vehicle(const shared_ptr<Vehicle> other_vehicle)
    : Q(other_vehicle->Q)
    , customers(other_vehicle->customers.size())
    , cummulative_demands(other_vehicle->cummulative_demands)
    , load(other_vehicle->load)
    , currentCustomer(other_vehicle->currentCustomer)
    , cost(other_vehicle->cost)
{
    for(int i = 0; i < other_vehicle->customers.size(); i++) {
        customers[i] = make_shared<Customer>(other_vehicle->customers[i]);
    }
    this->vehicle_address = this;
}

Vehicle::Vehicle(int Q)
{
    /**
     *  Vehicle object that carries the goods
     *  @param Q maximum number of items carried by the vehicle.
     */

    this->Q = Q;
    this->customers = {};
    this->cummulative_demands = {};
    this->load = 0;
    this->currentCustomer = 0; // starting location is depot
    this->polar_angle_barycenter = 0.0;
    this->cost = 0;
    this->vehicle_address = this;
}

Vehicle::~Vehicle()
{
    /*
     * for(int i = this->customers.size() - 1; i > -1; i--){
            Customer* c = this->customers[i];
            delete c;
            this->customers[i] = nullptr;
    }
    */
}

void Vehicle::add_customer(shared_ptr <Customer> customer)
{
    /**
     * Add customer to vehicle route
     * @param customer Object of class Node
     */

    this->customers.push_back(customer);
    this->load += customer->demand;
    this->currentCustomer = customer->id;
}

void Vehicle::insert_customer(shared_ptr <Customer> customer, int position)
{
    this->customers.insert(this->customers.begin() + position, customer); // IMPORTANT, COMPLEXITY
    this->load += customer->demand;
    customer->isRouted = true;
}

void Vehicle::remove_customer(int position)
{
    this->load -= this->customers[position]->demand;
    this->customers[position]->isRouted = false;
    this->customers.erase(this->customers.begin() + position);
}

void Vehicle::temporary_remove_customer(int position)
{
    this->load -= this->customers[position]->demand;
    this->customers.erase(this->customers.begin() + position);
}

bool Vehicle::check_constraint(int newLoad)
{
    return (this->load + newLoad) <= this->Q;
}

bool Vehicle::current_load_below_Q()
{
    return (this->load <= this->Q);
}

void Vehicle::compute_cummulative_demands()
{
    this->cummulative_demands.clear(); // Complexity?
    int idx = 0;
    for(const auto& customer : this->customers) {
        this->cummulative_demands.push_back(customer->demand);
        if(idx > 0)
            this->cummulative_demands[idx] += cummulative_demands[idx - 1];
        idx++;
    }
}

int Vehicle::get_demand_between(int left, int right)
{
    int demand = this->cummulative_demands[right];
    if(left > 0)
        demand -= this->cummulative_demands[left - 1];
    return demand;
}

void Vehicle::update_load()
{
    this->load = 0;
    bool first_it = true;
    for(const auto& customer : this->customers) {
        this->load += customer->demand;
        if(customer->id != 0) {
            if(first_it) {
                this->sector.initialize(customer->polar_angle);
                first_it = false;
            } else {
                this->sector.extend(customer->polar_angle);
            }
        }
    }
}

float Vehicle::compute_cost(const Instance& instance)
{
    this->cost = 0;
    int prev = -1;
    int curr = -1;
    bool first_it = true;
    for(int i = 0; i < this->customers.size(); i++) {
        curr = this->customers[i]->id;
        if(prev != -1) {
            this->cost += instance.dist_matrix[prev][curr];
        }
        // cout << "Customer " << i << "Angle " << this->customers[i]->polar_angle << endl;
        if(this->customers[i]->id != 0) {
            if(first_it) {
                this->sector.initialize(this->customers[i]->polar_angle);
                first_it = false;
            } else {
                this->sector.extend(this->customers[i]->polar_angle);
            }
        }
        prev = curr;
    }
    return this->cost;
}

Customer::Customer()
{
}

Customer::Customer(shared_ptr <Customer>other_customer)
    : id(other_customer->id)
    , demand(other_customer->demand)
    , vehicle_route(other_customer->vehicle_route)
    , isRouted(other_customer->isRouted)
    , x(other_customer->x)
    , y(other_customer->y)
    , pos(other_customer->pos)
    , dist_to_origin(other_customer->dist_to_origin)
    , dist_to_prev(other_customer->dist_to_prev)
    , dist_to_forw(other_customer->dist_to_forw)
    , cumu_demand(other_customer->cumu_demand)
    , cumu_dist(other_customer->cumu_dist)
    , polar_angle(other_customer->polar_angle)
{
    this->customer_address = this;
}

Customer::Customer(int id, int demand, vector<vector<int>> cust_locations)
{
    this->id = id;
    this->demand = demand;
    this->vehicle_route = -1;
    this->isRouted = false;
    this->x = cust_locations[0][id];
    this->y = cust_locations[1][id];
    this->pos = { x, y };
    // new attributes
    this->dist_to_origin = NULL;
    this->dist_to_prev = NULL;
    this->dist_to_forw = NULL;
    this->cumu_demand = NULL;
    this->cumu_dist = NULL;
    this->polar_angle = CircleSector::positive_mod(
        32768. * atan2(this->y - cust_locations[1][0], this->x - cust_locations[0][0]) / 3.14159265359);
    this->customer_address = this;
}

Customer::~Customer()
{
}

Solution::Solution()
{
}

Solution::Solution(const Instance& instance, int n_vehicles)
{
    this->vehicles = vector<shared_ptr<Vehicle>>();
    //this->vehicles = vector<Vehicle*>();
    this->vehicles.reserve(n_vehicles);
    for(int n = 0; n < n_vehicles; n++) {
        //this->vehicles.push_back(new Vehicle(instance.Q));
        //this->vehicles.push_back(make_unique(Vehicle(instance.Q)));
        this->vehicles.push_back(std::make_shared<Vehicle>(instance.Q));
        //this->vehicles[n] = new Vehicle(instance.Q);
    }

    //this->customers = vector<Customer*>();
    this->customers = vector<shared_ptr<Customer>>();
    this->customers.reserve(instance.n_cust);
    for(int n = 0; n < instance.n_cust; n++) {
        //this->customers.push_back(new Customer(n, instance.demands[n], instance.cust_locations));
        //this->customers[n] = new Customer(n, instance.demands[n], instance.cust_locations);
        this->customers.push_back(std::make_shared<Customer>(n, instance.demands[n], instance.cust_locations));
    }

    this->cost = 0;
    this->n_vehicles = n_vehicles;
    this->solution_address = this;
}

Solution::Solution(shared_ptr<Solution> otherSolution)
    : vehicles(otherSolution->vehicles.size())
    , customers(otherSolution->customers.size())
    , cost(otherSolution->cost)
    , n_vehicles(otherSolution->n_vehicles)
{
    for(int i = 0; i < otherSolution->vehicles.size(); i++) {
        //vehicles[i] = new Vehicle(*otherSolution->vehicles[i]);
        vehicles[i] = make_shared <Vehicle>(otherSolution->vehicles[i]);
    }

    for(int i = 0; i < otherSolution->customers.size(); i++) {
        //customers[i] = new Customer(*otherSolution->customers[i]);
        customers[i] = make_shared <Customer>(otherSolution->customers[i]);
    }
    this->solution_address = this;
}

void Solution::cloneSolution(std::shared_ptr<Solution> otherSolution) {
    for(int i = 0; i < otherSolution->customers.size(); i++) {
        this->customers[i]->vehicle_route = otherSolution->customers[i]->vehicle_route;
    }
    if (this->vehicles.size() < otherSolution->vehicles.size()) {
        while (this->vehicles.size() < otherSolution->vehicles.size()){
            this->vehicles.push_back(std::make_shared<Vehicle>(otherSolution->vehicles[0]->Q));
        }
    }
    else if (this->vehicles.size() > otherSolution->vehicles.size()){
        while (this->vehicles.size() > otherSolution->vehicles.size()){
            this->vehicles.pop_back();
        }
    }
        for (int i = 0; i < otherSolution->vehicles.size(); i++) {
            if (this->vehicles[i]->customers.size() < otherSolution->vehicles[i]->customers.size()) {
                while (this->vehicles[i]->customers.size() < otherSolution->vehicles[i]->customers.size()){
                    this->vehicles[i]->add_customer(this->customers[otherSolution->vehicles[i]->customers[0]->id]);
                }
            }
            else if (this->vehicles[i]->customers.size() > otherSolution->vehicles[i]->customers.size()){
                while (this->vehicles[i]->customers.size() > otherSolution->vehicles[i]->customers.size()){
                    this->vehicles[i]->customers.pop_back();
                }
            }
            this->vehicles[i]->load = 0;
            for (int j = 0; j < otherSolution->vehicles[i]->customers.size(); j++) {
                this->vehicles[i]->customers[j] = this->customers[otherSolution->vehicles[i]->customers[j]->id];
                this->vehicles[i]->load += this->vehicles[i]->customers[j]->demand;
            }
            /*
            this->vehicles[i]->customers.clear();
            this->vehicles[i]->customers.reserve(otherSolution->vehicles[i]->customers.size());
            this->vehicles[i]->load = 0;
            for (int j = 0; j < otherSolution->vehicles[i]->customers.size(); j++) {
                this->vehicles[i]->add_customer(
                        this->customers[otherSolution->vehicles[i]->customers[j]->id]);
            }
            */
        }

    this->cost = otherSolution->cost;
}





Solution::~Solution()
{
/*
    for(int i = this->customers.size() - 1; i > -1; i--) {
        Customer* c = this->customers[i];
        delete c;
        //delete this->customers[i]->customer_address;
        this->customers[i] = nullptr;
    }
    customers.clear();
    for(int i = this->vehicles.size() - 1; i > -1; i--) {
        Vehicle* v = this->vehicles[i];
        delete v;
        //delete this->vehicles[i]->vehicle_address;
        this->vehicles[i] = nullptr;
    }
    vehicles.clear();
*/
/*
    for(auto customer: this->customers){
        delete customer;
    }
    for(auto vehicle: this->vehicles){
        delete vehicle;
    }
*/

}

bool Solution::is_feasible()
{
    for(auto& vehicle : this->vehicles) {
        if(!vehicle->current_load_below_Q()) {
            cout << "unfeasible current load " << endl;
            return false;
        }
    }

    for(auto& customer : this->customers){
        if(!customer->isRouted){
            cout << "Customers not routed "<< endl;
            return false;
        }
    }

    return true;
}

float Solution::compute_cost(const Instance& instance)
{
    this->cost = 0;
    int prev = -1;
    int curr = -1;
    bool first_it;
    float vehicle_cost;
    for(int j=0; j < this->vehicles.size(); j++) {
        vehicle_cost = 0;
        prev = -1;
        curr = -1;
        first_it = true;
        //this->cost += this->vehicles[j]->compute_cost(instance);
        
        for(int i = 0; i < this->vehicles[j]->customers.size(); i++) {
            curr = this->vehicles[j]->customers[i]->id;
            if(prev != -1) {
                this->cost += instance.dist_matrix[prev][curr];
                vehicle_cost += instance.dist_matrix[prev][curr];
            }
            if(this->vehicles[j]->customers[i]->id != 0){
            if(first_it) {
                this->vehicles[j]->sector.initialize( this->vehicles[j]->customers[i]->polar_angle);
                first_it = false;
            } else {
                this->vehicles[j]->sector.extend( this->vehicles[j]->customers[i]->polar_angle);
            }
			this->vehicles[j]->customers[i]->vehicle_route = j;
			
        }
			this->customers[this->vehicles[j]->customers[i]->id]->vehicle_route = j;
            prev = curr;
        }
		
        this->vehicles[j]->cost = vehicle_cost;
        
    }
    return this->cost;
}

float Solution::update_cost(const Instance& instance)
{
    this->cost = 0;
    int prev = -1;
    int curr = -1;
    bool first_it;
    //float vehicle_cost;
    for(int j=0; j < this->vehicles.size(); j++) {
        //vehicle_cost = 0;
        prev = -1;
        curr = -1;
        first_it = true;
        //this->cost += this->vehicles[j]->compute_cost(instance);

        for(int i = 0; i < this->vehicles[j]->customers.size(); i++) {
            curr = this->vehicles[j]->customers[i]->id;
            if(prev != -1) {
                this->cost += instance.dist_matrix[prev][curr];
                //vehicle_cost += instance.dist_matrix[prev][curr];
            }
            this->customers[this->vehicles[j]->customers[i]->id]->vehicle_route = j;
            prev = curr;
        }
        //this->vehicles[j]->cost = vehicle_cost;
    }
    return this->cost;
}

/*
float Solution::compute_cost(const Instance& instance) {
    this->cost = 0;
    int prev = -1;
    int curr = -1;
    //for(auto& vehicle : this->vehicles){
    for (int j= 0; j < this->vehicles.size() ; j++){
        prev = -1;
        curr = -1;
        for(int i = 0; i < this->vehicles[j]->customers.size() ; i++){
            curr = this->vehicles[j]->customers[i]->id;
            if(prev != -1){
                this->cost += instance.dist_matrix[prev][curr];
                        }
            prev = curr;
        }
    }
    return this->cost;
}
*/
void Solution::delete_empty_vehicles()
{
    //vector<Vehicle*> vehicles2;
    vector<shared_ptr<Vehicle>> vehicles2;
	bool found_empty_vehicle = false;
    for(auto& vehicle : this->vehicles) {
        if(vehicle->customers.size() <= 2)
		{
			found_empty_vehicle = true;
		}
		else{
            vehicles2.push_back(vehicle);
		}
    }
	this->vehicles = vehicles2;
	if (found_empty_vehicle == true){
	    //Update all routes and customers when a vehicle is removed from the solution
		for (int w=0; w < this->vehicles.size(); w++){
			for (int v=0; v < this->vehicles[w]->customers.size(); v++){
				this->vehicles[w]->customers[v]->vehicle_route = w;
				this->customers[this->vehicles[w]->customers[v]->id]->vehicle_route = w;
			}
		}
	}
}

void Solution::compute_cummulative_demands()
{
    for(auto& vehicle : this->vehicles) {
        vehicle->compute_cummulative_demands();
    }
}

void Solution::append_route(const int route, const int route_to_append, const Instance& instance) {

    const auto route_end = this->vehicles[route]->customers[this->vehicles[route]->customers.size()-2]->id;
    const auto route_to_append_start = this->vehicles[route_to_append]->customers[1]->id;

    assert(route_end != instance.depot);
    assert(route_to_append_start != instance.depot);

    const auto delta = +instance.dist_matrix[route_end][route_to_append_start]
                       - instance.dist_matrix[route_end][instance.depot]
                       - instance.dist_matrix[instance.depot][route_to_append_start];
    this->cost += delta;

    this->vehicles[route]->customers.pop_back();
    for (int i=1; i < this->vehicles[route_to_append]->customers.size()-1; i++) {
        this->vehicles[route]->add_customer(this->customers[this->vehicles[route_to_append]->customers[i]->id]);
        this->customers[this->vehicles[route_to_append]->customers[i]->id]->vehicle_route = route;
    }

    this->vehicles[route]->add_customer(this->customers[0]);
    int remove_customers = this->vehicles[route_to_append]->customers.size()-1;
    for (int i=1; i < remove_customers ; i++) {
        this->vehicles[route_to_append]->temporary_remove_customer(1);
    }
}
/*const string& Solution::str() const{ // Complejidad - Creación de multiples cadenas
    if(this->cost == 0)
        return "Optimization not run yet. Run by calling Greedy_VRP.run().";

    string string_ch = "=================================================================\n";
    int n_customers;
    for(int v_id = 0; v_id < this->vehicles.size(); v_id++){
        if(this->vehicles[v_id]->customers.size() != 0){
            string_ch += (boost::format("Vehicle %1%:\n\t") % v_id).str();
            n_customers = this->vehicles[v_id].customers.size();
            for(int i = 0; i < n_customers; i++){
                if(i == n_customers - 1){
                    string_ch += to_string(this->vehicles[v_id].customers[i].id);
                }else {
                    string_ch += (boost::format("%1%->") % this->vehicles[v_id].customers[i].id).str();
                }
            }
            string_ch += '\n';
        }
    }
    string_ch += (boost::format("Best Cost: %1%\n") % this->cost).str();
    return string_ch;
}*/

shared_ptr <Solution> Solution::create_giant_tour(const Instance& instance) {
    int n_vehicles = 1;
    shared_ptr <Solution> giant_tour = make_shared <Solution>(instance, n_vehicles);
    giant_tour->vehicles[0]->customers.push_back(this->vehicles[0]->customers[0]);
    for(auto& vehicle : this->vehicles) {
        for (auto &customer : vehicle->customers) {
            if (customer->id != 0) {
                giant_tour->vehicles[0]->customers.push_back(customer);
            }
        }
    }
    giant_tour->vehicles[0]->customers.push_back(this->vehicles[0]->customers[0]);
    giant_tour->compute_cost(instance);
    return giant_tour;
}

ostringstream Solution::local_solutions() const
{ // Complejidad - metodo join

    vector<string> lines = {};
    lines.push_back(to_string(this->vehicles.size()));
    lines.push_back(to_string(this->cost));

    string chain;
    vector<int> customers_ids;
    for(auto& vehicle : this->vehicles) {
        customers_ids = vector<int>();
        for(auto& customer : vehicle->customers) {
            // if(customer->id != 99)
            customers_ids.push_back(customer->id);
        }
        chain = ""; // revisar este metodo
        int cont = 0;
        for(auto& id : customers_ids) {
            if(cont == customers_ids.size() - 1) {
                chain += to_string(id);
            } else {
                chain += (boost::format("%1% ") % id).str();
            }
            cont++;
        }
        lines.push_back(chain);
    }

    ostringstream vts =  ostringstream();

    // Convert all but the last element to avoid a trailing ","
    copy(lines.begin(), lines.end() - 1, ostream_iterator<string>(vts, "\n"));
    // Now add the last element with no delimiter
    vts << lines.back();

    return vts;
}

const int Solution::dummy_vertex = -1;


int Solution::get_route_index(const int customer) const {
    return customers[customer]->vehicle_route;
}

int Solution::get_next_vertex(const int route, const int vertex) const {

    assert(contains_vertex(route, vertex));
    if (unlikely(this->vehicles[route]->customers[vertex]->id == 0)) {
        return 1;
    } else {
        return vertex+1;
    }

}

int Solution::get_prev_vertex(const int route, const int vertex) const {

    assert(contains_vertex(route, vertex));

    if (unlikely(this->vehicles[route]->customers[vertex]->id == 0)) {
        return vehicles[route]->customers.size()-1;
    } else {
        return vertex-1;
    }

}

int Solution::get_prev_vertex(const int customer) const {
    return customer-1;
}

int Solution::get_next_vertex(const int customer) const {
    return customer+1;
}

bool Solution::contains_vertex(const int route, const int vertex) const {
    return customers[vertex]->vehicle_route == route || customers[vertex]->id == 0;
}

bool Solution::is_customer_in_solution(const int customer) const {
    assert(customer != 0);
    return this->customers[customer]->isRouted;
}




