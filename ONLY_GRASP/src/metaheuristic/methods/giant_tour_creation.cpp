#include "./giant_tour_creation.h"
#include <algorithm>
#include <cstdlib>
#include <vector>

using namespace std;
GiantTour::GiantTour()
    : giant_tour_solution(make_shared <Solution>())
{
}

GiantTour::GiantTour(int seed)
    : giant_tour_solution(make_shared <Solution>())
{
    srand(seed);
}

GiantTour::~GiantTour()
{
}

bool GiantTour::unassignedCustomerExists()
{
    for(auto& customer : this->giant_tour_solution->customers) {
        if(!customer->isRouted) {
            return true;
        }
    }
    return false;
}

shared_ptr <Solution> GiantTour::run(const Instance& instance, int alpha /*= 3*/)
{
}

GiantTour_RNN::GiantTour_RNN()
{
}

GiantTour_RNN::GiantTour_RNN(int seed)
    : GiantTour(seed)
{
}


shared_ptr <Solution> GiantTour_RNN::run(const Instance& instance, int alpha /*= 3*/)
{
    int n_vehicles = 1;

    shared_ptr <Solution> solution = make_shared <Solution>(instance, n_vehicles);
    this->giant_tour_solution = solution;

    int vehicle_id = 0;
    int custId;
    Customer* candidate;
    bool isCandidateNull; // NULL value doesn't exists in c++
    // double minCost;
    // vector<int> RCL; // Complexity
    // vector<float> costList; // Complexity
    float candCost;
    // float c_min;
    // float c_max;
    // float endCost;
    int counter = 0;

    solution->vehicles[vehicle_id]->add_customer(solution->customers[0]);
    solution->customers[0]->isRouted = true;
    int randomIndex = (rand() % instance.n_cust);
    if(randomIndex == 0) {
        randomIndex = 1;
    }
    solution->vehicles[vehicle_id]->add_customer(solution->customers[randomIndex]);
    // cout << "before add rand customers Giant tour" << endl;
    solution->customers[randomIndex]->isRouted = true;
	solution->customers[randomIndex]->vehicle_route = vehicle_id; 

    int customers_without_route = instance.n_cust - 2;

    // while (this->unassignedCustomerExists()){

    while(customers_without_route > 0) {
        int number_of_nn = min(alpha, customers_without_route);
        auto current_c = solution->vehicles[vehicle_id]->customers.back();
        vector<unsigned int> neighbors = instance.sorted_dist_matrix[current_c->id];
        vector<unsigned int> near_neighbors(number_of_nn, 999999);
        int index = 1;
        int neighbors_index = 0;
        while(near_neighbors[number_of_nn - 1] == 999999) {
            if(!solution->customers[neighbors[index]]->isRouted) {
                near_neighbors[neighbors_index] = neighbors[index];
                neighbors_index++;
            }
            index++;
        }
        int randomNext = rand() % near_neighbors.size();
        custId = near_neighbors[randomNext];
        solution->vehicles[vehicle_id]->add_customer(solution->customers[custId]);
        solution->customers[custId]->isRouted = true;
		solution->customers[custId]->vehicle_route = vehicle_id;

        customers_without_route--;
    }
    solution->vehicles[vehicle_id]->add_customer(solution->customers[0]);
    // solution->cost = solution->compute_cost(instance);
    return solution;
}

Improved_GiantTour_RNN::Improved_GiantTour_RNN()
{
}

Improved_GiantTour_RNN::Improved_GiantTour_RNN(int seed)
    : GiantTour(seed)
{
}

shared_ptr <Solution> Improved_GiantTour_RNN::run(const Instance& instance, int alpha /*= 3*/)
{
    int n_vehicles = 1;

    shared_ptr <Solution> solution = make_shared <Solution>(instance, n_vehicles);
    this->giant_tour_solution = solution;


    int vehicle_id = 0;
    int custId;
    Customer* candidate;
    bool isCandidateNull; // NULL value doesn't exists in c++
    // double minCost;
    // vector<int> RCL; // Complexity
    // vector<float> costList; // Complexity
    float candCost;
    // float c_min;
    // float c_max;
    // float endCost;
    int counter = 0;

    solution->vehicles[vehicle_id]->add_customer(solution->customers[0]);
    solution->customers[0]->isRouted = true;
    int randomIndex = (rand() % instance.n_cust);
    if(randomIndex == 0) {
        randomIndex = 1;
    }
    solution->vehicles[vehicle_id]->add_customer(solution->customers[randomIndex]);
    // cout << "before add rand customers Giant tour" << endl;
    solution->customers[randomIndex]->isRouted = true;
	solution->customers[randomIndex]->vehicle_route = vehicle_id;

    int customers_without_route = instance.n_cust - 2;

    // while (this->unassignedCustomerExists()){

    while(customers_without_route > 0) {
        int number_of_nn = min(alpha, customers_without_route);
        auto current_c = solution->vehicles[vehicle_id]->customers.back();
        vector<unsigned int> neighbors = instance.sorted_dist_matrix[current_c->id];
        vector<unsigned int> near_neighbors(number_of_nn, 999999);
        int index = 1;
        int neighbors_index = 0;
        while(near_neighbors[number_of_nn - 1] == 999999) {
            if(!solution->customers[neighbors[index]]->isRouted) {
                near_neighbors[neighbors_index] = neighbors[index];
                neighbors_index++;
            }
            index++;
        }
        if(customers_without_route > instance.n_cust / 2) {
        //if(customers_without_route > 0) {
            int randomNext = rand() % near_neighbors.size();
            custId = near_neighbors[randomNext];
            solution->vehicles[vehicle_id]->add_customer(solution->customers[custId]);
            solution->customers[custId]->isRouted = true;
			solution->customers[custId]->vehicle_route = vehicle_id;
        } else {
            custId = near_neighbors[0];
            solution->vehicles[vehicle_id]->add_customer(solution->customers[custId]);
            solution->customers[custId]->isRouted = true;
			solution->customers[custId]->vehicle_route = vehicle_id;
        }
        customers_without_route--;
    }
    solution->vehicles[vehicle_id]->add_customer(solution->customers[0]);

    return solution;
}

GiantTour_RNI::GiantTour_RNI()
{
}

GiantTour_RNI::GiantTour_RNI(int seed)
    : GiantTour(seed)
{
}


shared_ptr <Solution> GiantTour_RNI::run(const Instance& instance, int alpha /*= 3*/)
{

    int n_vehicles = 1;

    shared_ptr <Solution> solution = make_shared<Solution>(instance, n_vehicles);

    this->giant_tour_solution = solution;

    int vehicle_id = 0;
    int custId;
    Customer a;
    Customer b;

    int nb_customers = instance.n_cust;
    vector<int> out_clients(nb_customers, NULL);
    vector<int> in_clients;
    vector<int> sorted_out_clients;
    for(int i = 0; i < nb_customers; i++) {
        out_clients[i] = i;
    }

    vector<int>::iterator position;

    solution->vehicles[vehicle_id]->add_customer(solution->customers[0]);

    solution->customers[0]->isRouted = true;
    int randomIndex = (rand() % instance.n_cust);
    if(randomIndex == 0) {
        randomIndex = 1;
    }
    solution->vehicles[vehicle_id]->add_customer(solution->customers[randomIndex]);
    // cout << "before add rand customers Giant tour" << endl;
    solution->customers[randomIndex]->isRouted = true;
	solution->customers[randomIndex]->vehicle_route = vehicle_id;
    // Round trip depot-1st customer-depot
    solution->vehicles[vehicle_id]->add_customer(solution->customers[0]);

    position = find(out_clients.begin(), out_clients.end(), 0);
    if(position != out_clients.end()) // == myVector.end() means the element was not found
        out_clients.erase(position);
    position = find(out_clients.begin(), out_clients.end(), randomIndex);
    if(position != out_clients.end()) // == myVector.end() means the element was not found
        out_clients.erase(position);

    in_clients.push_back(0);
    in_clients.push_back(randomIndex);

    int customers_without_route = instance.n_cust - 2;

    while(customers_without_route > 0) {
        int number_of_nn = min(alpha, customers_without_route);
        while(sorted_out_clients.size() < number_of_nn) {
            int distance = 0;
            int min_distance = 999999;
            int nearest_customer = 0;
            for(unsigned int i = 0; i < out_clients.size(); i++) {
                if(find(sorted_out_clients.begin(), sorted_out_clients.end(), out_clients[i]) !=
                    sorted_out_clients.end()) {
                    // do nothing
                } else {
                    for(unsigned int j = 0; j < in_clients.size(); j++) {
                        distance = instance.dist_matrix[out_clients[i]][in_clients[j]];
                        if(distance < min_distance) {
                            min_distance = distance;
                            nearest_customer = out_clients[i];
                        }
                    }
                }
            }
            sorted_out_clients.push_back(nearest_customer);
        }
        int randomNext = 0;
        if(customers_without_route > instance.n_cust / 2) {
        //if(customers_without_route > 0) {
            randomNext = rand() % sorted_out_clients.size();
        } else {
            randomNext = 0;
        }
        custId = sorted_out_clients[randomNext];
        int min_delta = 9999999;
        int delta = 0;
        int best_position = 0;
        for(int i = 0; i < instance.n_cust - customers_without_route; i++) {
            a = *solution->vehicles[vehicle_id]->customers[i];
            b = *solution->vehicles[vehicle_id]->customers[i + 1];
            delta = -instance.dist_matrix[a.id][b.id] + instance.dist_matrix[a.id][custId] + instance.dist_matrix[custId][b.id];
            if(delta < min_delta) {
                min_delta = delta;
                best_position = i + 1;
            }
        }
        solution->vehicles[vehicle_id]->insert_customer(solution->customers[custId], best_position);
        solution->customers[custId]->isRouted = true;
		solution->customers[custId]->vehicle_route = vehicle_id;
        in_clients.push_back(custId);
        position = find(out_clients.begin(), out_clients.end(), custId);
        if(position != out_clients.end()) // == myVector.end() means the element was not found
            out_clients.erase(position);
        sorted_out_clients.clear();
        customers_without_route--;
    }
    return solution;
}

GiantTour_RBI::GiantTour_RBI()
{
}

GiantTour_RBI::GiantTour_RBI(int seed)
    : GiantTour(seed)
{
}


shared_ptr <Solution> GiantTour_RBI::run(const Instance& instance, int alpha /*= 3*/)
{

    int n_vehicles = 1;

    shared_ptr <Solution> solution = make_shared<Solution>(instance, n_vehicles);
    this->giant_tour_solution = solution;

    int vehicle_id = 0;
    int custId;
    Customer a;
    Customer b;

    int nb_customers = instance.n_cust;
    vector<int> out_clients(nb_customers, NULL);
    vector<int> in_clients;
    vector<int> sorted_out_clients;
    for(int i = 0; i < nb_customers; i++) {
        out_clients[i] = i;
    }

    vector<int>::iterator position;

    solution->vehicles[vehicle_id]->add_customer(solution->customers[0]);

    solution->customers[0]->isRouted = true;
    int randomIndex = (rand() % instance.n_cust);
    if(randomIndex == 0) {
        randomIndex = 1;
    }
    solution->vehicles[vehicle_id]->add_customer(solution->customers[randomIndex]);
    // cout << "before add rand customers Giant tour" << endl;
    solution->customers[randomIndex]->isRouted = true;
	solution->customers[randomIndex]->vehicle_route = vehicle_id;
    // Round trip depot-1st customer-depot
    solution->vehicles[vehicle_id]->add_customer(solution->customers[0]);

    position = find(out_clients.begin(), out_clients.end(), 0);
    if(position != out_clients.end()) // == myVector.end() means the element was not found
        out_clients.erase(position);
    position = find(out_clients.begin(), out_clients.end(), randomIndex);
    if(position != out_clients.end()) // == myVector.end() means the element was not found
        out_clients.erase(position);

    in_clients.push_back(0);
    in_clients.push_back(randomIndex);

    int customers_without_route = instance.n_cust - 2;

    while(customers_without_route > 0) {
        int number_of_nn = min(alpha, customers_without_route);
        while(sorted_out_clients.size() < number_of_nn) {
            int delta_in = 0;
            int min_delta_in = 999999;
            int nearest_customer = 0;
            for(unsigned int i = 0; i < out_clients.size(); i++) {
                if(find(sorted_out_clients.begin(), sorted_out_clients.end(), out_clients[i]) !=
                    sorted_out_clients.end()) {
                    // do nothing
                } else {
                    for(unsigned int j = 0; j < in_clients.size() - 1; j++) {
                        //delta_in = instance.insertion_matrix[in_clients[j]][in_clients[j + 1]][out_clients[i]];
                        delta_in = -instance.dist_matrix[in_clients[j]][in_clients[j + 1]] + instance.dist_matrix[in_clients[j]][out_clients[i]] + instance.dist_matrix[out_clients[i]][in_clients[j + 1]];
                        if(delta_in < min_delta_in) {
                            min_delta_in = delta_in;
                            nearest_customer = out_clients[i];
                        }
                    }
                }
            }
            sorted_out_clients.push_back(nearest_customer);
        }
        int randomNext = 0;
        if(customers_without_route > instance.n_cust / 2) {
        //if(customers_without_route > 0) {
            randomNext = rand() % sorted_out_clients.size();
        } else {
            randomNext = 0;
        }
        custId = sorted_out_clients[randomNext];
        int min_delta = 9999999;
        int delta = 0;
        int best_position = 0;
        for(int i = 0; i < instance.n_cust - customers_without_route; i++) {
            a = *solution->vehicles[vehicle_id]->customers[i];
            b = *solution->vehicles[vehicle_id]->customers[i + 1];
            //delta = instance.insertion_matrix[a.id][b.id][custId];
            delta = -instance.dist_matrix[a.id][b.id] + instance.dist_matrix[a.id][custId] + instance.dist_matrix[custId][b.id];
            if(delta < min_delta) {
                min_delta = delta;
                best_position = i + 1;
            }
        }
        solution->vehicles[vehicle_id]->insert_customer(solution->customers[custId], best_position);
        solution->customers[custId]->isRouted = true;
		solution->customers[custId]->vehicle_route = vehicle_id;
        in_clients.push_back(custId);
        position = find(out_clients.begin(), out_clients.end(), custId);
        if(position != out_clients.end()) // == myVector.end() means the element was not found
            out_clients.erase(position);
        sorted_out_clients.clear();
        customers_without_route--;
    }
    return solution;
}