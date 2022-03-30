#include "./cw_savings.h"
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

using namespace std;

ParallelSavings::ParallelSavings()
{
}

ParallelSavings::ParallelSavings(int seed, int beta, vector<vector<int>> savings)
{
    srand(seed);
}


ParallelSavings::~ParallelSavings()
{
}


vector<vector<int>> clarke_wright_savings_function(vector<vector<float>> dist_matrix)
{
    vector<vector<int>> savings;
    int s = 0;
    for(int i = 1; i < dist_matrix.size(); i++) {
        for(int j = i + 1; j < dist_matrix.size(); j++) {
            s = dist_matrix[i][0] + dist_matrix[0][j] - dist_matrix[i][j];
            savings.push_back({ s, static_cast<int> (-dist_matrix[i][j]), i, j });
        }
    }
    std::sort(savings.begin(), savings.end(),
        [](const std::vector<int>& a, const std::vector<int>& b) { return a[0] > b[0]; });
    return savings;
}
vector<vector<int>> complete_randomized_clarke_wright_savings_function(vector<vector<float>> dist_matrix, float mean_demands, vector<int> demands, bool use_new_savings, int beta)
{
    vector<vector<int>> savings;
    int s = 0;
    float lambda = 1.3;
    float mu = 0.6;
    float v = 0.6;
    for(int i = 1; i < dist_matrix.size(); i++) {
        for(int j = i + 1; j < dist_matrix.size(); j++) {
            //s = dist_matrix[i][0] + dist_matrix[0][j] -  dist_matrix[i][j]; //Original savings calculation
            if (use_new_savings == true) {
                s = dist_matrix[i][0] + dist_matrix[0][j] - (lambda * dist_matrix[i][j]) \
 + (mu * abs(dist_matrix[i][0] - dist_matrix[0][j])) + (v * (demands[i] + demands[j]) / mean_demands);
            }
            else{
                s = dist_matrix[i][0] + dist_matrix[0][j] -  dist_matrix[i][j]; //Original savings calculation
            }
            savings.push_back({ s, static_cast<int> (-dist_matrix[i][j]), i, j });
        }
    }
    std::sort(savings.begin(), savings.end(),
        [](const std::vector<int>& a, const std::vector<int>& b) { return a[0] > b[0]; });
        
    vector<vector<int>> random_savings;
    vector<vector<int>> not_included_list;
    bool TrueFalse = (rand() % 100) < beta;
    for (int i = 0; i < savings.size() ; i++)
    {
        TrueFalse = (rand() % 100) < beta;
        if (TrueFalse){
            random_savings.push_back(savings[i]);
        }
        else {
            not_included_list.push_back(savings[i]);
        }
    }
    for (int j = 0; j < not_included_list.size(); j++){
        random_savings.push_back(not_included_list[j]);
    }
    
    return random_savings;
}
vector<vector<int>> randomized_clarke_wright_savings_function(vector<vector<int>> savings, int beta)
{
    vector<vector<int>> random_savings;
    vector<vector<int>> not_included_list;
    bool TrueFalse = (rand() % 100) < beta;
    for (int i = 0; i < savings.size() ; i++)
    {
        TrueFalse = (rand() % 100) < beta;
        if (TrueFalse){
            random_savings.push_back(savings[i]);
        }
        else {
            not_included_list.push_back(savings[i]);
        }
    }
    for (int j = 0; j < not_included_list.size(); j++){
        random_savings.push_back(not_included_list[j]);
    }
    
    return random_savings;
}

//Solution*
shared_ptr <Solution> ParallelSavings::run(shared_ptr <Solution> solution, const Instance& instance, int beta, vector<vector<int>> deterministic_savings)
{
    int n_vehicles = instance.n_cust;

    //Solution* solution = new Solution(instance, n_vehicles);
    //cout << "prlls address: "<< solution << endl;
    //solution->solution_address = solution;
    //this->created_solutions.push_back(solution);

    //this->savings_solution = *solution;
    //this->created_solutions.push_back(&savings_solution);
    //
    //cout<< "Before filling vehicles" << endl;
    for(int i = 1; i < solution->customers.size(); i++) {
        //solution->vehicles[i]->add_customer(solution->customers[0]);
        solution->customers[i]->vehicle_route = i;
    }
    solution->customers[0]->isRouted = true;
    solution->customers[0]->vehicle_route = 0;
    
    //cout<< "Before savings" << endl;
    vector<vector<int>> savings = randomized_clarke_wright_savings_function(deterministic_savings, beta);
    //cout<< "After savings" << endl;
    vector<int> top_saving;
    //Customer* a;
    //Customer* b;
    shared_ptr <Customer> a;
    shared_ptr <Customer> b;
    int current_vehicle = 0;
    int customer_a_pos = 0;
    int customer_b_pos = 0;
    int other_vehicle = 0;
    bool join_a_left = false;
    bool join_a_right = false;
    bool join_b_left = false;
    bool join_b_right = false;

    for(int i = 0; i < savings.size(); i++) {
        //cout<< "Inside savings loop" << endl;
        top_saving = savings[i];
        a = solution->customers[top_saving[2]];
        b = solution->customers[top_saving[3]];
        join_a_left = false;
        join_a_right = false;
        join_b_left = false;
        join_b_right = false;

        /*
        if (solution->vehicles[current_vehicle]->load + a->demand + b->demand > instance->Q)
        {
            continue;
        }
         */

        if((a->isRouted == false) && (b->isRouted == false)) {
            //cout<< "entered 1" << endl;
            current_vehicle = a->vehicle_route;
            //current_vehicle = solution->vehicles.size();
            if(solution->vehicles[current_vehicle]->load + a->demand + b->demand > instance.Q) {
                continue;
            }
            solution->vehicles[current_vehicle]->add_customer(solution->customers[0]);
            solution->vehicles[current_vehicle]->add_customer(a);
            solution->vehicles[current_vehicle]->add_customer(b);
            solution->vehicles[current_vehicle]->add_customer(solution->customers[0]);
            solution->customers[top_saving[2]]->vehicle_route = current_vehicle;
            solution->customers[top_saving[3]]->vehicle_route = current_vehicle;
            solution->customers[top_saving[2]]->isRouted = true;
            solution->customers[top_saving[3]]->isRouted = true;
            //cout<< "exit 1" << endl;
        }

        else if((a->isRouted == true) && (b->isRouted == false)) {
            //cout<< "entered 2" << endl;
            current_vehicle = a->vehicle_route;
            // Check demand constraint
            if(solution->vehicles[current_vehicle]->load + b->demand > instance.Q) {
                //cout<< "excess load 2" << endl;
                continue;
            }
            for(int j = 1; j < solution->vehicles[current_vehicle]->customers.size() - 1; j++) {
                if(a->id == solution->vehicles[current_vehicle]->customers[j]->id) {
                    customer_b_pos = j;
                    break;
                }
            }
            // Check if interior route
            if(solution->vehicles[current_vehicle]->customers[customer_b_pos + 1]->id == 0) {
                customer_b_pos = customer_b_pos + 1;
            } else if(solution->vehicles[current_vehicle]->customers[customer_b_pos - 1]->id == 0) {
                customer_b_pos = customer_b_pos;
            } else {
                //cout<< "Interior route 2" << endl;
                continue;
            }

            solution->vehicles[current_vehicle]->insert_customer(b, customer_b_pos);
            solution->vehicles[current_vehicle]->customers[customer_b_pos]->vehicle_route = current_vehicle;
            solution->customers[top_saving[3]]->isRouted = true;
            //cout<< "exit 2" << endl;
        }

        else if((a->isRouted == false) && (b->isRouted == true)) {
            //cout<< "entered 3" << endl;
            current_vehicle = b->vehicle_route;
            // Check demand constraint
            if(solution->vehicles[current_vehicle]->load + a->demand > instance.Q) {
                continue;
            }
            for(int j = 1; j < solution->vehicles[current_vehicle]->customers.size() - 1; j++) {
                if(b->id == solution->vehicles[current_vehicle]->customers[j]->id) {
                    customer_a_pos = j;
                    break;
                }
            }
            // Check if interior route
            if(solution->vehicles[current_vehicle]->customers[customer_a_pos - 1]->id == 0) {
                customer_a_pos = customer_a_pos;
            } else if(solution->vehicles[current_vehicle]->customers[customer_a_pos + 1]->id == 0) {
                customer_a_pos = customer_a_pos + 1;
            } else {
                continue;
            }

            solution->vehicles[current_vehicle]->insert_customer(a, customer_a_pos);
            solution->vehicles[current_vehicle]->customers[customer_a_pos]->vehicle_route = current_vehicle;
            solution->customers[top_saving[2]]->isRouted = true;
            //cout<< "exit 3" << endl;
        }
        
        else if((a->isRouted == true) && (b->isRouted == true)) {
            //cout<< "entered 4" << endl;
            current_vehicle = a->vehicle_route;
            other_vehicle = b->vehicle_route;
            
            //cout<<current_vehicle << " current vehicle size: " << solution->vehicles[current_vehicle]->customers.size() << endl;
            //cout<<other_vehicle << " other vehicle size: " << solution->vehicles[other_vehicle]->customers.size() << endl;
            
            if (current_vehicle == other_vehicle){
                //cout<<"same vehicle" << endl;
                continue;
            }
            // Check demand constraint
            if(solution->vehicles[current_vehicle]->load + solution->vehicles[other_vehicle]->load > instance.Q) {
                //cout<< "excess load 4" << endl;
                continue;
            }
            
            for(int j = 1; j < solution->vehicles[current_vehicle]->customers.size() - 1; j++) {
                if(a->id == solution->vehicles[current_vehicle]->customers[j]->id) {
                    customer_b_pos = j;
                    break;
                }
            }
            //cout<< customer_b_pos << endl;
            // Check if a is interior route
            if(solution->vehicles[current_vehicle]->customers[customer_b_pos + 1]->id == 0) {
                customer_b_pos = customer_b_pos + 1;
                join_a_right = true;
            } else if(solution->vehicles[current_vehicle]->customers[customer_b_pos - 1]->id == 0) {
                customer_b_pos = customer_b_pos;
                join_a_left = true;
            } else {
                //cout<< "a is interior " << endl;
                continue;
            }
            
            
            for(int j = 1; j < solution->vehicles[other_vehicle]->customers.size() - 1; j++) {
                if(b->id == solution->vehicles[other_vehicle]->customers[j]->id) {
                    customer_a_pos = j;
                    break;
                }
            }
            //cout<< customer_a_pos << endl;
            // Check if b is interior route
            if(solution->vehicles[other_vehicle]->customers[customer_a_pos - 1]->id == 0) {
                customer_a_pos = customer_a_pos;
                join_b_left = true;
            } else if(solution->vehicles[other_vehicle]->customers[customer_a_pos + 1]->id == 0) {
                customer_a_pos = customer_a_pos + 1;
                join_b_right = true;
            } else {
                //cout<< "b is interior " << endl;
                continue;
            }
            
            //Merge not interior routes
            if ((join_a_left == true) && (join_b_left==true)){
                //cout<< "left-left" << endl;
                int cont = 1;
                while (cont < solution->vehicles[other_vehicle]->customers.size() -1){
                    solution->vehicles[current_vehicle]->insert_customer(solution->vehicles[other_vehicle]->customers[1], 1);
                    solution->vehicles[other_vehicle]->customers[1]->vehicle_route = current_vehicle;
                    solution->vehicles[other_vehicle]->temporary_remove_customer(1);
                    cont++;
                }
            }
            
            else if ((join_a_left == true) && (join_b_right==true)){
                //cout<< "left-right" << endl;
                //cout<< current_vehicle << " size current: " << solution->vehicles[current_vehicle]->customers.size() << endl;
                //cout<< other_vehicle << " size other: " << solution->vehicles[other_vehicle]->customers.size() << endl;
                int cont = 1;
                while (cont < solution->vehicles[current_vehicle]->customers.size() -1){
                    auto customer_one = solution->vehicles[current_vehicle]->customers[1];
                    solution->vehicles[other_vehicle]->insert_customer(customer_one, solution->vehicles[other_vehicle]->customers.size() - 1);
                    solution->vehicles[current_vehicle]->customers[1]->vehicle_route = other_vehicle;
                    solution->vehicles[current_vehicle]->temporary_remove_customer(1);
                    cont++;
                }
            }
            
            else if ((join_a_right == true) && (join_b_left==true)){
                //cout<< "right-left" << endl;
                int cont = 1;
                 while ( cont < solution->vehicles[other_vehicle]->customers.size() -1){
                     solution->vehicles[current_vehicle]->insert_customer(solution->vehicles[other_vehicle]->customers[1], solution->vehicles[current_vehicle]->customers.size() - 1);
                    solution->vehicles[other_vehicle]->customers[1]->vehicle_route = current_vehicle;
                    solution->vehicles[other_vehicle]->temporary_remove_customer(1);
                    cont++;
                 }

            }
            
            else if ((join_a_right == true) && (join_b_right==true)){
                //cout<< "right-right" << endl;
                int cont = 1;
                while(cont <solution->vehicles[other_vehicle]->customers.size() -1){
                    solution->vehicles[current_vehicle]->insert_customer(solution->vehicles[other_vehicle]->customers[1], solution->vehicles[current_vehicle]->customers.size() - cont);
                    solution->vehicles[other_vehicle]->customers[1]->vehicle_route = current_vehicle;
                    solution->vehicles[other_vehicle]->temporary_remove_customer(1);
                    cont++;
                }
            }
            
            //cout<< "exit  4" << endl;
        }
    }
    solution->delete_empty_vehicles();
    //solution->cost = solution->compute_cost(instance);



    for (int k=1; k < solution->customers.size(); k++){
        if (solution->customers[k]->isRouted == false){
            if (solution->vehicles[solution->vehicles.size()-1]->load + solution->customers[k]->demand < instance.Q){
                solution->vehicles[solution->vehicles.size() - 1]->customers.pop_back();
                solution->vehicles[solution->vehicles.size() - 1]->add_customer(solution->customers[k]);
                solution->vehicles[solution->vehicles.size() - 1]->add_customer(solution->customers[0]);
                solution->customers[k]->vehicle_route = solution->vehicles.size() - 1;
                solution->customers[k]->isRouted = true;
            }
            else {
                //cout << "Unrouted customer: " << k <<endl;
                //solution->vehicles.push_back(new Vehicle(instance.Q));
                solution->vehicles.push_back(make_unique<Vehicle>(instance.Q));
                solution->vehicles[solution->vehicles.size() - 1]->add_customer(solution->customers[0]);
                solution->vehicles[solution->vehicles.size() - 1]->add_customer(solution->customers[k]);
                solution->vehicles[solution->vehicles.size() - 1]->add_customer(solution->customers[0]);
                solution->customers[k]->vehicle_route = solution->vehicles.size() - 1;
                solution->customers[k]->isRouted = true;
            }
        }
    }
    
    solution->cost = 0;
    for (int i=0; i < solution->vehicles.size(); i++){
        solution->vehicles[i]->compute_cost(instance);
        solution->cost += solution->vehicles[i]->cost;
    }
    //cout << "Cost with vehicles: " << solution->cost;
    
    //solution->cost = solution->compute_cost(instance);
    //cout << "Cost with old method: " << solution->cost;
    return solution;
}

shared_ptr <Solution> ParallelSavings::clarke_and_wright(const Instance& instance, shared_ptr <Solution> solution, float lambda, int neighbors_num, int beta, mt19937& rand_engine) {
    int rcl_list_size = 10;
    std::uniform_int_distribution<int> rand_uniform(0, rcl_list_size-1);

    for (auto i = instance.customers_begin; i < instance.customers_end; i++) {
        //solution.build_one_customer_route(i);
        //solution->vehicles.push_back(make_unique<Vehicle>(instance.Q));
        solution->vehicles[i-1]->add_customer(solution->customers[0]);
        solution->vehicles[i-1]->add_customer(solution->customers[i]);
        solution->vehicles[i-1]->add_customer(solution->customers[0]);
        solution->customers[i]->vehicle_route = i-1;
        solution->customers[i]->isRouted = true;
    }

    neighbors_num = std::min(instance.customers_num - 1, neighbors_num);

    const auto savings_num = instance.customers_num * neighbors_num;

    struct Saving {
        int i;
        int j;
        float value;
    };

    auto savings = std::vector<Saving>();
    savings.reserve(static_cast<unsigned long>(savings_num));

    for(int i = instance.customers_begin; i < instance.customers_end; i++) {

        for(auto n = 1u, added=0u; added < static_cast<unsigned int>(neighbors_num) && n < instance.sorted_dist_matrix[i].size(); n++) {

            const int j = instance.sorted_dist_matrix[i][n];

            if (i < j) {

                const auto value = +instance.get_cost(i, instance.depot) + instance.get_cost(instance.depot, j) - lambda*instance.get_cost(i, j);

                savings.push_back({i, j, value});

                added++;

            }
        }
    }


    std::sort(savings.begin(), savings.end(), [](const Saving &a, const Saving &b) { return a.value > b.value; });

    //Randomize a bit the savings order
    auto random_savings = std::vector<Saving>();
    auto not_included_list = std::vector<Saving>();
    random_savings.reserve(static_cast<unsigned long>(savings_num));

    /*
    bool TrueFalse = (rand() % 100) < beta;
    for (int i = 0; i < savings.size() ; i++)
    {
        TrueFalse = (rand() % 100) < beta;
        if (TrueFalse){
            random_savings.push_back(savings[i]);
        }
        else {
            not_included_list.push_back(savings[i]);
        }
    }
    for (int j = 0; j < not_included_list.size(); j++){
        random_savings.push_back(not_included_list[j]);
    }
    */

    int selected_candidate = 0;
    auto rcl_list = std::vector<Saving>();
    rcl_list.reserve(static_cast<unsigned long>(rcl_list_size));
    for (int i = 0; i < rcl_list_size ; i++)
    {
        rcl_list.push_back(savings[i]);
    }
    for (int i = rcl_list_size; i < savings.size() ; i++)
    {
        selected_candidate = rand_uniform(rand_engine);
        random_savings.push_back(rcl_list[selected_candidate]);
        rcl_list.erase(rcl_list.begin()+selected_candidate);
        rcl_list.push_back(savings[i]);
    }
    for (int i = 0; i < rcl_list.size() ; i++)
    {
        random_savings.push_back(rcl_list[i]);
    }


    for (auto &saving : random_savings) {

        const auto i = saving.i;
        const auto j = saving.j;

        const auto iRoute = solution->get_route_index(i);
        const auto jRoute = solution->get_route_index(j);

        if (iRoute == jRoute) { continue; }

        if (solution->vehicles[iRoute]->customers[solution->vehicles[iRoute]->customers.size()-2]->id == i && solution->vehicles[jRoute]->customers[1]->id == j &&
            solution->vehicles[iRoute]->load + solution->vehicles[jRoute]->load <= instance.Q) {

            solution->append_route(iRoute, jRoute, instance);


        } else if (solution->vehicles[jRoute]->customers[solution->vehicles[jRoute]->customers.size()-2]->id == j && solution->vehicles[iRoute]->customers[1]->id == i &&
                solution->vehicles[iRoute]->load + solution->vehicles[jRoute]->load <= instance.Q) {

            solution->append_route(jRoute, iRoute, instance);

        }
    }
    solution->delete_empty_vehicles();
    //solution->compute_cost(instance);
    //solution->is_feasible();
    return solution;
}

