#include "./neighborhood.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <numeric>
#include <queue>
#include <unordered_map>
#include <utility>

using namespace std;

Neighborhood::Neighborhood(bool first_improvement /*= false*/, const string& id /*= ""*/)
{
    this->first_improvement = first_improvement;
    this->id = id;
}

Neighborhood::~Neighborhood()
{
}

std::shared_ptr <Solution> Neighborhood::search_neighbors(std::shared_ptr <Solution> initial_solution, const Instance& instance)
{
}

SwapNeighborhood::SwapNeighborhood()
{
}

SwapNeighborhood::SwapNeighborhood(bool first_improvement /*= false*/)
        : Neighborhood(first_improvement, "A")
{
}

SwapNeighborhood::~SwapNeighborhood()
{
}

std::shared_ptr <Solution> SwapNeighborhood::search_neighbors(std::shared_ptr <Solution> initial_solution, const Instance& instance)
{
    int n_vehicles = initial_solution->vehicles.size();

    float min_delta = 0.0;
    // float initial_cost = initial_solution->cost;
    vector<int> best_move;
    best_move.clear();
    Vehicle vehicle_one;
    Vehicle vehicle_two;
    int initial_value_for;
    // customers
    Customer a;
    Customer b;
    Customer c;
    Customer d;

    Customer f;
    Customer g;
    Customer h;
    // demands
    int route1_new_demand;
    int route2_new_demand;
    // dist_matrix

    // delta
    float delta = 0.0;

    for(int vehicle_one_id = 0; vehicle_one_id < n_vehicles; vehicle_one_id++) {
        vehicle_one = *initial_solution->vehicles[vehicle_one_id];
        for(int customer_one_id = 1; customer_one_id < vehicle_one.customers.size() - 1; customer_one_id++) {
            for(int vehicle_two_id = vehicle_one_id; vehicle_two_id < n_vehicles; vehicle_two_id++) {
                vehicle_two = *initial_solution->vehicles[vehicle_two_id];
                if(vehicle_one_id != vehicle_two_id) {
                    initial_value_for = 1;
                    // continue; // With these line Swap is now only intraroute
                } else
                    initial_value_for = customer_one_id + 1;

                for(int customer_two_id = initial_value_for; customer_two_id < vehicle_two.customers.size() - 1;
                    customer_two_id++) {
                    /*
                        route1 =  [A B C D]
                        route2 =  [E F G H]
                        Swap(B, G) leaving us with:
                        new_route1 = [A G C D]
                        new_route2 = [E F B H]
                    */

                    a = *vehicle_one.customers[customer_one_id - 1];
                    b = *vehicle_one.customers[customer_one_id];
                    c = *vehicle_one.customers[customer_one_id + 1];

                    f = *vehicle_two.customers[customer_two_id - 1];
                    g = *vehicle_two.customers[customer_two_id];
                    h = *vehicle_two.customers[customer_two_id + 1];

                    // Check constraint of load violation
                    route1_new_demand = vehicle_one.load - b.demand + g.demand;
                    route2_new_demand = vehicle_two.load - g.demand + b.demand;

                    if((route1_new_demand > vehicle_one.Q) || (route2_new_demand > vehicle_two.Q) || (b.id == g.id)) {
                        continue;
                    }

                    //   Evaluate the improvement delta of both route1 and route2, using distance matrix
                    //   Delta = D(a,g) + D(g,c) + D(f,b) + D(b,h) - D(a,b) -D(b,c) - D(f,g) - D(g,h)
                    // delta = +dist[a.id][g.id] + dist[g.id][c.id] + dist[f.id][b.id] + dist[b.id][h.id]
                    // -dist[a.id][b.id] - dist[b.id][c.id] - dist[f.id][g.id] - dist[g.id][h.id];
                    /*
                     delta = instance.insertion_matrix[f.id][h.id][b.id] +
                                instance.insertion_matrix[a.id][c.id][g.id] + instance.remotion_matrix[f.id][h.id][g.id] +
                                instance.remotion_matrix[a.id][c.id][b.id];
                     */

                    if((vehicle_one_id == vehicle_two_id) && (abs(customer_one_id - customer_two_id) == 1)) {
                        d = *vehicle_one.customers[customer_one_id + 2];
                        delta = +instance.dist_matrix[a.id][g.id] + instance.dist_matrix[b.id][d.id] - instance.dist_matrix[a.id][b.id] - instance.dist_matrix[g.id][d.id];
                    } else if((vehicle_one_id == vehicle_two_id) && (abs(customer_one_id - customer_two_id) == 2)) {
                        // d = *vehicle_one.customers[customer_one_id + 2];
                        delta = +instance.dist_matrix[a.id][g.id] + instance.dist_matrix[b.id][h.id] - instance.dist_matrix[a.id][b.id] - instance.dist_matrix[g.id][h.id];
                    } else {
                        delta = +instance.dist_matrix[a.id][g.id] + instance.dist_matrix[g.id][c.id] + instance.dist_matrix[f.id][b.id] + instance.dist_matrix[b.id][h.id]
                         -instance.dist_matrix[a.id][b.id] - instance.dist_matrix[b.id][c.id] - instance.dist_matrix[f.id][g.id] - instance.dist_matrix[g.id][h.id];
                    }

                    // Store the current best possible move
                    if(delta < min_delta) {
                        min_delta = delta;
                        best_move.clear();
                        best_move.push_back(vehicle_one_id);
                        best_move.push_back(vehicle_two_id);
                        best_move.push_back(customer_one_id);
                        best_move.push_back(customer_two_id);
                        // best_move = { vehicle_one_id, vehicle_two_id, customer_one_id, customer_two_id,};
                        if(this->first_improvement) {
                            initial_solution->cost = initial_solution->cost + min_delta;
                            return this->apply_best_move(initial_solution, best_move);
                        }
                    }
                }
            }
        }
    }
    if(min_delta < 0) {
        /*
        cout << "min_delta " << min_delta << " C1 " << best_move[2] << " C2 " << best_move[3] << endl;
        float veh_1_cost = initial_solution->vehicles[best_move[0]]->compute_cost(instance);
        float veh_2_cost = initial_solution->vehicles[best_move[1]]->compute_cost(instance);
        cout << "cost_before " << veh_1_cost << " " << veh_2_cost << endl;
*/
        initial_solution = apply_best_move(initial_solution, best_move);
        initial_solution->cost = initial_solution->cost + min_delta;
        initial_solution->vehicles[best_move[0]]->compute_cost(instance);
        initial_solution->vehicles[best_move[1]]->compute_cost(instance);
        /*
                veh_1_cost = initial_solution->vehicles[best_move[0]]->compute_cost(instance);
                veh_2_cost = initial_solution->vehicles[best_move[1]]->compute_cost(instance);
                cout << "cost_after " << veh_1_cost << " " << veh_2_cost << endl;
                 * */
    }
    return initial_solution;
}

/**
 * Apply the Swap (exchange) to the best possible move in the neighborhood
 * @param initial_solution Current solution of the VRP
 * @param best_move {vehicle_one_id, vehicle_two_id, customer_one_id, customer_two_id}
 *
 * @return local_optima: Neew VRP solution corresponding to the local optima of this neighborhood
 *          if there's no local optima, returns initial_solution instead
 */
std::shared_ptr <Solution> SwapNeighborhood::apply_best_move(std::shared_ptr <Solution> local_optima, vector<int>& best_move)
{
    if(best_move.size() != 0) {
        // std::shared_ptr <Solution> local_optima = initial_solution;
        int vehicle_one_id = best_move[0];
        int vehicle_two_id = best_move[1];
        int customer_one_id = best_move[2];
        int customer_two_id = best_move[3];

        if(vehicle_one_id == vehicle_two_id) {
            shared_ptr<Customer> customer_one = local_optima->vehicles[vehicle_one_id]->customers[customer_one_id];
            shared_ptr<Customer> customer_two = local_optima->vehicles[vehicle_two_id]->customers[customer_two_id];
            if(customer_two_id > customer_one_id) {
                if(customer_two_id == customer_one_id + 1) {
                    local_optima->vehicles[vehicle_one_id]->insert_customer(customer_two, customer_one_id);
                    local_optima->vehicles[vehicle_two_id]->remove_customer(customer_two_id + 1);
                } else {
                    local_optima->vehicles[vehicle_one_id]->insert_customer(customer_two, customer_one_id);
                    local_optima->vehicles[vehicle_two_id]->remove_customer(customer_two_id + 1);

                    local_optima->vehicles[vehicle_two_id]->insert_customer(customer_one, customer_two_id + 1);
                    local_optima->vehicles[vehicle_one_id]->remove_customer(customer_one_id + 1);
                }
            } else {
                if(customer_one_id == customer_two_id + 1) {
                    local_optima->vehicles[vehicle_two_id]->insert_customer(customer_one, customer_two_id);
                    local_optima->vehicles[vehicle_one_id]->remove_customer(customer_one_id + 1);
                } else {
                    local_optima->vehicles[vehicle_two_id]->insert_customer(customer_one, customer_two_id);
                    local_optima->vehicles[vehicle_one_id]->remove_customer(customer_one_id + 1);

                    local_optima->vehicles[vehicle_one_id]->insert_customer(customer_two, customer_one_id + 1);
                    local_optima->vehicles[vehicle_two_id]->remove_customer(customer_two_id + 1);
                }
            }

        } else {

            shared_ptr<Customer> customer_two = local_optima->vehicles[vehicle_two_id]->customers[customer_two_id];

            // local_optima->vehicles[vehicle_one_id]->load -= customer_one->demand;
            // local_optima->vehicles[vehicle_one_id]->load += customer_two->demand;
            // local_optima->vehicles[vehicle_one_id]->customers[customer_one_id] = customer_two;

            local_optima->vehicles[vehicle_one_id]->insert_customer(customer_two, customer_one_id);

            // local_optima->vehicles[vehicle_two_id]->load -= customer_two->demand;
            // local_optima->vehicles[vehicle_two_id]->load += customer_one->demand;
            // local_optima->vehicles[vehicle_two_id]->customers[customer_two_id] = customer_one;

            shared_ptr<Customer> customer_one = local_optima->vehicles[vehicle_one_id]->customers[customer_one_id + 1];
            local_optima->vehicles[vehicle_two_id]->insert_customer(customer_one, customer_two_id);

            local_optima->vehicles[vehicle_one_id]->remove_customer(customer_one_id + 1);
            local_optima->vehicles[vehicle_two_id]->remove_customer(customer_two_id + 1);

            local_optima->vehicles[vehicle_one_id]->customers[customer_one_id]->vehicle_route = vehicle_one_id;
			local_optima->customers[local_optima->vehicles[vehicle_one_id]->customers[customer_one_id]->id]->vehicle_route = vehicle_one_id;
            local_optima->vehicles[vehicle_two_id]->customers[customer_two_id]->vehicle_route = vehicle_two_id;
			local_optima->customers[local_optima->vehicles[vehicle_two_id]->customers[customer_two_id]->id]->vehicle_route = vehicle_two_id;

            // local_cost = local_optima->cost ;
        }
        return local_optima;
    } else {
        return local_optima;
    }
}

SwapTwoNeighborhood::SwapTwoNeighborhood()
{
}

SwapTwoNeighborhood::SwapTwoNeighborhood(bool first_improvement /*= false*/)
        : Neighborhood(first_improvement, "E")
{
}

SwapTwoNeighborhood::~SwapTwoNeighborhood()
{
}

std::shared_ptr <Solution> SwapTwoNeighborhood::search_neighbors(std::shared_ptr <Solution> initial_solution, const Instance& instance)
{
    int n_vehicles = initial_solution->vehicles.size();

    float min_delta = 0;
    vector<int> best_move;
    Vehicle vehicle_one;
    Vehicle vehicle_two;
    int initial_value_for;
    // customers
    Customer a;
    Customer b;
    Customer c;
    Customer d;

    Customer f;
    Customer g;
    Customer h;
    Customer i;
    // demands
    int route1_new_demand;
    int route2_new_demand;


    // delta
    float delta;
    for(int vehicle_one_id = 0; vehicle_one_id < n_vehicles; vehicle_one_id++) {
        vehicle_one = *initial_solution->vehicles[vehicle_one_id];
        for(int customer_one_id = 1; customer_one_id < vehicle_one.customers.size() - 2; customer_one_id++) {
            for(int vehicle_two_id = vehicle_one_id; vehicle_two_id < n_vehicles; vehicle_two_id++) {
                vehicle_two = *initial_solution->vehicles[vehicle_two_id];
                if(vehicle_one_id != vehicle_two_id)
                    initial_value_for = 1;
                else
                    initial_value_for = customer_one_id + 2;
                for(int customer_two_id = initial_value_for; customer_two_id < vehicle_two.customers.size() - 2;
                    customer_two_id++) {
                    /*
                        route1 =  [A B C D]
                        route2 =  [E F G H I]
                        Swap(B C, G H) leaving us with:
                        new_route1 = [A G H D]
                        new_route2 = [E F B C I]
                    */
                    if((vehicle_one.customers.size() <= 3) || (vehicle_one.customers.size() <= 3)) {
                        continue;
                    }

                    a = *vehicle_one.customers[customer_one_id - 1];
                    b = *vehicle_one.customers[customer_one_id];
                    c = *vehicle_one.customers[customer_one_id + 1];
                    d = *vehicle_one.customers[customer_one_id + 2];

                    f = *vehicle_two.customers[customer_two_id - 1];
                    g = *vehicle_two.customers[customer_two_id];
                    h = *vehicle_two.customers[customer_two_id + 1];
                    i = *vehicle_two.customers[customer_two_id + 2];

                    // Check constraint of load violation
                    route1_new_demand = vehicle_one.load - b.demand - c.demand + g.demand + h.demand;
                    route2_new_demand = vehicle_two.load - g.demand - h.demand + b.demand + c.demand;

                    if(route1_new_demand > vehicle_one.Q || route2_new_demand > vehicle_two.Q) {
                        continue;
                    }

                    /*
                        Evaluate the improvement delta of both route1 and route2, using distance matrix
                        Delta = D(a,g) + D(g,h) +D(h,d) + D(f,b) + D(b,c) +D(c,i) - D(a,b) -D(b,c) -D(c,d) -D(f,g)
                       -D(g,h) -D(h,i)
                    */

                    delta = instance.dist_matrix[a.id][g.id] + instance.dist_matrix[g.id][h.id] + instance.dist_matrix[h.id][d.id] + instance.dist_matrix[f.id][b.id] +
                            instance.dist_matrix[b.id][c.id] + instance.dist_matrix[c.id][i.id] - instance.dist_matrix[a.id][b.id] - instance.dist_matrix[b.id][c.id] - instance.dist_matrix[c.id][d.id] -
                            instance.dist_matrix[f.id][g.id] - instance.dist_matrix[g.id][h.id] - instance.dist_matrix[h.id][i.id];

                    // Store the current best possible move
                    if(delta < min_delta) {
                        min_delta = delta;
                        best_move = { vehicle_one_id, vehicle_two_id, customer_one_id, customer_two_id,
                                      customer_one_id + 1, customer_two_id + 1 };
                        if(this->first_improvement) {
                            return this->apply_best_move(initial_solution, best_move);
                        }
                    }
                }
            }
        }
    }

    return this->apply_best_move(initial_solution, best_move);
}

/**
 * Apply the Swap(2,2) (exchange) to the best possible move in the neighborhood
 * @param initial_solution Current solution of the VRP
 * @param best_move {vehicle_one_id, vehicle_two_id, customer_one_id, customer_two_id, customer_one_id + 1,
 * customer_two_id + 1}
 *
 * @return local_optima: Neew VRP solution corresponding to the local optima of this neighborhood
 *          if there's no local optima, returns initial_solution instead
 */
std::shared_ptr <Solution> SwapTwoNeighborhood::apply_best_move(std::shared_ptr <Solution> initial_solution, vector<int>& best_move)
{
    if(best_move.size() != 0) {
        std::shared_ptr <Solution> local_optima = initial_solution;
        int vehicle_one_id = best_move[0];
        int vehicle_two_id = best_move[1];
        int customer_one_id = best_move[2];
        int customer_two_id = best_move[3];
        int customer_three_id = best_move[4];
        int customer_four_id = best_move[5];

        shared_ptr<Customer> customer_one = local_optima->vehicles[vehicle_one_id]->customers[customer_one_id];
        shared_ptr<Customer> customer_two = local_optima->vehicles[vehicle_two_id]->customers[customer_two_id];
        shared_ptr<Customer> customer_three = local_optima->vehicles[vehicle_one_id]->customers[customer_three_id];
        shared_ptr<Customer> customer_four = local_optima->vehicles[vehicle_two_id]->customers[customer_four_id];

        local_optima->vehicles[vehicle_one_id]->load -= customer_three->demand;
        local_optima->vehicles[vehicle_one_id]->load += customer_four->demand;
        local_optima->vehicles[vehicle_one_id]->customers[customer_three_id] = customer_four;

        local_optima->vehicles[vehicle_two_id]->load -= customer_four->demand;
        local_optima->vehicles[vehicle_two_id]->load += customer_three->demand;
        local_optima->vehicles[vehicle_two_id]->customers[customer_four_id] = customer_three;

        local_optima->vehicles[vehicle_one_id]->load -= customer_one->demand;
        local_optima->vehicles[vehicle_one_id]->load += customer_two->demand;
        local_optima->vehicles[vehicle_one_id]->customers[customer_one_id] = customer_two;

        local_optima->vehicles[vehicle_two_id]->load -= customer_two->demand;
        local_optima->vehicles[vehicle_two_id]->load += customer_one->demand;
        local_optima->vehicles[vehicle_two_id]->customers[customer_two_id] = customer_one;

        local_optima->vehicles[vehicle_one_id]->customers[customer_three_id]->vehicle_route = vehicle_one_id;
		local_optima->customers[local_optima->vehicles[vehicle_one_id]->customers[customer_three_id]->id]->vehicle_route = vehicle_one_id;
        local_optima->vehicles[vehicle_two_id]->customers[customer_four_id]->vehicle_route = vehicle_two_id;
		local_optima->customers[local_optima->vehicles[vehicle_two_id]->customers[customer_four_id]->id]->vehicle_route = vehicle_two_id;
        local_optima->vehicles[vehicle_one_id]->customers[customer_one_id]->vehicle_route = vehicle_one_id;
		local_optima->customers[local_optima->vehicles[vehicle_one_id]->customers[customer_one_id]->id]->vehicle_route = vehicle_one_id;
        local_optima->vehicles[vehicle_two_id]->customers[customer_two_id]->vehicle_route = vehicle_two_id;
		local_optima->customers[local_optima->vehicles[vehicle_two_id]->customers[customer_two_id]->id]->vehicle_route = vehicle_two_id;

        return local_optima;
    } else {
        return initial_solution;
    }
}

RelocateNeighborhood::RelocateNeighborhood()
{
}

RelocateNeighborhood::RelocateNeighborhood(bool first_improvement /*= false*/)
        : Neighborhood(first_improvement, "B")
{
}

RelocateNeighborhood::~RelocateNeighborhood()
{
}

std::shared_ptr <Solution> RelocateNeighborhood::search_neighbors(std::shared_ptr <Solution> initial_solution, const Instance& instance)
{
    int n_vehicles = initial_solution->vehicles.size();

    float min_delta = 0.0;
    vector<int> best_move;
    Vehicle vehicle_one;
    Vehicle vehicle_two;
    // customers
    Customer a;
    Customer b;
    Customer c;

    Customer f;
    Customer g;
    // demands
    int route1_new_demand;
    int route2_new_demand;
    // dist_matrix
    // delta
    float delta = 0.0;


    for(int vehicle_one_id = 0; vehicle_one_id < n_vehicles; vehicle_one_id++) {
        vehicle_one = *initial_solution->vehicles[vehicle_one_id];
        for(int customer_one_id = 1; customer_one_id < vehicle_one.customers.size() - 1; customer_one_id++) {
            for(int vehicle_two_id = 0; vehicle_two_id < n_vehicles; vehicle_two_id++) {
                // if(vehicle_two_id == vehicle_one_id){
                vehicle_two = *initial_solution->vehicles[vehicle_two_id];
                for(int customer_new_pos = 1; customer_new_pos < vehicle_two.customers.size() - 1; customer_new_pos++) {
                    /*
                        route1 =  [A B C D]
                        route2 =  [E F G H]
                        Relocate(B, pos(G)) leaving us with:
                        new_route1 = [A C D]
                        new_route2 = [E F B G H]
                    */

                    a = *vehicle_one.customers[customer_one_id - 1];
                    b = *vehicle_one.customers[customer_one_id];
                    c = *vehicle_one.customers[customer_one_id + 1];

                    f = *vehicle_two.customers[customer_new_pos - 1];
                    g = *vehicle_two.customers[customer_new_pos];

                    // Check constraint of load violation
                    route2_new_demand = vehicle_two.load + b.demand;

                    if(((vehicle_one_id == vehicle_two_id) && (c.id == g.id)) ||
                       ((vehicle_one_id == vehicle_two_id) && (b.id == g.id)) || (route2_new_demand > vehicle_two.Q)) {
                        continue;
                    }

                    /*
                        Evaluate the improvement delta of both route1 and route2, using distance matrix
                        Delta = D(a,c) + D(f,b) + D(b,g) - D(a,b) -D(b,c) - D(f,g);
                         *
                    delta = dist[a.id][c.id] + dist[f.id][b.id] + dist[b.id][g.id] \
                            - dist[a.id][b.id] - dist[b.id][c.id] - dist[f.id][g.id];

                     delta = instance.insertion_matrix[f.id][g.id][b.id] + instance.remotion_matrix[a.id][c.id][b.id];
                    */

                    delta = instance.dist_matrix[a.id][c.id] + instance.dist_matrix[f.id][b.id] + instance.dist_matrix[b.id][g.id] \
                            - instance.dist_matrix[a.id][b.id] - instance.dist_matrix[b.id][c.id] - instance.dist_matrix[f.id][g.id];

                    // Store the current best possible move
                    if(delta < min_delta) {
                        min_delta = delta;
                        best_move.clear();
                        best_move = { vehicle_one_id, vehicle_two_id, customer_one_id, customer_new_pos };
                        if(this->first_improvement) {
                            initial_solution->cost = initial_solution->cost + min_delta;
                            return this->apply_best_move(initial_solution, best_move);
                        }
                    }
                }
                //}
            }
        }
    }
    if(min_delta < 0) {
        /*
        cout << "min_delta " << min_delta << " C1 " << best_move[2] << " C2 " << best_move[3] << endl;
        float veh_1_cost = initial_solution->vehicles[best_move[0]]->compute_cost(instance);
        float veh_2_cost = initial_solution->vehicles[best_move[1]]->compute_cost(instance);
        cout << "cost_before " << veh_1_cost << " " << veh_2_cost << endl;
*/
        initial_solution = apply_best_move(initial_solution, best_move);
        initial_solution->cost = initial_solution->cost + min_delta;
        initial_solution->vehicles[best_move[0]]->compute_cost(instance);
        initial_solution->vehicles[best_move[1]]->compute_cost(instance);
        /*
                veh_1_cost = initial_solution->vehicles[best_move[0]]->compute_cost(instance);
                veh_2_cost = initial_solution->vehicles[best_move[1]]->compute_cost(instance);
                cout << "cost_after " << veh_1_cost << " " << veh_2_cost << endl;
                 * */
    }
    // return this->apply_best_move(initial_solution, best_move);
    return initial_solution;
}

/**
 * Apply the relocate (insertion) to the best possible move in the neighborhood
 * @param initial_solution Current solution of the VRP
 * @param best_move {vehicle_one_id, vehicle_two_id, customer_one_id, customer_new_pos}
 *
 * @return local_optima: Neew VRP solution corresponding to the local optima of this neighborhood
 *          if there's no local optima, returns initial_solution instead
 */
std::shared_ptr <Solution> RelocateNeighborhood::apply_best_move(std::shared_ptr <Solution> initial_solution, vector<int>& best_move)
{
    if(best_move.size() != 0) {
        std::shared_ptr <Solution> local_optima = initial_solution;
        int vehicle_one_id = best_move[0];
        int vehicle_two_id = best_move[1];
        int customer_one_id = best_move[2];
        int customer_new_pos = best_move[3];

        shared_ptr<Customer> customer_one = local_optima->vehicles[vehicle_one_id]->customers[customer_one_id];

        if(vehicle_one_id == vehicle_two_id) {
            if(customer_new_pos > customer_one_id) {
                customer_new_pos--;
            }
        }

        local_optima->vehicles[vehicle_one_id]->remove_customer(customer_one_id);
        local_optima->vehicles[vehicle_two_id]->insert_customer(customer_one, customer_new_pos);

        local_optima->vehicles[vehicle_two_id]->customers[customer_new_pos]->vehicle_route = vehicle_two_id;
		local_optima->customers[local_optima->vehicles[vehicle_two_id]->customers[customer_new_pos]->id]->vehicle_route = vehicle_two_id;
		
        return local_optima;
    } else {
        return initial_solution;
    }
}

TwoOptIntraNeighborhood::TwoOptIntraNeighborhood()
{
}

TwoOptIntraNeighborhood::TwoOptIntraNeighborhood(bool first_improvement /*= false*/)
        : Neighborhood(first_improvement, "C")
{
}

TwoOptIntraNeighborhood::~TwoOptIntraNeighborhood()
{
}

std::shared_ptr <Solution> TwoOptIntraNeighborhood::search_neighbors(std::shared_ptr <Solution> solution, const Instance& instance)
{

    float min_delta = 0;
    vector<int> best_move;
    // customers
    Customer a;
    Customer b;
    Customer c;
    Customer d;
    // delta
    float delta;
    int vehicle_idx = 0;
    for(const auto& vehicle : solution->vehicles) {
        for(int i = 1; i < vehicle->customers.size() - 1; i++) {
            for(int j = i + 1; j < vehicle->customers.size() - 1; j++) {
                a = *vehicle->customers[i];
                b = *vehicle->customers[i + 1];
                c = *vehicle->customers[j];
                d = *vehicle->customers[j + 1];

                // No need to check for Load constraint violation, because it is intra-route

                delta = instance.dist_matrix[a.id][c.id] + instance.dist_matrix[b.id][d.id] -
                        instance.dist_matrix[a.id][b.id] - instance.dist_matrix[c.id][d.id];

                // Store the current best possible move

                if(delta < min_delta) {
                    min_delta = delta;
                    best_move.clear();
                    best_move = { vehicle_idx, i, j };
                    if(this->first_improvement)
                        return this->apply_best_move(solution, best_move);
                }
            }
        }
        vehicle_idx++;
    }
    if(min_delta < 0) {
        /*
        cout << "min_delta " << min_delta << " C1 " << best_move[2] << " C2 " << best_move[3] << endl;
        float veh_1_cost = solution->vehicles[best_move[0]]->compute_cost(instance);
        float veh_2_cost = solution->vehicles[best_move[1]]->compute_cost(instance);
        cout << "cost_before " << veh_1_cost << " " << veh_2_cost << endl;
        */
        solution = apply_best_move(solution, best_move);
        solution->cost = solution->cost + min_delta;
        /*
        veh_1_cost = solution->vehicles[best_move[0]]->compute_cost(instance);
        veh_2_cost = solution->vehicles[best_move[1]]->compute_cost(instance);
        cout << "cost_after " << veh_1_cost << " " << veh_2_cost << endl;
         * */
    }
    return solution;

    // return this->apply_best_move(solution, best_move);
}

/**
 * Apply the 2-opt move to the best possible move in the neighborhood
 * @param initial_solution Current solution of the VRP
 * @param best_move {vehicle_idx, i, j}
 *
 * @return local_optima: Neew VRP solution corresponding to the local optima of this neighborhood
 *          if there's no local optima, returns initial_solution instead
 */
std::shared_ptr <Solution> TwoOptIntraNeighborhood::apply_best_move(std::shared_ptr <Solution> solution, vector<int>& best_move)
{
    if(best_move.size() != 0) {
        std::shared_ptr <Solution> local_optima = solution;
        int vehicle_idx = best_move[0];
        int i = best_move[1];
        int j = best_move[2];
        Vehicle& vehicle = *local_optima->vehicles[vehicle_idx];
        vector<shared_ptr<Customer>> aux_vector;

        for(int h = 0; h < i + 1; h++)
            aux_vector.push_back(vehicle.customers[h]);
        for(int h = j; h > i; h--)
            aux_vector.push_back(vehicle.customers[h]);
        for(int h = j + 1; h < vehicle.customers.size(); h++)
            aux_vector.push_back(vehicle.customers[h]);

        vehicle.customers = aux_vector;

        return local_optima;
    } else {
        return solution;
    }
}

TwoOptInterNeighborhood::TwoOptInterNeighborhood()
{
}

TwoOptInterNeighborhood::TwoOptInterNeighborhood(bool first_improvement /*= false*/)
        : Neighborhood(first_improvement, "D")
{
}

TwoOptInterNeighborhood::~TwoOptInterNeighborhood()
{
}

std::shared_ptr <Solution> TwoOptInterNeighborhood::search_neighbors(std::shared_ptr <Solution> aux_solution, const Instance& instance)
{
    int x_nearest = aux_solution->customers.size() / 1;
    aux_solution->compute_cummulative_demands();
	
    float min_delta = 0;
    vector<int> best_move;
    // customers
    Customer a;
    Customer b;
    Customer c;
    Customer d;
    // demands
    int route1_new_demand;
    int route2_new_demand;
    // delta
    float delta;
    float delta_2 = 0;
    int vehicle1_idx = 0;
    int vehicle2_idx = 0;
    int demand_1 = 0;
    int demand_2 = 0;
    int demand_3 = 0;
    int demand_4 = 0;
    int c_id = 0;
    int customer2_id = -1;
    bool isFeasible = false;
	
    for(auto& vehicle1 : aux_solution->vehicles) {
		//error esta aqui, al cambiar de vehiculo, del 27 al 28 que no deberoa de existir
		
        for(int i = 1; i < vehicle1->customers.size() - 1; i++) {
			
            a = *vehicle1->customers[i];
			
            b = *vehicle1->customers[i + 1];
			
            demand_1 = vehicle1->get_demand_between(0, i);
			
            demand_2 = vehicle1->get_demand_between(i + 1, vehicle1->customers.size() - 1);
			
            for(int j = 0; j < x_nearest; j++) {
				
				
                c_id = instance.sorted_dist_matrix[a.id][j];
                
                vehicle2_idx = aux_solution->customers[c_id]->vehicle_route;
				
                if((vehicle1_idx == vehicle2_idx)||(vehicle2_idx == aux_solution->vehicles.size())) {
                    // vehicle2_idx++;
                    continue;
                }
				
                auto& vehicle2 = aux_solution->vehicles[vehicle2_idx];
				
                customer2_id = -1;
                for(int k = 1; k < vehicle2->customers.size() - 1; k++) {
                    // //cout << "k:  " << k <<  endl;
                    if(vehicle2->customers[k]->id == c_id) {
                        c = *vehicle2->customers[k];
                        d = *vehicle2->customers[k + 1];
                        customer2_id = k;
                    
                        break;
                    }
                }
                if(customer2_id == -1) {
                    continue;
                }
                isFeasible = false;
                // Parte nueva 2-opt
				//cout << "before demand 3" << endl;
                demand_3 = vehicle2->get_demand_between(0, customer2_id);
                demand_4 = vehicle2->get_demand_between(customer2_id + 1, vehicle2->customers.size() - 1);

                route1_new_demand = demand_1 + demand_3;
                route2_new_demand = demand_2 + demand_4;
                delta_2 = 1;

                if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {
                    delta_2 = instance.dist_matrix[a.id][c.id] + instance.dist_matrix[b.id][d.id] -
                              instance.dist_matrix[a.id][b.id] - instance.dist_matrix[c.id][d.id];
                    isFeasible = true;
                }

                 //cout << "before demand 1" << endl;
                route1_new_demand = demand_1 + demand_4;
                route2_new_demand = demand_2 + demand_3;

                delta = 1;
				
                if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {
                    delta = instance.dist_matrix[a.id][d.id] + instance.dist_matrix[b.id][c.id] -
                            instance.dist_matrix[a.id][b.id] - instance.dist_matrix[c.id][d.id];
                    isFeasible = true;
                }
				
                if(isFeasible == false) {
                    continue;
                }

                // Store the current best possible move
				
                if(delta_2 < delta) {
                    if(delta_2 < min_delta) {
                        min_delta = delta_2;
                        best_move.clear();
                        best_move = { vehicle1_idx, vehicle2_idx, i, customer2_id, 1 };
                        if(this->first_improvement)
                            return this->apply_best_move(aux_solution, best_move);
                    }
                } else {
                    if(delta < min_delta) {
                        min_delta = delta;
                        best_move.clear();
                        best_move = { vehicle1_idx, vehicle2_idx, i, customer2_id, 0 };
                        if(this->first_improvement)
                            return this->apply_best_move(aux_solution, best_move);
                    }
                }
            }
        }
		
        vehicle1_idx++;
    }
    if(min_delta < 0) {
        /*
        cout << "min_delta " << min_delta << " C1 " << best_move[2] << " C2 " << best_move[3] << endl;
        float veh_1_cost = aux_solution->vehicles[best_move[0]]->compute_cost(instance);
        float veh_2_cost = aux_solution->vehicles[best_move[1]]->compute_cost(instance);
        cout << "cost_before " << veh_1_cost << " " << veh_2_cost << endl;
        */
        /*
        for(int i = 0; i < aux_solution->vehicles[best_move[0]]->customers.size(); i++) {
            cout << aux_solution->vehicles[best_move[0]]->customers[i]->id << " ";
        }
        cout << " " << endl;
        for(int i = 0; i < aux_solution->vehicles[best_move[1]]->customers.size(); i++) {
            cout << aux_solution->vehicles[best_move[1]]->customers[i]->id << " ";
        }
        cout << " " << endl;
        */
        aux_solution = apply_best_move(aux_solution, best_move);
        aux_solution->cost = aux_solution->cost + min_delta;
        aux_solution->vehicles[best_move[0]]->compute_cost(instance);
        aux_solution->vehicles[best_move[1]]->compute_cost(instance);
        /*
        for(int i = 0; i < aux_solution->vehicles[best_move[0]]->customers.size(); i++) {
            cout << aux_solution->vehicles[best_move[0]]->customers[i]->id << " ";
        }
        cout << " " << endl;
        for(int i = 0; i < aux_solution->vehicles[best_move[1]]->customers.size(); i++) {
            cout << aux_solution->vehicles[best_move[1]]->customers[i]->id << " ";
        }
        cout << " " << endl << " ------------ " << endl;
         */
        /*
        veh_1_cost = aux_solution->vehicles[best_move[0]]->compute_cost(instance);
        veh_2_cost = aux_solution->vehicles[best_move[1]]->compute_cost(instance);
        cout << "cost_after " << veh_1_cost << " " << veh_2_cost << endl;
         * */
    }
    return aux_solution;
    // return this->apply_best_move(aux_solution, best_move);
}

std::shared_ptr <Solution> TwoOptInterNeighborhood::apply_best_move(std::shared_ptr <Solution> solution, vector<int>& best_move)
{
    if(best_move.size() != 0) {
        std::shared_ptr <Solution> local_optima = solution;
        int vehicle1_idx = best_move[0];
        int vehicle2_idx = best_move[1];
        int i = best_move[2];
        int j = best_move[3];
        int invert = best_move[4];
        if(invert == 0) {
            Vehicle& vehicle1 = *local_optima->vehicles[vehicle1_idx];
            Vehicle& vehicle2 = *local_optima->vehicles[vehicle2_idx];
            vector<shared_ptr<Customer>> aux_vector;
            vector<shared_ptr<Customer>> temp_vehicle1_customers(vehicle1.customers.size());

            for(int cust = 0; cust < vehicle1.customers.size(); cust++) {
                temp_vehicle1_customers[cust] = vehicle1.customers[cust];
            }

            for(int h = 0; h < i + 1; h++) {
                aux_vector.push_back(vehicle1.customers[h]);
            }
            for(int h = j + 1; h < vehicle2.customers.size(); h++) {
                aux_vector.push_back(vehicle2.customers[h]);
            }

            vehicle1.customers = vector<shared_ptr<Customer>>(aux_vector.size());
            for(int cust = 0; cust < aux_vector.size(); cust++) {
                vehicle1.customers[cust] = aux_vector[cust];
            }

            vehicle1.update_load();
            aux_vector.clear();

            for(int h = temp_vehicle1_customers.size() - 1; h > i; h--) {
                aux_vector.push_back(temp_vehicle1_customers[h]);
            }
            for(int h = j; h > -1; h--) {
                aux_vector.push_back(vehicle2.customers[h]);
            }

            vehicle2.customers = vector<shared_ptr<Customer>>(aux_vector.size());
            for(int cust = 0; cust < aux_vector.size(); cust++) {
                vehicle2.customers[cust] = aux_vector[cust];
            }
            vehicle2.update_load();

            for(int k = 1; k < vehicle1.customers.size() - 1; k++) {
                vehicle1.customers[k]->vehicle_route = vehicle1_idx;
				local_optima->customers[vehicle1.customers[k]->id]->vehicle_route = vehicle1_idx;
            }

            for(int k = 1; k < vehicle2.customers.size() - 1; k++) {
                vehicle2.customers[k]->vehicle_route = vehicle2_idx;
				local_optima->customers[vehicle2.customers[k]->id]->vehicle_route = vehicle2_idx;
            }
        } else {
            Vehicle& vehicle1 = *local_optima->vehicles[vehicle1_idx];
            Vehicle& vehicle2 = *local_optima->vehicles[vehicle2_idx];
            vector<shared_ptr<Customer>> aux_vector;
            vector<shared_ptr<Customer>> temp_vehicle1_customers(vehicle1.customers.size());

            for(int cust = 0; cust < vehicle1.customers.size(); cust++) {
                temp_vehicle1_customers[cust] = vehicle1.customers[cust];
            }

            for(int h = 0; h < i + 1; h++) {
                aux_vector.push_back(vehicle1.customers[h]);
            }
            for(int h = j; h > -1; h--) {
                aux_vector.push_back(vehicle2.customers[h]);
            }

            vehicle1.customers = vector<shared_ptr<Customer>>(aux_vector.size());
            for(int cust = 0; cust < aux_vector.size(); cust++) {
                vehicle1.customers[cust] = aux_vector[cust];
            }

            vehicle1.update_load();
            aux_vector.clear();

            for(int h = temp_vehicle1_customers.size() - 1; h > i; h--) {
                aux_vector.push_back(temp_vehicle1_customers[h]);
            }
            for(int h = j + 1; h < vehicle2.customers.size(); h++) {
                aux_vector.push_back(vehicle2.customers[h]);
            }

            vehicle2.customers = vector<shared_ptr<Customer>>(aux_vector.size());
            for(int cust = 0; cust < aux_vector.size(); cust++) {
                vehicle2.customers[cust] = aux_vector[cust];
            }
            vehicle2.update_load();

            for(int k = 1; k < vehicle1.customers.size() - 1; k++) {
                vehicle1.customers[k]->vehicle_route = vehicle1_idx;
				local_optima->customers[vehicle1.customers[k]->id]->vehicle_route = vehicle1_idx;
            }

            for(int k = 1; k < vehicle2.customers.size() - 1; k++) {
                vehicle2.customers[k]->vehicle_route = vehicle2_idx;
				local_optima->customers[vehicle2.customers[k]->id]->vehicle_route = vehicle2_idx;
            }
        }
        return local_optima;
    } else {
        return solution;
    }
}

ThreeOptInterNeighborhood::ThreeOptInterNeighborhood()
{
}

ThreeOptInterNeighborhood::ThreeOptInterNeighborhood(bool first_improvement /*= false*/)
        : Neighborhood(first_improvement, "E")
{
}

ThreeOptInterNeighborhood::~ThreeOptInterNeighborhood()
{
}

std::shared_ptr <Solution> ThreeOptInterNeighborhood::search_neighbors(std::shared_ptr <Solution> aux_solution, const Instance& instance)
{

    int x_nearest = aux_solution->customers.size() / 3;
    aux_solution->compute_cummulative_demands();

    float min_delta = 0;
    vector<int> best_move;
    // customers
    Customer a;
    Customer b;
    Customer c;
    Customer d;
    Customer e;
    Customer f;
    Customer g;
    Customer h;
    Customer r;
    Customer s;
    // demands
    int route1_new_demand;
    int route2_new_demand;
    int route3_new_demand;
    // delta
    float delta_1 = 0;
    float delta_2 = 0;
    float delta_3 = 0;
    float delta_4 = 0;
    float delta_5 = 0;
    float delta_6 = 0;
    float delta_7 = 0;
    float delta_8 = 0;
    float delta_9 = 0;
    float delta_10 = 0;
    float delta_11 = 0;
    float delta_12 = 0;
    float delta_13 = 0;
    float delta_14 = 0;
    float delta_15 = 0;
    float delta_16 = 0;
    float delta_17 = 0;
    float delta_18 = 0;
    float delta_21 = 0;
    float delta_22 = 0;
    float delta_23 = 0;
    float delta_24 = 0;
    float delta_25 = 0;
    float delta_26 = 0;
    string min_delta_name;
    int vehicle1_idx = 0;
    int vehicle2_idx = 0;
    int vehicle3_idx = 0;
    int demand_1 = 0;
    int demand_2 = 0;
    int demand_3 = 0;
    int demand_4 = 0;
    int demand_5 = 0;
    int demand_6 = 0;
    int demand_11 = 0;
    int demand_12 = 0;
    int demand_13 = 0;
    int demand_14 = 0;
    int demand_15 = 0;
    int demand_16 = 0;
    int c_id = 0;
    int e_id = 0;
    int customer2_id = -1;
    int customer3_id = -1;
    bool isFeasible = false;
    int doSwapOrRelocateOrInvert = 0; // 0 if none, 1 if swap, 2 if relocate

    for(auto& vehicle1 : aux_solution->vehicles) {
        for(int i = 1; i < vehicle1->customers.size() - 1; i++) {
            a = *vehicle1->customers[i];
            b = *vehicle1->customers[i + 1];
            demand_1 = vehicle1->get_demand_between(0, i);
            demand_2 = vehicle1->get_demand_between(i + 1, vehicle1->customers.size() - 1);
            // cout <<"Determina demandas 1 y 2" << endl;
            for(int j = 0; j < x_nearest; j++) {

                c_id = instance.sorted_dist_matrix[a.id][j];
                // cout << c_id << endl;
                vehicle2_idx = aux_solution->customers[c_id]->vehicle_route;

                if((vehicle1_idx == vehicle2_idx)||(vehicle2_idx == aux_solution->vehicles.size())) {
                    // vehicle2_idx++;
                    continue;
                }
                auto& vehicle2 = aux_solution->vehicles[vehicle2_idx];

                customer2_id = -1;
                for(int k = 1; k < vehicle2->customers.size() - 1; k++) {
                    if(vehicle2->customers[k]->id == c_id) {
                        c = *vehicle2->customers[k];
                        d = *vehicle2->customers[k + 1];
                        customer2_id = k;
                        break;
                    }
                }
                if(customer2_id == -1) {
                    continue;
                }

                demand_3 = vehicle2->get_demand_between(0, customer2_id);
                demand_4 = vehicle2->get_demand_between(customer2_id + 1, vehicle2->customers.size() - 1);

                // Aquí empiezan las combinaciones de 2-opt
                /*
                Route 1 = [* * A B * *]  demand_1 = [0 - A]  demand_2 = [B - 0]
                Route 2 = [* * C D * *]  demand_3 = [0 - C]  demand_4 = [D - 0]
                demand_1 = [0 - A]
                demand_2 = [B - 0]
                demand_3 = [0 - C]
                demand_4 = [D - 0]
                  */

                // Option 9
                // 1-2, 2-1, 1-2
                // 0-A-C-0, 0-D-B-0, 0-A-C-0
                route1_new_demand = demand_1 + demand_3;
                route2_new_demand = demand_4 + demand_2;
                delta_9 = 1;

                if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {
                    delta_9 = instance.dist_matrix[a.id][c.id] + instance.dist_matrix[d.id][b.id] -
                              instance.dist_matrix[a.id][b.id] - instance.dist_matrix[c.id][d.id];
                    isFeasible = true;
                    if(delta_9 < min_delta) {
                        min_delta = delta_9;
                        min_delta_name = "delta_9";
                        best_move.clear();
                        best_move = {
                                vehicle1_idx,
                                vehicle2_idx,
                                vehicle1_idx,
                                i,
                                0,
                                customer2_id,
                                0,
                                customer2_id + 1,
                                1,
                                i + 1,
                                1,
                        };
                        if(this->first_improvement)
                            return this->apply_best_move(aux_solution, best_move);
                    }
                }
                /*
                    Route 1 = [* * A B * *]  demand_1 = [0 - A]  demand_2 = [B - 0]
                    Route 2 = [* * C D * *]  demand_3 = [0 - C]  demand_4 = [D - 0]
                    demand_1 = [0 - A]
                    demand_2 = [B - 0]
                    demand_3 = [0 - C]
                    demand_4 = [D - 0]
                      */

                // Option 10
                // 1-2, 2-1, 1-2
                // 0-A-D-0, 0-C-B-0, 0-A-D-0
                route1_new_demand = demand_1 + demand_4;
                route2_new_demand = demand_3 + demand_2;
                delta_10 = 1;

                if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {
                    delta_10 = instance.dist_matrix[a.id][d.id] + instance.dist_matrix[c.id][b.id] -
                               instance.dist_matrix[a.id][b.id] - instance.dist_matrix[c.id][d.id];
                    isFeasible = true;
                    if(delta_10 < min_delta) {
                        min_delta = delta_10;
                        min_delta_name = "delta_10";
                        best_move.clear();
                        best_move = {
                                vehicle1_idx,
                                vehicle2_idx,
                                vehicle1_idx,
                                i,
                                0,
                                customer2_id + 1,
                                1,
                                customer2_id,
                                0,
                                i + 1,
                                1,
                        };
                        if(this->first_improvement)
                            return this->apply_best_move(aux_solution, best_move);
                    }
                }

                // cout <<"Determina demandas 3 y 4" << endl;

                for(int w = 0; w < x_nearest; w++) {

                    e_id = instance.sorted_dist_matrix[c.id][w];
                    vehicle3_idx = aux_solution->customers[e_id]->vehicle_route;

                    
					if((vehicle2_idx == vehicle3_idx)||(vehicle3_idx == aux_solution->vehicles.size())) {
                        continue;
                    }
                    if(vehicle1_idx == vehicle3_idx) {
						
                        // Aquí empiezan las combinaciones de 3-opt 2 rutas
                        auto& vehicle3 = aux_solution->vehicles[vehicle3_idx];

                        customer3_id = -1;
                        if(vehicle3->customers.size() >= i + 3) {
                            for(int z = i; z < vehicle3->customers.size() - 1; z++) {
                                if(vehicle3->customers[z]->id == e_id) {
                                    r = *vehicle3->customers[i - 1];
                                    h = *vehicle3->customers[z - 1];
                                    e = *vehicle3->customers[z];
                                    f = *vehicle3->customers[z + 1];
                                    if(f.id != 0) {
                                        s = *vehicle3->customers[z + 2];
                                    }
                                    g = *vehicle3->customers[i + 2];

                                    customer3_id = z;
                                    // cout << "llega hasta 4 " << endl;
                                    break;
                                }
                            }
                        }
                        if((customer3_id == -1) || (e.id == a.id) || (e.id == b.id) || (f.id == a.id) ||
                           (f.id == b.id)) {
                            continue;
                        }
                        isFeasible = false;

                        demand_11 = vehicle1->get_demand_between(0, i);
                        demand_12 = vehicle1->get_demand_between(i + 1, customer3_id);
                        demand_13 = vehicle3->get_demand_between(customer3_id + 1, vehicle3->customers.size() - 1);
                        demand_14 = vehicle2->get_demand_between(0, customer2_id);
                        demand_15 = vehicle2->get_demand_between(customer2_id + 1, vehicle2->customers.size() - 1);
                        // Option 21
                        // 1-2, 1-1, 2-1
                        // Swap 0-A-B-G-**-H-E-F-0,    0-D-B-G-**-H-E-C-0
                        //                          2-opt 0-A-F-0,
                        route1_new_demand = demand_15 + demand_12 + demand_14;
                        route2_new_demand = demand_11 + demand_13;
                        delta_21 = 1;

                        if((b.id == 0) || (e.id == 0)) {
                            continue; // No podemos hacer inversion con el Depot
                        }
                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {
                            delta_21 = instance.dist_matrix[d.id][b.id] + instance.dist_matrix[e.id][c.id] +
                                       instance.dist_matrix[a.id][f.id] - instance.dist_matrix[a.id][b.id] -
                                       instance.dist_matrix[e.id][f.id] - instance.dist_matrix[c.id][d.id];

                            doSwapOrRelocateOrInvert = 4; // Do different movement

                            isFeasible = true;
                            if(delta_21 < min_delta) {
                                min_delta = delta_21;
                                min_delta_name = "delta_21";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, customer2_id + 1, 1, i + 1, 1,
                                              customer3_id, 0, customer2_id, 0, i, 0, customer3_id + 1, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        // Option 22
                        // 1-2, 1-1, 2-1
                        // Swap 0-A-B-G-**-H-E-F-0,    0-C-B-G-**-H-E-D-0
                        //                          2-opt 0-A-F-0,
                        route1_new_demand = demand_15 + demand_12 + demand_14;
                        route2_new_demand = demand_11 + demand_13;
                        delta_22 = 1;

                        if((b.id == 0) || (e.id == 0)) {
                            continue; // No podemos hacer inversion con el Depot
                        }
                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {
                            delta_22 = instance.dist_matrix[c.id][b.id] + instance.dist_matrix[e.id][d.id] +
                                       instance.dist_matrix[a.id][f.id] - instance.dist_matrix[a.id][b.id] -
                                       instance.dist_matrix[e.id][f.id] - instance.dist_matrix[c.id][d.id];

                            doSwapOrRelocateOrInvert = 4; // Do different movement

                            isFeasible = true;
                            if(delta_22 < min_delta) {
                                min_delta = delta_22;
                                min_delta_name = "delta_22";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, customer2_id, 0, i + 1, 1,
                                              customer3_id, 0, customer2_id + 1, 1, i, 0, customer3_id + 1, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        // Option 23
                        // 1-2, 1-1, 2-1
                        // Swap 0-A-B-G-**-H-E-F-0,    0-C-E-H**-G-B-F-0
                        //                          2-opt 0-D-A-0,
                        route1_new_demand = demand_13 + demand_12 + demand_14;
                        route2_new_demand = demand_11 + demand_15;
                        delta_23 = 1;

                        if((b.id == 0) || (e.id == 0)) {
                            continue; // No podemos hacer inversion con el Depot
                        }
                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {
                            delta_23 = instance.dist_matrix[c.id][e.id] + instance.dist_matrix[b.id][f.id] +
                                       instance.dist_matrix[d.id][a.id] - instance.dist_matrix[a.id][b.id] -
                                       instance.dist_matrix[e.id][f.id] - instance.dist_matrix[c.id][d.id];
                            if(e.id == g.id) {
                                // Significa que e es adyacente a b, por lo tanto se realiza un relocate
                                doSwapOrRelocateOrInvert = 2; // Do relocate
                            } else {
                                doSwapOrRelocateOrInvert = 3; // Do inversion
                            }

                            isFeasible = true;
                            if(delta_23 < min_delta) {
                                min_delta = delta_23;
                                min_delta_name = "delta_23";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, i, 0, customer2_id + 1, 1,
                                              customer2_id, 0, i + 1, 1, i + 1, 1, customer3_id, 1, doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        // Option 24
                        // 1-2, 1-1, 2-1
                        // Swap 0-A-B-G-**-H-E-F-0,    0-D-E-H**-G-B-F-0
                        //                          2-opt 0-C-A-0,
                        route1_new_demand = demand_13 + demand_12 + demand_15;
                        route2_new_demand = demand_11 + demand_14;
                        delta_24 = 1;

                        if((b.id == 0) || (e.id == 0)) {
                            continue; // No podemos hacer inversion con el Depot
                        }
                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {
                            delta_24 = instance.dist_matrix[d.id][e.id] + instance.dist_matrix[b.id][f.id] +
                                       instance.dist_matrix[c.id][a.id] - instance.dist_matrix[a.id][b.id] -
                                       instance.dist_matrix[e.id][f.id] - instance.dist_matrix[c.id][d.id];
                            if(e.id == g.id) {
                                // Significa que e es adyacente a b, por lo tanto se realiza un relocate
                                doSwapOrRelocateOrInvert = 2; // Do relocate
                            } else {
                                doSwapOrRelocateOrInvert = 3; // Do inversion
                            }

                            isFeasible = true;
                            if(delta_24 < min_delta) {
                                min_delta = delta_24;
                                min_delta_name = "delta_24";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, i, 0, customer2_id, 0,
                                              customer2_id + 1, 1, i + 1, 1, i + 1, 1, customer3_id, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        // Option 25
                        // 1-2, 1-1, 2-1
                        // Swap 0-A-B-G-**-H-E-F-0,    0-A-E-H**-G-B-D-0
                        //                          2-opt 0-C-F-0,
                        route1_new_demand = demand_11 + demand_12 + demand_15;
                        route2_new_demand = demand_13 + demand_14;
                        delta_25 = 1;

                        if((b.id == 0) || (e.id == 0)) {
                            continue; // No podemos hacer inversion con el Depot
                        }

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {
                            if(e.id == g.id) {
                                // Significa que e es adyacente a b, por lo tanto se realiza un relocate
                                // 0-A-B-E-F-0    -> 0-A-E-B-D-0
                                delta_25 = instance.dist_matrix[a.id][e.id] + instance.dist_matrix[b.id][d.id] +
                                           instance.dist_matrix[c.id][f.id] - instance.dist_matrix[a.id][b.id] -
                                           instance.dist_matrix[e.id][f.id] - instance.dist_matrix[c.id][d.id];
                                doSwapOrRelocateOrInvert = 2; // Do relocate
                            } else {
                                delta_25 = instance.dist_matrix[a.id][e.id] + instance.dist_matrix[b.id][d.id] +
                                           instance.dist_matrix[c.id][f.id] - instance.dist_matrix[a.id][b.id] -
                                           instance.dist_matrix[e.id][f.id] - instance.dist_matrix[c.id][d.id];
                                doSwapOrRelocateOrInvert = 3; // Do inversion
                            }

                            isFeasible = true;
                            if(delta_25 < min_delta) {
                                min_delta = delta_25;
                                min_delta_name = "delta_25";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, customer3_id, 0,
                                              customer2_id + 1, 1, customer2_id, 0, customer3_id + 1, 1, i + 1, 1, customer3_id,
                                              1, doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        // Option 26
                        // 1-2, 1-1, 2-1
                        // Swap 0-A-B-G-**-H-E-F-0,    0-A-E-H**-G-B-C-0
                        //                          2-opt 0-D-F-0,
                        route1_new_demand = demand_11 + demand_12 + demand_14;
                        route2_new_demand = demand_13 + demand_15;
                        delta_26 = 1;

                        if((b.id == 0) || (e.id == 0)) {
                            continue; // No podemos hacer inversion con el Depot
                        }

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {
                            if(e.id == g.id) {
                                // Significa que e es adyacente a b, por lo tanto se realiza un relocate
                                // 0-A-B-E-F-0    -> 0-A-E-B-C-0
                                delta_26 = instance.dist_matrix[a.id][e.id] + instance.dist_matrix[b.id][c.id] +
                                           instance.dist_matrix[d.id][f.id] - instance.dist_matrix[a.id][b.id] -
                                           instance.dist_matrix[e.id][f.id] - instance.dist_matrix[c.id][d.id];
                                doSwapOrRelocateOrInvert = 2; // Do relocate
                            } else {
                                delta_26 = instance.dist_matrix[a.id][e.id] + instance.dist_matrix[b.id][c.id] +
                                           instance.dist_matrix[d.id][f.id] - instance.dist_matrix[a.id][b.id] -
                                           instance.dist_matrix[e.id][f.id] - instance.dist_matrix[c.id][d.id];
                                doSwapOrRelocateOrInvert = 3; // Do inversion
                            }

                            isFeasible = true;
                            if(delta_26 < min_delta) {
                                min_delta = delta_26;
                                min_delta_name = "delta_26";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, customer3_id, 0, customer2_id,
                                              0, customer2_id + 1, 1, customer3_id + 1, 1, i + 1, 1, customer3_id, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        // Option 11
                        // 1-2, 1-1, 2-1
                        // Swap 0-A-B-G-**-H-E-F-0,    0-A-E-G**-H-B-C-0
                        //                          2-opt 0-D-F-0,
                        route1_new_demand = demand_11 + demand_12 + demand_14;
                        route2_new_demand = demand_13 + demand_15;
                        delta_11 = 1;

                        if((b.id == 0) || (e.id == 0)) {
                            continue; // No podemos hacer swap con el Depot
                        }

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {
                            if(e.id == g.id) {
                                // Significa que e es adyacente a b, por lo tanto se realiza un relocate
                                // 0-A-B-E-F-0    -> 0-A-E-B-C-0
                                delta_11 = instance.dist_matrix[a.id][e.id] + instance.dist_matrix[b.id][c.id] +
                                           instance.dist_matrix[d.id][f.id] - instance.dist_matrix[a.id][b.id] -
                                           -instance.dist_matrix[e.id][f.id] - instance.dist_matrix[c.id][d.id];
                                doSwapOrRelocateOrInvert = 2; // Do relocate
                            } else {
                                delta_11 = instance.dist_matrix[a.id][e.id] + instance.dist_matrix[e.id][g.id] +
                                           instance.dist_matrix[h.id][b.id] + instance.dist_matrix[d.id][f.id] +
                                           instance.dist_matrix[b.id][c.id] - instance.dist_matrix[a.id][b.id] -
                                           instance.dist_matrix[b.id][g.id] - instance.dist_matrix[e.id][f.id] -
                                           instance.dist_matrix[c.id][d.id] - instance.dist_matrix[h.id][e.id];
                                doSwapOrRelocateOrInvert = 1; // Do Swap
                            }

                            isFeasible = true;
                            if(delta_11 < min_delta) {
                                min_delta = delta_11;
                                min_delta_name = "delta_11";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, customer3_id, 0, customer2_id,
                                              0, customer2_id + 1, 1, customer3_id + 1, 1, i + 1, 1, customer3_id, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        // Option 12
                        // 1-2, 1-1, 2-1
                        // Swap 0-R-A-B-G-**-H-E-F-0,    0-R-E-B-G**-H-A-C-0
                        //                          2-opt 0-D-F-0,
                        route1_new_demand = demand_11 + demand_12 + demand_14;
                        route2_new_demand = demand_13 + demand_15;
                        delta_12 = 1;

                        if((a.id == 0) || (e.id == 0)) {
                            continue; // No podemos hacer swap con el Depot
                        }

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {

                            delta_12 = instance.dist_matrix[e.id][b.id] + instance.dist_matrix[r.id][e.id] +
                                       instance.dist_matrix[h.id][a.id] + instance.dist_matrix[a.id][c.id] +
                                       instance.dist_matrix[d.id][f.id] - instance.dist_matrix[a.id][b.id] -
                                       instance.dist_matrix[r.id][a.id] - instance.dist_matrix[e.id][f.id] -
                                       instance.dist_matrix[c.id][d.id] - instance.dist_matrix[h.id][e.id];
                            doSwapOrRelocateOrInvert = 1; // Do Swap

                            isFeasible = true;
                            if(delta_12 < min_delta) {
                                min_delta = delta_12;
                                min_delta_name = "delta_12";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, customer3_id, 0, customer2_id,
                                              0, customer2_id + 1, 1, customer3_id + 1, 1, i, 1, customer3_id, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        // Option 13
                        // 1-2, 1-1, 2-1
                        // Swap 0-A-B-G-**-H-E-F-0,    0-C-E-G**-H-B-F-0
                        //                          2-opt 0-D-A-0,
                        route1_new_demand = demand_13 + demand_12 + demand_14;
                        route2_new_demand = demand_11 + demand_15;
                        delta_13 = 1;

                        if((b.id == 0) || (e.id == 0)) {
                            continue; // No podemos hacer swap con el Depot
                        }

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {
                            if(e.id == g.id) {
                                // Significa que e es adyacente a b, por lo tanto se realiza un relocate
                                // 0-A-B-E-F-0    -> 0-C-E-B-F-0, 0-A-D-0
                                delta_13 = instance.dist_matrix[c.id][e.id] + instance.dist_matrix[b.id][f.id] +
                                           instance.dist_matrix[d.id][a.id] - instance.dist_matrix[a.id][b.id] -
                                           instance.dist_matrix[e.id][f.id] - instance.dist_matrix[c.id][d.id];
                                doSwapOrRelocateOrInvert = 2; // Do relocate
                            } else {
                                delta_13 = instance.dist_matrix[c.id][e.id] + instance.dist_matrix[e.id][g.id] +
                                           instance.dist_matrix[h.id][b.id] + instance.dist_matrix[b.id][f.id] +
                                           instance.dist_matrix[d.id][a.id] - instance.dist_matrix[a.id][b.id] -
                                           instance.dist_matrix[b.id][g.id] - instance.dist_matrix[e.id][f.id] -
                                           instance.dist_matrix[c.id][d.id] - instance.dist_matrix[h.id][e.id];
                                doSwapOrRelocateOrInvert = 1; // Do Swap
                            }

                            isFeasible = true;
                            if(delta_13 < min_delta) {
                                min_delta = delta_13;
                                min_delta_name = "delta_13";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, i + 1, 1, customer2_id, 0,
                                              customer2_id + 1, 1, i, 0, i + 1, 1, customer3_id, 1, doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        // Option 14
                        // 1-2, 1-1, 2-1
                        // Swap 0-R-A-B-G-**-H-E-F-0,   0-R-E-B-G-**-H-A-F-0,   0-C-B-G**-H-A-F-0
                        //                          2-opt 0-D-E-R-0,

                        route1_new_demand = demand_13 + demand_12 + demand_14 + a.demand - e.demand;
                        route2_new_demand = demand_11 - a.demand + e.demand + demand_15;
                        delta_14 = 1;

                        if((a.id == 0) || (e.id == 0)) {
                            continue; // No podemos hacer swap con el Depot
                        }

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {

                            delta_14 = instance.dist_matrix[c.id][b.id] + instance.dist_matrix[h.id][a.id] +
                                       instance.dist_matrix[a.id][f.id] + instance.dist_matrix[d.id][e.id] +
                                       instance.dist_matrix[e.id][r.id] - instance.dist_matrix[a.id][b.id] -
                                       instance.dist_matrix[r.id][a.id] - instance.dist_matrix[e.id][f.id] -
                                       instance.dist_matrix[c.id][d.id] - instance.dist_matrix[h.id][e.id];
                            doSwapOrRelocateOrInvert = 1; // Do Swap

                            isFeasible = true;
                            if(delta_14 < min_delta) {
                                min_delta = delta_14;
                                min_delta_name = "delta_14";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, i + 1, 1, customer2_id, 0,
                                              customer2_id + 1, 1, i, 0, i, 1, customer3_id, 1, doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        // Option 15
                        // 1-2, 1-1, 2-1
                        // Swap 0-A-B-G-**-H-E-F-S-0,    0-A-F-G-**-H-E-B-S-0,     0-A-F-G**-H-E-C-0
                        //                          2-opt 0-D-B-S-0,
                        route1_new_demand = demand_11 + demand_12 + demand_14 - b.demand + f.demand;
                        route2_new_demand = demand_13 + demand_15 + b.demand - f.demand;
                        delta_15 = 1;

                        if((b.id == 0) || (f.id == 0)) {
                            continue; // No podemos hacer swap con el Depot
                        }

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {

                            delta_15 = instance.dist_matrix[a.id][f.id] + instance.dist_matrix[f.id][g.id] +
                                       instance.dist_matrix[e.id][c.id] + instance.dist_matrix[d.id][b.id] +
                                       instance.dist_matrix[b.id][s.id] - instance.dist_matrix[a.id][b.id] -
                                       instance.dist_matrix[b.id][g.id] - instance.dist_matrix[e.id][f.id] -
                                       instance.dist_matrix[c.id][d.id] - instance.dist_matrix[f.id][s.id];
                            doSwapOrRelocateOrInvert = 1; // Do Swap

                            isFeasible = true;
                            if(delta_15 < min_delta) {
                                min_delta = delta_15;
                                min_delta_name = "delta_15";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, customer3_id, 0, customer2_id,
                                              0, customer2_id + 1, 1, customer3_id + 1, 1, i + 1, 1, customer3_id + 1, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        // Option 16
                        // 1-2, 1-1, 2-1
                        // Swap 0-A-B-G-**-H-E-F-S-0,    0-A-F-G-**-H-E-B-S-0,     0-A-F-G**-H-E-D-0
                        //                          2-opt 0-C-B-S-0,
                        route1_new_demand = demand_11 + demand_12 + demand_15 - b.demand + f.demand;
                        route2_new_demand = demand_13 + demand_14 + b.demand - f.demand;
                        delta_16 = 1;

                        if((b.id == 0) || (f.id == 0)) {
                            continue; // No podemos hacer swap con el Depot
                        }

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {

                            delta_16 = instance.dist_matrix[a.id][f.id] + instance.dist_matrix[f.id][g.id] +
                                       instance.dist_matrix[e.id][d.id] + instance.dist_matrix[c.id][b.id] +
                                       instance.dist_matrix[b.id][s.id] - instance.dist_matrix[a.id][b.id] -
                                       instance.dist_matrix[b.id][g.id] - instance.dist_matrix[e.id][f.id] -
                                       instance.dist_matrix[c.id][d.id] - instance.dist_matrix[f.id][s.id];
                            doSwapOrRelocateOrInvert = 1; // Do Swap

                            isFeasible = true;
                            if(delta_16 < min_delta) {
                                min_delta = delta_16;
                                min_delta_name = "delta_16";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, customer3_id, 0,
                                              customer2_id + 1, 1, customer2_id, 0, customer3_id + 1, 1, i + 1, 1,
                                              customer3_id + 1, 1, doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        // Option 17
                        // 1-2, 1-1, 2-1
                        // Swap 0-R-A-B-G-**-H-E-F-S-0,    0-R-F-B-G-**-H-E-A-S-0,     0-R-F-B-G-**-H-E-D-0
                        //                          2-opt 0-C-A-S-0,

                        route1_new_demand = demand_11 + demand_12 + demand_15 - a.demand + f.demand;
                        route2_new_demand = demand_13 + demand_14 + a.demand - f.demand;
                        delta_17 = 1;

                        if((a.id == 0) || (f.id == 0)) {
                            continue; // No podemos hacer swap con el Depot
                        }

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {

                            delta_17 = instance.dist_matrix[r.id][f.id] + instance.dist_matrix[f.id][b.id] +
                                       instance.dist_matrix[e.id][d.id] + instance.dist_matrix[c.id][a.id] +
                                       instance.dist_matrix[a.id][s.id] - instance.dist_matrix[r.id][a.id] -
                                       instance.dist_matrix[a.id][b.id] - instance.dist_matrix[e.id][f.id] -
                                       instance.dist_matrix[c.id][d.id] - instance.dist_matrix[f.id][s.id];
                            doSwapOrRelocateOrInvert = 1; // Do Swap

                            isFeasible = true;
                            if(delta_17 < min_delta) {
                                min_delta = delta_17;
                                min_delta_name = "delta_17";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, customer3_id, 0,
                                              customer2_id + 1, 1, customer2_id, 0, customer3_id + 1, 1, i, 1, customer3_id + 1,
                                              1, doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        // Option 18
                        // 1-2, 1-1, 2-1
                        // Swap 0-R-A-B-G-**-H-E-F-S-0,    0-R-F-B-G-**-H-E-A-S-0,     0-R-F-B-G-**-H-E-C-0
                        //                          2-opt 0-D-A-S-0,
                        route1_new_demand = demand_11 + demand_12 + demand_14 - a.demand + f.demand;
                        route2_new_demand = demand_13 + demand_15 + a.demand - f.demand;
                        delta_18 = 1;

                        if((a.id == 0) || (f.id == 0)) {
                            continue; // No podemos hacer swap con el Depot
                        }

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q)) {

                            delta_18 = instance.dist_matrix[r.id][f.id] + instance.dist_matrix[f.id][b.id] +
                                       instance.dist_matrix[e.id][c.id] + instance.dist_matrix[d.id][a.id] +
                                       instance.dist_matrix[a.id][s.id] - instance.dist_matrix[r.id][a.id] -
                                       instance.dist_matrix[a.id][b.id] - instance.dist_matrix[e.id][f.id] -
                                       instance.dist_matrix[c.id][d.id] - instance.dist_matrix[f.id][s.id];
                            doSwapOrRelocateOrInvert = 1; // Do Swap

                            isFeasible = true;
                            if(delta_18 < min_delta) {
                                min_delta = delta_18;
                                min_delta_name = "delta_18";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, customer3_id, 0, customer2_id,
                                              0, customer2_id + 1, 1, customer3_id + 1, 1, i, 1, customer3_id + 1, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                    } else {

                        // Aquí empiezan las combinaciones de 3-opt
                        auto& vehicle3 = aux_solution->vehicles[vehicle3_idx];

                        customer3_id = -1;
                        for(int z = 1; z < vehicle3->customers.size() - 1; z++) {
                            if(vehicle3->customers[z]->id == e_id) {
                                e = *vehicle3->customers[z];
                                f = *vehicle3->customers[z + 1];
                                customer3_id = z;
                                // cout << "llega hasta 4 " << endl;
                                break;
                            }
                        }
                        if(customer3_id == -1) {
                            continue;
                        }
                        isFeasible = false;

                        demand_5 = vehicle3->get_demand_between(0, customer3_id);
                        demand_6 = vehicle3->get_demand_between(customer3_id + 1, vehicle3->customers.size() - 1);

                        // cout <<"Determina demandas 5 y 6" << endl;

                        // cout << a.id << " " << b.id << " " <<c.id << " " <<d.id << " " <<e.id << " " <<f.id << " "
                        // <<endl;

                        /*
                        Route 1 = [* * A B * *]  demand_1 = [0 - A]  demand_2 = [B - 0]
                        Route 2 = [* * C D * *]  demand_3 = [0 - C]  demand_4 = [D - 0]
                        Route 3 = [* * E F * *]  demand_5 = [0 - E]  demand_6 = [F - 0]
                        demand_1 = [0 - A]
                        demand_2 = [B - 0]
                        demand_3 = [0 - C]
                        demand_4 = [D - 0]
                        demand_5 = [0 - E]
                        demand_6 = [F - 0]
                          */
                        // Option 1
                        // 1-2, 2-3, 3-1
                        // 0-A-C-0, 0-D-E-0, 0-F-B-0
                        route1_new_demand = demand_1 + demand_3;
                        route2_new_demand = demand_4 + demand_5;
                        route3_new_demand = demand_6 + demand_2;
                        delta_1 = 1;

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q) &&
                           (route3_new_demand <= vehicle3->Q)) {
                            delta_1 = instance.dist_matrix[a.id][c.id] + instance.dist_matrix[d.id][e.id] +
                                      instance.dist_matrix[f.id][b.id] - instance.dist_matrix[a.id][b.id] -
                                      instance.dist_matrix[c.id][d.id] - instance.dist_matrix[e.id][f.id];
                            isFeasible = true;
                            doSwapOrRelocateOrInvert = 0;
                            if(delta_1 < min_delta) {
                                min_delta = delta_1;
                                min_delta_name = "delta_1";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, i, 0, customer2_id, 0,
                                              customer2_id + 1, 1, customer3_id, 0, customer3_id + 1, 1, i + 1, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        /*
                        Route 1 = [* * A B * *]  demand_1 = [0 - A]  demand_2 = [B - 0]
                        Route 2 = [* * C D * *]  demand_3 = [0 - C]  demand_4 = [D - 0]
                        Route 3 = [* * E F * *]  demand_5 = [0 - E]  demand_6 = [F - 0]
                        */
                        // Option 2
                        // 1-2, 2-3, 3-1
                        // 0-A-D-0, 0-C-E-0, 0-F-B-0
                        route1_new_demand = demand_1 + demand_4;
                        route2_new_demand = demand_3 + demand_5;
                        route3_new_demand = demand_6 + demand_2;
                        delta_2 = 1;

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q) &&
                           (route3_new_demand <= vehicle3->Q)) {
                            delta_2 = instance.dist_matrix[a.id][d.id] + instance.dist_matrix[c.id][e.id] +
                                      instance.dist_matrix[f.id][b.id] - instance.dist_matrix[a.id][b.id] -
                                      instance.dist_matrix[c.id][d.id] - instance.dist_matrix[e.id][f.id];
                            isFeasible = true;
                            doSwapOrRelocateOrInvert = 0;
                            if(delta_2 < min_delta) {
                                min_delta = delta_2;
                                min_delta_name = "delta_2";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, i, 0, customer2_id + 1, 1,
                                              customer2_id, 0, customer3_id, 0, customer3_id + 1, 1, i + 1, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        /*
                        Route 1 = [* * A B * *]  demand_1 = [0 - A]  demand_2 = [B - 0]
                        Route 2 = [* * C D * *]  demand_3 = [0 - C]  demand_4 = [D - 0]
                        Route 3 = [* * E F * *]  demand_5 = [0 - E]  demand_6 = [F - 0]
                        */
                        // Option 3
                        // 1-2, 2-3, 3-1
                        // 0-A-D-0, 0-C-F-0, 0-E-B-0
                        route1_new_demand = demand_1 + demand_4;
                        route2_new_demand = demand_3 + demand_6;
                        route3_new_demand = demand_5 + demand_2;
                        delta_3 = 1;

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q) &&
                           (route3_new_demand <= vehicle3->Q)) {
                            delta_3 = instance.dist_matrix[a.id][d.id] + instance.dist_matrix[c.id][f.id] +
                                      instance.dist_matrix[e.id][b.id] - instance.dist_matrix[a.id][b.id] -
                                      instance.dist_matrix[c.id][d.id] - instance.dist_matrix[e.id][f.id];
                            isFeasible = true;
                            doSwapOrRelocateOrInvert = 0;
                            if(delta_3 < min_delta) {
                                min_delta = delta_3;
                                min_delta_name = "delta_3";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, i, 0, customer2_id + 1, 1,
                                              customer2_id, 0, customer3_id + 1, 1, customer3_id, 0, i + 1, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }
                        /*
                        Route 1 = [* * A B * *]  demand_1 = [0 - A]  demand_2 = [B - 0]
                        Route 2 = [* * C D * *]  demand_3 = [0 - C]  demand_4 = [D - 0]
                        Route 3 = [* * E F * *]  demand_5 = [0 - E]  demand_6 = [F - 0]
                          */
                        // Option 4
                        // 1-2, 2-3, 3-1
                        // 0-A-C-0, 0-D-F-0, 0-E-B-0
                        route1_new_demand = demand_1 + demand_3;
                        route2_new_demand = demand_4 + demand_6;
                        route3_new_demand = demand_5 + demand_2;
                        delta_4 = 1;

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q) &&
                           (route3_new_demand <= vehicle3->Q)) {
                            delta_4 = instance.dist_matrix[a.id][c.id] + instance.dist_matrix[d.id][f.id] +
                                      instance.dist_matrix[e.id][b.id] - instance.dist_matrix[a.id][b.id] -
                                      instance.dist_matrix[c.id][d.id] - instance.dist_matrix[e.id][f.id];
                            isFeasible = true;
                            doSwapOrRelocateOrInvert = 0;
                            if(delta_4 < min_delta) {
                                min_delta = delta_4;
                                min_delta_name = "delta_4";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle2_idx, vehicle3_idx, i, 0, customer2_id, 0,
                                              customer2_id + 1, 1, customer3_id + 1, 1, customer3_id, 0, i + 1, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }
                        /*
                        Route 1 = [* * A B * *]  demand_1 = [0 - A]  demand_2 = [B - 0]
                        Route 2 = [* * C D * *]  demand_3 = [0 - C]  demand_4 = [D - 0]
                        Route 3 = [* * E F * *]  demand_5 = [0 - E]  demand_6 = [F - 0]
                          */
                        // Option 5
                        // 1-3, 3-2, 2-1
                        // 0-A-E-0, 0-F-C-0, 0-D-B-0
                        route1_new_demand = demand_1 + demand_5;
                        route2_new_demand = demand_6 + demand_3;
                        route3_new_demand = demand_4 + demand_2;
                        delta_5 = 1;

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q) &&
                           (route3_new_demand <= vehicle3->Q)) {
                            delta_5 = instance.dist_matrix[a.id][e.id] + instance.dist_matrix[f.id][c.id] +
                                      instance.dist_matrix[d.id][b.id] - instance.dist_matrix[a.id][b.id] -
                                      instance.dist_matrix[c.id][d.id] - instance.dist_matrix[e.id][f.id];
                            isFeasible = true;
                            doSwapOrRelocateOrInvert = 0;
                            if(delta_5 < min_delta) {
                                min_delta = delta_5;
                                min_delta_name = "delta_5";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle3_idx, vehicle2_idx, i, 0, customer3_id, 0,
                                              customer3_id + 1, 1, customer2_id, 0, customer2_id + 1, 1, i + 1, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }
                        /*
                                            Route 1 = [* * A B * *]  demand_1 = [0 - A]  demand_2 = [B - 0]
                                            Route 2 = [* * C D * *]  demand_3 = [0 - C]  demand_4 = [D - 0]
                                            Route 3 = [* * E F * *]  demand_5 = [0 - E]  demand_6 = [F - 0]
                                              */
                        // Option 6
                        // 1-3, 3-2, 2-1
                        // 0-A-E-0, 0-F-D-0, 0-C-B-0
                        route1_new_demand = demand_1 + demand_5;
                        route2_new_demand = demand_6 + demand_4;
                        route3_new_demand = demand_3 + demand_2;
                        delta_6 = 1;

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q) &&
                           (route3_new_demand <= vehicle3->Q)) {
                            delta_6 = instance.dist_matrix[a.id][e.id] + instance.dist_matrix[f.id][d.id] +
                                      instance.dist_matrix[c.id][b.id] - instance.dist_matrix[a.id][b.id] -
                                      instance.dist_matrix[c.id][d.id] - instance.dist_matrix[e.id][f.id];
                            isFeasible = true;
                            doSwapOrRelocateOrInvert = 0;
                            if(delta_6 < min_delta) {
                                min_delta = delta_6;
                                min_delta_name = "delta_6";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle3_idx, vehicle2_idx, i, 0, customer3_id, 0,
                                              customer3_id + 1, 1, customer2_id + 1, 1, customer2_id, 0, i + 1, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }
                        /*
                        Route 1 = [* * A B * *]  demand_1 = [0 - A]  demand_2 = [B - 0]
                        Route 2 = [* * C D * *]  demand_3 = [0 - C]  demand_4 = [D - 0]
                        Route 3 = [* * E F * *]  demand_5 = [0 - E]  demand_6 = [F - 0]
                          */
                        // Option 7
                        // 1-3, 3-2, 2-1
                        // 0-A-F-0, 0-E-C-0, 0-D-B-0
                        route1_new_demand = demand_1 + demand_6;
                        route2_new_demand = demand_5 + demand_3;
                        route3_new_demand = demand_4 + demand_2;
                        delta_7 = 1;

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q) &&
                           (route3_new_demand <= vehicle3->Q)) {
                            delta_7 = instance.dist_matrix[a.id][f.id] + instance.dist_matrix[e.id][c.id] +
                                      instance.dist_matrix[d.id][b.id] - instance.dist_matrix[a.id][b.id] -
                                      instance.dist_matrix[c.id][d.id] - instance.dist_matrix[e.id][f.id];
                            isFeasible = true;
                            doSwapOrRelocateOrInvert = 0;
                            if(delta_7 < min_delta) {
                                min_delta = delta_7;
                                min_delta_name = "delta_7";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle3_idx, vehicle2_idx, i, 0, customer3_id + 1, 1,
                                              customer3_id, 0, customer2_id, 0, customer2_id + 1, 1, i + 1, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }
                        /*
                       Route 1 = [* * A B * *]  demand_1 = [0 - A]  demand_2 = [B - 0]
                       Route 2 = [* * C D * *]  demand_3 = [0 - C]  demand_4 = [D - 0]
                       Route 3 = [* * E F * *]  demand_5 = [0 - E]  demand_6 = [F - 0]
                         */
                        // Option 8
                        // 1-3, 3-2, 2-1
                        // 0-A-F-0, 0-E-D-0, 0-C-B-0
                        route1_new_demand = demand_1 + demand_6;
                        route2_new_demand = demand_5 + demand_4;
                        route3_new_demand = demand_3 + demand_2;
                        delta_8 = 1;

                        if((route1_new_demand <= vehicle1->Q) && (route2_new_demand <= vehicle2->Q) &&
                           (route3_new_demand <= vehicle3->Q)) {
                            delta_8 = instance.dist_matrix[a.id][f.id] + instance.dist_matrix[e.id][d.id] +
                                      instance.dist_matrix[c.id][b.id] - instance.dist_matrix[a.id][b.id] -
                                      instance.dist_matrix[c.id][d.id] - instance.dist_matrix[e.id][f.id];
                            isFeasible = true;
                            doSwapOrRelocateOrInvert = 0;
                            if(delta_8 < min_delta) {
                                min_delta = delta_8;
                                min_delta_name = "delta_8";
                                best_move.clear();
                                best_move = { vehicle1_idx, vehicle3_idx, vehicle2_idx, i, 0, customer3_id + 1, 1,
                                              customer3_id, 0, customer2_id + 1, 1, customer2_id, 0, i + 1, 1,
                                              doSwapOrRelocateOrInvert };
                                if(this->first_improvement)
                                    return this->apply_best_move(aux_solution, best_move);
                            }
                        }

                        if(isFeasible == false) {
                            continue;
                        }
                    }
                }
            }
        }
        vehicle1_idx++;
    }
    if(min_delta < 0) {
        /*
        cout << "min_delta " << min_delta << " C1 " << best_move[2] << " C2 " << best_move[3] << endl;
        float veh_1_cost = aux_solution->vehicles[best_move[0]]->compute_cost(instance);
        float veh_2_cost = aux_solution->vehicles[best_move[1]]->compute_cost(instance);
        cout << "cost_before " << veh_1_cost << " " << veh_2_cost << endl;
        */
        /*
        for(int i = 0; i < aux_solution->vehicles[best_move[0]]->customers.size(); i++) {
            cout << aux_solution->vehicles[best_move[0]]->customers[i]->id << " ";
        }
        cout << "  -----" << best_move[3] << "-" << best_move[4] << " " << best_move[5] << "-" << best_move[6] << endl;
        for(int i = 0; i < aux_solution->vehicles[best_move[1]]->customers.size(); i++) {
            cout << aux_solution->vehicles[best_move[1]]->customers[i]->id << " ";
        }
        cout << "  -----" << best_move[7] << "-" << best_move[8] << " " << best_move[9] << "-" << best_move[10] << endl;
        for(int i = 0; i < aux_solution->vehicles[best_move[2]]->customers.size(); i++) {
            cout << aux_solution->vehicles[best_move[2]]->customers[i]->id << " ";
        }
        cout << "  -----" << best_move[11] << "-" << best_move[12] << " " << best_move[13] << "-" << best_move[14]
             << endl;
        */
        // cout << min_delta << " before " << aux_solution->cost << endl;
        /*
        all_deltas = {
               delta_1,
               delta_2,
               delta_3,
               delta_4,
               delta_5,
               delta_6,
               delta_7,
               delta_8,
               delta_9,
               delta_10,
               delta_11,
               delta_12,
               delta_13,
               delta_14,
               delta_15,
               delta_16,
               delta_17,
               delta_18,
               delta_21,
               delta_22,
                  };
       vector<string> all_deltas_names = {
            "delta_1 ", "delta_2 ", "delta_3 ", "delta_4 ", "delta_5 ", "delta_6 ", "delta_7 ", "delta_8 ", "delta_9 ", "delta_10 ",
            "delta_11 ", "delta_12 ", "delta_13 ", "delta_14 ", "delta_15 ", "delta_16 ", "delta_17 ", "delta_18 ", "delta_21 ",
            "delta_22 ", "delta_23 ", "delta_24 ", "delta_25 ", "delta_26 ",
           };
           for (int i=0; i<all_deltas.size(); i++){
               if (all_deltas[i] == min_delta){
                   cout << all_deltas_names[i] << endl;
               }
           }
           */
        //cout << min_delta_name <<" " << min_delta << " before " << aux_solution->cost << endl;
        aux_solution = apply_best_move(aux_solution, best_move);
        //aux_solution->cost = aux_solution->cost + min_delta;
        aux_solution->cost = aux_solution->compute_cost(instance);
        //cout << "after " << aux_solution->cost << endl;
        // aux_solution->is_feasible();
        /*
        for(int i = 0; i < aux_solution->vehicles[best_move[0]]->customers.size(); i++) {
            cout << aux_solution->vehicles[best_move[0]]->customers[i]->id << " ";
        }
        cout << " " << endl;
        for(int i = 0; i < aux_solution->vehicles[best_move[1]]->customers.size(); i++) {
            cout << aux_solution->vehicles[best_move[1]]->customers[i]->id << " ";
        }
        cout << " " << endl;
        for(int i = 0; i < aux_solution->vehicles[best_move[2]]->customers.size(); i++) {
            cout << aux_solution->vehicles[best_move[2]]->customers[i]->id << " ";
        }
        cout << " " << endl << " ------------ " << endl;
        */
        /*
        veh_1_cost = aux_solution->vehicles[best_move[0]]->compute_cost(instance);
        veh_2_cost = aux_solution->vehicles[best_move[1]]->compute_cost(instance);
        cout << "cost_after " << veh_1_cost << " " << veh_2_cost << endl;
         * */
    }
    return aux_solution;
}

std::shared_ptr <Solution> ThreeOptInterNeighborhood::apply_best_move(std::shared_ptr <Solution> solution, vector<int>& best_move)
{ /*
 * best_move = { vehicle1_idx, vehicle3_idx, vehicle2_idx,
    i, 0, customer3_id + 1, 1, customer3_id, 0, customer2_id + 1, 1, customer2_id, 0, i + 1, 1 };
 * */
    if(best_move.size() > 11) {
        std::shared_ptr <Solution> local_optima = solution;
        int vehicle1_idx = best_move[0];
        int vehicle2_idx = best_move[1];
        int vehicle3_idx = best_move[2];
        int custom_1 = best_move[3];
        int direction_custom_1 = best_move[4];
        int custom_2 = best_move[5];
        int direction_custom_2 = best_move[6];
        int custom_3 = best_move[7];
        int direction_custom_3 = best_move[8];
        int custom_4 = best_move[9];
        int direction_custom_4 = best_move[10];
        int custom_5 = best_move[11];
        int direction_custom_5 = best_move[12];
        int custom_6 = best_move[13];
        int direction_custom_6 = best_move[14];
        int doSwapOrRelocateOrInvert = best_move[15];

        Vehicle& vehicle1 = *local_optima->vehicles[vehicle1_idx];
        Vehicle& vehicle2 = *local_optima->vehicles[vehicle2_idx];
        Vehicle& vehicle3 = *local_optima->vehicles[vehicle3_idx];

        vector<shared_ptr<Customer>> aux_vehicle1;
        vector<shared_ptr<Customer>> aux_vehicle2;
        vector<shared_ptr<Customer>> aux_vehicle3;
        vector<shared_ptr<Customer>> temp_vehicle1_customers(vehicle1.customers.size());
        if(doSwapOrRelocateOrInvert == 0) {
            // Start neew vehicle 1
            if(direction_custom_1 == 0) {
                for(int cust = 0; cust < custom_1 + 1; cust++) {
                    aux_vehicle1.push_back(vehicle1.customers[cust]);
                }
            } else {
                for(int cust = vehicle1.customers.size() - 1; cust > custom_1 - 1; cust--) {
                    aux_vehicle1.push_back(vehicle1.customers[cust]);
                }
            }
            if(direction_custom_2 == 0) {
                for(int cust = custom_2; cust > -1; cust--) {
                    aux_vehicle1.push_back(vehicle2.customers[cust]);
                }
            } else {
                for(int cust = custom_2; cust < vehicle2.customers.size(); cust++) {
                    aux_vehicle1.push_back(vehicle2.customers[cust]);
                }
            }

            // Start vehicle 2
            if(direction_custom_3 == 0) {
                for(int cust = 0; cust < custom_3 + 1; cust++) {
                    aux_vehicle2.push_back(vehicle2.customers[cust]);
                }
            } else {
                for(int cust = vehicle2.customers.size() - 1; cust > custom_3 - 1; cust--) {
                    aux_vehicle2.push_back(vehicle2.customers[cust]);
                }
            }
            if(direction_custom_4 == 0) {
                for(int cust = custom_4; cust > -1; cust--) {
                    aux_vehicle2.push_back(vehicle3.customers[cust]);
                }
            } else {
                for(int cust = custom_4; cust < vehicle3.customers.size(); cust++) {
                    aux_vehicle2.push_back(vehicle3.customers[cust]);
                }
            }

            // Start vehicle 3
            if(direction_custom_5 == 0) {
                for(int cust = 0; cust < custom_5 + 1; cust++) {
                    aux_vehicle3.push_back(vehicle3.customers[cust]);
                }
            } else {
                for(int cust = vehicle3.customers.size() - 1; cust > custom_5 - 1; cust--) {
                    aux_vehicle3.push_back(vehicle3.customers[cust]);
                }
            }
            if(direction_custom_6 == 0) {
                for(int cust = custom_6; cust > -1; cust--) {
                    aux_vehicle3.push_back(vehicle1.customers[cust]);
                }
            } else {
                for(int cust = custom_6; cust < vehicle1.customers.size(); cust++) {
                    aux_vehicle3.push_back(vehicle1.customers[cust]);
                }
            }

            // Actualizar vehiculo 1
            vehicle1.customers = vector<shared_ptr<Customer>>(aux_vehicle1.size());
            for(int cust = 0; cust < aux_vehicle1.size(); cust++) {
                vehicle1.customers[cust] = aux_vehicle1[cust];
                vehicle1.customers[cust]->vehicle_route = vehicle1_idx;
				local_optima->customers[vehicle1.customers[cust]->id]->vehicle_route = vehicle1_idx;
            }
            // Actualizar vehiculo 2
            vehicle2.customers = vector<shared_ptr<Customer>>(aux_vehicle2.size());
            for(int cust = 0; cust < aux_vehicle2.size(); cust++) {
                vehicle2.customers[cust] = aux_vehicle2[cust];
                vehicle2.customers[cust]->vehicle_route = vehicle2_idx;
				local_optima->customers[vehicle2.customers[cust]->id]->vehicle_route = vehicle2_idx;
            }
            // Actualizar vehiculo 3
            vehicle3.customers = vector<shared_ptr<Customer>>(aux_vehicle3.size());
            for(int cust = 0; cust < aux_vehicle3.size(); cust++) {
                vehicle3.customers[cust] = aux_vehicle3[cust];
                vehicle3.customers[cust]->vehicle_route = vehicle3_idx;
				local_optima->customers[vehicle3.customers[cust]->id]->vehicle_route = vehicle3_idx;
            }

            vehicle1.update_load();
            vehicle2.update_load();
            vehicle3.update_load();
        } else if(doSwapOrRelocateOrInvert == 1) {
            // Perform Swap first on Vehicle 1, before 2 opt
            // Swap first
            shared_ptr<Customer> customer_one = local_optima->vehicles[vehicle1_idx]->customers[custom_5];
            shared_ptr<Customer> customer_two = local_optima->vehicles[vehicle1_idx]->customers[custom_6];
            if(custom_6 > custom_5) {
                if(custom_6 == custom_5 + 1) {
                    local_optima->vehicles[vehicle1_idx]->insert_customer(customer_two, custom_5);
                    local_optima->vehicles[vehicle1_idx]->remove_customer(custom_6 + 1);
                } else {
                    local_optima->vehicles[vehicle1_idx]->insert_customer(customer_two, custom_5);
                    local_optima->vehicles[vehicle1_idx]->remove_customer(custom_6 + 1);

                    local_optima->vehicles[vehicle1_idx]->insert_customer(customer_one, custom_6 + 1);
                    local_optima->vehicles[vehicle1_idx]->remove_customer(custom_5 + 1);
                }
            } else {
                if(custom_5 == custom_6 + 1) {
                    local_optima->vehicles[vehicle1_idx]->insert_customer(customer_one, custom_6);
                    local_optima->vehicles[vehicle1_idx]->remove_customer(custom_5 + 1);
                } else {
                    local_optima->vehicles[vehicle1_idx]->insert_customer(customer_one, custom_6);
                    local_optima->vehicles[vehicle1_idx]->remove_customer(custom_5 + 1);

                    local_optima->vehicles[vehicle1_idx]->insert_customer(customer_two, custom_5 + 1);
                    local_optima->vehicles[vehicle1_idx]->remove_customer(custom_6 + 1);
                }
            }
            // End Swap. Now 2 opt.
            // Start neew vehicle 1
            vehicle1 = *local_optima->vehicles[vehicle1_idx];
            if(direction_custom_1 == 0) {
                for(int cust = 0; cust < custom_1 + 1; cust++) {
                    aux_vehicle1.push_back(vehicle1.customers[cust]);
                }
            } else {
                for(int cust = vehicle1.customers.size() - 1; cust > custom_1 - 1; cust--) {
                    aux_vehicle1.push_back(vehicle1.customers[cust]);
                }
            }
            if(direction_custom_2 == 0) {
                for(int cust = custom_2; cust > -1; cust--) {
                    aux_vehicle1.push_back(vehicle2.customers[cust]);
                }
            } else {
                for(int cust = custom_2; cust < vehicle2.customers.size(); cust++) {
                    aux_vehicle1.push_back(vehicle2.customers[cust]);
                }
            }
            // Start vehicle 2
            if(direction_custom_3 == 0) {
                for(int cust = 0; cust < custom_3 + 1; cust++) {
                    aux_vehicle2.push_back(vehicle2.customers[cust]);
                }
            } else {
                for(int cust = vehicle2.customers.size() - 1; cust > custom_3 - 1; cust--) {
                    aux_vehicle2.push_back(vehicle2.customers[cust]);
                }
            }
            if(direction_custom_4 == 0) {
                for(int cust = custom_4; cust > -1; cust--) {
                    aux_vehicle2.push_back(vehicle1.customers[cust]);
                }
            } else {
                for(int cust = custom_4; cust < vehicle1.customers.size(); cust++) {
                    aux_vehicle2.push_back(vehicle1.customers[cust]);
                }
            }
            // Actualizar vehiculo 1
            vehicle1.customers = vector<shared_ptr<Customer>>(aux_vehicle1.size());
            for(int cust = 0; cust < aux_vehicle1.size(); cust++) {
                vehicle1.customers[cust] = aux_vehicle1[cust];
                vehicle1.customers[cust]->vehicle_route = vehicle1_idx;
				local_optima->customers[vehicle1.customers[cust]->id]->vehicle_route = vehicle1_idx;
            }
            // Actualizar vehiculo 2
            vehicle2.customers = vector<shared_ptr<Customer>>(aux_vehicle2.size());
            for(int cust = 0; cust < aux_vehicle2.size(); cust++) {
                vehicle2.customers[cust] = aux_vehicle2[cust];
                vehicle2.customers[cust]->vehicle_route = vehicle2_idx;
				local_optima->customers[vehicle2.customers[cust]->id]->vehicle_route = vehicle2_idx;
            }
            vehicle1.update_load();
            vehicle2.update_load();
            //
        } else if(doSwapOrRelocateOrInvert == 2) {

            // Perform relocate first on Vehicle 1, before 2 opt
            shared_ptr<Customer> customer_one = local_optima->vehicles[vehicle1_idx]->customers[custom_5];
            shared_ptr<Customer> customer_two = local_optima->vehicles[vehicle1_idx]->customers[custom_6];
            if(custom_6 > custom_5) {
                if(custom_6 == custom_5 + 1) {
                    local_optima->vehicles[vehicle1_idx]->insert_customer(customer_two, custom_5);
                    local_optima->vehicles[vehicle1_idx]->remove_customer(custom_6 + 1);
                }
            } else {
                if(custom_5 == custom_6 + 1) {
                    local_optima->vehicles[vehicle1_idx]->insert_customer(customer_one, custom_6);
                    local_optima->vehicles[vehicle1_idx]->remove_customer(custom_5 + 1);
                }
            }
            // End Swap. Now 2 opt.
            // Start neew vehicle 1
            vehicle1 = *local_optima->vehicles[vehicle1_idx];
            if(direction_custom_1 == 0) {
                for(int cust = 0; cust < custom_1 + 1; cust++) {
                    aux_vehicle1.push_back(vehicle1.customers[cust]);
                }
            } else {
                for(int cust = vehicle1.customers.size() - 1; cust > custom_1 - 1; cust--) {
                    aux_vehicle1.push_back(vehicle1.customers[cust]);
                }
            }
            if(direction_custom_2 == 0) {
                for(int cust = custom_2; cust > -1; cust--) {
                    aux_vehicle1.push_back(vehicle2.customers[cust]);
                }
            } else {
                for(int cust = custom_2; cust < vehicle2.customers.size(); cust++) {
                    aux_vehicle1.push_back(vehicle2.customers[cust]);
                }
            }
            // Start vehicle 2
            if(direction_custom_3 == 0) {
                for(int cust = 0; cust < custom_3 + 1; cust++) {
                    aux_vehicle2.push_back(vehicle2.customers[cust]);
                }
            } else {
                for(int cust = vehicle2.customers.size() - 1; cust > custom_3 - 1; cust--) {
                    aux_vehicle2.push_back(vehicle2.customers[cust]);
                }
            }
            if(direction_custom_4 == 0) {
                for(int cust = custom_4; cust > -1; cust--) {
                    aux_vehicle2.push_back(vehicle1.customers[cust]);
                }
            } else {
                for(int cust = custom_4; cust < vehicle1.customers.size(); cust++) {
                    aux_vehicle2.push_back(vehicle1.customers[cust]);
                }
            }
            // Actualizar vehiculo 1
            vehicle1.customers = vector<shared_ptr<Customer>>(aux_vehicle1.size());
            for(int cust = 0; cust < aux_vehicle1.size(); cust++) {
                vehicle1.customers[cust] = aux_vehicle1[cust];
                vehicle1.customers[cust]->vehicle_route = vehicle1_idx;
				local_optima->customers[vehicle1.customers[cust]->id]->vehicle_route = vehicle1_idx;
            }
            // Actualizar vehiculo 2
            vehicle2.customers = vector<shared_ptr<Customer>>(aux_vehicle2.size());
            for(int cust = 0; cust < aux_vehicle2.size(); cust++) {
                vehicle2.customers[cust] = aux_vehicle2[cust];
                vehicle2.customers[cust]->vehicle_route = vehicle2_idx;
				local_optima->customers[vehicle2.customers[cust]->id]->vehicle_route = vehicle2_idx;
            }
            vehicle1.update_load();
            vehicle2.update_load();
        } else if(doSwapOrRelocateOrInvert == 3) {

            // Hacer inversion y luego 2-opt
            // Inversion
            for(int cust = 0; cust < custom_5; cust++) {
                aux_vehicle3.push_back(vehicle1.customers[cust]);
            }
            for(int cust = custom_6; cust > custom_5 - 1; cust--) {
                aux_vehicle3.push_back(vehicle1.customers[cust]);
            }
            for(int cust = custom_6 + 1; cust < vehicle1.customers.size(); cust++) {
                aux_vehicle3.push_back(vehicle1.customers[cust]);
            }
            // Actualizar vehiculo 1
            vehicle1.customers = vector<shared_ptr<Customer>>(aux_vehicle3.size());
            for(int cust = 0; cust < aux_vehicle3.size(); cust++) {
                vehicle1.customers[cust] = aux_vehicle3[cust];
                vehicle1.customers[cust]->vehicle_route = vehicle1_idx;
				local_optima->customers[vehicle1.customers[cust]->id]->vehicle_route = vehicle1_idx;
            }
            // Termina inversion
            // Start neew vehicle 1
            // vehicle1 = *local_optima->vehicles[vehicle1_idx];
            if(direction_custom_1 == 0) {
                for(int cust = 0; cust < custom_1 + 1; cust++) {
                    aux_vehicle1.push_back(vehicle1.customers[cust]);
                }
            } else {
                for(int cust = vehicle1.customers.size() - 1; cust > custom_1 - 1; cust--) {
                    aux_vehicle1.push_back(vehicle1.customers[cust]);
                }
            }
            if(direction_custom_2 == 0) {
                for(int cust = custom_2; cust > -1; cust--) {
                    aux_vehicle1.push_back(vehicle2.customers[cust]);
                }
            } else {
                for(int cust = custom_2; cust < vehicle2.customers.size(); cust++) {
                    aux_vehicle1.push_back(vehicle2.customers[cust]);
                }
            }
            // Start vehicle 2
            if(direction_custom_3 == 0) {
                for(int cust = 0; cust < custom_3 + 1; cust++) {
                    aux_vehicle2.push_back(vehicle2.customers[cust]);
                }
            } else {
                for(int cust = vehicle2.customers.size() - 1; cust > custom_3 - 1; cust--) {
                    aux_vehicle2.push_back(vehicle2.customers[cust]);
                }
            }
            if(direction_custom_4 == 0) {
                for(int cust = custom_4; cust > -1; cust--) {
                    aux_vehicle2.push_back(vehicle1.customers[cust]);
                }
            } else {
                for(int cust = custom_4; cust < vehicle1.customers.size(); cust++) {
                    aux_vehicle2.push_back(vehicle1.customers[cust]);
                }
            }
            // Actualizar vehiculo 1
            vehicle1.customers = vector<shared_ptr<Customer>>(aux_vehicle1.size());
            for(int cust = 0; cust < aux_vehicle1.size(); cust++) {
                vehicle1.customers[cust] = aux_vehicle1[cust];
                vehicle1.customers[cust]->vehicle_route = vehicle1_idx;
				local_optima->customers[vehicle1.customers[cust]->id]->vehicle_route = vehicle1_idx;
            }
            // Actualizar vehiculo 2
            vehicle2.customers = vector<shared_ptr<Customer>>(aux_vehicle2.size());
            for(int cust = 0; cust < aux_vehicle2.size(); cust++) {
                vehicle2.customers[cust] = aux_vehicle2[cust];
                vehicle2.customers[cust]->vehicle_route = vehicle2_idx;
				local_optima->customers[vehicle2.customers[cust]->id]->vehicle_route = vehicle2_idx;
            }
            vehicle1.update_load();
            vehicle2.update_load();
        } else if(doSwapOrRelocateOrInvert == 4) {
            // cout << "-----Hacer movimiento diferente -------" << endl;
            // Start neew vehicle 1
            // vehicle1 = *local_optima->vehicles[vehicle1_idx];
            if(direction_custom_1 == 0) {
                for(int cust = 0; cust < custom_1 + 1; cust++) {
                    aux_vehicle1.push_back(vehicle2.customers[cust]);
                }
            } else {
                for(int cust = vehicle2.customers.size() - 1; cust > custom_1 - 1; cust--) {
                    aux_vehicle1.push_back(vehicle2.customers[cust]);
                }
            }
            for(int cust = custom_2; cust < custom_3 + 1; cust++) {
                aux_vehicle1.push_back(vehicle1.customers[cust]);
            }
            if(direction_custom_4 == 0) {
                for(int cust = custom_4; cust > -1; cust--) {
                    aux_vehicle1.push_back(vehicle2.customers[cust]);
                }
            } else {
                for(int cust = custom_4; cust < vehicle2.customers.size(); cust++) {
                    aux_vehicle1.push_back(vehicle2.customers[cust]);
                }
            }
            // Start vehicle 2
            for(int cust = 0; cust < custom_5 + 1; cust++) {
                aux_vehicle2.push_back(vehicle1.customers[cust]);
            }
            for(int cust = custom_6; cust < vehicle1.customers.size(); cust++) {
                aux_vehicle2.push_back(vehicle1.customers[cust]);
            }
            // Actualizar vehiculo 1
            vehicle1.customers = vector<shared_ptr<Customer>>(aux_vehicle1.size());
            for(int cust = 0; cust < aux_vehicle1.size(); cust++) {
                vehicle1.customers[cust] = aux_vehicle1[cust];
                vehicle1.customers[cust]->vehicle_route = vehicle1_idx;
				local_optima->customers[vehicle1.customers[cust]->id]->vehicle_route = vehicle1_idx;
            }
            // Actualizar vehiculo 2
            vehicle2.customers = vector<shared_ptr<Customer>>(aux_vehicle2.size());
            for(int cust = 0; cust < aux_vehicle2.size(); cust++) {
                vehicle2.customers[cust] = aux_vehicle2[cust];
                vehicle2.customers[cust]->vehicle_route = vehicle2_idx;
				local_optima->customers[vehicle2.customers[cust]->id]->vehicle_route = vehicle2_idx;
            }
            vehicle1.update_load();
            vehicle2.update_load();
        }
        return local_optima;
    } else if(best_move.size() == 11) {
        std::shared_ptr <Solution> local_optima = solution;
        int vehicle1_idx = best_move[0];
        int vehicle2_idx = best_move[1];
        int vehicle3_idx = best_move[2];
        int custom_1 = best_move[3];
        int direction_custom_1 = best_move[4];
        int custom_2 = best_move[5];
        int direction_custom_2 = best_move[6];
        int custom_3 = best_move[7];
        int direction_custom_3 = best_move[8];
        int custom_4 = best_move[9];
        int direction_custom_4 = best_move[10];

        Vehicle& vehicle1 = *local_optima->vehicles[vehicle1_idx];
        Vehicle& vehicle2 = *local_optima->vehicles[vehicle2_idx];
        // Vehicle& vehicle3 = *local_optima->vehicles[vehicle3_idx];

        vector<shared_ptr<Customer>> aux_vehicle1;
        vector<shared_ptr<Customer>> aux_vehicle2;
        // vector<shared_ptr<Customer>> aux_vehicle3;
        vector<shared_ptr<Customer>> temp_vehicle1_customers(vehicle1.customers.size());

        // Start neew vehicle 1
        if(direction_custom_1 == 0) {
            for(int cust = 0; cust < custom_1 + 1; cust++) {
                aux_vehicle1.push_back(vehicle1.customers[cust]);
            }
        } else {
            for(int cust = vehicle1.customers.size() - 1; cust > custom_1 - 1; cust--) {
                aux_vehicle1.push_back(vehicle1.customers[cust]);
            }
        }
        if(direction_custom_2 == 0) {
            for(int cust = custom_2; cust > -1; cust--) {
                aux_vehicle1.push_back(vehicle2.customers[cust]);
            }
        } else {
            for(int cust = custom_2; cust < vehicle2.customers.size(); cust++) {
                aux_vehicle1.push_back(vehicle2.customers[cust]);
            }
        }

        // Start vehicle 2
        if(direction_custom_3 == 0) {
            for(int cust = 0; cust < custom_3 + 1; cust++) {
                aux_vehicle2.push_back(vehicle2.customers[cust]);
            }
        } else {
            for(int cust = vehicle2.customers.size() - 1; cust > custom_3 - 1; cust--) {
                aux_vehicle2.push_back(vehicle2.customers[cust]);
            }
        }
        if(direction_custom_4 == 0) {
            for(int cust = custom_4; cust > -1; cust--) {
                aux_vehicle2.push_back(vehicle1.customers[cust]);
            }
        } else {
            for(int cust = custom_4; cust < vehicle1.customers.size(); cust++) {
                aux_vehicle2.push_back(vehicle1.customers[cust]);
            }
        }
        // Actualizar vehiculo 1
        vehicle1.customers = vector<shared_ptr<Customer>>(aux_vehicle1.size());
        for(int cust = 0; cust < aux_vehicle1.size(); cust++) {
            vehicle1.customers[cust] = aux_vehicle1[cust];
            vehicle1.customers[cust]->vehicle_route = vehicle1_idx;
			local_optima->customers[vehicle1.customers[cust]->id]->vehicle_route = vehicle1_idx;
        }
        // Actualizar vehiculo 2
        vehicle2.customers = vector<shared_ptr<Customer>>(aux_vehicle2.size());
        for(int cust = 0; cust < aux_vehicle2.size(); cust++) {
            vehicle2.customers[cust] = aux_vehicle2[cust];
            vehicle2.customers[cust]->vehicle_route = vehicle2_idx;
			local_optima->customers[vehicle2.customers[cust]->id]->vehicle_route = vehicle2_idx;
        }
        vehicle1.update_load();
        vehicle2.update_load();
        return local_optima;
    } else {
        return solution;
    }
}

struct my_cross_exchange_comparator {
    // queue elements are vectors so we need to compare those
    bool operator()(vector<int> const& a, vector<int> const& b) const
    {
        // reverse sort puts the lowest value at the top
        return a[10] > b[10];
    }
};

using cross_exchange_priority_queue = priority_queue<vector<int>, vector<vector<int>>, my_cross_exchange_comparator>;

CrossExchangeNeighborhood::CrossExchangeNeighborhood()
{
}

CrossExchangeNeighborhood::CrossExchangeNeighborhood(bool first_improvement /*= false*/)
        : Neighborhood(first_improvement, "F")
{
}

CrossExchangeNeighborhood::~CrossExchangeNeighborhood()
{
}

std::shared_ptr <Solution> CrossExchangeNeighborhood::search_neighbors(std::shared_ptr <Solution> aux_solution, const Instance& instance)
{

    // //cout << "Entra cross exchange" << endl;

    // int x_nearest = aux_solution->customers.size() / 30;
    int x_nearest = aux_solution->customers.size() / 1;

    aux_solution->compute_cummulative_demands();

    int min_delta = 0;
    vector<int> best_move;
    vector<vector<int>> possible_moves_c1;
    cross_exchange_priority_queue possible_moves;
    vector<int> next_possible_move;
    vector<int> c1_possible;
    // customers
    Customer a;
    Customer b;
    Customer c;
    Customer d;
    Customer e;
    // demands
    int route1_new_demand;
    int route2_new_demand;
    // delta
    float delta;
    int vehicle1_idx = 0;
    int vehicle2_idx = 0;
    int cost_1_1 = 1;
    int cost_1_2 = 1;
    int cost_1_star = 1;
    // Vehicle& vehicle2 = *aux_solution->vehicles[vehicle2_idx];

    for(auto& vehicle1 : aux_solution->vehicles) {
        for(int i = 1; i < vehicle1->customers.size() - 1; i++) {
            a = *vehicle1->customers[i];
            b = *vehicle1->customers[i + 1];
            bool is_there = false;
            // //cout << "llega hasta 1 " << endl;
            for(int j = 0; j < x_nearest; j++) {

                int c_id = instance.sorted_dist_matrix[a.id][j];
                vehicle2_idx = aux_solution->customers[c_id]->vehicle_route;
				if((vehicle1_idx == vehicle2_idx)||(vehicle2_idx == aux_solution->vehicles.size())) {
					continue;
				}
                // Vehicle& vehicle2 = *aux_solution->vehicles[vehicle2_idx];
                auto& vehicle2 = aux_solution->vehicles[vehicle2_idx];

                if(vehicle1_idx != vehicle2_idx) {

                    int customer2_id = -1;
                    if(vehicle2->customers.size() <= 2) {
                        continue;
                    }
                    // //cout << "vehicle2->customers.size  " << vehicle2->customers.size() <<  endl;
                    for(int k = 1; k < vehicle2->customers.size() - 1; k++) {
                        // //cout << "k:  " << k <<  endl;
                        if(vehicle2->customers[k]->id == c_id) {
                            c = *vehicle2->customers[k];
                            d = *vehicle2->customers[k + 1];
                            e = *vehicle2->customers[k - 1];
                            customer2_id = k;
                            // //cout << "llega hasta 4 " << endl;
                            break;
                        }
                    }
                    if(customer2_id == -1) {
                        continue;
                    }

                    // Eval both c1 possible start moves joining ik+1 with  jl+1 or jl-1.
                    cost_1_1 = instance.dist_matrix[a.id][c.id] + instance.dist_matrix[b.id][e.id] -
                               instance.dist_matrix[a.id][b.id] - instance.dist_matrix[c.id][e.id];

                    cost_1_2 = instance.dist_matrix[a.id][c.id] + instance.dist_matrix[b.id][d.id] -
                               instance.dist_matrix[a.id][b.id] - instance.dist_matrix[c.id][d.id];

                    cost_1_star = instance.dist_matrix[a.id][c.id] - instance.dist_matrix[a.id][b.id];

                    if(cost_1_1 <= 0) {
                        // Check if route constraints are not violated before storing candidate move
                        route1_new_demand = vehicle1->get_demand_between(0, i) +
                                            vehicle2->get_demand_between(customer2_id, vehicle2->customers.size() - 1);
                        // otra opcion que puede ser mas rapida (menos operaciones de consulta) para calcular demand de
                        // ruta 2  route2_new_demand = vehicle1.load + vehicle2->load - route1_new_demand;
                        route2_new_demand = vehicle1->get_demand_between(i + 1, vehicle1->customers.size() - 1) +
                                            vehicle2->get_demand_between(0, customer2_id - 1);

                        if(route1_new_demand > vehicle1->Q || route2_new_demand > vehicle2->Q) {
                            continue;
                        }
                        c1_possible.reserve(8);
                        c1_possible.push_back(cost_1_1);
                        c1_possible.push_back(vehicle1_idx);
                        c1_possible.push_back(vehicle2_idx);
                        c1_possible.push_back(i);                // Ik position
                        c1_possible.push_back(customer2_id);     // Jl position
                        c1_possible.push_back(i + 1);            // Ik+1 position
                        c1_possible.push_back(customer2_id - 1); // Jl-1 position
                        c1_possible.push_back(1);
                        // c1_possible.push_back(a.id);
                        // c1_possible.push_back(b.id);
                        // c1_possible.push_back(c.id);
                        // c1_possible.push_back(e.id);
                        possible_moves_c1.push_back(c1_possible);
                        c1_possible.clear();
                    }

                    if(cost_1_2 <= 0) {
                        // Check if route constraints are not violated before storing candidate move
                        route1_new_demand =
                                vehicle1->get_demand_between(0, i) + vehicle2->get_demand_between(0, customer2_id);
                        route2_new_demand = vehicle1->get_demand_between(i + 1, vehicle1->customers.size() - 1) +
                                            vehicle2->get_demand_between(customer2_id + 1, vehicle2->customers.size() - 1);

                        if(route1_new_demand > vehicle1->Q || route2_new_demand > vehicle2->Q) {
                            continue;
                        }
                        c1_possible.reserve(8);
                        c1_possible.push_back(cost_1_2);
                        c1_possible.push_back(vehicle1_idx);
                        c1_possible.push_back(vehicle2_idx);
                        c1_possible.push_back(i);                // Ik position
                        c1_possible.push_back(customer2_id);     // Jl position
                        c1_possible.push_back(i + 1);            // Ik+1 position
                        c1_possible.push_back(customer2_id + 1); // Jl+1 position
                        c1_possible.push_back(2);
                        // c1_possible.push_back(a.id);
                        // c1_possible.push_back(b.id);
                        // c1_possible.push_back(c.id);
                        // c1_possible.push_back(d.id);
                        possible_moves_c1.push_back(c1_possible);
                        c1_possible.clear();
                    }

                    if(cost_1_star <= 0) {
                        // Check if route constraints are not violated before storing candidate move

                        route1_new_demand = vehicle1->load + c.demand;
                        route2_new_demand = vehicle2->load - c.demand;

                        if(route1_new_demand > vehicle1->Q || route2_new_demand > vehicle2->Q) {
                            continue;
                        }

                        c1_possible.reserve(8);
                        c1_possible.push_back(cost_1_star);
                        c1_possible.push_back(vehicle1_idx);
                        c1_possible.push_back(vehicle2_idx);
                        c1_possible.push_back(i);                // Ik position
                        c1_possible.push_back(customer2_id);     // Jl position
                        c1_possible.push_back(i + 1);            // Ik+1 position
                        c1_possible.push_back(customer2_id + 1); // Jl+1 position
                        c1_possible.push_back(-1);
                        // c1_possible.push_back(a.id);
                        // c1_possible.push_back(b.id);
                        // c1_possible.push_back(c.id);
                        // c1_possible.push_back(d.id);
                        possible_moves_c1.push_back(c1_possible);
                        c1_possible.clear();
                    }

                }
            }

            // empieza analisis second cross
            //
            // //cout << "llega hasta 1 " << endl;
            // cout << "Possible moves c1 size: "<< possible_moves_c1.size() << endl;

            for(int move = 0; move < possible_moves_c1.size(); move++) {

                best_move.clear();
                best_move = possible_moves_c1[move];
                if(best_move[7] == -1) {
                    int vehicle1_idx = best_move[1];
                    int vehicle2_idx = best_move[2];
                    int i_k = best_move[3];
                    int j_l = best_move[4];
                    int i_kp1 = best_move[5];
                    // int j_lmp1 = best_move[6];

                    Vehicle& vehicle_1 = *aux_solution->vehicles[vehicle1_idx];
                    Vehicle& vehicle_2 = *aux_solution->vehicles[vehicle2_idx];

                    //for(int f = i_k; f < vehicle_1.customers.size() - 2; f++) {
                    // //cout << "llega hasta 3 " << endl;
                    for(int z = j_l; z < vehicle_2.customers.size() - 1; z++) {
                        // //cout << "llega hasta 4 " << endl;
                        int ik = vehicle_1.customers[i_k]->id; // Ik

                        int ikp1 = vehicle_1.customers[i_k + 1]->id; // Ik+1

                        int jl = vehicle_2.customers[j_l]->id; // Jl

                        int jlp1 = vehicle_2.customers[z]->id; // Jl+1

                        int jlp2 = vehicle_2.customers[z + 1]->id; // Jl+2

                        int jlm1 = vehicle_2.customers[j_l - 1]->id; // Jl-1

                        // if ((ikp1!=0) && (jlp1!=0) && (ikp2!=0) && (jl!=0) && (ik!=0) && (jlm1!=0) ){
                        // //cout << "llega hasta 5 " << endl;

                        //r_j = {j_l}
                        int c2_1 = instance.dist_matrix[ikp1][jlp1] + instance.dist_matrix[jlm1][jlp2] -
                                   instance.dist_matrix[jl][jlm1] - instance.dist_matrix[jlp1][jlp2];

                        //int c2_1 = instance.dist_matrix[ikp1][jlp1] + instance.dist_matrix[jl][ikp2] - \
                                instance.dist_matrix[ikp1][ikp2] - instance.dist_matrix[jl][jlp1];

                        //int c2_2 = instance.dist_matrix[ik][jl] + instance.dist_matrix[jlm1][ikp1] -\
                                instance.dist_matrix[ik][ikp1] - instance.dist_matrix[jlm1][jl];

                        int cost_to_check_1 = c2_1 + best_move[0];
                        //int cost_to_check_2 = c2_2 + best_move[0];

                        if(cost_to_check_1 < 0) {

                            // best_cost = cost_to_check_1;
                            // possible_move = { ID vehiculo 1, ID vehiculo 2, First cross Ik, First cross Jl,
                            // First cross Ik+1, First cross Jl+1 o Jl-1, Second cross Ik+1, Second cross
                            // jl+1,Second cross Ik+1, Second cross Jl }

                            possible_moves.push({ best_move[1], best_move[2], best_move[3], best_move[4],
                                                  best_move[5], best_move[6], ik, jl, jlp1, jlm1, cost_to_check_1,
                                                  best_move[7] }); // Add values to apply best_move
                            /*
                            if(this->first_improvement) {
                                return this->apply_best_move(aux_solution, possible_moves.top());
                            }
                             * */
                            //}
                        }
                    }
                    //}
                } else if(best_move[7] != -1) {
                    // std::shared_ptr <Solution> local_optima = aux_solution;
                    int vehicle1_idx = best_move[1];
                    int vehicle2_idx = best_move[2];
                    int i_k = best_move[3];
                    int j_l = best_move[4];
                    int i_kp1 = best_move[5];
                    int j_lmp1 = best_move[6];

                    Vehicle& vehicle_1 = *aux_solution->vehicles[vehicle1_idx];
                    Vehicle& vehicle_2 = *aux_solution->vehicles[vehicle2_idx];
                    vector<shared_ptr<Customer>> aux_vector;
                    vector<shared_ptr<Customer>> temp_vehicle1_customers(vehicle_1.customers.size());
                    vector<shared_ptr<Customer>> temp_vehicle1_1_customers(vehicle_1.customers.size());
                    vector<shared_ptr<Customer>> temp_vehicle2_customers(vehicle_2.customers.size());

                    // //cout << "vehicles before first cross " << endl;

                    for(int cust = 0; cust < vehicle_1.customers.size(); cust++) {
                        // //cout << vehicle_1.customers[cust]->id << " ";

                        temp_vehicle1_customers[cust] = vehicle_1.customers[cust];
                        temp_vehicle1_1_customers[cust] = vehicle_1.customers[cust];
                    }

                    // //cout << " " << endl;

                    for(int cust = 0; cust < vehicle_2.customers.size(); cust++) {
                        // //cout << vehicle_2.customers[cust]->id << " ";

                        temp_vehicle2_customers[cust] = vehicle_2.customers[cust];
                    }

                    // //cout << " " << endl;

                    // Starting first cross
                    for(int h = 0; h < i_k + 1; h++) {
                        aux_vector.push_back(vehicle_1.customers[h]);
                    }

                    if(j_l > j_lmp1) {
                        for(int h = j_l; h < vehicle_2.customers.size(); h++) {
                            aux_vector.push_back(vehicle_2.customers[h]);
                        }
                    } else {
                        for(int h = j_l; h > -1; h--) {
                            aux_vector.push_back(vehicle_2.customers[h]);
                        }
                    }

                    vehicle_1.customers = vector<shared_ptr<Customer>>(aux_vector.size());
                    for(int cust = 0; cust < aux_vector.size(); cust++) {

                        vehicle_1.customers[cust] = aux_vector[cust];
                    }

                    vehicle_1.update_load();
                    aux_vector.clear();

                    if(j_l > j_lmp1) {
                        for(int h = 0; h < j_lmp1 + 1; h++) {
                            aux_vector.push_back(vehicle_2.customers[h]);
                        }
                        for(int h = i_kp1; h < temp_vehicle1_customers.size(); h++) {
                            aux_vector.push_back(temp_vehicle1_customers[h]);
                        }
                    } else {
                        for(int h = temp_vehicle1_customers.size() - 1; h > i_k; h--) {
                            aux_vector.push_back(temp_vehicle1_customers[h]);
                        }
                        for(int h = j_lmp1; h < vehicle_2.customers.size(); h++) {
                            aux_vector.push_back(vehicle_2.customers[h]);
                        }
                    }

                    vehicle_2.customers = vector<shared_ptr<Customer>>(aux_vector.size());
                    for(int cust = 0; cust < aux_vector.size(); cust++) {

                        vehicle_2.customers[cust] = aux_vector[cust];
                    }
                    vehicle_2.update_load();
                    aux_vector.clear();

                    // Finishing First cross
                    // //cout << "Finishing First cross" << endl;

                    int best_cost = 0;
                    vector<int> c1_values = possible_moves_c1[move];
                    // if(c1_values[7] == 1) {
                    // Vehicle& vehicle_2 = vehicle2;
                    // Vehicle& vehicle_1 = vehicle1;
                    // //cout << "llega hasta 2 " << endl;
                    // if ((c1_values[3] != 0) && (c1_values[4] != 0)){
                    for(int f = c1_values[3]; f < vehicle_1.customers.size() - 3; f++) {
                        // //cout << "llega hasta 3 " << endl;
                        for(int z = c1_values[4]; z < vehicle_2.customers.size() - 2; z++) {
                            // //cout << "llega hasta 4 " << endl;
                            int ikp1 = vehicle_1.customers[f + 1]->id; // Ik+1

                            int jlp1 = vehicle_2.customers[z + 1]->id; // Jl+1

                            int jl = vehicle_2.customers[z]->id; // Jl

                            int ikp2 = vehicle_1.customers[f + 2]->id; // Ik+2

                            int ik = vehicle_1.customers[f]->id; // Ik

                            int jlm1 = vehicle_2.customers[z - 1]->id; // Jl-1

                            // if ((ikp1!=0) && (jlp1!=0) && (ikp2!=0) && (jl!=0) && (ik!=0) && (jlm1!=0) ){
                            // //cout << "llega hasta 5 " << endl;

                            int c2_1 = instance.dist_matrix[ikp1][jlp1] + instance.dist_matrix[jl][ikp2] -
                                       instance.dist_matrix[ikp1][ikp2] - instance.dist_matrix[jl][jlp1];

                            int c2_2 = instance.dist_matrix[ik][jl] + instance.dist_matrix[jlm1][ikp1] -
                                       instance.dist_matrix[ik][ikp1] - instance.dist_matrix[jlm1][jl];

                            int cost_to_check_1 = c2_1 + c1_values[0];
                            int cost_to_check_2 = c2_2 + c1_values[0];

                            if(cost_to_check_1 < 0) {

                                // best_cost = cost_to_check_1;
                                // possible_move = { ID vehiculo 1, ID vehiculo 2, First cross Ik, First cross
                                // Jl, First cross Ik+1, First cross Jl+1 o Jl-1, Second cross Ik+1, Second
                                // cross jl+1,Second cross Ik+1, Second cross Jl }

                                possible_moves.push({ c1_values[1], c1_values[2], c1_values[3], c1_values[4],
                                                      c1_values[5], c1_values[6], ikp1, jlp1, ikp2, jl, cost_to_check_1,
                                                      c1_values[7] }); // Add values to apply best_move
                                /*
                                if(this->first_improvement) {
                                    return this->apply_best_move(aux_solution, possible_moves.top());
                                }
                                 * */
                                //}
                            }

                            if(cost_to_check_2 < 0) {
                                // if(cost_to_check_2 < best_cost) {
                                // best_cost = cost_to_check_2;
                                // possible_move = { ID vehiculo 1, ID vehiculo 2, First cross Ik, First cross
                                // Jl, First cross Ik+1, First cross Jl+1 o Jl-1, Second cross Ik, Second cross
                                // jl,Second cross Ik+1, Second cross Jl-1 }

                                possible_moves.push({ c1_values[1], c1_values[2], c1_values[3], c1_values[4],
                                                      c1_values[5], c1_values[6], ik, jl, ikp1, jlm1, cost_to_check_2,
                                                      c1_values[7] }); // Add values to apply best_move
                                /*
                                if(this->first_improvement) {
                                    return this->apply_best_move(aux_solution, possible_moves.top());
                                }
                                */
                                //}
                            }
                            //}
                        }
                    }
                    //}
                    //} if c1_values[7]==1
                    vehicle_1.customers = vector<shared_ptr<Customer>>(temp_vehicle1_1_customers.size());
                    for(int cust = 0; cust < temp_vehicle1_1_customers.size(); cust++) {

                        vehicle_1.customers[cust] = temp_vehicle1_1_customers[cust];
                    }
                    vehicle_2.customers = vector<shared_ptr<Customer>>(temp_vehicle2_customers.size());
                    for(int cust = 0; cust < temp_vehicle2_customers.size(); cust++) {

                        vehicle_2.customers[cust] = temp_vehicle2_customers[cust];
                    }

                    vehicle_1.update_load();
                    vehicle_2.update_load();
                    // //cout << "termina primeros possibles second cross" << endl;

                    vector<shared_ptr<Customer>>().swap(aux_vector);
                    vector<shared_ptr<Customer>>().swap(temp_vehicle1_customers);
                    vector<shared_ptr<Customer>>().swap(temp_vehicle1_1_customers);
                    vector<shared_ptr<Customer>>().swap(temp_vehicle2_customers);

                }
            }
        }
        vehicle1_idx++;
    }
    // cout << "--------------- Final Possible moves size " << possible_moves.size() << endl;
    std::shared_ptr <Solution> final_solution;
    /*
    std::sort(possible_moves.begin(), possible_moves.end(),
        [](const std::vector<int>& a, const std::vector<int>& b) { return a[10] < b[10]; });
    */
    bool improved = false;
    // final_solution = aux_solution;
    for(int move = 0; move < possible_moves.size(); move++) {
        try {
            next_possible_move = possible_moves.top();
            possible_moves.pop();
            final_solution = this->apply_best_move(aux_solution, next_possible_move);
            // cout << "-------------Solution improved correctly----------------" << endl;
            // cout << " Applied C" << possible_moves[move][11] << endl;
            improved = true;
            vector<vector<int>>().swap(possible_moves_c1);
            // vector<vector<int>>().swap(possible_moves);
            vector<int>().swap(c1_possible);
            break;
        } catch(int error) {
            // cout << "saco throw" << endl;
        }
    }
    if(improved == true) {
        // cout << "entrega final solution" << endl;
        //cout << "----------- Cost before: " << final_solution->cost << endl;
        final_solution->cost = final_solution->compute_cost(instance);
        //cout << "-------------- n cost: " << final_solution->cost << endl;
        return final_solution;

    } else {
        // cout << "entrega aux solution" << endl;
        return aux_solution;
    }
}

std::shared_ptr <Solution> CrossExchangeNeighborhood::apply_best_move(std::shared_ptr <Solution> local_optima, vector<int>& best_move)
{
    if(best_move.size() != 0) {
        if (best_move[11] != -1) {
            // std::shared_ptr <Solution> local_optima = solution;
            int vehicle1_idx = best_move[0];
            int vehicle2_idx = best_move[1];
            int i_k = best_move[2];
            int j_l = best_move[3];
            int i_kp1 = best_move[4];
            int j_lmp1 = best_move[5];

            Vehicle& vehicle1 = *local_optima->vehicles[vehicle1_idx];
            Vehicle& vehicle2 = *local_optima->vehicles[vehicle2_idx];
            vector<shared_ptr<Customer>> aux_vector;
            vector<shared_ptr<Customer>> temp_vehicle1_customers(vehicle1.customers.size());
            vector<shared_ptr<Customer>> temp_vehicle1_1_customers(vehicle1.customers.size());
            vector<shared_ptr<Customer>> temp_vehicle2_customers(vehicle2.customers.size());

            // cout << "vehicles before first cross " << endl;

            for(int cust = 0; cust < vehicle1.customers.size(); cust++) {

                temp_vehicle1_customers[cust] = vehicle1.customers[cust];
                temp_vehicle1_1_customers[cust] = vehicle1.customers[cust];
            }

            // cout << " " << endl;

            for(int cust = 0; cust < vehicle2.customers.size(); cust++) {

                temp_vehicle2_customers[cust] = vehicle2.customers[cust];
            }

            // cout << " " << endl;
            /*//cout << "Customers de vehicle 1: " << endl;
            for(int cust = 0; cust < vehicle1.customers.size(); cust++){
                    //cout << vehicle1.customers[cust]->id << " ";
            }
            //cout << endl;
            //cout << "Customers de vehicle 2: " << endl;
            for(int cust = 0; cust < vehicle2.customers.size(); cust++){
                    //cout << vehicle2.customers[cust]->id << " ";
            }
            //cout << endl;*/

            // Starting first cross
            for(int h = 0; h < i_k + 1; h++) {
                aux_vector.push_back(vehicle1.customers[h]);
            }

            if(j_l > j_lmp1) {
                for(int h = j_l; h < vehicle2.customers.size(); h++) {
                    aux_vector.push_back(vehicle2.customers[h]);
                }
            } else {
                for(int h = j_l; h > -1; h--) {
                    aux_vector.push_back(vehicle2.customers[h]);
                }
            }

            vehicle1.customers = vector<shared_ptr<Customer>>(aux_vector.size());
            for(int cust = 0; cust < aux_vector.size(); cust++) {

                vehicle1.customers[cust] = aux_vector[cust];
            }

            vehicle1.update_load();
            aux_vector.clear();

            if(j_l > j_lmp1) {
                for(int h = 0; h < j_lmp1 + 1; h++) {
                    aux_vector.push_back(vehicle2.customers[h]);
                }
                for(int h = i_kp1; h < temp_vehicle1_customers.size(); h++) {
                    aux_vector.push_back(temp_vehicle1_customers[h]);
                }
            } else {
                for(int h = temp_vehicle1_customers.size() - 1; h > i_k; h--) {
                    aux_vector.push_back(temp_vehicle1_customers[h]);
                }
                for(int h = j_lmp1; h < vehicle2.customers.size(); h++) {
                    aux_vector.push_back(vehicle2.customers[h]);
                }
            }

            vehicle2.customers = vector<shared_ptr<Customer>>(aux_vector.size());
            for(int cust = 0; cust < aux_vector.size(); cust++) {

                vehicle2.customers[cust] = aux_vector[cust];
            }
            vehicle2.update_load();
            aux_vector.clear();

            // Finishing First cross

            // Starting Second cross

            int ik_11 = best_move[6];
            int jl_11 = best_move[7];
            int ik_21 = best_move[8];
            int jl_21 = best_move[9];
            int best_delta = best_move[10];
            int start = best_move[11];
            int ik_1 = -1;
            int jl_1 = -1;
            int ik_2 = -1;
            int jl_2 = -1;

            // cout << "Movements in first cross" << endl;
            // cout << i_k << " " << j_l << " " << i_kp1 << " " << j_lmp1 << " " << endl;
            // cout << "Starting second cross" << endl;
            // cout << "start = " << start << endl;
            // cout << ik_11 << " " << jl_11 << " " << ik_21 << " " << jl_21 << " " << endl;

            // cout << " before second cross " << endl;
            // temp_vehicle1_customers.clear();
            vector<shared_ptr<Customer>> temp_vehicle3_customers(vehicle1.customers.size());
            // cout << "Temp vehicle 1 " << endl;
            for(int cust = 0; cust < vehicle1.customers.size(); cust++) {

                temp_vehicle3_customers[cust] = vehicle1.customers[cust];
                // cout << temp_vehicle1_customers[cust]->id << " ";
            }

            // cout << " " << endl;
            // Esta seccion cambia, segun si el first cross era c1 o c2
            // if start == 2
            // cout << "old vehicle 1 " << endl;
            for(int f = 0; f < vehicle1.customers.size() - 1; f++) {
                // cout << vehicle1.customers[f]->id << " ";
                if(vehicle1.customers[f]->id == ik_11) {
                    ik_1 = f;
                }
                if(vehicle1.customers[f]->id == ik_21) {
                    ik_2 = f;
                }
            }
            // cout << " " << endl;
            // cout << "old vehicle 2 " << endl;
            for(int f = 0; f < vehicle2.customers.size() - 1; f++) {
                // cout << vehicle2.customers[f]->id << " ";
                if(vehicle2.customers[f]->id == jl_11) {
                    jl_1 = f;
                }
                if(vehicle2.customers[f]->id == jl_21) {
                    jl_2 = f;
                }
            }

            // cout << " " << endl;
            // cout << ik_1 << " " << jl_1 << " " << ik_2 << " " << jl_2 << " " << endl;

            if((ik_1 == -1) || (ik_2 == -1) || (jl_1 == -1) || (jl_2 == -1)) {
                // cout << "error" << endl;
                // return solution;
            } else {
                for(int h = 0; h < ik_1 + 1; h++) {
                    aux_vector.push_back(vehicle1.customers[h]);
                }
                for(int h = jl_1; h < vehicle2.customers.size(); h++) {
                    aux_vector.push_back(vehicle2.customers[h]);
                }
                // cout << "vehicles after second cross" << endl;
                vehicle1.customers = vector<shared_ptr<Customer>>(aux_vector.size());
                for(int cust = 0; cust < aux_vector.size(); cust++) {

                    vehicle1.customers[cust] = aux_vector[cust];
                    // cout << vehicle1.customers[cust]->id << " ";
                }
                // cout << " " << endl;

                vehicle1.update_load();
                aux_vector.clear();

                // cout << "aux vector 2_1" << endl;
                for(int h = 0; h < jl_2 + 1; h++) {
                    aux_vector.push_back(vehicle2.customers[h]);
                    // cout << aux_vector[h]->id << " ";
                }
                // //cout << "llega hasta 0" << endl;
                // cout << " " << endl;
                // cout << "aux vector 2_2" << endl;
                for(int h = ik_2; h < temp_vehicle3_customers.size(); h++) {
                    aux_vector.push_back(temp_vehicle3_customers[h]);
                    // cout << aux_vector[h]->id << " ";
                }
                // cout << " " << endl;

                vehicle2.customers = vector<shared_ptr<Customer>>(aux_vector.size());

                for(int cust = 0; cust < aux_vector.size(); cust++) {

                    vehicle2.customers[cust] = aux_vector[cust];

                }

                // cout << " " << endl;
                // //cout << "llega hasta 2" << endl;
                vehicle2.update_load();
                aux_vector.clear();
                // ends second cross
                // //cout << "llega hasta 3" << endl;

                if((vehicle1.load > vehicle1.Q) || (vehicle2.load > vehicle2.Q)) {
                    // cout << "Throw" << endl;
                    vehicle1.customers = vector<shared_ptr<Customer>>(temp_vehicle1_1_customers.size());
                    for(int cust = 0; cust < temp_vehicle1_1_customers.size(); cust++) {
                        vehicle1.customers[cust] = temp_vehicle1_1_customers[cust];
                    }
                    vehicle2.customers = vector<shared_ptr<Customer>>(temp_vehicle2_customers.size());
                    for(int cust = 0; cust < temp_vehicle2_customers.size(); cust++) {
                        vehicle2.customers[cust] = temp_vehicle2_customers[cust];
                    }
                    vehicle1.update_load();
                    vehicle2.update_load();
                    vector<shared_ptr<Customer>>().swap(aux_vector);
                    vector<shared_ptr<Customer>>().swap(temp_vehicle1_customers);
                    vector<shared_ptr<Customer>>().swap(temp_vehicle1_1_customers);
                    vector<shared_ptr<Customer>>().swap(temp_vehicle2_customers);
                    vector<shared_ptr<Customer>>().swap(temp_vehicle3_customers);

                    throw 1;
                }
                // cout << "entrega local optima" << endl;
                /*
                for(int cust = 0; cust < temp_vehicle1_1_customers.size(); cust++) {
                    cout << temp_vehicle1_1_customers[cust]->id << " ";
                }
                cout << " " << endl;
                for(int cust = 0; cust < temp_vehicle2_customers.size(); cust++) {
                    cout << temp_vehicle2_customers[cust]->id << " ";
                }
                cout << " " << endl;
                */
                vector<shared_ptr<Customer>>().swap(aux_vector);
                vector<shared_ptr<Customer>>().swap(temp_vehicle1_customers);
                vector<shared_ptr<Customer>>().swap(temp_vehicle1_1_customers);
                vector<shared_ptr<Customer>>().swap(temp_vehicle2_customers);
                vector<shared_ptr<Customer>>().swap(temp_vehicle3_customers);
                /*
                cout << "Expected improvement: " << best_delta << " Initial cost " << local_optima->cost << endl;
                ;
				*/
                for(int cust = 0; cust < vehicle1.customers.size(); cust++) {
                    vehicle1.customers[cust]->vehicle_route = vehicle1_idx;
					local_optima->customers[vehicle1.customers[cust]->id]->vehicle_route = vehicle1_idx;
                }            
                for(int cust = 0; cust < vehicle2.customers.size(); cust++) {
                    vehicle2.customers[cust]->vehicle_route = vehicle2_idx;
					local_optima->customers[vehicle2.customers[cust]->id]->vehicle_route = vehicle2_idx;
                }
              
                return local_optima;
            }
        }
        else if (best_move[11] == -1){
            //cout << "------------Maybe perform Or-exchange-------------------" << endl;
            int vehicle1_idx = best_move[0];
            int vehicle2_idx = best_move[1];
            int i_k = best_move[2];
            int j_l = best_move[3];
            int i_kp1 = best_move[4];
            int j_lmp1 = best_move[5];

            Vehicle& vehicle1 = *local_optima->vehicles[vehicle1_idx];
            Vehicle& vehicle2 = *local_optima->vehicles[vehicle2_idx];
            vector<shared_ptr<Customer>> aux_vector;
            vector<shared_ptr<Customer>> temp_vehicle1_customers(vehicle1.customers.size());
            vector<shared_ptr<Customer>> temp_vehicle1_1_customers(vehicle1.customers.size());
            vector<shared_ptr<Customer>> temp_vehicle2_customers(vehicle2.customers.size());

            // cout << "vehicles before first cross " << endl;

            for(int cust = 0; cust < vehicle1.customers.size(); cust++) {
                //      cout << vehicle1.customers[cust]->id << " ";
                temp_vehicle1_customers[cust] = vehicle1.customers[cust];
                temp_vehicle1_1_customers[cust] = vehicle1.customers[cust];
            }

            //  cout << " " << endl;

            for(int cust = 0; cust < vehicle2.customers.size(); cust++) {
                //      cout << vehicle2.customers[cust]->id << " ";
                temp_vehicle2_customers[cust] = vehicle2.customers[cust];
            }

            //  cout << " " << endl;
            /*//cout << "Customers de vehicle 1: " << endl;
            for(int cust = 0; cust < vehicle1.customers.size(); cust++){
                    //cout << vehicle1.customers[cust]->id << " ";
            }
            //cout << endl;
            //cout << "Customers de vehicle 2: " << endl;
            for(int cust = 0; cust < vehicle2.customers.size(); cust++){
                    //cout << vehicle2.customers[cust]->id << " ";
            }
            //cout << endl;*/
            // cout << "Starting first cross " << endl;
            // Starting first cross
            /*
            for(int h = 0; h < i_k + 1; h++) {
                aux_vector.push_back(vehicle1.customers[h]);
            }

            if(j_l > j_lmp1) {
                for(int h = j_l; h < vehicle2.customers.size(); h++) {
                    aux_vector.push_back(vehicle2.customers[h]);
                }
            } else {
                for(int h = j_l; h > -1; h--) {
                    aux_vector.push_back(vehicle2.customers[h]);
                }
            }

            vehicle1.customers = vector<shared_ptr<Customer>>(aux_vector.size());
            for(int cust = 0; cust < aux_vector.size(); cust++) {
                vehicle1.customers[cust] = aux_vector[cust];
            }

            vehicle1.update_load();
            aux_vector.clear();

            if(j_l > j_lmp1) {
                for(int h = 0; h < j_lmp1 + 1; h++) {
                    aux_vector.push_back(vehicle2.customers[h]);
                }
                for(int h = i_kp1; h < temp_vehicle1_customers.size(); h++) {
                    aux_vector.push_back(temp_vehicle1_customers[h]);
                }
            } else {
                for(int h = temp_vehicle1_customers.size() - 1; h > i_k; h--) {
                    aux_vector.push_back(temp_vehicle1_customers[h]);
                }
                for(int h = j_lmp1; h < vehicle2.customers.size(); h++) {
                    aux_vector.push_back(vehicle2.customers[h]);
                }
            }

            vehicle2.customers = vector<shared_ptr<Customer>>(aux_vector.size());
            for(int cust = 0; cust < aux_vector.size(); cust++) {
                vehicle2.customers[cust] = aux_vector[cust];
            }
            vehicle2.update_load();
            aux_vector.clear();
            cout << "Finishing first cross " << endl;
            */
            // Finishing First cross

            // Starting Second cross

            int ik_11 = best_move[6];
            int jl_11 = best_move[7];
            int ik_21 = best_move[8];
            int jl_21 = best_move[9];
            int best_delta = best_move[10];
            int start = best_move[11];
            int ik_1 = -1;
            int jl = -1;
            int jlp1 = -1;
            int jlm1 = -1;

            // cout << "Movements in first cross" << endl;
            // cout << i_k << " " << j_l << " " << i_kp1 << " " << j_lmp1 << " " << endl;
            // cout << "Starting second cross" << endl;
            // cout << "start = " << start << endl;
            // cout << ik_11 << " " << jl_11 << " " << ik_21 << " " << jl_21 << " " << endl;

            // cout << " before second cross " << endl;
            // temp_vehicle1_customers.clear();
            vector<shared_ptr<Customer>> temp_vehicle3_customers(vehicle1.customers.size());
            // cout << "Temp vehicle 1 " << endl;
            for(int cust = 0; cust < vehicle1.customers.size(); cust++) {
                temp_vehicle3_customers[cust] = vehicle1.customers[cust];
                //cout << temp_vehicle1_customers[cust]->id << " ";
            }

            // cout << " " << endl;
            // Esta seccion cambia, segun si el first cross era c1 o c2
            // if start == 2
            // cout << "old vehicle 1 " << endl;
            for(int f = 0; f < vehicle1.customers.size() - 1; f++) {
                // cout << vehicle1.customers[f]->id << " ";
                if(vehicle1.customers[f]->id == ik_11) {
                    ik_1 = f;
                }
            }
            // cout << " " << endl;
            // cout << "old vehicle 2 " << endl;
            for(int f = 0; f < vehicle2.customers.size() - 1; f++) {
                // cout << vehicle2.customers[f]->id << " ";
                if(vehicle2.customers[f]->id == ik_21) {
                    jlp1 = f;
                }
                if(vehicle2.customers[f]->id == jl_11) {
                    jl = f;
                }
                if(vehicle2.customers[f]->id == jl_21) {
                    jlm1 = f;
                }
            }

            // cout << " " << endl;
            // cout << "If one of these is -1 we get error: " << ik_1 << " " << jl << " " << jlp1 << " " << jlm1 << " " << endl;

            if((ik_1 == -1) || (jlp1 == -1) || (jl == -1) || (jlm1 == -1)) {
                // cout << "error" << endl;
            } else {
                for(int h = 0; h < ik_1 + 1; h++) {
                    aux_vector.push_back(vehicle1.customers[h]);
                }
                for(int h = jl; h < jlp1 + 1; h++) {
                    aux_vector.push_back(vehicle2.customers[h]);
                }
                for(int h = ik_1 +1; h < vehicle1.customers.size(); h++) {
                    aux_vector.push_back(vehicle1.customers[h]);
                }
                //  cout << "vehicles after second cross" << endl;
                vehicle1.customers = vector<shared_ptr<Customer>>(aux_vector.size());
                for(int cust = 0; cust < aux_vector.size(); cust++) {
                    vehicle1.customers[cust] = aux_vector[cust];
                    //  cout << vehicle1.customers[cust]->id << " ";
                }
                //  cout << " " << endl;

                vehicle1.update_load();
                aux_vector.clear();



                // cout << "aux vector 2_1" << endl;
                for(int h = 0; h < jlm1 + 1; h++) {
                    aux_vector.push_back(vehicle2.customers[h]);
                }
                for(int h = jlp1+1; h < vehicle2.customers.size(); h++) {
                    aux_vector.push_back(vehicle2.customers[h]);
                }
                // cout << " " << endl;

                vehicle2.customers = vector<shared_ptr<Customer>>(aux_vector.size());

                for(int cust = 0; cust < aux_vector.size(); cust++) {
                    vehicle2.customers[cust] = aux_vector[cust];
                    //  cout << vehicle2.customers[cust]->id << " ";
                }

                // cout << " " << endl;
                // //cout << "llega hasta 2" << endl;
                vehicle2.update_load();
                aux_vector.clear();
                // ends second cross
                // //cout << "llega hasta 3" << endl;

                if((vehicle1.load > vehicle1.Q) || (vehicle2.load > vehicle2.Q)) {
                    // cout << "Throw" << endl;
                    vehicle1.customers = vector<shared_ptr<Customer>>(temp_vehicle1_1_customers.size());
                    for(int cust = 0; cust < temp_vehicle1_1_customers.size(); cust++) {
                        vehicle1.customers[cust] = temp_vehicle1_1_customers[cust];
                    }
                    vehicle2.customers = vector<shared_ptr<Customer>>(temp_vehicle2_customers.size());
                    for(int cust = 0; cust < temp_vehicle2_customers.size(); cust++) {
                        vehicle2.customers[cust] = temp_vehicle2_customers[cust];
                    }
                    vehicle1.update_load();
                    vehicle2.update_load();
                    vector<shared_ptr<Customer>>().swap(aux_vector);
                    vector<shared_ptr<Customer>>().swap(temp_vehicle1_customers);
                    vector<shared_ptr<Customer>>().swap(temp_vehicle1_1_customers);
                    vector<shared_ptr<Customer>>().swap(temp_vehicle2_customers);
                    vector<shared_ptr<Customer>>().swap(temp_vehicle3_customers);
                    throw 1;
                }
                // cout << "entrega local optima" << endl;
                /*
                for(int cust = 0; cust < temp_vehicle1_1_customers.size(); cust++) {
                    cout << temp_vehicle1_1_customers[cust]->id << " ";
                }
                cout << " " << endl;
                for(int cust = 0; cust < temp_vehicle2_customers.size(); cust++) {
                    cout << temp_vehicle2_customers[cust]->id << " ";
                }
                cout << " " << endl;
                */
                vector<shared_ptr<Customer>>().swap(aux_vector);
                vector<shared_ptr<Customer>>().swap(temp_vehicle1_customers);
                vector<shared_ptr<Customer>>().swap(temp_vehicle1_1_customers);
                vector<shared_ptr<Customer>>().swap(temp_vehicle2_customers);
                vector<shared_ptr<Customer>>().swap(temp_vehicle3_customers);
                for(int cust = 0; cust < vehicle1.customers.size(); cust++) {
                    vehicle1.customers[cust]->vehicle_route = vehicle1_idx;
					local_optima->customers[vehicle1.customers[cust]->id]->vehicle_route = vehicle1_idx;
                }            
                for(int cust = 0; cust < vehicle2.customers.size(); cust++) {
                    vehicle2.customers[cust]->vehicle_route = vehicle2_idx;
					local_optima->customers[vehicle2.customers[cust]->id]->vehicle_route = vehicle2_idx;
                }
                return local_optima;
            }
        }
    } else {
        throw 1;
        return local_optima;
    }
}

RelocationChainNeighborhood_LimitNodes::RelocationChainNeighborhood_LimitNodes()
{
}

RelocationChainNeighborhood_LimitNodes::RelocationChainNeighborhood_LimitNodes(bool first_improvement /*= false*/)
        : Neighborhood(first_improvement, "B")
{
}

RelocationChainNeighborhood_LimitNodes::~RelocationChainNeighborhood_LimitNodes()
{
}


struct my_priority_queue_comparator {
    // queue elements are vectors so we need to compare those
    bool operator()(vector<int> const& a, vector<int> const& b) const
    {
        // reverse sort puts the lowest value at the top
        return a[0] > b[0];
    }
};

using my_priority_queue = priority_queue<vector<int>, vector<vector<int>>, my_priority_queue_comparator>;

my_priority_queue FirstRelocationFunction(std::shared_ptr <Solution> initial_solution, const Instance& instance, int initial_vehicle)
{
    int n_vehicles = initial_solution->vehicles.size();

    vector<int> best_move;
    vector<vector<int>> possible_moves;
    my_priority_queue possible_moves_q;
    vector<vector<int>> best_possible_moves;
    Vehicle vehicle_one;
    Vehicle vehicle_two;
    // customers
    Customer a;
    Customer b;
    Customer c;

    Customer f;
    Customer g;
    // demands
    int route1_new_demand;
    int route2_new_demand;
    int route1_old_demand;
    int route2_old_demand;
    // dist_matrix
    // delta
    int delta;
    int is_feasible = 0;

    //for(int vehicle_one_id = 0; vehicle_one_id < n_vehicles; vehicle_one_id++) {
    vehicle_one = *initial_solution->vehicles[initial_vehicle];
    int vehicle_one_id = initial_vehicle;
    int min_delta = 0;
    for(int customer_one_id = 1; customer_one_id < vehicle_one.customers.size() - 1; customer_one_id++) {
        for(int vehicle_two_id = 0; vehicle_two_id < n_vehicles; vehicle_two_id++) {
            // Avoid INTRAROUTE relocation
            if(vehicle_one_id == vehicle_two_id) {
                continue;
            }

            vehicle_two = *initial_solution->vehicles[vehicle_two_id];

            for(int customer_new_pos = 1; customer_new_pos < vehicle_two.customers.size() - 1; customer_new_pos++) {
                /*
                    route1 =  [A B C D]
                    route2 =  [E F G H]
                    Relocate(B, pos(G)) leaving us with:
                    new_route1 = [A C D]
                    new_route2 = [E F B G H]
                */

                a = *vehicle_one.customers[customer_one_id - 1];
                b = *vehicle_one.customers[customer_one_id];
                c = *vehicle_one.customers[customer_one_id + 1];

                f = *vehicle_two.customers[customer_new_pos - 1];
                g = *vehicle_two.customers[customer_new_pos];

                /*
                    Evaluate the improvement delta of both route1 and route2, using distance matrix
                    Delta = D(a,c) + D(f,b) + D(b,g) - D(a,b) -D(b,c) - D(f,g)

                delta = dist[a.id][c.id] + dist[f.id][b.id] + dist[b.id][g.id] \
                        - dist[a.id][b.id] - dist[b.id][c.id] - dist[f.id][g.id];
                delta = instance.insertion_matrix[f.id][g.id][b.id] + instance.remotion_matrix[a.id][c.id][b.id];
                */


                delta = instance.dist_matrix[a.id][c.id] + instance.dist_matrix[f.id][b.id] + instance.dist_matrix[b.id][g.id] \
                            - instance.dist_matrix[a.id][b.id] - instance.dist_matrix[b.id][c.id] - instance.dist_matrix[f.id][g.id];

                // Check if movement fixes the previous constraint of load violation, from first ejection

                route1_new_demand = vehicle_one.load - b.demand;

                if(route1_new_demand > vehicle_one.Q) {
                    continue;
                }

                // Store the best possible moves
                if(delta < 0) {
                    // min_delta = delta;
                    route2_new_demand = vehicle_two.load + b.demand;

                    if(route2_new_demand <= vehicle_two.Q) {
                        is_feasible = 1;
                    } else {
                        is_feasible = 0;
                    }

                    // possible_moves.push_back({ delta, vehicle_one_id, vehicle_two_id, b.id,
                    // customer_new_pos, is_feasible,
                    //      route1_new_demand, route2_new_demand, a.id, b.id, c.id, f.id, g.id, b.id, c.id,
                    //      g.id
                    //      });
                    possible_moves_q.push(
                            { delta, vehicle_one_id, vehicle_two_id, b.id, customer_new_pos, is_feasible,
                              route1_new_demand, route2_new_demand, a.id, b.id, c.id, f.id, g.id, b.id, c.id, g.id });
                    // Pensar en modificar aqui directamente f_suc, f_pred y veh_loads, para no repetir
                    // despues y ganar cpu time
                }
            }
        }
    }
    //}
    // std::sort(possible_moves.begin(), possible_moves.end(),
    //    [](const std::vector<int>& a, const std::vector<int>& b) { return a[0] < b[0]; });

    // cout << "Listo Sorting possible moves" << endl;
    /*
    if(possible_moves_q.size() > 0) {
        for(int i = 0; i < possible_moves_q.size(); i++) {
            //best_possible_moves.push_back(possible_moves[i]);
            vector<int> delta_q = possible_moves_q.top();
            best_possible_moves.push_back(delta_q);
            cout <<  " PQ delta: " << delta_q[0] << endl;
            possible_moves_q.pop();
            if(i >= 20) {
                break;
            }
        }
        */
    /*for(int i = 0; i < possible_moves.size(); i++) {
        if(possible_moves[i][5] == 1) {
            best_possible_moves.push_back(possible_moves[i]);
            break;
        }
    }*/
    /*
    int j = 0;
    int z = 2;
    for(int i = 1; i < possible_moves.size(); i++) {
        if((possible_moves[i][1] != best_possible_moves[j][1]) &&
            (possible_moves[i][2] != best_possible_moves[j][2])) {
            best_possible_moves.push_back(possible_moves[i]);
            j = j + 1;
            if(j >= 1) {
                z = i;
                break;
            }
        }
    }
    for(int i = z; i < possible_moves.size(); i++) {
        if((possible_moves[i][1] != best_possible_moves[j][1]) &&
            (possible_moves[i][2] != best_possible_moves[j][2]) &&
            (possible_moves[i][1] != best_possible_moves[j - 1][1]) &&
            (possible_moves[i][2] != best_possible_moves[j - 1][2])) {
            best_possible_moves.push_back(possible_moves[i]);
            j = j + 1;
            if(j >= 25) {
                break;
            }
        }
    }
}
*/
    return possible_moves_q;
}

struct Node {

    int delta;
    int is_feasible;
    vector<int> f_suc;
    vector<int> f_pred;
    unordered_map<int, int> vehicle_loads;
    vector<int> relocate_move;
    vector<Node*> child;
    Node* parent;
};

Node* newNode(vector<int> data, vector<int> n_f_suc, vector<int> n_f_pred, unordered_map<int, int> veh_loads)
{

    Node* temp = new Node;
    temp->delta = data[0];
    temp->relocate_move = { data[1], data[2], data[3], data[4], data[6] };
    temp->is_feasible = data[5];
    temp->f_suc = std::move(n_f_suc);
    temp->f_pred = std::move(n_f_pred);
    temp->vehicle_loads = std::move(veh_loads);
    return temp;
}

Node* InsertNode(Node* parent,
                 Node* ChildNode,
                 const vector<int>& data,
                 const vector<int>& n_f_suc,
                 const vector<int>& n_f_pred,
                 const unordered_map<int, int>& veh_loads)
{
    if(parent == NULL)
        ChildNode = newNode(data, n_f_suc, n_f_pred, veh_loads);
    else {
        Node* childNode = newNode(data, n_f_suc, n_f_pred, veh_loads);
        parent->child.push_back(childNode);
        childNode->parent = parent;
        return childNode;
    }
    return ChildNode;
}

void deleteNode(Node* node)
{
    for(int i = 0; i < node->child.size(); i++) {
        deleteNode(node->child[i]);
        delete node->child[i];
    }
}

my_priority_queue CompleteRelocationFunction(std::shared_ptr <Solution> initial_solution, const Instance& instance, Node* node)
{
    int n_vehicles = initial_solution->vehicles.size();

    vector<int> best_move;
    // vector<vector<int>> possible_moves;
    vector<vector<int>> best_possible_moves;
    my_priority_queue possible_moves;
    Vehicle vehicle_one;
    Vehicle vehicle_two;
    // customers
    Customer a;
    Customer b;
    Customer c;

    Customer f;
    Customer g;
    // demands
    int route1_new_demand;
    int route2_new_demand;
    int route1_old_demand;
    int route2_old_demand;
    // dist_matrix
    // delta
    int delta;
    int is_feasible = 0;
    int min_delta = 0;
    int vehicle_one_id = node->relocate_move[1];
    vector<int> f_suc = node->f_suc;
    vector<int> f_pred = node->f_pred;
    unordered_map<int, int> veh_loads = node->vehicle_loads;

    vehicle_one = *initial_solution->vehicles[vehicle_one_id];
    route1_old_demand = veh_loads[vehicle_one_id];

    for(int customer_one_id = 1; customer_one_id < vehicle_one.customers.size() - 1; customer_one_id++) {
        for(int vehicle_two_id = 0; vehicle_two_id < n_vehicles; vehicle_two_id++) {
            // Avoid INTRAROUTE relocation
            if(vehicle_one_id == vehicle_two_id) {
                continue;
            }

            vehicle_two = *initial_solution->vehicles[vehicle_two_id];

            // possible_moves.clear();
            for(int customer_new_pos = 1; customer_new_pos < vehicle_two.customers.size() - 1; customer_new_pos++) {
                /*
                    route1 =  [A B C D]
                    route2 =  [E F G H]
                    Relocate(B, pos(G)) leaving us with:
                    new_route1 = [A C D]
                    new_route2 = [E F B G H]
                */

                a = *vehicle_one.customers[customer_one_id - 1];
                b = *vehicle_one.customers[customer_one_id];
                c = *vehicle_one.customers[customer_one_id + 1];

                f = *vehicle_two.customers[customer_new_pos - 1];
                g = *vehicle_two.customers[customer_new_pos];

                // Check if B cant be relocated
                if(find(f_suc.begin(), f_suc.end(), b.id) != f_suc.end()) {
                    continue;
                }

                // Check if G cant be relocation destination
                // if(g.id != 0) {
                if(find(f_pred.begin(), f_pred.end(), g.id) != f_pred.end()) {
                    continue;
                }
                //}
                /*
                    Evaluate the improvement delta of both route1 and route2, using distance matrix
                    Delta = D(a,c) + D(f,b) + D(b,g) - D(a,b) -D(b,c) - D(f,g)

                delta = dist[a.id][c.id] + dist[f.id][b.id] + dist[b.id][g.id] \
                        - dist[a.id][b.id] - dist[b.id][c.id] - dist[f.id][g.id];
                delta = instance.insertion_matrix[f.id][g.id][b.id] + instance.remotion_matrix[a.id][c.id][b.id] + node->delta;
                */


                delta = instance.dist_matrix[a.id][c.id] + instance.dist_matrix[f.id][b.id] + instance.dist_matrix[b.id][g.id] \
                - instance.dist_matrix[a.id][b.id] - instance.dist_matrix[b.id][c.id] - instance.dist_matrix[f.id][g.id] +       node->delta;

                // Check if movement fixes the previous constraint of load violation, from first ejection

                route1_new_demand = route1_old_demand - b.demand;

                if(route1_new_demand > vehicle_one.Q) {
                    continue;
                }

                // Store the best possible moves
                if(delta < 0) {
                    // min_delta = delta;
                    route2_old_demand = veh_loads[vehicle_two_id];
                    route2_new_demand = route2_old_demand + b.demand;
                    if(route2_new_demand <= vehicle_two.Q) {
                        is_feasible = 1;
                    } else {
                        is_feasible = 0;
                    }

                    possible_moves.push({ delta, vehicle_one_id, vehicle_two_id, b.id, customer_new_pos, is_feasible,
                                          route1_new_demand, route2_new_demand, a.id, b.id, c.id, f.id, g.id, b.id, c.id, g.id });

                    // Pensar en modificar aqui directamente f_suc, f_pred y veh_loads, para no repetir despues
                    // y ganar cpu time
                }
            }
        }
    }
    /*
    std::sort(possible_moves.begin(), possible_moves.end(),
        [](const std::vector<int>& a, const std::vector<int>& b) { return a[0] < b[0]; });

    // cout << "Listo Sorting possible moves" << endl;

    if(possible_moves.size() > 0) {
        for(int i = 0; i < possible_moves.size(); i++) {
            best_possible_moves.push_back(possible_moves[i]);
            // if(i >= 15) {
            //   break;
            //}
        }
    }
    */
    // return best_possible_moves;
    return possible_moves;
}

std::shared_ptr <Solution> RelocationChainNeighborhood_LimitNodes::search_neighbors(std::shared_ptr <Solution> initial_solution, const Instance& instance)
{
    for(int k = 0; k < initial_solution->vehicles.size(); k++) {
        std::shared_ptr <Solution> local_optima = initial_solution;
        int n_vehicles = initial_solution->vehicles.size();

        vector<int> select_best_move;
        vector<int> move_data;
        vector<vector<int>> best_move;
        vector<int> possible_move;
        vector<int> first_best_move;
        vector<int> counters;
        // vector<vector<int>> first_possible_moves;
        my_priority_queue first_possible_moves;
        // vector<vector<int>> next_possible_moves;
        my_priority_queue next_possible_moves;
        vector<vector<int>> best_possible_moves;
        vector<unordered_map<int, int>> vehicle_loads;
        vector<vector<int>> f_suc;
        vector<vector<int>> f_pred;
        vector<vector<int>> delta_variation;
        vector<int> current_f_suc;
        vector<int> current_f_pred;
        unordered_map<int, int> initial_vehicle_loads;
        unordered_map<int, int> current_vehicle_loads;
        Node* best_move_node;
        Node* root;

        int level = 0;
        int final_move_feasible = 0;
        bool improved_solution = false;
        int cont_ext = 0;
        int current_counter = 0;
        int explored_nodes = 0;
        int max_first_nodes = 0;
        int max_explored_nodes = 100;


        for(int i = 0; i < initial_solution->vehicles.size(); i++) {
            initial_vehicle_loads[i] = initial_solution->vehicles[i]->load;
        }

        //cout << " " << k << " ";
        first_possible_moves = FirstRelocationFunction(initial_solution, instance, k);
        max_first_nodes = first_possible_moves.size();
        //cout << "First possible moves size: " << first_possible_moves.size() << endl;

        while ((improved_solution == false) && (cont_ext < max_first_nodes) &&
               (explored_nodes < max_explored_nodes + 1)) {
            explored_nodes++;
            // cout << "First possible move: " << cont_ext << endl;

            possible_move.clear();
            // possible_move = first_possible_moves[cont_ext];
            possible_move = first_possible_moves.top();
            first_possible_moves.pop();

            move_data.clear();

            move_data.push_back(possible_move[0]);  // Insert delta
            move_data.push_back(possible_move[1]);  // Insert vehicle_one_id
            move_data.push_back(possible_move[2]);  // Insert vehicle_two_id
            move_data.push_back(possible_move[3]);  // Insert client to relocate
            move_data.push_back(possible_move[4]);  // Insert relocation destination
            move_data.push_back(possible_move[5]);  // Insert is_feasible
            move_data.push_back(possible_move[12]); // Insert relocation destination id of client in position

            current_vehicle_loads = initial_vehicle_loads;
            current_vehicle_loads[possible_move[1]] = possible_move[6];
            current_vehicle_loads[possible_move[2]] = possible_move[7];

            current_f_suc.clear();
            current_f_pred.clear();

            current_f_suc.push_back(possible_move[8]);
            current_f_suc.push_back(possible_move[9]);
            current_f_suc.push_back(possible_move[10]);
            current_f_suc.push_back(possible_move[11]);
            current_f_suc.push_back(possible_move[12]);

            current_f_pred.push_back(possible_move[13]);
            current_f_pred.push_back(possible_move[14]);
            current_f_pred.push_back(possible_move[15]);

            root = newNode(move_data, current_f_suc, current_f_pred, current_vehicle_loads);
            root->parent = NULL;

            // cout <<"Logro crear root node" << endl;
            // cout << cont_ext << endl;
            next_possible_moves = CompleteRelocationFunction(initial_solution, instance, root);
            // possible_move[2], current_f_suc, current_f_pred, current_vehicle_loads);
            // cout << "Next possible move size: " << next_possible_moves.size() << endl;

            if (next_possible_moves.size() > 0) {
                for (int i = 0; i < next_possible_moves.size(); i++) {

                    possible_move.clear();
                    // possible_move = next_possible_moves[i];
                    possible_move = next_possible_moves.top();
                    next_possible_moves.pop();

                    move_data.clear();

                    move_data.push_back(possible_move[0]);  // Insert delta
                    move_data.push_back(possible_move[1]);  // Insert vehicle_one_id
                    move_data.push_back(possible_move[2]);  // Insert vehicle_two_id
                    move_data.push_back(possible_move[3]);  // Insert client to relocate
                    move_data.push_back(possible_move[4]);  // Insert relocation destination
                    move_data.push_back(possible_move[5]);  // Insert is_feasible
                    move_data.push_back(possible_move[12]); // Insert relocation destination id of client in position

                    if (possible_move[5] == 1) { // feasible solution
                        best_move_node = root;
                        first_best_move = {possible_move[1], possible_move[2], possible_move[3], possible_move[4],
                                           possible_move[12], possible_move[0]};
                        improved_solution = true;
                        // cout << "Improved at first move " << endl;
                        break;
                    }

                    current_vehicle_loads = root->vehicle_loads;
                    current_vehicle_loads[possible_move[1]] = possible_move[6]; // Update vehicle_one load
                    current_vehicle_loads[possible_move[2]] = possible_move[7]; // Update vehicle_two load

                    current_f_suc = root->f_suc;
                    current_f_pred = root->f_pred;

                    // Update forbidden clients moves

                    current_f_suc.push_back(possible_move[8]);
                    current_f_suc.push_back(possible_move[9]);
                    current_f_suc.push_back(possible_move[10]);
                    current_f_suc.push_back(possible_move[11]);
                    current_f_suc.push_back(possible_move[12]);

                    current_f_pred.push_back(possible_move[13]);
                    current_f_pred.push_back(possible_move[14]);
                    current_f_pred.push_back(possible_move[15]);

                    InsertNode(root, NULL, move_data, current_f_suc, current_f_pred, current_vehicle_loads);
                    //(root->child).push_back(newNode(move_data, current_f_suc, current_f_pred,
                    // current_vehicle_loads));
                }

                int cont = 0;

                // next_possible_moves.clear();
                next_possible_moves = my_priority_queue();
                bool went_back = false;
                bool no_more_childs = false;

                if ((root->child.size() > cont) && (improved_solution == false) && (no_more_childs == false)) {
                    while ((root->child.size() > cont) && (improved_solution == false) && (no_more_childs == false) &&
                           (explored_nodes < max_explored_nodes)) {

                        explored_nodes++;
                        if (explored_nodes >= max_explored_nodes) {
                            // counters[current_counter - 1] = 0;
                            current_counter = 0;
                            // cout << "Reached max nodes 3" << endl;
                        }

                        counters.clear();
                        counters.push_back(0);
                        Node *current_node = root->child[cont];
                        Node *parent_node = root;

                        while ((no_more_childs == false) && (improved_solution == false) &&
                               (explored_nodes < max_explored_nodes)) {
                            // current_node = parent_node->child[cont];
                            // next_possible_moves.clear();
                            next_possible_moves = my_priority_queue();
                            next_possible_moves = CompleteRelocationFunction(initial_solution, instance, current_node);
                            explored_nodes++;

                            if (next_possible_moves.size() > 0) {
                                counters.push_back(0);
                                current_counter += 1;
                                // cout << "now on level: " << current_counter << endl;
                                // cout << "next possible moves size: " << next_possible_moves.size() << endl;
                                for (int i = 0; i < next_possible_moves.size(); i++) {
                                    // cout << "Analizando " << i << endl;
                                    possible_move.clear();
                                    // possible_move = next_possible_moves[i];
                                    possible_move = next_possible_moves.top();
                                    next_possible_moves.pop();

                                    move_data.clear();

                                    move_data.push_back(possible_move[0]); // Insert delta
                                    move_data.push_back(possible_move[1]); // Insert vehicle_one_id
                                    move_data.push_back(possible_move[2]); // Insert vehicle_two_id
                                    move_data.push_back(possible_move[3]); // Insert client to relocate
                                    move_data.push_back(possible_move[4]); // Insert relocation destination
                                    move_data.push_back(possible_move[5]); // Insert is_feasible
                                    move_data.push_back(
                                            possible_move[12]); // Insert relocation destination id of client in position

                                    if (possible_move[5] == 1) { // feasible solution
                                        best_move_node = current_node;
                                        first_best_move = {possible_move[1], possible_move[2], possible_move[3],
                                                           possible_move[4], possible_move[12], possible_move[0]};
                                        improved_solution = true;
                                        // cout << "improved interior" << endl;

                                        break;
                                    }

                                    current_vehicle_loads = current_node->vehicle_loads;
                                    current_vehicle_loads[possible_move[1]] = possible_move[6]; // Update vehicle_one load
                                    current_vehicle_loads[possible_move[2]] = possible_move[7]; // Update vehicle_two load

                                    current_f_suc = current_node->f_suc;
                                    current_f_pred = current_node->f_pred;

                                    // Update forbidden clients moves

                                    current_f_suc.push_back(possible_move[8]);
                                    current_f_suc.push_back(possible_move[9]);
                                    current_f_suc.push_back(possible_move[10]);
                                    current_f_suc.push_back(possible_move[11]);
                                    current_f_suc.push_back(possible_move[12]);

                                    current_f_pred.push_back(possible_move[13]);
                                    current_f_pred.push_back(possible_move[14]);
                                    current_f_pred.push_back(possible_move[15]);

                                    InsertNode(current_node, NULL, move_data, current_f_suc, current_f_pred,
                                               current_vehicle_loads);
                                    // cout << "LLego aca " << i << endl;
                                }
                                if (improved_solution == false) {
                                    // cout << "Parent node child size " <<parent_node->child.size()<< endl;
                                    // cout << "Counters[current counter] " << counters[current_counter] << endl;
                                    parent_node = current_node;
                                    current_node = parent_node->child[counters[current_counter]];
                                    // cout << "llego aca" << endl;
                                }

                            } else {
                                went_back = false;
                                while ((current_counter >= 0) && (went_back == false)) {

                                    counters[current_counter] += 1;
                                    if (parent_node->child.size() > counters[current_counter]) {
                                        current_node = parent_node->child[counters[current_counter]];
                                        went_back = true;
                                    } else {
                                        if (parent_node->parent != NULL) {
                                            counters[current_counter] = 0;
                                            current_counter -= 1;
                                            // cout << "parent node va a devolverse" << endl;
                                            // if(parent_node->parent == NULL){
                                            //     cout <<"dpto error" <<endl;
                                            //}
                                            parent_node = parent_node->parent;
                                            // cout << "se devolvio correctamente" <<endl;
                                            // current_node = parent_node->child[counters[current_counter]];
                                        } else {
                                            counters[current_counter] = 0;
                                            no_more_childs = true;
                                            // cout << "No more childs" << endl;
                                            current_counter = -1;
                                        }
                                    }
                                }
                                if (current_counter == -1) {
                                    // cout << "entro aqui" << endl;
                                    current_counter = 0;
                                    // current_node = parent_node->child[counters[current_counter]];
                                }
                            }
                        }
                        cont = cont + 1;
                    }
                }
            }
            cont_ext += 1;
            if (improved_solution == false) {
                deleteNode(root);
                delete root;
            }
        }
        if (improved_solution == true) {

            Node *actual_node = best_move_node;
            best_move.push_back(first_best_move);
            best_move.push_back(best_move_node->relocate_move);
            //cout << "Expected delta: " << first_best_move[5] << endl;

            bool stop_now = false;
            while (stop_now == false) {
                if (actual_node->parent == NULL) {
                    // cout << "llegue a null" << endl;
                    // best_move.push_back(actual_node->relocate_move);
                    stop_now = true;
                    break;
                } else {
                    actual_node = actual_node->parent;
                    // cout << "Vehicles: " << actual_node->relocate_move[0] << " " << actual_node->relocate_move[1] <<
                    // endl;
                    best_move.push_back(actual_node->relocate_move);
                    // cout << "hice push_back" << endl;
                }
            }
            // cout << "Vehicles: " << best_move_node->relocate_move[0] << " " << best_move_node->relocate_move[1] << endl;
            // cout << "Vehicles: " << first_best_move[0] << " " << first_best_move[1] << endl;

            // best_move.push_back(first_best_move);

            // cout << "improved solution" << endl;
            deleteNode(root);
            // delete root;
            //cout << "Initial cost: " <<initial_solution->cost;
            initial_solution = apply_best_move(initial_solution, best_move);
            initial_solution->cost = initial_solution->compute_cost(instance);
            //cout << "Final cost: " <<initial_solution->cost <<endl;
            // initial_solution->is_feasible();
            //return initial_solution;
            delete root;
        }
    }
    return initial_solution;
}

/**
 * Apply the relocate (insertion) to the best possible move in the neighborhood
 * @param initial_solution Current solution of the VRP
 * @param best_move {vehicle_one_id, vehicle_two_id, customer_one_id, customer_new_pos}
 *
 * @return local_optima: VRP solution corresponding to the local optima of this neighborhood
 *          if there's no local optima, returns initial_solution instead
 */
std::shared_ptr <Solution> RelocationChainNeighborhood_LimitNodes::apply_best_move(std::shared_ptr <Solution> initial_solution,
                                                                  vector<vector<int>>& complete_best_move)
{
    if(complete_best_move.size() != 0) {
        for(int i = complete_best_move.size() - 1; i >= 0; i--) {
            vector<int> best_move = complete_best_move[i];
            // std::shared_ptr <Solution> local_optima = initial_solution;
            // cout << "Vehicles: " << best_move[0] << " " << best_move[1] << endl;
            int vehicle_one_id = best_move[0];
            int vehicle_two_id = best_move[1];
            int customer_to_relocate = best_move[2];
            int customer_new_pos = best_move[3];
            int customer_new_pos_id = best_move[4];
            int customer_one_id;

            for(int j = 1; j < initial_solution->vehicles[vehicle_one_id]->customers.size() - 1; j++) {
                if(customer_to_relocate == initial_solution->vehicles[vehicle_one_id]->customers[j]->id) {
                    customer_one_id = j;
                    break;
                }
            }
            for(int k = 1; k < initial_solution->vehicles[vehicle_two_id]->customers.size() - 1; k++) {
                if(customer_new_pos_id == initial_solution->vehicles[vehicle_two_id]->customers[k]->id) {
                    customer_new_pos = k;
                    break;
                }
            }

            shared_ptr<Customer> customer_one = initial_solution->vehicles[vehicle_one_id]->customers[customer_one_id];
            // shared_ptr<Customer> customer_one = initial_solution->customers[customer_to_relocate];

            initial_solution->vehicles[vehicle_one_id]->remove_customer(customer_one_id);
            initial_solution->vehicles[vehicle_two_id]->insert_customer(customer_one, customer_new_pos);

            initial_solution->vehicles[vehicle_two_id]->customers[customer_new_pos]->vehicle_route = vehicle_two_id;
			initial_solution->customers[initial_solution->vehicles[vehicle_two_id]->customers[customer_new_pos]->id]->vehicle_route = vehicle_two_id;

            /*
            cout << "-----------loads-------------" << endl;
            for(int i = 0; i < initial_solution->vehicles.size(); i++) {
                cout << initial_solution->vehicles[i]->load << " ";
            }
            cout << " " << endl;
             * */
            // cout << "vehicle_one " << vehicle_one_id << endl;
            // cout << "vehicle_two " << vehicle_two_id << endl;
        }
        // initial_solution->is_feasible();
        return initial_solution;
    } else {
        return initial_solution;
    }
}

RelocationChainNeighborhood_LimitDepth::RelocationChainNeighborhood_LimitDepth()
{
}

RelocationChainNeighborhood_LimitDepth::RelocationChainNeighborhood_LimitDepth(bool first_improvement /*= false*/)
        : Neighborhood(first_improvement, "B")
{
}

RelocationChainNeighborhood_LimitDepth::~RelocationChainNeighborhood_LimitDepth()
{
}

std::shared_ptr <Solution> RelocationChainNeighborhood_LimitDepth::search_neighbors(std::shared_ptr <Solution> initial_solution, const Instance& instance)
{
    for(int k = 0; k < initial_solution->vehicles.size(); k++) {
    std::shared_ptr <Solution> local_optima = initial_solution;
    int n_vehicles = initial_solution->vehicles.size();

    vector<int> select_best_move;
    vector<int> move_data;
    vector<vector<int>> best_move;
    vector<int> possible_move;
    vector<int> first_best_move;
    vector<int> counters;
    // vector<vector<int>> first_possible_moves;
    my_priority_queue first_possible_moves;
    // vector<vector<int>> next_possible_moves;
    my_priority_queue next_possible_moves;
    vector<vector<int>> best_possible_moves;
    vector<unordered_map<int, int>> vehicle_loads;
    vector<vector<int>> f_suc;
    vector<vector<int>> f_pred;
    vector<vector<int>> delta_variation;
    vector<int> current_f_suc;
    vector<int> current_f_pred;
    unordered_map<int, int> initial_vehicle_loads;
    unordered_map<int, int> current_vehicle_loads;
    Node* best_move_node;
    Node* root;


    int level = 0;
    int final_move_feasible = 0;
    bool improved_solution = false;
    int cont_ext = 0;
    int current_counter = 0;
    int explored_nodes = 0;
    int max_first_nodes = 0;
    int max_explored_nodes = 100;

    for(int i = 0; i < initial_solution->vehicles.size(); i++) {
        initial_vehicle_loads[i] = initial_solution->vehicles[i]->load;
    }

        first_possible_moves = FirstRelocationFunction(initial_solution, instance, k);
        max_first_nodes = first_possible_moves.size();
        // cout << "First possible moves size: " << first_possible_moves.size() << endl;

        while ((improved_solution == false) && (cont_ext < max_first_nodes) &&
               (explored_nodes < max_explored_nodes + 3)) {
            explored_nodes++;
            // cout << "First possible move: " << cont_ext << endl;

            /*
            if (explored_nodes >= max_explored_nodes){
                explored_nodes = 0;
                current_counter = 0;
            }
            */

            possible_move.clear();
            // possible_move = first_possible_moves[cont_ext];
            possible_move = first_possible_moves.top();
            first_possible_moves.pop();

            move_data.clear();

            move_data.push_back(possible_move[0]);  // Insert delta
            move_data.push_back(possible_move[1]);  // Insert vehicle_one_id
            move_data.push_back(possible_move[2]);  // Insert vehicle_two_id
            move_data.push_back(possible_move[3]);  // Insert client to relocate
            move_data.push_back(possible_move[4]);  // Insert relocation destination
            move_data.push_back(possible_move[5]);  // Insert is_feasible
            move_data.push_back(possible_move[12]); // Insert relocation destination id of client in position

            current_vehicle_loads = initial_vehicle_loads;
            current_vehicle_loads[possible_move[1]] = possible_move[6];
            current_vehicle_loads[possible_move[2]] = possible_move[7];

            current_f_suc.clear();
            current_f_pred.clear();

            current_f_suc.push_back(possible_move[8]);
            current_f_suc.push_back(possible_move[9]);
            current_f_suc.push_back(possible_move[10]);
            current_f_suc.push_back(possible_move[11]);
            current_f_suc.push_back(possible_move[12]);

            current_f_pred.push_back(possible_move[13]);
            current_f_pred.push_back(possible_move[14]);
            current_f_pred.push_back(possible_move[15]);

            root = newNode(move_data, current_f_suc, current_f_pred, current_vehicle_loads);
            root->parent = NULL;

            // cout <<"Logro crear root node" << endl;
            // cout << cont_ext << endl;
            next_possible_moves = CompleteRelocationFunction(initial_solution, instance, root);
            // possible_move[2], current_f_suc, current_f_pred, current_vehicle_loads);
            // cout << "Next possible move size: " << next_possible_moves.size() << endl;

            if (next_possible_moves.size() > 0) {
                for (int i = 0; i < next_possible_moves.size(); i++) {

                    possible_move.clear();
                    // possible_move = next_possible_moves[i];
                    possible_move = next_possible_moves.top();
                    next_possible_moves.pop();

                    move_data.clear();

                    move_data.push_back(possible_move[0]);  // Insert delta
                    move_data.push_back(possible_move[1]);  // Insert vehicle_one_id
                    move_data.push_back(possible_move[2]);  // Insert vehicle_two_id
                    move_data.push_back(possible_move[3]);  // Insert client to relocate
                    move_data.push_back(possible_move[4]);  // Insert relocation destination
                    move_data.push_back(possible_move[5]);  // Insert is_feasible
                    move_data.push_back(possible_move[12]); // Insert relocation destination id of client in position

                    if (possible_move[5] == 1) { // feasible solution
                        best_move_node = root;
                        first_best_move = {possible_move[1], possible_move[2], possible_move[3], possible_move[4],
                                           possible_move[12], possible_move[0]};
                        improved_solution = true;
                        // cout << "Improved at first move " << endl;
                        break;
                    }

                    current_vehicle_loads = root->vehicle_loads;
                    current_vehicle_loads[possible_move[1]] = possible_move[6]; // Update vehicle_one load
                    current_vehicle_loads[possible_move[2]] = possible_move[7]; // Update vehicle_two load

                    current_f_suc = root->f_suc;
                    current_f_pred = root->f_pred;

                    // Update forbidden clients moves

                    current_f_suc.push_back(possible_move[8]);
                    current_f_suc.push_back(possible_move[9]);
                    current_f_suc.push_back(possible_move[10]);
                    current_f_suc.push_back(possible_move[11]);
                    current_f_suc.push_back(possible_move[12]);

                    current_f_pred.push_back(possible_move[13]);
                    current_f_pred.push_back(possible_move[14]);
                    current_f_pred.push_back(possible_move[15]);

                    InsertNode(root, NULL, move_data, current_f_suc, current_f_pred, current_vehicle_loads);
                    //(root->child).push_back(newNode(move_data, current_f_suc, current_f_pred,
                    // current_vehicle_loads));
                }

                int cont = 0;

                // next_possible_moves.clear();
                next_possible_moves = my_priority_queue();
                bool went_back = false;
                bool no_more_childs = false;

                if ((root->child.size() > cont) && (improved_solution == false) && (no_more_childs == false)) {
                    while ((root->child.size() > cont) && (improved_solution == false) && (no_more_childs == false) &&
                           (explored_nodes < max_explored_nodes)) {

                        explored_nodes++;
                        if (explored_nodes >= max_explored_nodes) {
                            current_counter = 0;
                        }
                        if (counters.size() >= 3) {
                            current_counter = 0;
                        }

                        counters.clear();
                        counters.push_back(0);
                        Node *current_node = root->child[cont];
                        Node *parent_node = root;

                        while ((no_more_childs == false) && (improved_solution == false) &&
                               (explored_nodes < max_explored_nodes) && (counters.size() < 3)) {
                            // current_node = parent_node->child[cont];
                            // next_possible_moves.clear();
                            next_possible_moves = my_priority_queue();
                            next_possible_moves = CompleteRelocationFunction(initial_solution, instance, current_node);
                            explored_nodes++;

                            if (next_possible_moves.size() > 0) {
                                counters.push_back(0);
                                current_counter += 1;
                                // cout << "now on level: " << current_counter << endl;
                                // cout << "next possible moves size: " << next_possible_moves.size() << endl;
                                for (int i = 0; i < next_possible_moves.size(); i++) {
                                    // cout << "Analizando " << i << endl;
                                    possible_move.clear();
                                    // possible_move = next_possible_moves[i];
                                    possible_move = next_possible_moves.top();
                                    next_possible_moves.pop();

                                    move_data.clear();

                                    move_data.push_back(possible_move[0]); // Insert delta
                                    move_data.push_back(possible_move[1]); // Insert vehicle_one_id
                                    move_data.push_back(possible_move[2]); // Insert vehicle_two_id
                                    move_data.push_back(possible_move[3]); // Insert client to relocate
                                    move_data.push_back(possible_move[4]); // Insert relocation destination
                                    move_data.push_back(possible_move[5]); // Insert is_feasible
                                    move_data.push_back(
                                            possible_move[12]); // Insert relocation destination id of client in position

                                    if (possible_move[5] == 1) { // feasible solution
                                        best_move_node = current_node;
                                        first_best_move = {possible_move[1], possible_move[2], possible_move[3],
                                                           possible_move[4], possible_move[12], possible_move[0]};
                                        improved_solution = true;
                                        // cout << "improved interior" << endl;

                                        break;
                                    }

                                    current_vehicle_loads = current_node->vehicle_loads;
                                    current_vehicle_loads[possible_move[1]] = possible_move[6]; // Update vehicle_one load
                                    current_vehicle_loads[possible_move[2]] = possible_move[7]; // Update vehicle_two load

                                    current_f_suc = current_node->f_suc;
                                    current_f_pred = current_node->f_pred;

                                    // Update forbidden clients moves

                                    current_f_suc.push_back(possible_move[8]);
                                    current_f_suc.push_back(possible_move[9]);
                                    current_f_suc.push_back(possible_move[10]);
                                    current_f_suc.push_back(possible_move[11]);
                                    current_f_suc.push_back(possible_move[12]);

                                    current_f_pred.push_back(possible_move[13]);
                                    current_f_pred.push_back(possible_move[14]);
                                    current_f_pred.push_back(possible_move[15]);

                                    InsertNode(current_node, NULL, move_data, current_f_suc, current_f_pred,
                                               current_vehicle_loads);
                                    // cout << "LLego aca " << i << endl;
                                }
                                if (improved_solution == false) {
                                    // cout << "Parent node child size " <<parent_node->child.size()<< endl;
                                    // cout << "Counters[current counter] " << counters[current_counter] << endl;
                                    parent_node = current_node;
                                    current_node = parent_node->child[counters[current_counter]];
                                    // cout << "llego aca" << endl;
                                }

                            } else {
                                went_back = false;
                                while ((current_counter >= 0) && (went_back == false)) {

                                    counters[current_counter] += 1;
                                    if (parent_node->child.size() > counters[current_counter]) {
                                        current_node = parent_node->child[counters[current_counter]];
                                        went_back = true;
                                    } else {
                                        if (parent_node->parent != NULL) {
                                            counters[current_counter] = 0;
                                            current_counter -= 1;
                                            // cout << "parent node va a devolverse" << endl;
                                            // if(parent_node->parent == NULL){
                                            //     cout <<"dpto error" <<endl;
                                            //}
                                            parent_node = parent_node->parent;
                                            // cout << "se devolvio correctamente" <<endl;
                                            // current_node = parent_node->child[counters[current_counter]];
                                        } else {
                                            counters[current_counter] = 0;
                                            no_more_childs = true;
                                            // cout << "No more childs" << endl;
                                            current_counter = -1;
                                        }
                                    }
                                }
                                if (current_counter == -1) {
                                    // cout << "entro aqui" << endl;
                                    current_counter = 0;
                                    // current_node = parent_node->child[counters[current_counter]];
                                }
                            }
                        }
                        cont = cont + 1;
                    }
                }
            }
            cont_ext += 1;
            if (improved_solution == false) {
                deleteNode(root);
                delete root;
            }
        }
        if (improved_solution == true) {

            Node *actual_node = best_move_node;
            best_move.push_back(first_best_move);
            best_move.push_back(best_move_node->relocate_move);
            // cout << "Expected delta: " << first_best_move[5] << endl;

            bool stop_now = false;
            while (stop_now == false) {
                if (actual_node->parent == NULL) {
                    // cout << "llegue a null" << endl;
                    // best_move.push_back(actual_node->relocate_move);
                    stop_now = true;
                    break;
                } else {
                    actual_node = actual_node->parent;
                    // cout << "Vehicles: " << actual_node->relocate_move[0] << " " << actual_node->relocate_move[1] <<
                    // endl;
                    best_move.push_back(actual_node->relocate_move);
                    // cout << "hice push_back" << endl;
                }
            }
            // cout << "Vehicles: " << best_move_node->relocate_move[0] << " " << best_move_node->relocate_move[1] << endl;
            // cout << "Vehicles: " << first_best_move[0] << " " << first_best_move[1] << endl;

            // best_move.push_back(first_best_move);

            // cout << "improved solution" << endl;
            deleteNode(root);
            // delete root;
            // cout << "Initial cost: " <<initial_solution->cost;
            initial_solution = apply_best_move(initial_solution, best_move);
            initial_solution->cost = initial_solution->compute_cost(instance);
            // cout << "Final cost: " <<initial_solution->cost <<endl;
            // initial_solution->is_feasible();
            //return initial_solution;
            delete root;
        }
    }
    return initial_solution;

}

/**
 * Apply the relocate (insertion) to the best possible move in the neighborhood
 * @param initial_solution Current solution of the VRP
 * @param best_move {vehicle_one_id, vehicle_two_id, customer_one_id, customer_new_pos}
 *
 * @return local_optima: VRP solution corresponding to the local optima of this neighborhood
 *          if there's no local optima, returns initial_solution instead
 */
std::shared_ptr <Solution> RelocationChainNeighborhood_LimitDepth::apply_best_move(std::shared_ptr <Solution> initial_solution,
                                                                  vector<vector<int>>& complete_best_move)
{
    if(complete_best_move.size() != 0) {
        for(int i = complete_best_move.size() - 1; i >= 0; i--) {
            vector<int> best_move = complete_best_move[i];
            // std::shared_ptr <Solution> local_optima = initial_solution;
            // cout << "Vehicles: " << best_move[0] << " " << best_move[1] << endl;
            int vehicle_one_id = best_move[0];
            int vehicle_two_id = best_move[1];
            int customer_to_relocate = best_move[2];
            int customer_new_pos = best_move[3];
            int customer_new_pos_id = best_move[4];
            int customer_one_id;

            for(int i = 1; i < initial_solution->vehicles[vehicle_one_id]->customers.size() - 1; i++) {
                if(customer_to_relocate == initial_solution->vehicles[vehicle_one_id]->customers[i]->id) {
                    customer_one_id = i;
                    break;
                }
            }
            for(int i = 1; i < initial_solution->vehicles[vehicle_two_id]->customers.size() - 1; i++) {
                if(customer_new_pos_id == initial_solution->vehicles[vehicle_two_id]->customers[i]->id) {
                    customer_new_pos = i;
                    break;
                }
            }

            shared_ptr<Customer> customer_one = initial_solution->vehicles[vehicle_one_id]->customers[customer_one_id];
            // shared_ptr<Customer> customer_one = initial_solution->customers[customer_to_relocate];

            initial_solution->vehicles[vehicle_one_id]->remove_customer(customer_one_id);
            initial_solution->vehicles[vehicle_two_id]->insert_customer(customer_one, customer_new_pos);

            initial_solution->vehicles[vehicle_two_id]->customers[customer_new_pos]->vehicle_route = vehicle_two_id;
			initial_solution->customers[initial_solution->vehicles[vehicle_two_id]->customers[customer_new_pos]->id]->vehicle_route = vehicle_two_id;

            /*
            cout << "-----------loads-------------" << endl;
            for(int i = 0; i < initial_solution->vehicles.size(); i++) {
                cout << initial_solution->vehicles[i]->load << " ";
            }
            cout << " " << endl;
             * */
            // cout << "vehicle_one " << vehicle_one_id << endl;
            // cout << "vehicle_two " << vehicle_two_id << endl;
        }
        // initial_solution->is_feasible();
        return initial_solution;
    } else {
        return initial_solution;
    }
}


SwapStarNeighborhood::SwapStarNeighborhood()
{
}

SwapStarNeighborhood::SwapStarNeighborhood(bool first_improvement /*= false*/)
        : Neighborhood(first_improvement, "A")
{
}

SwapStarNeighborhood::~SwapStarNeighborhood()
{
}

vector<vector<int>>
PreprocessInsertions(std::shared_ptr <Solution> initial_solution, const Instance& instance, int vehicle_one_id, int vehicle_two_id)
{
    int n_vehicles = initial_solution->vehicles.size();

    int min_delta = 0;
    vector<int> best_move;
    vector<vector<int>> possible_moves;
    vector<vector<int>> best_possible_moves;
    vector<int> temp_move;
    Vehicle vehicle_one;
    Vehicle vehicle_two;
    // customers
    Customer a;
    Customer b;
    Customer c;

    Customer f;
    Customer g;

    // demands
    int route1_new_demand;
    int route2_new_demand;
    // dist_matrix
    // delta
    float delta;
    int temp_delta;

    best_possible_moves.clear();
    possible_moves.clear();
    best_move.clear();
    temp_move.clear();

    vehicle_one = *initial_solution->vehicles[vehicle_one_id];
    vehicle_two = *initial_solution->vehicles[vehicle_two_id];

    for(int customer_one_id = 1; customer_one_id < vehicle_one.customers.size() - 1; customer_one_id++) {
        a = *vehicle_one.customers[customer_one_id - 1];
        b = *vehicle_one.customers[customer_one_id];
        c = *vehicle_one.customers[customer_one_id + 1];

        // temp_delta = instance.insertion_matrix[f.id][g.id][b.id] +
        // instance.remotion_matrix[h.id][j.id][i.id];
        for(int customer_new_pos = 1; customer_new_pos < vehicle_two.customers.size() - 1; customer_new_pos++) {
            /*
                route1 =  [A B C D]
                route2 =  [E F G H]
                Relocate(B, pos(G)) leaving us with:
                new_route1 = [A C D]
                new_route2 = [E F B G H]
            */

            f = *vehicle_two.customers[customer_new_pos - 1];
            g = *vehicle_two.customers[customer_new_pos];

            /*
             delta = instance.insertion_matrix[f.id][g.id][b.id] + instance.remotion_matrix[a.id][c.id][b.id];
            */


            delta = instance.dist_matrix[a.id][c.id] + instance.dist_matrix[f.id][b.id] + instance.dist_matrix[b.id][g.id] \
                - instance.dist_matrix[a.id][b.id] - instance.dist_matrix[b.id][c.id] - instance.dist_matrix[f.id][g.id] ;

            // Store the current best possible move
            if(delta < 0) {
                // min_delta = delta;
                best_move.clear();
                best_move = {
                        static_cast<int> (delta),
                        vehicle_one_id,
                        vehicle_two_id,
                        b.id,
                        g.id,
                        f.id,
                        customer_one_id,
                        customer_new_pos,
                };
                possible_moves.push_back(best_move);
            }
        }
        std::sort(possible_moves.begin(), possible_moves.end(),
                  [](const std::vector<int>& a, const std::vector<int>& b) { return a[0] < b[0]; });

        // cout << "Listo Sorting possible moves" << endl;

        if(possible_moves.size() == 1) {
            best_possible_moves.push_back({possible_moves[0],});
        } else if(possible_moves.size() == 2) {
            for(int j = 0; j < 2; j++) {
                for(int i = 0; i < 8; i++) {
                    temp_move.push_back(possible_moves[j][i]);
                }
            }
            best_possible_moves.push_back(temp_move);

        } else if(possible_moves.size() >= 3) {
            for(int j = 0; j < 3; j++) {
                for(int i = 0; i < 8; i++) {
                    temp_move.push_back(possible_moves[j][i]);
                }
            }
            best_possible_moves.push_back(temp_move);
        }
        possible_moves.clear();
        temp_move.clear();
    }

    return best_possible_moves;
}

vector<int> EvalSwapStar(std::shared_ptr <Solution> initial_solution,
                         const Instance& instance,
                         vector<int> possible_moves,
                         int destination_id,
                         int destination)
{
    int n_vehicles = initial_solution->vehicles.size();

    int min_delta = 0;
    vector<int> best_move;
    int best_delta = 0;
    Vehicle vehicle_one;
    Vehicle vehicle_two;

    // customers
    Customer a;
    Customer b;
    Customer c;

    Customer f;
    Customer g;
    Customer h;
    // demands
    int route1_new_demand;
    int route2_new_demand;
    // dist_matrix
    // delta
    float delta;

    /*
                     Possible moves structure = {  delta, vehicle_one_id, vehicle_two_id, b.id, g.id, f.id,
       customer_one_id, customer_new_pos,
                      *  delta, vehicle_one_id, vehicle_two_id, b.id, g.id, f.id, customer_one_id,
       customer_new_pos,
                      *  delta, vehicle_one_id, vehicle_two_id, b.id, g.id, f.id, customer_one_id,
       customer_new_pos, }
                     */
    // cout << destination_id << " " << possible_moves[4] << " " << possible_moves[5] << endl;
    best_move.clear();

    if((possible_moves[4] != destination_id) && (possible_moves[5] != destination_id)) {
        for(int i = 0; i < 8; i++) {
            best_move.push_back(possible_moves[i]);
        }
        best_delta = best_move[0];
    } /*else {
        if(possible_moves.size() > 8) {
            // cout << destination_id << " " << possible_moves[12] << " " << possible_moves[13] << endl;
            if((possible_moves[12] != destination_id) && (possible_moves[13] != destination_id)) {
                for(int i = 8; i < 16; i++) {
                    best_move.push_back(possible_moves[i]);
                }
                best_delta = best_move[0];
            } else {
                if(possible_moves.size() > 16) {
                    // cout << destination_id << " " << possible_moves[20] << " " << possible_moves[21] << endl;
                    if((possible_moves[20] != destination_id) && (possible_moves[21] != destination_id)) {
                        for(int i = 16; i < 24; i++) {
                            best_move.push_back(possible_moves[i]);
                        }

                        best_delta = best_move[0];
                    }
                }
            }
        }
    }*/

    vehicle_one = *initial_solution->vehicles[possible_moves[1]];
    vehicle_two = *initial_solution->vehicles[possible_moves[2]];

    int customer_one_id = possible_moves[6];
    int customer_two_id = destination;
    /*
        route1 =  [A B C D]
        route2 =  [E F G H]
        Swap(B, G) leaving us with:
        new_route1 = [G A C D]
        new_route2 = [E F B H]
    */

    a = *vehicle_one.customers[customer_one_id - 1];
    b = *vehicle_one.customers[customer_one_id];
    c = *vehicle_one.customers[customer_one_id + 1];

    f = *vehicle_two.customers[customer_two_id - 1];
    g = *vehicle_two.customers[customer_two_id];
    h = *vehicle_two.customers[customer_two_id + 1];

    // Check constraint of load violation
    route1_new_demand = vehicle_one.load - b.demand + g.demand;
    route2_new_demand = vehicle_two.load - g.demand + b.demand;
	
    /*
        Evaluate the improvement delta of both route1 and route2, using distance matrix
        Delta = D(a,g) + D(g,c) + D(f,b) + D(b,h) - D(a,b) -D(b,c) - D(f,g) - D(g,h)
    delta = dist[a.id][g.id] + dist[g.id][c.id] + dist[f.id][b.id] \
                        + dist[b.id][h.id] - dist[a.id][b.id] - dist[b.id][c.id] \
                        - dist[f.id][g.id] - dist[g.id][h.id];
    delta = instance.insertion_matrix[f.id][h.id][b.id] + instance.insertion_matrix[a.id][c.id][g.id] +
            instance.remotion_matrix[f.id][h.id][g.id] + instance.remotion_matrix[a.id][c.id][b.id];
    */


    delta = +instance.dist_matrix[a.id][g.id] + instance.dist_matrix[g.id][c.id] + instance.dist_matrix[f.id][b.id] + instance.dist_matrix[b.id][h.id]
            -instance.dist_matrix[a.id][b.id] - instance.dist_matrix[b.id][c.id] - instance.dist_matrix[f.id][g.id] - instance.dist_matrix[g.id][h.id];

    // cout << "delta" << delta << endl;
    // Store the current best possible move
	if((route1_new_demand <= vehicle_one.Q) && (route2_new_demand <= vehicle_two.Q)) {
    if(delta < best_delta) {
        // cout << "Swap move " << b.id <<" "<< g.id<< " delta "<< delta/2 <<endl;
        best_delta = delta;
        // best_move = { vehicle_one_id, vehicle_two_id, customer_one_id, customer_two_id };
        best_move.clear();
        best_move = {
                best_delta / 2,
                possible_moves[1],
                possible_moves[2],
                b.id,
                g.id,
                f.id,
                customer_one_id,
                customer_two_id,
				route1_new_demand,
				route2_new_demand
        };
    }
	}
	else
	{
		best_move.clear();
	}

    return best_move;
}

std::shared_ptr <Solution> SwapStarNeighborhood::search_neighbors(std::shared_ptr <Solution> initial_solution, const Instance& instance)
{
    int n_vehicles = initial_solution->vehicles.size();

    float min_delta = 0;
    vector<vector<int>> final_best_moves;
    vector<vector<int>> possible_moves_1;
    vector<vector<int>> possible_moves_2;
    vector<int> best_possible_move_1;
    vector<int> best_possible_move_2;
    Vehicle vehicle_one;
    Vehicle vehicle_two;
    int initial_value_for;
    // customers
    Customer a;
    Customer b;
    Customer c;

    Customer f;
    Customer g;
    Customer h;
    // demands
    int route1_new_demand;
    int route2_new_demand;
    // dist_matrix
    // delta
    float delta;

    // for(int vehicle_one_id = 0; vehicle_one_id < n_vehicles - 1; vehicle_one_id++) {
    for(int vehicle_one_id = 0; vehicle_one_id < n_vehicles - 1; vehicle_one_id++) {
        // cout << "Start vehicle 2.." << endl;
        vehicle_one = *initial_solution->vehicles[vehicle_one_id];
        // Cambiar este if por rutina de coordenadas polares y mismo espacio polar
        // for(int vehicle_two_id = vehicle_one_id - 2; vehicle_two_id < vehicle_one_id + 2; vehicle_two_id++) {
        for(int vehicle_two_id = vehicle_one_id + 1; vehicle_two_id < n_vehicles; vehicle_two_id++) {
            // cout << "Start possible moves" << endl;
            /*
            if(vehicle_two_id < 0) {
                vehicle_two_id = n_vehicles - 1;
            }
            if(vehicle_two_id >= n_vehicles) {
                vehicle_two_id = 0;
            }
            */
            if(vehicle_one_id == vehicle_two_id) {
                continue;
            }
            vehicle_two = *initial_solution->vehicles[vehicle_two_id];

            if(CircleSector::overlap(vehicle_one.sector, vehicle_two.sector)) {
                // if (true){
                //    cout << "Overlapping sectors " <<  vehicle_one_id << " " << vehicle_two_id << endl;
                possible_moves_1.clear();
                possible_moves_2.clear();
                possible_moves_1 = PreprocessInsertions(initial_solution, instance, vehicle_one_id, vehicle_two_id);
                possible_moves_2 = PreprocessInsertions(initial_solution, instance, vehicle_two_id, vehicle_one_id);

                // cout << "Done finding possible moves" << endl;

                /*
                 Possible moves structure = {  delta, vehicle_one_id, vehicle_two_id, b.id, g.id, f.id,
                 customer_one_id, customer_new_pos,
                  *  delta, vehicle_one_id, vehicle_two_id, b.id, g.id, f.id, customer_one_id, customer_new_pos,
                  *  delta, vehicle_one_id, vehicle_two_id, b.id, g.id, f.id, customer_one_id,
                 customer_new_pos, }
                 */

                int delta_best = 0;
                final_best_moves.clear();
                final_best_moves.push_back({ 0 });
                final_best_moves.push_back({ 0 });

                if((possible_moves_1.size() > 0) && (possible_moves_2.size() > 0)) {
                    for(int i = 0; i < possible_moves_1.size(); i++) {
                        // cout << "Arrived here 1" << endl;
                        b = *vehicle_one.customers[possible_moves_1[i][6]];
                        // cout << "Arrived here 2" << endl;
                        for(int j = 0; j < possible_moves_2.size(); j++) {
                            // cout << "Arrived here 3" << endl;
                            g = *vehicle_two.customers[possible_moves_2[j][6]];
                            // cout << "Arrived here 4" << endl;

                            route1_new_demand = vehicle_one.load - b.demand + g.demand;
                            route2_new_demand = vehicle_two.load - g.demand + b.demand;

                            if((route1_new_demand <= vehicle_one.Q) && (route2_new_demand <= vehicle_two.Q)) {

                                // cout << "Start Eval swap" << endl;
                                best_possible_move_1.clear();
                                best_possible_move_2.clear();
                                best_possible_move_1 = EvalSwapStar(initial_solution, instance, possible_moves_1[i],
                                                                    possible_moves_2[j][3], possible_moves_2[j][6]);
                                best_possible_move_2 = EvalSwapStar(initial_solution, instance, possible_moves_2[j],
                                                                    possible_moves_1[i][3], possible_moves_1[i][6]);

                                // cout << "Done with EvalSwap. " << best_possible_move_1.size() << " " <<
                                // best_possible_move_2.size() <<endl;

                                if((best_possible_move_1.size() > 0) && (best_possible_move_2.size() > 0)) {
                                    if(best_possible_move_1[0] + best_possible_move_2[0] <
                                       final_best_moves[0][0] + final_best_moves[1][0]) {
                                        // cout << "Found better move" << endl;
                                        final_best_moves.clear();
                                        final_best_moves.push_back(best_possible_move_1);
                                        final_best_moves.push_back(best_possible_move_2);
                                        // cout << "Done" << endl;
                                    }
                                }
                            }
                        }
                    }
                }
                if(final_best_moves[0][0] + final_best_moves[1][0] < 0) {
                    // cout << "Applying best possible moves" << endl;
                    /*
                    cout << "min_delta " << final_best_moves[0][0] + final_best_moves[1][0] << " C1 " <<
                    final_best_moves[0][3] << " C2 " << final_best_moves[0][4] << endl; float veh_1_cost =
                    initial_solution->vehicles[final_best_moves[0][1]]->compute_cost(instance); float
                    veh_2_cost = initial_solution->vehicles[final_best_moves[0][2]]->compute_cost(instance);
                    cout << "cost_before
                    "
                    << veh_1_cost << " " << veh_2_cost << endl;
                    */
                    initial_solution = apply_best_move(initial_solution, final_best_moves);
                    // initial_solution->cost = initial_solution->cost + final_best_moves[0][0] +
                    // final_best_moves[1][0];
                    initial_solution->cost = initial_solution->compute_cost(instance);
                    initial_solution->vehicles[final_best_moves[0][1]]->compute_cost(instance);
                    initial_solution->vehicles[final_best_moves[1][2]]->compute_cost(instance);

                    // cout << initial_solution->cost << endl;
                    /*
                    veh_1_cost = initial_solution->vehicles[final_best_moves[0][1]]->compute_cost(instance);
                    veh_2_cost = initial_solution->vehicles[final_best_moves[0][2]]->compute_cost(instance);
                    cout << "cost_after " << veh_1_cost << " " << veh_2_cost << endl;
                    */
                    // cout << "Done applying best possible moves" << endl;
                }
            }
        }
    }
    return initial_solution;
}

/**
 * Apply the Swap (exchange) to the best possible move in the neighborhood
 * @param initial_solution Current solution of the VRP
 * @param best_move {vehicle_one_id, vehicle_two_id, customer_one_id, customer_two_id}
 *
 * @return local_optima: VRP solution corresponding to the local optima of this neighborhood
 *          if there's no local optima, returns initial_solution instead
 */
std::shared_ptr <Solution> SwapStarNeighborhood::apply_best_move(std::shared_ptr <Solution> initial_solution, vector<vector<int>>& complete_best_move)
{
    /*
			best_move_1 structure= {
				best_delta / 2,
                vehicle_one_id,
                vehicle_two_id,
                b.id,
                g.id,
                f.id,
                customer_one_id,
                customer_two_id,
				route1_new_demand,
				route2_new_demand
			}
             */

    if(complete_best_move.size() != 0) {
        vector<int> best_move_1 = complete_best_move[0];
        vector<int> best_move_2 = complete_best_move[1];
        // std::shared_ptr <Solution> local_optima = initial_solution;
        int vehicle_one_id = best_move_1[1];
        int vehicle_two_id = best_move_1[2];
        int customer_to_relocate_1 = best_move_1[6];
        int customer_to_relocate_2 = best_move_2[3];
        int customer_new_pos_1 = best_move_1[7];
        int customer_new_pos_2 = best_move_2[7];
        int customer_two_id = -1;
        int customer_one_id = -1;
		int route1_new_demand = best_move_1[8];
		int route2_new_demand = best_move_1[9];

        // cout << initial_solution->vehicles[vehicle_one_id]->load << " " <<
        // initial_solution->vehicles[vehicle_two_id]->load <<endl;
		/*
        cout << "Customers to move: " << best_move_1[3] << " " << customer_to_relocate_2 << endl;
		cout <<"Pos of customers :  " << customer_to_relocate_1 << " " << best_move_2[6] << endl;
		
        for (int i = 0; i < initial_solution->vehicles[vehicle_one_id]->customers.size(); i++){
            cout << initial_solution->vehicles[vehicle_one_id]->customers[i]->id << " ";
        }
        cout << " " << endl;
        for (int i = 0; i < initial_solution->vehicles[vehicle_two_id]->customers.size(); i++){
            cout << initial_solution->vehicles[vehicle_two_id]->customers[i]->id << " ";
        }
        cout << " " << endl;
        */
		
		//cout << "Previous demands: " << initial_solution->vehicles[vehicle_one_id]->load << " " << initial_solution->vehicles[vehicle_two_id]->load << endl;
		
		
        shared_ptr<Customer> customer_one = initial_solution->vehicles[vehicle_one_id]->customers[customer_to_relocate_1];
        initial_solution->vehicles[vehicle_two_id]->insert_customer(customer_one, customer_new_pos_1);
        initial_solution->vehicles[vehicle_two_id]->customers[customer_new_pos_1]->vehicle_route = vehicle_two_id;
		initial_solution->customers[initial_solution->vehicles[vehicle_two_id]->customers[customer_new_pos_1]->id]->vehicle_route = vehicle_two_id;

        for(int i = 1; i < initial_solution->vehicles[vehicle_two_id]->customers.size() - 1; i++) {
            if(customer_to_relocate_2 == initial_solution->vehicles[vehicle_two_id]->customers[i]->id) {
                customer_two_id = i;
                break;
            }
        }

        shared_ptr<Customer> customer_two = initial_solution->vehicles[vehicle_two_id]->customers[customer_two_id];
        initial_solution->vehicles[vehicle_one_id]->insert_customer(customer_two, customer_new_pos_2);
        initial_solution->vehicles[vehicle_one_id]->customers[customer_new_pos_2]->vehicle_route = vehicle_one_id;
		initial_solution->customers[initial_solution->vehicles[vehicle_one_id]->customers[customer_new_pos_2]->id]->vehicle_route = vehicle_one_id;

        initial_solution->vehicles[vehicle_two_id]->remove_customer(customer_two_id);

        for(int i = 1; i < initial_solution->vehicles[vehicle_one_id]->customers.size() - 1; i++) {
            if(best_move_1[3] == initial_solution->vehicles[vehicle_one_id]->customers[i]->id) {
                customer_one_id = i;
                break;
            }
        }

        initial_solution->vehicles[vehicle_one_id]->remove_customer(customer_one_id);

        /*
        for (int i = 0; i < initial_solution->vehicles[vehicle_one_id]->customers.size(); i++){
            cout << initial_solution->vehicles[vehicle_one_id]->customers[i]->id << " ";
        }
        cout << " " << endl;
        for (int i = 0; i < initial_solution->vehicles[vehicle_two_id]->customers.size(); i++){
            cout << initial_solution->vehicles[vehicle_two_id]->customers[i]->id << " ";
        }
        cout << " " << endl;
        */
        // cout << initial_solution->vehicles[vehicle_one_id]->load << " " <<
        // initial_solution->vehicles[vehicle_two_id]->load <<endl;

        // initial_solution->is_feasible();

        return initial_solution;
    } else {
        return initial_solution;
    }
}



