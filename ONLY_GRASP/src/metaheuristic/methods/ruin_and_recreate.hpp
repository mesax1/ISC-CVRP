//
// Created by mesar on 8/12/2021.
// Adapted from FILO metaheuristic by acco93
//
//
#include "../model/solution.h"
#include <algorithm>
#include <cmath>
#include <list>
#include <queue>
#include <unordered_map>
#include <random>
#include <unordered_set>
#include <assert.h>

class RuinAndRecreate {

    const Instance& instance;
    std::mt19937& rand_engine;
    std::uniform_int_distribution<int> boolean_dist;
    std::uniform_int_distribution<int> customers_distribution;
    std::uniform_int_distribution<int> rand_uniform;

public:

    RuinAndRecreate( const Instance& instance_, std::mt19937& rand_engine_) : instance(instance_),
                                                                              rand_engine(rand_engine_),
                                                                              boolean_dist(std::uniform_int_distribution(0, 1)),
                                                                              customers_distribution(instance.customers_begin, instance.customers_end - 1),
                                                                              rand_uniform(0, 3){

    }

    std::vector<int> removed = std::vector<int>();
    int apply(std::shared_ptr <Solution>& solution, std::vector<int>& omega) {

        //auto removed = std::vector<int>();


        std::uniform_int_distribution<int> routes_distribution (0, solution->vehicles.size()-1);

        removed.clear();
        auto routes = std::unordered_set<int>();
        auto i_route = routes_distribution(rand_engine);
        std::uniform_int_distribution<int> initial_customer_distribution (1, solution->vehicles[i_route]->customers.size()-2);
        auto seed = initial_customer_distribution(rand_engine);
        const auto N =  omega[seed];

        auto curr = seed;
        auto route = i_route;
        auto next_route = route;

        for(auto n = 0; n < N; n++) {


            assert(curr != instance.depot);
            assert(solution->vehicles[route]->customers[curr]->id != instance.depot);

            auto next = Solution::dummy_vertex;
            next_route = route;


            routes.insert(route);

            for (int i=0; i < solution->vehicles[route]->customers.size(); i++){
                //cout << solution->vehicles[route]->customers[i]->id << " ";
            }
            //cout << endl;

            if(solution->vehicles[route]->customers.size() > 3 && boolean_dist(rand_engine)) {
                // move within the current route

                if(boolean_dist(rand_engine)) {
                    //next = solution->get_next_vertex(curr);
                    next = curr;
                    //cout << "1_next_id: " << solution->vehicles[route]->customers[next+1]->id << endl;
                    //assert(solution->vehicles[route]->customers[next+1]->id != instance.depot);

                    if(next+1 >= solution->vehicles[route]->customers.size() -1) {
                        next=1;
                        //cout << "2_next_id: " << solution->vehicles[route]->customers[next]->id << endl;
                        assert(solution->vehicles[route]->customers[next]->id != instance.depot);
                    }
                } else {
                    next = solution->get_prev_vertex(curr);
                    //cout << "3_next_id: " << solution->vehicles[route]->customers[next]->id << " next was: " << next<< " Size will be: " << solution->vehicles[route]->customers.size()-1 << endl;


                    if(next <= 0) {
                        //next = solution->get_prev_vertex(route, next);
                        next = solution->vehicles[route]->customers.size()-3;
                        //cout << "4_next_id: " << solution->vehicles[route]->customers[next]->id << endl;
                        assert(solution->vehicles[route]->customers[next]->id != instance.depot);
                    }
                }


            } else {
                // jump to neighbor route

                if(boolean_dist(rand_engine)) {
                    if(boolean_dist(rand_engine)) {
                        for (auto m = 1u; m < instance.sorted_dist_matrix[curr].size(); m++) {
                            const auto neighbor = instance.sorted_dist_matrix[curr][m];
                            auto route_i = solution->get_route_index(neighbor);
                            if (neighbor == instance.depot || !solution->is_customer_in_solution(neighbor) ||
                                routes.count(solution->get_route_index(neighbor))) { continue; }
                            //next = neighbor;

                            for (int k = 1; k < solution->vehicles[route_i]->customers.size() - 1; k++) {

                                if (solution->vehicles[route_i]->customers[k]->id == neighbor) {
                                    //cout << neighbor << endl;
                                    next = k;
                                    next_route = route_i;
                                    //cout << "5_Changed route, next_id: " << solution->vehicles[route_i]->customers[next]->id<< endl;
                                    assert(solution->vehicles[route_i]->customers[next]->id != instance.depot);
                                    break;
                                }
                            }

                            break;
                        }
                    }
                    else {
                        for (auto m = 1u; m < instance.sorted_dist_matrix[curr].size(); m++) {

                            int route_i;
                            if (route < solution->vehicles.size() - 1) {
                                if (boolean_dist(rand_engine)) {
                                    route_i = route + 1;
                                } else if (route > 0) {
                                    route_i = route - 1;
                                } else {
                                    route_i = route + 1;
                                }
                            } else {
                                route_i = route - 1;
                            }

                            std::uniform_int_distribution<int> route_customer_distribution(1,
                                                                                           solution->vehicles[route_i]->customers.size() -
                                                                                           2);
                            const auto neighbor = route_customer_distribution(rand_engine);
                            if (neighbor == instance.depot || !solution->is_customer_in_solution(neighbor) ||
                                routes.count(route_i)) { continue; }
                            next = neighbor;
                            next_route = route_i;
                            assert(solution->vehicles[route_i]->customers[next]->id != instance.depot);
                            break;

                        }
                    }
                }
                else {
                    if(boolean_dist(rand_engine)){
                        for (auto m = 1u; m < instance.sorted_dist_matrix[curr].size(); m++) {
                            const auto neighbor = instance.sorted_dist_matrix[curr][m];
                            auto route_i = solution->get_route_index(neighbor);

                            if (neighbor == instance.depot || !solution->is_customer_in_solution(neighbor) ||
                                route_i == route) { continue; }
                            for (int k = 1; k < solution->vehicles[route_i]->customers.size() - 1; k++) {

                                if (solution->vehicles[route_i]->customers[k]->id == neighbor) {
                                    //cout << neighbor << endl;
                                    next = k;
                                    next_route = route_i;
                                    //cout << "6_Changed route, next_id: " << solution->vehicles[route_i]->customers[next]->id<< endl;
                                    assert(solution->vehicles[route_i]->customers[next]->id != instance.depot);

                                    break;
                                }
                            }
                            break;
                        }
                    }
                    else{
                        for (auto m = 1u; m < instance.sorted_dist_matrix[curr].size(); m++) {

                            int route_i;
                            if (route < solution->vehicles.size() - 1) {
                                if (boolean_dist(rand_engine)) {
                                    route_i = route + 1;
                                } else if (route > 0) {
                                    route_i = route - 1;
                                } else {
                                    route_i = route + 1;
                                }
                            } else {
                                route_i = route - 1;
                            }

                            std::uniform_int_distribution<int> route_customer_distribution(1,
                                                                                           solution->vehicles[route_i]->customers.size() -
                                                                                           2);
                            const auto neighbor = route_customer_distribution(rand_engine);
                            if (neighbor == instance.depot || !solution->is_customer_in_solution(neighbor) ||
                                routes.count(route_i == route)) { continue; }
                            next = neighbor;
                            next_route = route_i;
                            assert(solution->vehicles[route_i]->customers[next]->id != instance.depot);
                            break;

                        }
                    }

                }


            }

            assert(solution->vehicles[route]->customers[curr]->id != instance.depot);
            removed.push_back(solution->vehicles[route]->customers[curr]->id);
            solution->vehicles[route]->remove_customer(curr);

            if(solution->vehicles[route]->customers.size() <= 2) {
                solution->delete_empty_vehicles();
                if (route < next_route){
                    next_route--;
                }
            }

            if(next == Solution::dummy_vertex) {
                break;
            }


            curr = next;
            route = next_route;

        }
        //cout << "rr 1" << endl;

        std::shuffle(removed.begin(), removed.end(), rand_engine);

        //cout << "rr 2" << endl;

        //cout << endl;

        for (auto customer : removed) {

            assert(customer != instance.depot);

            std::shared_ptr<Vehicle> best_route;
            std::shared_ptr<Vehicle> dummy_route;
            auto best_where = Solution::dummy_vertex;
            auto best_cost = std::numeric_limits<float>::max();
            int mejor_id_ruta = -1;

            //cout << "rr 3" << endl;

            for(int w=0; w < solution->vehicles.size(); w++) {
                auto ruta = solution->vehicles[w];

                if (ruta->load + instance.get_demand(customer) > instance.Q) { continue; }
                //cout << "rr 4" << endl;
                for (int v=1; v < ruta->customers.size(); v++)  {

                    const auto prev = ruta->customers[v-1]->id;
                    const auto where = ruta->customers[v]->id;

                    const auto cost = -instance.get_cost(prev, where) + instance.get_cost(prev, customer) + instance.get_cost(customer, where);

                    if (cost < best_cost) {

                        best_cost = cost;
                        best_route = ruta;
                        best_where = v;
                        mejor_id_ruta = w;

                    }
                }
                //cout << "rr 5" << endl;


                const auto cost = -instance.get_cost(instance.depot, ruta->customers[ruta->customers.size()-2]->id)
                                  + instance.get_cost(ruta->customers[ruta->customers.size()-2]->id, customer)
                                  + instance.get_cost(customer, instance.depot);

                if (cost < best_cost) {

                    best_cost = cost;
                    best_route = ruta;
                    best_where = ruta->customers.size()-1;
                    mejor_id_ruta = w;

                }


            }
            //cout << "rr 3" << endl;
            if (best_route == dummy_route) {
                //solution.build_one_customer_route(customer);
                solution->vehicles.push_back(std::make_unique<Vehicle>(instance.Q));
                solution->vehicles[solution->vehicles.size()-1]->add_customer(solution->customers[0]);
                solution->vehicles[solution->vehicles.size()-1]->add_customer(solution->customers[customer]);
                solution->vehicles[solution->vehicles.size()-1]->add_customer(solution->customers[0]);
                solution->customers[customer]->vehicle_route = solution->vehicles.size()-1;
                solution->customers[customer]->isRouted = true;
            } else {
                //solution.insert_vertex_before(best_route, best_where, customer);
                //best_route->insert_customer(solution->customers[customer], best_where);
                solution->vehicles[mejor_id_ruta]->insert_customer(solution->customers[customer], best_where);
                solution->customers[customer]->vehicle_route = mejor_id_ruta;
            }

        }
        solution->compute_cost(instance);
        /*
        cout << "solution->cost: "<< solution->cost << endl;
        for (int i=0; i<solution->vehicles.size(); i++){
            //cout << "Vehicle: " << i << " cost: " << solution->vehicles[i]->cost << endl;
            for (int j=0; j<solution->vehicles[i]->customers.size(); j++){
                cout <<  solution->vehicles[i]->customers[j]->id << " ";
            }
            cout << endl;
        }
         */
        return seed;

    }



};


