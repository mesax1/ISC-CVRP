#include "./vnd.h"
#include <iostream>

using namespace std;

shared_ptr <Solution> run_vnd(shared_ptr <Solution> initial_solution,
    const Instance& instance,
    vector<Neighborhood*>& neighborhoods,
    int& relocation_applications)
{	
	
    shared_ptr <Solution> best_solution = make_shared<Solution>(initial_solution);
	
    //float previous_cost = best_solution->compute_cost(instance);
	float previous_cost = initial_solution->cost;
    int neighborhood = 0;
    float current_cost = 0;
    bool perform_relocation = true;
    //Neighborhood* perturbation = new CrossPerturbationNeighborhood(false);
    //Neighborhood* relocation_chain = new RelocationChainNeighborhood(false);
    /*
    while(perform_relocation) {
        Solution* local_optima = relocation_chain->search_neighbors(best_solution, instance);
        // current_cost = local_optima->compute_cost(instance);
        current_cost = local_optima->cost;
        if(current_cost < previous_cost) {
            best_solution = local_optima;
            previous_cost = current_cost;
            neighborhood = 0;
        } else {
            perform_relocation = false;
            best_solution = local_optima;
        }
    }
    */
	
    while(neighborhood < neighborhoods.size()) {
		
        shared_ptr <Solution> local_optima = neighborhoods[neighborhood]->search_neighbors(best_solution, instance);
		
        // current_cost = local_optima->compute_cost(instance);
        current_cost = local_optima->cost;
        if(current_cost < previous_cost) {
            best_solution = local_optima;
            previous_cost = current_cost;
            neighborhood = 0;
        } else {
            neighborhood++;
            if(neighborhood == neighborhoods.size() - 1) {
                relocation_applications++;
            }
        }
        /*
        if ((neighborhood >= neighborhoods.size()) && (perturbation_count<100))
        {
            perturbation_count++;
            local_optima = perturbation->search_neighbors(best_solution, instance);
            current_cost = local_optima->compute_cost(instance);
            best_solution = local_optima;
                        previous_cost = current_cost;
            neighborhood = 0;
        }

        if ((neighborhood >= neighborhoods.size()) && (perturbation_count>=100))
        {
            local_optima = relocation_chain->search_neighbors(best_solution, instance);
            current_cost = local_optima->compute_cost(instance);
            if(current_cost < previous_cost){
            best_solution = local_optima;
                        previous_cost = current_cost;
            neighborhood = 0;
            perturbation_count = 0;
        }
        }
        */
    }
	best_solution->cost = best_solution->compute_cost(instance);
    best_solution->delete_empty_vehicles();
    /*
    cout << "local optima cost: " << best_solution->cost << endl;
    best_solution->cost= best_solution->compute_cost(instance);
    cout << "calculated_cost: " << best_solution->cost << endl;
     */
    return best_solution;
}

/*run_one_vnd_return run_one_vnd(Solution& initial_solution, vector<Neighborhood*>& neighborhoods){
    Solution best_solution = Solution(initial_solution);
    int neighborhood = 0;
    unordered_map<string, int> vnd_improvements;
    for(auto& item : neighborhoods){
        vnd_improvements[item->id] = 0;
    }
    while (neighborhood < neighborhoods.size()){
        Solution local_optima = neighborhoods[neighborhood]->search_neighbors(best_solution);
        if(local_optima.compute_cost() < best_solution.compute_cost()){
            best_solution = local_optima;
            vnd_improvements[neighborhoods[neighborhood]->id]++;
            neighborhood = 0;
        }else{
            neighborhood++;
        }
    }
    return {best_solution, vnd_improvements};
}*/