#pragma once
#include "./neighborhood.h"
#include "../model/solution.h"
#include <unordered_map>

struct run_one_vnd_return{
    Solution best_solution;
    std::unordered_map<std::string, int> vnd_improvements;
};

std::shared_ptr <Solution> run_vnd(std::shared_ptr <Solution> initial_solution, const Instance& instance, std::vector<Neighborhood*>& neighborhoods, int& relocation_applications);
run_one_vnd_return run_one_vnd(Solution& initial_solution, std::vector<Neighborhood*>& neighborhoods);