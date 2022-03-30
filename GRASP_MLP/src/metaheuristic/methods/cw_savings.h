#pragma once
#include "../model/solution.h"
#include "../model/instance.h"
#include <random>

class ParallelSavings {
    public:
    Instance instance;
    
    ParallelSavings();
    ParallelSavings(int seed, int beta, std::vector<std::vector<int>> savings);
    ~ParallelSavings();

    std::shared_ptr <Solution> run(std::shared_ptr <Solution> solution, const Instance& instance, int beta, std::vector<std::vector<int>> deterministic_savings);
    std::shared_ptr <Solution> clarke_and_wright(const Instance& instance, std::shared_ptr <Solution> solution, float lambda, int neighbors_num, int beta, std::mt19937& rand_engine);
};