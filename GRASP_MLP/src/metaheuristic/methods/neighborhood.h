#pragma once
#include "../model/solution.h"

class Neighborhood{
    public:
    bool first_improvement;
    std::string id;

    Neighborhood(bool first_improvement = false, const std::string& id = "");
    ~Neighborhood();

    virtual std::shared_ptr <Solution> search_neighbors(std::shared_ptr <Solution> initial_solution, const Instance& instance);
};

/**
    Perform inter-route and intra-route Swap (exchange operation)

    Currently it searches the entire neighborhood and applies the
    best possible movement.
*/
class SwapNeighborhood: public Neighborhood{
    public:
    SwapNeighborhood();

    SwapNeighborhood(bool first_improvement = false);

    ~SwapNeighborhood();

    std::shared_ptr <Solution> search_neighbors(std::shared_ptr <Solution> initial_solution, const Instance& instance);
    std::shared_ptr <Solution> apply_best_move(std::shared_ptr <Solution> initial_solution, std::vector<int>& best_move);
};

/**
    Perform inter-route and intra-route Swap(2,2) (exchange operation)

    Currently it searches the entire neighborhood and applies the
    best possible movement.
*/
class SwapTwoNeighborhood: public Neighborhood{
    public:
    SwapTwoNeighborhood();

    SwapTwoNeighborhood(bool first_improvement = false);

    ~SwapTwoNeighborhood();

    std::shared_ptr <Solution> search_neighbors(std::shared_ptr <Solution> initial_solution, const Instance& instance);
    std::shared_ptr <Solution> apply_best_move(std::shared_ptr <Solution> initial_solution, std::vector<int>& best_move);
};

/**
    Perform inter-route and intra-route Relocate (insertionoperation)

    Currently it searches the entire neighborhood and applies the
    best possible movement.
*/
class RelocateNeighborhood: public Neighborhood{
    public:
    RelocateNeighborhood();

    RelocateNeighborhood(bool first_improvement = false);

    ~RelocateNeighborhood();

    std::shared_ptr <Solution> search_neighbors(std::shared_ptr <Solution> initial_solution, const Instance& instance);
    std::shared_ptr <Solution> apply_best_move(std::shared_ptr <Solution> initial_solution, std::vector<int>& best_move);
};

/**
    Perform intra-route 2-opt move

    Currently it searches the entire neighborhood and applies the
    best possible movement.
*/
class TwoOptIntraNeighborhood: public Neighborhood{
    public:
    TwoOptIntraNeighborhood();

    TwoOptIntraNeighborhood(bool first_improvement = false);

    ~TwoOptIntraNeighborhood();

    std::shared_ptr <Solution> search_neighbors(std::shared_ptr <Solution> solution, const Instance& instance);
    std::shared_ptr <Solution> apply_best_move(std::shared_ptr <Solution> solution, std::vector<int>& best_move);
};

/**
    r1 = w w a - b x x
    r2 = y y c - d z z

    r1 = w w a d z z
    r2 = reverse(b x x) + reverse(c y y)
    r2 = x x b c y y

    r1 = r1[None: i + 1] + r2[j+1: None]
    r2 = r1[None:i:-1] + r2[j:None:-1]
    r1_sgm2[0] -> None
    r1_sgm2[1] -> i
    r1_sgm2[2] -> -1

    r2_sgm2[0] -> j
    r2_sgm2[0] -> None
    r2_sgm2[0] -> -1
*/
class TwoOptInterNeighborhood: public Neighborhood{
    public:
    TwoOptInterNeighborhood();

    TwoOptInterNeighborhood(bool first_improvement = false);

    ~TwoOptInterNeighborhood();

    std::shared_ptr <Solution> search_neighbors(std::shared_ptr <Solution> solution, const Instance& instance);
    std::shared_ptr <Solution> apply_best_move(std::shared_ptr <Solution> solution, std::vector<int>& best_move);
};

class ThreeOptInterNeighborhood: public Neighborhood{
    public:
    ThreeOptInterNeighborhood();

    ThreeOptInterNeighborhood(bool first_improvement = false);

    ~ThreeOptInterNeighborhood();

    std::shared_ptr <Solution> search_neighbors(std::shared_ptr <Solution> solution, const Instance& instance);
    std::shared_ptr <Solution> apply_best_move(std::shared_ptr <Solution> solution, std::vector<int>& best_move);
};

class CrossExchangeNeighborhood: public Neighborhood{
	public:
	CrossExchangeNeighborhood();
	
	CrossExchangeNeighborhood(bool first_improvement = false);

    ~CrossExchangeNeighborhood();

    std::shared_ptr <Solution> search_neighbors(std::shared_ptr <Solution> solution, const Instance& instance);
    std::shared_ptr <Solution> apply_best_move(std::shared_ptr <Solution> solution, std::vector<int>& best_move);
};



/**
    Perform inter-route Relocation chain

    Currently it searches the X nearest neighbors and applies the
    best possible movement.
*/
class RelocationChainNeighborhood_LimitNodes: public Neighborhood{
    public:
    RelocationChainNeighborhood_LimitNodes();

    RelocationChainNeighborhood_LimitNodes(bool first_improvement = false);

    ~RelocationChainNeighborhood_LimitNodes();

    std::shared_ptr <Solution> search_neighbors(std::shared_ptr <Solution> initial_solution, const Instance& instance);
    std::shared_ptr <Solution> apply_best_move(std::shared_ptr <Solution> initial_solution, std::vector<std::vector<int>>& best_moves);
};

class RelocationChainNeighborhood_LimitDepth: public Neighborhood{
    public:
    RelocationChainNeighborhood_LimitDepth();

    RelocationChainNeighborhood_LimitDepth(bool first_improvement = false);

    ~RelocationChainNeighborhood_LimitDepth();

    std::shared_ptr <Solution> search_neighbors(std::shared_ptr <Solution> initial_solution, const Instance& instance);
    std::shared_ptr <Solution> apply_best_move(std::shared_ptr <Solution> initial_solution, std::vector<std::vector<int>>& best_moves);
};

/**
    Perform inter-route and intra-route SwapStar (exchange operation)

    Currently it searches the entire neighborhood and applies the
    best possible movement.
*/
class SwapStarNeighborhood: public Neighborhood{
    public:
    SwapStarNeighborhood();

    SwapStarNeighborhood(bool first_improvement = false);

    ~SwapStarNeighborhood();

    std::shared_ptr <Solution> search_neighbors(std::shared_ptr <Solution> initial_solution, const Instance& instance);
    std::shared_ptr <Solution> apply_best_move(std::shared_ptr <Solution> initial_solution, std::vector<std::vector<int>>& best_move);
};
