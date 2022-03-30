#pragma once
#include "../model/solution.h"
#include "../model/instance.h"
#include <vector>

class GiantTour {
    public:
    std::shared_ptr <Solution> giant_tour_solution;
    Instance instance;
    
    GiantTour();
    GiantTour(int seed);
    ~GiantTour();

    bool unassignedCustomerExists();
    virtual std::shared_ptr <Solution> run(const Instance& instance, int alpha);
};

class GiantTour_RNN:public GiantTour {
    public: 
    
    GiantTour_RNN();
    GiantTour_RNN(int seed);
    ~GiantTour_RNN();

    std::shared_ptr <Solution> run(const Instance& instance, int alpha);
};

class Improved_GiantTour_RNN:public GiantTour {
    public: 
    
    Improved_GiantTour_RNN();
    Improved_GiantTour_RNN(int seed);
    ~Improved_GiantTour_RNN();

    std::shared_ptr <Solution> run(const Instance& instance, int alpha);
};

class GiantTour_RNI:public GiantTour {
    public: 
    
    GiantTour_RNI();
    GiantTour_RNI(int seed);
    ~GiantTour_RNI();

    std::shared_ptr <Solution> run(const Instance& instance, int alpha);
};

class GiantTour_RBI:public GiantTour {
    public: 
    
    GiantTour_RBI();
    GiantTour_RBI(int seed);
    ~GiantTour_RBI();

    std::shared_ptr <Solution> run(const Instance& instance, int alpha);
};