#ifndef DATA_H
#define DATA_H

#include "../metaheuristic/model/instance.h"
#include "../metaheuristic/model/solution.h"

#include <stdlib.h>
#include <stdio.h> 
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream> 
#include <time.h>
#include <math.h>

using namespace std ;

// Solver types:
// BELLMAN --> Classic O(nB) approach -- B is the average size of a route
// BELLMAN_SOFT --> Classic O(nB') approach with possible penalized capacity excess and unlimited fleet -- B' is the average size of a route (need a limit on capacity excess otherwise B' = n)
// BELLMAN_BOUNDED --> Classic O(nBK) approach, with limited fleet
// LINEAR --> Efficient O(n) approach for Split
// LINEAR_SOFT --> Efficient O(n) approach, with possible penalized capacity excess
// LINEAR_BOUNDED --> Efficient O(n) approach, with bounded fleet

struct Client
{
	int index ;
	double demand ;
	double dreturn ;
	double dnext ;
};

class Pb_Data
{
public:

/* PROBLEM DATA */

// number of nodes
int nbNodes ;

// vehicle capacity
double vehCapacity ;

// vector of clients, i.e., locations to be visited in the TSP
vector < Client > cli ;

// Starting time of the optimization
clock_t time_StartComput ;

// End time of the optimization
clock_t time_EndComput ;

/* SOLUTION STRUCTURE */

// list of indices of clients where a route starts, -1 if not filled
vector < int > solution ;

// cost of the solution
double solutionCost ;

// number of routes in the solution
int solutionNbRoutes ;

/* METHODS TO TEST AND EXPORT THE FINAL SOLUTION */

// Method to display a solution
void printSolution()
{ 
	cout << endl ;
	cout << "------------------------------------" << endl ;
	cout << "SOLUTION COST : " << std::setprecision(12) << solutionCost << endl ;
	cout << "NB ROUTES : " << solutionNbRoutes << endl ;
	cout << "SOLUTION : [ " ;
	for (int i=0 ; i < solutionNbRoutes ; i++)
		cout << solution[i] << " " ; 
	cout << "]" << endl ;
	cout << "------------------------------------" << endl ;
	cout << endl ;
}

// Method to test a solution, verify the solution cost 
void checkSolution()
{ 
	double costTotal = 0 ;
	int begin, end ;
	double load, distance ;

	if (solution[0] != 1)
		cout << " ERROR : First route should start with customer 0" << endl ;

	for (int i=0 ; i < solutionNbRoutes ; i++)
	{
		begin = solution[i] ;
		if (i < solutionNbRoutes-1)
			end = solution[i+1]-1 ;
		else
			end = nbNodes ;

		load = 0 ;
		for (int j = begin ; j <= end ; j++)
			load += cli[j].demand ;

		distance = cli[begin].dreturn + cli[end].dreturn ;
		for (int j = begin ; j < end ; j++)
			distance += cli[j].dnext ;

		costTotal += distance;
	}

	//cout << "Cost reported by the Split algorithm : " << solutionCost << endl ;
	//cout << "Cost evaluated by the solution checker : " << costTotal << endl ;

	if (costTotal > solutionCost + 0.0001 || costTotal < solutionCost - 0.0001)
		cout << "ERROR : Solution checker does not find the same solution cost" << endl ;
}

// Constructor
Pb_Data(Instance& instance, vector<shared_ptr<Customer>>& customers);

~Pb_Data(void);
};

#endif

