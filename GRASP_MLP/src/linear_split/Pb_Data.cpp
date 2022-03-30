
//--------------------------------------------------------
//LIBRARY OF SPLIT ALGORITHM FOR VEHICLE ROUTING PROBLEMS
//Author : Thibaut VIDAL
//Date   : August 15th, 2015. 
//E-mail : vidalt@inf.puc-rio.br
//
//This code is distributed for research purposes.
//All rights reserved.
//--------------------------------------------------------

#include "Pb_Data.h"

Pb_Data::Pb_Data(Instance& instance, vector<shared_ptr<Customer>>& customers)
{
	// For now it's hard coded, but should be a variable of the problem or a commandline input

	// to generate different demands for each instance
	this->nbNodes = instance.n_cust;
	this->vehCapacity = instance.Q;
	this->cli = vector<Client> (this->nbNodes + 1);
	this->cli[0].demand = 0;
	this->cli[0].dreturn = 0;
	this->cli[0].dnext = 0;
	this->cli[0].index = 0;
	
	for(int i = 1; i < customers.size() - 1; i++){
		this->cli[i].index = customers[i]->id;
		this->cli[i].demand = customers[i]->demand;
		this->cli[i].dreturn = instance.dist_warehouse[customers[i]->id];
		if(i == this->nbNodes-1){
			this->cli[i].dnext = -1;
		}else{
			this->cli[i].dnext = instance.dist_matrix[customers[i]->id][customers[i+1]->id];
		}
	}
	
	this->solution = vector<int> (this->nbNodes);
	for(int i = 0; i < this->nbNodes; i++){
		this->solution[i] = -1;
	}
}

Pb_Data::~Pb_Data(void)
{}