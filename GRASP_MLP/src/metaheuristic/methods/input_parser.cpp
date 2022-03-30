#include "../model/instance.h"
#include "../model/solution.h"
#include "input_parser.h"
#include <boost/algorithm/string.hpp>

using namespace std;

pair<string, string> line_to_field_value(string& line, const string& sep=":"){
	vector<string> result;
	boost::split(result, line, boost::is_any_of(":"));
	// Falta el strip?
	string field_name = result[0];
	string value_str = result[1];
	
	return make_pair(field_name, value_str);
}

shared_ptr <Solution> read_solution(vector<string>& lines, int curr_line, Instance& instance){
	int n_vehicles = stoi(lines[curr_line]);

    std::shared_ptr <Solution> solution = make_shared<Solution>(instance, n_vehicles);
	vector<int> customers_served;
	vector<string> result;
	
	for(int i = 2; i < n_vehicles + 2; i++){
		boost::split(result, lines[i + curr_line], [](unsigned char x){return std::isspace(x);});
		customers_served = vector<int>(result.size());
		//Insert depot at start of route
		//customers_served[0] = 0;
		
		for(int j = 0; j < result.size(); j++){
			if (result[j] != ""){
			
				customers_served[j] = stoi(result[j]);
			}
		}

		//Insert depot at end of route
		//customers_served[customers_served.size()-1] = 0;
		
		for(const auto& customer_id : customers_served){
            std::shared_ptr <Customer> customer = solution->customers[customer_id];
			customer->isRouted = true;
			customer->vehicle_route = i-2;
			solution->vehicles[i - 2]->add_customer(customer);
		}
	}
	
	solution->compute_cost(instance);
	return solution;
}

vector<SampleInfo*> read_dataset(Instance& instance, string& dataset_file){
	ifstream f{dataset_file};
	string line;
	vector<string> lines;
	
	if(!f){
		cout << "File not found" << endl;
	}else{
		while(!f.eof()){
			getline(f, line);
			// Falta el strip?
			lines.push_back(line);
		}
	}
	
	f.close();
	
	vector<SampleInfo*> sample_infos;
	SampleInfo* sample_info;
	int i = 0;
	
	string neighborhood;
	int sample_number;
    std::shared_ptr <Solution> initial_solution;
    std::shared_ptr <Solution> tsp_solution;
	double final_cost;
	
	while(i < lines.size() && lines[i] != "" && lines[i][0] != '-'){
		neighborhood = line_to_field_value(lines[i]).second;
		sample_number = stoi(line_to_field_value(lines[i + 1]).second);
		i += 4;
		
		
		tsp_solution = read_solution(lines, i, instance);
		i += tsp_solution->n_vehicles + 3;
		
		
		initial_solution = read_solution(lines, i, instance);
		i += initial_solution->n_vehicles + 3;
		
		final_cost = stod(lines[i]);
		i++;
		
		sample_info = new SampleInfo();
		sample_info->neighborhood = neighborhood;
		sample_info->sample_number = sample_number;
		sample_info->tsp_solution = tsp_solution;
		sample_info->initial_solution = initial_solution;
		sample_info->final_cost = final_cost;
		sample_infos.push_back(sample_info);
	}
	
	return sample_infos;
}

vector<SampleInfo*> read_select_vnd_dataset(Instance& instance, string& dataset_file, int number_of_vnds){
    ifstream f{dataset_file};
    string line;
    vector<string> lines;

    if(!f){
        cout << "File not found" << endl;
    }else{
        while(!f.eof()){
            getline(f, line);
            // Falta el strip?
            lines.push_back(line);
        }
    }

    f.close();

    vector<SampleInfo*> sample_infos;
    SampleInfo* sample_info;
    int i = 0;

    int neighborhood;
    int sample_number;
    std::shared_ptr <Solution> initial_solution;
    std::shared_ptr <Solution> tsp_solution;
    float final_cost;
    vector<int> neighborhoods = {};
    vector<float> final_costs = {};

    for (int j=0; j < number_of_vnds; j++){
        neighborhoods.push_back(j);
        final_costs.push_back(0.0);
    }

    while(i < lines.size() && lines[i] != "" && lines[i][0] != '-') {
        sample_number = stoi(line_to_field_value(lines[i]).second);
        i += 2;


        tsp_solution = read_solution(lines, i, instance);
        i += tsp_solution->n_vehicles + 3;


        initial_solution = read_solution(lines, i, instance);
        i += initial_solution->n_vehicles + 2;

        sample_info = new SampleInfo();
        for (int j=0; j < number_of_vnds; j++){
            neighborhood = stoi(line_to_field_value(lines[i]).second);
            i += 2;
            final_cost = stof(lines[i]);
            i++;
            //sample_info->neighborhoods[neighborhood] = neighborhood;
            final_costs[neighborhood]= final_cost;
        }
        sample_info->final_costs= final_costs;

        sample_info->sample_number = sample_number;
        sample_info->tsp_solution = tsp_solution;
        sample_info->initial_solution = initial_solution;
        sample_infos.push_back(sample_info);
    }

    return sample_infos;
}


vector<SampleInfo*> read_validation_dataset(Instance& instance, string& dataset_file){
	ifstream f{dataset_file};
	string line;
	vector<string> lines;
	
	if(!f){
		cout << "File not found" << endl;
	}else{
		while(!f.eof()){
			getline(f, line);
			// Falta el strip?
			lines.push_back(line);
		}
	}
	
	f.close();
	
	vector<SampleInfo*> sample_infos;
	SampleInfo* sample_info;
	int i = 0;
	
	string neighborhood;
	int sample_number;
	shared_ptr <Solution> initial_solution;
    shared_ptr <Solution> tsp_solution;
	double final_cost;
	
	while(i < lines.size() && lines[i] != "" && lines[i][0] != '-'){
		sample_number = stoi(line_to_field_value(lines[i]).second);
		i += 2;
		
		tsp_solution = read_solution(lines, i, instance);
		i += tsp_solution->n_vehicles + 3;
		
		initial_solution = read_solution(lines, i, instance);
		i += initial_solution->n_vehicles + 2;
		
		final_cost = 0.0;
        neighborhood = "0";
		
		sample_info = new SampleInfo();
		sample_info->neighborhood = neighborhood;
		sample_info->sample_number = sample_number;
		sample_info->tsp_solution = tsp_solution;
		sample_info->initial_solution = initial_solution;
		sample_info->final_cost = final_cost;
		sample_infos.push_back(sample_info);
	}
	
	return sample_infos;
}

vector<int> read_classification(string& dataset_file){
	ifstream f{dataset_file};
	string line;
	vector<string> lines;
	
	if(!f){
		cout << "File not found" << endl;
	}else{
		while(!f.eof()){
			getline(f, line);
			// Falta el strip?
			lines.push_back(line);
		}
	}
	
	f.close();
	
	int i = 0;
	int current_class;
    vector<int> classification;
	
	while(i < lines.size() && lines[i] != "" && lines[i][0] != '-'){
		
		current_class = stoi(lines[i]);
		i++;
        classification.push_back(current_class);
	}
	
	return classification;
}
