from gnn.input_parser import *
from gnn.preprocessing import *

import numpy as np
import pandas as pd
from os import listdir
from os.path import isfile, join

from datetime import date



def split_into_deltas_quantiles(final_costs, good_quantile, bad_quantile):
    """
        input: list of tuples of the form
        (neighborhood, iteration, initial_solution, final_solution)
    """
    Z = np.array(final_costs)
    labels = np.zeros(len(Z))
    bks = Z.min()
    
    good_sol = np.quantile(Z, good_quantile)
    bad_sol = np.quantile(Z, bad_quantile)
    
    
    for i in range(len(Z)):
        if Z[i] <= good_sol:
            labels[i] = 1
        elif Z[i] > bad_sol:
            labels[i] = 0
        else:
            labels[i] = 2  
    
    good_delta =  (good_sol - bks) / bks
    bad_delta =  (bad_sol - bks) / bks
    
    return  labels, bks, good_delta, bad_delta

def split_into_good_bad(dataset, delta, fit_encoder = True):
    """
        input: list of tuples of the form
        (neighborhood, iteration, initial_solution, final_solution)
    """
    X = [x[3] for x in dataset]
    Z = [x[4] for x in dataset]
    
    bks = np.array(Z).min()
    
   
    for i in range(len(Z)):
        if Z[i] < bks + bks*delta:
            Z[i] = 1
        else:
            Z[i] = 0  
            
    return Z

def classify_solutions(n_samples, train_samples, file_seed, instance_full_name, feature_type):
    """
    quantile_good = 0.18
    quantile_bad = 0.85
    depth_trees = 4
    select_thresh = 10
    class_weights=None
    limit_train_size = 500    
    #class_weights={0:1, 1:3}
    """
    instance_name = instance_full_name.split('.')[0]
    delta_good = 0
    seed_name = str(file_seed)
    
    
    print("reading dataset")
    
    all_feature_types = ["tsp", "vrp", "client-pairs", "tsp-vrp", "tsp-client-pairs", "vrp-client-pairs", "tsp-vrp-client-pairs"]
    features = all_feature_types[feature_type]
    try:
        grasp_mlp = pd.read_csv('results/evolution/rr_grasp_mlp_with_final_evolution_' + instance_name + '_s' + seed_name + '.log', header=None)
        only_grasp = pd.read_csv('results/evolution/rr_complete_grasp_with_final_evolution_' + instance_name + '_s' + seed_name + '.log', header=None)
    except:
        try:
            grasp_mlp = pd.read_csv('results/evolution/grasp_mlp_with_final_evolution_' + instance_name + '_s' + seed_name + '.log', header=None)
            only_grasp = pd.read_csv('results/evolution/complete_grasp_with_final_evolution_' + instance_name + '_s' + seed_name + '.log', header=None)
        except:
            grasp_mlp = pd.read_csv('results/evolution/grasp_mlp_evolution_' + instance_name + '_s' + seed_name + '.log', header=None)
            only_grasp = pd.read_csv('results/evolution/complete_grasp_evolution_' + instance_name + '_s' + seed_name + '.log', header=None)

    
    print('results/evolution/grasp_mlp_evolution_' + instance_name + '_s' + seed_name + '.log')

    grasp_mlp_times = []
    grasp_mlp_costs = []    
    for w in range(len(grasp_mlp)):
        a = grasp_mlp.iloc[w, 0].split(": ")
        if a[0] == 'time':
            grasp_mlp_times.append(int(a[1]))
        elif a[0] == "cost":
            grasp_mlp_costs.append(int(a[1]))
        
    #print(grasp_mlp_times)
    #print(grasp_mlp_costs)
    
    only_grasp_times = []
    only_grasp_costs = []    
    for w in range(len(only_grasp)):
        a = only_grasp.iloc[w, 0].split(": ")
        if a[0] == 'time':
            only_grasp_times.append(int(a[1]))
        elif a[0] == "cost":
            only_grasp_costs.append(int(a[1]))
        
    #print(only_grasp_times)
    #print(only_grasp_costs)
    
    
    df_final = pd.DataFrame([[0, 0, 0, 0, 0, 0]], columns = ['instancia', 'seed', 'only_grasp_best_solution', 'grasp_mlp_best_solution', 'gap',
            'execution_time' ])
    
    df_final['instancia'] = instance_name
    df_final['seed'] = file_seed
    df_final['only_grasp_best_solution'] = only_grasp_costs[-1]
    df_final['grasp_mlp_best_solution'] = grasp_mlp_costs[-1]
    df_final['execution_time'] = only_grasp_times[-1]
    df_final['gap'] = (df_final['grasp_mlp_best_solution'] - df_final['only_grasp_best_solution']) / df_final['only_grasp_best_solution']
    
    print(instance_name + "_s" + seed_name + " gap: " + str(df_final['gap'][0]))
    return df_final, only_grasp_times, only_grasp_costs, grasp_mlp_times, grasp_mlp_costs
    
    
if __name__ == "__main__":
    
    """
    n_samples = int(sys.argv[1])
    train_samples = int(sys.argv[2])
    seed = int(sys.argv[3])
    instance_name = sys.argv[4]
    feature_type = int(sys.argv[5])
    classify_solutions(n_samples, train_samples, seed, instance_name, feature_type)
    """
    quantile_good = 0.15
    quantile_bad = 0.50
    limit_train_size = 300
    depth_trees = 4
    select_thresh = 10
    threshold = 0.5
    class_weights=None
    use_undersample = False
    use_smote = False
    #df_final = classify_solutions(3000, 500, 2, 'X-n449-k29.vrp', feature_type=3)
    df3 = pd.DataFrame(None)
    df4 = pd.DataFrame(None)
    df_best = pd.DataFrame(None)
    df_worst = pd.DataFrame(None)
    df_complete_grasp_best = pd.DataFrame(None)
    df_importances = pd.DataFrame(None)
    df_imp_mean = pd.DataFrame(None)
    df_relocations = pd.DataFrame(None)
    
    df_only_grasp_best = pd.DataFrame(None)
    df_ml_only_grasp_best = pd.DataFrame(None)
    df_imp_grasp_best = pd.DataFrame(None)
    df_imp_ml_best = pd.DataFrame(None)
    
    instances_path = 'instances/'
    all_instances = [f for f in listdir(instances_path) if isfile(join(instances_path, f))]
    all_instances.sort()
    instances_names = all_instances[107:]
    #for j in range(6, 24):
    #instances_names = instances_names[1:]
    for instance_name in instances_names:
        #try:
    
        
            #instance_name = all_instances[j]
            
            
            
            
            #class_weights={0:1.2, 1:1}
            
            #class_weights = "balanced"
            df2 = pd.DataFrame(None)
            temp_importances = pd.DataFrame(None)
            
            for i in range(1, 11):
                try:
                    df, only_grasp_times, only_grasp_costs, grasp_mlp_times, grasp_mlp_costs = classify_solutions(1000, 300, i, instance_name, feature_type=3)
                    df2 = df2.append(df, ignore_index=True)
                    
                except:
                    print("Error with instance: " + instance_name + " seed: " + str(i))
            
            
            df3 = df3.append(df2, ignore_index=True)
            
            columns = ['instancia', 'seed', 'only_grasp_best_solution', 'grasp_mlp_best_solution', 'gap', 'execution_time' ]
            
            df_best = df_best.append(df2[df2.grasp_mlp_best_solution == min(df2.grasp_mlp_best_solution)])
            df_worst = df_worst.append(df2[df2.grasp_mlp_best_solution == max(df2.grasp_mlp_best_solution)])
            df_complete_grasp_best = df_complete_grasp_best.append(df2[df2.only_grasp_best_solution == min(df2.only_grasp_best_solution)])
            df_complete_grasp_worst = df_complete_grasp_best.append(df2[df2.only_grasp_best_solution == max(df2.only_grasp_best_solution)])
            
            
            df2.loc['mean'] = df2.mean()
            
            df2['instancia']['mean'] = instance_name
            df4 = df4.append(df2.loc['mean'], ignore_index=True)
            
            
            
        #except:
            print("Error in file: "+ instance_name)
        
    
    today = date.today()
    dt_string = today.strftime("%d-%m-%Y")
    #df3.to_csv('resultados_metricas/resultados_completos_vrp_tsp_1000_'+dt_string+'.csv', index=False)            
    df4.to_csv('results/consolidated_results/resultados_promedio_grasp_limited_time_'+dt_string+'.csv', index=False)        
    df_best.drop_duplicates(subset='instancia', keep="last", inplace=True)
    df_best.to_csv('results/consolidated_results/resultados_grasp_mlp_best_'+dt_string+'.csv', index=False)        
    df_worst.drop_duplicates(subset='instancia', keep="last", inplace=True)
    df_worst.to_csv('results/consolidated_results/resultados_grasp_mlp_worst_'+dt_string+'.csv', index=False)        
    df_complete_grasp_best.drop_duplicates(subset='instancia', keep="last", inplace=True)
    df_complete_grasp_best.to_csv('results/consolidated_results/resultados_only_grasp_best_'+dt_string+'.csv', index=False)        
    df_complete_grasp_worst.drop_duplicates(subset='instancia', keep="last", inplace=True)
    df_complete_grasp_worst.to_csv('results/consolidated_results/resultados_only_grasp_worst_'+dt_string+'.csv', index=False)        
    
           
    
    df_consolidated = pd.DataFrame(None)
    
    
    df_consolidated['instancia'] = df4['instancia']
    df_consolidated['only_grasp_best'] = df_complete_grasp_best['only_grasp_best_solution'].to_list()
    df_consolidated['grasp_mlp_best'] = df_best['grasp_mlp_best_solution'].to_list()
    df_consolidated['only_grasp_average'] = df4['only_grasp_best_solution'].to_list()
    df_consolidated['grasp_mlp_average'] = df4['grasp_mlp_best_solution'].to_list()
    df_consolidated['average_gap'] = df4['gap'].to_list()
    df_consolidated['execution_time'] = df4['execution_time'].to_list()
    
    
    df_consolidated.to_csv('results/consolidated_results/resultados_consolidados_'+dt_string+'.csv', index=False) 
    
    print("Done consolidating results")
    
    