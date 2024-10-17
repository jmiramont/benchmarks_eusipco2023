if __name__ == "__main__":
    # from unittest import result
    import importlib
    from src.methods import *
    from mcsm_benchs.benchmark_utils import MethodTemplate as MethodTemplate
    import time
    import inspect
    
    # Collects the methods in the folder/ module "methods" and make a global list of them.
    print('Collecting methods to benchmark...')
    modules = dir()
    modules = [mod_name for mod_name in modules if mod_name.startswith('method_')]
    global list_of_methods # Use with caution.

    list_of_methods = list()    
    for mod_name in modules:
        mod = importlib.import_module('methods.' + mod_name)
        classes_in_mod = inspect.getmembers(mod, inspect.isclass)
        for a_class in classes_in_mod:
            method_class = getattr(mod, a_class[0])
            class_parent = method_class.__bases__[0]
            if class_parent == MethodTemplate:
                method_name = method_class().id
                print(method_name)
                list_of_methods.append(method_class())


    from mcsm_benchs.Benchmark import *
    import numpy as np
    from mcsm_benchs.ResultsInterpreter import ResultsInterpreter
    import yaml


    dictionary_of_methods = dict()
    dictionary_of_parameters = dict()

    # Select only methods for denoising.
    for method_instance in list_of_methods:
        if method_instance.task == 'denoising':
            method_id = method_instance.id
            dictionary_of_methods[method_id] = method_instance.method
            dictionary_of_parameters[method_id] = method_instance.get_parameters()

    # Parameters of the benchmark:
    # Load parameters from configuration file.
    with open("config_benchmarks.yaml", "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    
    config['methods'] = dictionary_of_methods
    config['parameters'] = dictionary_of_parameters
    config['task'] = 'denoising'

    if 'add_new_methods' in config.keys():
        if config['add_new_methods']:
            
            filename = 'results\last_benchmark_denoising'
            with open(filename + '.pkl', 'rb') as f:
                benchmark = pickle.load(f)
            benchmark.add_new_method(config['methods'],config['parameters']) 
        else:
            config.pop('add_new_methods') 
            benchmark = Benchmark(**config)    
    else:
        benchmark = Benchmark(**config)
         
    np.random.seed(0)
    start = time.time()
    my_results = benchmark.run_test() # Run the test. my_results is a nested dictionary with the results for each of the variables of the simulation.
    end = time.time()
    print("The time of execution:", end-start)

    # Save the benchmark to a file. Notice that only the methods_ids are saved.
    benchmark.save_to_file(filename = 'results/last_benchmark_denoising')