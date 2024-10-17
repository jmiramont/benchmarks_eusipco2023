from mcsm_benchs.benchmark_utils import MethodTemplate
from mcsm_benchs.MatlabInterface import MatlabInterface
import os

# import sys
# sys.path.append("methods")
# You must import the MethodTemplate abstract class and the MatlabInterface class.

# Create an interface with the matlab engine by passing the name of the function file 
# (without the .m extension). Then get the matlab function as:

# Paths to additional code for the method to add to Matlab path variable.
paths = [   os.path.join('src','methods','EM_method_utils'),
            os.path.join('..','src','methods','EM_method_utils')
        ]

mlint = MatlabInterface('em_method', add2path=paths, matlab_warnings=False) 
matlab_function = mlint.matlab_function # A python function handler to the method.

# Parameters of the benchmark:
# Load parameters from configuration file.
import yaml
try:
    with open('src\methods\config_tasks.yaml', "r") as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    task =config['task']
except:
    task = 'denoising'


class NewMethod(MethodTemplate):

    def __init__(self):
        self.id = 'em_method'
        self.task = task
        

    def method(self, signal, nc=[], *params, **kwargs):
        """_summary_

        Args:
            signals (_type_): _description_
            params (_type_): _description_

        Returns:
            _type_: _description_
        """
        if nc==[]:
            nc = signal.total_comps

        # Matlab function: xr = 
        signal_output = matlab_function(signal, nc, *params, **kwargs) # Only positional args.
        return signal_output
        
    # xr = em_method(x,Ncomp,M,L,c,return_comps, return_freq)
    # xr = em_method(x,Ncomp,M,L,c, step_r, step_v, return_comps, return_freq) 
    # def get_parameters(self):            # Use it to parametrize your method.
    #     return (([], [], [], [], [], [], [], [], True,),)  # # Use this to return all
    def get_parameters(self):            
        if self.task == 'component_denoising':            
            return (([],[],[],[],[],[],True,),
                    # ([],[],[],[],1,100,True,),
                    )  # Return components

        if self.task == 'inst_frequency':            
            return (([],[],[],[],[],[],[],True,),
                    # ([],[],[],[],1,100,[],True,),
                    )  # Return inst freq.            

        if self.task == 'denoising':
            return ((([],[],[],[],[],[]),{}),
                    # (([],[],[],[],1,100),{}),
                    )        