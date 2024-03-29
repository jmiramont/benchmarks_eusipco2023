from benchmark_demo.benchmark_utils import MethodTemplate, MatlabInterface
# import sys
# sys.path.append("methods")
# You must import the MethodTemplate abstract class and the MatlabInterface class.

# Create an interface with the matlab engine by passing the name of the function file 
# (without the .m extension). Then get the matlab function as:

# Paths to additional code for the method to add to Matlab path variable.
paths = ['src\methods\FRI_method_utils',
        '..\src\methods\FRI_method_utils'
        ]

mlint = MatlabInterface('fri_method', add2path=paths) 
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
        self.id = 'fri_method'
        
        # self.task = 'denoising'
        # self.task = 'component_denoising'
        # self.task = 'inst_frequency'
        self.task = task
        

    def method(self, signal, nc=[], *params):
        """_summary_

        Args:
            signals (_type_): _description_
            params (_type_): _description_

        Returns:
            _type_: _description_
        """
        if nc == []:
            nc = signal.total_comps

        # Matlab function: xr = fri_method(x,Ncomp,M,L,Pnei,M0,Method)
        signal_output = matlab_function(signal, nc, *params) # Only positional args.
        return signal_output

    def get_parameters(self):            
        if self.task == 'component_denoising':            
            return (([], [], [], [], [], [], True,),)  # Return components

        if self.task == 'inst_frequency':            
            return (([], [], [], [], [], [], [], True,),)  # Return inst. freq.            

        if self.task == 'denoising':
            return (((),{}),)     

