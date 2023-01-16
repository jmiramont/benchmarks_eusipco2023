from benchmark_demo.benchmark_utils import MethodTemplate, MatlabInterface
# import sys
# sys.path.append("methods")
# You must import the MethodTemplate abstract class and the MatlabInterface class.

# Create an interface with the matlab engine by passing the name of the function file 
# (without the .m extension). Then get the matlab function as:

# Paths to additional code for the method to add to Matlab path variable.
paths = ['src\methods\ssa_decomp_utils',
        '..\src\methods\ssa_decomp_utils'
        ]

mlint = MatlabInterface('slide_ssa_method', add2path=paths) 
matlab_function = mlint.matlab_function # A python function handler to the method.

class NewMethod(MethodTemplate):

    def __init__(self):
        self.id = 'sliding_ssa'
        self.task = 'denoising'
        
    def method(self, signal, nc = [], *params):
        if nc == []:
            nc = signal.ncomps

        signal_output = matlab_function(signal, nc, *params) # Only positional args.
        return signal_output

    # def get_parameters(self):            # Use it to parametrize your method.
    #     return (((25,),{}),)    
    # [signal_r, Y] = method_ssa_decomp(signal, L, n_components, epsilon)        
