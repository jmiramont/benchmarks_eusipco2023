from benchmark_demo.benchmark_utils import MethodTemplate, MatlabInterface
# import sys
# sys.path.append("methods")
# You must import the MethodTemplate abstract class and the MatlabInterface class.

# Create an interface with the matlab engine by passing the name of the function file 
# (without the .m extension). Then get the matlab function as:

# Paths to additional code for the method to add to Matlab path variable.
paths = ['src\methods\EM_method_utils',
        '..\src\methods\EM_method_utils'
        ]

mlint = MatlabInterface('em_method', add2path=paths, matlab_warnings=False) 
matlab_function = mlint.matlab_function # A python function handler to the method.

class NewMethod(MethodTemplate):

    def __init__(self):
        self.id = 'em_method'
        self.task = 'denoising'
        

    def method(self, signal, nc=[], *params):
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
        signal_output = matlab_function(signal, nc, *params) # Only positional args.
        return signal_output
        
    # xr = em_method(x,Ncomp,M,L,c,cl,step_Nx,stepg,seuil,return_comps)    
    # def get_parameters(self):            # Use it to parametrize your method.
    #     return (([], [], [], [], [], [], [], [], True,),)  # # Use this to return all
        