from abc import ABC, abstractmethod
import sys
import importlib.util 
import numpy as np

matlab_is_present = importlib.util.find_spec('matlab')
if matlab_is_present:
    import matlab.engine


class MethodTemplate(ABC):
    """_summary_

    Args:
        ABC (_type_): _description_

    Returns:
        _type_: _description_
    """

    @abstractmethod
    def __init__(self):
        ...

    @abstractmethod
    def method(self):
        ...

    @property
    def id(self):
        return self._id
    
    @id.setter
    def id(self, id):
        self._id = id
    
    @property
    def task(self):
        return self._task
    
    @task.setter
    def task(self, task):
        self._task = task

    def get_parameters(self):
        return (((),{}),) 

class MatlabInterface():
    """_summary_
    """
    def __init__(self, matlab_function_name, add2path=[], matlab_warnings=False):
        """_summary_

        Args:
            matlab_function_name (_type_): _description_
        """
        self.matlab_function_name = matlab_function_name
        self.eng = matlab.engine.start_matlab()
        
        if not matlab_warnings:
            self.eng.eval("warning('off','all');", nargout=0)

        self.eng.eval("addpath('../src/methods')")
        self.eng.eval("addpath('src/methods')")
        
        for path in add2path:
            self.eng.eval("addpath('"+path+"')")
            self.eng.eval("addpath(genpath('"+path+"'))")
        # sys.path.insert(0, os.path.abspath('../src/methods/'))

    def matlab_function(self, signal, *params,**kwargs):
        """_summary_

        Args:
            signal (_type_): _description_

        Returns:
            _type_: _description_
        """
        all_params = list((signal.copy(),*params))
        params = self.pre_parameters(*all_params)
        fun_handler = getattr(self.eng, self.matlab_function_name)
        outputs = fun_handler(*params,**kwargs)
        if len(outputs)==1:
            outputs = outputs[0].toarray()
        else:
            outputs = [output.toarray() for output in outputs]
        return np.array(outputs)
        
    def pre_parameters(self, *params):
        """_summary_

        Returns:
            _type_: _description_
        """
        params_matlab = list()
        for param in params:
            if isinstance(param,np.ndarray):
                params_matlab.append(matlab.double(vector=param.tolist()))
            if isinstance(param,list) or isinstance(param,tuple):
                params_matlab.append(matlab.double(vector=list(param)))
            if isinstance(param,float):
                params_matlab.append(matlab.double(param))
            if isinstance(param,int):
                params_matlab.append(matlab.double(float(param)))
                
        return params_matlab    

# def sigmerge(x1,x2,ratio,return_noise=False):
#     """_summary_

#     Args:
#         x1 (_type_): _description_
#         x2 (_type_): _description_
#         ratio (_type_): _description_
#         return_noise (bool, optional): _description_. Defaults to False.

#     Returns:
#         _type_: _description_
#     """
#     ex1=np.mean(np.abs(x1)**2)
#     ex2=np.mean(np.abs(x2)**2)
#     h=np.sqrt(ex1/(ex2*10**(ratio/10)))
#     sig=x1+h*x2

#     if return_noise:
#         return sig, h*x2
#     else:
#         return sig