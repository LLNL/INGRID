#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 10:04:20 2020

@author: jguterl
"""
import types
import numpy as np
def DistribFunc(F:str or types.FunctionType='x,x',**kwargs)->types.FunctionType:
        
        def get_lambdafunc( func: str) -> types.FunctionType:
            """
            Create a function from a string input.
    
            Will be used to generate a poloidal or radial transformation
            function.
    
            Parameters
            ----------
            func : str
                An expression to generate a function from.
    
            Returns
            -------
                A function generated from the str input.
    
            Examples
            --------
            When calling method ``get_func`` the user must provide a **string**
            with the following general format:
    
            .. math::
    
                x, f(x)
    
            That is, the dependent variable and expression to evaluate are
            separated with a comma character.
    
            The `Sympy` library supports most expressions that can be generated with
            Python natively. See the `Sympy` documentation for advanced features.
    
            Defining a function representing f(x) = x ^ 2:
    
            >>> func = 'x, x ** 2'
            >>> f = MyTopologyUtils.get_func(func)
            >>> f(np.array([0, 1, 2, 3]))
            array([0, 1, 4, 9])
    
            Defining a function representing f(x) = exp(x)
    
            >>> func = 'x, exp(x)'
            >>> f = MyTopologyUtils.get_func(func)
            >>> f(np.array([0, 1, 2, 3]))
            array([ 1.        ,  2.71828183,  7.3890561 , 20.08553692])
    
            """
    
            def make_sympy_func(var, expression):
                import sympy as sp
                _f = sp.lambdify(var, expression, 'numpy')
                return _f
    
            f_str_raw = func
    
            f_str_raw = f_str_raw.replace(' ', '')
            delim = f_str_raw.index(',')
    
            var = f_str_raw[0: delim]
            expression = f_str_raw[delim + 1:]
    
            func = make_sympy_func(var, expression)
            # TODO: Check Normalization of the function to 1
            return func
    
        def CheckLambdaFunction(expression: str, Verbose: bool = False) -> bool:
            """
            Check if a str is in the correct format for method ``get_func``
    
            Parameters
            ----------
            expression : str
                Expression to check.
    
            Verbose : bool, optional
                Print full output to terminal. Default to False
    
            Returns
            -------
                True if expression is valid. False otherwise.
            """
            ExpValid = False
            try:
                com = 'lambda {}:{}'.format(expression.split(',')[0], expression.split(',')[1])
                if Verbose:
                    print(com)
                eval(com)
                ExpValid = True
            except:
                ExpValid=False
            return ExpValid
        
        
        if type(F)==str:
            if F=='exp':
                if kwargs.get('alpha') is None:
                    alpha=1
                else:
                    alpha=kwargs.get('alpha')
                def Func(x):
                    return (1-np.exp(-(x)/alpha))/(1-np.exp(-1/alpha))
            elif F=='mexp':
                if kwargs.get('alpha') is None:
                    alpha=1
                else:
                    alpha=kwargs.get('alpha')
                def Func(x):
                    return 1-((1-np.exp(-((1-x)/alpha)))/(1-np.exp(-1/alpha)))    
        
            elif F=='pow':
                if kwargs.get('alpha') is None:
                    alpha=1
                else:
                    alpha=kwargs.get('alpha')
                def Func(x):
                    return x**alpha
            elif F=='mpow':
                if kwargs.get('alpha') is None:
                    alpha=1
                else:
                    alpha=kwargs.get('alpha')
                def Func(x):
                    return 1-x**alpha
            else:
                if CheckLambdaFunction(F):    
                    Func=get_lambdafunc(F)
                else:
                    raise ValueError('Unable to parse expression entry "{}".'.format(F))
                    Func=None
        else:
            Func=F
        
        if not callable(Func):
            Func=None
            
        return Func 