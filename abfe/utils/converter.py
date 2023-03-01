import numpy as np
from scipy import constants as const
from typing import Union
from numbers import Number

boltzman_const = const.Boltzmann*const.Avogadro # J/(K)^-1 * mol^-1 = J/(K*mol)^-1

def IC50_to_dG(ic50s:Union[Number, np.array], temperature:float=298)->np.array:
    """_summary_

    Args:
        ic50s (Union[Number, np.array]): concentration in mol
        temperature (float, optional): _description_. Defaults to 298.

    Returns:
        np.array: _description_
    """
    beta = 1/(temperature*boltzman_const/1000) # 1/kJ/mol
    dGs = -(1/beta)*np.log(ic50s)
    
    return dGs
    

print(IC50_to_dG(np.array([0, 0.1*10**-9, 1*10**-9,2*10**-9,3*10**-9])))