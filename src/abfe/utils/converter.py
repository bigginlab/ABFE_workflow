from numbers import Number
from typing import Union

import numpy as np
from scipy import constants as const

boltzman_const = const.Boltzmann * const.Avogadro  # J/(K)^-1 * mol^-1 = J/(K*mol)^-1

def IC50_to_dG(ic50s: Union[Number, np.array], temperature: float = 298) -> np.array:
    """_summary_

    Args:
        ic50s (Union[Number, np.array]): concentration in mol
        temperature (float, optional): _description_. Defaults to 298.

    Returns:
        np.array: kJ/mol
    """
    beta = 1 / (temperature * boltzman_const / 1000)  # 1/kJ/mol
    dGs = (1 / beta) * np.log(ic50s)

    return dGs


print(IC50_to_dG(np.array([100/
                           /0, 100, 10, 1, 0.46, 0.01, 0.001]))/4.18)


