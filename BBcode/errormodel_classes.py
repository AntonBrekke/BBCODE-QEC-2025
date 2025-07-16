import functools
from typing import Tuple, Optional
import numpy as np
from panqec.codes import StabilizerCode
from panqec.error_models import BaseErrorModel
from panqec.bpauli import pauli_to_bsf
import random

from scipy.integrate import quad 
from scipy.optimize import newton


th = np.sqrt(np.pi)/2

def gaussian(x, sigma):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-x**2/(2*sigma)**2)

def root_func(sigma, p):
    return quad(gaussian, -th, th, args=(sigma,))[0] + p - 1

def fast_choice(options, probs, rng=None):
    mu = 0 
    # sigma_x = 1
    p_i, p_x, p_y, p_z = probs
    sigma_x = newton(root_func, x0=1, args=(p_x,))
    # sigma_y = find_sigma(p_y)
    # sigma_z = find_sigma(p_z)
    """Found on stack overflow to accelerate np.random.choice"""
    if rng is None:
        x = random.random()
    else:
        # x = rng.random()
        x = rng.normal(0, sigma_x)
        # y = rng.normal(0, sigma_y)
        # z = rng.normal(0, sigma_z)

    if abs(x) < th: return options[0], x, sigma_x
    else: return options[1], x, sigma_x
    cum = 0
    for i, p in enumerate(probs):
        cum += p
        if x < cum:
            return options[i]
    return options[-1]


class GaussianPauliErrorModel(BaseErrorModel):

    def __init__(self,
                 r_x: float, r_y: float, r_z: float, 
                 deformation_name: Optional[str] = None, 
                 deformation_kwargs: Optional[dict] = None):

        if not np.isclose(r_x + r_y + r_z, 1):
            raise ValueError(
                f'Noise direction ({r_x}, {r_y}, {r_z}) does not sum to 1.0'
            )
        self._direction = r_x, r_y, r_z
        self._deformation_name = deformation_name

        if deformation_kwargs is not None:
            self._deformation_kwargs = deformation_kwargs
        else:
            self._deformation_kwargs = {}

    @property
    def direction(self) -> Tuple[float, float, float]:
        """Rate of X, Y and Z errors, as given when initializing the
        error model

        Returns
        -------
        (r_x, r_y, r_z): Tuple[float]
            Rate of X, Y and Z errors
        """
        return self._direction

    @property
    def label(self):
        label = 'Pauli X{:.4f}Y{:.4f}Z{:.4f}'.format(*self.direction)
        if self._deformation_name:
            label = 'Deformed ' + self._deformation_name + ' ' + label

        return label

    @property
    def params(self) -> dict:
        """List of class arguments (as a dictionary), that can be saved
        and reused to instantiate the same code"""
        return {
            'r_x': self.direction[0],
            'r_y': self.direction[1],
            'r_z': self.direction[2],
            'deformation_name': self._deformation_name,
            'deformation_kwargs': self._deformation_kwargs
        }

    def generate(self, code: StabilizerCode, error_rate: float, rng=None):
        rng = np.random.default_rng() if rng is None else rng

        p_i, p_x, p_y, p_z = self.probability_distribution(code, error_rate)

        error_pauli = ''
        self.delta = []
        for i in range(code.n):
            pauli_str, delta_x, sigma_x = fast_choice(('I', 'X', 'Y', 'Z'),[p_i[i], p_x[i], p_y[i], p_z[i]],rng=rng)
            error_pauli += pauli_str
            self.delta.append(delta_x)
        
        # self.delta = np.random.normal(0, sigma_x, self.code.n)

        # error_pauli = ''.join([fast_choice(
        #     ('I', 'X', 'Y', 'Z'),
        #     [p_i[i], p_x[i], p_y[i], p_z[i]],
        #     rng=rng
        # ) for i in range(code.n)])


        error = pauli_to_bsf(error_pauli)

        return error

    @functools.lru_cache()
    def probability_distribution(
        self, code: StabilizerCode, error_rate: float
        ) -> Tuple:
        n = code.n
        r_x, r_y, r_z = self.direction

        p: dict = {}
        p['I'] = (1 - error_rate) * np.ones(n)
        p['X'] = (r_x * error_rate) * np.ones(n)
        p['Y'] = (r_y * error_rate) * np.ones(n)
        p['Z'] = (r_z * error_rate) * np.ones(n)

        if self._deformation_name is not None:
            for i in range(code.n):
                deformation = code.get_deformation(
                    code.qubit_coordinates[i], self._deformation_name,
                    **self._deformation_kwargs
                )
                previous_p = {pauli: p[pauli][i] for pauli in ['X', 'Y', 'Z']}
                for pauli in ['X', 'Y', 'Z']:
                    p[pauli][i] = previous_p[deformation[pauli]]

        return p['I'], p['X'], p['Y'], p['Z']