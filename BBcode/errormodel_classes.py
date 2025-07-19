import functools
from typing import Tuple, Optional
import numpy as np
from panqec.codes import StabilizerCode
from panqec.error_models import BaseErrorModel
from panqec.bpauli import pauli_to_bsf
import random

from scipy.integrate import quad 
from scipy.optimize import newton
from scipy.special import erfinv


th = np.sqrt(np.pi)/2

def gaussian(x, sigma):
    return 1/(sigma*np.sqrt(2*np.pi))*np.exp(-x**2/(2*sigma**2)) if sigma != 0.0 else 0.0

def root_func(sigma, p):
    return quad(gaussian, -th, th, args=(sigma,))[0] + p - 1

def fast_choice(options, probs, rng=None):
    """Found on stack overflow to accelerate np.random.choice"""
    if rng is None:
        x = random.random()
    else:
        x = rng.random()
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
        px = p_x[0]
        py = p_y[0]
        pz = p_z[0]
        # sigma_x = newton(root_func, x0=1, args=(p,))
        sigma_x = np.sqrt(np.pi/8) * 1/(erfinv(1-px))
        sigma_y = np.sqrt(np.pi/8) * 1/(erfinv(1-py))
        sigma_z = np.sqrt(np.pi/8) * 1/(erfinv(1-pz))

        error_pauli = ''
        self.gaussian_error_data = []
        options = ('I', 'X', 'Y', 'Z')
        for i in range(code.n):
            delta_x = abs(rng.normal(0, sigma_x))
            delta_y = abs(rng.normal(0, sigma_y))
            delta_z = abs(rng.normal(0, sigma_z))
            if (delta_x > th and delta_y < th and delta_z < th) or (delta_x < th and delta_y > th and delta_z > th):
                error_pauli += 'X'
                delta = delta_x
                sigma = sigma_x
            elif (delta_y > th and delta_x < th and delta_z < th) or (delta_y < th and delta_x > th and delta_z > th):
                error_pauli += 'Y'
                delta = delta_y
                sigma = sigma_y
            elif (delta_z > th and delta_x < th and delta_y < th) or (delta_z < th and delta_x > th and delta_y > th):
                error_pauli += 'Z'
                delta = delta_z
                sigma = sigma_z
            else: 
                error_pauli += 'I'
                delta = np.sqrt(np.pi) - max(delta_x, delta_y, delta_z)
                if delta == np.sqrt(np.pi)-delta_x: sigma = sigma_x
                elif delta == np.sqrt(np.pi)-delta_y: sigma = sigma_y
                else: sigma = sigma_z

            # if delta_x > th and delta_x > delta_y and delta_x > delta_z:
            #     error_pauli += 'X'
            #     delta = delta_x
            #     sigma = sigma_x
            # elif delta_y > th and delta_y > delta_x and delta_y > delta_z:
            #     error_pauli += 'Y'
            #     delta = delta_y
            #     sigma = sigma_y
            # elif delta_z > th and delta_z > delta_x and delta_z > delta_y:
            #     error_pauli += 'Z'
            #     delta = delta_z
            #     sigma = sigma_z
            # else: 
            #     error_pauli += 'I'
            #     delta = np.sqrt(np.pi) - max(delta_x, delta_y, delta_z)
            #     if delta == np.sqrt(np.pi)-delta_x: sigma = sigma_x
            #     elif delta == np.sqrt(np.pi)-delta_y: sigma = sigma_y
            #     else: sigma = sigma_z
                # sigma = min(sigma_x, sigma_y, sigma_z)

            # if delta_x < th: 
            #     error_pauli += options[0]
            #     delta = delta_x
            # else: 
            #     error_pauli += options[1]
            #     delta = np.sqrt(np.pi) - delta_x

            likelihood_wrong = gaussian(np.sqrt(np.pi)-delta, sigma)
            likelihood_correct = gaussian(delta, sigma)
            normalization = likelihood_wrong + likelihood_correct
            # If sigma --> 0, delta --> 0 and hence f_inc = 0 and f_corr = infinity. Then p_correct = 1, and p_wrong = 0
            p_wrong = likelihood_wrong/normalization if normalization != 0.0 else 0.0
            self.gaussian_error_data.append(p_wrong)

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