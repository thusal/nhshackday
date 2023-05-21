import pandas as pd
import numpy as np


class STI:
    def __init__(self, base_data):
        self.population = base_data  # read from the file

        self.reciprocal_phi = 14  # days
        self.phi = 1 / self.reciprocal_phi

        self.nu = 62.5  # percent

        self.reciprocal_gamma_a = 300  # days
        self.gamma_a = 1 / self.reciprocal_gamma_a

        self.reciprocal_gamma_s = 35  # days
        self.gamma_s = 1 / self.reciprocal_gamma_s

        self.reciprocal_xi = 90  # days
        self.xi = 1 / self.reciprocal_xi

        self.sigma = 3  # power law exponent
        self.e_g = 0.7  # age group mixing
        self.e_h = 0.3  # sex group mixing
        self.p = 0.0375  # transmission probability
        self.tau = 180  # days in partnership
        self.h_a = 1  # fraction of asymptomatic with protective immunity
        self.h_s = 0  # fraction of symptomatic with protective immunity
        self.m_ia = 1  # coital act frequency asymptomatic
        self.m_is = 0.645  # coital act frequency symptomatic

        self.alpha = 0.5  # fixme find the value

        self.C = None  # FIXME
        self.d = {'age': [16, 17, 18, 19,
                          20, 21, 22, 23, 24,
                          25, 26, 27, 28, 29,
                          30, 31, 32, 33, 34,
                          35, 36, 37, 38, 39,
                          40, 41, 42, 43, 44,
                          45, 46, 47, 48, 49,
                          50, 51, 52, 53, 54,
                          55, 56, 57, 58, 59,
                          60, 61, 62, 63, 64,
                          65, 66, 67, 68, 69,
                          70, 71, 72, 73, 74,
                          75, 76, 77, 78, 79,
                          80, 81, 82, 83, 84,
                          85, 86, 87, 88, 89,
                          90, 91, 92, 93, 94,
                          95, 96, 97, 98],
                  'spa_rate': [0.85, 0.85, 0.85, 0.85,
                               1.2, 1.2, 1.2, 1.2, 1.2,
                               0.61, 0.61, 0.61, 0.61, 0.61,
                               0.33, 0.33, 0.33, 0.33, 0.33,
                               0.25, 0.25, 0.25, 0.25, 0.25,
                               0.19, 0.19, 0.19, 0.19, 0.19,
                               0.14, 0.14, 0.14, 0.14, 0.14,
                               0.095, 0.095, 0.095, 0.095, 0.095,
                               0.065, 0.065, 0.065, 0.065, 0.065,
                               0.045, 0.045, 0.045, 0.045, 0.045,
                               0.033, 0.033, 0.033, 0.033, 0.033,
                               0.025, 0.025, 0.025, 0.025, 0.025]}
        self.l = pd.DataFrame(data=self.d).set_index('age').fillna(0)
        self.N = self.population.copy()

        self.population['S'] = self.population['count'] - 100

        self.population['E'] = 0

        self.population['Ia'] = 50
        self.population['Is'] = 50

        self.population['W'] = 0
        self.population['R'] = 0
        self.population['E2'] = 0
        self.population['Ia2'] = 0
        self.population['Is2'] = 0
        self.population['W2'] = 0

        _n = {'age': [16, 17, 18, 19,
                      20, 21, 22, 23, 24,
                      25, 26, 27, 28, 29,
                      30, 31, 32, 33, 34,
                      35, 36, 37, 38, 39,
                      40, 41, 42, 43, 44,
                      45, 46, 47, 48, 49,
                      50, 51, 52, 53, 54,
                      55, 56, 57, 58, 59,
                      60, 61, 62, 63, 64,
                      65, 66, 67, 68, 69,
                      70, 71, 72, 73, 74,
                      75, 76, 77, 78, 79,
                      80, 81, 82, 83, 84,
                      85, 86, 87, 88, 89,
                      90, 91, 92, 93, 94,
                      95, 96, 97, 98],
              'coital_frequency': [3.9, 3.9, 3.9, 3.9,
                                   3.5, 3.5, 3.5, 3.5, 3.5,
                                   3.2, 3.2, 3.2, 3.2, 3.2,
                                   2.9, 2.9, 2.9, 2.9, 2.9,
                                   2.5, 2.5, 2.5, 2.5, 2.5,
                                   2.2, 2.2, 2.2, 2.2, 2.2,
                                   1.9, 1.9, 1.9, 1.9, 1.9,
                                   1.5, 1.5, 1.5, 1.5, 1.5,
                                   1.2, 1.2, 1.2, 1.2, 1.2,
                                   0.87, 0.87, 0.87, 0.87, 0.87,
                                   0.53, 0.53, 0.53, 0.53, 0.53,
                                   0.2, 0.2, 0.2, 0.2, 0.2]}
        self.n = pd.DataFrame(data=_n).set_index('age').fillna(0)

    # TODO replace with classes?
    def delta_S(self, age: int, sexual_risk: int):
        return -self.get_lambda(age, sexual_risk) * self.get_S(age,
                                                               sexual_risk) + self.gamma_a * (
                1 - self.h_a) * self.get_Ia(age,
                                            sexual_risk) + self.gamma_s * (
                1 - self.h_s) * self.get_Is(age, sexual_risk)

    def get_S(self, age: int, sexual_risk: int) -> float:
        return self.population.loc[age, ['S']].values[0]

    def update_S(self, age: int, sexual_risk: int):
        self.population.loc[age, ['S']] += int(self.delta_S(age, sexual_risk))

    def delta_E(self, age: int, sexual_risk: int) -> float:
        return self.get_lambda(age, sexual_risk) * self.get_S(age,
                                                              sexual_risk) - self.phi * self.get_E(
            age, sexual_risk)
        # TODO check for mu with no indices

    def get_E(self, age: int, sexual_risk: int) -> float:
        return self.population.loc[age, ['E']].values[0]

    def update_E(self, age: int, sexual_risk: int):
        self.population.loc[age, ['E']] += int(self.delta_E(age, sexual_risk))

    def delta_Ia(self, age: int, sexual_risk: int):
        return self.nu * self.phi * self.get_E(age,
                                               sexual_risk) - self.gamma_s * self.get_Ia(
            age, sexual_risk)

    def get_Ia(self, age: int, sexual_risk: int) -> float:
        return self.population.loc[age, ['Ia']].values[0]

    def update_Ia(self, age: int, sexual_risk: int):
        self.population.loc[age, ['Ia']] += int(self.delta_Ia(age, sexual_risk))

    def delta_Is(self, age: int, sexual_risk: int):
        return (1 - self.nu) * self.phi * self.get_E(age,
                                                     sexual_risk) - self.gamma_s * self.get_Is(
            age, sexual_risk)

    def get_Is(self, age: int, sexual_risk: int) -> float:
        return self.population.loc[age, ['Is']].values[0]

    def update_Is(self, age: int, sexual_risk: int) -> None:
        self.population.loc[age, ['Is']] += int(self.delta_Is(age, sexual_risk))

    def delta_W(self, age: int, sexual_risk: int) -> float:
        return self.gamma_a * self.h_a * self.get_Ia(age,
                                                     sexual_risk) + self.gamma_s * self.h_s * self.get_Is(
            age, sexual_risk) - self.xi * self.get_W(age, sexual_risk)

    def get_W(self, age: int, sexual_risk: int) -> float:
        return self.population.loc[age, ['W']].values[0]

    def update_W(self, age: int, sexual_risk: int):
        self.population.loc[age, ['W']] += int(self.delta_W(age, sexual_risk))

    def delta_R(self, age: int, sexual_risk: int) -> float:
        return self.xi * (self.get_W(age, sexual_risk) + self.get_W2(age,
                                                                     sexual_risk)) - (
                (1 - self.alpha) * self.get_lambda(age,
                                                   sexual_risk)) * self.get_R(
            age, sexual_risk) + self.gamma_a * (1 - self.h_a) * self.get_Ia2(
            age, sexual_risk) + self.gamma_s * (1 - self.h_s) * self.get_Is2(
            age, sexual_risk)

    def get_R(self, age: int, sexual_risk: int) -> float:
        return self.population.loc[age, ['R']].values[0]

    def update_R(self, age: int, sexual_risk: int):
        self.population.loc[age, ['R']] += int(self.delta_R(age, sexual_risk))

    def delta_E2(self, age: int, sexual_risk: int) -> float:
        return (1 - self.alpha) * self.get_lambda(age,
                                                  sexual_risk) * self.get_R(age,
                                                                            sexual_risk) - self.phi * self.get_E2(
            age, sexual_risk)

    def get_E2(self, age: int, sexual_risk: int) -> float:
        return self.population.loc[age, ['E2']].values[0]

    def update_E2(self, age: int, sexual_risk: int):
        self.population.loc[age, ['E2']] += int(self.delta_E2(age, sexual_risk))

    def delta_Ia2(self, age: int, sexual_risk: int) -> float:
        return self.nu * self.phi * self.get_E2(age,
                                                sexual_risk) - self.gamma_a * self.get_Ia2(
            age, sexual_risk)

    def get_Ia2(self, age: int, sexual_risk: int) -> float:
        return self.population.loc[age, ['Ia2']].values[0]

    def update_Ia2(self, age: int, sexual_risk: int) -> None:
        self.population.loc[age, ['Ia2']] += int(self.delta_Ia2(age, sexual_risk))

    def delta_Is2(self, age: int, sexual_risk: int) -> float:
        return (1 - self.nu) * self.phi * self.get_E2(age,
                                                      sexual_risk) - self.gamma_s * self.get_Is2(
            age, sexual_risk)

    def get_Is2(self, age: int, sexual_risk: int) -> float:
        return self.population.loc[age, ['Is2']].values[0]

    def update_Is2(self, age: int, sexual_risk: int) -> None:
        self.population.loc[age, ['Is2']] += int(self.delta_Is2(age, sexual_risk))

    def delta_W2(self, age: int, sexual_risk: int) -> float:
        return self.gamma_a * self.h_a * self.get_Ia2(age,
                                                      sexual_risk) + self.gamma_s * self.h_s * self.get_Is2(
            age, sexual_risk) - self.xi * self.get_W2(age, sexual_risk)

    def get_W2(self, age: int, sexual_risk: int) -> float:
        return self.population.loc[age, ['W2']].values[0]

    def update_W2(self, age: int, sexual_risk: int) -> None:
        self.population.loc[age, ['W2']] += int(self.delta_W2(age, sexual_risk))

    def get_l(self, age: int) -> float:
        return self.l.loc[age, 'spa_rate']

    def get_rho(self, age: int, sexual_risk: int) -> float:
        return self.C * self.get_l(age) * np.power(sexual_risk, self.sigma)

    def get_G(self, x: int, j: int) -> float:
        k_ = 1  # todo only one sex risk group atm
        t = self.get_rho(j, k_) * self.get_N(j, k_)
        t1 = np.sum(
            [self.get_rho(j, k_) * self.get_N(j, k_) for j in self.d['age']])
        return self.e_g * np.identity(len(self.d['age']))[
            self.d['age'].index(x), self.d['age'].index(j)] + (
                    1 - self.e_g) * t / t1

    def get_H(self, y: int, k: int) -> float:
        return 1

    def get_N(self, age: int, sexual_risk: int) -> float:
        return self.N.loc[age, :].values[0]

    def get_lambda(self, age: int, sexual_risk: int) -> float:
        # todo check s7 formula, there is some odd division
        t = [self.get_G(age, j) * ((self.get_rho(j, 1) * self.get_qia(age,
                                                                      sexual_risk) * (
                                                self.get_Ia(j,
                                                            1) + self.get_Ia2(j,
                                                                              1)) + self.get_rho(
            j, 1) * self.get_qis(age, sexual_risk) * (self.get_Is(j,
                                                                  1) + self.get_Is2(
            j, 1))) / (self.get_rho(j, 1) * self.get_N(j, 1))) for j in
             self.d['age']]
        return self.get_rho(age, sexual_risk) * np.sum(t)

    def get_n(self, age: int, sexual_risk: int):
        return self.n.loc[age]

    def get_qia(self, age: int, sexual_risk: int):
        return 1 - np.power(1 - self.p,
                            self.m_ia * self.get_n(age, sexual_risk) * self.tau)

    def get_qis(self, age: int, sexual_risk: int):
        return 1 - np.power(1 - self.p,
                            self.m_is * self.get_n(age, sexual_risk) * self.tau)

