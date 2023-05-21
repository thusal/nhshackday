import pandas as pd
import numpy as np

class STI:
    def __init__(self):
        self.population = None # read from the file
        # add respective columns
        self.initial_infected = 10_000  # randomly assigned

        self.reciprocal_phi = 14  # days
        self.nu = 62.5  # percent
        self.reciprocal_gamma_a = 300  # days
        self.reciprocal_gamma_s = 35  #days
        self.reciprocal_xi = 90  # days
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

        self.C = None # FIXME
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
        self.N = None # total people per age group

    # TODO replace with classes?
    def delta_S(self, age: int, sexual_risk: int):
        return -self.get_lambda(age, sexual_risk) * self.get_S(age, sexual_risk) + (1/self.reciprocal_gamma_a) * (1 - self.h_a) * self.get_Ia(age, sexual_risk) + (1/self.reciprocal_gamma_s) * (1 - self.h_s) * self.get_Is(age, sexual_risk)

    def get_S(self, age: int, sexual_risk: int) -> float:
        pass

    def update_S(self, age: int, sexual_risk: int):
        pass

    def delta_E(self, age: int, sexual_risk: int) -> float:
        return self.get_lambda(age, sexual_risk) * self.get_S(age, sexual_risk) - (1 / self.reciprocal_phi) * self.get_E(age, sexual_risk)
        # TODO check for mu with no indices

    def get_E(self, age: int, sexual_risk: int) -> float:
        pass

    def update_E(self, age: int, sexual_risk: int):
        pass

    def delta_Ia(self, age: int, sexual_risk: int):
        return self.nu * (1 / self.reciprocal_phi) * self.get_E(age, sexual_risk) - self.get_Ia(age, sexual_risk)/self.reciprocal_gamma_s

    def get_Ia(self, age: int, sexual_risk: int) -> float:
        pass

    def update_Ia(self, age: int, sexual_risk: int):
        pass

    def delta_Is(self, age: int, sexual_risk: int):
        return (1 - self.nu) * (1/self.reciprocal_phi) * self.get_E(age, sexual_risk) - self.get_Is(age, sexual_risk) / self.reciprocal_gamma_s

    def get_Is(self, age: int, sexual_risk: int) -> float:
        pass

    def update_Is(self, age: int, sexual_risk: int) -> None:
        pass

    def delta_W(self, age: int, sexual_risk: int) -> float:
        return  (1/self.reciprocal_gamma_a) * self.h_a * self.get_Ia(age, sexual_risk) + (1 / self.reciprocal_gamma_s) * self.h_s * self.get_Is(age, sexual_risk) - (1 / self.reciprocal_xi) * self.get_W(age, sexual_risk)

    def get_W(self, age: int, sexual_risk: int) -> float:
        pass

    def update_W(self, age: int, sexual_risk: int) -> float:
        pass

    def delta_R(self, age: int, sexual_risk: int) -> float:
        return (1 / self.reciprocal_xi) * (self.get_W(age, sexual_risk) + self.get_W2(age, sexual_risk)) - ((1 - self.alpha) * self.get_lambda(age, sexual_risk)) * self.get_R(age, sexual_risk) + (1 / self.reciprocal_gamma_a) * (1 - self.h_a) * self.get_Ia2(age, sexual_risk) + (1 / self.reciprocal_gamma_s) * (1 - self.h_s) * self.get_Is2(age, sexual_risk)

    def get_R(self, age: int, sexual_risk: int) -> float:
        pass

    def update_R(self, age: int, sexual_risk: int):
        pass

    def delta_E2(self, age: int, sexual_risk: int) -> float:
        return (1 - self.alpha) * self.get_lambda(age, sexual_risk) * self.get_R(age, sexual_risk) - (1 / self.reciprocal_phi) * self.get_E2(age, sexual_risk)

    def get_E2(self, age: int, sexual_risk: int) -> float:
        pass

    def update_E2(self, age: int, sexual_risk: int):
        pass

    def delta_Ia2(self, age: int, sexual_risk: int) -> float:
        return self.nu * (1/self.reciprocal_phi) * self.get_E2(age, sexual_risk) - (1 / self.reciprocal_gamma_a) * self.get_Ia2(age, sexual_risk)

    def get_Ia2(self, age: int, sexual_risk: int) -> float:
        pass

    def update_Ia2(self, age: int, sexual_risk: int) -> None:
        pass

    def delta_Is2(self, age: int, sexual_risk: int) -> float:
        return (1 - self.nu) * (1 / self.reciprocal_phi) * self.get_E2(age, sexual_risk) - (1 / self.reciprocal_gamma_s) * self.get_Is2(age, sexual_risk)

    def get_Is2(self, age: int, sexual_risk: int) -> float:
        pass

    def update_Is2(self, age: int, sexual_risk: int) -> None:
        pass

    def delta_W2(self, age: int, sexual_risk: int) -> float:
        return (1 / self.reciprocal_gamma_a) * self.h_a * self.get_Ia2(age, sexual_risk) + (1 / self.reciprocal_gamma_s) * self.h_s * self.get_Is2(age, sexual_risk) - (1 / self.reciprocal_xi) * self.get_W2(age, sexual_risk)

    def get_W2(self, age: int, sexual_risk: int) -> float:
        pass

    def update_W2(self, age: int, sexual_risk: int) -> None:
        pass

    def get_l(self, age:int) -> float:
        return self.l.loc[age, 'spa_rate']

    def get_rho(self, age: int, sexual_risk: int) -> float:
        return self.C * self.get_l(age) * np.power(sexual_risk, self.sigma)

    def get_G(self, x: int, age: int) -> float:
        pass

    def get_N(self, age:int) -> float:
        if self.N is None:
            self.N = self.population.groupby('age').count().to_frame()
        return self.N.loc[age, :].values[0]  # fixme not correct?

    def get_lambda(self,  age: int, sexual_risk: int) -> float:
        pass


