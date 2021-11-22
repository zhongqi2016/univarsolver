"""
Shekel's function generator
"""
import random

import pandas as pd
import random as rnd

class ShekelGen:
    """
    Shekel test functions generator
    Source: Sergeyev, Y. D., Nasso, M. C., Mukhametzhanov, M. S., Kvasov, D. E. (2021).
        Novel local tuning techniques for speeding up one-dimensional algorithms in expensive global optimization using Lipschitz derivatives.
        Journal of Computational and Applied Mathematics, 383, 113134.
    Attributes:
        k_range: the range for k parameter
        a_range: the range for a parameter
        c_range: the range for c parameter
        N: number of terms in the sum
        digits: number of rounded digits in function coefficients
    """

    def __init__(self):
        self.k_range = (1,3)
        self.a_range = (0,10)
        self.c_range = (0.1, 0.3)
        self.N = 10
        self.digits = 3

    def gen_one_problem(self, n, reverse = False):
        """
        Generates Schekel's or reverse Shekel's function
        Args:
            n: the number of a problem
            reverse: if True - generate reverse Shekel's function, if False - normal one

        Returns:
            data frame for a new problem

        """
        name = "shekel_" + str(n)
        formula = ""
        for i in range(0,self.N):
            a = random.uniform(self.a_range[0], self.a_range[1])
            k = random.uniform(self.k_range[0], self.k_range[1])
            c = random.uniform(self.c_range[0], self.c_range[1])
            term = "1./(" + str(round(k*k, self.digits)) + " * (10. * x - " + str(round(a, self.digits)) + ")^2 + " + str(round(c, self.digits)) + ")"
            formula += term if i == 0 else " + " + term

        # print("rand = ", a, k, c)
        dct = dict(name=name, objective=formula, a=0., b=1., min_f=0.11, mins_x=str([1,2]))
        return pd.DataFrame(dct, index = [0])

random.seed(1)
df = None
sk_gen = ShekelGen()
for i in range(0,5):
    dfn = sk_gen.gen_one_problem(i)
    if df is None:
        df = dfn
    else:
        df = pd.concat([df, dfn])

print(df)
df.to_csv('/tmp/shek.csv', index = False)
