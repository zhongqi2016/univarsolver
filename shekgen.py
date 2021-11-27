"""
Shekel's function generator
"""
import random
import pandas as pd
import gridsearch as gs
import uvarprob as uvpr


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
        gs_step: the grid step size
    """

    def __init__(self):
        self.k_range = (1,3)
        self.a_range = (0,10)
        self.c_range = (0.1, 0.3)
        self.N = 10
        self.digits = 5
        self.gs_step = 1e-5

    def gen_one_problem(self, n, reverse = False):
        """
        Generates Schekel's or reverse Shekel's function
        Args:
            n: the number of a problem
            reverse: if True - generate reverse Shekel's function, if False - normal one

        Returns:
            data frame for a new problem

        """
        name = ("rshekel_" if reverse else "shekel_") + str(n)
        print("Generating problem ", name, " ... ", end='')
        formula = ""
        for i in range(0,self.N):
            a = random.uniform(self.a_range[0], self.a_range[1])
            k = random.uniform(self.k_range[0], self.k_range[1])
            c = random.uniform(self.c_range[0], self.c_range[1])
            term = ("1./(" if reverse else "-1/(") + str(round(k*k, self.digits)) + " * (10. * x - " + str(round(a, self.digits)) + ")^2 + " + str(round(c, self.digits)) + ")"
            formula += term if i == 0 else " + " + term


        # print("rand = ", a, k, c)

        prob = uvpr.UniVarProblem("name", formula, 0., 1., 0.0, 0.0)

        true_min = gs.grid_search(prob.objective, prob.a, prob.b, self.gs_step)
         # print("true_min = ", round(true_min[0], self.digits), round(true_min[1], self.digits))
        dct = dict(name=name, objective=formula, a=0., b=1., min_f=round(true_min[1], self.digits), min_x=round(true_min[0], self.digits))
        print(" done.")
        return pd.DataFrame(dct, index = [0])

# Tuning:

# Seeding
random.seed(1)

# Number of problems
np = 10

# Reverse Shekel (True) or normal Shekel (False) functions
rshek = False

# File name to write the results
fname = "/tmp/shek.csv"

df = None
sk_gen = ShekelGen()
for i in range(0, np):
    dfn = sk_gen.gen_one_problem(i, rshek)
    if df is None:
        df = dfn
    else:
        df = pd.concat([df, dfn])

print(df)
df.to_csv(fname, index = False)
