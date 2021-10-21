import sympy as sym

class UniVarProblem:
    """ A class to store univariate optimization problem instance """

    def __init__(self, name, objective, a, b, min_f, mins_x):
        self.name = name
        self.sym_objective = sym.sympify(objective)
        self.sym_df = self.sym_objective.diff()
        self.sym_ddf = self.sym_df.diff()
        x = sym.symbols('x')
        self.objective = sym.lambdify(x, self.sym_objective)
        self.a = a
        self.b = b
        self.min_f = min_f
        self.mins_x = mins_x

    def __repr__(self):
        return self.name + ": " + str(self.sym_objective) + " -> min, " + str(self.a) + " <= x <= " + str(self.b) + ", f* = " + str(self.min_f) + ", " + str(self.mins_x)