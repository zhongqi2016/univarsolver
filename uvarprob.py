import sympy as sym

class UniVarProblem:
    """ A class to store univariate optimization problem instance """

    def __init__(self, name, objective, a, b):
        self.name = name
        self.sym_objective = sym.sympify(objective)
        x = sym.symbols('x')
        self.objective = sym.lambdify(x, self.sym_objective)
        self.a = a
        self.b = b

    def __repr__(self):
        return self.name + ": " + str(self.sym_objective) + " -> min, " + str(self.a) + " <= x <= " + str(self.b)