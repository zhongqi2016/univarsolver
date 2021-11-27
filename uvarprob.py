import sympy as sym

class UniVarProblem:
    """ Univariate problem

    A class for storing information about a univariate optimization problem

    Attributes:
        name: name of a test problem
        sym_objective: objective in sympy formal
        objective: objective as python function
        a: left interval end
        b: right interval end
        min_f: global minumum function's value
        min_x: global minimum point
    """
    def __init__(self, name, objective, a, b, min_f, min_x):
        """ Constructor
        Args:
            name: name of a test example
            objective: objective in sympy format
            a: left end of a feasible interval
            b: right end of a feasible interval
            min_f: global minimum (best known) value
            min_x: global minimum (best known) point
        """
        self.name = name
        self.sym_objective = sym.sympify(objective)
        self.sym_df = self.sym_objective.diff()
        self.sym_ddf = self.sym_df.diff()
        x = sym.symbols('x')
        self.objective = sym.lambdify(x, self.sym_objective)
        self.a = a
        self.b = b
        self.min_f = min_f
        self.min_x = min_x

    def __repr__(self):
        return self.name + ": " + str(self.sym_objective) + " -> min, " + str(self.a) + " <= x <= " + str(self.b) + ", f* = " + str(self.min_f) + ", x* = " + str(self.min_x)