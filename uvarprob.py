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
        mins_x: global minimum points
    """
    def __init__(self, name, objective, a, b, min_f, mins_x):
        """ Constructor
        Args:
            name: name of a test example
            objective: objective in sympy format
            a: left end of a feasible interval
            b: right end of a feasible interval
            min_f: global minimum value
            mins_x: global minimum (known) points
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
        self.mins_x = mins_x

    def __repr__(self):
        return self.name + ": " + str(self.sym_objective) + " -> min, " + str(self.a) + " <= x <= " + str(self.b) + ", f* = " + str(self.min_f) + ", " + str(self.mins_x)