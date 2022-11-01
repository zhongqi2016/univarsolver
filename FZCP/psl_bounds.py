"""
The piece-wise linear underestimator from Casado, L. G., MartÍnez, J. A., GarcÍa, I., & Sergeyev, Y. D. (2003). New interval analysis support functions using gradient information in a global minimization algorithm. Journal of Global Optimization, 25(4), 345-362.
"""


class PSL_Bounds:
    """
    Piecewise linear underestimator
    """

    def __init__(self, a, b, alp, bet, f, under):
        """
        The piecewise linear estimator constiructor
        Args:
            a: left interval end
            b: right interval end
            alp: lower end of the Lipschitzian interval
            bet: upper end of the Lipschitzian interval
            f: objective
            under: if True compute the under bound
        """
        self.a = a
        self.b = b
        self.f = f
        self.under = under
        if under:
            self.alp = alp
            self.bet = bet
            self.fa = f(a)
            self.fb = f(b)
        else:
            self.alp = -bet
            self.bet = -alp
            self.fa = -f(a)
            self.fb = -f(b)

        self.c = (self.fa - self.fb + self.bet * self.b - self.alp * self.a) / (self.bet - self.alp)

    def __repr__(self):
        return "Piecewise linear estimator " + "a = " + str(self.a) + ", b = " + str(self.b) + ", c = " + str(
            self.c) + ", alp = " + str(self.alp) + ", bet = " + str(self.bet) \
               + ", fa = " + str(self.fa) + ", fb = " + str(self.fb)

    def estimator(self, x):
        """
        The piecewise linear underestimator
        Args:
            x: argument within [a,b]

        Returns: underestimator's value
        """
        if x <= self.c:
            res = self.fa + self.alp * (x - self.a)
        else:
            res = self.fb + self.bet * (x - self.b)
        if self.under:
            return res
        else:
            return -res

    def nestimator(self, x):
        return -self.estimator(x)

    def lower_bound_and_point(self):
        """
        Returns: Tuple (point where it is achieved, lower bound on interval [a,b])
        """
        record_x = None
        record_v = None
        if self.alp >= 0:
            record_x = self.a
            record_v = self.fa
        elif self.bet <= 0:
            record_x = self.b
            record_v = self.fb
        else:
            record_x = self.c
            record_v = self.estimator(self.c)

        return record_x, record_v

    def record_and_point(self):
        """
        Returns: Tuple (point c where the best value of objective is achieved (a or b), f(c))
        """
        return (self.a, self.fa) if self.fa <= self.fb else (self.b, self.fb)
