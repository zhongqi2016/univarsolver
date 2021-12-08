#Piecewise smooth quadratic estimators



class PSQE:
    """
    Piecewise quadratic estimators
    """

    def __init__(self, a, b, alp, bet, f, df):
        """
        The smooth piecewise quadratic estimator constiructor
        Args:
            a: left interval end
            b: right interval end
            alp: lower end of the Lipschitzian interval
            bet: upper end of the Lipschitzian interval
            f: objective function
            df: objective function derivative
        """
        self.a = a
        self.b = b
        self.alp = alp
        self.bet = bet
        self.fa = f(a)
        self.fb = f(b)
        self.dfa = df(a)
        self.dfb = df(b)

        delt = (self.dfb - self.dfa - alp * (b - a))/(bet - alp)
        print("delt = ", delt)
        self.c = ((delt - a) * self.dfa + (b - delt) * self.dfb + 0.5 * delt**2 * (bet - alp) + alp * delt * (b - a) + 0.5 * alp * (a**2 - b**2) + self.fa - self.fb)/(delt * (bet - alp))
        self.d = self.c + delt

    def __repr__(self):
        return "Estimator " + "a = " + str(self.a) + ", b = " + str(self.b) + ", c = " + str(self.c) + ", d = " + str(self.d) + ", alp = " + str(self.alp) + ", bet = " + str(self.bet) + ", fa = " + str(self.fa) + ", fb = " + str(self.fb) + ", dfa = " + str(self.dfa) + ", dfb = " + str(self.dfb)

    def underestimator(self, x):
        """
        The piecewise quadratic underestimator
        Args:
            x: argument

        Returns: underestimator's value
        """
        if x < self.c:
            return self.fa + self.dfa * (x - self.a) + 0.5 * self.alp * (x - self.a)**2
        elif x < self.d:
            return self.fa + self.dfa * (self.c - self.a) + 0.5 * self.alp * (self.c - self.a)**2 + (self.dfa + self.alp * (self.c - self.a))*(x - self.c) + 0.5 * self.bet * (x - self.c)**2
        else:
            return self.fb + self.dfb * (x - self.b) + 0.5 * self.alp * (x - self.b)**2
