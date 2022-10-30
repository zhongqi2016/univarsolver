# Piecewise smooth quadratic estimators


class PSQE_Bounds:
    """
    Piecewise quadratic underestimator
    """

    def __init__(self, a, b, alp, bet, f, df, under):
        """
        The smooth piecewise quadratic estimator constiructor
        Args:
            a: left interval end
            b: right interval end
            alp: lower end of the Lipschitzian interval for derivative
            bet: upper end of the Lipschitzian interval for derivative
            f: objective
            df: objective's derivative
            under: if True compute the under bound
        """
        self.a = a
        self.b = b
        self.under = under
        if under:
            self.alp = alp
            self.bet = bet
            self.fa = f(a)
            self.fb = f(b)
            self.dfa = df(a)
            self.dfb = df(b)
        else:
            self.alp = -bet
            self.bet = -alp
            self.fa = -f(a)
            self.fb = -f(b)
            self.dfa = -df(a)
            self.dfb = -df(b)
        self.f = f
        self.df = df

        delt = (self.dfb - self.dfa - alp * (b - a)) / (bet - alp)
        # print("delt = ", delt)
        self.c = ((delt - a) * self.dfa + (b - delt) * self.dfb + 0.5 * delt ** 2 * (bet - alp) + alp * delt * (
                b - a) + 0.5 * alp * (a ** 2 - b ** 2) + self.fa - self.fb) / (delt * (bet - alp))
        self.d = self.c + delt

    def __repr__(self):
        return "Estimator " + "a = " + str(self.a) + ", b = " + str(self.b) + ", c = " + str(self.c) + ", d = " + str(
            self.d) + ", alp = " + str(self.alp) + ", bet = " + str(self.bet) + ", fa = " + str(
            self.fa) + ", fb = " + str(self.fb) + ", dfa = " + str(self.dfa) + ", dfb = " + str(self.dfb)

    def estimator(self, x):
        """
        The piecewise quadratic underestimator
        Args:
            x: argument

        Returns: underestimator's value
        """
        if x < self.c:
            return self.fa + self.dfa * (x - self.a) + 0.5 * self.alp * (x - self.a) ** 2
        elif x < self.d:
            return self.fa + self.dfa * (self.c - self.a) + 0.5 * self.alp * (self.c - self.a) ** 2 + (
                    self.dfa + self.alp * (self.c - self.a)) * (x - self.c) + 0.5 * self.bet * (x - self.c) ** 2
        else:
            return self.fb + self.dfb * (x - self.b) + 0.5 * self.alp * (x - self.b) ** 2

    def nestimator(self, x):
        return -self.estimator(x)

    def estimators_derivative(self, x):
        """
        The piecewise linear underestimator's derivative
        Args:
            x: argument

        Returns: underestimator's derivative value
        """
        if x < self.c:
            return self.dfa + self.alp * (x - self.a)
        elif x < self.d:
            return self.dfa + self.alp * (self.c - self.a) + self.bet * (x - self.c)
        else:
            return self.dfb + self.alp * (x - self.b)

    def lower_bound_and_point(self):
        """
        Returns: Tuple (point where it is achieved, lower bound on interval [a,b])
        """
        x_list = [self.a, self.c, self.d, self.b]
        df_list = [self.estimators_derivative(x) for x in x_list]
        check_list = [self.a, self.b]
        record_x = None
        record_v = None
        ln = len(x_list)
        for i in range(ln - 1):
            x = self.find_argmin(x_list[i], df_list[i], x_list[i + 1], df_list[i + 1])
            if not (x is None):
                check_list.append(x)
        # print(check_list)
        for x in check_list:
            v = self.estimator(x)
            if record_v is None or v < record_v:
                record_v = v
                record_x = x
        if not self.under:
            record_v = -record_v
        return record_x, record_v

    def record_and_point(self):
        """
        Returns: Tuple (point c where the best value of objective is achieved (a or b), f(c))
        """
        return (self.a, self.fa) if self.fa <= self.fb else (self.b, self.fb)

    def find_argmin(self, x1, df1, x2, df2):
        if df1 == 0 and df2 == 0:
            xs = 0.5 * (x1 + x2)
        elif df1 <= 0 <= df2:
            xs = x1 + (-df1) * (x2 - x1) / (df2 - df1)
        else:
            xs = None
        return xs
