# Piecewise smooth quadratic estimators


class PSQE_Under:
    """
    Piecewise quadratic underestimator
    """

    def __init__(self, a, b, alp, bet, f, df):
        """
        The smooth piecewise quadratic estimator constiructor
        Args:
            a: left interval end
            b: right interval end
            alp: lower end of the Lipschitzian interval for derivative
            bet: upper end of the Lipschitzian interval for derivative
            f: objective
            df: objective's derivative
        """
        self.a = a
        self.b = b
        self.alp = alp
        self.bet = bet
        self.fa = f(a)
        self.fb = f(b)
        self.dfa = df(a)
        self.dfb = df(b)
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
        record = (None, None)
        ln = len(x_list)
        for i in range(ln - 1):
            x = self.find_argmin(x_list[i], df_list[i], x_list[i + 1], df_list[i + 1])
            if not (x is None):
                check_list.append(x)
        # print(check_list)
        for x in check_list:
            v = self.estimator(x)
            if record[1] is None or v < record[1]:
                record = (x, v)
        return record

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

    def delta_first(self):
        return self.dfa ** 2 - 2 * self.alp * self.fa

    def delta_second(self):
        c = self.fa + self.dfa * (self.c - self.a) + self.alp / 2 * (self.c - self.a) ** 2
        b = self.dfa + self.alp * (self.c - self.a)
        return b ** 2 - 2 * self.bet * c

    def delta_third(self):
        return self.dfb ** 2 - 2 * self.alp * self.fb

    def root_first_left(self, d1):
        return self.a + (-self.dfa - d1 ** 0.5) / self.alp

    def root_second_left(self, d2):
        return self.c + (-self.dfa - self.alp * (self.c - self.a) - d2 ** 0.5) / self.bet

    def root_third_left(self, d3):
        return self.b + (-self.dfb - d3 ** 0.5) / self.alp

    def root_first_right(self, d1):
        return self.a + (-self.dfa + d1 ** 0.5) / self.alp

    def root_second_right(self, d2):
        return self.c + (-self.dfa - self.alp * (self.c - self.a) + d2 ** 0.5) / self.bet

    def root_third_right(self, d3):
        return self.b + (-self.dfb + d3 ** 0.5) / self.alp

    def get_left_end(self):
        # if (flag1)==true, in first section of estimator have root.
        if self.estimator(self.c) <= 0:
            d1 = self.delta_first()
            return self.root_first_left(d1)
        else:
            d1 = self.delta_first()
            if self.under_est_der_le_0(self.a, self.c) and d1 >= 0:
                return self.root_first_left(d1)

        if self.estimator(self.d) <= 0:
            d2 = self.delta_second()
            return self.root_second_left(d2)
        else:
            d2 = self.delta_second()
            if self.under_est_der_le_0(self.c, self.d) and d2 >= 0:
                return self.root_second_left(d2)

        if self.estimator(self.b) <= 0:
            d3 = self.delta_third()
            return self.root_third_left(d3)
        else:
            d3 = self.delta_third()
            if self.under_est_der_le_0(self.d, self.b) and d3 >= 0:
                return self.root_third_left(d3)

        # d3 = self.delta_third()
        # if d3 >= 0:
        #     return self.root_third_left(d3)

        return None

    def get_right_end2(self):
        # if (flag1)==true, in first section of estimator have root.

        if self.estimator(self.d) <= 0:
            d1 = self.delta_third()
            return self.root_third_right(d1)
        else:
            d1 = self.delta_third()
            if self.under_est_der_le_0(self.d, self.b) and d1 >= 0:
                return self.root_third_right(d1)
        if self.estimator(self.c) <= 0:
            d2 = self.delta_second()
            return self.root_second_right(d2)
        else:
            d2 = self.delta_second()
            if self.under_est_der_le_0(self.c, self.d) and d2 >= 0:
                return self.root_second_right(d2)
        if self.estimator(self.a) <= 0:
            d3 = self.delta_first()
            return self.root_first_right(d3)
        else:
            d3 = self.delta_third()
            if self.under_est_der_le_0(self.a, self.c) and d3 >= 0:
                return self.root_first_right(d3)

        # d3 = self.delta_third()
        # if d3 >= 0:
        #     return self.root_third_left(d3)
        return None
