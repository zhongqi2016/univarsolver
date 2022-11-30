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

    def find_min_under_zero(self, x1, df1, x2, df2):
        if df1 > 0 and df2 > 0:
            if self.estimator(x1) <= 0:
                return True
            else:
                return False
        elif df1 < 0 and df2 < 0:
            if self.estimator(x2) <= 0:
                return True
            else:
                return False
        elif df1 < 0:
            if self.estimator(x1 + (-df1) * (x2 - x1) / (df2 - df1)) <= 0:
                return True
            else:
                return False
        else:
            if self.estimator(x1) <= 0:
                return True
            else:
                if self.estimator(x2) <= 0:
                    return True
                else:
                    return False

    def find_max_above_zero(self, x1, df1, x2, df2):
        if df1 > 0 and df2 > 0:
            if self.estimator(x2) >= 0:
                return True
            else:
                return False
        elif df1 < 0 and df2 < 0:
            if self.estimator(x1) >= 0:
                return True
            else:
                return False
        elif df1 > 0:
            max_x = x1 + (-df1) * (x2 - x1) / (df2 - df1)
            if self.estimator(max_x) >= 0:
                return True
            else:
                return False
        else:
            if self.estimator(x1) >= 0:
                return True
            else:
                if self.estimator(x2) >= 0:
                    return True
                else:
                    return False

    def right_root_first(self):
        # print("first right")
        if self.alp < 0:
            return self.a + (-self.dfa - (self.dfa ** 2 - 2 * self.alp * self.fa) ** 0.5) / self.alp
        else:
            return self.a + (-self.dfa + (self.dfa ** 2 - 2 * self.alp * self.fa) ** 0.5) / self.alp

    def left_root_first(self):
        # print("first left")
        if self.alp < 0:
            return self.a + (-self.dfa + (self.dfa ** 2 - 2 * self.alp * self.fa) ** 0.5) / self.alp
        else:
            return self.a + (-self.dfa - (self.dfa ** 2 - 2 * self.alp * self.fa) ** 0.5) / self.alp

    def right_root_third(self):
        # print("third right")
        if self.alp < 0:
            return self.b + (-self.dfb - (self.dfb ** 2 - 2 * self.alp * self.fb) ** 0.5) / self.alp
        else:
            return self.b + (-self.dfb + (self.dfb ** 2 - 2 * self.alp * self.fb) ** 0.5) / self.alp

    def left_root_third(self):
        # print("third left")
        if self.alp < 0:
            return self.b + (-self.dfb + (self.dfb ** 2 - 2 * self.alp * self.fb) ** 0.5) / self.alp
        else:
            return self.b + (-self.dfb - (self.dfb ** 2 - 2 * self.alp * self.fb) ** 0.5) / self.alp

    def right_root_second(self):
        # print("second right")
        c = self.fa + self.dfa * (self.c - self.a) + self.alp / 2 * (self.c - self.a) ** 2
        b = self.dfa + self.alp * (self.c - self.a)
        if self.bet>0:
            res = self.c + (-b + (b ** 2 - 2 * self.bet * c) ** 0.5) / self.bet
        else:
            res = self.c + (-b - (b ** 2 - 2 * self.bet * c) ** 0.5) / self.bet
        return res

    def left_root_second(self):
        # print("second left")
        c = self.fa + self.dfa * (self.c - self.a) + self.alp / 2 * (self.c - self.a) ** 2
        b = self.dfa + self.alp * (self.c - self.a)
        if self.bet > 0:
            res = self.c + (-b - (b ** 2 - 2 * self.bet * c) ** 0.5) / self.bet
        else:
            res = self.c + (-b + (b ** 2 - 2 * self.bet * c) ** 0.5) / self.bet
        return res

    def getNewTrialPoint(self):
        if self.under:
            return self.get_trial_point_under()
        else:
            return self.get_trial_point_upper()

    def get_trial_point_under(self):
        if self.estimator(self.c) <= 0:
            new_point = self.right_root_first()
        else:
            est_der_c = self.estimators_derivative(self.c)
            est_der_d = self.estimators_derivative(self.d)
            if (est_der_c > 0 and est_der_d < 0) or (est_der_c < 0 and est_der_d > 0):
                if self.find_min_under_zero(self.c, est_der_c, self.d, est_der_d):
                    new_point = self.left_root_second()
                else:
                    new_point = self.right_root_third()
            else:
                if self.estimator(self.d) > 0:
                    new_point = self.right_root_third()
                else:
                    new_point = self.left_root_second()
        return new_point

    def get_trial_point_upper(self):
        est_der_a = self.estimators_derivative(self.a)
        est_der_c = self.estimators_derivative(self.c)
        if self.find_max_above_zero(self.a, est_der_a, self.c, est_der_c):
            new_point = self.left_root_first()
        else:
            if self.estimator(self.d) >= 0:
                new_point = self.right_root_second()
            else:
                est_der_d = self.estimators_derivative(self.d)
                est_der_b = self.estimators_derivative(self.b)
                if self.find_max_above_zero(self.d, est_der_d, self.b, est_der_b):
                    new_point = self.left_root_third()
                else:
                    new_point = self.b
        return new_point
