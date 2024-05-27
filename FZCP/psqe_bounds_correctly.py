# Piecewise smooth quadratic estimators
import interval_arithmetics as int_arith
import decimal as dec


class PSQE_Bounds:
    """
    Piecewise quadratic underestimator
    """

    def __init__(self, a: dec.Decimal, b: dec.Decimal, alp: dec.Decimal, bet: dec.Decimal, ival_fa: int_arith.Interval,
                 ival_fb: int_arith.Interval, ival_dfa: int_arith.Interval, ival_dfb: int_arith.Interval, under: bool):
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
            self.ival_fa = ival_fa
            self.ival_fb = ival_fb
            self.ival_dfa = ival_dfa
            self.ival_dfb = ival_dfb
            self.fa = ival_fa.a
            self.fb = ival_fb.a
            self.dfa = ival_dfa.a
            self.dfb = ival_dfb.a
        else:
            self.alp = -bet
            self.bet = -alp
            self.ival_fa = -ival_fa
            self.ival_fb = -ival_fb
            self.ival_dfa = -ival_dfa
            self.ival_dfb = -ival_dfb
            self.fa = -ival_fa.a
            self.fb = -ival_fb.a
            self.dfa = -ival_dfa.a
            self.dfb = -ival_dfb.b

        self.ival_a = int_arith.Interval(self.a, self.a)
        self.ival_b = int_arith.Interval(self.b, self.b)

        self.ival_alp = int_arith.Interval(self.alp, self.alp)
        self.ival_bet = int_arith.Interval(self.bet, self.bet)

        delt = (self.ival_dfb - self.ival_dfa - self.ival_alp * (self.ival_b - self.ival_a)) / (
                self.ival_bet - self.ival_alp)
        # print("delt = ", delt)
        # self.c = ((delt - a) * self.dfa + (b - delt) * self.dfb + 0.5 * delt ** 2 * (
        #         self.bet - self.alp) + self.alp * delt * (
        #                   b - a) + 0.5 * self.alp * (a ** 2 - b ** 2) + self.fa - self.fb) / (
        #                  delt * (self.bet - self.alp))
        self.ival_c = (self.ival_alp * (self.ival_b - self.ival_a) + self.ival_dfa - self.ival_dfb) / (
                2 * (self.ival_bet - self.ival_alp)) + (
                              self.ival_fa - self.ival_fb + self.ival_b * self.ival_dfb - self.ival_a * self.ival_dfa + self.ival_alp * (
                              self.ival_a ** 2 - self.ival_b ** 2) / 2) / (
                              self.ival_dfb - self.ival_dfa - self.ival_alp * (self.ival_b - self.ival_a))
        self.ival_d = self.ival_c + delt
        # if self.c > self.b:
        #     print('error c>b')
        #     print('a=', self.a, 'b=', self.b, 'c=', self.c, 'delt=', delt, 'fa=', self.fa, 'fb=', self.fb, 'alp=',
        #           self.alp, 'beta=', self.bet)
        if self.ival_d.a > self.b:
            print('error d>b')
            print('a=', self.a, 'b=', self.b, 'c=', self.ival_c, 'd=', self.ival_d, 'delt=', delt, 'fa=', self.fa,
                  'fb=', self.fb, 'alp=',
                  self.alp, 'beta=', self.bet)

    def __repr__(self):
        return "Estimator " + "a = " + str(self.a) + ", b = " + str(self.b) + ", c = " + str(
            self.ival_c) + ", d = " + str(
            self.ival_d) + ", alp = " + str(self.alp) + ", bet = " + str(self.bet) + ", fa = " + str(
            self.fa) + ", fb = " + str(self.fb) + ", dfa = " + str(self.dfa) + ", dfb = " + str(self.dfb)

    def normal_cd(self) -> bool:
        if self.ival_c.a > self.a and self.ival_d.b < self.b:
            return True
        else:
            return False

    def estimator_q1(self, ival_x: int_arith.Interval):
        """ first part of quadratic estimator """
        return self.ival_fa + self.ival_dfa * (ival_x - self.ival_a) + 0.5 * self.ival_alp * (ival_x - self.ival_a) ** 2

    def estimator_q2(self, ival_x: int_arith.Interval):
        """ second part of quadratic estimator """
        return self.ival_fa + self.ival_dfa * (self.ival_c - self.ival_a) + 0.5 * self.ival_alp * (
                self.ival_c - self.ival_a) ** 2 + (
                self.ival_dfa + self.ival_alp * (self.ival_c - self.ival_a)) * (
                ival_x - self.ival_c) + 0.5 * self.ival_bet * (ival_x - self.ival_c) ** 2

    def estimator_q3(self, ival_x: int_arith.Interval):
        """ third part of quadratic estimator """
        return self.ival_fb + self.ival_dfb * (ival_x - self.ival_b) + 0.5 * self.ival_alp * (ival_x - self.ival_b) ** 2

    def estimator(self, x: dec.Decimal):
        """
        The piecewise quadratic underestimator
        Args:
            x: argument

        Returns: underestimator's value
        """
        ival_x = int_arith.Interval(x, x)
        if x <= self.ival_c.a:
            res = self.estimator_q1(ival_x)
        elif x < self.ival_d:
            res = self.estimator_q2(ival_x)
        else:
            res = self.estimator_q3(ival_x)

        # if self.under:
        return res
        # else:
        #     return -res

    def get_fb(self):
        return self.fb

    def nestimator(self, x):
        return -self.estimator(x)

    def estd_q1(self, ival_x: int_arith.Interval):
        return self.ival_dfa + self.ival_alp * (ival_x - self.ival_a)

    def estd_q2(self, ival_x: int_arith.Interval):
        return self.ival_dfa + self.ival_alp * (self.ival_c - self.ival_a) + self.ival_bet * (ival_x - self.ival_c)

    def estd_q3(self, ival_x: int_arith.Interval):
        return self.ival_dfb + self.ival_alp * (ival_x - self.ival_b)

    def estimators_derivative(self, ival_x: int_arith.Interval):
        """
        The piecewise linear underestimator's derivative
        Args:
            x: argument

        Returns: underestimator's derivative value
        """
        if ival_x.b < self.ival_c.a:
            res = self.estd_q1(ival_x)
        elif ival_x.b < self.ival_d.a:
            res = self.estd_q2(ival_x)
        else:
            res = self.estd_q3(ival_x)
        # if self.under:
        return res
        # else:
        #     return -res

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
        for x in check_list:
            v = self.estimator(x)
            if record_v is None or v < record_v:
                record_v = v
                record_x = x
        if not self.under:
            record_v = -record_v
        return record_x, record_v

    def lower_bound_and_point2(self):
        if self.fb <= 0:
            return 0
        x_list = [self.a, self.c, self.d, self.b]
        df_list = [self.estimators_derivative(x) for x in x_list]
        ln = len(x_list)

        for i in range(ln - 1):
            x = self.find_argmin(x_list[i], df_list[i], x_list[i + 1], df_list[i + 1])
            if (x is not None) and self.estimator(x) <= 0:
                return i + 1

        return -1

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
        """
            Return True if the maximum value of f(x)>0 (x in [x1,x2]), otherwise return False
        """
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

    def under_est_der_le_0(self, num1: int_arith.Interval, num2: int_arith.Interval):
        """
        if(\phi'(num1) * \phi'(num2)<0) return true, else return false.
        """
        est_der_num1 = self.estimators_derivative(num1)
        est_der_num2 = self.estimators_derivative(num2)
        if est_der_num1.a < 0 and est_der_num2.b > 0:
            return True
        else:
            return False

    def upper_est_der_le_0(self, num1, num2):
        """
        if(\phi'(num1) * \phi'(num2)<0) return true, else return false.
        """
        est_der_num1 = self.estimators_derivative(num1)
        est_der_num2 = self.estimators_derivative(num2)
        if est_der_num1 > 0 and est_der_num2 < 0:
            return True
        else:
            return False

    def delta_first(self) -> int_arith.Interval:
        return self.ival_dfa ** 2 - 2 * self.ival_alp * self.ival_fa

    def delta_second(self) -> int_arith.Interval:
        c = self.ival_fa + self.ival_dfa * (self.ival_c - self.ival_a) + self.ival_alp / 2 * (
                self.ival_c - self.ival_a) ** 2
        b = self.ival_dfa + self.ival_alp * (self.ival_c - self.ival_a)
        return b ** 2 - 2 * self.ival_bet * c

    def delta_third(self) -> int_arith.Interval:
        return self.ival_dfb ** 2 - 2 * self.ival_alp * self.ival_fb

    def root_first_left(self, d1):
        return self.ival_a + (-self.ival_dfa - d1.sqrt()) / self.ival_alp

    def root_second_left(self, d2):
        return self.ival_c + (-self.ival_dfa - self.ival_alp * (self.ival_c - self.ival_a) - d2.sqrt()) / self.ival_bet

    def root_third_left(self, d3):
        return self.ival_b + (-self.ival_dfb - d3.sqrt()) / self.ival_alp

    def root_first_right(self, d1):
        return self.ival_a + (-self.ival_dfa + d1.sqrt()) / self.ival_alp

    def root_second_right(self, d2):
        return self.ival_c + (-self.ival_dfa - self.ival_alp * (self.ival_c - self.ival_a) + d2.sqrt()) / self.ival_bet

    def root_third_right(self, d3):
        return self.ival_b + (-self.ival_dfb + d3.sqrt()) / self.ival_alp

    def get_left_end(self) -> dec.Decimal | None:
        """get the first root of under bound"""
        if self.ival_c.a < self.a:
            print('warning: c*<a')
            return self.a
        est1 = self.estimator_q1(self.ival_c)

        if est1.a <= 0:
            d1 = self.delta_first()
            if d1.a < 0: d1.a = int_arith.c_zero
            res = self.root_first_left(d1).a
            if res < self.a:
                print('warning: left end q1 < self.a')
                res = self.a
            return res
        else:
            d1 = self.delta_first()
            if self.under_est_der_le_0(self.ival_a, self.ival_c) and d1.a >= 0:
                res = self.root_first_left(d1).a
                if self.a <= res <= self.ival_c.b:
                    return res
        est2 = self.estimator_q3(self.ival_d)
        if est2.a <= 0:
            d2 = self.delta_second()
            if d2.a < 0: d2.a = int_arith.c_zero
            res = self.root_second_left(d2).a
            if res < self.ival_c.a:
                print('warning: left end q2 < self.c')
                return self.ival_c.a
            return res
        else:
            d2 = self.delta_second()
            if self.under_est_der_le_0(self.ival_c, self.ival_d) and d2.a >= 0:
                res = self.root_second_left(d2).a
                if self.ival_c.a <= res <= self.ival_d.b:
                    return res
        est3 = self.estimator_q3(self.ival_b)
        if est3.a <= 0:
            d3 = self.delta_third()
            if d3.a < 0: d3.a = int_arith.c_zero
            res = self.root_third_left(d3).a
            if res < self.ival_d.a:
                print('warning: left end q3 < self.d')
                return self.ival_d.a
            assert self.ival_d.a <= res <= self.b
            return res
        else:
            d3 = self.delta_third()
            if self.under_est_der_le_0(self.ival_d, self.ival_b) and d3.a >= 0:
                res = self.root_third_left(d3).a
                if self.ival_d.a <= res <= self.b:
                    return res

        # d3 = self.delta_third()
        # if d3 >= 0:
        #     return self.root_third_left(d3)
        return None

    def get_right_end_upper_bound(self) -> dec.Decimal | None:
        """get the first root of the upper bound"""
        if self.ival_c.a < self.a:
            print('warning: c*<a')
            return self.b
        if self.estimator_q1(self.ival_c).a >= 0:
            d1 = self.delta_first()
            if d1.a < 0: d1.a = int_arith.c_zero
            res = self.root_first_right(d1).b
            if res > self.ival_c.b:
                print('warning: right end upper bound q1 < self.c')
                return self.ival_c.b
            assert self.a <= res <= self.ival_c.b
            return res
        else:
            d1 = self.delta_first()
            if self.under_est_der_le_0(self.ival_a, self.ival_c) and d1.a >= 0:
                res = self.root_first_right(d1).b
                if self.a <= res <= self.ival_c.b:
                    return res
        if self.estimator_q3(self.ival_d).a >= 0:
            d2 = self.delta_second()
            if not d2.a >= 0:
                return self.ival_d.b
            res = self.root_second_right(d2).b
            if res > self.ival_d.b:
                print('warning: right end upper bound q2 < self.d')
                return self.ival_d.b
            assert self.ival_c.a <= res <= self.ival_d.b
            return res
        else:
            d2 = self.delta_second()
            if self.under_est_der_le_0(self.ival_c, self.ival_d) and d2.a >= 0:
                res = self.root_second_right(d2).b
                if self.ival_c.a <= res <= self.ival_d.b:
                    return res
        if self.estimator_q3(self.ival_b).a >= 0:
            d3 = self.delta_third()
            res = self.root_third_right(d3).b
            if res > self.b:
                print('warning: right end upper bound q3 < self.b')
                return self.b
            assert self.ival_d.a <= res <= self.b
            return res
        else:
            d3 = self.delta_third()
            if self.under_est_der_le_0(self.ival_d, self.ival_b) and d3.a >= 0:
                res = self.root_third_right(d3).b
                if self.ival_d.a <= res <= self.b:
                    return res

        return self.b

    def get_right_end_under_bound(self):
        """get the last root of the under bound"""
        if self.ival_d.b > self.b:
            print('warning: d*>b')
            return self.b
        est3 = self.estimator_q3(self.ival_d)
        if est3.a <= 0:
            d3 = self.delta_third()
            if d3.a < 0: d3.a = int_arith.c_zero
            res = self.root_third_right(d3).b
            if res > self.b:
                print('warning: right end under bound q3 < self.b')
                return self.b
            if self.ival_d.a > res:
                print(self.root_third_right(d3))
                print(self.root_third_left(d3))
                print('err')
            assert self.ival_d.a <= res <= self.b
            return res
        else:
            d3 = self.delta_third()
            if d3.a < 0: d3.a = int_arith.c_zero
            if self.under_est_der_le_0(self.ival_d, self.ival_b) and d3.a >= 0:
                res = self.root_third_right(d3).b
                if self.ival_d.a <= res <= self.b:
                    return res
        est2 = self.estimator_q1(self.ival_c)
        if est2.a <= 0:
            d2 = self.delta_second()
            if not d2.a >= 0:
                return self.ival_d.b
            res = self.root_second_right(d2).b
            if res > self.ival_d.b:
                print('warning: right end under bound q2 < self.d')
                return self.ival_d.b
            assert self.ival_c.a <= res <= self.ival_d.b
            return res
        else:
            d2 = self.delta_second()
            if d2.a < 0: d2.a = int_arith.c_zero
            if self.under_est_der_le_0(self.ival_c, self.ival_d) and d2.a >= 0:
                res = self.root_second_right(d2).b
                if self.ival_c.a <= res <= self.ival_d.b:
                    return res
        est1 = self.estimator_q1(self.ival_a)
        if est1.a <= 0:
            d1 = self.delta_first()
            if d1.a < 0: d1.a = int_arith.c_zero
            res = self.root_first_right(d1).b
            if res > self.ival_c.b:
                print('warning: right end under bound q1 < self.c')
                return self.ival_c.b
            assert self.a <= res <= self.ival_c.b
            return res
        else:
            d1 = self.delta_first()
            if d1.a < 0: d1.a = int_arith.c_zero
            if self.under_est_der_le_0(self.ival_a, self.ival_c) and d1.a >= 0:
                res = self.root_first_right(d1).b
                if self.a <= res <= self.ival_c.b:
                    return res

        return None
