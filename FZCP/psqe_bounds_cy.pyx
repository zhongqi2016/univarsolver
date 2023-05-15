# Piecewise smooth quadratic estimators


cdef class PSQE_Bounds:
    """
    Piecewise quadratic underestimator
    """
    cdef double a, b, alp, bet, fa, fb, dfa, dfb, c, d
    cdef bint under

    def __init__(self, a, b, alp, bet, fa, fb, dfa, dfb, under: bool):
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
            self.fa = fa
            self.fb = fb
            self.dfa = dfa
            self.dfb = dfb
        else:
            self.alp = -bet
            self.bet = -alp
            self.fa = -fa
            self.fb = -fb
            self.dfa = -dfa
            self.dfb = -dfb
        cdef delt
        delt = (self.dfb - self.dfa - self.alp * (b - a)) / (bet - alp)
        # print("delt = ", delt)
        # self.c = ((delt - a) * self.dfa + (b - delt) * self.dfb + 0.5 * delt ** 2 * (
        #         self.bet - self.alp) + self.alp * delt * (
        #                   b - a) + 0.5 * self.alp * (a ** 2 - b ** 2) + self.fa - self.fb) / (
        #                  delt * (self.bet - self.alp))
        self.c = (self.alp * (b - a) + self.dfa - self.dfb) / (2 * (self.bet - self.alp)) + (
                self.fa - self.fb + b * self.dfb - a * self.dfa + self.alp * (a ** 2 - b ** 2) / 2) / (
                         self.dfb - self.dfa - self.alp * (b - a))
        self.d = self.c + delt
        # if self.c > self.b:
        #     print('error c>b')
        #     print('a=', self.a, 'b=', self.b, 'c=', self.c, 'delt=', delt, 'fa=', self.fa, 'fb=', self.fb, 'alp=',
        #           self.alp, 'beta=', self.bet)
        if self.d > self.b:
            print('error d>b')
            print('a=', self.a, 'b=', self.b, 'c=', self.c, 'delt=', delt, 'fa=', self.fa, 'fb=', self.fb, 'alp=',
                  self.alp, 'beta=', self.bet)

    def __repr__(self):
        return "Estimator " + "a = " + str(self.a) + ", b = " + str(self.b) + ", c = " + str(self.c) + ", d = " + str(
            self.d) + ", alp = " + str(self.alp) + ", bet = " + str(self.bet) + ", fa = " + str(
            self.fa) + ", fb = " + str(self.fb) + ", dfa = " + str(self.dfa) + ", dfb = " + str(self.dfb)

    cpdef estimator(self, double x):
        """
        The piecewise quadratic underestimator
        Args:
            x: argument

        Returns: underestimator's value
        """
        cdef double res
        if x <= self.c:
            res = self.fa + self.dfa * (x - self.a) + 0.5 * self.alp * (x - self.a) ** 2
        elif x < self.d:
            res = self.fa + self.dfa * (self.c - self.a) + 0.5 * self.alp * (self.c - self.a) ** 2 + (
                    self.dfa + self.alp * (self.c - self.a)) * (x - self.c) + 0.5 * self.bet * (x - self.c) ** 2
        else:
            res = self.fb + self.dfb * (x - self.b) + 0.5 * self.alp * (x - self.b) ** 2

        # if self.under:
        return res
        # else:
        #     return -res

    cpdef nestimator(self, double x):
        return -self.estimator(x)

    cdef estimators_derivative(self, double x):
        """
        The piecewise linear underestimator's derivative
        Args:
            x: argument

        Returns: underestimator's derivative value
        """
        cdef double res
        if x < self.c:
            res = self.dfa + self.alp * (x - self.a)
        elif x < self.d:
            res = self.dfa + self.alp * (self.c - self.a) + self.bet * (x - self.c)
        else:
            res = self.dfb + self.alp * (x - self.b)
        # if self.under:
        return res
        # else:
        #     return -res

    # def lower_bound_and_point(self):
    #     """
    #     Returns: Tuple (point where it is achieved, lower bound on interval [a,b])
    #     """
    #     x_list = [self.a, self.c, self.d, self.b]
    #     df_list = [self.estimators_derivative(x) for x in x_list]
    #     check_list = [self.a, self.b]
    #     record_x = None
    #     record_v = None
    #     ln = len(x_list)
    #     for i in range(ln - 1):
    #         x = self.find_argmin(x_list[i], df_list[i], x_list[i + 1], df_list[i + 1])
    #         if not (x is None):
    #             check_list.append(x)
    #     # print(check_list)
    #     for x in check_list:
    #         v = self.estimator(x)
    #         if record_v is None or v < record_v:
    #             record_v = v
    #             record_x = x
    #     if not self.under:
    #         record_v = -record_v
    #     return record_x, record_v

    # def lower_bound_and_point2(self):
    #     """
    #     Returns: Tuple (point where it is achieved, lower bound on interval [a,b])
    #     """
    #     if self.fb <= 0:
    #         return 0
    #     x_list = [self.a, self.c, self.d, self.b]
    #     df_list = [self.estimators_derivative(x) for x in x_list]
    #     ln = len(x_list)
    #
    #     for i in range(ln - 1):
    #         x = self.find_argmin(x_list[i], df_list[i], x_list[i + 1], df_list[i + 1])
    #         if (x is not None) and self.estimator(x) <= 0:
    #             return i + 1
    #
    #     return -1

    def record_and_point(self):
        """
        Returns: Tuple (point c where the best value of objective is achieved (a or b), f(c))
        """
        return (self.a, self.fa) if self.fa <= self.fb else (self.b, self.fb)

    # cdef find_argmin(self, double x1, double df1, double x2, double df2):
    #     cdef double xs
    #     if df1 == 0 and df2 == 0:
    #         xs = 0.5 * (x1 + x2)
    #     elif df1 <= 0 <= df2:
    #         xs = x1 + (-df1) * (x2 - x1) / (df2 - df1)
    #     else:
    #         xs = None
    #     return xs

    cdef find_min_under_zero(self, double x1, double df1, double x2, double df2):
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

    cdef find_max_above_zero(self, double x1, double df1, double x2, double df2):
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

    cdef right_root_first(self):
        """
        right root of the first row of the estimate equation
        """
        # print("first right")
        if self.alp < 0:
            return self.a + (-self.dfa - (self.dfa ** 2 - 2 * self.alp * self.fa) ** 0.5) / self.alp
        else:
            return self.a + (-self.dfa + (self.dfa ** 2 - 2 * self.alp * self.fa) ** 0.5) / self.alp

    cdef left_root_first(self):
        """
        left root of the first row of the estimate equation
        """
        # print("first left")
        if self.alp < 0:
            return self.a + (-self.dfa + (self.dfa ** 2 - 2 * self.alp * self.fa) ** 0.5) / self.alp
        else:
            return self.a + (-self.dfa - (self.dfa ** 2 - 2 * self.alp * self.fa) ** 0.5) / self.alp

    cdef right_root_third(self):
        """
        right root of the third row of the estimate equation
        """
        # print("third right")
        if self.alp < 0:
            return self.b + (-self.dfb - (self.dfb ** 2 - 2 * self.alp * self.fb) ** 0.5) / self.alp
        else:
            return self.b + (-self.dfb + (self.dfb ** 2 - 2 * self.alp * self.fb) ** 0.5) / self.alp

    cdef left_root_third(self):
        """
        left root of the third row of the estimate equation
        """
        # print("third left")
        if self.alp < 0:
            return self.b + (-self.dfb + (self.dfb ** 2 - 2 * self.alp * self.fb) ** 0.5) / self.alp
        else:
            return self.b + (-self.dfb - (self.dfb ** 2 - 2 * self.alp * self.fb) ** 0.5) / self.alp

    cdef right_root_second(self):
        """
        right root of the second row of the estimate equation
        """
        # print("second right")
        cdef double c, d, res
        c = self.fa + self.dfa * (self.c - self.a) + self.alp / 2 * (self.c - self.a) ** 2
        b = self.dfa + self.alp * (self.c - self.a)
        if self.bet > 0:
            res = self.c + (-b + (b ** 2 - 2 * self.bet * c) ** 0.5) / self.bet
        else:
            res = self.c + (-b - (b ** 2 - 2 * self.bet * c) ** 0.5) / self.bet
        return res

    cdef left_root_second(self):
        """
        left root of the second row of the estimate equation
        """
        # print("second left")
        cdef double c, d, res
        c = self.fa + self.dfa * (self.c - self.a) + self.alp / 2 * (self.c - self.a) ** 2
        b = self.dfa + self.alp * (self.c - self.a)
        if self.bet > 0:
            res = self.c + (-b - (b ** 2 - 2 * self.bet * c) ** 0.5) / self.bet
        else:
            res = self.c + (-b + (b ** 2 - 2 * self.bet * c) ** 0.5) / self.bet
        return res

    cdef under_est_der_le_0(self, double num1, double num2):
        """
        if(\phi'(num1) * \phi'(num2)<0) return true, else return false.
        """

        if self.estimators_derivative(num1) < 0 and self.estimators_derivative(num2) > 0:
            return True
        else:
            return False

    cdef upper_est_der_le_0(self, double num1, double num2):
        """
        if(\phi'(num1) * \phi'(num2)<0) return true, else return false.
        """
        if self.estimators_derivative(num1) > 0 and self.estimators_derivative(num2) < 0:
            return True
        else:
            return False

    cdef delta_first(self):
        return self.dfa ** 2 - 2 * self.alp * self.fa

    cdef delta_second(self):
        c = self.fa + self.dfa * (self.c - self.a) + self.alp / 2 * (self.c - self.a) ** 2
        b = self.dfa + self.alp * (self.c - self.a)
        return b ** 2 - 2 * self.bet * c

    cdef delta_third(self):
        return self.dfb ** 2 - 2 * self.alp * self.fb

    cdef root_first_left(self, d1):
        return self.a + (-self.dfa - d1 ** 0.5) / self.alp

    cdef root_second_left(self, d2):
        return self.c + (-self.dfa - self.alp * (self.c - self.a) - d2 ** 0.5) / self.bet

    cdef root_third_left(self, d3):
        return self.b + (-self.dfb - d3 ** 0.5) / self.alp

    cdef root_first_right(self, d1):
        return self.a + (-self.dfa + d1 ** 0.5) / self.alp

    cdef root_second_right(self, d2):
        return self.c + (-self.dfa - self.alp * (self.c - self.a) + d2 ** 0.5) / self.bet

    cdef root_third_right(self, d3):
        return self.b + (-self.dfb + d3 ** 0.5) / self.alp

    cpdef get_left_end(self):
        # if (flag1)==true, in first section of estimator have root.
        cdef double d1
        if self.estimator(self.c) <= 0:
            d1 = self.delta_first()
            return self.root_first_left(d1)
        else:
            d1 = self.delta_first()
            if self.under_est_der_le_0(self.a, self.c) and d1 >= 0:
                return self.root_first_left(d1)
        cdef double d2
        if self.estimator(self.d) <= 0:
            d2 = self.delta_second()
            return self.root_second_left(d2)
        else:
            d2 = self.delta_second()
            if self.under_est_der_le_0(self.c, self.d) and d2 >= 0:
                return self.root_second_left(d2)
        cdef double d3
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

    cpdef get_right_end(self):
        cdef double d1
        # if (flag1)==true, in first section of estimator have root.
        if self.estimator(self.c) >= 0:
            d1 = self.delta_first()
            return self.root_first_right(d1)
        else:
            d1 = self.delta_first()
            if self.upper_est_der_le_0(self.a, self.c) and d1 >= 0:
                return self.root_first_right(d1)
        cdef double d2
        if self.estimator(self.d) >= 0:
            d2 = self.delta_second()
            return self.root_second_right(d2)
        else:
            d2 = self.delta_second()
            if self.upper_est_der_le_0(self.c, self.d) and d2 >= 0:
                return self.root_second_right(d2)
        cdef double d3
        if self.estimator(self.b) >= 0:
            d3 = self.delta_third()
            return self.root_third_right(d3)
        else:
            d3 = self.delta_third()
            if self.upper_est_der_le_0(self.d, self.b) and d3 >= 0:
                return self.root_third_right(d3)

        return self.b

# def get_trial_point_under(self):
#     if self.estimator(self.c) <= 0:
#         if self.alp > 0:
#             new_point = self.left_root_first()
#         else:
#             new_point = self.right_root_first()
#     else:
#         est_der_c = self.estimators_derivative(self.c)
#         est_der_d = self.estimators_derivative(self.d)
#         if self.find_min_under_zero(self.c, est_der_c, self.d, est_der_d):
#             if self.bet > 0:
#                 new_point = self.left_root_second()
#             else:
#                 new_point = self.right_root_second()
#         else:
#             if self.alp > 0:
#                 new_point = self.left_root_third()
#             else:
#                 new_point = self.right_root_third()
#     return new_point

# def get_trial_point_upper(self):
# """
# Get trial point (zero point of upper estimate)
# """
# if self.estimator(self.c) >= 0:
#     if self.alp > 0:
#         new_point = self.right_root_first()
#     else:
#         new_point = self.left_root_first()
# else:
#     est_der_c = self.estimators_derivative(self.c)
#     est_der_d = self.estimators_derivative(self.d)
#     if self.find_max_above_zero(self.c, est_der_c, self.d, est_der_d):
#         if self.bet > 0:
#             new_point = self.right_root_second()
#         else:
#             new_point = self.left_root_second()
#     else:
#         if self.fb > 0:
#             if self.alp > 0:
#                 new_point = self.right_root_third()
#             else:
#                 new_point = self.left_root_third()
#         else:
#             new_point = self.b
# return new_point
