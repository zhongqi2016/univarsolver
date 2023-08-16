"""
The piece-wise linear underestimator from Casado, L. G., MartÍnez, J. A., GarcÍa, I., & Sergeyev, Y. D. (2003). New interval analysis support functions using gradient information in a global minimization algorithm. Journal of Global Optimization, 25(4), 345-362.
"""

cdef class PSL_Bounds:
    """
    Piecewise linear underestimator
    """
    cdef double a, b, alp, bet, fa, fb, c
    cdef bint under
    def __init__(self, a, b, alp, bet, fa, fb, under):
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
        self.under = under
        if under:
            self.alp = alp
            self.bet = bet
            self.fa = fa
            self.fb = fb
        else:
            self.alp = -bet
            self.bet = -alp
            self.fa = -fa
            self.fb = -fb
        if alp != bet:
            self.c = (self.fa - self.fb + self.bet * self.b - self.alp * self.a) / (self.bet - self.alp)

    def __repr__(self):
        return "Piecewise linear estimator " + "a = " + str(self.a) + ", b = " + str(self.b) + ", c = " + str(
            self.c) + ", alp = " + str(self.alp) + ", bet = " + str(self.bet) \
            + ", fa = " + str(self.fa) + ", fb = " + str(self.fb)

    cpdef get_fb(self):
        return self.fb

    cpdef estimator(self, double x):
        """
        The piecewise linear underestimator
        Args:
            x: argument within [a,b]

        Returns: underestimator's value
        """
        cdef double res
        if x <= self.c:
            res = self.fa + self.alp * (x - self.a)
        else:
            res = self.fb + self.bet * (x - self.b)
        return res

    cpdef nestimator(self, double x):
        return -self.estimator(x)

    cpdef public lower_bound_and_point(self):
        """
        Returns: Tuple (point where it is achieved, lower bound on interval [a,b])
        """
        cdef double record_x
        cdef double record_v
        if self.alp >= 0:
            record_x = self.a
            record_v = self.fa
        elif self.bet <= 0:
            record_x = self.b
            record_v = self.fb
        else:
            record_x = self.c
            record_v = self.estimator(self.c)
        if not self.under:
            record_v = -record_v
        return record_x, record_v

    cpdef public lower_bound_and_point2(self):
        """
        Returns: Tuple (point where it is achieved, lower bound on interval [a,b])
        """
        cdef double record_v
        if self.alp >= 0:
            # record_x = self.a
            record_v = self.fa
        elif self.bet <= 0:
            # record_x = self.b
            record_v = self.fb
        else:
            # record_x = self.c
            record_v = self.estimator(self.c)

        if record_v <= 0:
            return 0
        else:
            return -1

    cpdef record_and_point(self):
        """
        Returns: Tuple (point c where the best value of objective is achieved (a or b), f(c))
        """
        return (self.a, self.fa) if self.fa <= self.fb else (self.b, self.fb)

    cpdef public get_left_end(self):
        if self.alp >= 0:
            return None
        cdef double root_of_left_part = self.a - self.fa / self.alp
        if self.alp == self.bet:
            print('Nan')
            if root_of_left_part <= self.b:
                return root_of_left_part
            else:
                return None
        if root_of_left_part <= self.c:
            return root_of_left_part
        cdef double root_of_right_part
        if self.bet < 0:
            root_of_right_part = self.b - self.fb / self.bet
            if root_of_right_part <= self.b:
                return root_of_right_part
        return None
    cpdef public get_right_end2(self):
        if self.bet < 0:
            return None
        cdef double root_of_right_part = self.b - self.fb / self.bet
        if self.alp == self.bet:
            print('Nan')
            if root_of_right_part >= self.a:
                return root_of_right_part
            else:
                return None
        if root_of_right_part >= self.c:
            return root_of_right_part
        cdef double root_of_left_part
        if self.bet < 0:
            root_of_left_part = self.a - self.fa / self.alp
            if root_of_left_part >= self.a:
                return root_of_left_part
        return None

    # cpdef public get_right_end(self):
    #     cdef double root_of_right_part = self.b - self.fb / self.bet
    #     cdef double root_of_left_part
    #     if root_of_right_part <= self.b:
    #         if root_of_right_part >= self.c:
    #             return root_of_right_part
    #         else:
    #             root_of_left_part = self.a - self.fa / self.alp
    #             return root_of_left_part
    #     else:
    #         return self.b
    cpdef public get_right_end(self):
        cdef double root_of_right_part = self.b - self.fb / self.bet
        cdef double root_of_left_part
        if root_of_right_part <= self.b:
            if root_of_right_part >= self.c:
                return root_of_right_part
            else:
                root_of_left_part = self.a - self.fa / self.alp
                return root_of_left_part
        else:
            return None
