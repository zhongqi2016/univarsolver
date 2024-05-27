"""
The piece-wise linear underestimator from Casado, L. G., MartÍnez, J. A., GarcÍa, I., & Sergeyev, Y. D. (2003). New interval analysis support functions using gradient information in a global minimization algorithm. Journal of Global Optimization, 25(4), 345-362.
"""
import interval_arithmetics as int_arith
import decimal as dec


def to_dec(a: float) -> int_arith.Interval:
    return int_arith.Interval(dec.Decimal(a), dec.Decimal(a))


class PSL_Bounds:
    """
    Piecewise linear underestimator
    """

    def __init__(self, a: dec.Decimal, b: dec.Decimal, alp: dec.Decimal, bet: dec.Decimal, ival_fa: int_arith.Interval,
                 ival_fb: int_arith.Interval, under: bool):
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
            self.ival_fa = ival_fa
            self.ival_fb = ival_fb

        else:
            self.alp = -bet
            self.bet = -alp
            self.ival_fa = -ival_fa
            self.ival_fb = -ival_fb
        self.ival_a = int_arith.Interval(self.a, self.a)
        self.ival_b = int_arith.Interval(self.b, self.b)
        self.fa = self.ival_fa.a
        self.fb = self.ival_fb.a
        self.ival_alp = int_arith.Interval(self.alp, self.alp)
        self.ival_bet = int_arith.Interval(self.bet, self.bet)

        self.ival_c = (self.ival_fa - self.ival_fb + self.ival_bet * self.ival_b - self.ival_alp * self.ival_a) / (
                self.ival_bet - self.ival_alp)

    def __repr__(self):
        return "Piecewise linear estimator " + "a = " + str(self.a) + ", b = " + str(self.b) + ", c = " + str(
            self.c) + ", alp = " + str(self.alp) + ", bet = " + str(self.bet) \
            + ", fa = " + str(self.fa) + ", fb = " + str(self.fb)

    def normal_cd(self) -> bool:
        return True
    def estimator_l1(self, x: dec.Decimal):
        ival_x = int_arith.Interval(x, x)
        """ left part of linear estimator"""
        return self.ival_fa + self.ival_alp * (ival_x - self.ival_a)

    def estimator_l2(self, x: dec.Decimal):
        ival_x = int_arith.Interval(x, x)
        """ right part of linear estimator"""
        return self.ival_fb + self.ival_bet * (ival_x - self.ival_b)

    def estimator(self, x) -> int_arith.Interval:
        """
        The piecewise linear underestimator
        Args:
            x: argument within [a,b]

        Returns: underestimator's value
        """
        if x < self.ival_c.a:
            res = self.estimator_l1(x)
        elif x > self.ival_c.b:
            res = self.estimator_l2(x)
        else:
            res = max(self.estimator_l1(x), self.estimator_l2(x))
        return res

    def get_fb(self):
        return self.fb

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
            record_x = self.ival_c.b
            record_v = self.estimator_l2(self.ival_c.b).a
        if not self.under:
            record_v = -record_v
        return record_x, record_v

    def bound_cross_zero(self):
        """
            Returns: if lower bound under zero return true,else return false.
            if upper bound above zero return true,else return false.
        """
        # record_x = None
        # record_v = None
        if self.alp >= 0:
            # record_x = self.a
            record_v = self.fa
        elif self.bet <= 0:
            # record_x = self.b
            record_v = self.fb
        else:
            # record_x = self.c
            record_v = self.estimator_l2(self.ival_c.b)
        if record_v.b <= 0:
            return True
        else:
            return False

    def record_and_point(self):
        """
        Returns: Tuple (point c where the best value of objective is achieved (a or b), f(c))
        """
        return (self.a, self.fa) if self.fa <= self.fb else (self.b, self.fb)

    def get_left_end(self):
        if self.alp >= 0:
            return None
        root_of_left_part = self.ival_a - self.ival_fa / self.ival_alp
        if root_of_left_part.a <= self.ival_c.a:
            return root_of_left_part.a
        if self.bet < 0:
            root_of_right_part = self.ival_b - self.ival_fb / self.ival_bet
            if root_of_right_part.a <= self.b:
                return root_of_right_part.a
        else:
            return None

    def get_right_end_upper_bound(self):
        if self.alp>0:
            root_of_left_part = self.ival_a - self.ival_fa / self.ival_alp
            if root_of_left_part.b<=self.ival_c.b:
                return root_of_left_part.b
        root_of_right_part = self.ival_b - self.ival_fb / self.ival_bet
        if root_of_right_part.b <= self.b:
            if root_of_right_part.b >= self.ival_c.a:
                return root_of_right_part.b
        else:
            return self.b

    def get_right_end_under_bound(self):
        if self.bet < 0:
            return None
        root_of_right_part = self.ival_b - self.ival_fb / self.ival_bet
        if root_of_right_part.b >= self.ival_c.b:
            return root_of_right_part.b
        root_of_left_part = self.ival_a - self.ival_fa / self.ival_alp
        if root_of_left_part.b <= self.ival_c.b:
            return root_of_left_part.b
        return None
