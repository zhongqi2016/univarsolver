import sub as sb
import interval as ival
import numpy as np

class IntervalProcessor:
    """    Simple Interval processor that uses interval bounds to prune subproblems
    """
    def __init__(self, rec_v, rec_x, problem, eps):
        """
        Initializes processor
        Args:
            rec_v: record value
            rec_x: record point
            problem: problem to solve
            eps: tolerance
        """
        self.rec_v = rec_v
        self.rec_x = rec_x
        self.problem = problem
        self.eps = eps

    def compute_bounds(self, sub):
        sub.bound = self.problem.objective(sub.data).x[0]
        c = sub.data.mid()
        v = self.problem.objective(c)
        if np.isnan(self.rec_v) or v < self.rec_v:
            self.rec_x = c
            self.rec_v = v

    def process(self, sub):
        """
        Process a subproblem
        Args:
            sub: subproblem to process

        Returns:
            list of generated subproblems

        """
        lst = []
        #TODO: separate discard and split routines
        if sub.bound < self.rec_v - self.eps:
            sub_1 = sb.Sub(sub.level + 1, 0, ival.Interval([sub.data.x[0], sub.data.mid()]))
            self.compute_bounds(sub_1)
            sub_2 = sb.Sub(sub.level + 1, 0, ival.Interval([sub.data.mid(), sub.data.x[1]]))
            self.compute_bounds(sub_2)
            if sub_1.bound < self.rec_v - self.eps:
                lst.append(sub_1)
            if sub_2.bound < self.rec_v - self.eps:
                lst.append(sub_2)
        return lst

