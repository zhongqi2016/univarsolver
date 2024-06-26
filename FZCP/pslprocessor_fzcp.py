import sys

sys.path.append("..")
import sub as sb
import interval as ival
import psl_bounds as ps


class PSLData:
    """
    The subproblem data used on PSL subproblems
    """

    def __init__(self, ival, split_point):
        """
        The constructor
        Args:
            ival: the subproblem's interval
            split_point: the point to split interval
        """
        self.ival = ival
        self.split_point = split_point


class PSLProcessor:
    """
    The processor based on piece-wise linear underestimator
    """

    def __init__(self, rec_v, rec_x, problem, eps, global_lipint=False, use_symm_lipint=False):
        """
        Initializes processor
        Args:
            rec_v: record value
            rec_x: record point
            problem: problem to solve
            eps: tolerance
            global_lipint: if True use global Lipschitz constant computed for the whole interval
            use_symm_lipint: if True use [-L,L], where L = max(|a|,|b|)
        """
        self.res_list = []
        self.use_symm_lipint = use_symm_lipint
        self.global_lipint = global_lipint
        self.rec_v = rec_v
        self.rec_x = rec_x
        self.problem = problem
        self.eps = eps
        self.di = problem.df(ival.Interval([problem.a, problem.b]))

    def compute_bounds(self, sub, under):
        if self.global_lipint:
            di = self.di
        else:
            di = self.problem.df(sub.data.ival)
        if self.use_symm_lipint:
            L = max(-di[0], di[1])
            di = ival.Interval([-L, L])

        return ps.PSL_Bounds(sub.data.ival[0], sub.data.ival[1], di[0], di[1], self.problem.objective, under)

    def updateSplitAndBounds(self, sub):
        psqe_upper = self.compute_bounds(sub, False)
        max_x, sub.bound[1] = psqe_upper.lower_bound_and_point()
        min_x, sub.bound[0] = self.compute_bounds(sub, True).lower_bound_and_point()
        widthX = sub.data.ival[1] - sub.data.ival[0]
        widthF = sub.bound[1] - sub.bound[0]
        beta = sub.bound[1] / widthF
        if beta <= 0.33:
            sub.data.split_point = sub.data.ival[0] + 0.33 * widthX
        elif beta <= 0.66:
            sub.data.split_point = sub.data.ival[0] + beta * widthX
        else:
            sub.data.split_point = sub.data.ival[0] + 0.66 * widthX

    def fzcp_process(self, sub):
        lst = []
        obj = self.problem.objective
        obj(1)
        if sub.bound[0] <= 0 <= sub.bound[1] and sub.data.ival[0] < self.rec_x:
            if sub.data.ival.x[1] - sub.data.ival.x[0] < self.eps and obj(sub.data.ival.x[1]) <= 0:
                self.res_list.append(sub.data.split_point)
            else:
                sub_1 = sb.Sub(sub.level + 1, [0, 0],
                               PSLData(ival.Interval([sub.data.ival.x[0], sub.data.split_point]), None))
                self.updateSplitAndBounds(sub_1)
                sub_2 = sb.Sub(sub.level + 1, [0, 0],
                               PSLData(ival.Interval([sub.data.split_point, sub.data.ival.x[1]]), None))
                self.updateSplitAndBounds(sub_2)
                lst.append(sub_2)
                if obj(sub_1.data.ival[1]) <= 0 and sub_1.data.ival[1] < self.rec_x:
                    self.rec_x = sub_1.data.ival[1]
                lst.append(sub_1)
        return lst
