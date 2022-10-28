import sub as sb
import interval as ival
import psqe_under as ps
import FZCP.psqe_bounds as fps

class PSQEData:
    """
    The subproblem data used on PSQE subproblems
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


class PSQEProcessor:
    """
    The processor based on piece-wise underestimator
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
        self.use_symm_lipint = use_symm_lipint
        self.global_lipint = global_lipint
        self.rec_v = rec_v
        self.rec_x = rec_x
        self.problem = problem
        self.eps = eps
        self.ddi = problem.ddf(ival.Interval([problem.a, problem.b]))
        self.nddi = problem.nddf(ival.Interval([problem.a, problem.b]))

    def compute_bounds(self, sub):
        psqe = self.compute_under(sub)
        sub.data.split_point, sub.bound = psqe.lower_bound_and_point()
        x, v = psqe.record_and_point()
        if v < self.rec_v:
            self.rec_x = x
            self.rec_v = v

    def compute_upper(self, sub):
        if self.global_lipint:
            ddi = self.nddi
        else:
            nddi = self.problem.nddf(sub.data.ival)
        if self.use_symm_lipint:
            L = max(-self.nddi.x[0], self.nddi.x[1])
            nddi = ival.Interval([-L, L])
        return ps.PSQE_Under(sub.data.ival[0], sub.data.ival[1], nddi[0], nddi[1], self.problem.nobj, self.problem.ndf)

    def compute_under(self, sub):
        if self.global_lipint:
            ddi = self.ddi
        else:
            ddi = self.problem.ddf(sub.data.ival)
        if self.use_symm_lipint:
            L = max(-self.ddi.x[0], self.ddi.x[1])
            ddi = ival.Interval([-L, L])
        return ps.PSQE_Under(sub.data.ival[0], sub.data.ival[1], ddi[0], ddi[1], self.problem.objective,
                             self.problem.df)

    def updateSplitAndBounds(self, sub):
        psqe_upper = self.compute_upper(sub)
        max_x, sub.bound[1] = psqe_upper.lower_bound_and_point()
        sub.bound[1] = -sub.bound[1]
        sub.data.upperOfB = -psqe_upper.estimator(sub.data.ival[1])
        min_x, sub.bound[0] = self.compute_under(sub).lower_bound_and_point()
        widthX = sub.data.ival[1] - sub.data.ival[0]
        widthF = sub.bound[1] - sub.bound[0]
        beta = sub.bound[1] / widthF
        if beta <= 0.33:
            sub.data.split_point = sub.data.ival[0] + 0.33 * widthX
        elif beta <= 0.66:
            sub.data.split_point = sub.data.ival[0] + beta * widthX
        else:
            sub.data.split_point = sub.data.ival[0] + 0.66 * widthX

    def process(self, sub):
        """
        Process a subproblem
        Args:
            sub: subproblem to process

        Returns:
            list of generated subproblems

        """
        lst = []
        if sub.bound < self.rec_v - self.eps:
            sub_1 = sb.Sub(sub.level + 1, 0, PSQEData(ival.Interval([sub.data.ival.x[0], sub.data.split_point]), None))
            self.compute_bounds(sub_1)
            sub_2 = sb.Sub(sub.level + 1, 0, PSQEData(ival.Interval([sub.data.split_point, sub.data.ival.x[1]]), None))
            self.compute_bounds(sub_2)
            if sub_1.bound < self.rec_v - self.eps:
                lst.append(sub_1)
            if sub_2.bound < self.rec_v - self.eps:
                lst.append(sub_2)
        return lst
