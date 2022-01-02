import sub as sb
import interval as ival
import psl_under as ps

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
    def __init__(self, rec_v, rec_x, problem, eps, global_lipint = False, use_symm_lipint = False):
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
        self.di = problem.df(ival.Interval([problem.a, problem.b]))

    def compute_bounds(self, sub):
        if self.global_lipint:
            di = self.di
        else:
            di = self.problem.df(sub.data.ival)
        if self.use_symm_lipint:
            L = max(-di[0], di[1])
            di = ival.Interval([-L,L])
        psl = ps.PSL_Under(sub.data.ival[0], sub.data.ival[1], di[0], di[1], self.problem.objective)
        sub.data.split_point, sub.bound = psl.lower_bound_and_point()
        x, v = psl.record_and_point()
        if v < self.rec_v:
            self.rec_x = x
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
        if sub.bound < self.rec_v - self.eps:
            sub_1 = sb.Sub(sub.level + 1, 0, PSLData(ival.Interval([sub.data.ival.x[0], sub.data.split_point]), None))
            self.compute_bounds(sub_1)
            sub_2 = sb.Sub(sub.level + 1, 0, PSLData(ival.Interval([sub.data.split_point, sub.data.ival.x[1]]), None))
            self.compute_bounds(sub_2)
            if sub_1.bound < self.rec_v - self.eps:
                lst.append(sub_1)
            if sub_2.bound < self.rec_v - self.eps:
                lst.append(sub_2)
        return lst
