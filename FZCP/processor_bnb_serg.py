import psqe_bounds as fps
import sys

sys.path.append("..")
import sub as sb
import interval as ival


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


class ProcessorBNBSerg:

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
        self.ddi = problem.ddf(ival.Interval([problem.a, problem.b]))

    def compute_bounds(self, sub, under):
        """
        Args:
            under: if True compute the under bound
        """
        if self.global_lipint:
            ddi = self.ddi
        else:
            ddi = self.problem.ddf(sub.data.ival)
        if self.use_symm_lipint:
            L = max(-self.ddi.x[0], self.ddi.x[1])
            ddi = ival.Interval([-L, L])
        return fps.PSQE_Bounds(sub.data.ival[0], sub.data.ival[1], ddi[0], ddi[1], self.problem.objective,
                               self.problem.df, under)

    def fzcp_process(self, sub):
        lst = []
        obj = self.problem.objective
        if sub.bound[0] <= 0 <= sub.bound[1] and sub.data.ival[0] < self.rec_x:
            if sub.data.ival.x[1] - sub.data.ival.x[0] < self.eps and obj(sub.data.ival.x[1]) <= 0:
                self.res_list.append(sub.data.split_point)
            else:
                sub_1 = sb.Sub(sub.level + 1, [0, 0],
                               PSQEData(ival.Interval([sub.data.ival.x[0], sub.data.split_point]), None))

                self.update_interval(sub_1)
                if obj(sub_1.data.ival[1]) <= 0 and sub_1.data.ival[1] < self.rec_x:
                    self.rec_x = sub_1.data.ival[1]
                else:
                    sub_2 = sb.Sub(sub.level + 1, [0, 0],
                                   PSQEData(ival.Interval([sub.data.split_point, sub.data.ival.x[1]]), None))

                    self.update_interval(sub_2)
                    lst.append(sub_2)
                lst.append(sub_1)

        return lst

    def update_interval(self, sub):
        psqe_upper = self.compute_bounds(sub, False)
        max_x, sub.bound[1] = psqe_upper.lower_bound_and_point()
        psqe_under = self.compute_bounds(sub, True)
        min_x, sub.bound[0] = psqe_under.lower_bound_and_point()
        if sub.bound[0] <= 0 <= sub.bound[1]:
            left_end = psqe_under.getNewTrialPoint()
            right_end = psqe_upper.getNewTrialPoint()
            sub.data.ival.x[0] = left_end
            if right_end < sub.data.ival.x[1]:
                sub.data.ival.x[1] = right_end
            sub.data.split_point = sub.data.ival.x[0] + (sub.data.ival.x[1] - sub.data.ival.x[0]) / 2
            # print("left:", self.problem.objective(sub.data.ival.x[0]))
            # print("[", sub.data.ival.x[0], ", ", sub.data.ival.x[1], "],", sub.data.split_point)