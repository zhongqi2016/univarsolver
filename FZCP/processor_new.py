import psqe_bounds as psqe
import psl_bounds as psl
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


class ProcessorNew:

    def __init__(self, rec_v, rec_x, problem, eps, global_lipint=False, use_symm_lipint=False, estimator=2):
        """
        Initializes processor
        Args:
            rec_v: record value
            rec_x: record point
            problem: problem to solve
            eps: tolerance
            global_lipint: if True use global Lipschitz constant computed for the whole interval
            use_symm_lipint: if True use [-L,L], where L = max(|a|,|b|)
            estimator:estimator==1 -> Algorithm Piyavksii, estimator==2 -> PSQE
        """
        self.res_list = []
        self.use_symm_lipint = use_symm_lipint
        self.global_lipint = global_lipint
        self.rec_v = rec_v
        self.rec_x = rec_x
        self.problem = problem
        self.eps = eps
        self.ddi = problem.ddf(ival.Interval([problem.a, problem.b]))
        self.di = problem.df(ival.Interval([problem.a, problem.b]))
        self.estimator = estimator

    def compute_bounds(self, sub, under):
        """
        Args:
            under: if True compute the under bound
        """

        if self.estimator == 1:
            if self.global_lipint:
                di = self.di
            else:
                di = self.problem.df(sub.data.ival)
            if self.use_symm_lipint:
                L = max(-di[0], di[1])
                di = ival.Interval([-L, L])
            return psl.PSL_Bounds(sub.data.ival[0], sub.data.ival[1], di[0], di[1], self.problem.objective,
                                  under)
        else:
            if self.global_lipint:
                ddi = self.ddi
            else:
                ddi = self.problem.ddf(sub.data.ival)
            if self.use_symm_lipint:
                L = max(-ddi.x[0], ddi.x[1])
                ddi = ival.Interval([-L, L])
            return psqe.PSQE_Bounds(sub.data.ival[0], sub.data.ival[1], ddi[0], ddi[1], self.problem.objective,
                                    self.problem.df, under)

    def fzcp_process(self, sub):
        """
        Process of branching
        """
        lst = []
        obj = self.problem.objective
        width_of_interval = sub.data.ival.x[1] - sub.data.ival.x[0]
        in_bounds = self.update_interval(sub)
        if in_bounds and sub.data.ival[0] <= self.rec_x:
            if sub.data.ival.x[1] - sub.data.ival.x[0] < self.eps and obj(sub.data.ival.x[1]) <= 0:
                # If width of the interval satisfies the precision requirement
                self.res_list.append(sub.data.split_point)
            else:
                new_width = sub.data.ival.x[1] - sub.data.ival.x[0]
                # print(new_width / width_of_interval)
                if new_width / width_of_interval > 0.7:
                    sub_1 = sb.Sub(sub.level + 1, [0, 0],
                                   PSQEData(ival.Interval([sub.data.ival.x[0], sub.data.split_point]), None))
                    if obj(sub_1.data.ival[1]) <= 0 and sub_1.data.ival[1] < self.rec_x:
                        self.rec_x = sub_1.data.ival[1]
                    else:
                        sub_2 = sb.Sub(sub.level + 1, [0, 0],
                                       PSQEData(ival.Interval([sub.data.split_point, sub.data.ival.x[1]]), None))

                        lst.append(sub_2)
                    lst.append(sub_1)
                else:
                    if obj(sub.data.ival[1]) <= 0 and sub.data.ival[1] < self.rec_x:
                        self.rec_x = sub.data.ival[1]
                    lst.append(sub)

        return lst

    def update_interval(self, sub):
        """
        Narrow the interval by the zeros of the upper and under bounds
        """
        upper_estimator = self.compute_bounds(sub, False)
        max_x, sub.bound[1] = upper_estimator.lower_bound_and_point()
        lower_estimator = self.compute_bounds(sub, True)
        min_x, sub.bound[0] = lower_estimator.lower_bound_and_point()
        if sub.bound[0] <= 0 <= sub.bound[1]:
            left_end = lower_estimator.getNewTrialPoint()
            right_end = upper_estimator.getNewTrialPoint()
            sub.data.ival.x[0] = left_end
            if right_end < sub.data.ival.x[1]:
                sub.data.ival.x[1] = right_end
            sub.data.split_point = sub.data.ival.x[0] + (sub.data.ival.x[1] - sub.data.ival.x[0]) / 2
            return True
            # print("left:", self.problem.objective(sub.data.ival.x[0]))
            # print("[", sub.data.ival.x[0], ", ", sub.data.ival.x[1], "],", sub.data.split_point)
        else:
            return False
