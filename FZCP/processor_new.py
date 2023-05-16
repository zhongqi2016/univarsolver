import psqe_bounds_cy as psqe
import psl_bounds_cy as psl
import sys

sys.path.append("..")
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

    def __init__(self, rec_v, rec_x, problem, eps, global_lipint=False, use_symm_lipint=False, estimator=2,
                 reduction=True):
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
        self.reduction = reduction
        self.running = True

    def compute_bounds(self, sub_interval: ival.Interval, under: bool):
        """
        Args:
            sub_interval:
            under: if True compute the under bound
        """

        if self.estimator == 1:
            if self.global_lipint:
                di = self.di
            else:
                di = self.problem.df(sub_interval)
            if self.use_symm_lipint:
                L = max(-di[0], di[1])
                di = ival.Interval([-L, L])
            return self.Bounds(sub_interval.x[0], sub_interval.x[1], di[0], di[1], under)
        else:
            if self.global_lipint:
                ddi = self.ddi
            else:
                ddi = self.problem.ddf(sub_interval)
            if self.use_symm_lipint:
                L = max(-ddi.x[0], ddi.x[1])
                ddi = ival.Interval([-L, L])
            return self.Bounds(sub_interval.x[0], sub_interval.x[1], ddi[0], ddi[1], under)

    def Bounds(self, a: float, b: float, lip_alp: float, lip_bet: float, under: bool):
        if self.estimator == 1:
            return psl.PSL_Bounds(a, b, lip_alp, lip_bet, self.problem.objective(a), self.problem.objective(b),
                                  under)
        else:
            return psqe.PSQE_Bounds(a, b, lip_alp, lip_bet, self.problem.objective(a), self.problem.objective(b),
                                    self.problem.df(a), self.problem.df(b), under)

    def fzcp_process(self, sub_interval: ival.Interval):
        """
        Process of branching
        """

        lst = []
        obj = self.problem.objective
        # if self.rec_x - sub.data.ival.x[0] <= self.eps:
        #     self.res_list.append(sub.data.ival.x[0])
        #     self.running = False
        #     return lst
        if sub_interval.x[1] - sub_interval.x[0] <= self.eps:
            if obj(sub_interval.x[1]) <= 1e-14:
                self.res_list.append(sub_interval.x[0])
                self.running = False
                return lst
            else:
                return lst
        width_of_interval = sub_interval.x[1] - sub_interval.x[0]
        # in_bounds = self.update_interval(sub_interval)

        lower_estimator = self.compute_bounds(sub_interval, under=True)
        left_end = lower_estimator.get_left_end()
        upper_estimator = self.compute_bounds(sub_interval, under=False)
        right_end = upper_estimator.get_right_end()
        if right_end > sub_interval.x[1]:
            right_end = sub_interval.x[1]
        if left_end is not None and sub_interval.x[0] <= self.rec_x:
            if self.reduction is False:
                left_end = sub_interval.x[0]
                right_end = sub_interval.x[1]
            split_point = left_end + (right_end - left_end) / 2
            if right_end - left_end < self.eps and obj(right_end) <= 0:
                # If width of the interval satisfies the precision requirement
                self.res_list.append(split_point)
                self.running = False
            else:
                new_width = right_end - left_end
                # print(new_width / width_of_interval)
                if new_width / width_of_interval > 0.7:
                    sub_1 = ival.Interval([left_end, split_point])
                    if obj(sub_1.x[1]) <= 0 and sub_1.x[1] < self.rec_x:
                        self.rec_x = sub_1.x[1]
                    else:
                        sub_2 = ival.Interval([split_point, right_end])
                        # if obj(sub_2.data.ival[1]) <= 0:
                        #     self.rec_x = sub_2.data.ival[1]
                        lst.append(sub_2)
                    lst.append(sub_1)
                else:
                    # if obj(sub.data.ival[1]) <= 0 and sub.data.ival[1] < self.rec_x:
                    #     self.rec_x = sub.data.ival[1]
                    sub_interval[0] = left_end
                    sub_interval[1] = right_end
                    lst.append(sub_interval)

        return lst

    # def fzcp_process(self, sub):
    #     """
    #     Process of branching
    #     """
    #
    #     lst = []
    #     obj = self.problem.objective
    #     if self.rec_x - sub.data.ival.x[0] <= self.eps:
    #         self.res_list.append(sub.data.ival.x[0])
    #         self.running = False
    #         return lst
    #     if sub.data.ival.x[1] - sub.data.ival.x[0] <= self.eps:
    #         if obj(sub.data.ival.x[1]) <= 1e-14:
    #             self.res_list.append(sub.data.ival.x[0])
    #             self.running = False
    #             return lst
    #     width_of_interval = sub.data.ival.x[1] - sub.data.ival.x[0]
    #     in_bounds = self.update_interval(sub)
    #
    #     if in_bounds and sub.data.ival[0] <= self.rec_x:
    #         sub.data.split_point = sub.data.ival.x[0] + (sub.data.ival.x[1] - sub.data.ival.x[0]) / 2
    #
    #         new_width = sub.data.ival.x[1] - sub.data.ival.x[0]
    #         # print(new_width / width_of_interval)
    #         if new_width / width_of_interval > 0.7:
    #             sub_1 = sb.Sub(sub.level + 1, [0, 0],
    #                            PSQEData(ival.Interval([sub.data.ival.x[0], sub.data.split_point]), None))
    #             if obj(sub_1.data.ival[1]) <= 0 and sub_1.data.ival[1] < self.rec_x:
    #                 self.rec_x = sub_1.data.ival[1]
    #             else:
    #                 sub_2 = sb.Sub(sub.level + 1, [0, 0],
    #                                PSQEData(ival.Interval([sub.data.split_point, sub.data.ival.x[1]]), None))
    #                 if obj(sub_2.data.ival[1]) <= 0:
    #                     self.rec_x = sub_2.data.ival[1]
    #                 lst.append(sub_2)
    #             lst.append(sub_1)
    #         else:
    #             if obj(sub.data.ival[1]) <= 0 and sub.data.ival[1] < self.rec_x:
    #                 self.rec_x = sub.data.ival[1]
    #             lst.append(sub)
    #
    #     return lst

    def update_interval(self, sub_interval: ival.Interval):
        """
        Narrow the interval by the zeros of the upper and under bounds
        """
        lower_estimator = self.compute_bounds(sub_interval, under=True)
        left_end = lower_estimator.get_left_end()
        upper_estimator = self.compute_bounds(sub_interval, under=False)
        right_end = upper_estimator.get_right_end()
        le = sub_interval.x[0]
        re = sub_interval.x[1]
        if left_end is not None:
            if self.reduction:
                # upper_estimator = self.compute_bounds(sub, under=False)
                # right_end = upper_estimator.get_right_end()
                if left_end > sub_interval.x[0]:
                    le = left_end
                if right_end < sub_interval.x[1]:
                    re = right_end
            return True
        else:
            return False
