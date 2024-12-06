import math
import psqe_bounds_correctly as psqe
import psl_bounds_correctly as psl
import sys
import copy
import decimal as dec

sys.path.append("..")
import interval_arithmetics as ival


class ProcData:
    """
    The subproblem data used in algorithm
    """

    def __init__(self, sub_interval: ival.Interval, lip: ival.Interval, quadratic: bool, counter=0, period_comp_lip=0):
        """
        The constructor
        Args:
            sub_interval: the subproblem's interval
            split_point: the point to split interval
        """
        self.sub_interval = sub_interval
        self.lip = lip
        self.counter = counter
        self.period_comp_lip = period_comp_lip
        self.quadratic = quadratic


class ProcessorNew:

    def __init__(self, rec_v, rec_x, problem, eps, global_lipint=False, use_symm_lipint=False, estimator=2,
                 reduction=1, adaptive=False):
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
        self.ddi = problem.ddf(ival.Interval(dec.Decimal(problem.a), dec.Decimal(problem.b)))
        self.di = problem.df(ival.Interval(dec.Decimal(problem.a), dec.Decimal(problem.b)))
        self.estimator = estimator
        self.reduction = reduction
        self.running = True
        self.adaptive = adaptive
        self.sign_b = True if problem.objective(
            ival.Interval(dec.Decimal(problem.b), dec.Decimal(problem.b))).b > 0 else False

    def update_lipschitz(self, data: ProcData):
        """
        Args:
            data: data of sub_interval
        """
        if self.estimator == 1:
            data.lip = self.problem.df(data.sub_interval)
            if self.use_symm_lipint:
                L = max(-data.lip.a, data.lip.b)
                data.lip = ival.Interval(-L, L)
        else:
            data.lip = self.problem.ddf(data.sub_interval)
            if self.use_symm_lipint:
                L = max(-data.lip.a, data.lip.b)
                data.lip = ival.Interval(-L, L)

    def compute_bounds(self, data: ProcData, under: bool):
        a = data.sub_interval.a
        b = data.sub_interval.b
        ival_a = ival.Interval(a, a)
        ival_b = ival.Interval(b, b)
        if self.estimator == 1:
            return psl.PSL_Bounds(a, b, data.lip.a, data.lip.b, self.problem.objective(ival_a),
                                  self.problem.objective(ival_b), under)
        else:
            assert self.estimator == 2
            return psqe.PSQE_Bounds(a=a, b=b, alp=data.lip.a, bet=data.lip.b,
                                    ival_fa=self.problem.objective(ival_a),
                                    ival_fb=self.problem.objective(ival_b),
                                    ival_dfa=self.problem.df(ival_a), ival_dfb=self.problem.df(ival_b), under=under)

    def fzcp_process(self, data: ProcData):
        """
        Process of branching
        """
        sub_interval = data.sub_interval
        lst = []
        obj = self.problem.objective
        # if self.rec_x - sub.data.ival.x[0] <= self.eps:
        #     self.res_list.append(sub.data.ival.x[0])
        #     self.running = False
        #     return lst
        # if sub_interval.x[0] > self.rec_x:
        #     return []
        if sub_interval.b > self.rec_x: return lst
        if self.rec_x - sub_interval.a <= self.eps:
            self.res_list.append(ival.Interval(sub_interval.a, self.rec_x))
            # print("(%lf,%lf),f(x_r)=%lf" % (sub_interval.a, self.rec_x, obj(self.rec_x)))
            self.running = False
            return lst
        width_of_interval = sub_interval.b - sub_interval.a

        if not self.global_lipint:
            if self.adaptive:
                if data.counter < data.period_comp_lip:
                    data.counter += 1
                else:
                    self.update_lipschitz(data)
                    data.period_comp_lip = int(512 / math.sqrt(data.lip.b - data.lip.a))
                    # print(data.lip.b - data.lip.a, data.period_comp_lip)
                    data.counter = 0
            else:
                self.update_lipschitz(data)
        if self.estimator == 1 and data.quadratic:
            data.quadratic = False
            self.update_lipschitz(data)
        lower_estimator = self.compute_bounds(data, under=True)
        if not lower_estimator.normal_cd():
            self.estimator = 1
            data.quadratic = False
            self.update_lipschitz(data)
            lower_estimator = self.compute_bounds(data, under=True)
        left_end = lower_estimator.get_left_end()

        if left_end is not None:
            if self.reduction > 0:
                # if self.reduction == 2:
                #     upper_estimator = self.compute_bounds(data, under=False)
                #     right_end = upper_estimator.get_right_end()
                #     if right_end is None:
                #         right_end = lower_estimator.get_right_end2()
                if self.reduction == 1:
                    # if lower_estimator.get_fb() > 0:
                    if ((sub_interval.b < dec.Decimal(self.problem.b) and sub_interval.b < self.rec_x) or
                            (sub_interval.b == dec.Decimal(self.problem.b) and self.sign_b)):
                        right_end = lower_estimator.get_right_end_under_bound()
                        if right_end is None:
                            print('err')
                    elif sub_interval.b == self.rec_x:
                        upper_estimator = self.compute_bounds(data, under=False)
                        right_end = upper_estimator.get_right_end_upper_bound()
                        if right_end > sub_interval.b:
                            right_end = sub_interval.b
                        self.rec_x = right_end
                else:
                    right_end = None

            elif self.reduction == 0:
                left_end = sub_interval.a
                right_end = sub_interval.b
            split_point = left_end + (right_end - left_end) / 2
            if right_end - left_end < self.eps:
                # If width of the interval satisfies the precision requirement
                res_ival = ival.Interval(left_end, right_end)
                self.res_list.append(res_ival)
                # print("(%lf,%lf),lb=%lf,f(x2)=%lf" % (
                #     left_end, right_end, lower_estimator.lower_bound_and_point()[1], obj(right_end)))
            else:
                new_width = right_end - left_end
                # print(new_width / width_of_interval)
                if new_width / width_of_interval > 0.7:
                    sub_1 = ival.Interval(left_end, split_point)
                    if obj(ival.Interval(sub_1.b, sub_1.b)).b <= 0:
                        self.rec_x = sub_1.b
                    else:
                        data2 = ProcData(sub_interval=ival.Interval(split_point, right_end),
                                         lip=copy.deepcopy(data.lip), counter=data.counter,
                                         quadratic=data.quadratic,
                                         period_comp_lip=data.period_comp_lip)
                        lst.append(data2)

                    data.sub_interval = sub_1
                    lst.append(data)
                else:
                    data.sub_interval.a = left_end
                    data.sub_interval.b = right_end
                    lst.append(data)
        return lst
