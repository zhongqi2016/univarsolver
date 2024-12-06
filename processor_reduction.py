import math
import FZCP.psqe_bounds as psqe
import FZCP.psl_bounds as psl
import copy

import interval as ival


class ProcData:
    """
    The subproblem data used in algorithm
    """

    def __init__(self, sub_interval: ival.Interval, lip: ival.Interval, counter: int, period_comp_lip):
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


class ProcessorReduction:

    def __init__(self, rec_v, rec_x, problem, eps, global_lipint=False, use_symm_lipint=False, estimator=2,
                 reduction=1, adaptive=False, period_comp_lip=0):
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
        self.adaptive = adaptive

    def update_lipschitz(self, data: ProcData):
        """
        Args:
            data: data of sub_interval
        """
        if self.estimator == 1:
            data.lip = self.problem.df(data.sub_interval)
            if self.use_symm_lipint:
                L = max(-data.lip.x[0], data.lip.x[1])
                data.lip = ival.Interval([-L, L])
        else:
            data.lip = self.problem.ddf(data.sub_interval)
            if self.use_symm_lipint:
                L = max(-data.lip.x[0], data.lip.x[1])
                data.lip = ival.Interval([-L, L])

    def compute_bounds(self, data: ProcData, under: bool, dfa: float, dfb: float):
        a = data.sub_interval.x[0]
        b = data.sub_interval.x[1]
        if self.estimator == 1:
            return psl.PSL_Bounds(a=a, b=b, alp=data.lip.x[0], bet=data.lip.x[1],
                                  fa=self.problem.objective(a) - self.rec_v,
                                  fb=self.problem.objective(b) - self.rec_v, under=under)
        else:
            return psqe.PSQE_Bounds(a=a, b=b, alp=data.lip.x[0], bet=data.lip.x[1],
                                    fa=self.problem.objective(a) - self.rec_v,
                                    fb=self.problem.objective(b) - self.rec_v,
                                    dfa=dfa, dfb=dfb, under=under)

    def reduction2(self, data: ProcData, dfa: float, dfb: float) -> bool:
        """
        Return: interval will be reduced or not
        """
        if self.estimator == 1:
            if data.lip.include_zero():
                return False
            else:
                return True
        else:
            lip_left = data.lip.x[0]
            lip_right = data.lip.x[1]
            if dfa < 0:
                dfa = -dfa
                dfb = -dfb
                lip_left = -data.lip.x[1]
                lip_right = -data.lip.x[0]
            lower_estimator = psl.PSL_Bounds(a=data.sub_interval.x[0], b=data.sub_interval.x[1],
                                             alp=lip_left,
                                             bet=lip_right, fa=dfa, fb=dfb, under=True)
            left_end = lower_estimator.get_left_end()
            if left_end:
                if dfb > 0:
                    right_end = lower_estimator.get_right_end_under_bound()
                else:
                    upper_estimator = psl.PSL_Bounds(a=data.sub_interval.x[0], b=data.sub_interval.x[1],
                                                     alp=lip_left,
                                                     bet=lip_right, fa=dfa, fb=dfb, under=False)
                    right_end = upper_estimator.get_right_end_upper_bound()
                data.sub_interval.x[0] = left_end
                data.sub_interval.x[1] = right_end
                return False
            else:
                return True

    def fzcp_process(self, data: ProcData):
        """
        Process of branching
        """
        sub_interval = data.sub_interval
        # print(sub_interval.x)
        lst = []
        obj = self.problem.objective
        width_of_interval = sub_interval.x[1] - sub_interval.x[0]
        # in_bounds = self.update_interval(sub_interval)

        if not self.global_lipint:
            if self.adaptive:
                if data.counter < data.period_comp_lip:
                    data.counter += 1
                else:
                    self.update_lipschitz(data)
                    data.period_comp_lip = int(512 / math.sqrt(data.lip.x[1] - data.lip.x[0]))
                    print(data.lip.x[1] - data.lip[0], data.period_comp_lip)
                    data.counter = 0
            else:
                self.update_lipschitz(data)
        dfa = self.problem.df(data.sub_interval.x[0])
        dfb = self.problem.df(data.sub_interval.x[1])
        if self.reduction == 2 and self.reduction2(data=data, dfa=dfa, dfb=dfb):
            print('reduced!')
            return []
        lower_estimator = self.compute_bounds(data=data, under=True, dfa=dfa, dfb=dfb)
        (split_point, bound_y) = lower_estimator.lower_bound_and_point()
        if bound_y > -self.eps:
            return []
        f_split = obj(split_point)
        if f_split < self.rec_v:
            self.rec_v = f_split
            self.rec_x = split_point

        if self.reduction >= 1:
            right_end = lower_estimator.get_right_end_under_bound()
            left_end = lower_estimator.get_left_end()

            new_width = right_end - left_end
            if new_width / width_of_interval > 0.7:
                sub_1 = ival.Interval([left_end, split_point])
                data2 = ProcData(sub_interval=ival.Interval([split_point, right_end]),
                                 lip=copy.deepcopy(data.lip), counter=data.counter,
                                 period_comp_lip=data.period_comp_lip)
                lst.append(data2)
                data.sub_interval = sub_1
                lst.append(data)
            else:
                data.sub_interval = ival.Interval([left_end, right_end])
                lst.append(data)
        else:
            sub_1 = ival.Interval([sub_interval.x[0], split_point])
            data2 = ProcData(sub_interval=ival.Interval([split_point, sub_interval.x[1]]),
                             lip=copy.deepcopy(data.lip), counter=data.counter, period_comp_lip=data.period_comp_lip)
            lst.append(data2)
            data.sub_interval = sub_1
            lst.append(data)
        return lst
