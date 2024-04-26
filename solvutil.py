import math
import pslprocessor as pslproc
import psqeprocessor as psqproc
import interval as ival
import bnb as bnb
import sub as sub
from sortedcontainers import SortedKeyList
from collections import namedtuple
import sys
import processor_reduction

TestResult = namedtuple('TestResult', ['nsteps', 'record_value'])


def get_initial_recval(prob, known_record):
    if known_record:
        return prob.min_f
    else:
        return float('inf')


def psl(prob, sym=True, max_steps=sys.maxsize, epsilon=1e-2, global_lipschitz_interval=True, known_record=False):
    """
    Runs quadratic minorant solver
    :param prob: problem to solver
    :param sym: if use symmetric Lipschitz constant (default True)
    :param max_steps: maximal number of steps (default infinity)
    :param epsilon: tolerance (default 1e-2)
    :param global_lipschitz_interval: if True - use global Lipshitz constant or interval computed for the initial interval (default True)
    :param known_record: if True, use the known optimum as a record (default False)
    :return: the testing results 
    :rtype: TestResult
    """
    psp = pslproc.PSLProcessor(rec_v=get_initial_recval(prob, known_record), rec_x=None, problem=prob, eps=epsilon,
                               global_lipint=global_lipschitz_interval, use_symm_lipint=sym)
    sl = SortedKeyList(key=lambda s: s.level)
    subp = sub.Sub(0, 0, pslproc.PSLData(ival.Interval([prob.a, prob.b]), 0))
    psp.compute_bounds(subp)
    sl.add(subp)
    cnt = max_steps
    steps = bnb.bnb(sl, cnt, psp)
    return TestResult(nsteps=steps, record_value=psp.rec_v)


def psqe(prob, sym=True, max_steps=sys.maxsize, epsilon=1e-2, global_lipschitz_interval=True, known_record=False):
    """
    Runs quadratic minorant solver
    :param prob: problem to solver
    :param sym: if use symmetric Lipschitz constant (default True)
    :param max_steps: maximal number of steps (default infinity)
    :param epsilon: tolerance (default 1e-2)
    :param global_lipschitz_interval: if True - use global Lipshitz constant or interval computed for the initial interval (default True)
    :param known_record: if True, use the known optimum as a record (default False)
    :return: the testing results 
    :rtype: TestResult
    """
    psp = psqproc.PSQEProcessor(rec_v=get_initial_recval(prob, known_record), rec_x=None, problem=prob, eps=epsilon,
                                global_lipint=global_lipschitz_interval, use_symm_lipint=sym)
    if isinstance(psp.ddi, int):
        min_x = -psp.problem.df(0) / psp.ddi
        if psp.problem.a <= min_x <= psp.problem.b:
            psp.rec_x = min_x

        else:
            if psp.problem.objective(psp.problem.a) < psp.problem.objective(psp.problem.b):
                psp.rec_x = psp.problem.a
            else:
                psp.rec_x = psp.problem.b
        psp.rec_v = psp.problem.objective(psp.rec_x)
        return TestResult(nsteps=1, record_value=psp.rec_v)
    sl = SortedKeyList(key=lambda s: s.level)
    subp = sub.Sub(0, 0, psqproc.PSQEData(ival.Interval([prob.a, prob.b]), 0))
    psp.compute_bounds(subp)
    sl.add(subp)
    cnt = max_steps
    steps = bnb.bnb(sl, cnt, psp)
    return TestResult(nsteps=steps, record_value=psp.rec_v)


def new_proc(prob, symm=True, max_steps=sys.maxsize, epsilon=1e-2, global_lipschitz_interval=False,
             known_record=False, estimator=2, reduction=1, adaptive=False):
    fa = prob.objective(prob.a)
    fb = prob.objective(prob.b)
    if fa < fb:
        rec_x = prob.a
        rec_v = fa
    else:
        rec_x = prob.b
        rec_v = fb

    psp = processor_reduction.ProcessorReduction(rec_v=rec_v, rec_x=rec_x,
                                                 problem=prob,
                                                 eps=epsilon, global_lipint=global_lipschitz_interval,
                                                 use_symm_lipint=symm,
                                                 estimator=estimator, reduction=reduction, adaptive=adaptive)
    sl = []
    interval = ival.Interval([prob.a, prob.b])
    data = processor_reduction.ProcData(sub_interval=interval, lip=ival.Interval([0, 0]), counter=0, period_comp_lip=0)
    psp.update_lipschitz(data)
    if adaptive:
        data.period_comp_lip = int(128 / math.sqrt(data.lip.x[1] - data.lip.x[0]))
        if data.period_comp_lip == 0:
            data.period_comp_lip = 1
    print(data.lip.x[1] - data.lip.x[0], data.period_comp_lip)
    # psp.update_interval(interval)
    sl.append(data)
    cnt = max_steps
    steps = bnb.bnb_fzcp(sl, cnt, psp)
    return TestResult(nsteps=steps, record_value=psp.rec_v)
