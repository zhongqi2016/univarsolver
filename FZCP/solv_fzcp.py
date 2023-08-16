import sys
# from sortedcontainers import SortedKeyList
from collections import namedtuple
import psqeprocessor_fzcp as psqproc
import pslprocessor_fzcp as pslproc
import processor_new as sergproc
import processor_Casado as casproc

sys.path.append("..")
import interval as ival
import bnb as bnb
import sub as sub

TestResult = namedtuple('TestResult', ['nsteps', 'first_crossing_zero_point'])


def get_initial_recval(prob, known_record):
    if known_record:
        return prob.min_f
    else:
        return float('inf')


def psl(prob, sym=True, max_steps=sys.maxsize, epsilon=1e-2, global_lipschitz_interval=True, known_record=False):
    psp = pslproc.PSLProcessor(rec_v=get_initial_recval(prob, known_record), rec_x=prob.b, problem=prob, eps=epsilon,
                               global_lipint=global_lipschitz_interval, use_symm_lipint=sym)
    sl = []
    subp = sub.Sub(0, [0, 0], pslproc.PSLData(ival.Interval([prob.a, prob.b]), 0))
    psp.updateSplitAndBounds(subp)
    sl.append(subp)
    cnt = max_steps
    steps = bnb.bnb_fzcp(sl, cnt, psp)
    return TestResult(nsteps=steps, first_crossing_zero_point=psp.res_list[0])


def psqe(prob, sym=True, max_steps=sys.maxsize, epsilon=1e-2, global_lipschitz_interval=True, known_record=False):
    psp = psqproc.PSQEProcessor_FZCP(rec_v=get_initial_recval(prob, known_record), rec_x=prob.b, problem=prob,
                                     eps=epsilon,
                                     global_lipint=global_lipschitz_interval, use_symm_lipint=sym)
    sl = []
    subp = sub.Sub(0, [0, 0], psqproc.PSQEData(ival.Interval([prob.a, prob.b]), 0))
    psp.updateSplitAndBounds(subp)
    sl.append(subp)
    cnt = max_steps
    steps = bnb.bnb_fzcp(sl, cnt, psp)
    return TestResult(nsteps=steps, first_crossing_zero_point=psp.res_list[0])


def cas(prob, sym=True, max_steps=sys.maxsize, epsilon=1e-2, known_record=False):
    psp = casproc.CasProcessor(rec_v=get_initial_recval(prob, known_record), rec_x=prob.b, problem=prob,
                               eps=epsilon)
    sl = []
    subp = sub.Sub(0, [0, 0], psqproc.PSQEData(ival.Interval([prob.a, prob.b]), 0))
    psp.updateSplitAndBounds(subp)
    sl.append(subp)
    cnt = max_steps
    steps = bnb.bnb_fzcp(sl, cnt, psp)
    return TestResult(nsteps=steps, first_crossing_zero_point=psp.res_list[0])


def method2(prob, sym=True, max_steps=sys.maxsize, epsilon=1e-2, global_lipschitz_interval=True, known_record=False):
    psp = psqproc.PSQEProcessor_FZCP(rec_v=get_initial_recval(prob, known_record), rec_x=prob.b, problem=prob,
                                     eps=epsilon,
                                     global_lipint=global_lipschitz_interval, use_symm_lipint=sym)
    # sl = SortedKeyList(key=lambda s: s.level)
    subp = sub.Sub(0, [0, 0], psqproc.PSQEData(ival.Interval([prob.a, prob.b]), 0))
    steps = 0
    obj = psp.problem.objective
    while steps <= max_steps:
        psp.update_interval(subp)
        steps = steps + 1
        if subp.data.ival.x[1] - subp.data.ival.x[0] < epsilon:
            break
    return TestResult(nsteps=steps, first_crossing_zero_point=subp.data.ival.x[0])


def new_method(prob, symm=True, max_steps=sys.maxsize, epsilon=1e-2, global_lipschitz_interval=False,
               known_record=False, estimator=2, reduction=1, ):
    psp = sergproc.ProcessorNew(rec_v=get_initial_recval(prob, known_record), rec_x=prob.b, problem=prob,
                                eps=epsilon,
                                global_lipint=global_lipschitz_interval, use_symm_lipint=symm, estimator=estimator,
                                reduction=reduction)
    sl = []
    interval = ival.Interval([prob.a, prob.b])
    # psp.update_interval(interval)
    sl.append(interval)
    cnt = max_steps
    steps = bnb.bnb_fzcp(sl, cnt, psp)
    return TestResult(nsteps=steps, first_crossing_zero_point=psp.res_list)
