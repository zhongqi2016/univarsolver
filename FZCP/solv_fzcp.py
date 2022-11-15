import sys
from sortedcontainers import SortedKeyList
from collections import namedtuple
import psqeprocessor_fzcp as psqproc
import pslprocessor_fzcp as pslproc

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
    sl = SortedKeyList(key=lambda s: s.level)
    subp = sub.Sub(0, [0, 0], pslproc.PSLData(ival.Interval([prob.a, prob.b]), 0))
    psp.updateSplitAndBounds(subp)
    sl.add(subp)
    cnt = max_steps
    steps = bnb.bnb_fzcp(sl, cnt, psp)
    return TestResult(nsteps=steps, first_crossing_zero_point=psp.res_list[0])


def psqe(prob, sym=True, max_steps=sys.maxsize, epsilon=1e-2, global_lipschitz_interval=True, known_record=False):
    psp = psqproc.PSQEProcessor_FZCP(rec_v=get_initial_recval(prob, known_record), rec_x=prob.b, problem=prob,
                                     eps=epsilon,
                                     global_lipint=global_lipschitz_interval, use_symm_lipint=sym)
    sl = SortedKeyList(key=lambda s: s.level)
    subp = sub.Sub(0, [0, 0], psqproc.PSQEData(ival.Interval([prob.a, prob.b]), 0))
    psp.updateSplitAndBounds2(subp)
    sl.add(subp)
    cnt = max_steps
    steps = bnb.bnb_fzcp(sl, cnt, psp)
    return TestResult(nsteps=steps, first_crossing_zero_point=psp.res_list[0])
