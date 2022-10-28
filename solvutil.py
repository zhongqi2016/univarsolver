import pslprocessor as pslproc
import psqeprocessor as psqproc
import interval as ival
import bnb as bnb
import sub as sub
from sortedcontainers import SortedKeyList
from collections import namedtuple
import sys

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
