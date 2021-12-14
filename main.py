import math
import sympy as sym
import pandas as pd
import numpy as np
import interval as ival
import uvarprob as uvpr
import bnb as bnb
import sub as sub
import ivalprocessor as ivproc
import psqeprocessor as psproc
import gridsearch as gs

# from sortedcontainers import SortedList
from sortedcontainers import SortedKeyList


def read_problems(fname):
    data = pd.read_csv(fname, index_col = 'name', comment = '#')
    return data

# problems = read_problems("problems.csv")
problems = read_problems("/tmp/shek.csv")
print(problems)
name = 'rshekel_3'
prob = uvpr.UniVarProblem(name, problems.loc[name,'objective'], problems.loc[name,'a'], problems.loc[name,'b'], problems.loc[name,'min_f'], problems.loc[name,'min_x'])
print("Parsed objective = ", problems.loc[name,'objective'])
print("Parsed problem = ", prob)

dp = ivproc.IntervalProcessor(np.nan, np.nan, prob, 1e-5)
sl = SortedKeyList(key = lambda s : s.level)
subp = sub.Sub(0, 0, ival.Interval([prob.a, prob.b]))
dp.compute_bounds(subp)
sl.add(subp)
print(sl)
cnt = 100000
print("steps performed: ", bnb.bnb(sl, cnt, dp))
print("Record value = ", dp.rec_v, " at ", dp.rec_x);
print(sl)

psp = psproc.PSQEProcessor(float('inf'), None, prob, 1e-5, False, True)
sl = SortedKeyList(key = lambda s : s.level)
subp = sub.Sub(0, 0, psproc.PSQEData(ival.Interval([prob.a, prob.b]),0))
psp.compute_bounds(subp)
sl.add(subp)
# print(sl)
cnt = 100000
print("steps performed: ", bnb.bnb(sl, cnt, psp))
print("Record value = ", psp.rec_v, " at ", psp.rec_x);
print(sl)

true_min = gs.grid_search(prob.objective, prob.a, prob.b, 1e-4)
print("Grid search:", true_min)

# while len(sl) > 0:
#     print(sl.pop())
