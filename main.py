import math
import sympy as sym
import pandas as pd
import interval as ival
import uvarprob as uvpr
import bnb as bnb
import ivalprocessor as ivproc

# from sortedcontainers import SortedList
from sortedcontainers import SortedKeyList


def read_problems(fname):
    data = pd.read_csv(fname, index_col = 'name', comment = '#')
    return data




problems = read_problems("problems.csv")
print(problems)
prob = uvpr.UniVarProblem('prob1', problems.loc['prob1','objective'], problems.loc['prob1','a'], problems.loc['prob1','b'], problems.loc['prob1','min_f'], problems.loc['prob1','mins_x'])
print(prob)
# problems.to_csv('/tmp/prob.csv')
dp = ivproc.IntervalProcessor(0, prob.objective(0), prob, 0.001)

sl = SortedKeyList(key = lambda s : s.level)
sub = bnb.Sub(0, 0, ival.Interval([prob.a, prob.b]))
dp.compute_bounds(sub)
sl.add(sub)
print(sl)
cnt = 1000
print("steps performed: ", bnb.bnb(sl, cnt, dp))
print(sl)

# while len(sl) > 0:
#     print(sl.pop())
