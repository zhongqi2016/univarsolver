import math
import sympy as sym
import pandas as pd
import interval as ival
import uvarprob as uvpr
import bnb as bnb

# from sortedcontainers import SortedList
from sortedcontainers import SortedKeyList


def read_problems(fname):
    data = pd.read_csv(fname, index_col = 'name', comment = '#')
    return data



class DummyProcessor:
    def __init__(self, rec_v, rec_x, problem, eps):
        self.rec_v = rec_v
        self.rec_x = rec_x
        self.problem = problem
        self.eps = eps

    def compute_bounds(self, sub):
        sub.bound = self.problem.objective(sub.data).x[0]
        c = sub.data.mid()
        v = self.problem.objective(c)
        if v < self.rec_v:
            self.rec_x = c
            self.rec_v = v

    def process(self, sub):
        print(sub)
        sub_1 = bnb.Sub(sub.level + 1, 0, ival.Interval([sub.data.x[0], sub.data.mid()]))
        self.compute_bounds(sub_1)
        print(sub_1)
        sub_2 = bnb.Sub(sub.level + 1, 0, ival.Interval([sub.data.mid(), sub.data.x[1]]))
        self.compute_bounds(sub_2)
        print(sub_2)
        print(self.rec_v)
        lst = []
        if sub_1.bound < self.rec_v - self.eps:
            lst.append(sub_1)
        if sub_2.bound < self.rec_v - self.eps:
            lst.append(sub_2)

        return lst


problems = read_problems("problems.csv")
print(problems)
prob = uvpr.UniVarProblem('prob1', problems.loc['prob1','objective'], problems.loc['prob1','a'], problems.loc['prob1','b'], problems.loc['prob1','min_f'], problems.loc['prob1','mins_x'])
print(prob)
# problems.to_csv('/tmp/prob.csv')
dp = DummyProcessor(0, prob.objective(0), prob, 0.001)

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
