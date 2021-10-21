import math
import sympy as sym
import pandas as pd
import interval as ival
import  uvarprob as uvpr

from sortedcontainers import SortedList
from sortedcontainers import SortedKeyList


def read_problems(fname):
    data = pd.read_csv(fname, index_col = 'name', comment = '#')
    return data


def make_objective(func):
    sf = sym.sympify(func)
    print(sf)
    x = sym.symbols('x')
    f = sym.lambdify(x, sf)
    return f



class Sub:
    def __init__(self, level, bound, data):
        self.level = level
        self.bound = bound
        self.data = data

    def __repr__(self):
        return "[ level = " + str(self.level) + ", bound = " + str(self.bound) + ", data =  " + str(self.data) + "]"

I1 = ival.Interval([-1,1])
s1 = Sub(1,-1,I1)
print(s1)
I2 = ival.Interval([-1,2])
s2 = Sub(3,-1.1,I2)
subs = [s1, s2]
print (subs)
# sl = SortedKeyList(key = lambda s : -s.bound)
sl = SortedKeyList(key = lambda s : s.level)
sl.update(subs)
print(sl)
sl.add(Sub(2, 2, ival.Interval([-2,1])))
print(sl)
while len(sl) > 0:
    print(sl.pop())
problems = read_problems("problems.csv")
print(problems)
obj = make_objective(problems.loc['prob1','objective'])
print(obj(1))
prob = uvpr.UniVarProblem('prob1', problems.loc['prob1','objective'], problems.loc['prob1','a'], problems.loc['prob1','b'])
print(prob)
print(prob.objective(prob.b))
problems.to_csv('/tmp/prob.csv')
# print(problems.objective)
# print(problems.loc[lambda df : df['name'] == 'prob1'])
# df1 = problems.loc[lambda df : df['name'] == 'prob1']
# print(df1.loc[0,'objective'])
# make_objective(df1.loc[0,'objective'])
# print ("Hello")
