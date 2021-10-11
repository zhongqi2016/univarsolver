import math
import sympy
import interval as ival

from sortedcontainers import SortedList
from sortedcontainers import SortedKeyList

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
# print ("Hello")
