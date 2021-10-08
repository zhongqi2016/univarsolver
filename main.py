import math
import sympy
import interval as ival

from sortedcontainers import SortedList
from sortedcontainers import SortedKeyList

subs = [{"chr" : - 1.0, "dat" : 1}, {"chr" : - 1.2, "dat" : 2}]
print (subs)
sl = SortedKeyList(key = lambda x : -x["chr"])
sl.update(subs)
print(sl)
sl.add({"chr" : - 1.5, "dat" : 3})
print(sl)
while len(sl) > 0:
    print(sl.pop())
# print ("Hello")
