"""
Shekel's function generator
"""
import pandas as pd

def gen_one_problem(n, reverse = False):
    """
    Generates Schekel's or reverse Shekel's function
    Args:
        n: the number of a problem
        reverse: if True - generate reverse Shekel's function, if False - normal one

    Returns:
        data frame for a new problem

    """
    name = "shekel_" + str(n)
    dct = dict(name=name, objective="x^2 - 1", a=0., b=1., min_f=0.11, mins_x=str([1,2]))
    return pd.DataFrame(dct, index = [0])

df = None
for i in range(0,5):
    dfn = gen_one_problem(i)
    if df is None:
        df = dfn
    else:
        df = pd.concat([df, dfn])

print(df)
df.to_csv('/tmp/shek.csv', index = False)
