import matplotlib.pyplot as plt
import numpy as np

default_colormap = {'function' : 'b-'}

def plot_problem(prob, colormap = default_colormap, npoints = 1000):
    """
       Plots a problem
       
       Parameters
       ----------
       prob : UniVarProblem
           The problem to plot
       npoints : int
           number of points used to draw a plot
       """

#     colors = ['r-', 'b-', 'g-', 'y-', 'm-', 'c-']
    step = (prob.b-prob.a)/npoints
    ta = np.arange(prob.a, prob.b + step, step)
    num_points = len(ta)
    fta = np.empty([num_points])
       
    for i in range(num_points):
        fta[i] = prob.objective(ta[i])
    lb = np.amin(fta)
    ub = np.amax(fta)
    d = (ub - lb) * 0.1
    plt.plot(ta, fta, colormap['function'])
    plt.ylim([lb - d,ub + d])
#     plt.show() 

    

def plot_points(x, y, colormap = default_colormap):
    plt.scatter(x, y, s = 0.01)
    
    