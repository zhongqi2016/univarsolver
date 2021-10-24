#TODO: transfer sub to another file
class Sub:
    """ Generic subproblem

    Class for storing information regarding the subproblem.

    Attributes:
        level(int): level in the BnB-tree
        bound(double): the bound (can be upper or lower depending on the problem)
        data(any): the problem-specific data associated with the subproblem
    """

    def __init__(self, level, bound, data):
        """
        Initializes the subproblem.

        Args:
            level: level in the BnB-tree
            bound: the bound (can be upper or lower depending on the problem)
            data: the problem-specific data associated with the subproblem
        """
        self.level = level
        self.bound = bound
        self.data = data

    def __repr__(self):
        return "[ level = " + str(self.level) + ", bound = " + str(self.bound) + ", data =  " + str(self.data) + "]"


def bnb(subs, max_steps, processor):
    """
    Main Branch-and-Bound driver

    Args:
        subs: the list of subproblems
        max_steps: mximal number of steps to perform
        processor: the processor

    Returns:
        number of actually performed steps
    """
    steps = 0
    while len(subs) > 0 and steps <= max_steps:
        s = subs.pop()
        new_subs = processor.process(s)
        subs.update(new_subs)
        steps = steps + 1
    return steps
