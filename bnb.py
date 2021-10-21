
class Sub:
    """ Generic sub-problem representation """

    def __init__(self, level, bound, data):
        self.level = level
        self.bound = bound
        self.data = data

    def __repr__(self):
        return "[ level = " + str(self.level) + ", bound = " + str(self.bound) + ", data =  " + str(self.data) + "]"


def bnb(subs, max_steps, processor):
    """ Branch-and-bound method generic driver """
    steps = 0
    while len(subs) > 0 and steps <= max_steps:
        s = subs.pop()
        new_subs = processor.process(s)
        subs.update(new_subs)
        steps = steps + 1
    return steps
