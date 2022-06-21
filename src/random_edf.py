import random
import bisect


class Edfreq(object):

        def __init__(self, histogram, name=''):
                """Makes a edf from an unsorted sequence of (value, frequency) pairs Dictionary type."""

                self.name = name
                self.items = histogram.iteritems()
                runsum = 0
                xs = []
                cs = []

                for value, count in sorted(self.items):
                        runsum += count
                        xs.append(value)
                        cs.append(runsum)

                total = float(runsum)
                ps = [c/total for c in cs]

    
                self.xs = [] if xs is None else xs
                self.ps = [] if ps is None else ps
        

        def prob(self, x):
                """Returns edf(x), the probability that corresponds to value x.

                Args:
            x: number

                Returns:
            float probability
                """
                if x < self.xs[0]: return 0.0
                index = bisect.bisect(self.xs, x)
                p = self.ps[index-1]
                return p

        def value(self, p):
                """Returns Inverseedf(p), the value that corresponds to probability p.

        Args:
            p: number in the range [0, 1]

        Returns:
            number value
                """
                if p < 0 or p > 1:
                        raise ValueError('Probability p must be in range [0, 1]')

                if p == 0: return self.xs[0]
                if p == 1: return self.xs[-1]
                index = bisect.bisect(self.ps, p)
                if p == self.ps[index-1]:
                        return self.xs[index-1]
                else:
                        return self.xs[index]

        def render(self,filename):
            """Generates a sequence of points suitable for plotting.

        An empirical edf is a step function; linear interpolation
        can be misleading."""

            content = '0'+'\t'+'0'+'\n'
            add2file(filename,content)
            
            top = self.value(1)
            
            for a in range(0,top,2):
                content = str(a)+'\t'+str(self.prob(a))+'\n'
                add2file(filename,content)


            content = str(top)+'\t'+str(self.prob(top))+'\n'
            add2file(filename,content)

            
            return


def add2file(filename, content):
    
    filename ='%s' % filename
    f = open(filename, 'a')
    f.write( content )
    f.close()
    
    return

numbers = []
for i in range(0,9000):
    numbers.append(random.randrange(2,200))
for i in range(0,1000):
    numbers.append(random.randrange(200,500))
for i in range(0,400):
    numbers.append(random.randrange(600,1000))


histogram = {}
for number in numbers:
    histogram[number] = histogram.get(number, 0) + 1

edf = Edfreq(histogram,'random')
edf.render('random.dat')
