import random
import bisect
import math


numbers = []
interval = 0
histogram = {}

def add2file(filename, content):
    
    filename ='%s' % filename
    f = open(filename, 'a')
    f.write( content )
    f.close()
    
    return


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



def meanstdv(x):
    """ Calculate mean and standard deviation of data x[]:
        mean = {\sum_i x_i \over n}
        std = sqrt(\sum_i (x_i - mean)^2 \over n-1) """

    if len(x) == 1:
        return x[0],0
    
    n, mean, std = len(x), 0, 0
    
    for a in x:
        mean = mean + a
    
    mean = mean / float(n)
    
    for a in x:
        std = std + (a - mean)**2
    
    std = math.sqrt(std / float(n-1))
    
    return mean, std

def contrb(datos,prob):

	cluster = []
	centroids = [10,20,30,40,50,60,70,80,90,100,200,500,1000]

	for i in range(13):

		cluster.append([])

		#(c,db) = z_ddbb.link_db(cfg)
	#c.execute('SELECT occ_time FROM water_intervals WHERE wetspot_id=?' ,(wetspot_key,))

	for data in datos:

		#if int(data)<5:
			#cluster[0].append(int(data))
			#continue
		if int(data)<11:
			cluster[0].append(int(data))
			continue
		if int(data)<21:
			cluster[1].append(int(data))
			continue
		if int(data)<31:
			cluster[2].append(int(data))
			continue
		if int(data)<41:
			cluster[3].append(int(data))
			continue
		if int(data)<51:
			cluster[4].append(int(data))
			continue
		if int(data)<61:
			cluster[5].append(int(data))
			continue
		if int(data)<71:
			cluster[6].append(int(data))
			continue
		if int(data)<81:
			cluster[7].append(int(data))
			continue
		if int(data)<91:
			cluster[8].append(int(data))
			continue
		if int(data)<101:
			cluster[9].append(int(data))
			continue
		if int(data)<201:
			cluster[10].append(int(data))
			continue	
		if int(data)<501:
			cluster[11].append(int(data))
			continue

		else:
			cluster[12].append(int(data))
			continue

	wet_time = sum(datos)

	filename = 'random'+'_'+prob
	content= '#######    '+prob+'\n'
	add2file(filename,content)
	for i in range(13):
		try:
            
			#N = float(wet_timdata)/float(centroids[i])
			content=str(centroids[i])+'  '+str(float(sum(cluster[i]))/float(wet_time))+'\n'
			add2file(filename,content)
		except:
			content= '#  0.0'+'\n'
			add2file(filename,content)
         
	return

for i in range(0,4800000):
	if random.random() >= 0.35:
		interval += 1

	else:
		numbers.append(interval)
		interval = 0

lista = []
for e in numbers:
	lista.append(e*2)

for e in lista:
    histogram[e] = histogram.get(e,0) + 1
    edf = Edfreq(histogram,'prob')

contrb(lista,'0.35')
edf.render('random_WATprob35_edf.dat')

#print 'sum, Tmax, avg,std,porciento,numero de eventos '
avg , std = meanstdv(lista)
print sum(lista),';', max(lista),';',avg,';',std,';', float(sum(lista))/float(24000)*100,';', len(lista)


