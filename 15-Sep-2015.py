
# coding: utf-8

# With this code I capture the first 5k stages of the simulation for both resources and pop-size in the sexual model with the 'pop' condition, i.e. with resources and population size depending on each other. The run I am capturing is pop5 from 17-Jun-2015   
# The output file is ready to be used as an input file for R, using [this](https://github.com/dvalenzano/R-sessions/blob/10-Jul-2015.R/15-Jul-2015.R) script

# In[68]:

import numpy as np
from scipy import interpolate
#import matplotlib.pyplot as plt
import cPickle
from time import sleep


# In[69]:

class cp(object):
    def __init__(self, inp):
        self.inp = inp
        self.o = open(self.inp, 'rb')
        self.pop_in = cPickle.load(self.o)
        self.n_stage = len(self.pop_in)-1
        self.res_in = cPickle.load(self.o)
        self.age_distr_in = cPickle.load(self.o)
        #print np.shape(age_distr_in[0])                                                                                                                                                                                         
        self.repr_rate_in = cPickle.load(self.o)
        self.repr_rate_sd_in = cPickle.load(self.o)
        for i in range(len(self.repr_rate_sd_in)):
            self.repr_rate_sd_in[i] = np.array(self.repr_rate_sd_in[i])
            self.repr_rate_sd_in[i].shape = (71,1)
        self.repr_rate_junk_in = cPickle.load(self.o)
        self.surv_rate_in = cPickle.load(self.o)
        self.surv_rate_sd_in = cPickle.load(self.o)
        for i in range(len(self.surv_rate_sd_in)):
            self.surv_rate_sd_in[i] = np.array(self.surv_rate_sd_in[i])
            self.surv_rate_sd_in[i].shape = (71,1)
        self.surv_rate_junk_in = cPickle.load(self.o)
        self.repr_fit_in = cPickle.load(self.o)
        self.repr_fit_junk_in = cPickle.load(self.o)
        self.surv_fit_in = cPickle.load(self.o)
        self.surv_fit_junk_in = cPickle.load(self.o)
        self.fit_in = np.array(self.repr_fit_in)*np.array(self.surv_fit_in)
        self.fit_junk_in = np.array(self.repr_fit_junk_in)*np.array(self.surv_fit_junk_in)
        self.dens_surv_in = cPickle.load(self.o)
        self.dens_repr_in = cPickle.load(self.o)
        self.hetrz_mea = cPickle.load(self.o)
        self.hetrz_mea_sd = cPickle.load(self.o) # when simul version > 0.6                                                                                                                                                                
        self.males_females_ages = cPickle.load(self.o)
        self.o.close()
        
    # actual survival series                                                                                                                                                                                                 
    def compute_actual_surv_rate(self, p):
        #self.p = p
        """                                                                                                                                                                                                                  
        Takes age distribution of two consecutive stages and computes the                                                                                                                                                    
        fractions of those survived from age x to age x+1. The cumulative product                                                                                                                                            
        of those values builds the final result.                                                                                                                                                                             
        Returns a numpy array.                                                                                                                                                                                               
        """
        div = self.age_distr_in[p]*self.pop_in[p]
        div[div == 0] = 1
        stage2 = np.array(list((self.age_distr_in[p+1]*self.pop_in[p+1]))[1:]+[0])

        res = stage2 / div
        for i in range(1,len(res)):
            res[i] = res[i-1] * res[i]
        return res        

    def avr_actual_surv_rate(self, m):
        """Averages actual survival rate over 100 stages."""
        if m <= 50:
            res2 = self.compute_actual_surv_rate(m+100)
            for i in range(m,m+100):
                res2 += self.compute_actual_surv_rate(i)
            return res2/100
        if m >= self.n_stage-50:
            res2 = self.compute_actual_surv_rate(self.n_stage-101)
            for i in range(self.n_stage-100,self.n_stage-1):
                res2 += self.compute_actual_surv_rate(i)
            return res2/100
        else:
            res2 = self.compute_actual_surv_rate(m+50)
            for i in range(m-50,m+50):
                res2 += self.compute_actual_surv_rate(i)
            return res2/100


# In[70]:

class het(object):
    
    """ This class allows to averaging 10 consecutive intervals to plot the genome evolution scatterplot """
    
    def __init__(self, inp):
        self.inp = inp
        self.rng = range(len(self.inp[0]))[::10]

    def loop1(self, ind):
        self.ind = int(ind)
        return [np.average(self.inp[self.ind][i:i+10]) for i in self.rng]
    
    def loop2(self):
        return [self.loop1(i) for i in range(len(self.inp))]


# In[71]:

s = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/plot_values_run1.txt'


# In[72]:

spop = cp(s)
sp = spop.pop_in
sr = spop.res_in


# In[79]:

pop_5k = sp[:5000]
res_5k = sr[:5000]
time = range(1,5001)


# In[80]:

p = ['population']+[str(i) for i in pop_5k]
r = ['resources']+[str(i) for i in res_5k]
t = ['time']+[str(i) for i in time]

l = [p,r,t]
lt = zip(*l)
#ltl = [list(i) for i in lt]
ltl = ','.join([','.join(list(i))+'\n' for i in lt]).replace('\n,','\n')
z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/popres2.csv','w')
z.write(ltl)
z.close()


# Below I analyse the pop run analysis to generate the figure for the genome evolution - started on 15-Sep-2015

# In[83]:

s1 = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/first-run/plot_values_run1.txt'
s2 = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/second-run/plot_values_run1.txt'
sz = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/last-run/plot_values_run1.txt'

s1_c = cp(s1)
s2_c = cp(s2)
sz_c = cp(sz)


# In[93]:

hcs1 = het(s1_c.hetrz_mea)
red1 = hcs1.loop2()
hcsz = het(sz_c.hetrz_mea)
redz = hcsz.loop2()


# In the following block I generate the het progression for stage 0, 1.6k, 5k, 10k and 60k

# In[96]:

a = ','.join([str(i) for i in red1[0]]).replace(',',',0\n')+',0\n'
b = ','.join([str(i) for i in red1[2]]).replace(',',',1\n')+',1\n'
c = ','.join([str(i) for i in red1[5]]).replace(',',',3\n')+',3\n'
d = ','.join([str(i) for i in red1[11]]).replace(',',',6\n')+',6\n'
z = ','.join([str(i) for i in redz[-9]]).replace(',',',z\n')+',z\n'
abcdz = 'het,group\n'+a+b+c+d+z


# In[97]:

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/het-freq.csv', 'w')
z.write(abcdz)
z.close()


# In[ ]:



