
# coding: utf-8

# In[7]:

import numpy as np
from scipy import interpolate
#import matplotlib.pyplot as plt
import cPickle
from time import sleep


# In[8]:

class cp2(object):
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
                
        
    def compute_actual_death_rate(self, s):
        """                                                                                                                                              
        Takes age distribution of two consecutive stages. Computes the fraction                                                                          
        of those died from age x to age x+1.                                                                                                             
        Returns a numpy array.                                                                                                                           
        """
        stage1 = self.age_distr_in[s]*self.pop_in[s]
        stage2 = np.array(list((self.age_distr_in[s+1]*self.pop_in[s+1]))[1:]+[0])
        div = stage1
        div[div == 0] = 1

        return (stage1 - stage2) / div
    
    
   
    def avr_actual_death_rate(self, s):
        """Averages actual death rate over 100 stages."""
        if s <= 50:
            res = self.compute_actual_death_rate(s)
            for i in range(s+1,s+101):
                res += self.compute_actual_death_rate(i)
            return res/100
        if s >= self.n_stage-120:
            res = self.compute_actual_death_rate(self.n_stage-70)
            for i in range(self.n_stage-170,self.n_stage-70):
                res += self.compute_actual_death_rate(i)
            return res/100
        else:
            res = self.compute_actual_death_rate(s+50)
            for i in range(s-50,s+50):
                res += self.compute_actual_death_rate(i)
            return res/100


# In[9]:

lc_s = '/Volumes/group_dv/personal/DValenzano/month-by-month/Jun2015/simul/sex/17-Jun-2015/plot_values_run1.txt'
lp_s = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/last-run/plot_values_run1.txt'
sc_s = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/small-popsize/const-last.txt'
sp_s = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/small-popsize/pop-last.txt'

lc_as = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/lcl/plot_values_run1.txt'
lp_as = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/lpl/plot_values_run1.txt'
sc_as = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/scl/plot_values_run1.txt'
sp_as = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/spl/plot_values_run1.txt'


# In[10]:

Clc_s = cp2(lc_s)
Clp_s = cp2(lp_s)
Csc_s = cp2(sc_s)
Csp_s = cp2(sp_s)
Clc_as = cp2(lc_as)
Clp_as = cp2(lp_as)
Csc_as = cp2(sc_as)
Csp_as = cp2(sp_as)


# In[69]:

Clc_sdr = 'Clc_s,'+','.join([str(i) for i in list(Clc_s.avr_actual_death_rate(Clc_s.n_stage))])
Clp_sdr = 'Clp_s,'+','.join([str(i) for i in list(Clp_s.avr_actual_death_rate(Clp_s.n_stage))])
Csc_sdr = 'Csc_s,'+','.join([str(i) for i in list(Csc_s.avr_actual_death_rate(Csc_s.n_stage))])
Csp_sdr = 'Csp_s,'+','.join([str(i) for i in list(Csp_s.avr_actual_death_rate(Csp_s.n_stage))])
Clc_asdr = 'Clc_as,'+','.join([str(i) for i in list(Clc_as.avr_actual_death_rate(Clc_as.n_stage))])
Clp_asdr = 'Clp_as,'+','.join([str(i) for i in list(Clp_as.avr_actual_death_rate(Clp_as.n_stage))])
Csc_asdr = 'Csc_as,'+','.join([str(i) for i in list(Csc_as.avr_actual_death_rate(Csc_as.n_stage))])
Csp_asdr = 'Csp_as,'+','.join([str(i) for i in list(Csp_as.avr_actual_death_rate(Csp_as.n_stage))])
xdr = 'Age,'+','.join([str(i) for i in range(71)])

drinput = zip(xdr.split(','), Clc_sdr.split(','), Clp_sdr.split(','), Csc_sdr.split(','), Csp_sdr.split(','), Clc_asdr.split(','), Clp_asdr.split(','), Csc_asdr.split(','), Csp_asdr.split(','))
drinput1 = ','.join([','.join(list(i))+'\n' for i in drinput ]).replace('\n,','\n')


# In[70]:

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/death-rate_rinput.csv', 'w')
z.write(drinput1)
z.close()


# In[ ]:



