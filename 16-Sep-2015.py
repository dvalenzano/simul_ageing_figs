
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


# In the block below I generate the sd progression for the sex-pop evol

# In[100]:

sdcsz = het(sz_c.hetrz_mea_sd)
sd_red = sdcsz.loop2()
scz_sd_60k = sd_red[-9]

sd_cs1 = het(s1_c.hetrz_mea_sd)
sd_red1 = sd_cs1.loop2()

sd_a = ','.join([str(i) for i in sd_red1[0]]).replace(',',',0\n')+',0\n'
sd_b = ','.join([str(i) for i in sd_red1[2]]).replace(',',',1\n')+',1\n'
sd_c = ','.join([str(i) for i in sd_red1[5]]).replace(',',',3\n')+',3\n'
sd_d = ','.join([str(i) for i in sd_red1[11]]).replace(',',',6\n')+',6\n'
sd_z = ','.join([str(i) for i in scz_sd_60k]).replace(',',',z\n')+',z\n'
sd_abcdz = 'het,group\n'+sd_a+sd_b+sd_c+sd_d+sd_z

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/het-sd.csv', 'w')
z.write(sd_abcdz)
z.close()


# Here below I generate a figure where I compare - in the sexual model - const vs pop with low and high resources. 

# In[171]:

large_pop_early = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/first-run/plot_values_run1.txt'
large_pop_late = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/last-run/plot_values_run1.txt'
large_const_early = '/Volumes/group_dv/personal/DValenzano/papers/simulation_arXiv/Figure3/first-run/plot_values_run1.txt'
large_const_late = '/Volumes/group_dv/personal/DValenzano/month-by-month/Jun2015/simul/sex/17-Jun-2015/plot_values_run1.txt'

small_pop_early = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/small-popsize/pop-first.txt'
small_pop_late = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/small-popsize/pop-last.txt'
small_const_early = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/small-popsize/const-first.txt'
small_const_late = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/small-popsize/const-last.txt'


# In[172]:

lpe = cp(large_pop_early)
lpl = cp(large_pop_late)
lce = cp(large_const_early)
lcl = cp(large_const_late)
spe = cp(small_pop_early)
spl = cp(small_pop_late)    
sce = cp(small_const_early)
scl = cp(small_const_late)


# I need to types of plots: one where I show early and late pop-res oscillations in the small and large populations,  
# both constant and pop

# In[175]:

# This gives the arrays that enable me to compute Si and Ri values in large and small early and late stages for 
# "const" and "pop" runs

hlpe = het(lpe.hetrz_mea)
red_lpe = hlpe.loop2()

hlpl = het(lpl.hetrz_mea)
red_lpl = hlpl.loop2()

hlce = het(lce.hetrz_mea)
red_lce = hlce.loop2()

hlcl = het(lcl.hetrz_mea)
red_lcl = hlcl.loop2()

hspe = het(spe.hetrz_mea)
red_spe = hspe.loop2()

hspl = het(spl.hetrz_mea)
red_spl = hspl.loop2()

hsce = het(sce.hetrz_mea)
red_sce = hsce.loop2()

hscl = het(scl.hetrz_mea)
red_scl = hscl.loop2()


# In[158]:

# This gives the arrays that enable me to compute pop and resources values in large and small early and late 
# stages for "const" and "pop" runs

lpe_pop = lpe.pop_in
lpe_res = lpe.res_in

lpl_pop = lpl.pop_in
lpl_res = lpl.res_in

lce_pop = lce.pop_in
lce_res = lce.res_in

lcl_pop = lcl.pop_in
lcl_res = lcl.res_in

spe_pop = spe.pop_in
spe_res = spe.res_in

spl_pop = spl.pop_in
spl_res = spl.res_in

sce_pop = sce.pop_in
sce_res = sce.res_in

scl_pop = scl.pop_in
scl_res = scl.res_in


# First, I do the population-resources plot

# In[132]:

# For large populations - sex-pop:
lpe_p = lpe_pop[:6000]
lpe_r = lpe_res[:6000]
lpl_p = lpl_pop[9000:15000]
lpl_r = lpl_res[9000:15000]
lp_time = range(1,6001)+range(54001,60000)

lp_p = ['population']+[str(i) for i in lpe_p]+[str(i) for i in lpl_p]
lp_r = ['resources']+[str(i) for i in lpe_r]+[str(i) for i in lpl_r]
lp_t = ['time']+[str(i) for i in lp_time]

lp = [lp_p,lp_r,lp_t]
lptr = zip(*lp)
#ltl = [list(i) for i in lt]
lptr_l = ','.join([','.join(list(i))+'\n' for i in lptr]).replace('\n,','\n')


# In[133]:

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/large-sex-popres.csv','w')
z.write(lptr_l)
z.close()


# In[163]:

# For large populations - sex-const:
lce_p = lce_pop[:6000]
lce_r = lce_res[:6000]
lcl_p = lcl_pop[-6000:]
lcl_r = lcl_res[-6000:]
lc_time = range(1,6001)+range(54001,60001)

lc_p = ['population']+[str(i) for i in lce_p]+[str(i) for i in lcl_p]
lc_r = ['resources']+[str(i) for i in lce_r]+[str(i) for i in lcl_r]
lc_t = ['time']+[str(i) for i in lc_time]

lc = [lc_p,lc_r,lc_t]
lctr = zip(*lc)
lctr_l = ','.join([','.join(list(i))+'\n' for i in lctr]).replace('\n,','\n')

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/large-sex-constres.csv','w')
z.write(lctr_l)
z.close()


# In[155]:

# For small populations - sex-pop:
spe_p = spe_pop[:6000]
spe_r = spe_res[:6000]
spl_p = spl_pop[-6000:]
spl_r = spl_res[-6000:]
sp_time = range(1,6001)+range(54001,60001)

sp_p = ['population']+[str(i) for i in spe_p]+[str(i) for i in spl_p]
sp_r = ['resources']+[str(i) for i in spe_r]+[str(i) for i in spl_r]
sp_t = ['time']+[str(i) for i in sp_time]

sp = [sp_p,sp_r,sp_t]
sptr = zip(*sp)
#ltl = [list(i) for i in lt]
sptr_l = ','.join([','.join(list(i))+'\n' for i in sptr]).replace('\n,','\n')

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/small-sex-popres.csv','w')
z.write(sptr_l)
z.close()


# In[151]:

# For small populations - const-pop:
sce_p = sce_pop[:6000]
sce_r = sce_res[:6000]
scl_p = scl_pop[-6000:]
scl_r = scl_res[-6000:]
sc_time = range(1,6001)+range(54001,60001)

sc_p = ['population']+[str(i) for i in sce_p]+[str(i) for i in scl_p]
sc_r = ['resources']+[str(i) for i in sce_r]+[str(i) for i in scl_r]
sc_t = ['time']+[str(i) for i in sc_time]

sc = [sc_p,sc_r,sc_t]
sctr = zip(*sc)
#ltl = [list(i) for i in lt]
sctr_l = ','.join([','.join(list(i))+'\n' for i in sctr]).replace('\n,','\n')

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/small-sex-constres.csv','w')
z.write(sctr_l)
z.close()


# 16-Sep-2015 is from here  
# Below, I extract, for the 60k's stage, and compare the trend of the surv-repr plot in large and small populations

# In[176]:

lps_z = red_lpl[-9]
lcs_z = red_lcl[-1]
sps_z = red_spl[-1]
scs_z = red_scl[-1]


# In[177]:

a = ','.join([str(i) for i in lps_z]).replace(',',',lps\n')+',lps\n'
b = ','.join([str(i) for i in lcs_z]).replace(',',',lcs\n')+',lcs\n'
c = ','.join([str(i) for i in sps_z]).replace(',',',sps\n')+',sps\n'
d = ','.join([str(i) for i in scs_z]).replace(',',',scs\n')+',scs\n'
abcd = 'surv,group\n'+a+b+c+d

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/pop/small-large_surv.csv', 'w')
z.write(abcd)
z.close()


# In[ ]:



