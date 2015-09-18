
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

# In[229]:

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


# From here it's 17-Sep-2015  
# Below, I extract the population-resources oscillations for the asexual model

# In[260]:

# First, large population (5k)

lp_e_as = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/lpe/plot_values_run1.txt'
lp_l_as = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/lpl/plot_values_run1.txt'
lc_e_as = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/lce/plot_values_run1.txt'
lc_l_as = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/lcl/plot_values_run1.txt'

sp_e_as = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/spe/plot_values_run1.txt'
sp_l_as = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/spl/plot_values_run1.txt'
sc_e_as = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/sce/plot_values_run1.txt'
sc_l_as = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/scl/plot_values_run1.txt'


# In[197]:

# This gives the arrays that enable me to compute pop and resources values in large populations, for early and late 
# stages in "const" and "pop" runs

lpe_as = cp(lp_e_as)
lpl_as = cp(lp_l_as)
lce_as = cp(lc_e_as)
lcl_as = cp(lc_l_as)

lpe_as_p = lpe_as.pop_in
lpe_as_r = lpe_as.res_in

lpl_as_p = lpl_as.pop_in
lpl_as_r = lpl_as.res_in

lce_as_p = lce_as.pop_in
lce_as_r = lce_as.res_in

lcl_as_p = lcl_as.pop_in
lcl_as_r = lcl_as.res_in


# In[261]:

# This gives the arrays that enable me to compute pop and resources values in small populations, for early and late 
# stages in "const" and "pop" runs

spe_as = cp(sp_e_as)
spl_as = cp(sp_l_as)
sce_as = cp(sc_e_as)
scl_as = cp(sc_l_as)

spe_as_p = spe_as.pop_in
spe_as_r = spe_as.res_in

spl_as_p = spl_as.pop_in
spl_as_r = spl_as.res_in

sce_as_p = sce_as.pop_in
sce_as_r = sce_as.res_in

scl_as_p = scl_as.pop_in
scl_as_r = scl_as.res_in


# In[183]:

# For large populations - asex-const:
lce_asp = lce_as_p[:6000]
lce_asr = lce_as_r[:6000]
lcl_asp = lcl_as_p[-6000:]
lcl_asr = lcl_as_r[-6000:]
lc_astime = range(1,6001)+range(54001,60001)

lc_asp = ['population']+[str(i) for i in lce_asp]+[str(i) for i in lcl_asp]
lc_asr = ['resources']+[str(i) for i in lce_asr]+[str(i) for i in lcl_asr]
lc_ast = ['time']+[str(i) for i in lc_astime]

lc_as = [lc_asp,lc_asr,lc_ast]
lctr_as = zip(*lc_as)
lctr_as_l = ','.join([','.join(list(i))+'\n' for i in lctr_as]).replace('\n,','\n')

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/large-sex-constres.csv','w')
z.write(lctr_as_l)
z.close()


# In[184]:

# For large populations - asex-pop:
lpe_asp = lpe_as_p[:6000]
lpe_asr = lpe_as_r[:6000]
lpl_asp = lpl_as_p[-6000:]
lpl_asr = lpl_as_r[-6000:]
lp_astime = range(1,6001)+range(54001,60001)

lp_asp = ['population']+[str(i) for i in lpe_asp]+[str(i) for i in lpl_asp]
lp_asr = ['resources']+[str(i) for i in lpe_asr]+[str(i) for i in lpl_asr]
lp_ast = ['time']+[str(i) for i in lp_astime]

lp_as = [lp_asp,lp_asr,lp_ast]
lptr_as = zip(*lp_as)
lptr_as_l = ','.join([','.join(list(i))+'\n' for i in lptr_as]).replace('\n,','\n')

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/large-sex-popres.csv','w')
z.write(lptr_as_l)
z.close()


# In[262]:

# For small populations - asex-const:
sce_asp = sce_as_p[:6000]
sce_asr = sce_as_r[:6000]
scl_asp = scl_as_p[-6000:]
scl_asr = scl_as_r[-6000:]
sc_astime = range(1,6001)+range(54001,60001)

sc_asp = ['population']+[str(i) for i in sce_asp]+[str(i) for i in scl_asp]
sc_asr = ['resources']+[str(i) for i in sce_asr]+[str(i) for i in scl_asr]
sc_ast = ['time']+[str(i) for i in sc_astime]

sc_as = [sc_asp,sc_asr,sc_ast]
sctr_as = zip(*sc_as)
sctr_as_l = ','.join([','.join(list(i))+'\n' for i in sctr_as]).replace('\n,','\n')

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/small-sex-constres.csv','w')
z.write(sctr_as_l)
z.close()

# For small populations - asex-pop:
spe_asp = spe_as_p[:6000]
spe_asr = spe_as_r[:6000]
spl_asp = spl_as_p[-6000:]
spl_asr = spl_as_r[-6000:]
sp_astime = range(1,6001)+range(54001,60001)

sp_asp = ['population']+[str(i) for i in spe_asp]+[str(i) for i in spl_asp]
sp_asr = ['resources']+[str(i) for i in spe_asr]+[str(i) for i in spl_asr]
sp_ast = ['time']+[str(i) for i in sp_astime]

sp_as = [sp_asp,sp_asr,sp_ast]
sptr_as = zip(*sp_as)
sptr_as_l = ','.join([','.join(list(i))+'\n' for i in sptr_as]).replace('\n,','\n')

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/small-sex-popres.csv','w')
z.write(sptr_as_l)
z.close()


# In[185]:

# For asex large populations, 60kth stage - surv-repr

hlpl_as = het(lpl_as.hetrz_mea)
red_lpl_as = hlpl_as.loop2()

hlcl_as = het(lcl_as.hetrz_mea)
red_lcl_as = hlcl_as.loop2()

lps_as_z = red_lpl_as[-1]
lcs_as_z = red_lcl_as[-1]

a_as = ','.join([str(i) for i in lps_as_z]).replace(',',',lps\n')+',lps\n'
b_as = ','.join([str(i) for i in lcs_as_z]).replace(',',',lcs\n')+',lcs\n'
ab_as = 'surv,group\n'+a_as+b_as

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/large_surv_repr.csv', 'w')
z.write(ab_as)
z.close()


# In[263]:

# For asex small populations, 60kth stage - surv-repr

#hspl_as = het(spl_as.hetrz_mea)
#red_spl_as = hspl_as.loop2()

#hscl_as = het(scl_as.hetrz_mea)
#red_scl_as = hscl_as.loop2()

#sps_as_z = red_spl_as[-1]
#scs_as_z = red_scl_as[-1]

#c_as = ','.join([str(i) for i in sps_as_z]).replace(',',',sps\n')+',sps\n'
#d_as = ','.join([str(i) for i in scs_as_z]).replace(',',',scs\n')+',scs\n'
#cd_as = 'surv,group\n'+c_as+d_as

#z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/small_surv_repr.csv', 'w')
#z.write(cd_as)
#z.close()


# In[187]:

# Now off to doing the asex-pop very large population - 25k individuals - pop-resources
llc_e_as = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/llce/plot_values_run1.txt'
llc_l_as = '/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/llcl/plot_values_run1.txt'

llce_as = cp(llc_e_as)
llcl_as = cp(llc_l_as)

llce_as_p = llce_as.pop_in
llce_as_r = llce_as.res_in

llcl_as_p = llcl_as.pop_in
llcl_as_r = llcl_as.res_in

# For very large populations - asex-const:
llce_asp = llce_as_p[:6000]
llce_asr = llce_as_r[:6000]
llcl_asp = llcl_as_p[-6000:]
llcl_asr = llcl_as_r[-6000:]
llc_astime = range(1,6001)+range(54001,60001)

llc_asp = ['population']+[str(i) for i in llce_asp]+[str(i) for i in llcl_asp]
llc_asr = ['resources']+[str(i) for i in llce_asr]+[str(i) for i in llcl_asr]
llc_ast = ['time']+[str(i) for i in llc_astime]

llc_as = [llc_asp,llc_asr,llc_ast]
llctr_as = zip(*llc_as)
llctr_as_l = ','.join([','.join(list(i))+'\n' for i in llctr_as]).replace('\n,','\n')

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/very-large-sex-constres.csv','w')
z.write(llctr_as_l)
z.close()


# In[265]:

# This generates the 60kth stage survival-reproduction plot for small, large, and very large populations, in the 
# constant resources condition

hllcl_as = het(llcl_as.hetrz_mea)
red_llcl_as = hllcl_as.loop2()

#lps_as_z = red_lpl_as[-1]
lcs_as_z = red_lcl_as[-1]
#sps_as_z = red_spl_as[-1]
scs_as_z = red_scl_as[-1]
llcs_as_z = red_llcl_as[-1]

#a_as = ','.join([str(i) for i in lps_as_z]).replace(',',',lps\n')+',lps\n'
b_as = ','.join([str(i) for i in lcs_as_z]).replace(',',',lcs\n')+',lcs\n'
#c_as = ','.join([str(i) for i in sps_as_z]).replace(',',',sps\n')+',sps\n'
d_as = ','.join([str(i) for i in scs_as_z]).replace(',',',scs\n')+',scs\n'
e_as = ','.join([str(i) for i in llcs_as_z]).replace(',',',llcs\n')+',llcs\n'
bde_as = 'surv,group\n'+b_as+d_as+e_as

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/const-small-large_surv_repr.csv', 'w')
z.write(bde_as)
z.close()


# In[267]:

# This box computes the expected survival and reproduction value, given the junk DNA string

# First, for the sex reproduction

Jlcr = float(lcl.repr_rate_junk_in[-1])/0.4
Jlcs = (float(lcl.surv_rate_junk_in[-1])-0.98)/0.02
Jlpr = float(lpl.repr_rate_junk_in[-1])/0.4
Jlps = (float(lpl.surv_rate_junk_in[-1])-0.98)/0.02
Jscr = float(scl.repr_rate_junk_in[-1])/0.4
Jscs = (float(scl.surv_rate_junk_in[-1])-0.98)/0.02
Jspr = float(spl.repr_rate_junk_in[-1])/0.4
Jsps = (float(spl.surv_rate_junk_in[-1])-0.98)/0.02

Jlc = ((str(Jlcs)+',')*71+(str(Jlcr)+',')*55).split(',')
Jlp = ((str(Jlps)+',')*71+(str(Jlpr)+',')*55).split(',')
Jsc = ((str(Jscs)+',')*71+(str(Jscr)+',')*55).split(',')
Jsp = ((str(Jsps)+',')*71+(str(Jspr)+',')*55).split(',')

Jnk_s = 'Jlc,Jlp,Jsc,Jsp\n'+','.join([','.join(list(i))+'\n' for i in zip(Jlc, Jlp, Jsc, Jsp)]).replace('\n,','\n')

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/Jnk_sex.csv', 'w')
z.write(Jnk_s)
z.close()


# In[268]:

# This box computes the expected survival and reproduction value, given the junk DNA string
# Asex reproduction

Jlcr_as = float(lcl_as.repr_rate_junk_in[-1])/0.4
Jlcs_as = (float(lcl_as.surv_rate_junk_in[-1])-0.98)/0.02
Jscr_as = float(scl_as.repr_rate_junk_in[-1])/0.4
Jscs_as = (float(scl_as.surv_rate_junk_in[-1])-0.98)/0.02
Jllcr_as = float(llcl_as.repr_rate_junk_in[-1])/0.4
Jllcs_as = (float(llcl_as.surv_rate_junk_in[-1])-0.98)/0.02


Jlc_as = ((str(Jlcs_as)+',')*71+(str(Jlcr_as)+',')*55).split(',')
Jsc_as = ((str(Jscs_as)+',')*71+(str(Jscr_as)+',')*55).split(',')
Jllc_as = ((str(Jllcs_as)+',')*71+(str(Jllcr_as)+',')*55).split(',')

Jnk_s_as = 'Jlc_as,Jsc_as,Jllc_as\n'+','.join([','.join(list(i))+'\n' for i in zip(Jlc_as, Jsc_as, Jllc_as)]).replace('\n,','\n')

z = open('/Volumes/group_dv/personal/DValenzano/month-by-month/Sep2015/simul-paper/asex/Jnk_sex-repr_as.csv', 'w')
z.write(Jnk_s_as)
z.close()


# In[ ]:



