from __future__ import division
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import time
import pickle


class Populations(object):

    def __init__(self,nReps,nInds,nLoci,ploidy=45, genomic_mu=0.1, selcoef=0.1, ctrl_nLoci_lower = 10, ctrl_nLoci_upper = 10,\
        ctrl_mu_lower = 0.01, ctrl_mu_upper = 0.1, ctrl_eff_lower = 0.9, ctrl_eff_upper = 0.9):

        self.nReps = nReps
        self.nInds = nInds
        self.nLoci = nLoci
        
        self.soma = np.zeros((nReps,nInds,nLoci),dtype='int')
        self.germ = np.zeros((nReps,nInds,nLoci),dtype='int')
        
        self.ploidy = ploidy
        self.mu = genomic_mu/(ploidy*nLoci)
        self.selcoef = selcoef

        # Mu rate controller in Soma
        self.ctrl_nLoci_lower = ctrl_nLoci_lower
        self.ctrl_nLoci_upper = ctrl_nLoci_upper

        self.contruct_mu_ctrl_loci()

        # Mu rate controller appearing in Germ
  
        self.ctrl_mu_lower_initial = ctrl_mu_lower/(ploidy*ctrl_nLoci_lower)
        self.ctrl_mu_upper_initial = ctrl_mu_upper/(ploidy*ctrl_nLoci_upper)

        self.ctrl_mu_lower = np.ones((nReps, nInds, ctrl_nLoci_lower))*(ctrl_mu_lower/(ploidy*ctrl_nLoci_lower))
        self.ctrl_mu_upper = np.ones((nReps, nInds, ctrl_nLoci_upper))*(ctrl_mu_upper/(ploidy*ctrl_nLoci_upper))

        self.ctrl_eff_lower = ctrl_eff_lower/ploidy
        self.ctrl_eff_upper = ctrl_eff_upper/ploidy

        self.sm_mu = np.ones((nReps, nInds, nLoci))*(genomic_mu/(ploidy*nLoci))
        self.gm_mu = np.ones((nReps, nInds, nLoci))*(genomic_mu/(ploidy*nLoci))

        self.total_sm_mu = []
        self.total_gm_mu = []
        self.total_pop_mfit = []

        self.current_step = 0



    def contruct_mu_ctrl_loci(self):

        # Mu rate controller in Soma

        self.sm_soma_ctrl_upper = np.zeros((self.nReps, self.nInds, self.ctrl_nLoci_upper), dtype = 'int') # corresponding to self.sm_soma_ctrl 
        self.sm_soma_ctrl_lower = np.zeros((self.nReps, self.nInds, self.ctrl_nLoci_lower), dtype = 'int') # corresponding to self.sm_soma_ctrl 

        self.sm_germ_ctrl_upper = np.zeros((self.nReps, self.nInds, self.ctrl_nLoci_upper), dtype = 'int') # corresponding to self.sm_germ_ctrl 
        self.sm_germ_ctrl_lower = np.zeros((self.nReps, self.nInds, self.ctrl_nLoci_lower), dtype = 'int') # corresponding to self.sm_germ_ctrl 

        # Mu rate controller appearing in Germ

        self.gm_soma_ctrl_upper = np.zeros((self.nReps, self.nInds, self.ctrl_nLoci_upper), dtype = 'int') # corresponding to self.sm_soma_ctrl 
        self.gm_soma_ctrl_lower = np.zeros((self.nReps, self.nInds, self.ctrl_nLoci_lower), dtype = 'int') # corresponding to self.sm_soma_ctrl 

        self.gm_germ_ctrl_upper = np.zeros((self.nReps, self.nInds, self.ctrl_nLoci_upper), dtype = 'int') # corresponding to self.sm_germ_ctrl 
        self.gm_germ_ctrl_lower = np.zeros((self.nReps, self.nInds, self.ctrl_nLoci_lower), dtype = 'int') # corresponding to self.sm_germ_ctrl  



    def fitness(self):
        """return a numpy array containing the fitness for each individual in each population"""
        ws = (self.soma.astype('float')/self.ploidy)*self.selcoef

        each_locus = 1-ws
        fitnesses = np.prod(each_locus, axis =2)
            
        fitnesses[fitnesses <= 0] = 0  

        return fitnesses

            
    def relative_fitness(self):
        """return a numpy array containing each individual's relative fitness"""
        w = self.fitness()
        total_w = np.sum(w,axis=1)
        totalw = np.expand_dims(total_w,axis=1)
        relfit = w/totalw
        return relfit
                   


    @staticmethod
    def calculate_mu(ctrl_upper, ctrl_lower, upper_eff, lower_eff, mu):

        sm_ws_upper = ctrl_upper.astype('float')*upper_eff
        sm_ws_lower = ctrl_lower.astype('float')*lower_eff

        sm_ctrl_locus_upper = 1+ sm_ws_upper
        sm_ctrl_locus_lower = 1- sm_ws_lower

        sm_mu = np.prod(sm_ctrl_locus_upper, axis =2)*np.prod(sm_ctrl_locus_lower, axis =2)*mu

        return sm_mu  


    def get_mu_rate(self):

        # Calculate the Mu in both Soma and Germ

        sm_mu = self.calculate_mu(self.sm_soma_ctrl_upper, self.sm_soma_ctrl_lower, self.ctrl_eff_upper, self.ctrl_eff_lower, self.mu)
        self.sm_mu = np.repeat(sm_mu[:, :, np.newaxis], self.nLoci, axis=2)

        gm_mu = self.calculate_mu(self.sm_germ_ctrl_upper, self.sm_germ_ctrl_lower, self.ctrl_eff_upper, self.ctrl_eff_lower, self.mu)    
        self.gm_mu = np.repeat(gm_mu[:, :, np.newaxis], self.nLoci, axis=2)


        sm_ctrl_mu_upper = self.calculate_mu(self.sm_soma_ctrl_upper, self.sm_soma_ctrl_lower, self.ctrl_eff_upper, self.ctrl_eff_lower, \
            self.ctrl_mu_upper_initial)
        self.sm_ctrl_mu_upper = np.repeat(sm_ctrl_mu_upper[:, :, np.newaxis], self.ctrl_nLoci_upper, axis=2)


        sm_ctrl_mu_lower = self.calculate_mu(self.sm_soma_ctrl_upper, self.sm_soma_ctrl_lower, self.ctrl_eff_upper, self.ctrl_eff_lower, \
            self.ctrl_mu_lower_initial)
        self.sm_ctrl_mu_lower = np.repeat(sm_ctrl_mu_lower[:, :, np.newaxis], self.ctrl_nLoci_lower, axis=2)


        gm_ctrl_mu_upper = self.calculate_mu(self.sm_germ_ctrl_upper, self.sm_germ_ctrl_lower, self.ctrl_eff_upper, self.ctrl_eff_lower, \
            self.ctrl_mu_upper_initial)
        self.gm_ctrl_mu_upper = np.repeat(gm_ctrl_mu_upper[:, :, np.newaxis], self.ctrl_nLoci_upper, axis=2)


        gm_ctrl_mu_lower = self.calculate_mu(self.sm_germ_ctrl_upper, self.sm_germ_ctrl_lower, self.ctrl_eff_upper, self.ctrl_eff_lower, \
            self.ctrl_mu_lower_initial)
        self.gm_ctrl_mu_lower = np.repeat(gm_ctrl_mu_lower[:, :, np.newaxis], self.ctrl_nLoci_lower, axis=2)



    @staticmethod
    def get_mutations(ctrl, ctrl_mu, ploidy):

        wt_ctrl = ploidy - ctrl
        mut_ctrl = np.random.binomial(wt_ctrl, ctrl_mu)
        return mut_ctrl



    def mutate_all_before(self):
        # Mutate both Soma, Germ and also the Mu-controller loci

        self.get_mu_rate()

        # Mutations in Soma and Germ
        soma_mutations = self.get_mutations(self.soma, self.sm_mu, self.ploidy)
        self.soma += soma_mutations

        germ_mutations = self.get_mutations(self.germ, self.gm_mu, 2)
        self.germ += germ_mutations


        # Mu controller loci
        mut_sm_soma_ctrl_upper = self.get_mutations(self.sm_soma_ctrl_upper, self.sm_ctrl_mu_upper, self.ploidy)
        self.sm_soma_ctrl_upper += mut_sm_soma_ctrl_upper

        mut_sm_soma_ctrl_lower = self.get_mutations(self.sm_soma_ctrl_lower, self.sm_ctrl_mu_lower, self.ploidy)
        self.sm_soma_ctrl_lower += mut_sm_soma_ctrl_lower


        mut_sm_germ_ctrl_upper = self.get_mutations(self.sm_germ_ctrl_upper, self.sm_ctrl_mu_upper, self.ploidy)
        self.sm_germ_ctrl_upper += mut_sm_germ_ctrl_upper

        mut_sm_germ_ctrl_lower = self.get_mutations(self.sm_germ_ctrl_lower, self.sm_ctrl_mu_lower, self.ploidy)
        self.sm_germ_ctrl_lower += mut_sm_germ_ctrl_lower


        # Then deal with the Mu controller in Germ
        mut_gm_soma_ctrl_upper = self.get_mutations(self.gm_soma_ctrl_upper, self.gm_ctrl_mu_upper, 2)
        self.gm_soma_ctrl_upper += mut_gm_soma_ctrl_upper

        mut_gm_soma_ctrl_lower = self.get_mutations(self.gm_soma_ctrl_lower, self.gm_ctrl_mu_lower, 2)
        self.gm_soma_ctrl_lower += mut_gm_soma_ctrl_lower


        mut_gm_germ_ctrl_upper = self.get_mutations(self.gm_germ_ctrl_upper, self.gm_ctrl_mu_upper, 2)
        self.gm_germ_ctrl_upper += mut_gm_germ_ctrl_upper

        mut_gm_germ_ctrl_lower = self.get_mutations(self.gm_germ_ctrl_lower, self.gm_ctrl_mu_lower, 2)
        self.gm_germ_ctrl_lower += mut_gm_germ_ctrl_lower




    
    def pick_parents_all(self):
        """Randomly choose N parents to produce offspring to populate the next generation.
        Each individual's probability of being chosen is weighted by its relative fitness.
        Wright-Fisher model method."""
        nReps,nInds = self.soma.shape[0:2]        
        relfit = self.relative_fitness()
        csrelfit = np.cumsum(relfit,axis=1)
        randvals = np.random.random((nReps,nInds))
        parents = map(np.searchsorted,csrelfit,randvals)
        return parents
    
    
    
    def mitosis_all(self,parents):
        """Generate mitotic offspring from all individuals selected as parents.
        Wright-Fisher model method."""
        nReps,nInds = self.soma.shape[0:2]        
        rReps = np.ones((nReps,nInds),dtype='int')*np.expand_dims(np.arange(nReps),axis=1)
        self.soma = self.soma[rReps,parents,]
        self.germ = self.germ[rReps,parents,]

        # Then deal with the soma and germ mu control
        sm_soma_ctrl_upper_good = (self.ploidy-self.sm_soma_ctrl_upper[rReps,parents,])*2
        sm_soma_ctrl_upper_bad = self.sm_soma_ctrl_upper[rReps,parents,]*2
        self.sm_soma_ctrl_upper = np.random.hypergeometric(sm_soma_ctrl_upper_bad, sm_soma_ctrl_upper_good, self.ploidy)

        sm_soma_ctrl_lower_good = (self.ploidy-self.sm_soma_ctrl_lower[rReps,parents,])*2
        sm_soma_ctrl_lower_bad = self.sm_soma_ctrl_lower[rReps,parents,]*2
        self.sm_soma_ctrl_lower = np.random.hypergeometric(sm_soma_ctrl_lower_bad, sm_soma_ctrl_lower_good, self.ploidy)


        sm_germ_ctrl_upper_good = (self.ploidy-self.sm_germ_ctrl_upper[rReps,parents,])*2
        sm_germ_ctrl_upper_bad = self.sm_germ_ctrl_upper[rReps,parents,]*2
        self.sm_germ_ctrl_upper = np.random.hypergeometric(sm_germ_ctrl_upper_bad, sm_germ_ctrl_upper_good, self.ploidy)

        sm_germ_ctrl_lower_good = (self.ploidy-self.sm_germ_ctrl_lower[rReps,parents,])*2
        sm_germ_ctrl_lower_bad = self.sm_germ_ctrl_lower[rReps,parents,]*2
        self.sm_germ_ctrl_lower = np.random.hypergeometric(sm_germ_ctrl_lower_bad, sm_germ_ctrl_lower_good, self.ploidy)


        self.gm_soma_ctrl_upper = self.gm_soma_ctrl_upper[rReps, parents,]
        self.gm_soma_ctrl_lower = self.gm_soma_ctrl_lower[rReps, parents,]

        self.gm_germ_ctrl_upper = self.gm_germ_ctrl_upper[rReps, parents,]
        self.gm_germ_ctrl_lower = self.gm_germ_ctrl_lower[rReps, parents,]  


        return
    


    def amitosis_all(self,parents):
        """Generate amitotic offspring from all individuals selected as parents. Only one
        amitotic offspring is generated from each parent, so this method does not reflect
        the reciprocity of amitosis.
        Wright-Fisher model method."""
        nReps,nInds = self.soma.shape[0:2]        
        rReps = np.ones((nReps,nInds),dtype='int')*np.expand_dims(np.arange(nReps),axis=1)
        good = (self.ploidy-self.soma[rReps,parents,])*2
        bad = self.soma[rReps,parents,]*2
        self.soma = np.random.hypergeometric(bad,good,self.ploidy)
        self.germ = self.germ[rReps,parents,]

        # Then deal with the soma and germ mu control
        self.sm_soma_ctrl_upper = self.sm_soma_ctrl_upper[rReps, parents,]
        self.sm_soma_ctrl_lower = self.sm_soma_ctrl_lower[rReps, parents,]

        self.sm_germ_ctrl_upper = self.sm_germ_ctrl_upper[rReps, parents,]
        self.sm_germ_ctrl_lower = self.sm_germ_ctrl_lower[rReps, parents,]         

        self.gm_soma_ctrl_upper = self.gm_soma_ctrl_upper[rReps, parents,]
        self.gm_soma_ctrl_lower = self.gm_soma_ctrl_lower[rReps, parents,]

        self.gm_germ_ctrl_upper = self.gm_germ_ctrl_upper[rReps, parents,]
        self.gm_germ_ctrl_lower = self.gm_germ_ctrl_lower[rReps, parents,] 
    
        return
    
    
    def wright_fisher_step(self,asex_type='amitosis'):
        """Take populations through a single Wright-Fisher model generation"""

        self.mutate_all_before()
        parents = self.pick_parents_all()
        if asex_type == 'amitosis':
            self.amitosis_all(parents)
        elif asex_type == 'mitosis':
            self.mitosis_all(parents)
        self.current_step += 1
        return
    
    
    #Asexul reproduction_WF model and Moran model
    
    def step(self,model='M',asex_type='amitosis',sex_freq=None):
        """Take populations through one time step if model='M', or one generation if model='WF'"""
        if model == 'M':
            self.moran_step(asex_type)
        elif model == 'WF':
            self.wright_fisher_step(asex_type)
        return 
    


    def make_gametes(self, parents):

        rReps = np.ones((self.nReps,self.nInds),dtype='int')*np.expand_dims(np.arange(self.nReps),axis=1)
        gametes = np.random.hypergeometric(self.germ[rReps,parents,],2-self.germ[rReps,parents,],1)

        gametes_soma_ctrl_upper = np.random.hypergeometric(self.gm_soma_ctrl_upper[rReps,parents,],2-self.gm_soma_ctrl_upper[rReps,parents,],1)
        gametes_soma_ctrl_lower = np.random.hypergeometric(self.gm_soma_ctrl_lower[rReps,parents,],2-self.gm_soma_ctrl_lower[rReps,parents,],1)

        gametes_germ_ctrl_upper = np.random.hypergeometric(self.gm_germ_ctrl_upper[rReps,parents,],2-self.gm_germ_ctrl_upper[rReps,parents,],1)
        gametes_germ_ctrl_lower = np.random.hypergeometric(self.gm_germ_ctrl_lower[rReps,parents,],2-self.gm_germ_ctrl_lower[rReps,parents,],1)

        return gametes, gametes_soma_ctrl_upper, gametes_soma_ctrl_lower, gametes_germ_ctrl_upper, gametes_germ_ctrl_lower



    def make_zygotes(self, random_mating=False):

        self.mutate_all_before()

        if random_mating == False:
            parents = self.pick_parents_all()
            these_g1 = self.make_gametes(parents)
            those_g1 = self.make_gametes(parents)

        else: 
            first_parents = self.pick_parents_all()
            second_parents = self.pick_parents_all()

            these_g1 = self.make_gametes(first_parents)
            those_g1 = self.make_gametes(second_parents)           

        these_gametes = these_g1[0]
        these_g_soma_ctrl_upper = these_g1[1]
        these_g_soma_ctrl_lower = these_g1[2]
        these_g_germ_ctrl_upper = these_g1[3]
        these_g_germ_ctrl_lower = these_g1[4]


        those_gametes = those_g1[0]
        those_g_soma_ctrl_upper = those_g1[1]
        those_g_soma_ctrl_lower = those_g1[2]
        those_g_germ_ctrl_upper = those_g1[3]
        those_g_germ_ctrl_lower = those_g1[4]


        zygotes = these_gametes + those_gametes
        zy_soma_ctrl_upper = these_g_soma_ctrl_upper + those_g_soma_ctrl_upper
        zy_soma_ctrl_lower = these_g_soma_ctrl_lower + those_g_soma_ctrl_lower

        zy_germ_ctrl_upper = these_g_germ_ctrl_upper + those_g_germ_ctrl_upper
        zy_germ_ctrl_lower = these_g_germ_ctrl_lower + those_g_germ_ctrl_lower

        return zygotes, zy_soma_ctrl_upper, zy_soma_ctrl_lower, zy_germ_ctrl_upper, zy_germ_ctrl_lower



    @staticmethod
    def get_soma_allele_after_sex(germ, ploidy):

        soma = np.zeros_like(germ)

        hetcount = len(germ[germ == 1])
        hetvals = int(ploidy/2) + np.random.binomial(ploidy-2*int(ploidy/2),0.5, hetcount)

        soma[germ == 0] = 0
        soma[germ == 1] = hetvals
        soma[germ == 2] = ploidy

        return soma


    def sex(self,random_mating=False):
        
        zy = self.make_zygotes(random_mating)  # make zygotes. Here the random_mating is just an argument of self.make
                                                    # _zygotes, doesn't mean that here is the random_mating
        self.germ = zy[0]  # the germline genome will be the same as the zygotes.
        self.gm_soma_ctrl_upper = zy[1]
        self.gm_soma_ctrl_lower = zy[2]
        self.gm_germ_ctrl_upper = zy[3]
        self.gm_germ_ctrl_lower = zy[4]
        
        self.soma = self.get_soma_allele_after_sex(self.germ, self.ploidy)
        
        self.sm_soma_ctrl_upper = self.get_soma_allele_after_sex(self.gm_soma_ctrl_upper, self.ploidy)
        self.sm_soma_ctrl_lower = self.get_soma_allele_after_sex(self.gm_soma_ctrl_lower, self.ploidy)        

        self.sm_germ_ctrl_upper = self.get_soma_allele_after_sex(self.gm_germ_ctrl_upper, self.ploidy)
        self.sm_germ_ctrl_lower = self.get_soma_allele_after_sex(self.gm_germ_ctrl_lower, self.ploidy)  

        self.current_step +=1
                             
        return


    
    def get_results(self):
        """calculate stuff like mean fitness, Gst, and number of fixed mutations"""

        W = self.fitness()
        mW = np.nanmean(W)
        log_mW = np.log(mW)

        pop_mean_fit = np.nanmean(W, axis =1)
        self.total_pop_mfit.append(pop_mean_fit)
        mean_fit = np.nanmean(pop_mean_fit)
        pop_mfit_std = np.nanstd(pop_mean_fit) 


        total_fit_var = []
        for i in range(self.nReps):
            fit_var = np.var(W[i])
            total_fit_var.append(fit_var)

        varW_mean = np.nanmean(total_fit_var)
        varW_std = np.nanstd(total_fit_var)

        # Then calculate the mean mu in Soma and Germ per locus
        self.get_mu_rate()
        
        sm_mu_per_ind = np.nanmean(self.sm_mu, axis =2)
        sm_mu_per_pop = np.nanmean(sm_mu_per_ind, axis =1)

        self.total_sm_mu.append(sm_mu_per_pop)

        sm_mu_mean = np.nanmean(sm_mu_per_pop)
        sm_mu_std = np.nanstd(sm_mu_per_pop)


        gm_mu_per_ind = np.nanmean(self.gm_mu, axis =2)
        gm_mu_per_pop = np.nanmean(gm_mu_per_ind, axis =1)

        self.total_gm_mu.append(gm_mu_per_pop)
        
        gm_mu_mean = np.nanmean(gm_mu_per_pop)
        gm_mu_std = np.nanstd(gm_mu_per_pop)

        # Calculate the mean genomic Mu in Soma and Germ

        sm_genomic_mu_per_ind = np.sum(self.sm_mu*self.ploidy, axis =2)
        sm_genomic_mu_per_pop = np.nanmean(sm_genomic_mu_per_ind, axis =1)

        sm_genomic_mu_mean = np.nanmean(sm_genomic_mu_per_pop)
        sm_genomic_mu_std = np.nanstd(sm_genomic_mu_per_pop)


        gm_genomic_mu_per_ind = np.sum(self.gm_mu*2, axis =2)
        gm_genomic_mu_per_pop = np.nanmean(gm_genomic_mu_per_ind, axis =1)

        gm_genomic_mu_mean = np.nanmean(gm_genomic_mu_per_pop)
        gm_genomic_mu_std = np.nanstd(gm_genomic_mu_per_pop)        



        return [mW, log_mW, mean_fit, pop_mfit_std, varW_mean, varW_std, sm_mu_mean, sm_mu_std, gm_mu_mean, gm_mu_std, \
        sm_genomic_mu_mean, sm_genomic_mu_std,  gm_genomic_mu_mean, gm_genomic_mu_std]
    
    


    def simulate(self,stepcount,model='M',asex_type='amitosis',sex_freq=None,random_mating=False,strides=10):

        results = [self.get_results()]
        

        start = time.time()
    
        while self.current_step <= stepcount:
            if (sex_freq != None) and (self.current_step%sex_freq == 0):
                self.sex(random_mating)
                
            else:
                self.step(model,asex_type,sex_freq)

            
            if self.current_step%strides == 0:
                results.append(self.get_results())

        colnames = ['meanFit','ln_meanFit','PopMeanFit_Mean', 'PopMeanFit_STD', 'PopVar_Mean', 'PopVar_STD',\
        'SomaMu_Mean', 'SomaMu_Std', 'GermMu_Mean', 'GermMu_Std', \
        'SomaGenomicMu_Mean', 'SomaGenomicMu_Std', 'GermGenomicMu_Mean', 'GermGenomicMu_Std']

        results = pd.DataFrame(np.array(results),columns=colnames)
        


        end = time.time()
        print 'TOTAL TIME: ',end-start
        return results
        

    def simulateNsave(self,outfile,stepcount,model='M',asex_type='amitosis',sex_freq=None,random_mating=False,strides=10):
        """same as simulate method except results are written to the file specified by outfile"""
        results = self.simulate(stepcount,model,asex_type,sex_freq,random_mating,strides)
        results.to_csv(outfile,index=False)
        return
            
         

        
def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)




amito_50 =Populations(nReps = 100,nInds = 2000, nLoci = 100, ploidy= 45, genomic_mu= 0.1, selcoef=0.1, ctrl_nLoci_lower = 10, ctrl_nLoci_upper = 10,\
    ctrl_mu_lower = 0.002, ctrl_mu_upper = 0.01, ctrl_eff_lower = 0.9, ctrl_eff_upper = 0.9)

amito_50.simulateNsave('Fit_FitAmito_CtrlMito_N2K_DeleOnly_UP01LW002_MuEvo_191002R2.csv', 10*1000, model ='WF', asex_type='amitosis',sex_freq=None,random_mating=False,strides=1)

save_object(amito_50.total_sm_mu, 'FitAmito_CtrlMito_N2K_DeleOnly_UP01LW002_Total_SM_MU_191002R2')
save_object(amito_50.total_gm_mu, 'FitAmito_CtrlMito_N2K_DeleOnly_UP01LW002_Total_GM_MU_191002R2')
save_object(amito_50.total_pop_mfit, 'FitAmito_CtrlMito_N2K_DeleOnly_UP01LW002_Total_Pop_MFit_191002R2')

save_object(amito_50.sm_soma_ctrl_upper, 'FitAmito_CtrlMito_N2K_DeleOnly_UP01LW002_SM_Soma_Ctrl_Upper_191002R2')
save_object(amito_50.sm_soma_ctrl_lower, 'FitAmito_CtrlMito_N2K_DeleOnly_UP01LW002_SM_Soma_Ctrl_Lower_191002R2')

save_object(amito_50.sm_germ_ctrl_upper, 'FitAmito_CtrlMito_N2K_DeleOnly_UP01LW002_SM_Germ_Ctrl_Upper_191002R2')
save_object(amito_50.sm_germ_ctrl_lower, 'FitAmito_CtrlMito_N2K_DeleOnly_UP01LW002_SM_Germ_Ctrl_Lower_191002R2')

save_object(amito_50.gm_soma_ctrl_upper, 'FitAmito_CtrlMito_N2K_DeleOnly_UP01LW002_GM_Soma_Ctrl_Upper_191002R2')
save_object(amito_50.gm_soma_ctrl_lower, 'FitAmito_CtrlMito_N2K_DeleOnly_UP01LW002_GM_Soma_Ctrl_Lower_191002R2')

save_object(amito_50.gm_germ_ctrl_upper, 'FitAmito_CtrlMito_N2K_DeleOnly_UP01LW002_GM_Germ_Ctrl_Upper_191002R2')
save_object(amito_50.gm_germ_ctrl_lower, 'FitAmito_CtrlMito_N2K_DeleOnly_UP01LW002_GM_Germ_Ctrl_Lower_191002R2')

save_object(amito_50.soma,'FitAmito_CtrlMito_N2K_DeleOnly_UP01LW002_Soma_191002R2')
save_object(amito_50.germ,'FitAmito_CtrlMito_N2K_DeleOnly_UP01LW002_Germ_191002R2')