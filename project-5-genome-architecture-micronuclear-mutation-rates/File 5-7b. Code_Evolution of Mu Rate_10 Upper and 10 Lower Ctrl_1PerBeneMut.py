from __future__ import division
import numpy as np
from scipy import stats
import pandas as pd
import time
import pickle
import multiprocessing as mp


class Populations(object):

    def __init__(self,nReps,nInds,nLoci,ploidy=45, dele_genomic_mu=0.1, bene_genomic_mu= 0.01*0.1, beneficial_selcoef=0.1, deleterious_selcoef = 0.1, \
        ctrl_nLoci_lower = 10, ctrl_nLoci_upper = 10, ctrl_mu_lower = 0.01, ctrl_mu_upper = 0.1, ctrl_eff_lower = 0.9, ctrl_eff_upper = 0.9):

        self.nReps = nReps
        self.nInds = nInds
        self.nLoci = nLoci
        
        self.soma_bene = np.zeros((nReps,nInds,nLoci),dtype='int') # store the beneficial mutations
        self.soma_dele = np.zeros((nReps,nInds,nLoci),dtype='int') # store the deleterious mutations
        
        self.germ_bene = np.zeros((nReps,nInds,nLoci),dtype='int')
        self.germ_dele = np.zeros((nReps,nInds,nLoci),dtype='int')        
        
        self.ploidy = ploidy

        self.mu_bene = bene_genomic_mu/(nLoci*ploidy)
        self.mu_dele = dele_genomic_mu/(nLoci*ploidy)

        self.bene_selcoef = beneficial_selcoef
        self.dele_selcoef = deleterious_selcoef

        # Mu control loci in Soma and Germ
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

        self.sm_mu_bene = np.ones((nReps, nInds, nLoci))*(bene_genomic_mu/(nLoci*ploidy))
        self.sm_mu_dele = np.ones((nReps, nInds, nLoci))*(dele_genomic_mu/(nLoci*ploidy))

        self.gm_mu_bene = np.ones((nReps, nInds, nLoci))*(bene_genomic_mu/(nLoci*ploidy))
        self.gm_mu_dele = np.ones((nReps, nInds, nLoci))*(dele_genomic_mu/(nLoci*ploidy))

        self.total_sm_mu_bene = []
        self.total_sm_mu_dele = []

        self.total_gm_mu_bene = []
        self.total_gm_mu_dele = []

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
        self.gm_germ_ctrl_lower = np.zeros((self.nReps, self.nInds, self.ctrl_nLoci_lower), dtype = 'int') # corr
        


    def fitness(self):
        """return a numpy array containing the fitness for each individual in each population"""
        
        # need to add beneficial mutations
        bene_ws = (self.soma_bene.astype('float')/self.ploidy)*self.bene_selcoef
        dele_ws = (self.soma_dele.astype('float')/self.ploidy)*self.dele_selcoef

        bene_locus = 1 + bene_ws
        dele_locus = 1- dele_ws

        fitnesses = np.prod(bene_locus, axis =2)*np.prod(dele_locus, axis =2)
            
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

        sm_mu_bene = self.calculate_mu(self.sm_soma_ctrl_upper, self.sm_soma_ctrl_lower, self.ctrl_eff_upper, self.ctrl_eff_lower, \
            self.mu_bene)
        self.sm_mu_bene = np.repeat(sm_mu_bene[:, :, np.newaxis], self.nLoci, axis=2)

        sm_mu_dele = self.calculate_mu(self.sm_soma_ctrl_upper, self.sm_soma_ctrl_lower, self.ctrl_eff_upper, self.ctrl_eff_lower, \
            self.mu_dele)
        self.sm_mu_dele = np.repeat(sm_mu_dele[:, :, np.newaxis], self.nLoci, axis=2)


        gm_mu_bene = self.calculate_mu(self.sm_germ_ctrl_upper, self.sm_germ_ctrl_lower, self.ctrl_eff_upper, self.ctrl_eff_lower, \
            self.mu_bene)    
        self.gm_mu_bene = np.repeat(gm_mu_bene[:, :, np.newaxis], self.nLoci, axis=2)

        gm_mu_dele = self.calculate_mu(self.sm_germ_ctrl_upper, self.sm_germ_ctrl_lower, self.ctrl_eff_upper, self.ctrl_eff_lower, \
            self.mu_dele)    
        self.gm_mu_dele = np.repeat(gm_mu_dele[:, :, np.newaxis], self.nLoci, axis=2)


        # Calculate the Mu in controller loci
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
        soma_bene_mutations = self.get_mutations(self.soma_bene, self.sm_mu_bene, self.ploidy)
        self.soma_bene += soma_bene_mutations

        soma_dele_mutations = self.get_mutations(self.soma_dele, self.sm_mu_dele, self.ploidy)
        self.soma_dele += soma_dele_mutations

        germ_bene_mutations = self.get_mutations(self.germ_bene, self.gm_mu_bene, 2)
        self.germ_bene += germ_bene_mutations

        germ_dele_mutations = self.get_mutations(self.germ_dele, self.gm_mu_dele, 2)
        self.germ_dele += germ_dele_mutations


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
        nReps,nInds = self.soma_bene.shape[0:2]        
        relfit = self.relative_fitness()

        csrelfit = np.cumsum(relfit,axis=1)
        randvals = np.random.random((nReps,nInds))
        parents = map(np.searchsorted,csrelfit,randvals)
        return parents    
    
    

    def mitosis_all(self,parents):
        """Generate mitotic offspring from all individuals selected as parents.
        Wright-Fisher model method."""
        nReps,nInds = self.soma_bene.shape[0:2]        
        rReps = np.ones((nReps,nInds),dtype='int')*np.expand_dims(np.arange(nReps),axis=1)
       
        self.soma_bene = self.soma_bene[rReps,parents,]
        self.germ_bene = self.germ_bene[rReps,parents,]        

        self.soma_dele = self.soma_dele[rReps,parents,]
        self.germ_dele = self.germ_dele[rReps,parents,] 

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
        
        # may consider multi hypergeometric distribution
        
        nReps,nInds = self.soma_bene.shape[0:2]        
        rReps = np.ones((nReps,nInds),dtype='int')*np.expand_dims(np.arange(nReps),axis=1)
        
        soma_bene_good = (self.ploidy-self.soma_bene[rReps,parents,])*2
        soma_bene_bad = self.soma_bene[rReps,parents,]*2
        self.soma_bene = np.random.hypergeometric(soma_bene_bad, soma_bene_good, self.ploidy)

        soma_dele_good = (self.ploidy-self.soma_dele[rReps,parents,])*2
        soma_dele_bad = self.soma_dele[rReps,parents,]*2
        self.soma_dele = np.random.hypergeometric(soma_dele_bad, soma_dele_good, self.ploidy)

        self.germ_bene = self.germ_bene[rReps, parents,]
        self.germ_dele = self.germ_dele[rReps, parents,]


        # Mu control
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
    
    
    # Asexul reproduction_WF model and Moran model
    
    def step(self,model='M',asex_type='amitosis',sex_freq=None):
        """Take populations through one time step if model='M', or one generation if model='WF'"""
        if model == 'M':
            self.moran_step(asex_type)
        elif model == 'WF':
            self.wright_fisher_step(asex_type)
        return 

    
    
    def make_gametes(self, parents):

        rReps = np.ones((self.nReps,self.nInds),dtype='int')*np.expand_dims(np.arange(self.nReps),axis=1)

        gametes_bene = np.random.hypergeometric(self.germ_bene[rReps,parents,],2-self.germ_bene[rReps,parents,],1)
        gametes_dele = np.random.hypergeometric(self.germ_dele[rReps,parents,],2-self.germ_dele[rReps,parents,],1)        

        # Mu controller
        gametes_soma_ctrl_upper = np.random.hypergeometric(self.gm_soma_ctrl_upper[rReps,parents,],2-self.gm_soma_ctrl_upper[rReps,parents,],1)
        gametes_soma_ctrl_lower = np.random.hypergeometric(self.gm_soma_ctrl_lower[rReps,parents,],2-self.gm_soma_ctrl_lower[rReps,parents,],1)

        gametes_germ_ctrl_upper = np.random.hypergeometric(self.gm_germ_ctrl_upper[rReps,parents,],2-self.gm_germ_ctrl_upper[rReps,parents,],1)
        gametes_germ_ctrl_lower = np.random.hypergeometric(self.gm_germ_ctrl_lower[rReps,parents,],2-self.gm_germ_ctrl_lower[rReps,parents,],1)

        return gametes_bene, gametes_dele, gametes_soma_ctrl_upper, gametes_soma_ctrl_lower, \
        gametes_germ_ctrl_upper, gametes_germ_ctrl_lower


    

    def make_zygotes(self, random_mating=False):

        self.mutate_all_before()

        if random_mating == False:
            parents = self.pick_parents_all()

            these = self.make_gametes(parents)
            those = self.make_gametes(parents) 

        else: 
            first_parents = self.pick_parents_all()
            second_parents = self.pick_parents_all()

            these = self.make_gametes(first_parents)
            those = self.make_gametes(second_parents)           

        these_gametes_bene = these[0]
        these_gametes_dele = these[1]
        these_gametes_soma_ctrl_upper = these[2]
        these_gametes_soma_ctrl_lower = these[3]
        these_gametes_germ_ctrl_upper = these[4]
        these_gametes_germ_ctrl_lower = these[5]

        those_gametes_bene = those[0]
        those_gametes_dele = those[1]
        those_gametes_soma_ctrl_upper = those[2]
        those_gametes_soma_ctrl_lower = those[3]
        those_gametes_germ_ctrl_upper = those[4]
        those_gametes_germ_ctrl_lower = those[5]

        zy_bene = these_gametes_bene + those_gametes_bene
        zy_dele = these_gametes_dele + those_gametes_dele

        zy_soma_ctrl_upper = these_gametes_soma_ctrl_upper + those_gametes_soma_ctrl_upper
        zy_soma_ctrl_lower = these_gametes_soma_ctrl_lower + those_gametes_soma_ctrl_lower

        zy_germ_ctrl_upper = these_gametes_germ_ctrl_upper + those_gametes_germ_ctrl_upper
        zy_germ_ctrl_lower = these_gametes_germ_ctrl_lower + those_gametes_germ_ctrl_lower
        
        return zy_bene, zy_dele, zy_soma_ctrl_upper, zy_soma_ctrl_lower, zy_germ_ctrl_upper, zy_germ_ctrl_lower
    


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
    
        zy = self.make_zygotes(random_mating)  

        self.germ_bene = zy[0]
        self.germ_dele = zy[1]
        self.gm_soma_ctrl_upper = zy[2]
        self.gm_soma_ctrl_lower = zy[3]
        self.gm_germ_ctrl_upper = zy[4]   
        self.gm_germ_ctrl_lower = zy[5]


        self.soma_bene = self.get_soma_allele_after_sex(self.germ_bene, self.ploidy)
        self.soma_dele = self.get_soma_allele_after_sex(self.germ_dele, self.ploidy)


        # Then deal with the soma_mu_control
        self.sm_soma_ctrl_upper = self.get_soma_allele_after_sex(self.gm_soma_ctrl_upper, self.ploidy)
        self.sm_soma_ctrl_lower = self.get_soma_allele_after_sex(self.gm_soma_ctrl_lower, self.ploidy)        

        self.sm_germ_ctrl_upper = self.get_soma_allele_after_sex(self.gm_germ_ctrl_upper, self.ploidy)
        self.sm_germ_ctrl_lower = self.get_soma_allele_after_sex(self.gm_germ_ctrl_lower, self.ploidy) 

        self.current_step +=1
        return     
    
    

    def get_results(self):
        """calculate stuff like mean fitness, Gst, and number of fixed mutations"""
        """calculate stuff like mean fitness, Gst, and number of fixed mutations"""

        W = self.fitness()
        mW = np.nanmean(W)

        pop_mean_fit = np.nanmean(W, axis =1)
        self.total_pop_mfit.append(pop_mean_fit)

        mean_fit = np.nanmean(pop_mean_fit)

        total_fit_var = []
        for i in range(self.nReps):
            fit_var = np.var(W[i])
            total_fit_var.append(fit_var)

        varW_mean = np.nanmean(total_fit_var)

        # Then calculate the mean mu in Soma and Germ per locus
        self.get_mu_rate()
        
        sm_mu_bene_per_ind = np.nanmean(self.sm_mu_bene, axis =2)
        sm_mu_bene_per_pop = np.nanmean(sm_mu_bene_per_ind, axis =1)

        self.total_sm_mu_bene.append(sm_mu_bene_per_pop)
        sm_mu_bene_mean = np.nanmean(sm_mu_bene_per_pop)


        sm_mu_dele_per_ind = np.nanmean(self.sm_mu_dele, axis =2)
        sm_mu_dele_per_pop = np.nanmean(sm_mu_dele_per_ind, axis =1)

        self.total_sm_mu_dele.append(sm_mu_dele_per_pop)
        sm_mu_dele_mean = np.nanmean(sm_mu_dele_per_pop)


        gm_mu_bene_per_ind = np.nanmean(self.gm_mu_bene, axis =2)
        gm_mu_bene_per_pop = np.nanmean(gm_mu_bene_per_ind, axis =1)

        self.total_gm_mu_bene.append(gm_mu_bene_per_pop)
        gm_mu_bene_mean = np.nanmean(gm_mu_bene_per_pop)


        gm_mu_dele_per_ind = np.nanmean(self.gm_mu_dele, axis =2)
        gm_mu_dele_per_pop = np.nanmean(gm_mu_dele_per_ind, axis =1)

        self.total_gm_mu_dele.append(gm_mu_dele_per_pop)
        gm_mu_dele_mean = np.nanmean(gm_mu_dele_per_pop)

            
        return [mW, mean_fit, varW_mean, sm_mu_bene_mean, sm_mu_dele_mean, gm_mu_bene_mean,  gm_mu_dele_mean]
    
    

    def simulate(self,stepcount,model='M',asex_type='amitosis',sex_freq=None,random_mating=False,strides=10):
        results = [self.get_results()]
        
        start = time.time()
   
        while self.current_step <= stepcount:
            if (sex_freq != None) and (self.current_step%sex_freq == 0):
                self.sex(random_mating)
                # self.sex_gen.append(self.current_step)
              
            else:
                self.step(model,asex_type,sex_freq)
                
            if self.current_step%strides == 0:
#                 print 'STEP', self.current_step
                results.append(self.get_results())


                
        colnames = ['meanFit','PopMeanFit_Mean',  'PopVar_Mean', \
        'SomaMu_Bene_Mean', 'SomaMu_Dele_Mean','GermMu_Bene_Mean', 'GermMu_Dele_Mean']


        results = pd.DataFrame(np.array(results),columns=colnames)

        total_sm_mu_bene = self.total_sm_mu_bene
        total_sm_mu_dele = self.total_sm_mu_dele
        total_gm_mu_bene = self.total_gm_mu_bene
        total_gm_mu_dele = self.total_gm_mu_dele
        total_pop_mfit = self.total_pop_mfit

        sm_soma_ctrl_upper = self.sm_soma_ctrl_upper
        sm_soma_ctrl_lower = self.sm_soma_ctrl_lower
        sm_germ_ctrl_upper = self.sm_germ_ctrl_upper
        sm_germ_ctrl_lower = self.sm_germ_ctrl_lower

        gm_soma_ctrl_upper = self.gm_soma_ctrl_upper
        gm_soma_ctrl_lower = self.gm_soma_ctrl_lower
        gm_germ_ctrl_upper = self.gm_germ_ctrl_upper
        gm_germ_ctrl_lower = self.gm_germ_ctrl_lower   

        soma_bene = self.soma_bene
        soma_dele = self.soma_dele
        germ_bene = self.germ_bene
        germ_dele = self.germ_dele    
    
        end = time.time()
        print 'TOTAL TIME: ',end-start
        
        return results, total_sm_mu_bene, total_sm_mu_dele, total_gm_mu_bene, total_gm_mu_dele, total_pop_mfit, \
        sm_soma_ctrl_upper, sm_soma_ctrl_lower, sm_germ_ctrl_upper, sm_germ_ctrl_lower, \
        gm_soma_ctrl_upper, gm_soma_ctrl_lower, gm_germ_ctrl_upper, gm_germ_ctrl_lower, \
        soma_bene, soma_dele, germ_bene, germ_dele

    
    def simulateNsave(self,outfile,stepcount,model='M',asex_type='amitosis',sex_freq=None,random_mating=False,strides=10):
        """same as simulate method except results are written to the file specified by outfile"""
        results = self.simulate(stepcount,model,asex_type,sex_freq,random_mating,strides)
        
        results.to_csv(outfile,index=False)
             
        return


def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)



def get_result(data_list, data_name):

    '''Get the data from the stored data_list according to their name'''

    n = len(data_list)
    mean_data = []

    for i in range(n):
        mean_data.append(list(data_list[i][data_name]))

    m = len(list(data_list[0][data_name]))
    # print 'M', m

    each_gen_mean_data = []

    each_gen_std_data = [] # Need to check

    for j in range(m):
        each_gen_data = []
        for s in range(n):
            each_gen_data.append(mean_data[s][j])

        each_gen_mean_data.append(np.nanmean(each_gen_data))
        each_gen_std_data.append(np.nanstd(each_gen_data))  # Need check

        if j == 0:
            print 'Initial', data_name, each_gen_data
            print data_name, len(list(set(each_gen_data)))

        elif j == m-1:
            print 'Final', data_name, each_gen_data
            print data_name, len(list(set(each_gen_data)))
  
    return each_gen_mean_data, each_gen_std_data



def transform_to_gen_data(replicate_data):

    rep = len(replicate_data)
    generation = len(replicate_data[0])

    generation_data = []

    for i in range(generation):
        each_gen_data = []
        for j in range(rep):
            each_gen_data.append(replicate_data[j][i])

        generation_data.append(each_gen_data)

    return generation_data



def worker(n):
    '''define a function worker to conduct a single task. Here the worker was used to run a single replicate for the
    whole simulation process. Here for I am using object-oriented programming, to let the code run, the function worker 
    is used to run the object. The function worker should have parameters otherwise it will raise error when running. But 
    the parameter n here won't be used later (here n is just to make sure that the function has parameters)'''

    np.random.seed()

    a =Populations( nReps = 1,nInds = 2000,nLoci = 100,ploidy=45, dele_genomic_mu=0.1, bene_genomic_mu= 0.01*0.1, beneficial_selcoef=0.1, deleterious_selcoef = 0.1, \
        ctrl_nLoci_lower = 10, ctrl_nLoci_upper = 10, ctrl_mu_lower = 0.002, ctrl_mu_upper = 0.01, ctrl_eff_lower = 0.9, ctrl_eff_upper = 0.9)

    results = a.simulate(10*1000, model ='WF', asex_type='amitosis',sex_freq=None,random_mating=False,strides=1)

    return results




if __name__ == "__main__":  # add this line can make sure that multiprocess can be tested in Windows PC. 
                    

    start = time.time()

    pool = mp.Pool(5)   # creat a multiprocessing pool with 4 CPU (should be identical to the number of requested
                           # CPU in the .sh file, or = requested CPU-1)

    print 'CPU COUNT', mp.cpu_count()  
        
    results = pool.map_async(worker, list(range(100))) # most time consuming part, use parallel computing
                                                # Here map_async is asynchronical function, and what this function do
                                                # is to map asynchronical function to our defined function worker.
                                                # list(range(10)) is the parameters used for worker, and here it means run
                                                # worker for 10 times (i.e., the number of replicates)

    results = results.get()  # get the results
 
    fitness = [i[0] for i in results]

    total_sm_mu_bene = [i[1] for i in results]
    total_sm_mu_dele = [i[2] for i in results]
    total_gm_mu_bene = [i[3] for i in results]
    total_gm_mu_dele = [i[4] for i in results]
    total_pop_mfit = [i[5] for i in results]


    sm_soma_ctrl_upper = [i[6] for i in results]
    sm_soma_ctrl_lower = [i[7] for i in results]
    sm_germ_ctrl_upper = [i[8] for i in results]
    sm_germ_ctrl_lower = [i[9] for i in results]

    gm_soma_ctrl_upper = [i[10] for i in results]
    gm_soma_ctrl_lower = [i[11] for i in results]
    gm_germ_ctrl_upper = [i[12] for i in results]
    gm_germ_ctrl_lower = [i[13] for i in results] 

    soma_bene = [i[14] for i in results]
    soma_dele = [i[15] for i in results]
    germ_bene = [i[16] for i in results]
    germ_dele = [i[17] for i in results] 



    gen_mfit = get_result(fitness, 'meanFit')
    gen_mfit_mean = gen_mfit[0]
    gen_mfit_std = gen_mfit[1]

    gen_pop_mfit = get_result(fitness, 'PopMeanFit_Mean') 
    gen_pop_mfit_mean = gen_pop_mfit[0]
    gen_pop_mfit_std = gen_pop_mfit[1]  

    gen_varfit = get_result(fitness, 'PopVar_Mean')
    gen_varfit_mean = gen_varfit[0]
    gen_varfit_std = gen_varfit[1]

    gen_sm_mu_bene = get_result(fitness, 'SomaMu_Bene_Mean')
    gen_sm_mu_bene_mean = gen_sm_mu_bene[0]
    gen_sm_mu_bene_std = gen_sm_mu_bene[1]

    gen_sm_mu_dele = get_result(fitness, 'SomaMu_Dele_Mean')
    gen_sm_mu_dele_mean = gen_sm_mu_dele[0]
    gen_sm_mu_dele_std = gen_sm_mu_dele[1]

    gen_gm_mu_bene = get_result(fitness, 'GermMu_Bene_Mean')
    gen_gm_mu_bene_mean = gen_gm_mu_bene[0]
    gen_gm_mu_bene_std = gen_gm_mu_bene[1]

    gen_gm_mu_dele = get_result(fitness, 'GermMu_Dele_Mean')
    gen_gm_mu_dele_mean = gen_gm_mu_dele[0]
    gen_gm_mu_dele_std = gen_gm_mu_dele[1]

     
    fitness_result = []

    for j in range(len(gen_mfit_mean)):

        fitness_result.append([gen_mfit_mean[j], gen_mfit_std[j], gen_pop_mfit_mean[j], gen_pop_mfit_std[j], gen_varfit_mean[j], gen_varfit_std[j], \
            gen_sm_mu_bene_mean[j], gen_sm_mu_bene_std[j], gen_sm_mu_dele_mean[j], gen_sm_mu_dele_std[j], \
            gen_gm_mu_bene_mean[j], gen_gm_mu_bene_std[j], gen_gm_mu_dele_mean[j], gen_gm_mu_dele_std[j]])

    colnames = ['mFit_Mean','mFit_Std', 'PopMFit_Mean','PopMFit_Std','varFit_Mean', 'varFit_Std', \
    'SomaMu_Bene_Mean', 'SomaMu_Bene_Std', 'SomaMu_Dele_Mean', 'SomaMu_Dele_Std', \
    'GermMu_Bene_Mean', 'GermMu_Bene_Std', 'GermMu_Dele_Mean', 'GermMu_Dele_Std']

    fitness_result = pd.DataFrame(np.array(fitness_result),columns=colnames)
    fitness_result.to_csv('Fit_FACM_N2K_Bene01_UP01LW002_MuEvo_190910R_MP2.csv', index =False)


    total_sm_mu_bene_eg = transform_to_gen_data(total_sm_mu_bene)
    total_sm_mu_dele_eg = transform_to_gen_data(total_sm_mu_dele)    
    total_gm_mu_bene_eg = transform_to_gen_data(total_gm_mu_bene)
    total_gm_mu_dele_eg = transform_to_gen_data(total_gm_mu_dele)   
    total_pop_mfit_eg = transform_to_gen_data(total_pop_mfit)


    save_object(total_sm_mu_bene_eg, 'FACM_N2K_Bene01_UP01LW002_Total_SM_MU_Bene_190910R_MP2')
    save_object(total_sm_mu_dele_eg,'FACM_N2K_Bene01_UP01LW002_Total_SM_MU_Dele_190910R_MP2')
    save_object(total_gm_mu_bene_eg,'FACM_N2K_Bene01_UP01LW002_Total_GM_MU_Bene_190910R_MP2')
    save_object(total_gm_mu_dele_eg,'FACM_N2K_Bene01_UP01LW002_Total_GM_MU_Dele_190910R_MP2')
    save_object(total_pop_mfit_eg,'FACM_N2K_Bene01_UP01LW002_Total_Pop_MFit_190910R_MP2')

    save_object(sm_soma_ctrl_upper, 'FACM_N2K_Bene01_UP01LW002_SM_Soma_Ctrl_Upper_190910R_MP2')
    save_object(sm_soma_ctrl_lower, 'FACM_N2K_Bene01_UP01LW002_SM_Soma_Ctrl_Lower_190910R_MP2')
    save_object(sm_germ_ctrl_upper, 'FACM_N2K_Bene01_UP01LW002_SM_Germ_Ctrl_Upper_190910R_MP2')
    save_object(sm_germ_ctrl_lower, 'FACM_N2K_Bene01_UP01LW002_SM_Germ_Ctrl_Lower_190910R_MP2')   

    save_object(gm_soma_ctrl_upper, 'FACM_N2K_Bene01_UP01LW002_GM_Soma_Ctrl_Upper_190910R_MP2')
    save_object(gm_soma_ctrl_lower, 'FACM_N2K_Bene01_UP01LW002_GM_Soma_Ctrl_Lower_190910R_MP2')
    save_object(gm_germ_ctrl_upper, 'FACM_N2K_Bene01_UP01LW002_GM_Germ_Ctrl_Upper_190910R_MP2')
    save_object(gm_germ_ctrl_lower, 'FACM_N2K_Bene01_UP01LW002_GM_Germ_Ctrl_Lower_190910R_MP2') 

    save_object(soma_bene,'FACM_N2K_Bene01_UP01LW002_Soma_Bene_190910R_MP2')
    save_object(soma_dele,'FACM_N2K_Bene01_UP01LW002_Soma_Dele_190910R_MP2')
    save_object(germ_bene, 'FACM_N2K_Bene01_UP01LW002_Germ_Bene_190910R_MP2')
    save_object(germ_dele, 'FACM_N2K_Bene01_UP01LW002_Germ_Dele_190910R_MP2')    



    end = time.time()
    print 'TOTAL TIME:', end-start

