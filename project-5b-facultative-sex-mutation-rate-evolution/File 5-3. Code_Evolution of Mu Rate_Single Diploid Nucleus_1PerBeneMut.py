from __future__ import division
import numpy as np
from scipy import stats
import pandas as pd
import time
import pickle

class Populations(object):

    def __init__(self,nReps,nInds,nLoci,ploidy=45, dele_genomic_mu=0.1, bene_genomic_mu= 0.01*0.1, beneficial_selcoef=0.1, deleterious_selcoef = 0.1, \
        ctrl_nLoci_lower = 10, ctrl_nLoci_upper = 10, ctrl_mu_lower = 0.01, ctrl_mu_upper = 0.1, ctrl_eff_lower = 0.9, ctrl_eff_upper = 0.9):

        self.nReps = nReps
        self.nInds = nInds
        self.nLoci = nLoci
        
        self.soma_bene = np.zeros((nReps,nInds,nLoci),dtype='int') # store the beneficial mutations
        self.soma_dele = np.zeros((nReps,nInds,nLoci),dtype='int') # store the deleterious mutations
        
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

        # self.total_sm_mu_bene = []
        self.total_sm_mu_dele = []
        self.total_pop_mfit = []

        self.total_fit_var = []
        self.total_mu_dele_var = []
        self.total_dele_covariance = [] # store the covariance between fitness and mutation rate
        
        self.current_step = 0
        


    def contruct_mu_ctrl_loci(self):

        # Mu rate controller in Soma

        self.sm_ctrl_upper = np.zeros((self.nReps, self.nInds, self.ctrl_nLoci_upper), dtype = 'int') # corresponding to self.sm_soma_ctrl 
        self.sm_ctrl_lower = np.zeros((self.nReps, self.nInds, self.ctrl_nLoci_lower), dtype = 'int') # corresponding to self.sm_soma_ctrl 




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

        sm_mu_bene = self.calculate_mu(self.sm_ctrl_upper, self.sm_ctrl_lower, self.ctrl_eff_upper, self.ctrl_eff_lower, \
            self.mu_bene)
        self.sm_mu_bene = np.repeat(sm_mu_bene[:, :, np.newaxis], self.nLoci, axis=2)

        sm_mu_dele = self.calculate_mu(self.sm_ctrl_upper, self.sm_ctrl_lower, self.ctrl_eff_upper, self.ctrl_eff_lower, \
            self.mu_dele)
        self.sm_mu_dele = np.repeat(sm_mu_dele[:, :, np.newaxis], self.nLoci, axis=2)


        # Calculate the Mu in controller loci
        sm_ctrl_mu_upper = self.calculate_mu(self.sm_ctrl_upper, self.sm_ctrl_lower, self.ctrl_eff_upper, self.ctrl_eff_lower, \
            self.ctrl_mu_upper_initial)
        self.sm_ctrl_mu_upper = np.repeat(sm_ctrl_mu_upper[:, :, np.newaxis], self.ctrl_nLoci_upper, axis=2)


        sm_ctrl_mu_lower = self.calculate_mu(self.sm_ctrl_upper, self.sm_ctrl_lower, self.ctrl_eff_upper, self.ctrl_eff_lower, \
            self.ctrl_mu_lower_initial)
        self.sm_ctrl_mu_lower = np.repeat(sm_ctrl_mu_lower[:, :, np.newaxis], self.ctrl_nLoci_lower, axis=2)




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


        # Mu controller loci
        mut_sm_ctrl_upper = self.get_mutations(self.sm_ctrl_upper, self.sm_ctrl_mu_upper, self.ploidy)
        self.sm_ctrl_upper += mut_sm_ctrl_upper

        mut_sm_ctrl_lower = self.get_mutations(self.sm_ctrl_lower, self.sm_ctrl_mu_lower, self.ploidy)
        self.sm_ctrl_lower += mut_sm_ctrl_lower



    
    

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
        self.soma_dele = self.soma_dele[rReps,parents,]
 
        # Mu control
        self.sm_ctrl_upper = self.sm_ctrl_upper[rReps, parents,]
        self.sm_ctrl_lower = self.sm_ctrl_lower[rReps, parents,]

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


        # Then deal with the soma and germ mu control

        sm_ctrl_upper_good = (self.ploidy-self.sm_ctrl_upper[rReps,parents,])*2
        sm_ctrl_upper_bad = self.sm_ctrl_upper[rReps,parents,]*2
        self.sm_ctrl_upper = np.random.hypergeometric(sm_ctrl_upper_bad, sm_ctrl_upper_good, self.ploidy)

        sm_ctrl_lower_good = (self.ploidy-self.sm_ctrl_lower[rReps,parents,])*2
        sm_ctrl_lower_bad = self.sm_ctrl_lower[rReps,parents,]*2
        self.sm_ctrl_lower = np.random.hypergeometric(sm_ctrl_lower_bad, sm_ctrl_lower_good, self.ploidy)


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

        gametes_bene = np.random.hypergeometric(self.soma_bene[rReps,parents,],self.ploidy-self.soma_bene[rReps,parents,],1)
        gametes_dele = np.random.hypergeometric(self.soma_dele[rReps,parents,],self.ploidy-self.soma_dele[rReps,parents,],1)        

        # Mu controller
        gametes_ctrl_upper = np.random.hypergeometric(self.sm_ctrl_upper[rReps,parents,],self.ploidy-self.sm_ctrl_upper[rReps,parents,],1)
        gametes_ctrl_lower = np.random.hypergeometric(self.sm_ctrl_lower[rReps,parents,],self.ploidy-self.sm_ctrl_lower[rReps,parents,],1)

        return gametes_bene, gametes_dele, gametes_ctrl_upper, gametes_ctrl_lower


    

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
        these_gametes_ctrl_upper = these[2]
        these_gametes_ctrl_lower = these[3]

        those_gametes_bene = those[0]
        those_gametes_dele = those[1]
        those_gametes_ctrl_upper = those[2]
        those_gametes_ctrl_lower = those[3]


        zy_bene = these_gametes_bene + those_gametes_bene
        zy_dele = these_gametes_dele + those_gametes_dele

        zy_ctrl_upper = these_gametes_ctrl_upper + those_gametes_ctrl_upper
        zy_ctrl_lower = these_gametes_ctrl_lower + those_gametes_ctrl_lower

     
        return zy_bene, zy_dele, zy_ctrl_upper, zy_ctrl_lower
    


    
    def sex(self,random_mating=False):
    
        zy = self.make_zygotes(random_mating)  

        self.soma_bene = zy[0]
        self.soma_dele = zy[1]
        self.sm_ctrl_upper = zy[2]
        self.sm_ctrl_lower = zy[3]
    
        self.current_step +=1
        return     
    
    

    def get_results(self):
        """calculate stuff like mean fitness, Gst, and number of fixed mutations"""
        """calculate stuff like mean fitness, Gst, and number of fixed mutations"""

        W = self.fitness()
        mW = np.nanmean(W)
        log_mW = np.log(mW)

        pop_mean_fit = np.nanmean(W, axis =1)
        self.total_pop_mfit.append(pop_mean_fit)

        mean_fit = np.nanmean(pop_mean_fit)
        pop_mfit_std = np.nanstd(pop_mean_fit) 


        fit_var = np.var(W, axis =1)
        self.total_fit_var.append(fit_var)

        varW_mean = np.nanmean(fit_var)
        varW_std = np.nanstd(fit_var)

        # Then calculate the mean mu in Soma and Germ per locus
        self.get_mu_rate()
        
        sm_mu_bene_per_ind = np.nanmean(self.sm_mu_bene, axis =2)
        sm_mu_bene_per_pop = np.nanmean(sm_mu_bene_per_ind, axis =1)

        # self.total_sm_mu_bene.append(sm_mu_bene_per_pop)
        sm_mu_bene_mean = np.nanmean(sm_mu_bene_per_pop)
        sm_mu_bene_std = np.nanstd(sm_mu_bene_per_pop)

        sm_mu_dele_per_ind = np.nanmean(self.sm_mu_dele, axis =2)
        sm_mu_dele_per_pop = np.nanmean(sm_mu_dele_per_ind, axis =1)

        self.total_sm_mu_dele.append(sm_mu_dele_per_pop)
        sm_mu_dele_mean = np.nanmean(sm_mu_dele_per_pop)
        sm_mu_dele_std = np.nanstd(sm_mu_dele_per_pop)
    

        sm_mu_dele_var = np.var(sm_mu_dele_per_ind, axis =1)
        self.total_mu_dele_var.append(sm_mu_dele_var)

        sm_mu_dele_var_mean = np.nanmean(sm_mu_dele_var)
        sm_mu_dele_var_std = np.nanstd(sm_mu_dele_var)


        # Calculate the covariance between fitness and mutation rate
        cov_each_pop = []
        for i in range(self.nReps):
            covariance = np.cov(W[i], self.sm_mu_dele[:,:, 0][i])[0,1]
            cov_each_pop.append(covariance)

        self.total_dele_covariance.append(cov_each_pop)

        dele_cov_mean = np.nanmean(cov_each_pop)
        dele_cov_std = np.nanstd(cov_each_pop)

            
        return [mW, log_mW, mean_fit, pop_mfit_std, varW_mean, varW_std, sm_mu_bene_mean, sm_mu_bene_std, sm_mu_dele_mean, sm_mu_dele_std, \
        sm_mu_dele_var_mean, sm_mu_dele_var_std, dele_cov_mean, dele_cov_std]
    
    

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

                
        colnames = ['meanFit','ln_meanFit','PopMeanFit_Mean', 'PopMeanFit_STD', 'PopVar_Mean', 'PopVar_STD',\
        'SomaMu_Bene_Mean', 'SomaMu_Bene_Std', 'SomaMu_Dele_Mean', 'SomaMu_Dele_Std', \
        'SomaMu_Dele_Var_Mean', 'SomaMu_Dele_Var_Std', 'Cov_FitDeleMu_Mean', 'Cov_FitDeleMu_Std']


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


        

amito_50 =Populations(nReps = 100,nInds = 2000,nLoci = 100,ploidy=2, dele_genomic_mu=0.1*2/45, bene_genomic_mu= 0.01*0.1*2/45, beneficial_selcoef=0.1, deleterious_selcoef = 0.1, \
        ctrl_nLoci_lower = 10, ctrl_nLoci_upper = 10, ctrl_mu_lower = 0.002*2/45, ctrl_mu_upper = 0.01*2/45, ctrl_eff_lower = 0.9, ctrl_eff_upper = 0.9)


amito_50.simulateNsave('Fit_RME100M_P2SN_N2K_Bene01_MuEvo_200207.csv', 20*1000, model ='WF', asex_type='mitosis',sex_freq=100,random_mating=True,strides=1)


save_object(amito_50.total_sm_mu_dele, 'RME100M_P2SN_N2K_Bene01_Total_SM_MU_Dele_200207')
save_object(amito_50.total_pop_mfit, 'RME100M_P2SN_N2K_Bene01_Total_Pop_MFit_200207')

save_object(amito_50.total_fit_var, 'RME100M_P2SN_N2K_Bene01_Total_Pop_Fit_Var_200207')
save_object(amito_50.total_mu_dele_var, 'RME100M_P2SN_N2K_Bene01_Total_SM_DeleMu_Var_200207')
save_object(amito_50.total_dele_covariance, 'RME100M_P2SN_N2K_Bene01_Covariance_FitvsDeleMu_200207')

save_object(amito_50.sm_ctrl_upper, 'RME100M_P2SN_N2K_Bene01_SM_Ctrl_Upper_200207')
save_object(amito_50.sm_ctrl_lower, 'RME100M_P2SN_N2K_Bene01_SM_Ctrl_Lower_200207')

