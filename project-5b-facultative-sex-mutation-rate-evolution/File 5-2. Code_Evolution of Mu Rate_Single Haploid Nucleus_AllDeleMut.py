from __future__ import division
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import time
import pickle


class Populations(object):

    def __init__(self,nReps,nInds,nLoci,genomic_mu=0.1, selcoef=0.1, ctrl_nLoci_lower = 10, ctrl_nLoci_upper = 10,\
        ctrl_mu_lower = 0.01, ctrl_mu_upper = 0.1, ctrl_eff_lower = 0.9, ctrl_eff_upper = 0.9):

        self.nReps = nReps
        self.nInds = nInds
        self.nLoci = nLoci
        
        self.soma = np.zeros((nReps,nInds,nLoci),dtype='int')
        
        self.mu = genomic_mu/(1*nLoci)
        self.selcoef = selcoef

        # Mu rate controller in Soma
        self.ctrl_nLoci_lower = ctrl_nLoci_lower
        self.ctrl_nLoci_upper = ctrl_nLoci_upper

        self.contruct_mu_ctrl_loci()

        # Mu rate controller appearing in Germ
  
        self.ctrl_mu_lower_initial = ctrl_mu_lower/(1*ctrl_nLoci_lower)
        self.ctrl_mu_upper_initial = ctrl_mu_upper/(1*ctrl_nLoci_upper)

        self.ctrl_mu_lower = np.ones((nReps, nInds, ctrl_nLoci_lower))*(ctrl_mu_lower/(1*ctrl_nLoci_lower))
        self.ctrl_mu_upper = np.ones((nReps, nInds, ctrl_nLoci_upper))*(ctrl_mu_upper/(1*ctrl_nLoci_upper))

        self.ctrl_eff_lower = ctrl_eff_lower/1
        self.ctrl_eff_upper = ctrl_eff_upper/1

        self.sm_mu = np.ones((nReps, nInds, nLoci))*(genomic_mu/(1*nLoci))

        self.total_sm_mu = []
        self.total_pop_mfit = []

        self.total_fit_var = []
        self.total_mu_var = []
        self.total_covariance = [] # store the covariance between fitness and mutation rate

        self.current_step = 0



    def contruct_mu_ctrl_loci(self):

        # Mu rate controller in Soma

        self.sm_ctrl_upper = np.zeros((self.nReps, self.nInds, self.ctrl_nLoci_upper), dtype = 'int') # corresponding to self.sm_soma_ctrl 
        self.sm_ctrl_lower = np.zeros((self.nReps, self.nInds, self.ctrl_nLoci_lower), dtype = 'int') # corresponding to self.sm_soma_ctrl 
 




    def fitness(self):
        """return a numpy array containing the fitness for each individual in each population"""
        ws = (self.soma.astype('float')/1)*self.selcoef

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

        sm_mu = self.calculate_mu(self.sm_ctrl_upper, self.sm_ctrl_lower, self.ctrl_eff_upper, self.ctrl_eff_lower, self.mu)
        self.sm_mu = np.repeat(sm_mu[:, :, np.newaxis], self.nLoci, axis=2)


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
        soma_mutations = self.get_mutations(self.soma, self.sm_mu, 1)
        self.soma += soma_mutations


        # Mu controller loci
        mut_sm_ctrl_upper = self.get_mutations(self.sm_ctrl_upper, self.sm_ctrl_mu_upper, 1)
        self.sm_ctrl_upper += mut_sm_ctrl_upper

        mut_sm_ctrl_lower = self.get_mutations(self.sm_ctrl_lower, self.sm_ctrl_mu_lower, 1)
        self.sm_ctrl_lower += mut_sm_ctrl_lower




    
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
    

    
    def asexual_reproduction(self,parents):
        """Generate mitotic offspring from all individuals selected as parents.
        Wright-Fisher model method."""
        nReps,nInds = self.soma.shape[0:2]        
        rReps = np.ones((nReps,nInds),dtype='int')*np.expand_dims(np.arange(nReps),axis=1)

        self.soma = self.soma[rReps,parents,]
        self.sm_ctrl_upper = self.sm_ctrl_upper[rReps, parents,]
        self.sm_ctrl_lower = self.sm_ctrl_lower[rReps, parents,]     

        return

    

    def wright_fisher_step(self):
        """Take populations through a single Wright-Fisher model generation"""

        self.mutate_all_before()
        parents = self.pick_parents_all()
        self.asexual_reproduction(parents)

        self.current_step += 1
        return
       


    def make_gametes(self, parents):

        rReps = np.ones((self.nReps,self.nInds),dtype='int')*np.expand_dims(np.arange(self.nReps),axis=1)

        gametes = self.soma[rReps, parents,]
        gametes_ctrl_upper = self.sm_ctrl_upper[rReps, parents,]
        gametes_ctrl_lower = self.sm_ctrl_lower[rReps, parents,]

        return gametes, gametes_ctrl_upper, gametes_ctrl_lower



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
        these_g_ctrl_upper = these_g1[1]
        these_g_ctrl_lower = these_g1[2]

        those_gametes = those_g1[0]
        those_g_ctrl_upper = those_g1[1]
        those_g_ctrl_lower = those_g1[2]


        zygotes = these_gametes + those_gametes
        zy_ctrl_upper = these_g_ctrl_upper + those_g_ctrl_upper
        zy_ctrl_lower = these_g_ctrl_lower + those_g_ctrl_lower


        return zygotes, zy_ctrl_upper, zy_ctrl_lower





    def sex(self,random_mating=False):
        
        zy = self.make_zygotes(random_mating)  # make zygotes. Here the random_mating is just an argument of self.make
                                                    # _zygotes, doesn't mean that here is the random_mating

        self.soma = np.random.hypergeometric(zy[0], 2-zy[0], 1)
        self.sm_ctrl_mu_upper = np.random.hypergeometric(zy[1], 2-zy[1], 1)
        self.sm_ctrl_mu_lower = np.random.hypergeometric(zy[2], 2-zy[2], 1)        

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


        fit_var = np.var(W, axis =1)
        self.total_fit_var.append(fit_var)

        varW_mean = np.nanmean(fit_var)
        varW_std = np.nanstd(fit_var)

        # Then calculate the mean mu in Soma and Germ per locus
        self.get_mu_rate()
        
        sm_mu_per_ind = np.nanmean(self.sm_mu, axis =2)
        sm_mu_per_pop = np.nanmean(sm_mu_per_ind, axis =1)

        self.total_sm_mu.append(sm_mu_per_pop)

        sm_mu_mean = np.nanmean(sm_mu_per_pop)
        sm_mu_std = np.nanstd(sm_mu_per_pop)


        sm_mu_var = np.var(sm_mu_per_ind, axis =1)
        self.total_mu_var.append(sm_mu_var)

        sm_mu_var_mean = np.nanmean(sm_mu_var)
        sm_mu_var_std = np.nanstd(sm_mu_var)


        # Calculate the covariance between fitness and mutation rate
        cov_each_pop = []
        for i in range(self.nReps):
            covariance = np.cov(W[i], self.sm_mu[:,:, 0][i])[0,1]
            cov_each_pop.append(covariance)

        self.total_covariance.append(cov_each_pop)

        cov_mean = np.nanmean(cov_each_pop)
        cov_std = np.nanstd(cov_each_pop)


        return [mW, log_mW, mean_fit, pop_mfit_std, varW_mean, varW_std, sm_mu_mean, sm_mu_std, \
        sm_mu_var_mean, sm_mu_var_std, cov_mean, cov_std]
    
    


    def simulate(self,stepcount,sex_freq=None,random_mating=False,strides=10):

        results = [self.get_results()]
        

        start = time.time()
    
        while self.current_step <= stepcount:
            if (sex_freq != None) and (self.current_step%sex_freq == 0):
                self.sex(random_mating)
                
            else:
                self.wright_fisher_step()

            
            if self.current_step%strides == 0:
                results.append(self.get_results())

        colnames = ['meanFit','ln_meanFit','PopMeanFit_Mean', 'PopMeanFit_STD', 'PopVar_Mean', 'PopVar_STD',\
        'SomaMu_Mean', 'SomaMu_Std', 'SomaMuVar_Mean', 'SomaMuVar_Std', 'Cov_FitMu_Mean', 'Cov_FitMu_Std']

        results = pd.DataFrame(np.array(results),columns=colnames)
        
        end = time.time()
        print 'TOTAL TIME: ',end-start
        return results
        

    def simulateNsave(self,outfile,stepcount,sex_freq=None,random_mating=False,strides=10):
        """same as simulate method except results are written to the file specified by outfile"""
        results = self.simulate(stepcount,sex_freq,random_mating,strides)
        results.to_csv(outfile,index=False)
        return
            
         

        
def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)




amito_50 =Populations(nReps = 100,nInds = 2000, nLoci = 50, genomic_mu= 0.1*1/45*0.5, selcoef=0.1, ctrl_nLoci_lower = 5, ctrl_nLoci_upper = 5,\
    ctrl_mu_lower = 0.002*1/45*0.5, ctrl_mu_upper = 0.01*1/45*0.5, ctrl_eff_lower = 0.9, ctrl_eff_upper = 0.9)

amito_50.simulateNsave('Fit_RME100_SgP1Nucl_N2K_Dele_MuEvo_50Loci.csv', 10*1000, sex_freq=100,random_mating=True,strides=1)

save_object(amito_50.total_sm_mu, 'RME100_SgP1Nucl_N2K_Dele_Total_SM_MU_50Loci')
save_object(amito_50.total_pop_mfit, 'RME100_SgP1Nucl_N2K_Dele_Total_Pop_MFit_50Loci')
save_object(amito_50.total_fit_var, 'RME100_SgP1Nucl_N2K_Dele_Total_Pop_Fit_Var_50Loci')
save_object(amito_50.total_mu_var, 'RME100_SgP1Nucl_N2K_Dele_Total_SM_Mu_Var_50Loci')
save_object(amito_50.total_covariance, 'RME100_SgP1Nucl_N2K_Dele_Covariance_FitvsMu_50Loci')

save_object(amito_50.sm_ctrl_upper, 'RME100_SgP1Nucl_N2K_Dele_SM_Ctrl_Upper_50Loci')
save_object(amito_50.sm_ctrl_lower, 'RME100_SgP1Nucl_N2K_Dele_SM_Ctrl_Lower_50Loci')


