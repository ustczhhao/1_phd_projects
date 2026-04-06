from __future__ import division
import numpy as np
from scipy import stats
import pandas as pd
import time
import pickle


class Populations(object):

    def __init__(self,nReps,nInds,nLoci,ploidy=45,genomic_mu=0.1,selcoef=0.1, model_version ='ADD'):

        self.nReps = nReps
        self.nInds = nInds
        self.nLoci = nLoci
        
        self.soma = np.zeros((nReps,nInds,nLoci),dtype='int')
        self.germ = np.zeros((nReps,nInds,nLoci),dtype='int')
        self.ploidy = ploidy
        self.mu = genomic_mu/(nLoci*ploidy)
        self.selcoef = selcoef
        self.current_step = 0

        self. sex_gen = []
        self. sex_mean_fit = []  # The mean fitness of each population after they just have sex 

        self.model_version = model_version
        

        

    def fitness(self):
        """return a numpy array containing the fitness for each individual in each population"""
        ws = (self.soma.astype('float')/self.ploidy)*self.selcoef

        if self. model_version == 'ADD':
            fitnesses = 1 - np.sum(ws,axis=2)
            
        elif self. model_version == 'MUL':
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

                

        
    #Wright_Fisher model
    def mutate_all_before(self):
        """Mutate all individuals in the population with mutation rate mu.
        Wright-Fisher model method. Use this one instead of mutate_all if you want to mutate the population before selection."""
        wt_soma = self.ploidy - self.soma
        wt_germ = 2 - self.germ
        soma_mutations = np.random.binomial(wt_soma,self.mu)
        germ_mutations = np.random.binomial(wt_germ,self.mu)
        self.soma += soma_mutations
        self.germ += germ_mutations
        return
    
    
    
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

        # self.mutate_all_before()
        # parents = self.pick_parents_all(model_version)

        rReps = np.ones((self.nReps,self.nInds),dtype='int')*np.expand_dims(np.arange(self.nReps),axis=1)
        gametes = np.random.hypergeometric(self.germ[rReps,parents,],2-self.germ[rReps,parents,],1)

        return gametes




    def make_zygotes(self, random_mating=False):

        self.mutate_all_before()

        if random_mating == False:
            parents = self.pick_parents_all()

            these_gametes = self.make_gametes(parents)
            those_gametes = self.make_gametes(parents)

        else: 
            first_parents = self.pick_parents_all()
            second_parents = self.pick_parents_all()

            these_gametes = self.make_gametes(first_parents)
            those_gametes = self.make_gametes(second_parents)           


        zygotes = these_gametes + those_gametes

        return zygotes




    def sex(self,random_mating=False):
        
        zygotes = self.make_zygotes(random_mating)  # make zygotes. Here the random_mating is just an argument of self.make
                                                    # _zygotes, doesn't mean that here is the random_mating
        self.germ = zygotes  # the germline genome will be the same as the zygotes.
        

        
        hetcount = len(self.germ[self.germ==1])  # count the heterozygote loci. self.germ ==1 is just the heterozygote loci as
                                                # one of the two alleles at that locus is mutation and the other is wild type.
        hetvals = np.random.binomial(self.ploidy,0.5,hetcount) 
        self.soma[self.germ==1] = hetvals # replace the heterozygote sites (in somatic genome) with the hetvals.
        
        self.soma[self.germ==2] = self.ploidy # self.germ == 2 means the germline is homozygote of mutation. So the somatic 
                                            # developed from this germline should also be homozygote of mutation (beause there
                                            # is only mutation at the germline, so every time it only have mutation to pick to
                                            # form the somatic genome). Here replace these sites with self.ploidy (means all copies
                                            # at that locus are mutation).
        self.soma[self.germ==0] = 0  # self.germ == 0 means the germline is homozygote of wild type (0 means no mutation copy at 
                                    # that locus). Same reason as self.germ ==2, the somatic developed from this germline should 
                                    # also be homozygote of wild type. Here replace these sites with 0 (means all copies at that
                                    # locus are wild type).
        self.current_step +=1
                             
        return 




    
    def get_results(self):
        """calculate stuff like mean fitness, Gst, and number of fixed mutations"""

        
        W = self.fitness()

        mW = np.nanmean(W)
        log_mW = np.log(mW)


        pop_mean_fit = np.nanmean(W, axis =1)  # Calculate mean fitness from another way: first get the mean fitness of each population
        mean_fit = np.nanmean(pop_mean_fit)  # Then average across the populations

        pop_mfit_std = np.nanstd(pop_mean_fit)    # Get the STD of population mean fitness across each population
        pop_mfit_lower = np.percentile(pop_mean_fit, 2.5)   # The 2.5% lower bound of population mean fitness among populations
        pop_mfit_upper = np.percentile(pop_mean_fit, 97.5)   # The 97.5% upper bound of population mean fitness among populations


        total_fit_var = []     
        for i in range(self.nReps):  
            fit_var = np.var(W[i])   # Variance within each population
            total_fit_var.append(fit_var)


        varW_mean = np.nanmean(total_fit_var)  # Mean variance among populations

        varW_std = np.nanstd(total_fit_var)    # STD of variance among populations
        varW_lower = np.percentile(total_fit_var, 2.5)   # 2.5% lower bound of variance among populations
        varW_upper = np.percentile(total_fit_var, 97.5)  # 97.5% upper bound of variance among populations


        pop_fit_lower = np.percentile(W, 25, axis =1)   # The 25% lower bound of fitness within each population
        pop_fit_upper = np.percentile(W, 75, axis =1)   # The 75% upper bound of fitness within each population

        pop_fit_50_range = pop_fit_upper - pop_fit_lower  # The middle 50% fitness range within each population


        pop_fit_r50_mean = np.nanmean(pop_fit_50_range)   # Mean of middle 50% fitness range among populations
        pop_fit_r50_std = np.nanstd(pop_fit_50_range)    # STD of middle 50% fitness range among populations
        pop_fit_r50_lower = np.percentile(pop_fit_50_range, 2.5)  # 2.5% lower bound of middle 50% fitness range among populations
        pop_fit_r50_upper = np.percentile(pop_fit_50_range, 97.5) # 97.5% lower bound of middle 50% fitness range among populations

        
            
        return [mW, log_mW, mean_fit, pop_mfit_std, pop_mfit_lower, pop_mfit_upper, varW_mean, varW_std, varW_lower, varW_upper, \
        pop_fit_r50_mean, pop_fit_r50_std, pop_fit_r50_lower, pop_fit_r50_upper]



        
    def simulate(self,stepcount,model='M',asex_type='amitosis',sex_freq=None,random_mating=False,strides=10):

        results = [self.get_results()]
        
        start = time.time()


    
        while self.current_step <= stepcount:
            if (sex_freq != None) and (self.current_step%sex_freq == 0):
                self.sex(random_mating)

                self.sex_gen.append(self.current_step)

                Fit = self.fitness()
                sex_pop_mean_fit = np.nanmean(Fit, axis =1)
                self.sex_mean_fit.append(sex_pop_mean_fit)


            else:
                self.step(model,asex_type,sex_freq)

            if self.current_step%strides == 0:
                results.append(self.get_results())

                
        colnames = ['meanFit','ln_meanFit','PopMeanFit_Mean', 'PopMeanFit_STD', 'PopMeanFit_Lower', 'PopMeanFit_Upper', \
        'PopVar_Mean', 'PopVar_STD', 'PopVar_Lower', 'PopVar_Upper', 'PopR50Fit_Mean', 'PopR50Fit_STD', 'PopR50Fit_Lower', \
        'PopR50Fit_Upper']
        
        results = pd.DataFrame(np.array(results),columns=colnames)


        self.slope_list_1 = self.get_regression_slope(self.nReps, self.sex_mean_fit, sex_freq, -100) # linear slope of last 100G
        self.slope_list_2 = self.get_regression_slope(self.nReps, self.sex_mean_fit, sex_freq, -500) # linear slope of last 500G
        self.slope_list_3 = self.get_regression_slope(self.nReps, self.sex_mean_fit, sex_freq, -1000)
        self.slope_list_4 = self.get_regression_slope(self.nReps, self.sex_mean_fit, sex_freq, -2000)
        self.slope_list_5 = self.get_regression_slope(self.nReps, self.sex_mean_fit, sex_freq, -5000)
        self.slope_list_6 = self.get_regression_slope(self.nReps, self.sex_mean_fit, sex_freq, -10000)

        end = time.time()
        print 'TOTAL TIME: ',end-start
        return results
        

    def simulateNsave(self,outfile, stepcount,model='M',asex_type='amitosis',sex_freq=None,random_mating=False,strides=10):
        """same as simulate method except results are written to the file specified by outfile"""
        results = self.simulate(stepcount,model,asex_type,sex_freq,random_mating,strides)
 
        results.to_csv(outfile,index=False)

        return


         
    @staticmethod
    def get_regression_slope(nReps, sex_mean_fit, sex_freq, start_index, end_index = None):
        
        if end_index == None:
            aver_nReps_all_generation = sex_mean_fit[int(start_index/sex_freq): end_index]
        else:
            aver_nReps_all_generation = sex_mean_fit[int(start_index/sex_freq): int(end_index/sex_freq)]
        
            
        total_fitness_each_nReps = []
        
        for i in range(nReps):
            fitness_each_nReps = []
            for j in aver_nReps_all_generation:
                fitness_each_nReps.append(j[i])
            
            total_fitness_each_nReps.append(fitness_each_nReps)
    
        total_linear_regression_slope = [] 
        
        x = [i*sex_freq for i in list(range(len(total_fitness_each_nReps[0])))]  # stride = 1
        for i in total_fitness_each_nReps:
            slope = stats.linregress(x, i)[0]
            total_linear_regression_slope.append(slope)

        total_linear_regression_slope = np.array(total_linear_regression_slope)
            
        return total_linear_regression_slope


def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)   




sf_5 =Populations(nReps =100, nInds =1000, nLoci =100, ploidy =45, genomic_mu=0.1, selcoef =0.1, model_version ='MUL')
sf_5.simulateNsave('Fitness_SF_E1K_N1000_180510.csv', 20*1000, model ='WF', asex_type='amitosis',sex_freq=1000,random_mating=False,strides=1)


save_object(sf_5.sex_gen, 'SF_N1K_E1K_SEX_GEN')
save_object(sf_5.sex_mean_fit, 'SF_N1K_E1K_SEX_POP_MFIT')

save_object(sf_5.slope_list_1, 'SF_N1K_E1K_100G_SLOPE_0510')
save_object(sf_5.slope_list_2, 'SF_N1K_E1K_500G_SLOPE_0510')
save_object(sf_5.slope_list_3, 'SF_N1K_E1K_1KG_SLOPE_0510')
save_object(sf_5.slope_list_4, 'SF_N1K_E1K_2KG_SLOPE_0510')
save_object(sf_5.slope_list_5, 'SF_N1K_E1K_5KG_SLOPE_0510')
save_object(sf_5.slope_list_6, 'SF_N1K_E1K_10KG_SLOPE_0510')


save_object(sf_5.soma, 'SF_N1K_E1K_SOMA')
save_object(sf_5.germ, 'SF_N1K_E1K_GERM')
save_object(sf_5.current_step, 'SF_N1K_E1K_STEP')



rm_5 =Populations(nReps =100, nInds =1000, nLoci =100, ploidy =45, genomic_mu=0.1, selcoef =0.1, model_version ='MUL')
rm_5.simulateNsave('Fitness_RM_E1K_N1000_180510.csv', 20*1000, model ='WF', asex_type='amitosis',sex_freq=1000,random_mating=True,strides=1)


save_object(rm_5.sex_gen, 'RM_N1K_E1K_SEX_GEN')
save_object(rm_5.sex_mean_fit, 'RM_N1K_E1K_SEX_POP_MFIT')

save_object(rm_5.slope_list_1, 'RM_N1K_E1K_100G_SLOPE_0510')
save_object(rm_5.slope_list_2, 'RM_N1K_E1K_500G_SLOPE_0510')
save_object(rm_5.slope_list_3, 'RM_N1K_E1K_1KG_SLOPE_0510')
save_object(rm_5.slope_list_4, 'RM_N1K_E1K_2KG_SLOPE_0510')
save_object(rm_5.slope_list_5, 'RM_N1K_E1K_5KG_SLOPE_0510')
save_object(rm_5.slope_list_6, 'RM_N1K_E1K_10KG_SLOPE_0510')


save_object(rm_5.soma, 'RM_N1K_E1K_SOMA')
save_object(rm_5.germ, 'RM_N1K_E1K_GERM')
save_object(rm_5.current_step, 'RM_N1K_E1K_STEP')