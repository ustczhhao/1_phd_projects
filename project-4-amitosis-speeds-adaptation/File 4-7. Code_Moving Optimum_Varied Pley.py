from __future__ import division
import numpy as np
from scipy import stats
import scipy.spatial as spa
import numpy.random as rnd
import copy
import time
import pandas as pd
import math
import pickle
import multiprocessing as mp


class Population(object):
    '''Define a single population in FGM. Assuming the fixed optimum is at the orgin of the Coordinate System'''

    def __init__(self,size, nLoci, dim, pleiotropy_level, fit,mut_rate,sigma, ploidy =1):
        
        

        self.size = size  # The population size
        self.nLoci = nLoci
        self.dim = dim  # The dimension of the landscape

        self.pleiotropy_level = pleiotropy_level # How many dimenstions that a single locus correspond

        self.initial_fitness = fit  # The initial fitness of the population

        self.mutation_rate = mut_rate/(nLoci*ploidy) # mutation rate

        self.sigma = sigma/np.sqrt(pleiotropy_level)  # The STD of generating new mutations
        
        self.ploidy = ploidy  # The ploidy of the individual within the population
        
        self.soma_mut_moving_pos = np.zeros((size, nLoci, ploidy, dim)) # These two np.array will be used to store how much that each mutation in each locus in
        # self.germ_mut_moving_pos = np.zeros((size, nLoci, 2, dim))     # Soma and Germ will move along the corresponding dimension. Initially no mutations
                                                                # Thus the values within the array are 0.

        self.pleiotropy_matrix = self.generate_pleiotropy()
        self.mut_covariance = self.generate_mut_covariance()
        
        self.generate_ancestor()  # Generate ancestors with pre-defined initial fitness
        self.create_pop()  # Create populations composed of ancestors
        
        self.current_step =0  # Will be used to indicate the number of generations that have been simulated."
        
        
    def calculate_pop_fitness(self): 
        '''Calculate the fitness of each individual within the whole population simultaneously'''
        
        euclidean_distance = np.sqrt(np.sum((self.pop_point-np.zeros((self.size,self.dim)))**2, axis =1))      
        fitness = np.exp(-(euclidean_distance**2)/2)
        return fitness


    
    def generate_ancestor(self): 
        '''Generate ancestor with the predefined initial fitness'''
   
        if self.initial_fitness == 1: #special case if want ancestor at optimum
            self.anc_point = np.zeros(self.dim) #place ancestor at optimum
            self.anc_fit = 1
        else:
            self.anc_point = np.random.normal(size=self.dim,scale=.01) 
            self.anc_fit = np.exp(-(spa.distance.euclidean(self.anc_point,np.zeros(self.dim))**2)/2) 
            scaling = 1 


            if self.anc_fit < self.initial_fitness: 

                while self.anc_fit < self.initial_fitness: 
                    scaling = scaling + .0001 #increase scalar
                    self.anc_point = self.anc_point/scaling 
                    self.anc_fit = np.exp(-(spa.distance.euclidean(self.anc_point,np.zeros(self.dim))**2)/2) #test fitness

            elif self.anc_fit > self.initial_fitness: 

                while self.anc_fit > self.initial_fitness:  #keep moving ancestor away from optimum until desired initial fitness reached
                    scaling = scaling + .0001 #increase scalar
                    self.anc_point = self.anc_point * scaling 
                    self.anc_fit = np.exp(-(spa.distance.euclidean(self.anc_point,np.zeros(self.dim))**2)/2)            
            

            
    def create_pop(self):
        """ Create population - initially monomorphic. """

        self.anc_pop_point = np.repeat([self.anc_point], self.size, axis =0)
        
        self.pop_point = np.repeat([self.anc_point], self.size, axis =0)  
        self.pop_fit = np.repeat([float(self.anc_fit)], self.size, axis =0)
        
        self.soma_site = np.ones((self.size, self.nLoci, self.ploidy), dtype ='int') # Will be used to indicate whethter mutation occurs
        # self.germ_site = np.ones((self.size, self.nLoci, 2), dtype ='int')    
        
        

    def generate_pleiotropy(self):
        
        pleiotropy_matrix = np.zeros((self.nLoci, self.dim), dtype = 'int')
        
        pleiotropy_matrix[:,:self.pleiotropy_level] =1
                        
        sum_dim = np.sum(pleiotropy_matrix, axis = 0)
        
        while np.min(sum_dim) ==0:
            for i in range(self.nLoci):
                np.random.shuffle(pleiotropy_matrix[i])            
            sum_dim = np.sum(pleiotropy_matrix, axis = 0)    
        
        return pleiotropy_matrix
    

    def generate_mut_covariance(self):

        mut_matrix = np.identity(self.dim)*((self.sigma**2)) 

        total_mut_matrix = []
        for i in range(self.nLoci):
            total_mut_matrix.append(mut_matrix*self.pleiotropy_matrix[i])

        total_mut_matrix = np.array(total_mut_matrix)

        return total_mut_matrix


    
    def mutate(self):

        soma_mut_occur = np.random.binomial(self.soma_site, self.mutation_rate)
        # germ_mut_occur = np.random.binomial(self.germ_site, self.mutation_rate)  

        total_soma_mut = []
        # total_germ_mut = []
        for i in range(self.nLoci):
            total_soma_mut.append(np.random.multivariate_normal(mean=np.zeros(self.dim), cov= self.mut_covariance[i], size = self.size*self.ploidy))
            # total_germ_mut.append(np.random.multivariate_normal(mean=np.zeros(self.dim), cov= self.mut_covariance[i], size = self.size*2))
        
        total_soma_mut = np.array(total_soma_mut)
        # total_germ_mut = np.array(total_germ_mut)


        new_soma_mut = np.split(total_soma_mut, self.size, axis =1)
        new_soma_mut = np.array(new_soma_mut)

        # new_germ_mut = np.split(total_germ_mut, self.size, axis =1)
        # new_germ_mut = np.array(new_germ_mut)        

        self.soma_mut_moving_pos[soma_mut_occur ==1] = new_soma_mut[soma_mut_occur ==1]
        # self.germ_mut_moving_pos[germ_mut_occur ==1] = new_germ_mut[germ_mut_occur ==1]


        soma_moving_pos_each_loci = np.sum(self.soma_mut_moving_pos, axis =2)
        soma_moving_pos = np.sum(soma_moving_pos_each_loci, axis =1)

        self.pop_point = self.anc_pop_point + soma_moving_pos
        self.pop_fit = self.calculate_pop_fitness()          



    def selection(self):
        '''Select the parents for next generation based on their relative fitness. Select with replacement.'''

        assert self.pop_fit.min() >= 0 # ensure that the fitness of every one within the population is non-negative.
        
        total_w = np.sum(self.pop_fit)   # First caclulate the total fitness within the population
        relfit = self.pop_fit/total_w  # Calculate the relative fitness of each individual
        csrelfit = np.cumsum(relfit)
        randvals = np.random.random(self.size)
        parents = np.searchsorted(csrelfit,randvals)
        
        return parents


 
    def mitosis(self, parents):

        '''Asexual reproduction by mitosis. Just get an copy of the selected parents (inclduing the position, the fitness of the population, and
        the moving position of mutations in both Soma and Germ .'''

             
        self.soma_mut_moving_pos = self.soma_mut_moving_pos[parents]
        # self.germ_mut_moving_pos = self.germ_mut_moving_pos[parents]
        
        # Then based on the new soma_moving_pos, recalculate the pop_point and pop_fit

        soma_moving_pos_each_loci = np.sum(self.soma_mut_moving_pos, axis =2)
        soma_moving_pos = np.sum(soma_moving_pos_each_loci, axis =1)

        self.pop_point = self.anc_pop_point + soma_moving_pos
        self.pop_fit = self.calculate_pop_fitness()    



    def amitosis(self, parents):
        '''Asexual reproduction by amitosis.'''
        
        self.soma_mut_moving_pos = self.soma_mut_moving_pos[parents] # First get the moving positions of mutations in Soma of selected parents
        new_soma_mut_moving_pos = np.repeat(self.soma_mut_moving_pos, 2, axis =2) # Duplicate the mutations in Soma of selected parents
        
        # Make sure that shuffle can actually do the shuffling. Maybe try random.choice first.
        for i in range(self.size):
            for j in range(self.nLoci):
                np.random.shuffle(new_soma_mut_moving_pos[i][j])  # Shuffle the duplicated mutations in Soma
                
        self.soma_mut_moving_pos = new_soma_mut_moving_pos[:, :, :self.ploidy]  # Pick the first self.ploidy mutations as the mutations in offspring Soma
        
        soma_moving_pos_each_loci = np.sum(self.soma_mut_moving_pos, axis =2)
        soma_moving_pos = np.sum(soma_moving_pos_each_loci, axis =1)

        self.pop_point = self.anc_pop_point + soma_moving_pos
        self.pop_fit = self.calculate_pop_fitness() 
       
        # self.germ_mut_moving_pos = self.germ_mut_moving_pos[parents]   # Get a copy of the Germ mutations in selected parents



    def asexual_reproduction(self,asex_type='amitosis'):
        """Take populations through a single Wright-Fisher model generation"""

        self.mutate()
        parents = self.selection()
        if asex_type == 'amitosis':
            self.amitosis(parents)
        elif asex_type == 'mitosis':
            self.mitosis(parents)
        self.current_step += 1




    def get_results(self):
        '''Get soma data from the simulation'''

        W = self.pop_fit


        mW = np.nanmean(W)
        log_mW = np.log(mW)

        varW = np.var(W)  # Variance of the fitness within the population.


        return [mW, log_mW, varW]



    def simulate(self,stepcount, asex_type='amitosis',sex_freq=None,random_mating=False,strides=10):

        '''Run the simulation for certain generations with defined reproduction strategy.'''
        
        """stepcount - total number of generations or time steps to run
        asex_type - type of asexual reproduction: mitosis, amitosis
        sex_freq - frequency of sexual reproduction
        random_mating - selfing or random mating
        strides - how frequently to get data using the get_results method
        
        returns results as a list of lists"""


        results = [self.get_results()]

        start = time.time()

        while self.current_step <= stepcount:
            if (sex_freq != None) and (self.current_step%sex_freq == 0):
                self.sex(random_mating)
            else:
                self.asexual_reproduction(asex_type)


            if self.current_step%strides == 0:
                results.append(self.get_results())


        colnames = ['meanFit','ln_meanFit','varFit']

        results = pd.DataFrame(np.array(results),columns=colnames)


        end = time.time()
        print 'TOTAL TIME: ',end-start
        return results



    def simulateNsave(self,outfile,stepcount,asex_type='amitosis',sex_freq=None,random_mating=False,strides=10):
        """same as simulate method except results are written to the file specified by outfile"""
        results = self.simulate(stepcount, asex_type,sex_freq,random_mating,strides)
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
            print 'Initial', data_name, len(list(set(each_gen_data)))

        elif j == m-1:
            print 'Final', data_name, len(list(set(each_gen_data)))

    return each_gen_mean_data, each_gen_std_data



def worker(n):
    '''define a function worker to conduct a single task. Here the worker was used to run a single replicate for the
    whole simulation process. Here for I am using object-oriented programming, to let the code run, the function worker 
    is used to run the object. The function worker should have parameters otherwise it will raise error when running. But 
    the parameter n here won't be used later (here n is just to make sure that the function has parameters)'''

    # np.random.seed(n)
    np.random.seed()
    
    a =Population(size =500, nLoci = 10, dim =10, pleiotropy_level=2, fit =0.5, mut_rate =0.01, sigma =0.05, ploidy =45)

    results = a.simulate(2000, asex_type='amitosis',sex_freq=None,random_mating=False,strides=1)

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
 
    fitness = [i for i in results]



    gen_mfit = get_result(fitness, 'meanFit')
    gen_mfit_mean = gen_mfit[0]
    gen_mfit_std = gen_mfit[1]

    gen_mlogfit = get_result(fitness, 'ln_meanFit') 
    gen_mlogfit_mean = gen_mlogfit[0]
    gen_mlogfit_std = gen_mlogfit[1]  

    gen_varfit = get_result(fitness, 'varFit')
    gen_varfit_mean = gen_varfit[0]
    gen_varfit_std = gen_varfit[1]

     
    fitness_result = []

    for j in range(len(gen_mfit_mean)):

        fitness_result.append([gen_mfit_mean[j], gen_mfit_std[j], gen_mlogfit_mean[j], gen_mlogfit_std[j], \
            gen_varfit_mean[j], gen_varfit_std[j]])

    colnames = ['mFit_Mean','mFit_Std', 'mlogFit_Mean','mlogFit_Std','varFit_Mean', 'varFit_Std']


    fitness_result = pd.DataFrame(np.array(fitness_result),columns=colnames)

    fitness_result.to_csv('Rv_FO_N500_L10_D10_PL2_F05_R001_SM005_Amito_P45_Cor.csv', index =False)

    end = time.time()
    print 'TOTAL TIME:', end-start