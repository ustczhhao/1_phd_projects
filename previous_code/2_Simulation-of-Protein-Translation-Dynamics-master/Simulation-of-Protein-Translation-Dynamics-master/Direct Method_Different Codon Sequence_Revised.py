# -*- coding: utf-8 -*-
"""
Created on April 24 2016 

"""


import numpy as np
from numpy import random
import math
import pandas as pd


def propensity(k, *args):
    '''Calculate the propensity of a reaction 
    
    parameters
    ----------
    k:int
      reaction constant
    *args:int(s)
      molecular number of reactant(s)
    
    Returns
    -------
    out: int
       propensity of a reaction
    '''
    
    for i in args:
        k = k*i
    return k


def ribosomebinding(a, b, c):
    '''the binding of ribosome into RNA
       Calculate the molecular number change of this reaction
       
    parameters
    ----------
    a:int
      molecular number of first reactant (free ribosome)
    b:int
      molecular number of second reactant (RBS on mRNA)
    c:int
      molecular number of product (Ribosome-RNA_RBS complex)
      
    Returns
    -------
    out: int
       molecular numbers of reactants and product after the reaction has occurred once
    '''
    
    if a > 0 and b > 0:
        a -= 1
        b -= 1
        c += 1
    return a, b, c


def releaserbs(a, b, c):
    '''Ribosome read through the second codon, release the RBS of mRNA
    Calculate the molecular number change of this reaction type
    
    parameters
    ----------
    a:int
      molecular number of first reactant (Ribosome-RNA_RBS complex)
    b:int
      molecular number of second reactant (the second codon on mRNA)
    c:int
      molecular number of product (RBS on mRNA)
      
    Returns
    -------
    out: int
       molecular numbers of reactants and product after the reaction has occurred once
    
    '''
    
    if a > 0:
        a -= 1
        b += 1
        c += 1
    return a, b, c


def readthrough(a, b):
    '''ribosome read through mRNA, actually the transition of complex
    or encounter the stop codon, release the free ribosome
    Calculate the molecular number change of this reaction type
    
    parameters
    ----------
    a:int
      molecular number of reactant (Ribosome-RNA complex[i])
    b:int
      molecular number of product (Ribosome-RNA complex[i+1] or free ribosome)
    
    Returns
    -------
    out: int
       molecular numbers of reactant and product after the reaction has occurred once
    
    '''
    
    if a > 0:
        a -= 1
        b += 1
    return a,b


def slowcodon_index(codon_no, slowcodon_no):
    '''choose the index of certain number of slow translated codons among the total number of codons
    
    parameters
    ----------
    codon_no:int
      total number of codons 
    slowcodon_no:int
      the number of slow codons
    
    Returns
    -------
    out: tuple with two elements
       first element: index of the slow codons among the total number of codons
       second element: string that represents the arrangement of slow and fast codons (fast:'0', slwo:'1')
    '''
    
    position_list = list(range(codon_no))  # get the index of all the codons
    random.seed()
    primary_list = random.choice(position_list, slowcodon_no, replace=False) # select the index of slow codons in the primary list
    index_list = sorted(primary_list)  # sorted the index of slow codons and get the sorted index_list
    
    codon_str = ''
    for i in position_list: 
        if i in index_list:  # slow codons
            codon_str += '1'
        else:                # fast codons
            codon_str += '0'
    
    return index_list, codon_str


def binary2decimal(num_str):
    '''convert a binary number into a decimal number
    
    parameter
    ----------
    num_str:str
        string that represents a binary number
    
    Return
    -------
    out: int
        the decimal number converted from the binary number
       
    '''
    
    deci_no = int(num_str, 2)  
    
    return deci_no


def decimal2binary(deci_num, codon_no):
    '''Convert a decimal number into a binary number
    
    parameters
    ----------
    deci_num:int
        the decimal number to be converted
    codon_no: int
        the total number of fast and slow codons, used to add the missing leading '0'
        when converting the decimal no. to binary no.
    
    Return
    -------
    out: str
        the string that represents the binary number converted from the decimal number
       
    '''
    
    bin_no = np.binary_repr(deci_num)
    
    for i in range(codon_no - len(bin_no)): # add the missing leading '0' to the bin_no
        bin_no = '0' + bin_no
    
    return bin_no


def get_reconstant(bin_no):
    '''get the reaction constant with certain arrangement of slow codons
    
    parameter
    ----------
    bin_no:str
      string that represnts the arrangement of fast and slow codons 
      ('0' represents fast codon, '1' represents slow codon)
   
    Returns
    -------
    out: list
        the list that stores the reaction constants with  certain arrangement of slow codons
    '''
        
        
    
    bin_list = list(bin_no)    # convert the string into list
    reconstant_list = [0.002]  # first add the reaction constant of start codon into the reconstant_list
    
    for i in bin_list:
        if i == '0':    # i=='0' means fast codon, append 12 into reconstant_list
            reconstant_list.append(12)
        elif i == '1':  # i=='1' means slow codon, append 2 into reconstant_list
            reconstant_list.append(2)
  
    return reconstant_list


def directmethod(totalcodon_no, complex_list, reconstant_list, total_rib, rna, t1, occupied_pos):
    '''Use the Direct Method to simulate one reaction step of the translation reaction in Arkin et al 
    
    parameters
    ----------
    totalcodon_no: int
        the number of all the codons (including the start codon)
    complex_list: list
        list that stores the molecular number of the transitional Ribosome-RNA[i] complex
    reconstant_list: list
        list that stores reaction constants of all the reactions (values will not change as reaction occur)
    total_rib: int
        molecular number of the free ribosome
    rna: int
        molecular number of mRNA 
    t1: int/float
        simulation initial time
    occupied_pos: list
        list that store the indexes of codons that have been occupied by Ribosome-RNA[i] complex
        
    Returns
    -------
    out: pd.DataFrame with 5 elements
       first element: list
            update the No. of Ribosome-RNA[i] complex after one step
       second element: int
            update the No. of free ribosome after one step
       third element: int
            update the No. of free mRNA-RBS after one step
       forth element: float
            update the time after one step
       fifth element: int
            the selected codon position which the translation reaction will occur
    '''     
    
    n = totalcodon_no
    
    rib = total_rib
    rna_rbs = rna
    
    prop_list = [] 
    prop1 = propensity(reconstant_list[0], rib, rna_rbs)

    prop_list.append(prop1)
     
    for i in range(1, n):    # calculate the propensity of ribosome reading through the mRNA
        prop_list.append(propensity(reconstant_list[i], complex_list[i-1]))
        
    if occupied_pos!= []:  # occupied_pos stores the occupied codon positions, exclude them by resetting their propensities to 0

        new_prop = []  # will be used to store the propensities of reactions occurred in available codon positions

        for i, p in enumerate(prop_list):   
            if i not in occupied_pos:  # if the index of element is not in occupied_pos, it means that the position is available
                new_prop.append(p) # just store the value into new_prop list

            else:                # else means the index of element is in occupied_pos, it means that the position has been occupied
                new_prop.append(0)  # just append 0 into the new_prop (means the reaction cannnot occur)
    else:
        new_prop = prop_list
        
    sumprop = sum(new_prop)
    prob_list = []
    
    if sumprop > 0:
        for i in new_prop: 
            prob_list.append(i/float(sumprop))   # calculate the probablity of each reaction
    
    if sum(prob_list)!= 0:
        u = random.choice(np.arange(0, n), p=prob_list)  # choose u 
        tau = random.exponential(1/float(sumprop))     # choose time tau
        
        if u == 0:  # u==0, means ribosome binding into mRNA, update the number of ribosome, RNA_RBS and Ribosome-RNA_RBS complex
            rib, rna_rbs, complex_list[0] = ribosomebinding(rib, rna_rbs, complex_list[0]) 
            
        elif u == 1: # u==1, means ribosome read through the second codon and release the RBS
            complex_list[0], complex_list[1], rna_rbs = releaserbs(complex_list[0], complex_list[1], rna_rbs) 
            
        elif u == n-1: # u==n-1, means the ribosome encounters the last codon
            complex_list[u-1], rib = readthrough(complex_list[u-1], rib) 
            
        else: # else means the ribosome reads through other codon positions, regarded as the transition 
            # of Ribosome-RNA complex[i]-->Ribosome-RNA complex[i+1]
            complex_list[u-1], complex_list[u] = readthrough(complex_list[u-1], complex_list[u])
        
    newcomplex_list = []  
    for i in range(len(complex_list)): # get and store the updated number of Ribosome-RNA complex[i] into newcomplex_list
        newcomplex_list.append(complex_list[i])
        
    rib1 = rib  # get the updated number of free ribosome
    rna_rbs1 = rna_rbs  # get the updated number of RNA_RBS


    data_dict = {'newcomplex_list': pd.Series([newcomplex_list]), 'rib1': pd.Series([rib1]), \
             'rna_rbs1': pd.Series([rna_rbs1]), 'tau': pd.Series([tau]), 'u': pd.Series([u])}  
    
    data = pd.DataFrame.from_dict(data_dict, orient='columns')
    
    return data


def simulation(reconstant_list, total_rib, t1):
    '''Stimulate the translation reactions to translate 100 proteins
    
    parameters
    ----------
    reconstant_list: list
        the list that stores the reaction constant of all the codons (including the start codon)
    total_rib: int
        initial value of the free ribosome
    t1: int/float
        simulation initial time
    
    Returns
    -------
    out: pd. DataFrame with 4 elements
    
       first element:the dictionary that store the time and positions of each individual reaction  
       second element: the list that stores the generated u in each step
       third element: the time that needed to translate 100 proteins
       fourth element: how many stalling times in translating 100 proteins
    '''
    
    t0 = t1    # use t0 to store the initial value of t1 (stimulation starting time)
    t1_list = []  # will be used to store the updating reaction time t1
    
    n = len(reconstant_list)   # get the no. of codons (including the start codon)
        
    complex_list = [0]*(n-1)  # set the initial value of complex_list (all the elements are 0)
    
    rib = total_rib  # get the initial number of free ribosome
    
    rna_rbs = 1  # the initial number of mRNA-RBS is 1 as only 1 mRNA is involved in the simulation
    
    result = {}   # will be used to store the time and positions of each individual reaction
    
    ids = 0  # will be used as the ID of new reaction (also as the name of reaction)
    u_list = []   # will be used to store the generated u at each step
    
    occupied_pos = []  # will be used to store the codon postions that has been occupied (reaction cannot occur)

    number = 0   # the number of proteins that have been translated completely
    
    stall_count = 0   # the number of stalling times
    
    
    while number < 100:  

        y = directmethod(n, complex_list, reconstant_list, rib, rna_rbs, t1, occupied_pos)  # run direct method
        u = list(y.u)[0]    # get u
            
        if u_list ==  []:    # at first u_list is empty, 
            complex_list = list(y.newcomplex_list)[0]   # update the molecular number of Rib-RNA complex[i], free ribosome and RNA_RBS
            rib = list(y.rib1)[0]
            rna_rbs = list(y.rna_rbs1)[0]
            
            t1 += list(y.tau)[0]    # update the reaction time t1
            t1_list.append(t1)  # store the updated time   
            
            u_list.append(u)  
            
            ids += 1  # start a new reaction (ids will be used as the name of new reaction)
            result[ids] = [(t1,u)]  # store the tuple of (reaction time, u) into dictionary (the ids will be the key)
            
            occupied_pos = [u]  # update the occupied_pos
            
        else:  # then u_list is not empty 
            complex_list = list(y.newcomplex_list)[0]  # update number of Ribosome-RNA[i], free ribosome and RNA-RBS
            rib = list(y.rib1)[0]
            rna_rbs = list(y.rna_rbs1)[0]

            u_list.append(u)
            
            if u == 0:   # if u==0, means a new peptide translation reaction has occurred
                ids += 1   # use the new ids as the name of this new reaction
                t1 += list(y.tau)[0]  
                t1_list.append(t1)  
                result[ids] = [(t1,u)]  
                    
            else:  # if u!=0, means one of the current peptide reaction is continuing
                t1 += list(y.tau)[0]  # updating the reaction time
                t1_list.append(t1)  
                for k, v in result.items(): # then check which reaction is continuing,
                    if v[-1][-1] == u-1: 
                        result[k].append((t1, u))  
                        break                  
                        
            occupied_pos = [v[-1][-1] for v in result.values() if v[-1][-1]!= n-1]
            occupied_pos = [v - 9 for v in occupied_pos if v > 8]+[v - 8 for v in occupied_pos if v > 7]+\
            [v - 7 for v in occupied_pos if v > 6]+[v - 6 for v in occupied_pos if v > 5]+[v - 5 for v in occupied_pos if v > 4]+\
            [v - 4 for v in occupied_pos if v > 3]+[v - 3 for v in occupied_pos if v > 2]+[v - 2 for v in occupied_pos if v > 1]+\
            [v - 1 for v in occupied_pos if v > 0] + occupied_pos  
            
            if u == n-1:  # if u==n-1, it means the one protein has been translated completely
                number += 1  # update the number of proteins that have been translated completely.
            else:
                pass
            
            if u+1 in occupied_pos: # if u+1 in occupied_pos, it means the next position for current selected ribosome has been
                stall_count += 1   # occupied by another ribosome, current selected ribosome cannot move forward-->stalling
   

    data_dict = {'result': pd.Series([result]), 'u_list': pd.Series([u_list]), \
             't1': pd.Series([t1]), 'stall_count': pd.Series([stall_count])}
    
    data = pd.DataFrame.from_dict(data_dict, orient='columns')  


    return data



def calculate_oneprotein(reconstant_list, total_rib, exclude_pr, t1):
    '''Calculate the time needed to translate one protein
    
    parameters
    ----------
    
    reconstant_list: list
        the list that stores the reaction constant of all the codons (including the start codon)
    total_rib: int
        initial value of the free ribosome
    exclude_pr: int
        the number of proteins that will be excluded from calculation of mean translation time
        (they are regarded as not stably translated)
    t1: int/float
        simulation initial time
    
    Returns
    -------
    out: pd. DataFrame with 5 elements
        first element: dictionary that stores the time needed to translate one protein (key is the protein name/number)
        second element: list that stores the time needed to translate one protein (with the order of protein 0, 1, 2, 3,... ) 
        third element: list that stores the time needed to translate one protein from the 20th protein (based on previous result
                       from this protein the translation time seems stable--> stably translation)
        fourth element: the total time needed to translate 80 proteins from 20th protein(stably translation)
        fifth element: the total number of stalling times during the translation of 100 proteins
    '''
    
    y = simulation(reconstant_list, total_rib, t1)  # run the simulation
    
    stall_count = list(y.stall_count)[0]   # get the stalling times 
    
    protein_dic = {}  # creat a dictionary to store the time to completely translate a protein (only the protein that has been 
                    # completely translated will be added into this dictionary)

    for k, v in list(y.result)[0].items():
        if v[-1][-1] == len(reconstant_list)-1:     # v[-1][-1]==len(reconstant_list)-1 means this protein has been translated completely
            protein_dic[k] = v[-1][0]-v[0][0]    # get the time needed to completely translate this protein
         
    keylist = list(protein_dic.keys())   # as dictionary has no order, get the protein_dic keys into a list and sort it
    keylist.sort()
    
    time_list = []
    
    for k in keylist: # store the time needed to completely translate a protein according to the order of index(protein 1,2,3...)
        time_list.append(protein_dic[k])
    
    stable_list = time_list[exclude_pr:] # get the time to stably translate a complete protein
    
    sum_stable = sum(stable_list)


    data_dict = {'protein_dic': pd.Series([protein_dic]), 'time_list': pd.Series([time_list]), \
             'stable_list': pd.Series([stable_list]), 'sum_stable': pd.Series([sum_stable]),\
               'stall_count': pd.Series([stall_count])}
    
    data = pd.DataFrame.from_dict(data_dict, orient='columns') 
                
    return data
                
    
    
def repeat_calculation(reconstant_list, total_rib, exclude_pr, t1, repeat_time):
    '''Repeat the calculation of the time needed to translate one protein for certain times
    
    parameters
    ----------
    
    reconstant_list: list
        the list that stores the reaction constant of all the codons (including the start codon)
    total_rib: int
        initial value of the free ribosome
    exclude_pr: int
        the number of proteins that will be excluded from calculation of mean translation time
        (they are regarded as not stably translated)
    t1: int/float
        simulation initial time
    repeat_time: int
        how many times we would like to repeat for the calculation
    
    Returns
    -------
    out: pd. DataFrame with 6 elements
        first element: list that stores the mean time needed to translate one protein 
        second element: list that stores the standard deviation of time needed to translate one protein
        third element: list that stores the standard error of time 
        fourth element: list that stores the mean no. of stalling during translation of 100 proteins
        fifth element: list that stores the standard deviation of stalling during translation of 100 proteins
        sixth element: list that stores the standard error of stalling during translation of 100 proteins
    '''
    
    total_time_list = []  
    total_sum = []
    
    stall_list = []
    
    for i in range(repeat_time):   # for a certain reaction constant list, repeat running for repeat_time times
        y = calculate_oneprotein(reconstant_list, total_rib, exclude_pr, t1)
        total_time_list.append(list(y.stable_list)[0])
        total_sum.append(list(y.sum_stable)[0])
        stall_list.append(list(y.stall_count)[0])
    

    pro_num = 100-exclude_pr   # get the number of the stably translated protein

    s_list = []
    for j in range(pro_num): 
        for i in total_time_list:
            s_list.append(i[j])
         
    mean_time = sum(s_list)/(pro_num*repeat_time)
    std_time = np.std(s_list, ddof=1)
    se_time = std_time/math.sqrt(pro_num*repeat_time)
    
    mean_stall = sum(stall_list)/(repeat_time)
    std_stall = np.std(stall_list, ddof=1)
    se_stall = std_stall/math.sqrt(repeat_time)
            

    data_dict = {'mean_time': pd.Series([mean_time]), 'std_time': pd.Series([std_time]), \
             'se_time': pd.Series([se_time]), 'mean_stall': pd.Series([mean_stall]),\
               'std_stall': pd.Series([std_stall]), 'se_stall': pd.Series([se_stall])}
    
    data = pd.DataFrame.from_dict(data_dict, orient='columns') 


    return data
        

'''Below is an example: the total no. of fast and slow codons is 100 (20 of them are slow codons), we would like to simulate the 
translation pf 100 proteins with 50  ribosomes. For the arrangement of slow codons, we try to generate 5 different arrangement, and for each 
arrangement, the simulation process is repeated 3 times to get the mean and standard error of time needed to stably tranlated one
protein (the first 20 proteins are excluded as they are not stably translated). Also, the mean and standard error of stalling times during
the translation of 100 proteins are also recorded.'''

codon_no = 100      # set the initial value of the total no. of all the codons (exclude the start codon)
slowcodon_no = 20   # the no. of slow codons 
exclude_pr = 20    # how many proteins we would like to exclude as non-stably translated protein
rib = 50   # the no. of ribosomes 
t1 = 0   # the simulation starting time
repeat_time = 3   # how many times we would like to repeat for a certain arrangement 

total_deci_no = []
for i in range(5):  # generate 5 random slow codon arrangements
    y = slowcodon_index(codon_no, slowcodon_no)
    deci_no = binary2decimal(y[1])  # store the arrangements by converting them into decimal numbers
    total_deci_no.append(deci_no)


total_reconstant_list = []
for j in total_deci_no:   # get the reconstant_list by converting the decimal number
    z = decimal2binary(j, codon_no)     
    reconstant_list = get_reconstant(z)    
    total_reconstant_list.append(reconstant_list)

o = open('Sim_result_py2.txt', 'w')

o.write('FINAL TOTAL RECONS: {}\n'.format(total_reconstant_list))


total_mean_time = []
total_se_time = []
total_mean_stall = []
total_se_stall = []


for i in total_reconstant_list:
    z = repeat_calculation(i, rib, exclude_pr, t1, repeat_time)   # run the repeat_calculation with each arrangement of slow codons
    total_mean_time.append(list(z.mean_time)[0])
    total_se_time.append(list(z.se_time)[0])
    total_mean_stall.append(list(z.mean_stall)[0])
    total_se_stall.append(list(z.se_stall)[0])


    
o.write('TOTAL M TIME: {},{}\n'. format(total_mean_time, len(total_mean_time)))
o.write('TOTAL SE TIME: {},{}\n'. format(total_se_time, len(total_se_time)))
o.write('TOTAL M STALL: {},{}\n'. format(total_mean_stall, len(total_mean_stall)))
o.write('TOTAL SE STALL: {},{}\n'. format(total_se_stall, len(total_se_stall)))
    

min_time = min(total_mean_time)        # get the min and max mean time to stably translate one protein
max_time = max(total_mean_time)

min_time_index = total_mean_time.index(min_time)   
max_time_index = total_mean_time.index(max_time)

min_time_reconstant = total_reconstant_list[min_time_index]   # get the slow codon arrangement corresponding to the min and max mean time to 
max_time_reconstant = total_reconstant_list[max_time_index]   # stably translate one protein


min_stall = min(total_mean_stall)    # get the min and max mean stalling times to translate 100 proteins
max_stall = max(total_mean_stall)

min_stall_index = total_mean_stall.index(min_stall)
max_stall_index = total_mean_stall.index(max_stall)

min_stall_reconstant = total_reconstant_list[min_stall_index]  # get the slow codon arrangement corresponding to the min and max mean stalling
max_stall_reconstant = total_reconstant_list[max_stall_index]  # times to translate 100 proteins

o.write('MIN TIME: {}\n'. format(min_time))
o.write('MIN TIME INDEX: {}\n'. format(min_time_index))
o.write('MIN TIME RECON: {}\n'. format(min_time_reconstant))

o.write('MAX TIME: {}\n'. format(max_time))
o.write('MAX TIME INDEX: {}\n'. format(max_time_index))
o.write('MAX TIME RECON: {}\n'. format(max_time_reconstant))

o.write('MIN STALL: {}\n'. format(min_stall))
o.write('MIN STALL INDEX: {}\n'. format(min_stall_index))
o.write('MIN STALL RECON: {}\n'. format(min_stall_reconstant))

o.write('MAX STALL: {}\n'. format(max_stall))
o.write('MAX STALL INDEX: {}\n'. format(max_stall_index))
o.write('MAX STALL RECON: {}\n'. format(max_stall_reconstant))

o.close()