"""
# -*- coding: utf-8 -*-

Created on April 21, 2016 

Python 3

"""



import numpy as np
from numpy import random
import matplotlib.pyplot as plt
import math
%matplotlib inline


def propensity(k,*args):
    '''Calculate the propersity of a reaction 
    
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
        k=k*i
    return k


def ribosomebinding(a,b,c):
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
    
    if a>0 and b>0:
        a-=1
        b-=1
        c+=1
    return a,b,c


def releaserbs(a,b,c):
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
    
    if a>0:
        a-=1
        b+=1
        c+=1
    return a,b,c


def readthrough(a,b):
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
    if a>0:
        a-=1
        b+=1
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
    out: list
       index of the slow codons among the total number of codons
    '''
    
    position_list=list(range(codon_no)) # get the index of all the codons
    random.seed()
    primary_list=random.choice(position_list, slowcodon_no, replace=False) # select the index of slow codons in the primary list
    index_list=sorted(primary_list) # sorted the index of slow codons and get the sorted index_list
    
    return index_list


def get_reconstant(codon_no, slowcodon_no):
    '''get the reaction constant with certain arrangement of slow codons
    
    parameters
    ----------
    codon_no:int
      total number of codons (do not include the start codon)
    slowcodon_no:int
      the number of slow codons
    
    Returns
    -------
    out: tuple with two elements
        First element: the list that stores the reaction constants with  certain arrangement of slow codons
        Second element: the list that stores the primary index of slow codons among all the codons
    '''
    
    reconstant_list=[0.002]  # the first element in reaction constant list is the reaction constant for start codon
    for i in range(codon_no):
        reconstant_list.append(12) # first assume that all the codons are fast codons(k=12), and store the k into reconstant_list
    
    index_list=slowcodon_index(codon_no, slowcodon_no) # get the index of slow codons
    
    for j in index_list:  # change the k of slow codons to 2 according to the index_list (positions of slow codons)
        reconstant_list[j+1]=2 # use j+1 because the first element in reconstant_list is k of start codon
        
    return reconstant_list, index_list


def get_reconstant1(codon_no, slowcodon_no):
    '''get the reaction constant with certain arrangement of slow codons: The first 20 codons are slow codons
    (for the documentation, here we assume slowcodon_no=20)
    
    parameters
    ----------
    codon_no:int
      total number of codons (do not include the start codon)
    slowcodon_no:int
      the number of slow codons
    
    Returns
    -------
    out: list
        the list that stores the reaction constants with  certain arrangement of slow codons
    '''
    
    reconstant_list=[0.002]
    for i in range(codon_no):
        reconstant_list.append(12)
    
    for j in range(slowcodon_no):
        reconstant_list[j+1]=2
        
    return reconstant_list


def get_reconstant2(codon_no, slowcodon_no):
    '''get the reaction constant with certain arrangement of slow codons: The last 20 codons are slow codons
    (for the documentation, here we assume slowcodon_no=20)
    
    parameters
    ----------
    codon_no:int
      total number of codons (do not include the start codon)
    slowcodon_no:int
      the number of slow codons
    
    Returns
    -------
    out: list
        the list that stores the reaction constants with  certain arrangement of slow codons
    '''
    
    reconstant_list=[0.002]
    for i in range(codon_no):
        reconstant_list.append(12)
        
    for j in range(codon_no-slowcodon_no,codon_no):
        reconstant_list[j+1]=2
        
    return reconstant_list


def get_reconstant_n(codon_no, slowcodon_no,n):
    '''get the reaction constant with certain arrangement of slow codons: from position 'n' the following 
    slowcodon_no are slow codons
    
    parameters
    ----------
    codon_no:int
      total number of codons (do not include the start codon)
    slowcodon_no:int
      the number of slow codons
    n:int
      the position that the slowcodon_no-slow-codon-cluster begins (from this position the following slowcodon_no 
      are slow codons)  
    
    Returns
    -------
    out: list
        the list that stores the reaction constants with  certain arrangement of slow codons
    '''
    
    reconstant_list=[0.002]
    for i in range(codon_no):
        reconstant_list.append(12)
        
    for j in range(n,n+slowcodon_no):
        reconstant_list[j+1]=2
        
    return reconstant_list
    
    
def get_reconstant3(codon_no, slowcodon_no):
    '''get the reaction constant with certain arrangement of slow codons: 
    The 20 slow codons distribute evenly: 4 fast codons followed by 1 slow codon
    (for the documentation, here we assume slowcodon_no=20)
    
    parameters
    ----------
    codon_no:int
      total number of codons (do not include the start codon)
    slowcodon_no:int
      the number of slow codons
    
    Returns
    -------
    out: list
        the list that stores the reaction constants with  certain arrangement of slow codons
    '''
    
    reconstant_list=[0.002]
    for i in range(codon_no):
        reconstant_list.append(12)
    
    for j in range(1,21):
        reconstant_list[5*j]=2
        
    return reconstant_list


def get_reconstant4(codon_no, slowcodon_no):
    '''get the reaction constant with certain arrangement of slow codons: 
    The 20 slow codons distribute evenly: 3 fast codons followed by 1 slow codon then another fast codon
    (for the documentation, here we assume slowcodon_no=20)
    
    parameters
    ----------
    codon_no:int
      total number of codons (do not include the start codon)
    slowcodon_no:int
      the number of slow codons
    
    Returns
    -------
    out: list
        the list that stores the reaction constants with  certain arrangement of slow codons
    '''
    
    reconstant_list=[0.002]
    for i in range(codon_no):
        reconstant_list.append(12)
    
#     index_list=slowcodon_index(codon_no, slowcodon_no)
    
    for j in range(1,21):
        reconstant_list[5*j-1]=2
        
    return reconstant_list


def get_reconstant5(codon_no, slowcodon_no):
    '''get the reaction constant with certain arrangement of slow codons: 
    The 20 slow codons distribute evenly: 2 fast codons followed by 1 slow codon then another 2 fast codons
    (for the documentation, here we assume slowcodon_no=20)
    
    parameters
    ----------
    codon_no:int
      total number of codons (do not include the start codon)
    slowcodon_no:int
      the number of slow codons
    
    Returns
    -------
    out: list
        the list that stores the reaction constants with  certain arrangement of slow codons
    '''
    
    reconstant_list=[0.002]
    for i in range(codon_no):
        reconstant_list.append(12)
    
    for j in range(1,21):
        reconstant_list[5*j-2]=2
        
    return reconstant_list


def get_reconstant6(codon_no, slowcodon_no):
    '''get the reaction constant with certain arrangement of slow codons: 
    The 20 slow codons distribute evenly: 1 fast codons followed by 1 slow codon then another 3 fast codons
    (for the documentation, here we assume slowcodon_no=20)
    
    parameters
    ----------
    codon_no:int
      total number of codons (do not include the start codon)
    slowcodon_no:int
      the number of slow codons
    
    Returns
    -------
    out: list
        the list that stores the reaction constants with  certain arrangement of slow codons
    '''
    
    reconstant_list=[0.002]
    for i in range(codon_no):
        reconstant_list.append(12)
    
    for j in range(1,21):
        reconstant_list[5*j-3]=2
        
    return reconstant_list


def get_reconstant7(codon_no, slowcodon_no):
    '''get the reaction constant with certain arrangement of slow codons: 
    The 20 slow codons distribute evenly: 1 slow codon followed by 4 fast codons
    (for the documentation, here we assume slowcodon_no=20)
    
    parameters
    ----------
    codon_no:int
      total number of codons (do not include the start codon)
    slowcodon_no:int
      the number of slow codons
    
    Returns
    -------
    out: list
        the list that stores the reaction constants with  certain arrangement of slow codons
    '''
    
    reconstant_list=[0.002]
    for i in range(codon_no):
        reconstant_list.append(12)
    
    for j in range(1,21):
        reconstant_list[5*j-4]=2
        
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
    out: tuple with 5 elements
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
    n=totalcodon_no
    
    rib=total_rib
    rna_rbs=rna
    
    prop_list=[] 
    prop1=propensity(reconstant_list[0],rib,rna_rbs)

    prop_list.append(prop1)
     
    for i in range(1,n):   # calculate the propensity of ribosome reading through the mRNA
        prop_list.append(propensity(reconstant_list[i],complex_list[i-1]))
        
    if occupied_pos!=[]:  # occupied_pos stores the occupied codon positions, exclude them by resetting their propensities to 0
        new_prop = []   # will be used to store the propensities of reactions occurred in available codon positions
        for i, p in enumerate(prop_list):   # get the indexes and values of elements in prop_list
            if i not in occupied_pos:  # if the index of element is not in occupied_pos, it means that the position is available
                new_prop.append(p)  # just store the value into new_prop list
            else:   # else means the index of element is in occupied_pos, it means that the position has been occupied
                new_prop.append(0)  # just append 0 into the new_prop (means the reaction cannnot occur)
    else:
        new_prop = prop_list
        
    sumprop=sum(new_prop)
    prob_list=[]
    
    if sumprop > 0:
        for i in new_prop: 
            prob_list.append(i/sumprop)  # calculate the probablity of each reaction
    
    if sum(prob_list)!=0:
        u=random.choice(np.arange(0,n),p=prob_list)  # choose u 
        tau=random.exponential(1/sumprop)    # choose time tau
        
        if u==0: # u==0, means ribosome binding into mRNA, update the number of ribosome, RNA_RBS and Ribosome-RNA_RBS complex
            rib,rna_rbs,complex_list[0]=ribosomebinding(rib,rna_rbs,complex_list[0]) 
            
        elif u==1: # u==1, means ribosome read through the second codon and release the RBS
            complex_list[0],complex_list[1],rna_rbs=releaserbs(complex_list[0],complex_list[1],rna_rbs) 
            
        elif u==n-1: # u==n-1, means the ribosome encounters the STOP or the last codon
            complex_list[u-1],rib=readthrough(complex_list[u-1],rib) 
            
        else:# else means the ribosome reads through other codon positions, regarded as the transition 
            # of Ribosome-RNA complex[i]-->Ribosome-RNA complex[i+1]
            complex_list[u-1],complex_list[u]=readthrough(complex_list[u-1],complex_list[u])
        
    newcomplex_list=[]  
    for i in range(len(complex_list)): # get and store the updated number of Ribosome-RNA complex[i] into newcomplex_list
        newcomplex_list.append(complex_list[i])
        
    rib1=rib  # get the updated number of free ribosome
    rna_rbs1=rna_rbs # get the updated number of RNA_RBS
    
    return newcomplex_list, rib1, rna_rbs1, tau, u


def simulation(codon_no,slowcodon_no, reconstant_list, total_rib,t1):
    '''Stimulate the translation reactions to translate 100 proteins
    
    parameters
    ----------
    codon_no: int
        the total number of fast and slow codons (do not include the start codon)
    slowcodon_no: int
        the number of slow codons
    reconstant_list: list
        the list that stores the reaction constant of all the codons (including the start codon)
    total_rib: int
        initial value of the free ribosome
    t1: int/float
        simulation initial time
    
    Returns
    -------
    out: tuple with 4 elements
       first element:the dictionary that store the time and positions of each individual reaction  
       second element: the list that stores the generated u in each step
       third element: the time that needed to translate 100 proteins
       fourth element: how many stalling times in translating 100 proteins
    '''
    
    t0=t1  # use t0 to store the initial value of t1 (stimulation starting time)
    t1_list=[]  # will be used to store the updating reaction time t1
    
    n=codon_no+1
    
    complex_list=[0]*(n-1)  # set the initial value of complex_list (all the elements are 0)
    
    rib=total_rib # get the initial number of free ribosome
    
    rna_rbs=1  # the initial number of mRNA-RBS is 1 as only 1 mRNA is involved in the simulation
    
    result={}  # will be used to store the time and positions of each individual reaction
    
    ids=0  # will be used as the ID of new reaction (also as the name of reaction)
    u_list=[]   # will be used to store the generated u at each step
    
    occupied_pos=[]  # will be used to store the codon postions that has been occupied (reaction cannot occur)

    number=0   # the number of proteins that have been translated completely
    
    stall_count=0 # will be used to store the stalling time during translation
    
    
    while number<100:  

        y=directmethod(n, complex_list, reconstant_list, rib, rna_rbs, t1, occupied_pos)  # run direct method
        u=y[-1]  # get u
    
        if u_list==[]:    # at first u_list is empty, 
            complex_list=y[0]  # update the molecular number of Rib-RNA complex[i], free ribosome and RNA_RBS
            rib=y[1]
            rna_rbs=y[2]
            
            t1+=y[3]  # update the reaction time t1
            t1_list.append(t1)  # store the updated time 
            
            u_list.append(u)  # store the u generated from this step (then u_list is not empty)
            
            ids+=1  # start a new reaction (ids will be used as the name of new reaction)
            result[ids]=[(t1,u)]  # store the tuple of (reaction time, u) into dictionary (the ids will be the key)
            
            occupied_pos=[u]  # update the occupied_pos
            
        else:  # then u_list is not empty 
            complex_list=y[0]   # update number of Ribosome-RNA[i], free ribosome and RNA-RBS
            rib=y[1]
            rna_rbs=y[2]

            u_list.append(u)
            
            if u==0:   # if u==0, means a new peptide translation reaction has occurred
                ids+=1   # use the new ids as the name of this new reaction
                t1 += y[3]  # update the reaction time
                t1_list.append(t1)  # store the updated reaction time t1 into t1_list
                result[ids]=[(t1,u)]  # store the tuple of (reaction time, u) into dictionary (the new ids will be the key)
                    
            else:  # if u!=0, means one of the current peptide reaction is continuing
                t1 += y[3]  # updating the reaction time
                t1_list.append(t1)  # store the updated reaction time t1 into t1_list
                for k, v in result.items(): # then check which reaction is continuing,
                    if v[-1][-1]==u-1: # v[-1][-1] is the last codon position of each individual reaction, if v[-1][-1]==u-1,
                        result[k].append((t1, u)) # it means this peptide translation reaction is occuring, just append the 
                                                 # tuple of (reaction time, u) into result dictionary under this key, and also
                        break                   # append the translated amino acid into aa_result dictionary under this key
                        
            occupied_pos = [v[-1][-1] for v in result.values() if v[-1][-1]!=n-1]
            occupied_pos = [v - 9 for v in occupied_pos if v>8]+[v - 8 for v in occupied_pos if v>7]+\
            [v - 7 for v in occupied_pos if v>6]+[v - 6 for v in occupied_pos if v>5]+[v - 5 for v in occupied_pos if v>4]+\
            [v - 4 for v in occupied_pos if v>3]+[v - 3 for v in occupied_pos if v>2]+[v - 2 for v in occupied_pos if v>1]+\
            [v - 1 for v in occupied_pos if v>0] + occupied_pos # here we assumed the footprint is 10 codons 
                                                                #(from literature:28bp)
            
            if u==n-1:  # if u==n-1, it means the one protein has been translated completely
                number+=1  # update the number of proteins that have been translated completely.
            else:
                pass
            
            if u+1 in occupied_pos: # if u+1 in occupied_pos, it means the next position for current selected ribosome has been 
                stall_count+=1      # occupied by another ribosome, current selected ribosome cannot move forward-->stalling
   
    return result, u_list, t1, stall_count



def calculate_oneprotein(codon_no,slowcodon_no, reconstant_list, total_rib,t1):
    '''Calculate the time needed to translate one protein
    
    parameters
    ----------
    codon_no: int
        the total number of fast and slow codons (do not include the start codon)
    slowcodon_no: int
        the number of slow codons
    reconstant_list: list
        the list that stores the reaction constant of all the codons (including the start codon)
    total_rib: int
        initial value of the free ribosome
    t1: int/float
        simulation initial time
    
    Returns
    -------
    out: tuple with 5 elements
        first element: dictionary that stores the time needed to translate one protein (key is the protein name/number)
        second element: list that stores the time needed to translate one protein (with the order of protein 0, 1, 2, 3,... ) 
        third element: list that stores the time needed to translate one protein from the 20th protein (based on previous result
                       from this protein the translation time seems stable--> stably translation)
        fourth element: the total time needed to translate 80 proteins from 20th protein(stably translation)
        fifth element: the total number of stalling times during the translation of 100 proteins
    '''
    
    y=simulation(codon_no,slowcodon_no, reconstant_list, total_rib,t1)  # run the simulation
    
    stall_count=y[-1]  # get the stalling time during the translation of 100 proteins
    
    protein_dic={}  # creat a dictionary to store the time to completely translate a protein (only the protein that has been 
                    # completely translated will be added into this dictionary)

    for k, v in y[0].items():
        if v[-1][-1]==codon_no:  # this means that this protein has been translated completely
            protein_dic[k]=v[-1][0]-v[0][0] # get the time to completely translate the protein (finishing time- initial time)
         
    keylist=list(protein_dic.keys()) # as dictionary has no order, get the protein_dic keys into a list and sort it
    keylist.sort()  # it will appear according to the order of index (protein 1, 2, 3...)
    
    time_list=[]  # will be used to store the time needed to completely translate a protein
    
    for k in keylist: # store the time needed to completely translate a protein according to the order of index(protein 1,2,3...)
        time_list.append(protein_dic[k])
        
    stable_list=time_list[20:] # get the time to stably translate a complete protein (based on previous result, this will be
                               # from 20th protein the translation time seems stable--> stably translation)
    
    sum_stable=sum(stable_list)  # get the total time to stably translate 20th ~100th protein (totally 80 proteins)
                
    return protein_dic, time_list, stable_list, sum_stable, stall_count
                
    
    
def repeat_calculation(codon_no,slowcodon_no, reconstant_list, total_rib,t1, repeat_time):
    '''Repeat the calculation of the time needed to translate one protein for certain times
    
    parameters
    ----------
    codon_no: int
        the total number of fast and slow codons (do not include the start codon)
    slowcodon_no: int
        the number of slow codons
    reconstant_list: list
        the list that stores the reaction constant of all the codons (including the start codon)
    total_rib: int
        initial value of the free ribosome
    t1: int/float
        simulation initial time
    repeat_time: int
        how many times we would like to repeat for the calculation
    
    Returns
    -------
    out: tuple with 6 elements
        first element: list that stores the mean time needed to translate one protein 
        second element: list that stores the standard deviation of time needed to translate one protein
        third element: list that stores the standard error of time 
        fourth element: list that stores the mean no. of stalling during translation of 100 proteins
        fifth element: list that stores the standard deviation of stalling during translation of 100 proteins
        sixth element: list that stores the standard error of stalling during translation of 100 proteins
    '''
    
    total_time_list=[] 
    total_sum=[]
    
    stall_list=[]
    
    for i in range(repeat_time): # for a certain reaction constant list, repeat running for repeat_time times
        y=calculate_oneprotein(codon_no,slowcodon_no, reconstant_list, total_rib,t1) 
        total_time_list.append(y[2]) # store the time needed to stably translate one protein 
        total_sum.append(y[-2]) # store the total time needed to translate 80 proteins from 20th protein(stably translation)
        stall_list.append(y[-1]) # store the the total number of stalling times during the translation of 100 proteins in each 
                                    # repeat times
        
    s_list=[]  # will be used to store every element in total_time_list (each repeat time, get the time needed to             
    for j in range(80):  # translate one stably translated protein and store it into s_list) 
        for i in total_time_list:
            s_list.append(i[j]) 
            
    mean_time=sum(s_list)/(80*repeat_time)  # get the mean time needed to stably translate one protein
    std_time=np.std(s_list, ddof=1)  # get the STD of time needed to stably translate one protein
    se_time=std_time/math.sqrt(80*repeat_time) # get the standard error of time needed to stably translate one protein
    
    mean_stall=sum(stall_list)/(repeat_time) # get the mean no. of stalling times needed to translate 100 proteins
    std_stall=np.std(stall_list, ddof=1) # get the STD of stalling times needed to translate 100 proteins
    se_stall=std_stall/math.sqrt(repeat_time) # get the standard error of stalling times needed to translate 100 proteins
            
        
    return mean_time,std_time,se_time, mean_stall, std_stall, se_stall
        
