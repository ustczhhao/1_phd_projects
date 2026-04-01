"""
# -*- coding: utf-8 -*-

Created on April 12 2016 

Python 3

"""



import numpy as np
from numpy import random
import statistics as st
import matplotlib.pyplot as plt
%matplotlib inline


# footprint of ribosome is 2 codons, and the second codon occupied by the ribosome is translated.

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


def translation(sequence):  
    '''Given a mRNA sequence, translate the sequence into amino acid sequence
    
    parameters
    ----------
    sequence:str
      the given mRNA sequence
    
    Returns
    -------
    out: tuple with 2 elements
        first element: int
            how many amino acids would be translated from the given mRNA sequence
        second element: list
            store the amino acid sequence translated from the given mRNA sequence
    '''
    
    codon_dict={'UUU':'F', 'UUC':'F', 'UUA':'L', 'UUG':'L',     
                 'UCU':'S', 'UCC':'S', 'UCA':'S', 'UCG':'S',     
                 'UAU':'Y', 'UAC':'Y', 'UAA':'STOP', 'UAG':'STOP',     
                 'UGU':'C', 'UGC':'C', 'UGA':'STOP', 'UGG':'W',     
                 'CUU':'L', 'CUC':'L', 'CUA':'L', 'CUG':'L',     
                 'CCU':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',     
                 'CAU':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',     
                 'CGU':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',     
                 'AUU':'I', 'AUC':'I', 'AUA':'I', 'AUG':'M',     
                 'ACU':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',     
                 'AAU':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',     
                 'AGU':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',     
                 'GUU':'V', 'GUC':'V', 'GUA':'V', 'GUG':'V',     
                 'GCU':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',     
                 'GAU':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',    
                 'GGU':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
    
    aa_list=[]    # will be used to store the amino acid sequence translated from the mRNA sequence
    
    start_point=sequence.find('AUG')    # search the start codon in the mRNA sequence
    trimmed_sequence=sequence[start_point:]  # trim the mRNA sequence from the start codon
    
    codons_list = [trimmed_sequence[i:i+3] for i in range(0, len(trimmed_sequence), 3)]  # trim the sequence into codon list
    
    stop_list=['UAA', 'UAG','UGA']  # stop codon list
    
    n=0
    for i in codons_list:
        if len(i)==1 or len(i)==2:  # means the length of last element is less than 3, not a codon, just pass
            pass
        elif i not in stop_list:   # the codon is not stop codon, translate it into amino acid, and store into aa_list
            n+=1
            aa_list.append(codon_dict[i])
        else:     # the codon is stop codon, break
            n+=1
            aa_list.append(codon_dict[i])
            break
            
    return n,aa_list
    

def calculation_constant(sequence):
    '''Given a mRNA sequence, get the reaction constant of each codon in the sequence
    
    parameters
    ----------
    sequence:str
      the given mRNA sequence
    
    Returns
    -------
    out: list 
        list that store the reaction constant of each codon in the mRNA sequence
    '''
    
    codonconstant_dict={'UUU':33.3, 'UUC':33.3, 'UUA':33.3, 'UUG':33.3,     
                        'UCU':33.3, 'UCC':33.3, 'UCA':33.3, 'UCG':33.3,     
                        'UAU':33.3, 'UAC':33.3, 'UAA':33.3, 'UAG':33.3,     
                        'UGU':33.3, 'UGC':33.3, 'UGA':33.3, 'UGG':33.3,     
                        'CUU':33.3, 'CUC':33.3, 'CUA':33.3, 'CUG':33.3,     
                        'CCU':33.3, 'CCC':33.3, 'CCA':33.3, 'CCG':33.3,     
                        'CAU':33.3, 'CAC':33.3, 'CAA':33.3, 'CAG':33.3,     
                        'CGU':33.3, 'CGC':33.3, 'CGA':33.3, 'CGG':33.3,     
                        'AUU':33.3, 'AUC':33.3, 'AUA':33.3, 'AUG':0.002,     
                        'ACU':33.3, 'ACC':33.3, 'ACA':33.3, 'ACG':33.3,     
                        'AAU':33.3, 'AAC':33.3, 'AAA':33.3, 'AAG':33.3,     
                        'AGU':33.3, 'AGC':33.3, 'AGA':33.3, 'AGG':33.3,     
                        'GUU':33.3, 'GUC':33.3, 'GUA':33.3, 'GUG':33.3,     
                        'GCU':33.3, 'GCC':33.3, 'GCA':33.3, 'GCG':33.3,     
                        'GAU':33.3, 'GAC':33.3, 'GAA':33.3, 'GAG':33.3,    
                        'GGU':33.3, 'GGC':33.3, 'GGA':33.3, 'GGG':33.3}
    
    
    reconstant_list=[]   # will be used to store the amino acid sequence translated from the mRNA sequence
    
    start_point=sequence.find('AUG')     # search the start codon in the mRNA sequence
    trimmed_sequence=sequence[start_point:]   # trim the mRNA sequence from the start codon
    
    codons_list = [trimmed_sequence[i:i+3] for i in range(0, len(trimmed_sequence), 3)]  # trim the sequence into codon list
    
    stop_list=['UAA', 'UAG','UGA']   # stop codon list
    
    for i in codons_list:
        if len(i)==1 or len(i)==2:   # means the length of last element is less than 3, not a codon, just pass
            pass
        elif i not in stop_list: # the codon is not stop codon, get the corresponding reaction constant and store into 
            reconstant_list.append(codonconstant_dict[i])  # reconstant_list
        else:  # the codon is stop codon, break
            reconstant_list.append(codonconstant_dict[i]) 
            break
    
    return reconstant_list
            

def directmethod(codon_number,aa_list, complex_list,reconstant_list, total_rib, rna, t1, occupied_pos): 
    '''Use the Direct Method to simulate one reaction step of the translation reaction in Arkin et al 
    
    parameters
    ----------
    codon_number: int
        the number of codons in the given mRNA sequence
    aa_list:list
        list that stores the amino acids translated from the given mRNA sequence
    complex_list: list
        list that stores the molecular number of the transitional Ribosome-RNA[i] complex
    reconstant_list: list
        list that stores reaction constants of all the reactions (values will not change as reaction occur)
    total_rib: int
        molecular number of the free ribosome
    rna: int
        molecular number of mRNA 
    t1: int
        simulation initial time
    occupied_pos: list
        list that store the indexes of codons that have been occupied by Ribosome-RNA[i] complex
        
    Returns
    -------
    out: tuple with 6 elements
       first element: list
            update the No. of Ribosome-RNA[i] complex after one step
       second element: int
            update the No. of free ribosome after one step
       third element: int
            update the No. of free mRNA-RBS after one step
       forth element: float
            update the time after one step
       fifth element: str
            the translated amino acid by this reaction step
       sixth element: int
            the selected codon position which the reaction will occur
    '''        
    
    n=codon_number  # get the number of codons in the mRNA sequence
    
    rib=total_rib  # get the initial number of free ribosome and RNA-RBS 
    rna_rbs=rna
    
    prop_list=[] 
    prop1=propensity(reconstant_list[0],rib,rna_rbs)  # calculate the propensity of ribosome binding into mRNA-RBS
    
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
                          # else means occupied_pos is [](empty), no position has been occupied, new_prop is equal to prop_list
    
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
    
    aa=aa_list[u]  # get the amino acid translated from this step of reaction
               
    return newcomplex_list, rib1, rna_rbs1, tau, aa, u
        
        

def simulation(sequence,total_rib,t1):
    '''Stimulate the translation reactions to translate 100 proteins
    
    parameters
    ----------
    sequence: str
        the given mRNA sequence
    total_rib: int
        initial value of the free ribosome
    t1: int
        simulation initial time
    
    Returns
    -------
    out: tuple with 4 elements
       first element:the dictionary that store the time and positions of each individual reaction
       second element:the dictionary that store the translated amino acid sequence of each individual reaction   
       third element: the list that stores the generated u in each step
       forth element: the time that needed to translate 100 proteins
    '''
    
    
    t0=t1  # use t0 to store the initial value of t1 (stimulation starting time)
    t1_list=[]  # will be used to store the updating reaction time t1
    
    s=translation(sequence)  # translate the mRNA sequence
    n=s[0]  # the number of codons in the mRNA sequence
    aa_list=s[1]  # the amino acids translated from the given mRNA sequence
    
    reconstant_list=calculation_constant(sequence) # get the reaction constant of each codon in the given mRNA sequence
    
    complex_list=[0]*(n-1)  # set the initial value of complex_list (all the elements are 0)
     
    rib=total_rib # get the initial number of free ribosome
    
    rna_rbs=1  # the initial number of mRNA-RBS is 1 as only 1 mRNA is involved in the simulation
    
    result={}  # will be used to store the time and positions of each individual reaction
    aa_result={}  # will be used to store the translated amino acid sequence of each individual reaction
    
    ids=0  # will be used as the ID of new reaction (also as the name of reaction)
    u_list=[]   # will be used to store the generated u at each step
    
    occupied_pos=[]  # will be used to store the codon postions that has been occupied (reaction cannot occur)

    number=0   # the number of proteins that have been translated completely
    
    while number<100:  

        y=directmethod(n,aa_list, complex_list,reconstant_list, rib,rna_rbs,t1, occupied_pos)  # run direct method
        u=y[-1]  # get u
        aa=y[-2]  # get the translated amino acid from this step of reaction
        
        if u_list==[]:    # at first u_list is empty, 
            complex_list=y[0]  # update the molecular number of Rib-RNA complex[i], free ribosome and RNA_RBS
            rib=y[1]
            rna_rbs=y[2]

            t1+=y[3]  # update the reaction time t1
            t1_list.append(t1)  # store the updated time 
            
            u_list.append(u)  # store the u generated from this step (then u_list is not empty)
            
            ids+=1  # start a new reaction (ids will be used as the name of new reaction)
            result[ids]=[(t1,u)]  # store the tuple of (reaction time, u) into dictionary (the ids will be the key)
            aa_result[ids]=[aa]  # store the generated amino acid under the key of ids (ids is the name of reaction)
            
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
                aa_result[ids]=[aa] # store the generated amino acid under the key of ids (ids is the name of reaction)
                
            else:  # if u!=0, means one of the current peptide reaction is continuing
                t1 += y[3]  # updating the reaction time
                t1_list.append(t1)  # store the updated reaction time t1 into t1_list
                for k, v in result.items(): # then check which reaction is continuing,
                    if v[-1][-1]==u-1: # v[-1][-1] is the last codon position of each individual reaction, if v[-1][-1]==u-1,
                        result[k].append((t1, u)) # it means this peptide translation reaction is occuring, just append the 
                        aa_result[k].append(aa) # tuple of (reaction time, u) into result dictionary under this key, and also
                        break                   # append the translated amino acid into aa_result dictionary under this key

            occupied_pos = [v[-1][-1] for v in result.values() if v[-1][-1]!=n-1]
            occupied_pos = [v - 1 for v in occupied_pos if v!=0] + occupied_pos   
                                                   # first get the occupied codon positions and store them into steps (if
                                                  # v[-1][-1]==n-1, it means the ribosome has encountered the last codon in the 
                                                 # given sequence, and the ribosome should be released from the mRNA, so do not 
                                                 # need to store it into steps)
            
            if u==n-1:  # if u==n-1, it means the one protein has been translated completely
                number+=1  # update the number of proteins that have been translated completely.
            else:
                pass
    
#     for k,v in aa_result.items():
#         print('Reaction {}: {}'.format(k,''.join(v)))        
    
    return result,aa_result,u_list, t1
       
    

def variedribosome(sequence,t1, rib_number, repeat_times):
    '''Calculate the time needed to translate 100 proteins with different number of ribosomes
    
    parameters
    ----------
    sequence: str
        the given mRNA sequence
    t1: int
        simulation initial time
    rib_number:int
        the maximum number of ribosomes we would like to run the simulation
    repeat_times: int
        how many times we would to repeat the simulation for different number of ribosomes
    
    Returns
    -------
    out: tuple with 4 elements
       first element: the list that stores the number of ribosomes used in the simulation
       second element:the list that stores the time needed to translate 100 proteins with different number of ribosomes
       third element: the list that stores the mean time needed to translate 100 proteins with different number of ribosomes
       forth element: the list that stores the standard deviation of the time needed to translate 100 proteins with different
                      number of ribosomes.
    '''
    
    meantime_list=[]
    std_list=[]
    rib_list=[]
    
    for i in range(1, rib_number+1):
        rib_list.append(i)
        t_list=[]
        for j in range(repeat_times):
            x=simulation(sequence, i, t1)
            t_list.append(x[3])
            
        meantime_list.append(sum(t_list)/repeat_times)
        std_list.append(st.stdev(t_list))
        print('Meantime',meantime_list, len(meantime_list))
        print('STD',std_list, len(std_list))
            
    return rib_list, t_list, meantime_list, std_list
