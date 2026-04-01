"""
Created on Dec 18, 2015
Python 3
       
"""


import numpy as np
from numpy import random
import matplotlib.pyplot as plt
%matplotlib inline



def propensity1(k, x, y):
    '''Calculate the propersity of a reaction with two reactants
    
    parameters
    ----------
    k:int
      reaction constant
    x:int
      molecular number of first reactant
    y:int
      molecular number of first reactant
      
    Returns
    -------
    out: int
       propensity of a reaction
    '''
    return k*x*y

def propensity2(k, x):
    '''Calculate the propersity of a reaction with one reactant
    
    parameters
    ----------
    k:int
      reaction constant
    x:int
      molecular number of reactant
      
    Returns
    -------
    out: int
       propensity of a reaction
    '''
    return k*x


def reaction1(x,y,z):   #reaction 1,2,5: X+Y->Z
    '''Calculate the molecular number change of reaction type 1:X+Y->Z 
       (two reactants and one product)
    
    parameters
    ----------
    x:int
      molecular number of first reactant
    y:int
      molecular number of second reactant
    z:int
      molecular number of first product
      
    Returns
    -------
    out: int
       molecular numbers of reactants and product after the reaction has occurred once
    '''
    if x>0 and y>0:
        x, y, z = x-1,y-1,z+1
    return x, y, z


def reaction2(a, d, e=0):   #reaction 3: a+b->c+d
    '''Calculate the molecular number change of reaction type 2:A+E->E+D
       (two reactants and two products, but one of the reactants is also the product)
    
    parameters
    ----------
    a:int
      molecular number of first reactant
    d:int
      molecular number of second reactant
    e:
      molecular number of the chemical which both appears in reactant and product
      
    Returns
    -------
    out: int
       molecular numbers of reactants and products after the reaction has occurred once 
       (value will not change after reaction)
    '''
    if a>0:
        a, d, e = a-1, e, d+1
    return a, e, d


def reaction3(e,f,g):   #reaction 4: E->F+G
    '''Calculate the molecular number change of reaction type 3:E->F+G
       (one reactant and two products)
    
    parameters
    ----------
    e:int
      molecular number of reactant
    f:int
      molecular number of first product
    g:int
      molecular number of second product
      
    Returns
    -------
    out: int
       molecular numbers of reactants and products after the reaction has occurred once
    '''
    if e>0:
        e, f, g = e-1,f+1,g+1
    return e,f,g
    
    
def firstreaction(a,b,c,d,e,f,g,k1=0,k2=0,k3=0,k4=0,k5=0):
    '''Use the First Reaction Method to simulate one reaction step of the system with 5 reactions in Gibson et al 
    
    parameters
    ----------
    a~g:int
      molecular number of all the 7 chemicals (A~G) involved in the 5 reactions
    k1~k5:int
      reaction constants of 5 chemical reactions (values will not change as reaction occur)
      
    Returns
    -------
    out: tuple with 3 elements
       first element: list
            update the No. of each species after one step
       second element: float
            update the time after one step
       third element: list
            update the propensity of each reaction after one step
    '''        
    prop1=propensity1(k1,a,b)        
    prop2=propensity1(k2,b,c)
    prop3=propensity1(k3,d,e)
    prop4=propensity2(k4,f)
    prop5=propensity1(k5,e,g)      # calculate the propensity of reaction 1~5
    
    prop_list=[prop1,prop2,prop3,prop4,prop5]  # store the propensity of reaction 1~5
    tau_list=[]                                # will be used to store the putative reaction time tau
    
    for i in prop_list:                        # use propensity to calculate the putative time of each reaction
        if i == 0:                             # propensity==0, reaction will never occur, store a infinite no. in tau_list
            tau_list.append(np.inf)
        else:                                  # propensity!=0, calculate the putative time and store it in tau_list
            tau_list.append(((1/i)*np.log(1/(1-np.random.random()))))
    
    if prop_list:                    # if prop_list!=[0,0,0,0,0], which  means system does not reach equilibrium,  choose which
        tau= min(tau_list)           # reaction will occur according to tau (the reaction with minimal tau will occur). 
        u=np.argmin(tau_list)
    
        if u==0:
            a,b,c = reaction1(a,b,c)
        elif u==1:
            b,c,d = reaction1(b,c,d)
        elif u==2:
            d,e,f = reaction2(d,f,e)
        elif u==3:
            f,d,g = reaction3(f,d,g)
        else:
            e,g,a = reaction1(e,g,a)   # update the no. of moleculars to reflect the execution of selected reaction
    
        species_list=[a,b,c,d,e,f,g]
    
        return species_list, tau, prop_list
    return                            # return NoneType if prop_list==[0,0,0,0,0]
   

def window(t1,t2,t1_list, n):
    '''Calculate the index of first element in reaction occurrance time list 
        between each certain time intervals from t1 to t2
    
    parameters
    ----------
    t1: int
      initial time
    t2: int
      ending time
    t1_list: list
      store the reaction occurance time during the stimulation
    n: int
      how many even time intervals we would like to get from t1 to t2
    
    Returns
    -------
    out: list
       the index of the first element in reaction occurrance time list between each 
       certain time intervals from t1 to t2
    '''

    t_list=np.linspace(t1,t2,n)    # divide the time period from t1 to t2 into n time intervals evenly
    
    index_x = []                    
    for j in range(len(t_list)-1):
        elements_between = []
        for i in range(len(t1_list)):
            if t1_list[i]>=t_list[j] and t1_list[i]<t_list[j+1]: # find the indexes of all elements in t1_list between two time 
                elements_between.append(i)                      # intervals in t_list
        if len(elements_between)>0:             # if there are elements in t1_list between two time intervals in t_list, then 
            index_x.append(elements_between[0]) # put the index of the first element in index_x (i.e. the first element in
                                                # list elements_between)
    return index_x


    
