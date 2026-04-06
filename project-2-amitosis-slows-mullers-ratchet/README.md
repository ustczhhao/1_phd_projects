# Project 2: Amitosis Slows Muller’s Ratchet

## Overview

This project examines the **stochastic evolutionary benefit of amitosis** by asking whether amitotic division can reduce the accumulation of drift-induced mutation load relative to mitosis in finite populations. The analysis is based on the stochastic component of the accompanying *Tetrahymena* study and focuses on how amitosis alters long-term fitness trajectories under mutation, selection, and genetic drift.

In contrast to the deterministic equilibrium analysis in Project 1, this repository centers on **population-level stochastic simulations**. It documents computational experiments showing that amitosis can substantially slow, and in some settings effectively halt, the accumulation of deleterious mutations caused by drift. The project also includes comparisons with sexual reproduction under both deleterious-only and mixed beneficial/deleterious mutation regimes.

More broadly, this repository illustrates how a theoretical question in evolutionary genetics can be investigated through simulation-based modeling using multiple population-genetic frameworks.

## Scientific Motivation

The deterministic mutation-selection balance captures only part of the evolutionary picture. In finite populations, **genetic drift** can cause deleterious mutations to accumulate irreversibly even when selection acts against them. This process increases mutation load over time and lowers mean population fitness.

In asexual populations, this stochastic accumulation of deleterious mutations is known as **Muller’s ratchet**. Informally, Muller’s ratchet refers to the repeated stochastic loss of the least-mutated class of individuals in a finite asexual population. Once that class is lost, it cannot be perfectly regenerated in the absence of recombination, so the population “clicks” toward progressively higher deleterious load and lower fitness.

This project addresses the following question:

**Can amitosis slow Muller’s ratchet relative to mitosis by generating additional heritable variation on which selection can act?**

The broader significance of this question lies in understanding how unusual modes of inheritance may mitigate one of the classic long-term disadvantages of asexual reproduction.

## Project Goals

The main objectives of this repository are to:

- simulate mutation accumulation under finite population size using stochastic population-genetic models;
- compare fitness trajectories under mitosis, amitosis, and sex;
- quantify the extent to which amitosis slows the accumulation of drift load relative to mitosis;
- examine how the benefit depends on population size, ploidy, and mutation regime;
- reproduce and extend figure-generating workflows associated with the stochastic results in the accompanying study.

## Core Result

The central result of this project is that **amitosis makes populations less susceptible to Muller’s ratchet than mitosis**. In the accompanying study, finite mitotic populations accumulate drift load under deleterious mutation, whereas amitosis slows this process substantially and can even prevent sustained decline in some parameter settings. The benefit becomes stronger at higher ploidy and can be comparable to, or greater than, the benefit provided by facultative sex under biologically realistic reproductive schedules. :contentReference[oaicite:0]{index=0} :contentReference[oaicite:1]{index=1}



## Relationship to the Accompanying Manuscript

This repository corresponds to the **stochastic benefit** component of the accompanying manuscript `amitosis.pdf`. In that work, the deterministic and stochastic advantages of amitosis are analyzed separately. Project 1 focuses on the deterministic reduction in mutation load at equilibrium, whereas this project focuses on the stochastic slowing of drift-induced fitness decline in finite populations.

The manuscript explicitly describes this transition: after analyzing mutation-selection balance in large populations, it turns to the role of drift and evaluates the extent to which amitosis with chromosome copy-number control can slow the accumulation of drift load, i.e., Muller’s ratchet. :contentReference[oaicite:3]{index=3}

## Repository Contents

This folder contains simulation notebooks and figure-generation notebooks related to stochastic comparisons among mitosis, amitosis, and sex.

### Population-genetic simulation notebooks

- `Code_Moran Model_Both Bene and Dele_Amito and Sex.ipynb`  
  Moran-model simulations including both beneficial and deleterious mutations under amitosis and sex.

- `Code_Moran Model_Both Bene and Dele_Mitosis.ipynb`  
  Moran-model baseline simulations for mitosis under mixed mutation regimes.

- `Code_Moran Model_Dele Only_Amito and Sex.ipynb`  
  Moran-model simulations under deleterious-only mutation for amitosis and sex.

- `Code_Moran Model_Dele Only_Mito.ipynb`  
  Moran-model baseline simulations under deleterious-only mutation for mitosis.

- `Code_WF Model_Both Bene and Dele_Should be OK_20170920.ipynb`  
  Wright-Fisher simulations under mixed beneficial and deleterious mutation.

- `Code_WF Model_DeleOnly_Amito and Sex.ipynb`  
  Wright-Fisher simulations for deleterious-only scenarios under amitosis and sex.

- `Code_WF Model_DeleOnly_Amito and Sex_Revised_shouldbeOK.ipynb`  
  Revised Wright-Fisher implementation for deleterious-only simulations.

### Figure-generation and analysis notebooks

- `File 8-10-3 Plot Figure for Fit Comparison of Mitosis, Amitosis and Sex-20170107-MUL_10K_Stride =1.ipynb`  
  Plotting workflow for stochastic fitness comparisons among reproductive modes.

- `File 8-10-4 Test of Ver2 Code_Fitness Comparison_1224.ipynb`  
  Validation and comparison notebook for revised simulation outputs.

- `File 8-10-5_Plot Fig for Comparison of Variance of Amitosis, Mitosis, Selfing and RM_01022017.ipynb`  
  Analysis of variance generation under different reproductive strategies.

- `File 8-10-6_Plot Fig for Fitness Comparison_Amitosis for 200K_nInds =5, 10 and 20.ipynb`  
  Fitness comparison under amitosis across small population sizes.

- `File 8-10-7_Plot Fig for Fitness Comparison_Amitosis Vs Mitosis for 100K_For 20170505 EE Symposium.ipynb`  
  Amitosis-versus-mitosis comparison workflow.

- `File 8-10-8_Plot Fig for Fitness Comparison_Amitosis Vs Sex for 10K_N =20_For 20170505 EE Symposium.ipynb`  
  Amitosis-versus-sex comparison workflow.

- `File 8-10-10 Plot Figure for Both Beneficial and Deleterious Mutations (Same P and S)_20170920_N =20, 50 and 100_S =0.1_Fitness_Stride =1.ipynb`  
  Stochastic simulations and plotting under mixed mutation regimes across multiple population sizes.

- `File 8-10-11 Plot Figure for Both Beneficial and Deleterious Mutations (Same S_1% Beneficial)_20170929_N =20, 50 and 100_S =0.1_Fitness_First 1KG.ipynb`  
  Mixed-mutation analysis with rare beneficial mutations.

- `File 8-10-12 Plot Figure for Both Beneficial and Deleterious Mutations (Same S_1% Beneficial)_20171002_N =200, 500 and 1000_S =0.1_Fitness_1KG.ipynb`  
  Large-population extension of mixed-mutation simulations.

- `File 8-10-13 Fig 3 and Supp Fig 1. Stochastic benefit_P =2 and 45.ipynb`  
  Figure-generation notebook corresponding to the stochastic benefit analysis across ploidy, including manuscript Figure 3 and Supplementary Figure 1.

- `File 8-10-14 Figure 1.3. Fitness_Amito vs Mito_Different P_Amitosis_Diff Ploidy.ipynb`  
  Analysis of amitosis-versus-mitosis fitness dynamics across ploidy levels.

- `File 8-10-15 Figure 1.5 and Supp Figure 1A. Fitness comparison_Amito vs Sex_N =20_Dele only.ipynb`  
  Comparison of amitosis and sex in deleterious-only settings.

- `File 8-10-16 Figure 1.7B. 1% Beneficial Mutations_N =1K.ipynb`  
  Large-population mixed-mutation simulation with rare beneficial mutations.

- `File 8-10-17 Stoch benefit of amito over mito_N =10_Strategy 1_log fit.ipynb`  
  Focused analysis of stochastic advantage at small population size.

- `File 8-10-18 Stoch benefit of amito over mito_N =10_Strategy 1_log fit-Copy1.ipynb`  
  Alternate or backup version of the preceding analysis notebook.



## Methods and Approach

### Stochastic population-genetic modeling

This project studies mutation accumulation in **finite populations**, where random sampling effects can overpower selection and lead to progressive fitness decline. Unlike deterministic equilibrium analysis, the emphasis here is on evolutionary trajectories over time under repeated stochastic reproduction and mutation.

Two classic population-genetic frameworks appear in this repository:

- **Moran models**, which update populations in overlapping generations through sequential birth-death events;
- **Wright-Fisher models**, which update populations in discrete generations through multinomial sampling.

Using both frameworks provides a broader computational basis for evaluating whether the observed advantage of amitosis is robust to simulation architecture.

### Mutation regimes

The notebooks cover at least two major evolutionary settings:

- **deleterious-only mutation**, used to examine the operation of Muller’s ratchet directly;
- **mixed beneficial and deleterious mutation**, used to compare long-term adaptive responses under different reproductive strategies.

This allows the project to move beyond pure mutation accumulation and examine how amitosis performs in more general evolutionary scenarios.

### Comparative reproductive strategies

The simulations compare multiple reproductive modes, including:

- **mitosis**, as the baseline asexual mechanism;
- **amitosis**, with chromosome copy-number control;
- **sex**, including comparisons to selfing, random mating, and facultative versus more frequent sexual reproduction depending on the notebook.

These comparisons are scientifically important because the manuscript argues that some benefits of amitosis are analogous to classic mutational benefits of sex. :contentReference[oaicite:4]{index=4}

### Fitness trajectory analysis

Many notebooks generate time-series comparisons of mean population fitness, variance in fitness, and relative performance under different population sizes and ploidy levels. Several figure notebooks explicitly correspond to manuscript Figure 3, Figure 4, and supplementary stochastic analyses.

## Technical Highlights

This project demonstrates several strengths typical of quantitative evolutionary modeling:

- stochastic simulation of mutation accumulation in finite populations;
- implementation of both Moran and Wright-Fisher population-genetic models;
- comparative analysis of alternative reproductive strategies under matched mutation regimes;
- investigation of ploidy, population size, and mutation-spectrum effects on long-term fitness dynamics;
- integration of simulation output with figure-generation workflows suitable for scientific communication;
- organization of computational experiments around clearly defined evolutionary hypotheses.

## Biological Context

The biological setting is the evolution of amitosis in *Tetrahymena*, whose somatic nucleus divides amitotically while maintaining approximate chromosome copy-number control. The accompanying study proposes that this unusual mode of division generates more variation in fitness among offspring than mitosis, thereby increasing the efficiency of selection. In finite populations, that additional variance makes populations less vulnerable to the stochastic deterioration captured by Muller’s ratchet. :contentReference[oaicite:5]{index=5}

This mechanism is especially relevant in highly polyploid systems, where the stochastic advantage of amitosis over mitosis becomes more pronounced. :contentReference[oaicite:6]{index=6}

## Broader Relevance

Although this project is rooted in evolutionary theory, the workflow is broadly relevant to quantitative biology and computational research more generally. It demonstrates how to:

- formalize a biological hypothesis as a simulation problem;
- evaluate stochastic rather than purely deterministic system behavior;
- compare competing mechanisms under matched assumptions;
- interpret long-run population trajectories under uncertainty;
- connect raw simulation outputs to figure-ready scientific summaries.

The analytical habits reflected in this project are transferable to a wide range of problems involving evolutionary dynamics, probabilistic modeling, simulation-based inference, and computational biology.


## Note

This repository is presented as a research project archive documenting the simulation and analysis workflow for a theoretical evolutionary biology study. Some notebooks reflect iterative research-stage development and may require minor cleanup or environment-specific adjustment for re-execution in newer software versions. The data used to generate the figures can be reproduced by running the corresponding simulation and plotting notebooks included in this repository.