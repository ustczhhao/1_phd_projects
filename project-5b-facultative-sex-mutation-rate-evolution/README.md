# Project 5b: Facultative Sex and Mutation-Rate Evolution in Single-Nucleus Systems

## Overview

This project extends the mutation-rate evolution framework developed for *Tetrahymena* to a simpler class of organisms with a **single nucleus**. The goal is to isolate the effect of **facultative sex** on the evolution of mutation rate without the added complexity of nuclear dimorphism, amitotic division, or macronuclear copy-number control.

Using stochastic population-genetic simulations, this repository examines how mutation-rate modifiers evolve in **haploid** and **diploid** single-nucleus systems under **facultative sexual reproduction**, with comparisons across several important dimensions:

- deleterious-only versus mixed beneficial/deleterious mutation regimes;
- haploid versus diploid genome structure;
- additive, dominant, and recessive modifier effects;
- different frequencies of facultative sex;
- and different recombination assumptions between fitness loci and mutation-rate modifier loci. 

More broadly, this project asks a general evolutionary question:

**How does facultative sex affect the long-term evolution of mutation rate in single-nucleus organisms?**

## Scientific Motivation

Mutation rate is itself an evolvable trait. Although mutator and antimutator alleles do not usually affect fitness directly, they can change the rate at which beneficial and deleterious mutations arise across the genome. Their evolutionary fate is therefore shaped indirectly through association with the fitness consequences of the mutations they generate.

The related *Tetrahymena* project showed that life-cycle structure and genome architecture can strongly influence mutation-rate evolution. This repository removes the complications of dual nuclei and polyploid somatic inheritance in order to study a cleaner baseline system: a **single haploid or diploid nucleus** reproducing mostly asexually, but undergoing **facultative sex** at specified intervals. This makes it possible to distinguish which effects arise from facultative sex itself and which depend specifically on ciliate-like genome architecture. 

The project is especially motivated by three questions:

1. Can facultative sex reduce mutation rate in a standard single-nucleus system?
2. How do ploidy and dominance alter selection on mutator and antimutator modifiers?
3. How does recombination between fitness loci and modifier loci affect the strength of indirect selection on mutation rate? 

## Project Goals

The main objectives of this repository are to:

- model mutation-rate evolution in **single diploid** and **single haploid** nuclear systems;
- compare outcomes under **deleterious-only** and **mixed beneficial/deleterious** mutation regimes;
- test the effect of **dominant**, **recessive**, and **additive** mutation-rate modifiers;
- evaluate how different **facultative sex frequencies** affect the evolution of mutation rate;
- investigate how **recombination** between fitness loci and mutation-rate modifier loci changes evolutionary dynamics;
- and compare the simulation outputs with theoretical expectations, including tests related to **Raynes’ equation**. 

## Core Findings

This repository is designed around the idea that the effect of facultative sex on mutation-rate evolution depends on both **genetic architecture** and **modifier inheritance**.

### 1. Facultative sex can alter the evolution of mutation rate even in a single-nucleus system

By periodically reshuffling genotypes through sex, facultative sex changes the covariance structure between organismal fitness and mutation rate. The simulations explicitly track population mean fitness, mutation-rate means and variances, and the covariance between fitness and mutation rate, allowing the evolutionary response of mutator and antimutator alleles to be studied directly. 

### 2. Diploid and haploid systems behave differently

The repository includes dedicated implementations for both **single diploid** and **single haploid** nuclei, allowing direct comparison of how ploidy changes selection on mutation-rate modifiers. In diploids, the effect of modifier dominance can be studied explicitly, whereas haploid simulations provide a useful baseline without dominance masking. 

### 3. Modifier dominance changes mutation-rate evolution

Several scripts compare **additive**, **dominant**, and **recessive** mutator/antimutator models. These simulations make it possible to determine how the expression pattern of modifier loci affects the efficiency with which selection can reduce mutation rate. The dominant and recessive versions implement distinct genotype-to-mutation-rate mappings in the diploid state. 

### 4. Recombination structure matters

The project includes simulations in which recombination is restricted to different parts of the genome: one version allows recombination among **fitness loci**, while another allows recombination among **mutation-rate modifier loci**. These variants test how breaking or preserving linkage influences the buildup of covariance between fitness and mutation rate. 

## Relationship to Project 5

This repository is a conceptual extension of the preceding project on mutation-rate evolution in *Tetrahymena*. Project 5 focused on a dual-genome ciliate system with a diploid MIC and polyploid MAC, where facultative sex interacts with nuclear dimorphism and amitotic somatic reproduction.

By contrast, this project removes those ciliate-specific features and asks what remains in a more conventional **single-nucleus facultative-sex model**. It therefore serves as a useful comparative baseline for understanding which results are specific to *Tetrahymena* genome architecture and which reflect more general principles of mutation-rate evolution under facultative sex.

## Repository Contents

This folder contains both simulation scripts and summary notebooks related to mutation-rate evolution in single-nucleus systems.

### A. Deleterious-only mutation-rate evolution

- `File 5-1. Code_Evolution of Mu Rate_Single Diploid Nucleus_AllDeleMut.py`  
  Core diploid single-nucleus simulation under deleterious-only mutation with additive-style modifier effects. The script tracks population mean fitness, mutation rate, mutation-rate variance, and covariance between fitness and mutation rate. :contentReference[oaicite:8]{index=8}

- `File 5-1a. Code_Evolution of Mu Rate_Single Diploid Nucleus_AllDeleMut_DominantMut.py`  
  Diploid deleterious-only simulation with **dominant** modifier effects, where a nonzero modifier state is sufficient to alter mutation rate. :contentReference[oaicite:9]{index=9}

- `File 5-1b. Code_Evolution of Mu Rate_Single Diploid Nucleus_AllDeleMut_RecessiveMut.py`  
  Diploid deleterious-only simulation with **recessive** modifier effects, where modifier impact appears only in the homozygous state. :contentReference[oaicite:10]{index=10}

- `File 5-2. Code_Evolution of Mu Rate_Single Haploid Nucleus_AllDeleMut.py`  
  Haploid single-nucleus simulation under deleterious-only mutation, used to compare ploidy effects. The sexual cycle is modeled through temporary diploid zygote formation followed by restoration of the haploid state. :contentReference[oaicite:11]{index=11}

### B. Mixed beneficial and deleterious mutation regimes

- `File 5-3. Code_Evolution of Mu Rate_Single Diploid Nucleus_1PerBeneMut.py`  
  Diploid simulation with both beneficial and deleterious mutations, including explicit tracking of separate beneficial and deleterious mutation rates and their joint effect on fitness. :contentReference[oaicite:12]{index=12}

- `File 5-3a. Code_Evolution of Mu Rate_Single Diploid Nucleus_1PerBeneMut_DominantMut.py`  
  Mixed-mutation diploid simulation with **dominant** modifier effects. :contentReference[oaicite:13]{index=13}

- `File 5-3b. Code_Evolution of Mu Rate_Single Diploid Nucleus_1PerBeneMut_RecessiveMut.py`  
  Mixed-mutation diploid simulation with **recessive** modifier effects. :contentReference[oaicite:14]{index=14}

### C. Recombination-structure extensions

- `File 5-4a. Code_Evolution of Mu Rate_1PerBene_Part of Genome Recombine_FitLociRecom.py`  
  Mixed-mutation simulation in which recombination is introduced primarily among **fitness loci**, while modifier loci remain inherited with more limited reshuffling. This isolates the effect of recombination on the selected portion of the genome. :contentReference[oaicite:15]{index=15}

- `File 5-4b. Code_Evolution of Mu Rate_1PerBene_Part of Genome Recombine_MutModLociRecom.py`  
  Mixed-mutation simulation in which recombination is introduced for **mutation-rate modifier loci**, allowing direct comparison with the fitness-loci recombination scenario. :contentReference[oaicite:16]{index=16}

### D. Analysis and figure-generation notebooks

- `File 5-5. SN_1% Bene_N2K_Diploid_Facultative Sex_Mod Additive.ipynb`  
  Summary notebook for diploid facultative-sex simulations with 1% beneficial mutations and additive-style modifiers.

- `File 5-6. 0% Bene_Diploid SN_Comparing with Dele Only Model.ipynb`  
  Comparison notebook for deleterious-only settings and related baseline models.

- `File 5-7. SN_Bene_N2K_Diploid_Diff Bene Mut Prop.ipynb`  
  Analysis of how varying the proportion of beneficial mutations changes mutation-rate evolution.

- `File 5-8. SN_1% Bene_N2K_Diploid_Dominant and Recessive Mod.ipynb`  
  Comparison notebook for dominant versus recessive modifier behavior.

- `File 5-9. Test Raynes' Equation Using Our Parameters.ipynb`  
  Notebook testing the relationship between simulation outputs and theoretical predictions inspired by Raynes’ equation.

- `File 5-10. SN_ 1% Bene_Mito_Diff N_Test Raynes' Equation_Diff N_20KG.ipynb`  
  Parameter-sweep notebook evaluating theoretical expectations across different population sizes.

## Methods and Approach

### Stochastic population-genetic simulation

The simulations are implemented using Wright–Fisher-style population-genetic models in which each generation includes mutation, fitness-dependent parent sampling, and either asexual reproduction or sex at fixed intervals. Population state is stored explicitly across replicate populations, allowing the dynamics of fitness and mutation rate to be followed over time. 

### Mutation-rate modifier loci

Each simulation includes two sets of modifier loci:

- **upper controllers**, which increase mutation rate;
- **lower controllers**, which decrease mutation rate.

These loci do not directly determine organismal fitness. Instead, they modify the per-locus genomic mutation rate, and are therefore favored or disfavored only through the fitness effects of the mutations they help generate. The scripts track the mean mutation rate, variance in mutation rate, and the covariance between mutation rate and fitness over time. 

### Fitness loci and mutation regimes

The project includes both:

- **deleterious-only models**, where mutation pressure opposes selection directly; and
- **mixed beneficial/deleterious models**, where rare beneficial mutations introduce a countervailing force favoring higher mutation rates in some contexts. :contentReference[oaicite:19]{index=19}

In the mixed-mutation simulations, beneficial and deleterious mutations are stored separately and contribute multiplicatively to organismal fitness. :contentReference[oaicite:20]{index=20}

### Ploidy and dominance structure

The diploid scripts allow modifier effects to be modeled as:

- additive or dosage-like,
- dominant,
- or recessive.

This provides a direct way to test how genotype expression at modifier loci influences mutation-rate evolution. The haploid script removes this layer of complexity and serves as a simpler comparison model. 

### Facultative sex and recombination

Sex occurs at user-defined intervals through gamete formation, zygote construction, and restoration of the next-generation nuclear state. Additional scripts selectively alter recombination structure so that the contributions of fitness-locus recombination and modifier-locus recombination can be separated. 

## Scientific Interpretation

The key contribution of this project is to show how **facultative sex**, even in a standard single-nucleus organism, can reshape the evolution of mutation rate by modifying linkage, covariance, and the balance between short-term and long-term selection pressures.

This repository is especially useful for distinguishing several effects that are often confounded:

- the role of **ploidy** versus the role of **nuclear architecture**;
- the effect of **sex frequency** versus the effect of **modifier dominance**;
- and the contribution of **recombination pattern** versus the contribution of mutation regime.

By stripping away the unusual biology of ciliates, the project helps identify the more general evolutionary principles underlying mutation-rate evolution under facultative sex.

## Technical Highlights

This repository demonstrates several strengths characteristic of computational evolutionary genetics:

- implementation of stochastic mutation-rate evolution models in both haploid and diploid systems;
- explicit tracking of covariance between fitness and mutation rate;
- comparison of additive, dominant, and recessive modifier architectures;
- simulation of deleterious-only and mixed beneficial/deleterious mutation regimes;
- targeted manipulation of recombination structure to isolate linkage effects;
- and figure-oriented exploratory analysis linking simulations to theoretical expectations. 

## Broader Relevance

Although motivated by mutation-rate evolution, the framework used here is relevant to a wider class of questions in quantitative and evolutionary genetics, including:

- indirect selection on modifier alleles;
- the role of linkage in shaping long-term evolutionary outcomes;
- how sex changes covariance structure in adapting populations;
- and how genome organization influences the evolvability of fundamental parameters such as mutation rate.

More broadly, this repository illustrates a research workflow that combines mathematical thinking, simulation-based hypothesis testing, parameter sweeps, and theory-linked data analysis.

## Note

This repository is presented as a research project archive documenting simulation-based analyses of mutation-rate evolution under facultative sex in single-nucleus systems. Some scripts and notebooks reflect iterative research-stage development and may require minor cleanup or environment-specific adjustment for re-execution in newer software versions. The data used to generate the figures can be reproduced by running the corresponding simulation and plotting files included in this repository.