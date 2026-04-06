# Project: Effects of Genome Architecture on Tetrahymena Evolution

## Overview

This project explores how key features of **genome architecture** influence evolutionary dynamics in *Tetrahymena*. Building on the broader theoretical framework that amitosis can confer some of the benefits typically associated with sex, this repository investigates three related genomic factors:

1. **chromosome copy number control during amitosis**;
2. **multiple linked loci located on the same MIC chromosome**;
3. **cross-over events during sex**.

Together, these analyses examine how noncanonical genome structure and inheritance mechanisms shape mutation load, fitness dynamics, and the long-term evolutionary behavior of *Tetrahymena* populations.

The repository consists primarily of simulation and figure-generation notebooks that extend the core amitosis framework to more realistic or mechanistically enriched genomic settings.

## Scientific Motivation

*Tetrahymena* has an unusual nuclear architecture in which germline and somatic functions are separated into two nuclei within a single cell: the **micronucleus (MIC)** and the **macronucleus (MAC)**. During asexual growth, the MIC divides by mitosis, whereas the MAC divides by **amitosis**. A distinctive feature of *Tetrahymena* is that amitotic MAC division is accompanied by an unknown mechanism of **chromosome copy number control**, which maintains approximately stable ploidy across generations. This unusual architecture is thought to be central to the long-term persistence of asexual lineages in the genus. 

The analyses in this repository extend that idea by asking how additional genomic features influence evolutionary outcomes. In particular:

- How robust are the benefits of amitosis when chromosome copy number control is imperfect or unstable?
- How does the physical organization of multiple linked loci on the same MIC chromosome affect selection and long-term fitness?
- To what extent can cross-over during sex alter the relative advantage of amitosis, mitosis, or sex?

These questions connect genome structure directly to evolutionary performance.

## Project Goals

The main objectives of this repository are to:

- evaluate the evolutionary consequences of **unstable or imperfect chromosome copy number control** during amitosis;
- examine how **linkage among multiple loci within a MIC chromosome** influences fitness dynamics;
- assess the effect of introducing **cross-over** into linked-locus models;
- compare evolutionary outcomes under alternative reproductive or genomic scenarios;
- generate figure-ready outputs for interpreting how genome architecture shapes adaptation and mutation load.

## Scope of the Repository

This repository is organized around three related themes.

### 1. Chromosome Copy Number Control During Amitosis

Amitosis in *Tetrahymena* is biologically meaningful because it is associated with **roughly constant chromosome copy number across daughter cells** rather than unrestricted ploidy fluctuation. This project includes simulations that relax or perturb that assumption by introducing **unstable ploidy / unstable copy number control**. These analyses explore whether the evolutionary benefits of amitosis persist when copy number regulation is imperfect, and how sensitive long-term fitness is to such instability.

### 2. Multiple Linked Loci on the Same MIC Chromosome

The project also examines models in which **multiple fitness-affecting loci reside on the same MIC chromosome**. This introduces linkage structure into the framework and allows the analysis to move beyond single-locus or fully independent-locus approximations. These notebooks investigate how the number of loci per chromosome affects fitness outcomes and how chromosomal organization itself may alter the evolutionary consequences of amitosis and sex.

### 3. Cross-over as a Genomic Event

A further extension introduces **cross-over** into the linked-locus MIC chromosome models. By allowing recombination among loci, these analyses explore how the shuffling of linked alleles changes the relative dynamics of mutation accumulation, fitness recovery, and adaptation. This makes it possible to compare amitosis-driven variance generation with recombination-mediated reshuffling in a shared modeling framework.

## Repository Contents

### A. Unstable copy number / unstable ploidy during amitosis

- `File 8-11-1 Introduce Unstable Copy Number During Amitosis_20170312.ipynb`  
  Simulation notebook introducing instability in chromosome copy number during amitotic division.

- `File 8-11-2 Introduce Unstable Ploidy Copy Number during Amitosis_WF_20170204.ipynb`  
  Wright-Fisher implementation of unstable ploidy / copy number dynamics during amitosis.

- `File 8-11-3 Unstable Ploidy_Check Equilibrium_20170528.ipynb`  
  Equilibrium-focused analysis of unstable ploidy scenarios.

- `File 8-11-4 Plot Figures for Fitness Comparison of Unstable Ploidy During Amitosis_20170228.ipynb`  
  Figure-generation notebook comparing fitness outcomes under unstable ploidy.

- `File 8-11-14 Plot Figures for Unstable Ploidy During Amitosis_20170428.ipynb`  
  Additional plotting workflow for unstable ploidy analyses.

- `File 8-11-15 Plot Figures for Unstable Ploidy During Amitosis_20170428_Add PUC.ipynb`  
  Extended plotting notebook for unstable ploidy / copy-number-control comparisons.

### B. Multiple loci within a MIC chromosome

- `Code_Multiple loci within MIC Chromosome_No Crossover.ipynb`  
  Core simulation notebook for linked loci on a MIC chromosome without recombination.

- `Fitness comparison_Different # of loci in Germ Chrom_N =20_RM Every 100 Gens.ipynb`  
  Fitness comparison across different numbers of loci per germline chromosome under periodic random mating.

- `Fitness comparison_Different # of loci in Germ Chrom_N =20_RM Every Gen.ipynb`  
  Fitness comparison across different locus counts under more frequent random mating.

### C. Cross-over introduced into linked-locus MIC chromosome models

- `Code_Multiple loci within MIC Chromosome_Introduce Crossover.ipynb`  
  Core linked-locus simulation notebook with recombination introduced.

- `RM E100 Introduce CrossOver_L =20 loci and 1 CrossOver Per Chromosome.ipynb`  
  Simulation of cross-over in a 20-locus chromosome with periodic random mating.

- `RM E100 Introduce CrossOver_L =100 loci Per Chromosome_Different CrossOver level.ipynb`  
  Cross-over analysis for 100-locus chromosomes across multiple recombination levels.

- `RM Introduce CrossOver_L =100 loci Per Chromosome_Different CrossOver level.ipynb`  
  Additional simulation notebook for varying cross-over levels in high-locus linked models.

- `RM_Introduce CrossOver_L =20 loci and 1 CrossOver Per Chromosome.ipynb`  
  Related cross-over simulation workflow for the 20-locus setting.

## Methods and Approach

### Simulation-based evolutionary modeling

This project uses **population-genetic simulation notebooks** to study how genome architecture affects evolutionary outcomes under mutation, selection, and alternative reproductive mechanisms. Rather than treating all loci as fully independent, the notebooks progressively introduce biologically motivated structural constraints and genomic events.

### Comparative framework

Across the repository, the analyses compare how evolutionary outcomes change when one modifies:

- the stability of chromosome copy number during amitosis;
- the number of linked loci carried on a MIC chromosome;
- the frequency or presence of cross-over;
- the schedule of random mating or sexual reproduction in some scenarios.

This comparative design allows the effect of each genomic feature to be studied in relation to the baseline amitosis framework.

### Fitness-based evaluation

Many notebooks are organized around **fitness comparison**, including trajectory-based or equilibrium-based contrasts across alternative settings. This makes it possible to evaluate how specific genomic mechanisms influence mutation load, variance generation, and long-term mean fitness.

## Scientific Interpretation

A central idea motivating this repository is that the benefits of amitosis depend not only on the existence of amitotic segregation itself, but also on the **genomic architecture in which it operates**. In *Tetrahymena*, chromosome copy number control appears to be a crucial modifier of amitotic dynamics, because uncontrolled chromosome-number drift would otherwise lead to severe instability. Likewise, the arrangement of multiple loci on the same MIC chromosome and the possibility of cross-over during sexual processes may substantially change how variation is generated and how efficiently selection can act.

This project therefore moves beyond a minimal amitosis-versus-mitosis comparison and asks a broader question:

**How do specific genomic design features alter the evolutionary consequences of nonstandard inheritance?**

## Technical Highlights

This repository demonstrates several features of quantitative evolutionary research:

- simulation-based modeling of genome-architecture effects on evolution;
- extension of baseline inheritance models to structurally richer chromosomal settings;
- analysis of linked-locus systems rather than fully independent loci alone;
- explicit treatment of chromosome copy-number stability as an evolutionary parameter;
- incorporation of recombination-like genomic events into comparative evolutionary models;
- figure-generation workflows for interpreting simulation results in a research setting.

## Broader Relevance

Although the biological system is *Tetrahymena*, the questions addressed here are relevant more broadly to evolutionary genetics and quantitative biology. The repository illustrates how chromosome organization, ploidy regulation, and recombination structure can affect long-term evolutionary trajectories. It also demonstrates a general computational workflow for testing mechanistic hypotheses about genome structure using simulation-based approaches.

More broadly, the project reflects a style of research that links:

- biological mechanism,
- formal evolutionary reasoning,
- computational implementation,
- and interpretable quantitative outputs.



## Note

This repository is presented as a research project archive documenting simulation-based analyses of how genome architecture influences evolution in *Tetrahymena*. Some notebooks reflect iterative research-stage development and may require minor cleanup or environment-specific adjustment for re-execution in newer software versions. The data used to generate the figures can be reproduced by running the corresponding simulation and plotting notebooks included in this repository.