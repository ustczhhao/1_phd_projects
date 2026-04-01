# NK Fitness Landscape Simulation

This project implements a series of Python simulations based on the **NK model** of fitness landscapes, a classic framework for studying how epistasis shapes evolutionary dynamics. The code explores how landscape structure, interaction topology, and fitness-combination rules influence adaptive walks, random walks, local optima, and pairwise epistasis.

Originally developed during my PhD training, this project reflects an early effort to translate theoretical evolutionary concepts into executable computational models. It combines landscape generation, evolutionary trajectory simulation, and quantitative analysis into a unified Python-based workflow.

---

## Project Overview

The **NK model** is a widely used theoretical framework for studying rugged fitness landscapes. In this model:

- **N** represents the number of loci in a genotype
- **K** represents the number of interacting or epistatic loci affecting each site

As K increases, the landscape generally becomes more rugged, with more local optima and more complex adaptive dynamics.

This repository contains code for:

- generating complete NK fitness landscapes
- constructing landscape subsets for simulation
- comparing different epistatic network structures
- simulating **random walks** and **adaptive walks**
- identifying **local optima**
- computing **pairwise epistasis**
- summarizing walk length and fitness gain
- plotting results under different model settings

---

## Scientific Motivation

Evolution on rugged landscapes depends not only on the fitness of individual mutations, but also on how mutations interact with one another. The NK model provides a simple but powerful way to investigate questions such as:

- How does increasing epistasis change the ruggedness of the landscape?
- How many local optima emerge under different interaction structures?
- How far can an evolving genotype move before reaching a local peak?
- How do adaptive walks differ from random uphill walks?
- How does the distribution of pairwise epistasis differ between additive and multiplicative fitness models?

This project was built to explore such questions computationally.

---

## What This Repository Demonstrates

This repository demonstrates:

- implementation of the **NK fitness landscape model**
- simulation of genotype-to-fitness mappings over binary sequence space
- comparison of **random** versus **adjacent-neighbor** epistatic interactions
- analysis of **additive** and **multiplicative** fitness formulations
- simulation of **random walks** and **adaptive walks**
- identification of **local fitness peaks**
- calculation of **pairwise epistasis**
- statistical summarization of mean walk length and fitness outcomes
- plotting and interpretation of evolutionary simulation results

---

## Core Modeling Components

### 1. Whole-Landscape Construction
The project includes scripts that generate the full binary genotype space and assign fitness values based on NK-style epistatic interactions.

Key features include:

- enumeration of all possible genotypes of length `N`
- construction of site-specific epistatic neighborhoods
- assignment of fitness contributions for all allele combinations
- calculation of genotype-level fitness across the full landscape

### 2. Epistatic Network Design
Different interaction structures are explored, including:

- **random epistatic sites**
- **adjacent/flanking epistatic neighbors**

This allows comparison between different assumptions about how loci interact.

### 3. Random Walks and Adaptive Walks
The code simulates evolutionary trajectories across the landscape, including:

- **random uphill walk**: moving to a randomly chosen neighbor with fitness at least as high as the current genotype
- **adaptive walk**: moving according to more selective rules toward fitter neighbors

These simulations are used to characterize accessibility, walk length, and endpoint fitness.

### 4. Local Optima Search
The project identifies genotypes whose fitness exceeds that of all single-mutant neighbors, allowing analysis of landscape ruggedness and peak structure.

### 5. Pairwise Epistasis Calculation
Several notebooks analyze pairwise epistatic effects across genotype subsets or landscapes, helping quantify interaction patterns between mutations under different fitness models.

### 6. Figure Generation and Summary Analysis
The repository also contains plotting notebooks used to summarize:

- mean number of steps
- final fitness reached
- differences between walk types
- effects of additive versus multiplicative models
- effects of random versus structured epistasis

---

## Main Research Questions Explored

This project was designed to explore questions such as:

- How does the NK parameter `K` influence adaptive accessibility?
- How does the number of local optima change with interaction structure?
- How do additive and multiplicative fitness models differ in their evolutionary consequences?
- Do adjacent-neighbor interactions produce different walk behavior from random interaction networks?
- How does pairwise epistasis distribute across different regions of the landscape?

---

## Repository Structure

```text
.
├── 0_NK Model_Whole Landscape.py
├── 0_NK Model_Whole Landscape2.py
├── 0_NK Model_Flash in Landscape.py
├── 0_NK Model_Calculate Epistasis.py
├── Exercise 7-3_...
├── Exercise 7-4 No. of local opitima_0617.ipynb
├── Exercise 7-5_Adaptive walks_20160714.ipynb
├── Exercise 7-6_...
├── Exercise 7-9_...
├── Exercise 7_11_Calculate pairwise e_0815_Additive_OK.ipynb
├── Exercise 7_15_Mean Steps and fitness_0819-Use this one for ploting figure.ipynb
├── Plot Fig of NK_model-0924-...
└── Test_0921_Calculate all the e in a landscape subset-...