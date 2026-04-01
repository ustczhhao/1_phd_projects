# Simulation of Protein Translation Dynamics

This project implements a stochastic simulation framework for modeling protein translation dynamics along a single mRNA molecule. It combines codon-usage analysis, codon-dependent translation-rate modeling, and event-driven simulation of ribosome movement to study how sequence composition and codon arrangement influence translation efficiency and ribosome stalling.

Originally developed during my PhD training, this project reflects an early effort to connect biological sequence features with quantitative simulation, using Python to translate mechanistic assumptions into executable computational models.

---

## Project Overview

Protein translation is an inherently dynamic and stochastic process. The rate at which ribosomes move along an mRNA can vary depending on codon identity, codon usage bias, and local arrangement of slow-translating codons. This project was built to explore how such sequence-level variation affects overall translation behavior.

The repository contains several related components:

- **Codon usage analysis** of *E. coli* K12 MG1655
- **Direct-method stochastic simulation** of ribosome initiation, elongation, and termination
- **Sequence-dependent translation modeling** for specific mRNA sequences
- **Comparative simulation of different slow-codon arrangements**
- **Quantification of translation time and ribosome stalling**
- **Documented Jupyter notebooks** for exploratory analysis and simulation experiments

---

## Biological Motivation

Codon composition can shape translation kinetics by affecting ribosome progression along mRNA. In particular, clusters or distributions of slow-translating codons may alter elongation speed, induce ribosome queuing, and change the time required to complete protein synthesis.

This project explores these questions by simulating translation as a stochastic sequence of reaction events, including:

- ribosome binding to the mRNA ribosome-binding site
- release of the ribosome-binding site after initiation
- codon-by-codon ribosome readthrough
- polypeptide elongation
- ribosome release upon reaching the end of the coding region

The framework also accounts for ribosome footprint constraints, allowing the model to capture steric interference and stalling between neighboring ribosomes.

---

## What This Repository Demonstrates

This repository demonstrates:

- mechanistic modeling of a biological process at the sequence level
- stochastic simulation of biochemical events using a direct-event framework
- mapping from codon identity to reaction constants
- integration of sequence analysis with dynamic simulation
- comparison of alternative codon arrangements under a unified simulation setting
- extraction of quantitative outputs such as translation time and stalling counts
- early use of Python for computational biology and scientific programming

---

## Main Components

### 1. Codon Usage Analysis in *E. coli*
The project includes a script for analyzing codon usage in *E. coli* K12 MG1655 by parsing coding sequences and counting codon frequencies across the genome. This provides a biological basis for thinking about codon bias and translation efficiency.

### 2. Stochastic Simulation of Translation
The core simulation models translation on a single mRNA as a stochastic process. At each step, reaction propensities are computed, one event is sampled, and the system state is updated accordingly. The simulation tracks:

- free ribosomes
- mRNA ribosome-binding-site availability
- ribosome–mRNA complexes at different codon positions
- event time increments
- codon positions occupied by ribosomes

### 3. Sequence-Specific Translation Modeling
For a given mRNA sequence, the code can:

- identify the start codon
- translate codons into amino acids
- assign codon-specific reaction constants
- simulate peptide elongation along the sequence
- record amino-acid generation and translation trajectories

### 4. Slow-Codon Arrangement Analysis
A major focus of the project is to examine how different arrangements of slow codons affect translation outcomes. The code supports:

- random selection of slow-codon positions
- generation of predefined slow-codon patterns
- binary encoding/decoding of codon arrangements
- repeated simulations across alternative codon configurations
- comparison of mean translation time and ribosome stalling frequency

---

## Modeling Approach

The translation process is represented as a set of stochastic reactions, including:

1. **Ribosome binding** to the mRNA ribosome-binding site  
2. **Initiation progression** and release of the binding site  
3. **Elongation**, modeled as ribosome movement from one codon position to the next  
4. **Termination**, where the ribosome leaves the transcript after reaching the final codon  

Reaction selection is based on **propensity calculations**, and waiting times are sampled from an exponential distribution. This event-driven setup allows the model to capture both translation timing and competition among simultaneously translating ribosomes.

Some versions of the simulation explicitly incorporate a **ribosome footprint constraint**, preventing downstream movement when the next codon region is already occupied and thereby enabling measurement of **stalling events**.

---

## Example Questions Addressed

This project was designed to explore questions such as:

- How does codon usage relate to expected translation efficiency?
- How do different slow-codon arrangements affect the time needed to translate a protein?
- Does clustering slow codons produce more ribosome stalling than evenly distributed slow codons?
- How does ribosome number affect total protein output time?
- How do sequence-dependent elongation rates shape translation dynamics?

---

## Repository Structure

```text
.
├── Analysis of Codon Usage in E. coli K12 MG165.py
├── Direct Method_Protein Translation.py
├── Direct Method_Different Codon Sequence.py
├── Direct Method_Different Codon Sequence_Revised.py
├── Exercise 6-1 Simulation of Translation-...ipynb
├── Exercise 6-2 Simulation of Translation- OOP.ipynb
├── Exercise 6-3 ... slow codons ... footprint of Rib is 10 codons.ipynb
├── Exercise 6-4 ... random SC arrangements.ipynb
├── Exercise 6-5 ... random SC arrangements.ipynb
└── README.md
