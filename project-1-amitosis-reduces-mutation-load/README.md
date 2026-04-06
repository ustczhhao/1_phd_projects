# Project 1: Amitosis Reduces Mutation Load Compared to Mitosis

## Overview

This project investigates a central question in evolutionary genetics: whether **amitosis can confer some of the evolutionary benefits typically associated with sex, even in the absence of sexual reproduction**. Using a combination of **mathematical modeling, symbolic derivation, and computational analysis**, the project examines how amitotic division can reduce deleterious mutation load relative to mitosis under mutation-selection balance.

The biological motivation comes from *Tetrahymena*, a ciliate with separate germline and somatic nuclei, whose amitotic somatic division provides a distinctive system for studying the evolutionary consequences of noncanonical inheritance. This repository focuses specifically on the **deterministic fitness benefit** of amitosis and documents the analytical and numerical workflow used to derive, validate, and visualize that result.

More broadly, the project illustrates how a theoretical problem in evolutionary biology can be translated into a reproducible computational framework.

## Scientific Motivation

Sex is widely thought to persist because it improves the efficiency of selection, facilitates the removal of deleterious mutations, and enhances adaptive potential. Amitosis, by contrast, is often viewed as a less precise mode of nuclear division. In this system, however, amitotic segregation can generate heritable variation among offspring in a way that strengthens the action of natural selection.

This project addresses the following question:

**Can amitosis reduce mutation load compared with mitosis, even without recombination?**

The answer is relevant not only to the evolution of *Tetrahymena*, but also to broader questions concerning genome architecture, inheritance mode, mutation load, and long-term lineage persistence.

## Project Goals

The main objectives of this repository are to:

- derive and verify the equilibrium mean fitness under mitosis and amitosis;
- quantify the deterministic fitness advantage of amitosis across mutation and selection parameters;
- reproduce representative numerical results from the underlying study;
- connect analytical results to interpretable computational outputs and figures;
- preserve a transparent record of the modeling and analysis workflow.

## Core Result

Under the deterministic mutation-selection framework, amitosis reduces mutation load relative to mitosis by generating additional heritable variance among offspring. As a result, amitotic reproduction can maintain a higher equilibrium mean fitness than mitotic reproduction under the same deleterious mutation regime.

This provides a clear theoretical example of how an unusual mode of asexual division can partially mimic one of the classic advantages of sex: **more effective purging of deleterious mutations**.

## Repository Contents

This folder contains both symbolic and numerical analysis files used to study the deterministic benefit of amitosis:

- `direct_amitosis.nb`  
  Mathematica notebook for analytical treatment of the amitotic model.

- `direct_mitosis.nb`  
  Mathematica notebook for the corresponding mitotic baseline model.

- `mitosis_analyt.nb`  
  Mathematica notebook containing analytical derivations for the mitotic framework.

- `mitosis_analyt2.nb`  
  Additional Mathematica derivations and validation related to the mitotic analytical results.

- `Equilibrium fitness under amitosis and mitosis_Ud = 0.1_S = -0.1.ipynb`  
  Python notebook for numerical evaluation and visualization of equilibrium fitness under a representative parameter setting.

- `Deterministic benefit of amitosis over mitosis.ipynb`  
  Python notebook for plotting and exploring the deterministic selective advantage of amitosis over mitosis.



## Methods and Approach

### Mathematical Modeling

The project begins with a population-genetic model of asexual reproduction under deleterious mutation and selection. Separate formulations are developed for mitosis and amitosis, allowing direct comparison of genotype dynamics and equilibrium mean fitness.

The analytical framework focuses on deriving tractable equilibrium results under explicit assumptions about mutation, selection, ploidy, and inheritance mode.

### Symbolic Derivation

Mathematica notebooks are used to derive, simplify, and verify the expressions governing equilibrium behavior. This symbolic layer ensures that the numerical results are grounded in explicit mathematical reasoning rather than simulation alone.

### Numerical Analysis and Visualization

Python notebooks are used to evaluate representative parameter settings, compare equilibrium fitness under alternative reproductive modes, and generate visual summaries of the deterministic benefit of amitosis.

Together, the symbolic and numerical components form a workflow that is both rigorous and interpretable.

## Technical Highlights

This project demonstrates several strengths characteristic of quantitative and computational research:

- mathematical modeling of biological systems under explicit assumptions;
- symbolic computation to support analytical derivation and validation;
- scientific programming for numerical evaluation and visualization;
- theory-driven reasoning linking formal results to biological interpretation;
- reproducible organization of derivations, analyses, and outputs in notebook form.

## Tools and Environment

The project uses a mixed analytical and computational workflow built around:

- **Mathematica** for symbolic derivation and analytical verification;
- **Python / Jupyter Notebook** for numerical evaluation, exploration, and figure generation.

This structure is well suited to projects that require both formal theoretical development and computational implementation.

## Broader Relevance

Although rooted in evolutionary genetics, the framework developed here is broadly relevant to quantitative biology. The project demonstrates how to formulate a biological problem mathematically, derive interpretable results, validate them computationally, and present them in a form suitable for scientific communication.

The underlying workflow is transferable to a wide range of problems involving theoretical modeling, computational biology, population genetics, and quantitative analysis.

## Biological Context

The broader context of this work is the evolution of amitosis in *Tetrahymena*. Because amitotic division introduces variation in allele copy number among daughter nuclei while maintaining overall ploidy control, it can generate sufficient variance for selection to act more effectively than under strict mitosis.

This makes the system a particularly interesting natural example of how nonstandard inheritance mechanisms can shape mutation load, adaptation, and long-term evolutionary dynamics.

## Possible Extensions

Natural next steps for this project include:

- expanding the analysis to additional ploidy levels or broader parameter sweeps;
- comparing deterministic and stochastic effects within a unified framework;
- reproducing additional manuscript figures in a consistent plotting style;
- refactoring the notebooks into a more modular and reusable analysis pipeline.


## Note

This repository is presented as a research project archive documenting the analytical and computational workflow for a theoretical evolutionary biology study. Some notebooks reflect research-stage code organization and may require minor cleanup or environment-specific adjustment for re-execution in newer software versions.