# Project 5: Genome Architecture and the Evolution of Micronuclear Mutation Rates in *Tetrahymena*

## Overview

This project investigates the evolutionary causes of the exceptionally low **micronuclear (MIC) mutation rate** reported in *Tetrahymena thermophila*. The work combines stochastic population-genetic modeling with comparative simulations of **mitosis**, **amitosis**, **obligate sex**, and **facultative sex** to examine how *Tetrahymena*’s unusual genome architecture shapes both **mutation load** and **mutation-rate evolution**.

The central questions addressed in this repository are:

1. Why does the MIC carry a much lower mutation load than expected under a naive model in which MIC mutations remain invisible to selection until sex?
2. Can the distinctive reproductive strategies of *Tetrahymena*, especially **amitosis** and **facultative sex**, actively promote the evolution of lower mutation rates?
3. Does the polyploid macronuclear architecture of *Tetrahymena* enhance the reduction of mutation rate relative to a diploid system?

More broadly, this repository examines how **nuclear dimorphism, copy-number amplification, and life-cycle structure** interact with selection on mutation-rate modifiers.

## Scientific Motivation

*Tetrahymena* possesses two functionally distinct nuclei within a single cell: a transcriptionally active **macronucleus (MAC)** and a transcriptionally silent **micronucleus (MIC)**. During asexual reproduction, the MIC divides by mitosis while the MAC divides by amitosis; during sex, a diploid zygotic nucleus gives rise to both the new MIC and the new MAC. This architecture creates an unusual separation between the nucleus that determines phenotype during growth and the nucleus that transmits genetic information through sexual reproduction.

Previous work suggested that the extremely low MIC mutation rate might result from the fact that MIC mutations accumulate during long asexual intervals and are not exposed to selection until the next sexual cycle. This project tests that idea directly and shows that the explanation is incomplete: although MIC mutations are not expressed during asexual generations, they can still respond indirectly to selection because, after sex, the same mutations are present in the newly formed MAC, where selection acts on them during subsequent asexual growth.

The project then moves beyond this negative result and asks a broader question:

**Can the reproductive system of *Tetrahymena* itself favor the evolution of low mutation rates?**

## Project Goals

The main objectives of this repository are to:

- quantify the mutation load carried in the MIC under different frequencies of sexual reproduction;
- test whether MIC mutations respond indirectly to selection during asexual generations;
- investigate the evolution of mutation rate using explicit **mutator** and **antimutator** modifier loci;
- compare mutation-rate evolution under **mitosis**, **amitosis**, **obligate sex**, and **facultative sex**;
- distinguish the effects of **sexual processes** from those of **amitotic asexual cycles**;
- compare mutation-rate evolution in **45-ploid** and **diploid** MAC systems.

## Core Findings

This repository centers on four main conclusions.

### 1. The MIC carries a much lower mutation load than expected

A naive expectation assumes that MIC mutations accumulate during asexual generations without selection and are only exposed after sex. The simulations show that this prediction is too pessimistic: the equilibrium fitness immediately after sex is substantially higher than expected, indicating that the MIC carries a much lower mutation load than the naive theory predicts.

### 2. MIC mutations respond indirectly to selection during asexual generations

The explanation is that MIC mutations generated before sex are inherited into the newly formed MAC after sex. Because the MAC is expressed during subsequent asexual growth, selection acting on MAC fitness can indirectly change the frequencies of the corresponding mutations in the MIC. This indirect selection remains robust even when MAC mutations are also allowed to occur.

### 3. Amitosis and facultative sex can reduce mutation rates

When explicit mutation-rate modifiers are introduced, the simulations show that certain reproduction strategies can promote the evolution of lower mutation rates. Under deleterious-only mutation, **amitosis** lowers the MAC mutation rate, while **facultative sex** lowers mutation rates in both the MAC and the MIC. In contrast, mitosis does not reduce mutation rate effectively, and obligate sex does not produce the same overall reduction.

### 4. Polyploidy enhances the reduction of mutation rate

Under facultative sex, the **45-ploid MAC system** achieves a greater reduction in mutation rate than the corresponding diploid system. This suggests that the polyploid genome architecture of *Tetrahymena* strengthens the ability of selection to reduce mutation rate under the right reproductive regime.

## Biological Context

This project builds on the broader idea that *Tetrahymena*’s success may depend on its unusual combination of:

- **nuclear dimorphism** (separate MIC and MAC),
- **polyploid somatic genome structure**,
- **amitotic MAC division with copy-number control**,
- and **facultative sexual reproduction**.

Here, the focus is not only on short-term fitness consequences but also on how these features shape the **evolution of mutation rate itself**. In this sense, the project connects genome architecture to long-term evolutionary optimization.

## Repository Contents

This folder contains both simulation scripts and figure-generation notebooks related to mutation load and mutation-rate evolution.

### A. Mutation load in the MIC and effects of sexual frequency

- `File 5-1.  Code_Check the Equilibrium of Facultative Sex.py`  
  Simulates populations under different sexual frequencies and checks whether equilibrium is reached when fitness is evaluated immediately after sex.

- `File 5-2a. Code_Turn Off Mutation in MAC_Facultative Sex with Amito.py`  
  Tests whether MIC mutations respond indirectly to selection by turning off MAC mutation and allowing mutations to arise only in the MIC prior to sex, with amitotic asexual cycles.

- `File 5-2b. Code_Turn Off Mutation in MAC_Facultative Sex with Mito.py`  
  Parallel control with mitotic asexual cycles.

- `File 5-3.  Code_Same MT but different Selection Freq_Ploidy 45 FS with Amito.py`  
  Compares settings with matched mutation input but different opportunities for selection to act between sexual cycles.

- `File 5-4. Code_Same # of MT During Two Sex but Diff Sexual Freq_P45.py`  
  Additional comparison of sexual period while controlling mutation accumulation between two sexual events.

- `File 5-5. Code_Control_Comparison of Different MT Strategy But Same MT #.py`  
  Control simulations for alternative mutation-timing strategies.

- `File 5-6. Code_Turn on MAC mutation again_Soma MTU =22.5 During 2 Sexes.py`  
  Tests the robustness of indirect selection on MIC mutations when MAC mutations are reintroduced.

### B. Mutation-rate evolution with modifier loci

- `File 5-7a. Code_Evolution of Mu Rate_10 Upper and 10 Lower Ctrl_AllDeleMut.py`  
  Simulates mutation-rate evolution with 10 mutator and 10 antimutator loci under deleterious-only mutation.

- `File 5-7b. Code_Evolution of Mu Rate_10 Upper and 10 Lower Ctrl_1PerBeneMut.py`  
  Extends mutation-rate evolution to a regime including rare beneficial mutations.

These scripts implement the key Chapter 4 framework in which mutation-rate modifiers do not directly affect fitness, but instead alter mutation rates at fitness loci and are therefore indirectly selected through their covariance with fitness.

### C. Figure-generation and result-summary notebooks

- `File 5-8. Fig 4.1. Plot fig for N =500_Having FS_RM_Revised.ipynb`
- `File 5-9. Fig 4.1. Plot fig for N =500_Having FS_RM_Revised_New.ipynb`  
  Plotting notebooks for the comparison between simulated post-sex equilibrium fitness and naive mutation-load expectations.

- `File 5-10. Fig 4.2. Diff Sexual Freq Same Mut_P45 with Amito.ipynb`
- `File 5-11. Fig 4.2. Diff Sexual Freq Same Mut_P45 with Amito_New.ipynb`  
  Figure notebooks for testing how MIC mutations respond to selection under different sexual frequencies with matched mutation input.

- `File 5-12. Fig 4.3. Plot Figure for Same MT but Diff Sexual Freq_RM_MAC MT =45.ipynb`
- `File 5-13. Fig 4.3. Plot Figure for Same MT but Diff Sexual Freq_RM_MAC MT =45_New.ipynb`  
  Plotting notebooks for the robustness analysis with MAC mutation turned back on.

- `File 5-14. Fig 4.8. Test of Amitosis_N = 2K_Final Rate.ipynb`  
  Evaluates mutation-rate evolution under amitosis.

- `File 5-15. Fig 4.9. Test of Facultative Sex_N = 2K_Final Rate.ipynb`  
  Evaluates mutation-rate evolution under facultative sex.

- `File 5-16. Fig 4.11. Test of FS with Amito_N = 2K_Final Rate_P = 2 vs 45.ipynb`  
  Compares facultative sex in diploid and 45-ploid systems.

- `File 5-17. Fig 4.13 and 14. Test of FS with Amito_N = 2K_Bene_P =2 vs 45_New data.ipynb`  
  Polyploid-versus-diploid comparison in the presence of beneficial mutations.

- `File 5-18. Fig 4.13 and 14. Test of FS with Amito_N = 2K_Dele_P = 2 vs 45.ipynb`  
  Polyploid-versus-diploid comparison under deleterious-only mutation.

- `Chapter 4. Evolution of Mu_Ver3_Further Rv.docx`  
  Research chapter describing the biological motivation, methods, results, discussion, and theoretical interpretation of this project.

## Methods and Approach

### Stochastic population-genetic framework

The simulations extend a stochastic framework for mutation, selection, and reproduction in a dual-nucleus *Tetrahymena*-like system. Populations evolve with explicit handling of:

- a diploid **MIC**;
- an **n-ploid MAC** (usually 45-ploid);
- asexual cycles via **mitosis** or **amitosis**;
- and sexual cycles via random mating and zygote formation.

### Quantifying MIC mutation load

To assess mutation load in the MIC, the simulations monitor **population mean fitness immediately after sex**, when MIC-derived mutations become expressed in the newly formed MAC. This provides a direct numerical measure of the hidden burden carried by the MIC.

### Testing indirect selection on MIC mutations

A key design feature is the set of simulations in which mutation is turned off in the MAC and allowed only in the MIC immediately before sex. This isolates the effect of indirect selection on MIC mutations during subsequent asexual generations. Additional controls then reintroduce MAC mutation to test robustness.

### Mutation-rate modifier model

The mutation-rate evolution component introduces **mutator** and **antimutator** loci that alter the genomic mutation rate but do not directly contribute to organismal fitness. Their evolutionary dynamics therefore emerge indirectly through associations with the fitness consequences of mutations at other loci. The model includes both deleterious-only and mixed beneficial/deleterious settings.

### Polyploid versus diploid comparisons

The repository includes dedicated simulations comparing **45-ploid** and **diploid** MAC systems under facultative sex. These comparisons test whether polyploid architecture changes the efficiency with which mutation rates can evolve downward.

## Scientific Interpretation

A major conceptual contribution of this project is to challenge a simple narrative about hidden germline mutations in ciliates. The MIC is not merely a passive reservoir of mutations that escapes selection until the next round of sex. Instead, because MIC and MAC are jointly derived from the same zygotic nucleus, selection acting on MAC performance feeds back indirectly onto the evolutionary fate of MIC mutations.

The project also shows that mutation-rate evolution depends strongly on **life-cycle structure**. In this framework:

- **mitosis** does not efficiently reduce mutation rate;
- **amitosis** can reduce the **MAC** mutation rate;
- **facultative sex** can reduce mutation rates in both **MAC** and **MIC**;
- and **polyploidy** strengthens this reduction.

These results suggest that *Tetrahymena*’s unusual reproductive and genomic system may itself help explain why extremely low mutation rates are evolutionarily achievable.

## Technical Highlights

This repository demonstrates several strengths characteristic of computational evolutionary biology:

- explicit modeling of dual-nucleus genome architecture;
- simulation of mutation accumulation and mutation-load dynamics across sexual frequencies;
- implementation of evolvable mutation rates through mutator and antimutator loci;
- comparative analysis of mitosis, amitosis, obligate sex, and facultative sex;
- separation of mutation-input effects from selection-frequency effects through control simulations;
- direct comparison of diploid and polyploid systems under shared assumptions.

## Broader Relevance

Although centered on *Tetrahymena*, this project addresses broadly relevant questions in evolutionary genetics:

- how mutation load behaves in organisms with complex life cycles;
- how mutation-rate modifiers evolve under different inheritance systems;
- how hidden or delayed expression of mutations interacts with selection;
- and how genome architecture shapes the long-term evolution of genomic fidelity.

More generally, the repository illustrates a research workflow that integrates biological hypothesis testing, stochastic simulation, comparative experimental design in silico, and figure-oriented quantitative analysis.


## Note

This repository is presented as a research project archive documenting simulation-based analyses of mutation load and mutation-rate evolution in a dual-genome ciliate system. Some scripts and notebooks reflect iterative research-stage development and may require minor cleanup or environment-specific adjustment for re-execution in newer software versions. The data used to generate the figures can be reproduced by running the corresponding simulation and plotting files included in this repository.