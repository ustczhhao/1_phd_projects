# Stochastic Simulation of Chemical Reaction Systems  
### Reimplementation of Gibson & Bruck (2000) in Python

This project contains a Python-based reimplementation of core stochastic simulation algorithms described in Gibson and Bruck’s landmark paper on exact simulation of chemical systems with many species and reaction channels.

The repository was originally developed as a PhD-era coding project to simulate the stochastic dynamics of a multi-species chemical reaction network. It includes both procedural and object-oriented implementations of two classical stochastic simulation approaches:

- **Direct Method**
- **First Reaction Method**

In addition to reproducing the core simulation logic, the project also explores code modularization, reaction abstraction, equilibrium detection, and time-window-based sampling of system states.

---

## Project Background

Stochastic simulation algorithms (SSAs) are widely used to model chemical and biochemical systems in which random molecular interactions play an essential role. Rather than approximating concentrations deterministically, these methods explicitly simulate reaction events in time based on reaction propensities.

This repository was built to study and implement the simulation framework described in:

**Gibson, M. A., & Bruck, J. (2000).**  
*Efficient exact stochastic simulation of chemical systems with many species and many channels.*  
**The Journal of Physical Chemistry A, 104(9), 1876–1889.**

The code simulates the stochastic evolution of a reaction system involving multiple chemical species and reaction channels, with exact event timing determined by probabilistic rules.

---

## What This Repository Demonstrates

This project is more than a simple code exercise. It demonstrates:

- implementation of **exact stochastic simulation algorithms**
- decomposition of reaction systems into reusable computational components
- transition from **procedural scripting** to **object-oriented design**
- use of **Monte Carlo simulation logic** for event-driven systems
- tracking of **species counts over time**
- detection of **equilibrium states**
- extraction of system states at predefined **time windows**
- early scientific programming practice in **Python 3**, with numerical computing and visualization

---

## Implemented Methods

### 1. Direct Method
The Direct Method simulates one reaction event at a time by:

1. computing the propensity of each reaction,
2. sampling which reaction occurs next based on normalized propensities,
3. drawing the waiting time from an exponential distribution,
4. updating species counts accordingly,
5. iterating until the end of the simulation horizon or equilibrium is reached.

This implementation includes:
- propensity calculation
- reaction-state updates
- repeated simulation over time
- equilibrium detection
- optional window-based sampling of trajectories

---

### 2. First Reaction Method
The First Reaction Method independently samples a tentative firing time for each reaction channel and selects the reaction with the smallest putative reaction time.

This implementation includes:
- per-reaction tentative time generation
- minimum-time reaction selection
- exact state update logic
- repeated time evolution of the system
- extraction of discrete states across simulation intervals

---

### 3. Object-Oriented Refactoring
The repository also includes object-oriented versions of both methods. These versions introduce reusable abstractions for:

- **ChemicalSpecies**
- **Reaction**
- **System**

This refactoring improves modularity and readability, and makes it easier to extend the framework to larger or more general reaction systems.

---

## Repository Structure

```text
.
├── 1. Gibson Direct Method Code.py
├── 1.1. Gibson Direct Method_OOP.py
├── 2. Gibson First Reaction Method Code.py
├── 2.1. Gibson First Reaction Method_OOP.py
├── Exercise 5-1- Direct Method- New-01152016-Window.ipynb
├── Exercise 5-1-Direct Method- New-01152016-Equiltime- With Documentation.ipynb
├── Exercise 5-2 - First Reaction Method-New-12182015-Equiltime - With Documentation.ipynb
├── Exercise 5-2 - First Reaction Method-New-12182015-Window-With Documentation.ipynb
├── Exercise 5-3 -Object oriented programming- 01272016.ipynb
├── Exercise 5-3 -Object oriented programming- 01282016.ipynb
├── Exercise 5-3 -Object oriented programming- Directmethod- 02052016-With Documentation.ipynb
├── Exercise 5-3 -Object oriented programming- Directmethod- 02062016-With Documentation-Use Set.ipynb
├── Exercise 5-4 -Object oriented programming- First Reaction Method- 02052016-With Documentation.ipynb
└── Exercise5-1_Simulation of Gibson Paper_Store the dynamics and plot fig.ipynb
