# Stochastic simulation of gene expression

## Introduction
This package implements an algorithm to simulate the stochastic gene expression in an (optionally) growing cell.
The user can specify the desired reactions and the coupling to the cell growth.

If the cell is allowed to growth, the cell draws a growth rate form a gaussian distribution and grows until it has added a specific volume, also drwan from a gaussian. At division, the volume is halved and the daughter receives a binomial sample of the mother's molecules.


## Specify the reactions
The reactions are specified in a text file that must be passed to the program. You can find an example in reactions.txt.
For example, let's suppose we want to simulate a growing cell that produces mRNA, unfolded GFP and mature GFP.
Here are the settings to be used

* Every line beginning with # is a comment and is therefore ignored.

* !population: here we give the initial amount of each species. In our case we just write 1 0 0 to show that at the beginning we want 1 mRNA and zero proteins (in a latare section we can give names to the reactions).

* !matrix: Here we write the stochiometric matrix. The last element must be a semicolon.

* !rates: the reaction rates of the reactions (in a later section we can give name to the reactions)

* !reactants: here we define the reactant list. Every line is a different reaction. For example the first reaction is mRNA production and this doesn't depend on any amount (it is a source reaction). We can show this by writing -1 for a standard source reaction, or -2 if we want it to be proportional to the cell volume. In our example we suppose that the mRNA production rate depends on the cell volume, so we write -2 on the first line. The second line is mRNA decay which depends on the amount of mRNA. Since mRNA is the species number 0 (it's the first specie we wrote in the !population section), we write 0 on the second line. The last reaction is GFP bleaching, which depends on the amount of mature GFP, which is the species number 2, we write 2 on the last line.

* !time_step: usually experiments return data at discrete time step. If we want the results every 3 time units, we write 3.

* !species_names: here we have the possibility to write the names for each species. We call them m (mRNA), p (protein), g (folded GFP); so we write m p g

* !time_range: the time range we want to report. If we write 10000 30000 it means that we want to wait 10000 time units before starting outputting the results (for example because we want to reach steady state) and we want to simulate until we reach 30000 time units.

* !cell_growth: these are the parameters of the cell growth: initial volume, mean and standard deviation of the normal distribution from which to draw the growth rate and the amount of added volume before division.

