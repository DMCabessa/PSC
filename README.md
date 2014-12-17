# PSC

This is my personal *toolbox* for Particle Swarm Optimization in classification problmens. I have started from the [PSO toolbox](https://code.google.com/p/psomatlab/) and edited it to use the algorithm in classification problems. Feel free to use, edit and distribute this *toolbox* as log as you reference this source. Bear in mind that this is an algorithm **under construction**, so there will be changes.

## Contents

* The `/data` folder contains the datasets used in my tests. The file `pscdemo.m` loads files from this folder.
	+ Sources with training and test sets previously divided are listed as `.data` and `.test` respectively. Passing the argument `SEGMENTED` to the `psc.m` function loads both files separately.
	+ Sources with no segmentation are listed as `.data` only and are divided into training and test sets if you pass the argument `SINGULAR` to the `psc.m` function.
* The `/private` folder contains some functions, taken from the original source.
* On the main folder you will find many auxiliar functions as well and the 3 types of fitness functions (`pscfitnessfcn1.m`, `pscfitnessfcn2.m` and `pscfitnessfcn3.m`). The scripts are for testing purpouses, to clean the functions on RAM and loading changes on every run.
