# Multi-thermal and funnel corrections for free energy estimates
This is a modification of the original Invernizzi's [FES_from_Reweighting](https://github.com/invemichele/opes/blob/master/postprocessing/FES_from_Reweighting.py) script to add funnel correction and multi-thermal reweighting. For an Opes Expanded introduction, please read the original work<sup>1</sup>.

The main additions are
1. You can now re-weight Opes Expanded simulations in the case of multi-termal (multi-canonical) simulations[^1]
2. You can now correct the FES of a funnel-restrained simulation[^2]

To re-weight for multi-thermal runs, two flags are needed
* **`--multiT`** The target temperature (in K)
* **`--ene`** The potential energy of the system
Remember to add the multi-thermal bias to the other biases you are applying, in case your naming is different from the standard `.bias` of PLUMED.

To correct for funnel restraint, three flags are needed
* **`--rfunnel`** the radius of the funnel in nm
* **`--uat`** the value of the z CV above which the ligand is unbound
* **`--bat`** the value of the z CV above which the ligand is bound

A couple of other minor points have been addressed in the script
1. Additional output is reported in the headers of the output files (input file name, sigmas, etc.)
2. The flag `--blocks` associated with `--deltaFat` now reports the average deltaF and its associated error in the header of the main output file as weighted averages and standard deviations as in eqns. B3 and B4 of the original work[^1]

[^1]: Invernizzi, M., Piaggi, P. M., & Parrinello, M. (2020). Unified approach to enhanced sampling. Physical Review X, 10(4), 041034. [DOI:10.1103/PhysRevX.10.041034](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.041034)

[^2]: Rizzi, V., Bonati, L., Ansari, N., & Parrinello, M. (2021). The role of water in host-guest interaction. Nature Communications, 12(1), 93. [DOI:10.1038/s41467-020-20310-0](https://www.nature.com/articles/s41467-020-20310-0)
