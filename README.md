# Multi-thermal and funnel corrections for free energy estimates
This is a modification of the original Invernizzi's [FES_from_Reweighting](https://github.com/invemichele/opes/blob/master/postprocessing/FES_from_Reweighting.py) script to add funnel correction and multi-thermal reweighting. For an Opes Expanded introduction, please read the [original work](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.041034).
The main additions are
1. You can now re-weight Opes Expanded simulations in the case of multi-termal (multi-canonical) simulations
2. You can now correct the FES of a funnel-restrained simulation (from eq. (4) of the Supplementary Material of [this work](https://doi.org/10.1038/s41467-020-20310-0))
