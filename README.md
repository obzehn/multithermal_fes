# Multi-thermal and funnel corrections for free energy estimates
This is a modification of the original Invernizzi's [FES_from_Reweighting](https://github.com/invemichele/opes/blob/master/postprocessing/FES_from_Reweighting.py) script to add funnel correction and multi-thermal reweighting. For an Opes Expanded introduction, please read the [original work](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.041034).
The main additions are
1. You can now re-weight Opes Expanded simulations in the case of multi-termal (multi-canonical) simulations
2. You can now correct the FES of a funnel-restrained simulation



**References**
. For multi-thermal implementation
Invernizzi, M., Piaggi, P. M., & Parrinello, M. (2020). Unified approach to enhanced sampling. Physical Review X, 10(4), 041034. [DOI:10.1103/PhysRevX.10.041034](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.041034)
. For an implementation of the funnel correction (eq. (4) of the Supplementary information)
Rizzi, V., Bonati, L., Ansari, N., & Parrinello, M. (2021). The role of water in host-guest interaction. Nature Communications, 12(1), 93. [DOI:10.1038/s41467-020-20310-0](https://www.nature.com/articles/s41467-020-20310-0)
