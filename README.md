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

## Multi-thermal example
Let's suppose we are running a typical di-alanine simulation with the following `plumed.dat` file
```
MOLINFO STRUCTURE=ref_ala2.pdb
phi: TORSION ATOMS=@phi-2
psi: TORSION ATOMS=@psi-2
ene: ENERGY
ecv: ECV_MULTITHERMAL ARG=ene TEMP_MIN=300 TEMP_MAX=1000
opesX: OPES_EXPANDED ARG=ecv.* FILE=DELTAFS PACE=500
PRINT FMT=%g STRIDE=500 FILE=COLVAR ARG=phi,psi,opesX.*
```
where we are using Opes Expanded to sample the multi-canonical ensembles of dialanine associated to the tempereature range 300K-1000K, where 300K is the thermostat reference temperature. Then the `COLVAR` will look something like this
```
#! FIELDS time phi psi ecv.ene opesX.bias
#! SET min_phi -pi
#! SET max_phi pi
#! SET min_psi -pi
#! SET max_psi pi
 0.000000 -1.545129 3.008678 -24.172577 0.000000
 1.000000 -2.847320 2.473553 -38.454010 0.000000
 2.000000 -2.271588 3.045746 -32.297424 0.000000
 3.000000 -2.722175 3.057054 -44.255707 0.000000
 4.000000 -2.545126 2.809278 -29.646790 0.000000
[...]
```
Suppose we want to re-weight the CV `phi` at 500K with `sigma=0.06` and calculate the deltaF at 0.0, i.e. the deltaF between positive and negative `phi` values. Then we can run the following
```
python3 FES_from_Reweighting_multiT_funnel.py --colvar COLVAR --outfile fes_opesX_500K --sigma 0.06 --temp 300 --cv phi --bias opesX.bias --deltaFat 0.0 --multiT 500 --ene ecv.ene
```
Notice how `--temp` is the reference temperature of the thermostat (it doesn't change as we target different temperatures), while `--multiT 500` is the temperature we are targeting with the reweighting. Equivalently, you could have passed the columns numbers rather than the columns name.

The script will notify us of the columns used as data, check it to be sure about the output
```
 using cv "phi" found at column 2
 using bias "opesX.bias" found at column 5
 using "ecv.ene" as potential energy for multiT reweighting found at column 4
```
and the output header will contain the main information about the calculation along with the final fes.
```
#! FIELDS phi file.free
#! SET input_colvar_file COLVAR
#! SET sample_size 100001
#! SET effective_sample_size 40991.1
#! SET DeltaF 12.0847 kJ/mol
#! SET DeltaF_at 0
#! SET multiT_T0 300 K
#! SET multiT_Ttarget 500 K
#! SET sigma_phi 0.06
#! SET min_phi -3.14159
#! SET max_phi 3.14159
#! SET nbins_phi 100
#! SET periodic_phi true
   -3.141593     15.203143
   -3.078761     12.614765
[...]
```

## Multi-thermal + other biasing example
Let's suppose now that we are adding some other bias potential in our `plumed.dat` file, for example
```
MOLINFO STRUCTURE=ref_ala2.pdb
phi: TORSION ATOMS=@phi-2
psi: TORSION ATOMS=@psi-2
ene: ENERGY
ecv: ECV_MULTITHERMAL ARG=ene TEMP_MIN=300 TEMP_MAX=1000
opesX: OPES_EXPANDED ARG=ecv.* FILE=DELTAFS PACE=500
opes: OPES_METAD ...
  ARG=phi,psi
  PACE=5000
  BARRIER=40
...
PRINT FMT=%g STRIDE=50 FILE=COLVAR ARG=*
```
then the output will look something like this
```
#! FIELDS time phi psi ene ecv.ene opesX.bias opes.bias opes.rct opes.zed opes.neff opes.nker
#! SET min_phi -pi
#! SET max_phi pi
#! SET min_psi -pi
#! SET max_psi pi
 0.000000 -1.54513 3.00868 -24.1726 -24.1726 0 -40 -40 1 1 0
 0.100000 -1.55767 2.89311 -21.5896 -21.5896 0 -40 -40 1 1 0
 0.200000 -1.97153 2.56139 -37.4598 -37.4598 0 -40 -40 1 1 0
 0.300000 -2.09802 2.27693 -15.6814 -15.6814 0 -40 -40 1 1 0
 0.400000 -2.43291 2.09504 -4.2428 -4.2428 0 -40 -40 1 1 0
[...]
```
We can still re-weight for different temperature, but we have to consider *both* the bias added by OPES_METAD and the one added by OPES_EXPANDED. As such, you will have to pass both biases to the script, e.g.
```
python3 FES_from_Reweighting_multiT_funnel.py --colvar COLVAR --outfile fes_opesX_500K --sigma 0.06 --temp 300 --cv phi --bias opes.bias,opesX.bias --deltaFat 0.0 --multiT 500 --ene ecv.ene
```
Again, check the information from the script
```
 using cv "phi" found at column 2
 using bias "opes.bias" found at column 7
 using bias "opesX.bias" found at column 6
 using "ecv.ene" as potential energy for multiT reweighting found at column 5
```

[^1]: Invernizzi, M., Piaggi, P. M., & Parrinello, M. (2020). Unified approach to enhanced sampling. Physical Review X, 10(4), 041034. [DOI:10.1103/PhysRevX.10.041034](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.10.041034)

[^2]: Rizzi, V., Bonati, L., Ansari, N., & Parrinello, M. (2021). The role of water in host-guest interaction. Nature Communications, 12(1), 93. [DOI:10.1038/s41467-020-20310-0](https://www.nature.com/articles/s41467-020-20310-0)
