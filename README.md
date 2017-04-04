# mutated-godmd

Based on GOdMD (Sfriso, P, et al Bioinformatics 2013)

Code to find transition paths between two known structures mimicking the effect of mutations.

#### Requirements
* gfortran compiler

#### How to compile it

1. Modify the config.mk file. Change the variable `F77` variable to point to the gfortran compiler
2. Run `make clean` in the main directory
3. Run `make` in the main directory

#### Input files

 * -i        < parameter file >
 * -pdbin    < initial structure in PDB format >
 * -pdbtarg  < target structure in PDB format >
 * -p1       < alignment file, for initial structure>
 * -p2       < alignment file, for target structure>
 * -touch    < residues to mutate >
 * -ener     Output, < energy file >
 * -trj      Output, < trajectory file >
 * -o        Output, < log file >


#### How to run it

The exe `godmd` file is in the ./exe directory

##### In the example directory run:
$ ../exe/godmd -i param.in -pdbin 1ake_A.pdb  -pdbtarg 4ake_A.pdb -p1 1ake_A.aln -p2 4ake_A.aln -touch mutations.dat -trj test.trj.crd -ener test.ener.dat -o log_1ake.txt


##### Mutations format

33   50.

32   50.

34   50.

Should be interpreted as: Resiudes 33, 32, 34 will have they energy interactions 50 times more stable.
