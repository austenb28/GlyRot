# GlyRot
Takes an N-glycosylated glycoprotein in Gromacs generated with the Glycam forcefield and systematically generates a feasible configuration of N-glycans via rigid rotation about three bonds per N-glycan.

# Requirements
* Python3 or greater
* Gromacs 5.1

# Installation
Download and add the bin directory to `$PATH`.  For example, in Linux, add

`export PATH="/path/to/GlyRot/bin":$PATH`

in `~/.bashrc` and open a new terminal.

A precompiled binary using Gromacs 5.1 is provided in the bin directory. GlyRot must be executed while Gromacs 5.1 is installed.  If the GlyRot exectuable must be rebuilt, a Makefile is provided in the top directory, which can be built using `make`.

# Arguments
Execute `GlyRot.py -h` for a brief description of arguments.  Defaults values are provided  for optional arguments.

# Example
Starting in the top directory:

`cd example`

`GlyRot.py -gro BChEG_GMX.gro -top BChEG_GMX.top -dtheta 90`.

This assumes the gromacs executable is `gmx`. Otherwise, the `-gex` option should be used to specify a different gromacs executable name.  If GlyRot was successful, the file `Rotated.gro` is generated.  One can verify the glycans are sterically clashed in `BChEG_GMX.gro`, and not in `Rotated.gro`.  The rotation resolution `-dtheta 90` is used here to permit rapid execution on a personal computer.  For actual system preparation, the default `-dtheta 15` is advised, which generates glycan configurations at a rate of approximately 1 glycan per 15 minutes on a single CPU.  Additional CPUs will not significantly improve this speed.  BChEG_GMX.gro contains 9 glycans, and was used in [this](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5708630/) study.

An example script for executing `GlyRot.py` using Slurm on an HPC cluster is provided in srun_rot.sh. The script is executed using `sbatch srun_rot.sh`.

# Additional Information
GlyRot is executed on Gromacs coordinate and topology files using the Glycam forcefield.  These input files may be generated using [Acpype](https://github.com/alanwilter/acpype/) on a a system built using [Glycam Web](http://glycam.org/) and [AmberTools Leap](https://ambermd.org/AmberTools.ph). 

The three bonds that the N-glycans are rotated about are `ND2 - C1`, `CB - CG`, `CA - CB`.  Atom names provided using the Glycam forcefield.  All atoms correspond to the connected asparagine, except `C1`, which corresponds to the first attached sugar of the corresponding glycan.  Each glycan is sequentially rotated using all rotamers about the three bonds, totaling (360/`dtheta`)<sup>3</sup> rotamers per glycan.  The minimum in-vacuo energy rotamer for each glycan is selected at each stage.  This method has been shown to generate feasable glycan conformers for even dense glycoproteins, such as butyrylcholinesterase.
