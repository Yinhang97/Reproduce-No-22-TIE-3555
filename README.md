# Reproduce-No-22-TIE-3555
code and data for reproducing paper No-22-TIE-3555

The cls_machine.m is a class to define the basis parameters of machine

The force_struct.m is a struct which includes the radial and tagential components of forces, which are given by sampling sequences of force field.

The mode_struct.m is a struct which includes the radial and tagential components of modes, which are given by sampling sequences of mode field.

The original data for obtaining force field and mode field has been provided in "mode_vectors" folder and "nodal_forces" folder. The calculation method has been deduced in detail in the paper.

The vib_synthesis_for_submission.m is the code for vibration synthesis. The program can be started by running vib_synthesis_for_submission.m.

The program was tested in MATLAB R2021b.
