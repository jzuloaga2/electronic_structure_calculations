# electronic_structure_calculations
Fortran code to calculate electronic structure of and light absorption for a given geometrical configuration.



**************************************************
Files:
**************************************************

eqdavidson.f :
    Computes the equilibrium electron density
    by solving the Kohn-Sham equations of the system.
    We start with a plane-wave approximation and iterate
    until convergence.
    Atomic potentials are set in the jellium approximation.
    The output of
    this program, namely, the equilibrium electron density and the
    self-consistent potential of the system, are plotted by the matlab
    file equilibrium.m.

testrpa.f :
    Uses the output of eqdavidson.f to calculate the absorption of the
    system in the Random Phase Approximation.  It calls the
    subroutines: chartree.f, chi0.f, green.f, solverpa.f .


**************************************************
Subroutines: 
**************************************************

solverpa.f
    Subroutine called by testrpa.f, calculates
    absorption in the Random Phase Approximation.

chi0.f
    Subroutine called by testrpa.f and by solverpa.f, calculates chi0.

chartree.f
    Subroutine called by solverpa.f to calculate the Hartree potential.

green.f 
    Subroutine called by chi0.f to calculate the Green's function.



