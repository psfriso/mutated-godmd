!     
! File:   constants.f
! Author: gelpi
!
! Created on 20 de julio de 2012, 8:30
!

MODULE constants
    !sizes
    integer, parameter :: MAXATYPES = 1 ! Atom types
    integer, parameter :: MAXRES = 50, MAXATPERRES = 1 ! residue libraries
    integer, parameter :: MAXSTEPS = 4, MAXSPOTDEF = 3 ! Step Potentials

    !StepPotentials
    integer, parameter :: &
    WELL = 1   !Well centered in a position. 3 Param DELTA, and EMIN, TODO: infWall .true means always rebound ! 
    !FIT = 2, & !Fits curve A/r^n Param A, n, rcut (last step), nsteps. Step 0 es always rvdwij, use eref to scale
    !USER = 3 ! User provided: Nsteps, r and e arrays (1..nsteps) use eref to scale, 0 at step 0 means rvdwij
    integer, parameter :: &
    NULL = 0, & ! No interaction
    COVB = 1, & ! Actual Covalent bond
    SSEC = 2, &  ! Secondary structure restrain
    FAKE  = 3 ! Just steric term
!    COUL = 4, & ! Electrostatic
!    HFIL = 5, & ! Solvation Hydrophilic
!    HFOB = 6, & ! Solvation Hydrophobic
!    COVF = 7, & ! Covalent restrain forced
!    DIHH = 8    ! NB H-H Distance
    integer, parameter :: MAXTIPINT = 3
    ! energy
    real, parameter :: FACTE = 4184. ! Energy conversion
  !  real, parameter :: VDWBAR = 10. ! VDW Wall for automatic step potentials
    !!distances
  !  real, parameter :: RSSMAX = 2.5
  !  real, parameter :: RNOMAX = 4.1
  !  real, parameter :: RNOMIN = 2.5
  !  real, parameter :: RNCMAX = 5.
  !  real, parameter :: RNCMIN = 3.2
  !  real, parameter :: RCOMAX = 5.
  !  real, parameter :: RCOMIN = 3.2
   ! REAL, PARAMETER :: SIGMAGO = 0.05
  !  real, parameter :: RCUTNB2 = 30. * 30. !!!! CUIDADO SI SE USA PARA AGREGACION
  !  real, parameter :: RCUTGO2 = 20. * 20.
    !structure
    integer, parameter :: PROT = 1, NUC = 2, SMALL = 3, COMPLEX = 4 ! molType
 !   integer, parameter :: ALL = 0, HEAVY = 1, CAONLY = 2 ! CGTYpe 
 !   integer, parameter :: HELIX = 1, BETA = 2

  !  integer, parameter :: MD = 0, DOCKING = 1 ! tipCalc
  !  integer, parameter :: itopVersion = 026

    !    character(len = 10), parameter :: topVersion = 'v0.2.6'
    !
    character(len = 50), parameter :: pdbinputFmt = '(13X,A3,1X,A4,A1,I4,4X,3F8.3,2f6.2)'
    real, parameter :: MINREAL = 1.e-20
    real*8, parameter :: A2 = 1.d-20*1.d30   ! cambio de unidades?

CONTAINS

END MODULE constants
