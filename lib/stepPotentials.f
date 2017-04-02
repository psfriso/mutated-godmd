
MODULE stepPotentials
    use Constants
    
    type stepPot
        real r
        real e
        logical initial
    end type stepPot

    type stepPotInt
        integer nstep, tipInt
        type (stepPot) step(MAXSTEPS)
        logical active
        real dego
        real dmin
        integer saltos
    end type stepPotInt

    type stepPotIntDef
        character(len = 4) :: id
        integer nstep, stepTip
        real r(MAXSTEPS), e(MAXSTEPS)
    end type stepPotIntDef
    
    type intDataTop
        integer :: id, tipInt
        real :: rref, eref
    end type intDataTop    
    
    type stepPotIntDefList
        integer npots
        type(stepPotIntDef) list(MAXSPOTDEF)
    end type stepPotIntDefList
    
  
    
CONTAINS
    !===============================================================================
    function getStepPotFromDEF(GOENERGY,sigmago,dtol,wellwidth,rrefINI ,rrefTARG,tipINT) result (st)
        !
       ! type(stepPotIntDef), intent(IN) :: stDef
        REAL, INTENT(IN):: GOENERGY,sigmago,wellWidth
        INTEGER, INTENT(IN) :: tipINT
        REAL*8, INTENT (IN) :: dtol
        type(stepPotInt) :: st
        real*8, intent(IN) :: rrefINI, rrefTARG
        !
        logical overlaping
        real*8 :: dmin, dmax
        
        st=stepPotInt(0,0,(/stepPot(0.00000,0.00000,.FALSE.),stepPot(0.00000,0.00000,.FALSE.),&
        stepPot(0.00000,0.00000,.FALSE.),stepPot(0.00000,0.00000,.FALSE.)/),.FALSE.,0.0000,0.00000,0)
        
        ! armstrong
        !
        ! Buscamos cual es la distancia mas pequeña y la mas grande
        dmin = MIN( rrefINI, rrefTARG )
        dmax = MAX( rrefINI, rrefTARG )
        overlaping=.FALSE.
        !
        ! El posible solapamiento es con dmin(1+sigmago) y dmax( 1- sigmago)
        !
        IF ( ABS(dmin - dmax) < dtol) overlaping = .TRUE.
        !
!
      IF (  tipINT == COVB) THEN
          !do covb
          st % nstep = 0
          st % step(1) = stepPot( rrefINI - 0.02 , -0.75*GOENERGY * FACTE,.TRUE.)
          st % step(2) = stepPot( rrefINI + 0.02 ,  0.75*GOENERGY * FACTE,.TRUE.)
          !write(*,*) "COVB"
      ELSE
        IF ( overlaping ) THEN
            ! Pozo con un solo minimo
            st % nstep = 2
            st % step(1) = stepPot( dmin - wellWidth, -0.75*GOENERGY * FACTE,.TRUE.)
            st % step(2) = stepPot( dmax + wellWidth,  0.75*GOENERGY * FACTE,.TRUE.)
        ELSE
            ! Pozo con 4 steps, doble minimo
            st % nstep = 4
            IF (dmin == rrefINI) THEN
                !El primer pozo corresponde a la estructura inicial
            ! El primer pozo
            st % step(1) = stepPot( dmin - sigmago, -GOENERGY * FACTE,.TRUE.)
            st % step(2) = stepPot( dmin + sigmago,  GOENERGY * FACTE,.TRUE.)
            ! El segundo
            st % step(3) = stepPot( dmax - sigmago, -GOENERGY * FACTE,.FALSE.)
            st % step(4) = stepPot( dmax + sigmago,  GOENERGY * FACTE,.FALSE.)

            ELSE
                ! El segundo pozo corresponde a la estructura inicial
                ! El primer pozo
            st % step(1) = stepPot( dmin - sigmago, -GOENERGY * FACTE,.FALSE.)
            st % step(2) = stepPot( dmin + sigmago,  GOENERGY * FACTE,.FALSE.)
            ! El segundo
            st % step(3) = stepPot( dmax - sigmago, -GOENERGY * FACTE,.TRUE.)
            st % step(4) = stepPot( dmax + sigmago,  GOENERGY * FACTE,.TRUE.)
      !      Write(*,'(A4,10f8.4,X,3L2)')"stp ",st % step(1)%r,st % step(2)%r,st % step(3)%r,st % step(4)%r,&
      !      dtol,sigmago,rrefINI,rrefTARG,dmin,dmax,overlaping,st % step(1) %initial,st % step(3) %initial
      !      stop
            ENDIF
            !Write(*,'(A4,10f8.4,X,L1)')"stp ",st % step(1)%r,st % step(2)%r,st % step(3)%r,st % step(4)%r,&
            !dtol,sigmago,rrefINI,rrefTARG,dmin,dmax,overlaping
            !STOP
        ENDIF
      END IF
        st % tipInt = tipINT
        st % active = .FALSE.
        st % dego = 0.0
        st % dmin = 0.
        st % saltos = 0
!        call writeStPotDef(6,stDef)
!        call writeStPot(6,st)
    end function getStepPotFromDef
  !===============================================================================

      function getStepPotForDummy(GOENERGY,sigmaFake,rrefINI,tipInt) result (st)
        !
       ! type(stepPotIntDef), intent(IN) :: stDef
        REAL, INTENT(IN):: GOENERGY
        INTEGER, INTENT(IN) :: tipINT
        type(stepPotInt) :: st
        REAL*8, INTENT(IN):: rrefINI
        Real, Intent(IN) :: sigmaFake
        
        !
        
        st=stepPotInt(0,0,(/stepPot(0.00000,0.00000,.FALSE.),stepPot(0.00000,0.00000,.FALSE.),&
        stepPot(0.00000,0.00000,.FALSE.),stepPot(0.00000,0.00000,.FALSE.)/),.FALSE.,0.0000,0.00000,0)
        
            ! Pozo con un solo minimo
       st % nstep = 2
       st % step(1) = stepPot( rrefINI*(1.d0-DBLE(sigmaFake)), -0.75*GOENERGY * FACTE,.TRUE.)
       st % step(2) = stepPot( rrefINI*(1.d0+DBLE(sigmaFake)), 0.75*GOENERGY * FACTE,.TRUE.)
     
        st % tipInt = tipINT
        st % active = .FALSE.
        st % dego = 0.0
        st % dmin = 0.
        st % saltos = 0
!        call writeStPotDef(6,stDef)
!        call writeStPot(6,st)
    end function getStepPotForDummy
  !===============================================================================
 
    END MODULE stepPotentials
