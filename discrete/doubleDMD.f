!     
! File:   doubleDMD.f
! Author: psfriso
!
! Created on October 25, 2012, 6:17 PM
!
MODULE doubleDMD

    PRIVATE
    
    PUBLIC metaDMD,AreInB,updateXbeta,decideRCut4GO,fixDistances,modulateDEgo
 
  INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(p=13)! Double precision real number definition processor independent
  INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(p=6) ! Single Precision real number definition processor independent

    CONTAINS
 !==================================================================================   
    SUBROUTINE metaDMD(natom, r, stepPts,discardedEner)
        use stepPotentials
        use geometryDP
        !
        IMPLICIT NONE
        !
        type(pointDP), dimension(natom),INTENT(IN) :: r
        integer , intent (IN) :: natom
        type(stepPotInt),INTENT(INOUT),DIMENSION(natom,natom):: stepPts
        real,dimension(natom,natom),intent(INOUT):: discardedEner
        !
        integer i, j
        real :: newEner=0.0
      !  character(len=53) :: fmtsteppot='(2I4,I2,x,I1,4(2(f7.1,X),L1,X),L1,X,f9.4,f8.3,I2)'
        
    DO j=2,natom
        DO i =1,j-1
         
          IF(ThrowEner( stepPts(i,j),calcDistDP(r(i),r(j)) )) THEN
              
                 stepPts(i,j) % saltos = stepPts(i,j) % saltos-1
       
             IF (stepPts(i,j)%step(1)%initial) THEN
                 !
                 newEner=abs(stepPts(i,j)%step(1)%e) - stepPts(i,j)%dego * FACTE
                 
                 discardedEner(i,j)=discardedEner(i,j)+ stepPts(i,j)%dego * FACTE
                 ! ojo con el signo
                 stepPts(i,j)%step(1)%e= sign( newEner, stepPts(i,j)%step(1)%e)
                 stepPts(i,j)%step(2)%e= sign( newEner, stepPts(i,j)%step(2)%e)
                 
             ELSEIF (stepPts(i,j)%step(3)%initial ) THEN
                 !
                 newEner=abs(stepPts(i,j)%step(3)%e) - stepPts(i,j)%dego * FACTE
                 discardedEner(i,j)=discardedEner(i,j)+ stepPts(i,j)%dego * FACTE
                 ! ojo con el signo
                 stepPts(i,j)%step(3)%e= sign( newEner, stepPts(i,j)%step(3)%e)
                 stepPts(i,j)%step(4)%e= sign( newEner, stepPts(i,j)%step(4)%e)
             ENDIF
   !           
              stepPts(j,i)= stepPts(i,j)
            ENDIF
        END DO
    END DO
    
    END SUBROUTINE metaDMD 
    
    !===============================================================================
 logical function ThrowEner(stp,distIJ) result (enINITIAL)
     use stepPotentials
     !
     IMPLICIT NONE
     !
        type (stepPotInt), INTENT(IN) :: stp
        real*8, intent(IN) :: distij
        !
        real :: aux
        !
        enINITIAL=.FALSE.
        !
        aux=real(distIJ)
        !
        IF(  stp % saltos >0) THEN
        IF ( stp %step(1) %initial ) THEN
            enINITIAL=(stp %step(1) %r < aux .and. stp %step(2) %r > aux)
        ELSEIF ( stp %step(3) %initial)THEN
            enINITIAL=(stp %step(3) %r < aux .and. stp %step(4) %r > aux)
        ENDIF
        ELSE
            enINITIAL=.FALSE.
        ENDIF
!    
    end function ThrowEner
  !===============================================================================
  !===============================================================================
 function InFinalWell(stp,distIJ) result (enB)
     use stepPotentials
     !
     IMPLICIT NONE
     !
        type (stepPotInt), INTENT(IN) :: stp
        real, intent(IN) :: distij
        logical :: enB
        !
        enB=.FALSE.
        IF ( .not.stp %step(1) %initial ) THEN
            enB=(stp %step(1) %r < distIJ .and. stp %step(2) %r > distIJ)
        ELSEIF (.not.stp %step(3) %initial)THEN
            enB=(stp %step(3) %r < distIJ .and. stp %step(4) %r > distIJ)
        ENDIF
!    
 end function InFinalWell
  !===============================================================================
  !
    SUBROUTINE AreInB(natom,stepPts,r,numEnB)
        use geometryDP
        use stepPotentials
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: natom
        TYPE (stepPotInt),INTENT(IN),DIMENSION(natom,natom) :: StepPts
        TYPE (pointDP), INTENT(IN), DIMENSION (natom):: r
        INTEGER , INTENT(OUT) :: numEnB
        !
        integer :: i,j
!
        numEnB=0
   !     aux=0
    do j = 2, natom-1
        do i = 1, j - 1
            
          IF (StepPts(i,j) % nstep == 4 .and. StepPts(i,j) % active) then

               IF ( InFinalWell( stepPts(i,j), real(calcDistDP(r(i),r(j))))) numEnB=numEnB+1
          END IF

        END DO
    END DO
 !   
    END SUBROUTINE AreInB
!
!====================================================================================    
  function updateXbeta (xbeta,aux_acc,acceptance) result (xbeta_updated)
      !
      IMPLICIT NONE
      !
            real, intent(IN) :: xbeta
            real, intent (IN) :: aux_acc,acceptance
            !
            real :: xbeta_updated
            real :: a_xbeta, b_xbeta, c_xbeta , d_xbeta
            !
            real :: factor
            !
            ! pol(x)= a x**3 + b x**2 + c x +d , x es la aceptancion obtenida
            !
            d_xbeta=-0.27 * acceptance**2 +3.0 * acceptance**3
            c_xbeta= 0.54 * acceptance -9.0 * acceptance**2
            b_xbeta= -0.27  +9.0 * acceptance
            a_xbeta=-3.0
            factor=0.0
            xbeta_updated=0.0
            !
            factor= a_xbeta * aux_acc**3 + b_xbeta * aux_acc**2 + c_xbeta * aux_acc + d_xbeta
            !
            xbeta_updated=xbeta *(1+factor)
 
  end function updateXbeta
  !====================================================================================
    SUBROUTINE decideRCut4GO(natom,r,rtarg,xm,rgyrPercent,mxRcut4go,rgyr,rcut2GO4)
        use geometryDP
        !
        IMPLICIT NONE
        !
        INTEGER, INTENT(IN) :: natom
        TYPE(pointDP), INTENT (IN) :: r(natom),rtarg(natom)
        REAL, INTENT (IN) :: xm
        REAL, INTENT(IN) :: rgyrPercent
        REAL*8, INTENT(IN) :: mxRcut4go
        !
        REAL*8, INTENT(OUT) :: rgyr
        REAL*8, INTENT (INOUT) :: rcut2GO4
        !
        real :: rgyrINI, rgyrTARG
        type(pointDP) :: rcm, distanceRCM
        integer :: i
        
    !
    ! Radius de gyration computation ( o parecido)
    !
    ! Para la estructura inicial
    rcm = SPtoDP(calcCM(natom, r, xm))
    rgyrINI=0.0
    rgyrTARG=0.0
    distanceRCM=pointDP(0,0,0)
    
    Do i=1,natom
        distanceRCM=r(i)-rcm
        rgyrINI=rgyrINI+moduleDP(distanceRCM)
    END DO
   
   rcm = SPtoDP(calcCM(natom, rtarg, xm)) 
   distanceRCM=pointDP(0,0,0)
   
    Do i=1,natom
        distanceRCM=rtarg(i)-rcm
        rgyrTARG =rgyrTARG+moduleDP(distanceRCM)
    END DO
    
    
    rgyr=MAX(rgyrINI,rgyrTARG)/DBLE(natom)
    rcut2GO4=MIN((rgyrPercent*rgyr)**2,DBLE(mxRcut4go*mxRcut4go))
    
    END SUBROUTINE decideRCut4GO
 !======================================================================================   
!===============================================================================================    
    SUBROUTINE fixDistances(natom,natomCommon,oid,dist_INI,dist_TARG,ssectol,stepPts)
        use stepPotentials
        !
        IMPLICIT NONE
        ! Evaluar las distancias i -- i+4, i --i+5 para poner un enlace COVB si no cambian 
        !
        INTEGER , INTENT(IN) :: natom,natomCommon
        REAL*8 , INTENT(IN) ::dist_INI(natomCommon,natomCommon),dist_TARG(natomCommon,natomCommon)
        TYPE (stepPotInt),INTENT (INOUT) :: stepPts(natom,natom)
        REAL*8, INTENT(IN) :: ssectol
        INTEGER, DIMENSION(natomCommon),INTENT(IN):: oid
        !
        integer i, j
        
iplus2loop: DO i=1,natomCommon-4
  iplus2IF: IF ( ABS(dist_INI(i,i+2)-dist_TARG(i,i+2))/MIN(dist_INI(i,i+2),dist_TARG(i,i+2)) < ssectol) THEN
  iplus3IF:  IF (ABS(dist_INI(i,i+3)-dist_TARG(i,i+3))/MIN(dist_INI(i,i+3),dist_TARG(i,i+3)) < ssectol ) THEN
  iplus4IF:  IF (ABS(dist_INI(i,i+4)-dist_TARG(i,i+4))/MIN(dist_INI(i,i+4),dist_TARG(i,i+4)) < ssectol) THEN
                 j=i+4
                 IF ( stepPts(oid(i+1), oid(i) ) % tipInt == COVB) THEN
                ! Ajustar el i+1 
                  stepPts(oid(i+1) , oid(i) ) % step(1) %r = dist_INI(i,i+1) *DBLE(0.98)
                  stepPts(oid(i), oid(i+1) ) % step(1) %r = dist_INI(i,i+1) *DBLE(0.98)
                  stepPts(oid(i+1), oid(i) ) % step(2) %r = dist_INI(i,i+1) *DBLE(1.02)
                  stepPts(oid(i), oid(i+1) ) % step(2) %r = dist_INI(i,1+1) *DBLE(1.02) 
          !      Ajustar el i+4
                 stepPts(oid(j), oid(i) ) % tipInt = COVB
                 stepPts(oid(i), oid(j) ) % tipInt = COVB
                 stepPts(oid(j), oid(i) ) % step(1) %r = dist_INI(i,j) *DBLE(0.97)
                 stepPts(oid(i), oid(j) ) % step(1) %r = dist_INI(i,j) *DBLE(0.97)
                 stepPts(oid(j), oid(i) ) % step(2) %r = dist_INI(i,j) *DBLE(1.03)
                 stepPts(oid(i), oid(j) ) % step(2) %r = dist_INI(i,j) *DBLE(1.03)
                 stepPts(oid(i), oid(j) )  % nstep = 0
                 stepPts(oid(j), oid(i) )  % nstep = 0
                 END IF
                END IF iplus4IF
              END IF iplus3IF
            END IF iplus2IF
        END DO iplus2loop

  END SUBROUTINE fixDistances
 !==============================================================================================  
!======================================================================================  
    SUBROUTINE modulateDEgo(unit_o,natom,natomCommon,oid,evecs,evals,r,dist_INI,stepPts) 
       use paramSet
       use stepPotentials
       use geometryDP
       !
       IMPLICIT NONE
   
       INTEGER, INTENT(IN) :: unit_o,natom,natomCommon
       TYPE (pointDP), INTENT(IN), DIMENSION(nevecs,natom) :: evecs
       REAL, INTENT(IN), DIMENSION (nevecs) :: evals
       TYPE (pointDP), INTENT(IN) :: r(natom)
       REAL *8 , INTENT(IN) :: dist_INI(natomCommon,natomCommon)
       TYPE(stepPotInt), INTENT(INOUT), DIMENSION(natom,natom) :: stepPts
       INTEGER, DIMENSION(natomCommon),INTENT(IN):: oid
  
  !
  TYPE(pointDP), ALLOCATABLE, DIMENSION(:) :: nuevaCoord
  REAl, ALLOCATABLE, DIMENSION(:,:,:) :: cambioRelDist
  REAL*8 :: distanciaDespues
  REAL :: sumaEvals
  REAL, Allocatable, Dimension(:)::factorEval
  REAL,DIMENSION(nevecs) :: auxMax
  integer :: nSaltosOrig=0
  real :: originalRatio=0
  real :: factorDego
  real :: maxTemporal
  integer :: maxTrack
  integer :: i,j,k
  !
  ALLOCATE(factorEval(nevecs),cambioRelDist(nevecs,natom,natom),nuevaCoord(natom))
   
        
       
   sumaEvals=0.0
   auxMax=0.0
      !
    ! Vamos a evaluar como cambian las distancias entre los pozos cuando seguimos (doble pozo)
    ! un determinado eigenvector un terminado hstep, fijo h=0.1 
    ! si las distancias se acortan, exactamente 0.1 entonces el overlap del
    ! movimiento con el NM sera 100 %
 
   DO i=1,nevecs
   sumaEvals=sumaEvals+1.0/evals(i)
   END DO
   
   DO i=1,nevecs
       factorEval(i)= 1.0/evals(i)/sumaEvals
   WRITE (unit_o, '("PORCENTAJE ", i2, f8.4)')i, 1.0/evals(i)/sumaEvals
   END DO
   
   cambioRelDist=0.0000000000

   
   DO k=1,nevecs
       
     nuevaCoord=pointDP(0._DBL,0._DBL, 0._DBL) 
       
    DO i=1,natom
        
        nuevaCoord(i) = r(i) + hstep *evecs(k,i)
    
    END DO
    
    
    ! Nuevas distancias internas

    distanciaDespues=0.
 !   
    do j=2,natomCommon
  !
        do i = 1, j - 1
            distanciaDespues=abs(dist_INI(i,j) - calcDistDP(nuevaCoord(oid(i)), nuevaCoord(oid(j) )))
            cambioRelDist(k,i,j)= distanciaDespues/hstep 
            cambioRelDist(k,j,i) =  cambioRelDist(k,i,j)

     !       if(i<10 .and. j<10) write(unit_o,*) cambioRelDist(k,i,j)
       end do
    end do
    
    auxMax(k)=MAXVAL(cambioRelDist(k,:,:))
   END DO
 
  ! Ahora ya tenemos el tensor cambioRelDist( eigenvector , i , j) 
  ! vamos a normalizar las distancias respecto del maximo desplazamiento por eigenvector
  ! asi podemos tener una idea mas clara de que movimientos son importantes
  !
  !  Ajustar DEGO
   !
   ! Vamos a  buscar en cuales de los nevecs una dada distancia i-j sufre un desplazamiento 
   ! mayor
   !
     Do j=2,natom
         do i= 1, j-1
           maxTemporal=0.0         
             DO k=1,nevecs
                If( maxTemporal < cambioRelDist(k,i,j) ) THEN
                    maxTemporal =cambioRelDist(k,i,j)
                    maxTrack = k
                END IF
             END DO
            factorDego=DEGO_BASE*&
            (1._SGL+scalFactor*factorEval(maxTrack)*cambioRelDist(maxTrack,i,j)/auxMax(maxTrack))
            
            stepPts(i,j) % dego =factorDego* GOENER
             
         END DO
    END DO

   
   ! decidir cuantos saltos de energia puede hacer
    !
     do j = 2, natom
        do i = 1, j-1
  !          
     IsDoubleWell: IF ( stepPts(i,j) % nstep == 4 ) THEN
                
      maxEner:     IF (stepPts(i, j) % dego > 0.95*(GOENER-minEnerWell))  THEN
  
                    stepPts(i, j )  % saltos =1
                    stepPts(i, j ) % dego = GOENER-minEnerWell
  
      ELSE maxEner
                    
                  originalRatio=0._SGL
                  originalRatio=(GOENER-minEnerWell)/stepPts(i, j) % dego
 
                  nSaltosOrig=FLOOR(originalRatio)
                  
          checkSaltos:      IF (nSaltosOrig < 1 ) THEN
                      
                      Write(0,*) " Numero de saltos menor que 1 encontrado. Error"
                      STOP
                      
                  ELSE checkSaltos
                 
                   stepPts( i , j) % dego =  REAL (nSaltosOrig) *stepPts(i,j) % dego /real(originalRatio)
                   stepPts( i , j) % saltos=  nSaltosOrig+1
                   
                   END IF checkSaltos
                 
                ENDIF maxEner

                ELSEIF( stepPts(i,j) % nstep == 2) THEN !ISDoubleWell
                
                 stepPts(i, j ) % saltos =0
                 stepPts(i, j ) % dego=0.0
                 
            END IF IsDoubleWell
            stepPts(j, i )= stepPts(i,j)
            
        END DO
     END DO
            
   END SUBROUTINE modulateDEgo
 !============================================================================
             
END MODULE doubleDMD
