!     
! File:   ANM.f
! Author: psfriso
!
! Created on October 30, 2012, 12:45 PM
!

MODULE ANM
! Data dictionary and variable declaration
INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(p=13)! Double precision real number definition processor independent
INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(p=6) ! Single Precision real number definition processor independent

PRIVATE
PUBLIC ANM1
 CONTAINS
 !
       SUBROUTINE ANM1(nres,coord,nevecs, evec,evals)
       USE geometryDP
       USE geometry
       USE ls_rmsd

       IMPLICIT NONE  
       INTEGER ::RES3,i,k,J,ioerr
       INTEGER :: CUTOFF
       INTEGER, ALLOCATABLE :: INDX(:)
       REAL(SGL), DIMENSION(nevecs), INTENT(INOUT) :: evals
       REAL(DBL), ALLOCATABLE :: DIS(:,:),G(:,:),HESS(:,:)
       INTEGER ,intent(IN):: NRES ,nevecs  ! Residue number
       REAL(DBL), ALLOCATABLE :: W(:),V(:,:),WD(:),INVCONT(:,:)
       REAL(SGL), PARAMETER :: KBT=0.597,EIGENCUT=1E-5
       REAL(DBL), ALLOCATABLE :: rv1(:)
       REAL, PARAMETER :: SL=3, IN=2.8
       TYPE(pointDP), DIMENSION(nres), INTENT(IN) :: coord
       TYPE(pointDP), INTENT(INOUT), DIMENSION(nevecs,nres) :: evec
       REAL(DBL), dimension(nres, 3) :: r2
        r2(:,1)  =coord(:)%x
        r2(:,2)  =coord(:)%y
        r2(:,3)  =coord(:)%z
      RES3=NRES*3
!      WRITE (6,*)'Allocating arrays..'
      ALLOCATE (V(RES3,RES3) ,W(RES3), WD(RES3), INDX(RES3), INVCONT(NRES,NRES)&
      ,DIS(NRES,NRES),G(NRES,NRES),HESS(RES3,RES3),rv1(RES3),STAT=ioerr)
       call errorAllocmem (ioerr, "ANM ")
        V=0.0_DBL
        W=0.0_DBL
        WD=0.0_DBL
        INDX=0
        INVCONT=0.0_DBL
! Setting cutoff
            IF (NRES.LT.50) THEN
              CUTOFF=8
            ELSE
              CUTOFF=INT(SL*LOG(REAL(NRES))-IN)
            ENDIF
!
!        WRITE (*,*) 'Computing C-alpha distances...'

        CALL DISTANCES (NRES,R2,DIS)

!        WRITE (*,*) 'Setting force constants...'

        CALL SET_GAMMA(CUTOFF,NRES,DIS,G)
!
!        WRITE (*,*) 'Now we calculate hessian...'

        CALL HESSIAN(NRES,R2,DIS,G,HESS)
!
!        write(*,*) 'Diagonalizing....'
        CALL SVDCMP(RES3, HESS, RES3, RES3, W, V, rv1)
!        write(*,*) 'Putting the evals in order....'
        CALL INDEXX(RES3,W,INDX) 

   DO I=1,nevecs+6
           k=1
           IF(I.GT.6) THEN
              DO j=1,3*nres,3
!              write(*,*) "desde aqui" ,V(J, INDX(I))
              evec(I-6,k)%x=V(j,INDX(I))
              evec(I-6,k)%y=V(j+1,INDX(I))
              evec(I-6,k)%z=V(j+2,INDX(I))
              k=k+1
!              
              evals(I-6)=W(INDX(i))
!              
              END DO
           END IF
      ENDDO

  !    DO i=6,11
  !      WRITE(*,*)" INDX ", i,  W(INDX(i))
  !    END DO
      
!     WRITE(*,*) 'ANM Subroutine finished !!'
     DEALLOCATE(V ,W, WD, INDX, INVCONT,DIS,G,HESS,rv1)
       END SUBROUTINE ANM1


!********************************************************************************************
!*       SUBROUTINES TO SET GAMMA NRES x NRES ARRAY ACCORDING TO MODE...
!********************************************************************************************
        SUBROUTINE SET_GAMMA(CUTOFF,NRES,DIS,G)
        INTEGER NRES,S,CUTOFF
        REAL*4, PARAMETER :: Clin=10_SGL,Ckov=40_SGL,Ccart=6_SGL,Cseq=60_SGL
        INTEGER, PARAMETER :: Slim=3,Ex=6
        INTEGER:: K,J
        REAL(DBL), INTENT(INOUT) :: G(NRES,NRES)
        REAL(DBL),INTENT(IN) :: DIS(NRES,NRES)
!        WRITE (*,*)'Allocating gamma array..'
        G=0.0_DBL
!        HESS=0.0
!        WRITE (*,*) 'Computing with default values...'
 !Here we build g(j,k) terms that are a combination between gamma parameter and kirchhoff elements(-1/0) 
 DO J=1,nres
  DO K=1,nres
       IF(J.NE.K) THEN
         S=ABS(J-K)
         IF(S.LE.Slim) THEN
             G(J,K)=-(Cseq/(S**2))
         ELSE
             IF(DIS(J,K).LE.CUTOFF) THEN
                 G(J,K)=DBLE(-((Ccart/DIS(J,K))**Ex))
             ELSE
                 G(J,K)=0_DBL
            ENDIF
        ENDIF
       ENDIF
     ENDDO
   ENDDO
         END SUBROUTINE SET_GAMMA
!********************************************************************************************
!*       SUBROUTINES TO SET HESSIAN 3*NRES x 3*NRES ARRAY OF 2ND DERIVATIVES OF POTENTIAL
!********************************************************************************************
        SUBROUTINE HESSIAN (NRES,R,DIS,G,HESS)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NRES
        INTEGER RES3
        REAL(DBL) :: BX,BY,BZ,DIS2,GAMMA
        REAL(DBL), INTENT(IN) :: R(NRES,NRES)
        REAL(DBL), INTENT(IN) :: DIS(NRES,NRES),G(NRES,NRES)
        REAL(DBL), INTENT(INOUT) :: HESS(3*NRES,3*NRES)
        INTEGER :: J,K
        RES3=3*NRES 
!        WRITE (*,*)'Allocating hessian array, dimension', RES3
        HESS=0.0
!        write(*,*) 'Calculating the hessian..'
    DO J=1,NRES
     DO K=1,NRES
           bx=0.0
           by=0.0
           bz=0.0
           dis2=0.0
           gamma=0.0
          IF(J.NE.K)then
           BX=(r(J,1)-r(K,1))
           BY=(r(J,2)-r(K,2))
           BZ=(r(J,3)-r(K,3))
           DIS2=DIS(J,K)**2
           Gamma=G(J,K)
           HESS(3*J-2,3*K-2)=+GAMMA*BX*BX/DIS2
           HESS(3*J-2,3*K-1)=+GAMMA*BX*BY/DIS2
           HESS(3*J-2,3*K)=  +GAMMA*BX*BZ/DIS2
           HESS(3*J-1,3*K-2)=+GAMMA*BY*BX/DIS2
           HESS(3*J-1,3*K-1)=+GAMMA*BY*BY/DIS2
           HESS(3*J-1,3*K)=  +GAMMA*BY*BZ/DIS2
           HESS(3*J,3*K-2)=+GAMMA*BZ*BX/DIS2
           HESS(3*J,3*K-1)=+GAMMA*BZ*BY/DIS2
           HESS(3*J,3*K)=  +GAMMA*BZ*BZ/DIS2
!C	  SECOND: CREATION OF DIAGONAL Hii(Hii=-SUM(Hij/j#i) TERMS
           HESS(3*J-2,3*J-2)=HESS(3*J-2,3*J-2)-GAMMA*BX*BX/DIS2
           HESS(3*J-2,3*J-1)=HESS(3*J-2,3*J-1)-GAMMA*BX*BY/DIS2
           HESS(3*J-2,3*J)=HESS(3*J-2,3*J)  -GAMMA*BX*BZ/DIS2
           HESS(3*J-1,3*J-2)=HESS(3*J-1,3*J-2)-GAMMA*BY*BX/DIS2
           HESS(3*J-1,3*J-1)=HESS(3*J-1,3*J-1)-GAMMA*BY*BY/DIS2
           HESS(3*J-1,3*J)=HESS(3*J-1,3*J)  -GAMMA*BY*BZ/DIS2
           HESS(3*J,3*J-2)=HESS(3*J,3*J-2) -GAMMA*BZ*BX/DIS2
           HESS(3*J,3*J-1)=HESS(3*J,3*J-1) -GAMMA*BZ*BY/DIS2
           HESS(3*J,3*J)=HESS(3*J,3*J)  -GAMMA*BZ*BZ/DIS2
           ENDIF
          ENDDO
         ENDDO 
!         WRITE (*,*) 'Hessian done!!'
!         RETURN
         END SUBROUTINE HESSIAN

!*******************************************************************************
!* SUBROUTINE DISTANCES - COMPUTES NRES x NRES ARRAY OF CALPHA-CALPHA DISTANCES
!*******************************************************************************
SUBROUTINE DISTANCES (NRES,R,DIS)
      IMPLICIT NONE
      INTEGER :: NRES,I,J
      REAL(DBL) , INTENT(IN) :: R(NRES,3)
      REAL(DBL), INTENT(INOUT), DIMENSION(NRES,NRES) :: Dis
      REAL(DBL) :: x,y,z
!      WRITE(*,*) 'Now we compute c-ALPHA distances...'
     DO  I=1,NRES
        DO  J=1,NRES
        if(j.ne.i) then
          X=(R(I,1)-R(J,1))**2
          Y=(R(I,2)-R(J,2))**2
          Z=(R(I,3)-R(J,3))**2
          DIS(I,J)=DSQRT(X+Y+Z)
          endif
        ENDDO
      ENDDO
    END SUBROUTINE DISTANCES

!*******************************************************************************
!* SUBROUTINE INDEXX - ORDERS by eigenvalues
!*******************************************************************************
      SUBROUTINE INDEXX(N,ARRIN,INDX)
      INTEGER, INTENT(IN) :: N
      INTEGER :: J,L,IR,INDXT,I
      REAL(DBL), INTENT(IN), DIMENSION(N) :: ARRIN
      INTEGER, INTENT(INOUT) :: INDX(N)
      REAL(DBL) :: q
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      IF(N.EQ.1)RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10 
      END SUBROUTINE INDEXX
END MODULE ANM

