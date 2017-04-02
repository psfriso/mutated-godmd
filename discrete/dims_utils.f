!     
! File:   dims.f
! Author: psfriso
!
! Created on September 14, 2012, 11:35 AM
!

MODULE dims_utils
!
! Purpose:
! Set of tools regarding energy and sq-wells potential interactions.
!
! Record of revisions:
!      DATE            PROGRAMMER           DESCRIPTION OF CHANGE
!    =======          ============         =======================
!    14/09/11         Pedro                  Original Code
!    14/09/11         Pedro                  Defined as module. Comments added.

IMPLICIT NONE
!
PRIVATE
!

PUBLIC distanceCA, MCWeigth, MCCheck,ForceIT

!
! Data dictionary and variable declaration
INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(p=13)  ! Double precision real number definition processor independent
INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(p=6)   ! Single Precision real number definition processor independent
REAL, PARAMETER :: PI=3.141592654 ! PI number definition
!
CONTAINS

!===============================================================================
 function distanceCA(natom, oid,r,rtarg,w,natomCommon) result (score)
     ! distanceCA(natom, oid,r,rtarg,w,natomCommon)
 use geometryDP
   integer, intent(IN) :: natom,natomCommon
   integer, intent(in), dimension(natomCommon) :: oid
   type(pointDP), intent(IN) :: r(natom),rtarg(natomCommon)
   real, intent(in) :: w(natomCommon)
   real score, rij
   integer :: i
   score=0.
  do i=1,natomCommon
   rij = calcDistDP(r(oid(i)),rtarg(i))
   score=score+(rij)*w( i )
  enddo
 ! write(*,*) "SCORE ", score
 end function distanceCA
!===============================================================================

 function MCCheck (error,seed, xbeta, scoreprev, score) result (rej)
 ! If TRUE it makes trajectory rewind
     use ran_mod  
   integer, intent(INOUT) :: seed
   real, intent(IN) :: xbeta
   real(DBL), intent(in) :: score,error,scoreprev
   logical rej
   real sto,dscore
   real*8 :: fi

     dscore=scoreprev-score
   
    fi=ran1(seed)
     
    sto=exp(1.0/(xbeta)**2 * sign(1.,dscore)*(dscore)**2 /error**2)

    rej = (sto.lt.fi)

 end function MCCheck
!===============================================================================
SUBROUTINE MCWeigth(natom,oid,r, rtarg,w,natomCommon)
   !       MCWeigth(natom,oid,r, rtarg,w,natomCommon)
use geometryDP
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: natom,natomCommon
INTEGER,DIMENSION(natomCommon),INTENT(in)::oid
TYPE(pointDP),DIMENSION(natom), INTENT(IN) :: r
TYPE(pointDP),DIMENSION(natomCommon), INTENT(IN) :: rtarg
REAL(SGL) , INTENT(INOUT):: w(natomCommon)
INTEGER :: i
!integer, INTENT(IN):: natom
REAL :: temp_distance


do i=1,natomCommon
    temp_distance = calcDistDP(rtarg(i) ,r(oid(i)) ) 
    w(i)=(1.-switch_w(natom))*&
    (5000.0*Log(temp_Distance+1)/(2.5 *(20.0 - temp_Distance)**2 +500.0))+&
    (switch_w(natom))*&
    (5000.0*Log(temp_Distance+1)/(1.2*(15.0 - temp_Distance)**2 +500.0))
enddo

END SUBROUTINE MCWeigth
!======================================================================================
!===============================================================================
SUBROUTINE ForceIT(natom,natomCommon,oid,r, rtarg,w)
use geometryDP
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: natom,natomCommon
INTEGER, DIMENSION(natomCommon), INTENT(IN) :: oid
TYPE(pointDP),DIMENSION(natom), INTENT(IN) :: r
TYPE(pointDP),DIMENSION(natomCommon), INTENT(IN) :: rtarg
REAL(SGL) , INTENT(INOUT):: w(natomCommon)
INTEGER :: i
!integer, INTENT(IN):: natom

 do i=1,natomCommon
    w( i )= calcDistDP(rtarg( i ) ,r( oid(i) )  )*10.00
 enddo

END SUBROUTINE ForceIT
!======================================================================================
SUBROUTINE makeunit_3N(N, vec_3N)
USE geometry
USE geometryDP
IMPLICIT NONE
INTEGER, INTENT(IN) :: N
TYPE(pointDP), INTENT(INOUT), DIMENSION(N) :: vec_3N
REAL(DBL) :: norma
INTEGER ::i
!
norma=0._DBL
DO i=1,N
  norma=norma+vec_3N(i)%x**2 +vec_3N(i)%y**2 + vec_3N(i)%z**2
ENDDO
norma=SQRT(norma)
DO i=1,N
  vec_3N(i)%x=vec_3N(i)%x/norma
  vec_3N(i)%y=vec_3N(i)%y/norma
  vec_3N(i)%z=vec_3N(i)%z/norma
END DO
END SUBROUTINE makeunit_3N
!==========================================================
FUNCTION switch_w(natom) result (lambda)
    IMPLICIT NONE
    real lambda
    integer , intent(in) :: natom
    integer :: natomCutoff = 600
    real :: smoothing = 0.4_SGL
    real :: exponential=0
    !
    exponential=Exp(-smoothing*(natom-natomCutoff))
    lambda= 1._SGL/(1._SGL + exponential)
    !
END FUNCTION switch_w
!==========================================================
END MODULE dims_utils