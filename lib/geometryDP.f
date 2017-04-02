!
! $Date: 2012-06-15 18:58:02 $
! $Id: geometryDP.f,v 1.1 2012-06-15 18:58:02 gelpi Exp $
! $Revision: 1.1 $
!
 MODULE geometryDP

 TYPE pointDP 
   REAL*8 :: x,y,z 
 END TYPE pointDP

 interface operator (+) 
  module procedure sumavecDP
 end interface

 interface operator (-)
  module procedure negvecDP, restavecDP
 end interface

 interface operator (*)
  module procedure prodFactDP, dotDP
 end interface

CONTAINS

!===============================================================================   
pure function calcCM(natom, r, xm) result (rcm)
    use geometry
    type(point) rcm
    integer, intent(IN) :: natom
    type(pointDP), intent(IN) :: r(natom)
    real, intent(IN) :: xm
    integer i
    real xmassa
    !
    rcm = point(0., 0., 0.)
    xmassa = natom*xm
    do i = 1, natom
        rcm % x = rcm % x + xm * real(r(i) % x)
        rcm % y = rcm % y + xm * real(r(i) % y)
        rcm % z = rcm % z + xm * real(r(i) % z)
    enddo
    rcm = (1./xmassa) * rcm
    end function calcCM
!=============================================================================== 
pure function SPtoDP (rsp) result (rdp)
use geometry
 type(point), intent(IN) :: rsp
 type(pointDP) rdp
 rdp%x=rsp%x
 rdp%y=rsp%y
 rdp%z=rsp%z
end function SPtoDP

pure function DPtoSP (rdp) result (rsp)
use geometry
 type(point) rsp
 type(pointDP), intent(IN) :: rdp
 rsp%x=real(rdp%x)
 rsp%y=real(rdp%y)
 rsp%z=real(rdp%z)
end function DPtoSP

pure function sumavecDP (v1,v2)
type (pointDP):: sumavecDP
TYPE (pointDP), intent (IN):: v1,v2
sumavecDP%x=v1%x+v2%x
sumavecDP%y=v1%y+v2%y
sumavecDP%z=v1%z+v2%z
end function sumavecDP

pure function prodFactDP (f,v)
type (pointDP) :: prodFactDP
TYPE (pointDP), intent (IN):: v
real*8, intent (IN):: f
prodFactDP%x=v%x*f
prodFactDP%y=v%y*f
prodFactDP%z=v%z*f
end function prodFactDP

pure function negvecDP (v)
type (pointDP) :: negvecDP
type (pointDP), intent (in) :: v
negvecDP = prodFactDP (-1.d0, v)
end function negvecDP

pure function restaVecDP (v1,v2)
type (pointDP) :: restaVecDP
type (pointDP), intent(IN) :: v1, v2
restaVecDP = v1 + (-v2)
end function restaVecDP

pure function dotDP (v1,v2)
type (pointDP), intent (IN):: v1,v2
real*8 :: dotDP
dotDP=v1%x*v2%x+v1%y*v2%y+v1%z*v2%z
end function dotDP

pure function crossDP(v1,v2)
type (pointDP) crossDP
type (pointDP), intent (IN) :: v1,v2
 crossDP=pointDP(-v1%z*v2%y + v1%y*v2%z, v1%z*v2%x - v1%x*v2%z, -v1%y*v2%x + v1%x*v2%y)
end function crossDP

pure function moduleDP (v)
real*8 :: moduleDP
type (pointDP), intent (IN) :: v
moduleDP = dsqrt (dotDP(v,v))
end function moduleDP

pure function makeUnitDP (v)
type(pointDP):: makeUnitDP
type (pointDP), intent (IN) :: v
makeUnitDP = (1.e0/moduleDP(v))*v
end function makeUnitDP

pure function cosangDP (v1,v2)
real*8 :: cosangDP
type (pointDP), intent (IN) :: v1,v2
cosangDP=(v1*v2)/moduleDP(v1)/moduleDP(v2)
end function cosangDP

pure FUNCTION calcDist2DP (p1,p2)
  IMPLICIT NONE
  REAL*8 :: calcDist2DP
  TYPE (pointDP), INTENT(IN) :: p1,p2
  type (pointDP) :: v
  v = restaVecDP(p1,p2)
  calcDist2DP = v * v
END FUNCTION calcDist2DP

pure FUNCTION calcDistDP (p1,p2)
  IMPLICIT NONE
  REAL*8 :: calcDistDP
  TYPE (pointDP), INTENT(IN) :: p1,p2
  calcDistDP = DSQRT (calcDist2DP (p1,p2))
END FUNCTION calcDistDP

END MODULE geometryDP
