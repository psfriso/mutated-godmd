!     
! File:   ls_sq.f
! Author: psfriso
!
! Created on September 18, 2012, 4:43 PM
!
MODULE least_sq_fit
!
!  Purpose:
!    To perform a least-squares fit of an input data set
!    to a straight line, and print out the resulting slope
!    and intercept values.  The input data for this fit
!    comes from a user-specified input data file.
!
!
PRIVATE
PUBLIC lreg
!
CONTAINS
FUNCTION lreg(n,vector) RESULT (slope)
IMPLICIT NONE

! Data dictionary: declare variable types, definitions, & units  
! Note that cumulative variables are all initialized to zero.
INTEGER, INTENT(IN) :: n              ! Number of input data pairs (x,y)
REAL, INTENT(IN), DIMENSION(n,2) :: vector
REAL :: slope                 ! Slope of the line
REAL*8 :: sum_x             ! Sum of all input X values
REAL*8 :: sum_x2            ! Sum of all input X values squared
REAL*8 :: sum_xy            ! Sum of all input X*Y values
REAL*8 :: sum_y             ! Sum of all input Y values
REAL*8 :: x_bar                  ! Average X value
REAL*8 :: y_bar                 ! Average Y value
REAL*8 :: y_int                 ! Y-axis intercept of the line
INTEGER :: i                  ! Loop control index
!   Write(*,*) vector
   sum_x=0.0
   sum_y=0.0
   sum_x2=0.0
   sum_xy=0.0
   DO i=1,n
      sum_x  = sum_x + REAL(vector(i,1))               ! Calculate
     ! write(*,*) "x    ", sum_x
      sum_y  = sum_y + vector(i,2)               !   statistics
      sum_x2 = sum_x2 + REAL(vector(i,1))**2           !
      sum_xy = sum_xy + REAL(vector(i,1)) * vector(i,2)          !
   END DO
 slope=0.000
   ! Now calculate the slope and intercept. 
   x_bar=0.0
   y_bar=0.0
   y_int=0.0
   !
   x_bar = sum_x / REAL(n)
   y_bar = sum_y / REAL(n)
   slope = (sum_xy - sum_x * y_bar) / ( sum_x2 - sum_x * x_bar)
   y_int = y_bar - slope * x_bar 
!WRITE(*,'("SLOPE ",f5.2,1x, "xbar ", f6.2, " ybar ",f6.2, " sumx ",f8.2, " sumy ", f8.2, " sum_xy",f8.2, 1x, "yint ",f8.2)') &
!      slope, x_bar,y_bar , sum_x, sum_y , sum_xy,y_int
!
END FUNCTION lreg
END MODULE least_sq_fit
