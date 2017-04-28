! Modification of fmm.f90 from Public Domain Aeronautical Software
! Modification copyright (c) 2015, 2016 by Kernelyze LLC
! Author: Thomas A. Knox 
! This program is free software: you can redistribute it and/or modify 
! it under the terms of the GNU Affero General Public License as 
! published by the Free Software Foundation, either version 3 of the 
! License, or (at your option) any later version. 
! 
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! GNU Affero General Public License for more details. 
!  
! You should have received a copy of the GNU Affero General Public License 
! along with this program.  If not, see <http://www.gnu.org/licenses/>. 
!
! Modification:
! created on: 2015-02-07
! updated on: 2015-04-22
! updated on: 2016-04-24 (Avoid upstream use of internal procedures 
!             for thread safety, so use singlevar_func objects
!             rather than direct function interfaces to specify
!             the objective functions in these procedures)
!
! Information for fmm.f90:
! PURPOSE - A collection of Fortran procedures for mathematical computation
!   based on the procedures from the book "Computer Methods for Mathematical
!   Computations", by George E. Forsythe, Michael A. Malcolm, and 
!   Cleve B. Moler. Prentice-Hall, 1977.

!  COLLECTION AUTHORS - George E. Forsythe, Michael A. Malcolm, & 
!    Cleve B. Moler, Stanford University (1977)
!    Ralph L. Carmichael, Public Domain Aeronautical Software

! ORIGINAL AUTHORS - 
!  Decomp,Solve - Forsythe,Malcolm,Moler
!  Spline,Seval - Forsythe,Malcolm,Moler
!  Quanc8 - Forsythe,Malcolm,Moler
!  RKF45 - H.A. Watts AND L.F. Shampine  (Sandia)
!  ZEROIN - Richard Brent
!  FMIN - Van Winjngaarden,Decker,Brent, Et Al
!  SVD - Golub and Reinsch

MODULE brent_mod
      
use set_precision, only : wp
use constants_mod, only : err_msg_len
use singlevar_func_mod, only : singlevar_func
use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

REAL(wp), PARAMETER :: ZERO = 0.0E0_wp, HALF = 0.5E0_wp, ONE = 1.0E0_wp, &
  TWO = 2.0E0_wp, THREE = 3.0E0_wp

CONTAINS
  
!+
pure function zeroin (a,b,f,tol) RESULT(z)
! ---------------------------------------------------------------------------
! PURPOSE - Compute a zero of f in the interval (a,b)
  ! Arguments
  REAL(wp),INTENT(IN):: a,b   ! left and right endpoints of interval
  REAL(wp),INTENT(IN):: tol     ! desired interval of uncertainity 
  class(singlevar_func), intent(in) :: f ! Function to find a zero of
  ! Result
  REAL(wp):: z    
  ! Body
  INTEGER:: errCode
  INTEGER,PARAMETER:: max_iter=200
  INTEGER:: neval   ! not used
  REAL(wp):: xZero,fZero

!----------------------------------------------------------------------------
  CALL brentZero(a,b,f,tol,max_iter,neval,errCode,xZero,fZero)
  z=xZero
  RETURN
END Function zeroin   ! -----------------------------------------------------

!+
pure subroutine brentZero(ax,bx,f,tol,max_iter,neval,errCode,xZero,fZero)
! ---------------------------------------------------------------------------
! PURPOSE - Compute a zero of F in the interval (ax,bx)

  REAL(wp),INTENT(IN):: ax,bx   ! left and right endpoints of interval
  class(singlevar_func), intent(in) :: f ! Function to find a zero of
  REAL(wp),INTENT(IN):: tol     ! desired interval of uncertainity 
  INTEGER,INTENT(IN):: max_iter   ! max number of iterations allowed. 25 is good
  INTEGER,INTENT(OUT):: neval
  INTEGER,INTENT(OUT):: errCode   ! =0 is OK; =1 too many iterations
                                  ! =2 if F(ax) and F(bx) have the same sign
  REAL(wp),INTENT(OUT):: xZero,fZero ! the last and best value of the zero                                   

  REAL(wp):: a,b,c,d,e,eps
  REAL(wp):: fa,fb,fc,tol1
  INTEGER:: kIter
!  INTEGER:: method   ! =0 bisection; =1 linear; =2 inverse quadratic
  REAL(wp):: xm,p,q,r,s
!----------------------------------------------------------------------------
  eps=EPSILON(ax)
  tol1=ONE+eps

  a=ax   ! initialization
  b=bx
  fa=f%eval(a)
  fb=f%eval(b)
  neval=2
  IF (fa==ZERO) THEN
    xZero=a
    fZero=ZERO
    errCode=0
    RETURN
  END IF  
  IF (fb==ZERO) THEN
    xZero=b
    fZero=ZERO
    errCode=0
    RETURN
  END IF  
  IF (fa*fb > ZERO) THEN
    xZero = ieee_value(a, ieee_quiet_nan)
    fZero = ieee_value(a, ieee_quiet_nan)
    errCode=2
    RETURN
  END IF
! The trivial cases have now been dealt with. On to the real work...

  c=a
  fc=fa
  d=b-a
  e=d
    
  DO kIter=1,max_iter
    IF ( (fb>0 .AND. fc>0) .OR. (fb<0 .AND. fc<0) ) THEN
      c=a  ! we insist that b and c straddle the zero
      fc=fa
      d=b-a
      e=d
    END IF

    IF (ABS(fc) < ABS(fb)) THEN
      a=b    ! we insist that b be the better guess of b and c
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
    END IF

    tol1=TWO*eps*ABS(b)+HALF*tol   ! convergence test
    xm=HALF*(c-b)
    IF (ABS(xm) <= tol1 .OR. fb==ZERO) THEN
      xZero=b
      fZero=fb
      errCode=0   ! SUCCESS! The proper way to leave
      RETURN
    END IF

!    WRITE(DBG,*) "BrentZero, start of kIter=",kIter
!    WRITE(DBG,'(A,2ES25.16)') "a,fa=", a,fa
!    WRITE(DBG,'(A,2ES25.16)') "b,fb=", b,fb
!    WRITE(DBG,'(A,2ES25.16)') "c,fc=", c,fc
!    WRITE(DBG,'(A,2ES25.16)') "xm,tol1=", xm,tol1

    IF (ABS(e) < tol1 .OR. ABS(fa) <= ABS(fb) ) THEN
      d=xm   ! bisection
      e=d
!      method=0
    ELSE
      IF (a==c) THEN
        s=fb/fa   ! linear interpolation
        p=TWO*xm*s
        q=ONE-s
!        method=1
      ELSE
        q=fa/fc   ! inverse quadratic interpolation
        r=fb/fc
        s=fb/fa
        p=s*(TWO*xm*q*(q-r)-(b-a)*(r-ONE))
        q=(q-ONE)*(r-ONE)*(s-ONE)
!        method=2
      END IF
      IF (p > ZERO) q=-q   ! adjust signs
      p=ABS(p)
      IF (p+p >= (THREE*xm*q-ABS(tol1*q)) .OR. p+p >= ABS(e*q) ) THEN
        d=xm   ! don't interpolate. Use bisection
        e=d
 !       method=-1
      ELSE
        e=d   ! OK, use interpolation
        d=p/q
      END IF
    END IF  
   
    a=b   ! complete step. a becomes the previous iteration
    fa=fb
   IF (ABS(d) > tol1) THEN
     b=b+d
   ELSE  
     b=b+SIGN(tol1,xm)
   END IF
     
   fb=f%eval(b)   ! the newest and best value (we hope)
   neval=neval+1   ! keep count of the function evaluations
!!   IF ((fb*(fc/ABS(fc))) > ZERO) GO TO 20
!    WRITE(DBG,*) "BrentZero, end of kIter=",kIter, "   method=",method
!    WRITE(DBG,'(A,2ES25.16)') " new b,fb=", b,fb
!    WRITE(DBG,'(A,2ES25.16)') "d,e=", d,e
!    WRITE(DBG,'(A,2ES25.16)') "p,q=", p,q
  END DO
  
! The loop should never terminate. If it does, return a quiet NaN
!  and set errCode to 1
  xZero = ieee_value(a, ieee_quiet_nan)
  fZero = ieee_value(a, ieee_quiet_nan)
  errCode=1
  RETURN  

END Subroutine brentZero   ! ------------------------------------------------

!+
pure subroutine fmin( &
    ax, &
    bx, & 
    f, &
    tol, &
    xopt, &
    fopt, &
    err_stat, &
    err_msg)
! ---------------------------------------------------------------------------
! PURPOSE -

  REAL(wp),INTENT(IN)  :: ax      !  left endpoint of initial interval
  REAL(wp),INTENT(IN)  :: bx      !  right endpoint of initial interval
  class(singlevar_func), intent(in) :: f ! Function to find a min of
  REAL(wp),INTENT(IN)  :: tol     !  desired length of interval of uncertainity
  REAL(wp),INTENT(OUT) :: xopt    !  the x-coordinate of the minimum
  REAL(wp),INTENT(OUT) :: fopt    !  the func value at the minimum
  
  ! Optional arguments to pass back information on any errors
  integer, intent(out), optional                    :: err_stat
  character(len=err_msg_len), intent(out), optional :: err_msg
  
  INTEGER:: errCode

  INTEGER,PARAMETER:: MAX_ITER = 200
  INTEGER:: neval
!----------------------------------------------------------------------------
  
  ! Initialize the optional error arguments if present
  if (present(err_stat)) then
    err_stat = 0
  end if
  if (present(err_msg)) then
    err_msg = ''
  end if
  
  CALL brentmin(ax,bx,f,tol,MAX_ITER,neval,errCode,xopt,fopt)
  
  if (present(err_stat)) then
    err_stat = errCode
  end if
  
  if (errCode == 1) then
    if (present(err_msg)) then
      err_msg = 'Brent minimizer, maximum # iterations exceeded'
    end if
  end if
  
END subroutine fmin   ! -------------------------------------------------------

!+
pure subroutine brentmin(ax,bx,f,tol,max_iter,neval,errCode,xZero,fZero)
! ---------------------------------------------------------------------------
! PURPOSE -

  REAL(wp),INTENT(IN):: ax   !  left endpoint of initial interval
  REAL(wp),INTENT(IN):: bx   !  right endpoint of initial interval
  class(singlevar_func), intent(in) :: f ! Function to find a min of
  REAL(wp),INTENT(IN):: tol  !  desired length of interval of uncertainity
  INTEGER,INTENT(IN):: max_iter   ! max number of iterations allowed. 25 is good
  INTEGER,INTENT(OUT):: neval
  INTEGER,INTENT(OUT):: errCode   ! =0 is OK; =1 too many iterations
  REAL(wp),INTENT(OUT):: xZero,fZero ! the last and best value of the minimum                                  
      
  REAL(wp), PARAMETER :: C=0.381966011_wp ! (3-Sqrt(5))/2

  REAL(wp) :: a,b,d,e,eps,xm,p,q,r,tol1,tol2,u,v,w
  REAL(wp) :: fu,fv,fw,fx,x
  INTEGER :: iter
!----------------------------------------------------------------------------
  errCode=0
  eps=SQRT(EPSILON(ONE))

  a=ax ! initialization
  b=bx
  v=a+c*(b-a)
  w=v
  x=v
  e=ZERO
  fx=f%eval(x)
  fv=f%eval(v)
  neval=2
  fw=fx

  DO iter=1,max_iter   ! iterate until minimum is found
    xm=HALF*(a+b)
    tol1=eps*ABS(x)+tol/THREE
    tol2=tol1+tol1
    IF (ABS(x-xm)<=(tol2-HALF*(b-a)) ) EXIT    ! check stopping criterion.

    IF (ABS(e) > tol1) THEN   ! is golden section necessary???
      r=(x-w)*(fx-fv)   ! trial parabolic fit
      q=(x-v)*(fx-fw)
      p=(x-v)*q-(x-w)*r
      q=TWO*(q-r)
      IF (q > ZERO) p=-p
      q=ABS(q)
      r=e
      e=d
      IF (ABS(p)>=ABS(HALF*q*r) .OR. p<=q*(a-x) .OR. p>=q*(b-x)) GO TO 40

      d=p/q   ! a parabolic interpolation step
      u=x+d
      IF ((u-a)<tol2 .OR. (b-u)<tol2) D=SIGN(tol1,xm-x)
      GO TO 50
    END IF

!40       CONTINUE
 40 IF (x >= xm) THEN   ! a golden section step
      e=a-x
    ELSE
      e=b-x
    END IF
    d=c*e

50  IF (ABS(d) >= tol1) THEN   ! f must not be evaluated too close to x
      u=x+d
    ELSE
      u=x+SIGN(tol1,d)
    END IF
    fu=f%eval(u)
    neval=neval+1

    IF (fu <= fx) THEN   ! update a,b,v,w, and x
      IF (u >= x) THEN
        a=x
      ELSE
        b=x
      END IF
      v=w
      fv=fw
      w=x
      fw=fx
      x=u
      fx=fu
    ELSE
      IF (u<x) THEN
        a=u
      ELSE
        b=u
      END IF
      IF (fu<=fw .OR. w==x) THEN
        v=w
        fv=fw
        w=u
        fw=fu
      ELSE IF (fu<=fv .OR. v==x .OR. v==w) THEN
        v=u
        fv=fu
      END IF
    END IF
  END DO
  ! end of main loop
   
   ! If too many iterations have occurred,
   ! return NaN and an appropriate (nonzero)
   ! error code
   if (iter > max_iter) then
     xZero = ieee_value(x , ieee_quiet_nan)
     fZero = ieee_value(fx , ieee_quiet_nan)
     errCode = 1
   ! Otherwise, return with the values
   ! computed (note that errCode is
   ! already set to 0)
   else
     xZero=x   
     fZero=fx
   end if
  RETURN
END Subroutine brentmin   ! -------------------------------------------------

end module brent_mod