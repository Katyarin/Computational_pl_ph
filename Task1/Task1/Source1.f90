    use ivprk_int
    use umach_int
    use sset_int
    
    implicit none
    include 'link_fnl_shared_imsl.h'
    
    integer i, N, ido,istep, nout, mxparm
    parameter (mxparm=50, N = 10)
    real param(mxparm), tend, tol
    real fi, h, L, T, v(N)
    
    external FCN
    
    call umach (2, nout)
    
    !N = 10
    L = 1.0
    T = 0.0
    do i = 1, N
        v(i) =  3.0
    end do
    
    tol = 0.005
    
    CALL SSET (MXPARM, 0.0, PARAM, 1)

    param(10) = 1.0
    
    
    WRITE (nout, 99998) 
    
      ido = 1
      istep = 0
      WRITE (nout,'(I6,11F12.3)') istep, T, v

   10 CONTINUE

      istep = istep + 1

      tend = istep

      CALL IVPRK (ido, FCN, T, tend, v, TOL=tol, PARAM=param)

      IF (istep .LE. 10) THEN

         WRITE (nout,'(I6,11F12.3)') istep, T, v

!                                 Final call to release workspace

         IF (istep .EQ. 10) ido = 3

         GO TO 10
    END IF
    WRITE (NOUT,99999) PARAM(35)

99998 FORMAT (4X, 'ISTEP', 5X, 'Time', 9X, 'V')
99999 FORMAT (4X, 'Number of fcn calls with IVPRK =', F6.0)
      
      write(nout, *) 'This is the end'
      
      

      END

      SUBROUTINE FCN (N, T, v, vprime)

!                                 SPECIFICATIONS FOR ARGUMENTS

      INTEGER    N

      REAL       T, v(N), vprime(N)
      real p, k_test, q_test, f_test, nu(2)
      
      real k(N), f(N), q(N), h2(N)
      
      h = L/N
    
      h2(1) = h / 2
      h2(N) = h / 2
      do i = 2, N-1
          h2(i) = h
      end do
      
      !test 1
      k_test = 5.0
      q_test = 1.0
      f_test = 3.0
      nu(1) = 0.0
      nu(2) = 0.0
    
      p = 2 * h
      do i = 1, N
          k(i) = k_test
          q(i) = q_test
          f(i) = f_test
      end do
      
      
      
      vprime(1) = ((k(1) * ((v(2)-v(1)) / (h2(1) * p)) + (nu(1) / h2(1))) - q(1)*v(1)) + f(1)
      vprime(N) = (((- nu(2) / h2(N))  - (k(N-1) * ((v(N)-v(N-1)) / (h2(N) * p)))) - q(N)*v(N)) + f(N)

      do i = 2, N - 1
          vprime(i) = ((k(i+1) * ((v(i+1)-v(i)) / (h2(i) * p)) - k(i-1) * ((v(i)-v(i-1)) / (h2(i) * p))) - q(i)*v(i)) + f(i)
      end do

      RETURN
      
    end