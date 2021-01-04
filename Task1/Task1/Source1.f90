    use ivprk_int
    use ivpag_int
    use umach_int
    !use sset_int
    
    implicit none
    include 'link_fnl_shared_imsl.h'
    
    integer i, N, ido,istep, nout, mxparm, j
    parameter (mxparm=50, N = 20)
    
    !INTEGER MABSE, MBDF, MSOLVE
    !PARAMETER (MABSE=1, MBDF=2, MSOLVE=2)
    
    real A(N,N), param(mxparm), tend, tol, MABSE, MBDF, MSOLVE, param2(mxparm)
    real fi, h, L, T, v(N), error(N), error_ivprk(N), error_ivpag(N), max_err_ivprk, u(N, N), max_err_ivpag, T_max, h_time
    
    external SSET, FCN, FCNJ
    
    call umach (2, nout)
    
    MABSE=1.0
    MBDF=2.0
    MSOLVE=2.0
    
    L = 1.0
    T_max = 1.0
    T = 0.0
    
    h = L/N
    h_time = T_max / N
    
    !test 1
    !do i = 1, N
    !     do j = 1, N
    !        u(i, j) =  3
    !    end do
    !    v(i) =  3.0
    !end do
    !u = 3.0
    
    !test 2
    !do i = 1, N
    !    v(i) =  h * i
    !    do j = 1, N
    !        u(i, j) =  h * i * exp(j * h_time)
    !    end do
    !end do
    
    !test 3
    !do i = 1, N
    !    v(i) =  h * i
    !    do j = 1, N
    !        u(i, j) =  h * i * exp(- j * h_time)
    !    end do
    !end do
    
    !test 3
    !do i = 1, N
    !    v(i) =  0.0
    !    do j = 1, N
    !        u(i, j) =  h * i * j * h_time
    !    end do
    !end do
    
    !test 4
    !do i = 1, N
    !    v(i) =  0.0
    !    do j = 1, N
    !        u(i, j) =  j * 1.0
    !    end do
    !end do
    
    !test 5
    do i = 1, N
         v(i) =  sin(h * i)
       do j = 1, N
            u(i, j) =  sin(h * i) * exp(- j * h_time)
        end do
    end do
    
    tol = 0.005
    
    CALL SSET (MXPARM, 0.0, PARAM, 1)
    CALL SSET (MXPARM, 0.0, PARAM2, 1)

    param(10) = 1.0
    
    
    WRITE (nout, 99998) 
    
      ido = 1
      istep = 0
      WRITE (nout,'(I6,11F12.3)') istep, T, v

   10 CONTINUE

      istep = istep + 1

      tend = h_time * istep

      CALL IVPRK (ido, FCN, T, tend, v, TOL=tol, PARAM=param)

      IF (istep .LE. N) THEN

         WRITE (nout,'(I6,11F12.3)') istep, T, v
         do i = 1, N
             error(i) = abs(u(i, istep) - v(i))
        end do
         error_ivprk(istep) = maxval(error)

!                                 Final call to release workspace

         IF (istep .EQ. N) ido = 3
              

         GO TO 10
    END IF
    WRITE (nout,99999) param(35)
    max_err_ivprk = maxval(error_ivprk)
    WRITE (nout,*) 'Max approximate error = ', max_err_ivprk

99998 FORMAT (4X, 'ISTEP', 5X, 'Time', 9X, 'V')
99999 FORMAT (4X, 'Number of fcn calls =', F6.0)   
      
      write(nout, *) ''
      
      T = 0.0
      
    !test 1
    !do i = 1, N
    !     do j = 1, N
    !        u(i, j) =  3
    !    end do
    !    v(i) =  3.0
    !end do
    
    !test 2
    !do i = 1, N
    !    v(i) =  h * i
    !    do j = 1, N
    !        u(i, j) =  h * i * exp(j * h_time)
    !    end do
    !end do
      
    !test 3
    !do i = 1, N
    !    v(i) =  h * i
    !    do j = 1, N
    !        u(i, j) =  h * i * exp(- j * h_time)
    !    end do
    !end do
    
    !test 3
    !do i = 1, N
    !    v(i) =  0.0
    !    do j = 1, N
    !        u(i, j) =  h * i * j * h_time
    !    end do
    !end do
    
    !test 4
    !do i = 1, N
    !    v(i) =  0.0
    !    do j = 1, N
    !        u(i, j) =  j * 1.0
    !    end do
    !end do
      
    !test 5
    do i = 1, N
         v(i) =  sin(h * i)
       do j = 1, N
            u(i, j) =  sin(h * i) * exp(- j * h_time)
        end do
    end do
      
      param2(10) = MABSE
      param2(12) = 2
      param2(13) = 1
      param2(19) = 0
      
      WRITE (nout,99998)
      
      ido = 1
      istep = 0
      WRITE (nout,'(I6,11F12.3)') istep, T, v
   11 CONTINUE
      istep = istep + 1
      tend = h_time * istep
      
      CALL IVPAG (ido, FCN, FCNJ, T, tend, v, TOL=tol, PARAM=param2)
      IF (istep .LE. N) THEN
          WRITE (nout,'(I6,11F12.3)') istep, T, v
          do i = 1, N
             error(i) = abs(u(i, istep) - v(i))
        end do
         error_ivpag(istep) = maxval(error)
          IF (istep .EQ. N) ido = 3
          GO TO 11
    END IF
    WRITE (nout,99999) param2(35)
    max_err_ivpag = maxval(error_ivpag)
    WRITE (nout,*) 'Max approximate error = ', max_err_ivpag
      
      END

      SUBROUTINE FCN (N, T, v, vprime)

!                                 SPECIFICATIONS FOR ARGUMENTS

      INTEGER    N

      REAL       T, v(N), vprime(N)
      real p, k_test, q_test, f_test, nu(2)
      
      real k(N), f(N), q(N), h2(N)
      
      !external COS
      
      h = L/N
    
      h2(1) = h / 2
      h2(N) = h / 2
      do i = 2, N-1
          h2(i) = h
      end do
      
      !test 1
      !k_test = 5.0
      !q_test = 1.0
      !f_test = 3.0
      !nu(1) = 0.0
      !nu(2) = 0.0
      
      !test 2
      !k_test = 5.0
      !q_test = - 1.0
      !f_test = 0.0
      !nu(1) = 5.0 * exp(T)
      !nu(2) = 5.0 * exp(T)
      
      !test 3
      !k_test = 5.0
      !q_test = 1.0
      !f_test = 0.0
      !nu(1) = 5.0 * exp(- T)
      !nu(2) = 5.0 * exp(- T)
      
      
      !!test 3
      !  k_test = 5.0
      !  !q_test = -1.0 * (T ** (-1.0))
      !  q_test = 0.0
      !  f_test = 0.0
      !  nu(1) = 5.0 * T
      !  nu(2) = 5.0 * T
      
      !test 4
      !k_test = 5.0
      !q_test = - (1.0 / T)
      !f_test = 0.0
      !nu(1) = 0.0
      !nu(2) = 0.0
      
      !test 5
      k_test = 2.0
      q_test = - 1.0
      f_test = 0.0
      nu(1) = 2.0 * exp(-T)
      nu(2) = 0.0
    
      p = h
      do i = 1, N
          k(i) = k_test
          q(i) = q_test
          f(i) = f_test
          !f(i) = i * h
      end do
      
      
      vprime(1) = ((k(1) * ((v(2)-v(1)) / (h2(1) * p)) + (nu(1) / h2(1))) - q(1)*v(1)) + f(1)
      vprime(N) = (((- nu(2) / h2(N))  - (k(N-1) * ((v(N)-v(N-1)) / (h2(N) * p)))) - q(N)*v(N)) + f(N)

      do i = 2, N - 1
          vprime(i) = ((k(i+1) * ((v(i+1)-v(i)) / (h2(i) * p)) - k(i-1) * ((v(i)-v(i-1)) / (h2(i) * p))) - q(i)*v(i)) + f(i)
      end do

      RETURN
      
    end
    
    
    SUBROUTINE FCNJ (N, T, v, DYPDY)

        INTEGER N
        REAL T, v(N), DYPDY(N,*)
        
        real p, k_test, q_test, f_test, nu(2)
      
        real k(N), f(N), q(N), h2(N)
        
        !external COS
      
        h = L/N
    
        h2(1) = h / 2
        h2(N) = h / 2
        do i = 2, N-1
            h2(i) = h
        end do
      
        !test 1
        !k_test = 5.0
        !q_test = 1.0
        !f_test = 3.0
        !nu(1) = 0.0
        !nu(2) = 0.0
        
        !test 2
        !k_test = 5.0
        !q_test = - 1.0
        !f_test = 0.0
        !nu(1) = 5.0 * exp(T)
        !nu(2) = 5.0 * exp(T)
        
        !test 3
        !k_test = 5.0
        !q_test = 1.0
        !f_test = 0.0
        !nu(1) = 5.0 * exp(- T)
        !nu(2) = 5.0 * exp(- T)
        
        !test 3
        !k_test = 5.0
        !!q_test = -1.0 * (T ** (-1.0))
        !q_test = 0.0
        !f_test = 0.0
        !nu(1) = 5.0 * T
        !nu(2) = 5.0 * T
        
        !test 4
        !k_test = 5.0
        !q_test = - (1.0 / T)
        !f_test = 0.0
        !nu(1) = 0.0
        !nu(2) = 0.0
        
        !test 5
        k_test = 2.0
        q_test = - 1.0
        f_test = 0.0
        nu(1) = 2.0 * exp(-T)
        nu(2) = 0.0
        
        p = h
        do i = 1, N
            k(i) = k_test
            q(i) = q_test
            f(i) = f_test
            !f(i) = h * i
        end do
        
        CALL SSET (N**2, 0.0, DYPDY, 1)
        
        DYPDY(1,1) = (- k(1) / (h2(1) * p)) - q(1)
        DYPDY(1,2) = (k(1) / (h2(1) * p))
        
        do i = 2, N-1
            DYPDY(i, i-1) = k(i-1) / (h2(i) * p)
            DYPDY(i, i) = (- k(i+1) / (h2(i) * p)) - (k(i-1) / (h2(i) * p)) - q(i)
            DYPDY(i, i+1) = k(i+1) / (h2(i) * p)
            
        end do
        
        DYPDY(N,N-1) = k(N-1) / (h2(N) * p)
        DYPDY(N,N) = (-k(N-1) / (h2(N) * p)) - q(N)

        RETURN
    END