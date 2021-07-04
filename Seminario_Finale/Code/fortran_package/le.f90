      subroutine le(l)
          interface lv4d
              function lv4d(x, r, a) result(dx)
                  integer                 :: i, j
                  real*8                  :: s(4)
                  real*8                  :: dx(4)
                  real*8, intent(in)      :: a(4,4)
                  real*8, intent(in)      :: r(4)
                  real*8, intent(in)      :: x(4)
              end function
          end interface
          interface jac
              function jac(x, r, a)
                  integer                 :: i, j
                  real*8                  :: s(4)
                  real*8                  :: jac(4,4)
                  real*8                  :: dx(4)
                  real*8, intent(in)      :: a(4,4)
                  real*8, intent(in)      :: r(4)
                  real*8, intent(in)      :: x(4)
              end function
          end interface
          integer              :: i, j
          real*8               :: a(4,4), jmat(4, 4)
          real*8               :: r(4)
          real*8               :: x(4, 1000)
          real*8               :: dx(4)
          real*8, intent(inout):: l(4)
          data a  / 1d0, 1d0, 1d0, 1d0, &
                    1d0, 1d0, 1d0, 1d0, &
                    1d0, 1d0, 1d0, 1d0, &
                    1d0, 1d0, 1d0, 1d0 / 
          data r  / 1d0, 1d0, 1d0, 1d0 /
          data x(:, 1) / 0.1d0, 0.1d0, 0.1d0, 0.1d0 /
          x(:, 1) = x(:, 1) + lv4d(x(:, 1), r, a)
          jmat = jac(x(:, 1), r, a)
      end subroutine

      function lv4d(x, r, a) result(dx)
      !========================================
      ! Funzione che crea la un Lotka Volterra 4D 
      ! competitivo dati i parametri iniziali,
      ! restituisce gli array contenenti le x e le y.
      !========================================
          implicit none
          integer                 :: i, j
          real*8                  :: s(4)
          real*8                  :: dx(4)
          real*8, intent(in)      :: a(4,4)
          real*8, intent(in)      :: r(4)
          real*8, intent(in)      :: x(4)
 
          do j = 1, 4
              do i = 1, 4
                  s(j) = s(j) + a(j, i) * x(i)
              enddo
              dx(j) = r(j) * x(j) * (1d0 - s(j)) 
          enddo
      end function

      function jac(x, r, a) 
          implicit none
          integer                 :: i, j
          real*8                  :: s(4)
          real*8                  :: jac(4,4), ff
          real*8, intent(in)      :: a(4,4)
          real*8, intent(in)      :: r(4)
          real*8, intent(in)      :: x(4)
          do i = 1, 4
              do j = 1, 4
                  s(i) = s(i) + a(i, j) * x(j)
              enddo
              do j = 1, 4
                  ff = r(i) * x(j) * a(i,j)
                  if (i == j) then
                      jac(i, j) = r(i) * x(i) * (1 - s(i)) - ff
                  else
                      jac(i, j) = r(i) * x(i) * (1 - s(i)) 
                  endif
              enddo
          enddo
      end function

      function rg4(x, t1, t2, r, a) result(dx)
          interface
              function lv4d(x, r, a) result(dx)
                  integer                 :: i, j
                  real*8                  :: s(4)
                  real*8                  :: dx(4)
                  real*8, intent(in)      :: a(4,4)
                  real*8, intent(in)      :: r(4)
                  real*8, intent(in)      :: x(4)
              end function
          end interface
          integer                 :: i, j
          real*8                  :: dt, tmid
          real*8, dimension(4)    :: s, dx, k1, k2, k3, k4
          real*8, intent(in)      :: r(4), x(4), t1, t2, a(4,4)
          dt = t2 - t1
          tmid = (t2 + t1)/2
          k1 = lv4d(x, r, a)
          k2 = lv4d(x + dt * k1/2d0, r, a)
          k3 = lv4d(x + dt * k2/2d0, r, a)
          k4 = lv4d(x + dt * k3, r, a)
          dx = dt * (k1/2d0 + k2 + k3 + k4/2d0)/3d0
      end function

      function varrg4(x, t1, t2, r, a) result(dx)
          interface
              function jac(x, r, a) result(dx)
                  integer                 :: i, j
                  real*8                  :: s(4)
                  real*8                  :: dx(4)
                  real*8, intent(in)      :: a(4,4)
                  real*8, intent(in)      :: r(4)
                  real*8, intent(in)      :: x(4)
              end function
          end interface
          integer                 :: i, j
          real*8                  :: dt, tmid
          real*8, dimension(4)    :: s, dx, k1, k2, k3, k4
          real*8, intent(in)      :: r(4), x(4), t1, t2, a(4,4)
          dt = t2-t1
          tmid = (t2 + t1)/2
          k1 = jac(x, r, a)
          k2 = jac(x + dt * k1/2d0, r, a)
          k3 = jac(x + dt * k2/2d0, r, a)
          k4 = jac(x + dt * k3, r, a)
          dx = dt * (k1/2d0 + k2 + k3 + k4/2d0)/3d0
      end function


      function evalspectr(x0, t, r, a, trans, n, m) result(s)
          interface
              function rg4(x, t1, t2, r, a) result(dx)
                  integer                 :: i, j
                  real*8                  :: dt, tmid
                  real*8, dimension(4)    :: s, dx, k1, k2, k3, k4
                  real*8, intent(in)      :: r(4), x(4), t1, t2, a(4,4)
              end function
          end interface
          interface
              function lv4d(x, r, a) result(dx)
                  integer                 :: i, j
                  real*8                  :: s(4)
                  real*8                  :: dx(4)
                  real*8, intent(in)      :: a(4,4)
                  real*8, intent(in)      :: r(4)
                  real*8, intent(in)      :: x(4)
              end function
          end interface
          interface
              function jac(x, r, a) 
                  integer                 :: i, j
                  real*8                  :: s(4)
                  real*8                  :: jac(4,4), ff
                  real*8, intent(in)      :: a(4,4)
                  real*8, intent(in)      :: r(4)
                  real*8, intent(in)      :: x(4)
              end function
          end interface
          integer                 :: i, j
          integer, intent(in)     :: n, m
          real*8                  :: D = 4, rphi(4, 4), rdphi(4), dt, phi0(4,4)
          real*8, dimension(4)    :: s, dx, k1, k2, k3, k4
          real*8, intent(in)      :: r(4), x0(4), t(n), a(4,4), trans(m)
          dt = t(2)-t(1)
          ! Building Identity matrix phi0
          do i = 1, 4
              do j = 1, 4
                  if (i == j) then
                      phi0(i,j) = 1d0
                  else
                      phi0(i, j) = 0d0
                  endif
              enddo
          enddo
      end function
