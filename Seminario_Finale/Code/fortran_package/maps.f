      !=================================================
      ! Funzione di main ausiliaria per eventuale debug
      !=================================================
      ! program main
      !     integer              :: s = 100000
      !     real*8               :: a0 = 0.1
      !     real*8               :: x0 = 0.1
      !     real*8               :: y0 = 0.1
      !     call henonmap(s, x0, y0, a0)
      ! end

      ! Per rendere eseguibile in python lanciare il comando (unix):
      ! f2py -c -m maps maps.f --fcompiler=gnu95



      subroutine lv4d(n, x, y, z, w, r, a, h)
      !========================================
      ! Funzione che crea la un Lotka Volterra 4D 
      ! competitivo dati i parametri iniziali,
      ! restituisce gli array contenenti le x e le y.
      !========================================
          integer                 :: i, j
          integer, intent(in)     :: n
          real*8                  :: s(4)
          real*8, intent(in)      :: a(4,4)
          real*8, intent(in)      :: r(4), h
          real*8, intent(inout)   :: x(n)
          real*8, intent(inout)   :: y(n)
          real*8, intent(inout)   :: z(n)
          real*8, intent(inout)   :: w(n)
          real*8                  :: x1, y1, z1, w1
          real*8                  :: x2, y2, z2, w2, s1, s2
 
          do i = 1, size(x)-1
              do j = 1, 4
                  s1 = a(j, 1) * x(i) + a(j, 2) * y(i)
                  s2 = a(j, 3) * z(i) + a(j, 4) * w(i)
                  s(j) =  s1 + s2
              enddo
              x1 = x(i) + r(1) * x(i) * (1 - s(1)) * h
              y1 = y(i) + r(2) * y(i) * (1 - s(2)) * h
              z1 = z(i) + r(3) * z(i) * (1 - s(3)) * h
              w1 = w(i) + r(4) * w(i) * (1 - s(4)) * h

              do j = 1, 4
                  s1 = a(j, 1) * x1 + a(j, 2) * y1
                  s2 = a(j, 3) * z1 + a(j, 4) * w1
                  s(j) = s1 + s2
              enddo
              x2 = x(i) +r(1) * x1 * (1 - s(1)) * h
              y2 = y(i) +r(2) * y1 * (1 - s(2)) * h
              z2 = z(i) +r(3) * z1 * (1 - s(3)) * h
              w2 = w(i) +r(4) * w1 * (1 - s(4)) * h

              x(i+1) = 0.5d0 * (x1 + x2)
              y(i+1) = 0.5d0 * (y1 + y2)
              z(i+1) = 0.5d0 * (z1 + z2)
              w(i+1) = 0.5d0 * (w1 + w2)
          enddo
      end subroutine


      subroutine henonmap(n, x, y, a0)
      !========================================
      ! Funzione che crea la mappa di Henon
      ! dati i parametri iniziali, restituisce
      ! gli array contenenti le x e le y.
      !========================================
          integer                 :: i
          integer, intent(in)     :: n
          real*8, intent(in)      :: a0
          real*8, intent(inout)   :: x(n)
          real*8, intent(inout)   :: y(n)
          do i = 1, size(x)-1
              x(i+1) = x(i)*dcos(a0) - ( y(i) - x(i)*x(i) )*dsin(a0)
              y(i+1) = x(i)*dsin(a0) + ( y(i) - x(i)*x(i) )*dcos(a0)
          enddo
      end subroutine

      subroutine standardmap(n, q, p, k)
      !========================================
      ! Funzione che crea la mappa standard 
      ! dati i parametri iniziali, restituisce  
      ! gli array contenenti le x e le y.
      !========================================
          integer                 :: i
          ! Definisco un pi preciso
          real*8                  :: pi = 4.d0*datan(1.d0)
          integer, intent(in)     :: n ! Numero di iterate da fare
          real*8, intent(in)      :: k ! Costante della mappa
          real*8, intent(inout)   :: q(n)
          real*8, intent(inout)   :: p(n)
          do i = 1, size(q)-1
              p(i+1) = mod(p(i) + k/(2*pi) * dsin(2 * pi * q(i)), 1.d0 )
              q(i+1) = mod(q(i) + p(i+1), 1.d0 )
          enddo
      end subroutine

      subroutine logisticmap(n, x, y, lamda)
      !========================================
      ! Funzione che crea la mappa logistica
      ! dati i parametri iniziali, restituisce
      ! gli array contenenti le x e le y.
      !========================================
          integer                 :: i
          integer, intent(in)     :: n
          real*8, intent(in)      :: lamda
          real*8, intent(inout)   :: x(n)
          real*8, intent(out)     :: y(n)
          y(1) = x(1)
          do i = 1, size(x)-1
              y(i+1) = 4 * lamda * y(i) * ( 1 - y(i) )
          enddo
      end subroutine

      subroutine nonHhenonmap(n, x, y, a, b)
      !========================================
      ! Funzione che crea la mappa di Henon
      ! (non Hamiltoninana)
      ! dati i parametri iniziali, restituisce
      ! gli array contenenti le x e le y.
      !========================================
          integer                 :: i
          integer, intent(in)     :: n
          real*8, intent(in)      :: a, b
          real*8, intent(inout)   :: x(n)
          real*8, intent(inout)     :: y(n)
          do i = 1, size(x)-1
              x(i+1) = 1 - a * x(i)**2 + y(i)
              y(i+1) = b * x(i)
          enddo
      end subroutine

      subroutine lorentz(n, x, y, z, h, s, r, b)
      !========================================
      ! Funzione che crea la mappa di Lorentz
      ! dati i parametri iniziali, restituisce
      ! gli array contenenti le x, y, z.
      !========================================
          integer                 :: i
          integer, intent(in)     :: n
          real*8, intent(in)      :: h, s, r, b
          real*8                  :: x1, y1, z1
          real*8                  :: x2, y2, z2
          real*8, intent(inout)   :: x(n)
          real*8, intent(inout)   :: y(n)
          real*8, intent(inout)   :: z(n)
          do i = 1, size(x)-1
              x1 = x(i) + h * (s * ( y(i) - x(i) ))
              y1 = y(i) + h * (-x(i) * z(i) + r * x(i) - y(i))
              z1 = z(i) + h * (x(i) * y(i) - b * z(i))

              x2 = x(i) + h * (s * ( y1 - x1 ))
              y2 = y(i) + h * (-x1 * z1 + r * x1 - y1)
              z2 = z(i) + h * (x1 * y1 - b * z1)
              
              x(i+1) = 0.5d0 * (x1 + x2)
              y(i+1) = 0.5d0 * (y1 + y2)
              z(i+1) = 0.5d0 * (z1 + z2)
          enddo
      end subroutine
