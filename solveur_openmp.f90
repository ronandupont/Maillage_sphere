!==================
module solveur_openmp
!==================
!
contains
!

!=======================================================================
 subroutine conversionchainemorseter(lch,Aij,Aval,ligne)
    use deftype
    implicit none
    Type(listechaine), dimension(:), intent(in)       :: lch
    integer, dimension (:,:), intent(out)  :: Aij
    real (kind = 8), dimension(:), intent(out)  :: Aval
    integer, dimension (:), intent(inout)  :: ligne
    type (chainej), Pointer            :: chj
    integer :: i,l,k,tailleA
    

    k=0
    do i=1,size(lch)
       chj=>lch(i)%liste
       do while (associated(chj%prec))
          chj=>chj%prec
       end do

       ligne(i)=k+1
       do while (associated(chj))
          k=k+1
          Aij(k,2)=i
          Aij(k,1)=chj%j
          Aval(k)=chj%val

          chj=>chj%suiv
       end do
    end do
    ligne(size(lch)+1)=k+1
tailleA=k


!print*,'verification de la taille ter', tailleA
end subroutine conversionchainemorseter
!=======================================================================



!================================================================================================================!================================================================================================================
 subroutine solve_par(lch,u,b,niter,eps,taille)
!$USE OMP_LIB
use deftype
    implicit none

    Type(listechaine), dimension(:), intent(in)       :: lch
    real (kind = 8), dimension (:), intent(inout)  :: u
    real (kind = 8), dimension (:), intent(in)  :: b
    real (kind = 8), intent(in) :: eps
    integer, intent(in) :: niter,taille
    
    integer ::i,it,k
    real (kind = 8) :: r
    
    integer, dimension (:,:), allocatable  :: Aij
    real (kind = 8), dimension(:), allocatable  :: Aval
    integer, dimension (:), allocatable :: ligne
    


    real (kind = 8), dimension (:), allocatable  :: u_par,b_par    
    real (kind = 8), dimension (:), allocatable  :: u_inv      
    real (kind = 4),dimension(:),allocatable  :: Aval4
    
    allocate (u_par(size(u)),u_inv(size(u)),b_par(size(u)))
 
    allocate (Aval4(taille))

    
   
          allocate(Aij(taille,2),Aval(taille),ligne(size(u)+1))


!$OMP PARALLEL      

!$OMP DO
      do i=1,taille
         Aval(i)=0.d0
         Aij(i,1)=0
         Aij(i,2)=0
      end do
!$OMP END DO

!$OMP DO
      do i=1,size(u)+1
         ligne(i)=0
      end do
!$OMP END DO

!$OMP DO
      do i=1,size(u)
         u_par(i)=u(i)
      end do
!$OMP END DO

!$OMP END PARALLEL


    call conversionchainemorseter(lch,Aij,Aval,ligne)

!$OMP PARALLEL
!$OMP DO
       do k=1,size(u)
          u_inv(k)=1.d0
          do i=ligne(k),ligne(k+1)-1
             if(Aij(i,1).eq.k) then
                u_inv(k)=1.d0/Aval(i)
                exit
             endif
          enddo
	  b_par(k)=b(k)*u_inv(k)
          do i=ligne(k),ligne(k+1)-1
             Aval(i)=Aval(i)*u_inv(k)
             Aval4(i)=Aval(i)
          enddo
       end do
!$OMP END DO

!$OMP END PARALLEL

!print*,'apres 0'

    call bicgstab_par(Aij,Aval,Aval4,ligne,u_par,b_par,niter,eps)
!print*,'apres'
    u=u_par
    deallocate (u_par,u_inv,b_par,Aij,Aval,Aval4,ligne)

  end subroutine solve_par

!================================================================================================================


!================================================================================================================!================================================================================================================
 subroutine solve_cg_par(lch,u,b,niter,eps,taille)
!$USE OMP_LIB
use deftype
    implicit none

    Type(listechaine), dimension(:), intent(in)       :: lch
    real (kind = 8), dimension (:), intent(inout)  :: u
    real (kind = 8), dimension (:), intent(in)  :: b
    real (kind = 8), intent(in) :: eps
    integer, intent(in) :: niter,taille
    
    integer ::i,it,k
    real (kind = 8) :: r
    
    integer, dimension (:,:), allocatable  :: Aij
    real (kind = 8), dimension(:), allocatable  :: Aval
    integer, dimension (:), allocatable :: ligne
    


    real (kind = 8), dimension (:), allocatable  :: u_par,b_par    
    real (kind = 8), dimension (:), allocatable  :: u_inv      
    real (kind = 4),dimension(:),allocatable  :: Aval4
    real (kind = 8) ::res
    

    res=10.d0*eps

    allocate (u_par(size(u)),u_inv(size(u)),b_par(size(u)))
 
    allocate (Aval4(taille))

    
   
          allocate(Aij(taille,2),Aval(taille),ligne(size(u)+1))


!$OMP PARALLEL      

!$OMP DO
      do i=1,taille
         Aval(i)=0.d0
         Aij(i,1)=0
         Aij(i,2)=0
      end do
!$OMP END DO

!$OMP DO
      do i=1,size(u)+1
         ligne(i)=0
      end do
!$OMP END DO

!$OMP DO
      do i=1,size(u)
         u_par(i)=u(i)
      end do
!$OMP END DO

!$OMP END PARALLEL


    call conversionchainemorseter(lch,Aij,Aval,ligne)

!$OMP PARALLEL
!$OMP DO
       do k=1,size(u)
          u_inv(k)=1.d0
          do i=ligne(k),ligne(k+1)-1
             if(Aij(i,1).eq.k) then
                u_inv(k)=1.d0/Aval(i)
                exit
             endif
          enddo
	  b_par(k)=b(k)*u_inv(k)
          do i=ligne(k),ligne(k+1)-1
             Aval(i)=Aval(i)*u_inv(k)
             Aval4(i)=Aval(i)
          enddo
       end do
!$OMP END DO

!$OMP END PARALLEL

!print*,'apres 0'
       it=0
       do while (res>eps.and.it<=20)
          it=it+1
          call cg_par(Aij,Aval,Aval4,ligne,u_par,b_par,niter,eps,res)
       end do
!print*,'apres'
    u=u_par
    deallocate (u_par,u_inv,b_par,Aij,Aval,Aval4,ligne)

  end subroutine solve_cg_par

!================================================================================================================

subroutine bicgstab_par ( Aij,Aval,Aval4,ligne,u,b, niterxmax,epsx )
!================================================================================================================

!$use OMP_LIB

        implicit none

        integer, dimension (:,:), intent(in)  :: Aij
        real (kind = 8), dimension(:), intent(in)  :: Aval
        real (kind = 4), dimension(:), intent(in)  :: Aval4
        integer, dimension (:), intent(in)  :: ligne
        real (kind = 8), dimension (:), intent(in) :: b
        real (kind = 8), dimension (:), intent(inout)  :: u
        integer, intent(in) :: niterxmax
        real (kind = 8), intent(in) :: epsx

        real (kind = 8), dimension (size(b)) :: p,q,r,r0,t,y,z
        real (kind = 8) :: teta,gama,omega,rho,rho0,w,w1,w2,gamm,grandeur
        integer :: k, niterx, ll, m, l, nb

        grandeur=maxval(abs(b))
        nb = size(b)
        rho         = 1.d0
        teta        = 1.d0
        omega       = 1.d0
        gamm = 0.985d0

!$OMP PARALLEL

!$OMP DO
        do k=1,nb
           p(k)         = 0.d0
           q(k)         = 0.d0
        end do
!$OMP END DO

!$OMP DO PRIVATE (w)
        do ll=1,nb
           w = 0.
           do l=ligne(ll),ligne(ll+1)-1
              w = w + Aval(l)*u(Aij(l,1))
           enddo
           r(ll) = b(ll) - w
           r0(ll) = r(ll)
        enddo
!$OMP END DO

        do k = 1, niterxmax

!$OMP SINGLE
         niterx    = k
         rho0      = rho
         rho       = 0.
!$OMP END SINGLE

!$OMP DO REDUCTION (+: rho)
         Do l=1,nb
            rho = rho + r0(l) * r(l)
         end Do
!$OMP END DO

!$OMP SINGLE
         gama      = (rho/(rho0+1.d-25)) * (teta/(omega+1.d-25))
!         gama      = (rho/rho0) * (teta/omega)
         w1        = 0.
!$OMP END SINGLE

!$OMP DO
         Do l=1,nb
            p(l)         = r(l) + gama*(p(l)-omega*q(l))
         end Do
!$OMP END DO

!$OMP DO PRIVATE (w)
        do ll=1,nb
           w = 0.
           do l=ligne(ll),ligne(ll+1)-1
              w = w + Aval4(l)*p(Aij(l,1))
           enddo
           y(ll) = (1.+gamm)*p(ll)-gamm*w
        enddo
!$OMP END DO

!$OMP DO PRIVATE (w) REDUCTION (+: w1)
        do ll=1,nb
           w = 0.
           do l=ligne(ll),ligne(ll+1)-1
              w = w + Aval(l)*y(Aij(l,1))
           enddo
           q(ll) = w
           w1 = w1 + r0(ll) * q(ll)
        enddo
!$OMP END DO

!$OMP SINGLE
        teta = rho / (w1+1.d-25)
!        teta = rho / w1
        w1 = 0.
        w2 = 0.
!$OMP END SINGLE

!$OMP DO
        Do l=1,nb
           r(l)         = r(l) - teta*q(l)
        end do
!$OMP END DO

!$OMP DO PRIVATE (w)
        do ll=1,nb
           w = 0.
           do l=ligne(ll),ligne(ll+1)-1
              w = w + Aval4(l)*r(Aij(l,1))
           enddo
           z(ll) = (1.+gamm)*r(ll)-gamm*w
        enddo
!$OMP END DO

!$OMP DO PRIVATE (w) REDUCTION (+: w1, w2)
        do ll=1,nb
           w = 0.
           do l=ligne(ll),ligne(ll+1)-1
              w = w + Aval(l)*z(Aij(l,1))
           enddo
           t(ll) = w
           w1 = w1 + t(ll) * r(ll)
           w2 = w2 + t(ll) * t(ll)
        enddo
!$OMP END DO

!$OMP SINGLE
        omega = w1 / (w2+1.d-25)
        !omega = w1 / w2
        w1 = 0.
!$OMP END SINGLE

!$OMP DO REDUCTION (+: w1)
         Do l=1,nb
            u(l)      = u(l) + teta*y(l) + omega*z(l)
            r(l)   = r(l) - omega*t(l)
            w1 = w1 + r(l) * r(l)
         end do

!$OMP END DO

         if (sqrt(w1/nb) <=epsx*grandeur) exit

! !$OMP SINGLE
!	write(4,*) k,sqrt(w1),sqrt(w1/nb)
! !$OMP END SINGLE

        end do

!$OMP END PARALLEL

        print*,'norme residu iteration ',niterx,' du bicgstab = ',sqrt(w1/nb)

       return

!====================
end subroutine bicgstab_par
!====================
!
!
!
subroutine cg_par ( Aij,Aval,Aval4,ligne,u,b, niterxmax,epsx,res )
!================================================================================================================

!$use OMP_LIB

        implicit none

        integer, dimension (:,:), intent(in)  :: Aij
        real (kind = 8), dimension(:), intent(in)  :: Aval
        real (kind = 4), dimension(:), intent(in)  :: Aval4
        integer, dimension (:), intent(in)  :: ligne
        real (kind = 8), dimension (:), intent(in) :: b
        real (kind = 8), dimension (:), intent(inout)  :: u
        integer, intent(in) :: niterxmax
        real (kind = 8), intent(in) :: epsx
        real (kind = 8), intent(out) :: res

        real (kind = 8), dimension (size(b)) :: p,q,r,y
        real (kind = 8) :: teta,gama,rho,rho0,w,w1,gamm,grandeur
        integer :: k, niterx, ll, m, l, nb

        grandeur=maxval(abs(b))
        nb = size(b)
        rho         = 1.d0
        gamm = 0.97d0
       

!$OMP PARALLEL

!$OMP DO PRIVATE (w)
        do ll=1,nb
           w = 0.
           do l=ligne(ll),ligne(ll+1)-1
              w = w + Aval(l)*u(Aij(l,1))
           enddo
           r(ll) = b(ll) - w
           p(ll) = 0.d0
        enddo
!$OMP END DO

        do k = 1, niterxmax

!$OMP SINGLE
         niterx    = k
         rho0      = rho
         rho       = 0.
!$OMP END SINGLE

!$OMP DO PRIVATE (w) REDUCTION (+: rho)
         Do ll=1,nb
            w = 0.
            do l=ligne(ll),ligne(ll+1)-1
               w = w + Aval4(l)*r(Aij(l,1))
            enddo
            y(ll) = (1.+gamm)*r(ll)-gamm*w
            rho = rho + r(ll) * y(ll)
         end Do
!$OMP END DO

!$OMP SINGLE
         gama      = rho/(rho0+1.d-25)
         w1        = 0.
!$OMP END SINGLE

!$OMP DO
         Do l=1,nb
            p(l)         = y(l) + gama*p(l)
         end Do
!$OMP END DO

!$OMP DO PRIVATE (w) REDUCTION (+: w1)
        do ll=1,nb
           w = 0.
           do l=ligne(ll),ligne(ll+1)-1
              w = w + Aval(l)*p(Aij(l,1))
           enddo
           q(ll) = w
           w1 = w1 + p(ll) * q(ll)
        enddo
!$OMP END DO

!$OMP SINGLE
        teta = rho / (w1+1.d-25)
        w1 = 0.
!$OMP END SINGLE

!$OMP DO REDUCTION (+: w1)
         Do l=1,nb
            u(l)      = u(l) + teta*p(l)
            r(l)      = r(l) - teta*q(l)
            w1 = w1 + r(l) * r(l)
         end do

!$OMP END DO
         res=sqrt(w1/nb)
         if (res <=epsx*grandeur) exit

! !$OMP SINGLE
!  write(9,*) k,sqrt(w1),sqrt(w1/nb)
! !$OMP END SINGLE

        end do

!$OMP END PARALLEL

        print*,'norme residu iteration ',niterx,' du cg = ',res

       return

!====================
end subroutine cg_par
!====================









!======================
end module solveur_openmp
!======================
