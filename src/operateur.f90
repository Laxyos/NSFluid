module operateur

contains



!creation de la matrice de l'etape de prediction
subroutine prediction_u(lch, dx, dy, Nx, Ny, mu, dt, nk)
   use deftype
   implicit none
   type (listechaine), dimension(:), intent(inout) :: lch
   
   real (kind=8), intent(in) 	:: dx, dy, mu, dt
   integer, intent(in) 		:: Nx, Ny
   integer, intent(out)		:: nk
   integer :: i, j, l, C1, C2

!	initialisation liste chainee
   nk = Nx*Ny
   do i=1, nk
       allocate (lch(i)%liste)
        lch(i)%liste%j=i
        lch(i)%liste%val=1.0d0/dt
        nullify(lch(i)%liste%suiv)
        nullify(lch(i)%liste%prec)
    end do

   !parcours des faces verticales
   do j=1,Ny
      do i=1,Nx-1
         C1 = i+(j-1)*Nx
         C2 = C1+1
         call ajoutval(C1, mu*dt/dx**2, lch(C1)%liste, nk)
         call ajoutval(C2, mu*dt/dx**2, lch(C2)%liste, nk)
         call ajoutval(C1, -mu*dt/dx**2, lch(C2)%liste, nk)
         call ajoutval(C2, -mu*dt/dx**2, lch(C1)%liste, nk)
      enddo
    enddo

   !parcours des faces horizontales
   do i=1,Nx
      do j=1,Ny-1
         C1 = i+(j-1)*(Nx)
         C2=C1+Nx
         call ajoutval(C1, mu*dt/dy**2, lch(C1)%liste, nk)
         call ajoutval(C2, mu*dt/dy**2, lch(C2)%liste, nk)
         call ajoutval(C1, -mu*dt/dy**2, lch(C2)%liste, nk)
         call ajoutval(C2, -mu*dt/dy**2, lch(C1)%liste, nk)
      enddo
   enddo

   !parcours des faces horizontales du bas et du haut
   do i=1,Nx
      C1 = i
      C2 = i + (Ny-1)*Nx
      call ajoutval(C1, 2.0d0*mu*dt/dy**2, lch(C1)%liste, nk)
      call ajoutval(C2, 2.0d0*mu*dt/dy**2, lch(C2)%liste, nk)
      !call ajoutval(C1, mu*dt/dy**2, lch(C1)%liste, nk)
      !call ajoutval(C2, mu*dt/dy**2, lch(C2)%liste, nk)
   enddo


   !parcours des faces verticales de gauche et de droite
   !PENALISATION
   do j=1,Ny
      C1 = 1+(j-1)*Nx
      C2 = Nx+(j-1)*Nx
     ! call ajoutval(C1, mu*dt/dx**2, lch(C1)%liste, nk)
     ! call ajoutval(C2, mu*dt/dx**2, lch(C2)%liste, nk)
      call ajoutval(C1, 10.0d0**10.0d0, lch(C1)%liste, nk)
      call ajoutval(C2, 10.0d0**10.0d0, lch(C2)%liste, nk)
   enddo
end subroutine prediction_u



!creation de la matrice de l'etape de prediction
subroutine prediction_v(lch, dx, dy, Nx, Ny, mu, dt, nk)
   use deftype
   implicit none
   type (listechaine), dimension(:), intent(inout) :: lch
   
   real (kind=8), intent(in) 	:: dx, dy, mu, dt
   integer, intent(in) 		:: Nx, Ny
   integer, intent(out)		:: nk
   integer :: i, j, l, C1, C2

!	initialisation liste chainee
   nk = Nx*Ny
   do i=1, nk
       allocate (lch(i)%liste)
        lch(i)%liste%j=i
        lch(i)%liste%val=1.0d0/dt
        nullify(lch(i)%liste%suiv)
        nullify(lch(i)%liste%prec)
    end do

   !parcours des faces verticales
   do j=1,Ny
      do i=1,Nx-1
         C1 = i+(j-1)*Nx
         C2 = C1+1
         call ajoutval(C1, mu*dt/dx**2, lch(C1)%liste, nk)
         call ajoutval(C2, mu*dt/dx**2, lch(C2)%liste, nk)
         call ajoutval(C1, -mu*dt/dx**2, lch(C2)%liste, nk)
         call ajoutval(C2, -mu*dt/dx**2, lch(C1)%liste, nk)
      enddo
    enddo

   !parcours des faces horizontales
   do i=1,Nx
      do j=1,Ny-1
         C1 = i+(j-1)*(Nx)
         C2=C1+Nx
         call ajoutval(C1, mu*dt/dy**2, lch(C1)%liste, nk)
         call ajoutval(C2, mu*dt/dy**2, lch(C2)%liste, nk)
         call ajoutval(C1, -mu*dt/dy**2, lch(C2)%liste, nk)
         call ajoutval(C2, -mu*dt/dy**2, lch(C1)%liste, nk)
      enddo
   enddo

   !parcours des faces horizontales du bas et du haut
   !PENALISATION
   do i=1,Nx
      C1 = i
      C2 = i + (Ny-1)*Nx
      call ajoutval(C1, 10.0d0**10.0d0, lch(C1)%liste, nk)
      call ajoutval(C2, 10.0d0**10.0d0, lch(C2)%liste, nk)
   enddo

   !parcours des faces verticales de gauche et de droite
   do j=1,Ny
      C1 = 1+(j-1)*Nx
      C2 = Nx+(j-1)*Nx
      call ajoutval(C1, 2.0d0*mu*dt/dx**2, lch(C1)%liste, nk)
      call ajoutval(C2, 2.0d0*mu*dt/dx**2, lch(C2)%liste, nk)
   enddo
end subroutine prediction_v


!creation de la matrice de l'etape de projection
subroutine projection(lch, dx, dy, Nx, Ny, nk)
   use deftype
   implicit none
   type (listechaine), dimension(:), intent(inout) :: lch
   
   real (kind=8), intent(in) 	:: dx, dy
   integer, intent(in) 		:: Nx, Ny
   integer, intent(out)		:: nk
   integer :: i, j, l, C1, C2

!	initialisation liste chainee
   nk = Nx*Ny
   do i=1, nk
       allocate (lch(i)%liste)
        lch(i)%liste%j=i
        lch(i)%liste%val=-1.0d-6
        nullify(lch(i)%liste%suiv)
        nullify(lch(i)%liste%prec)
    end do

   !parcours des faces verticales
   do j=1,Ny
      do i=1,Nx-1
         C1 = i+(j-1)*Nx
         C2 = C1+1
         call ajoutval(C1, -1.0d0/dx**2, lch(C1)%liste, nk)
         call ajoutval(C2, -1.0d0/dx**2, lch(C2)%liste, nk)
         call ajoutval(C1,  1.0d0/dx**2, lch(C2)%liste, nk)
         call ajoutval(C2,  1.0d0/dx**2, lch(C1)%liste, nk)
      enddo
    enddo
   !parcours des faces horizontales
   do i=1,Nx
      do j=1,Ny-1
         C1 = i+(j-1)*Nx
         C2=C1+Nx
         call ajoutval(C1, -1.0d0/dy**2, lch(C1)%liste, nk)
         call ajoutval(C2, -1.0d0/dy**2, lch(C2)%liste, nk)
         call ajoutval(C1,  1.0d0/dy**2, lch(C2)%liste, nk)
         call ajoutval(C2,  1.0d0/dy**2, lch(C1)%liste, nk)
      enddo
   enddo


   !!parcours des faces horizontales du bas et du haut
   do i=1,Nx
      C1 = i
      C2 = i + (Ny-1)*Nx
      !call ajoutval(C1, -10.0d0**10.d0, lch(C1)%liste, nk)
      !call ajoutval(C2, -10.0d0**10.d0, lch(C2)%liste, nk)
      !call ajoutval(C1,  -1.0d0/dy**2, lch(C2)%liste, nk)
      !call ajoutval(C2,  -1.0d0/dy**2, lch(C1)%liste, nk)
      !call ajoutval(C2, -1.0d0/dy**2.0d0, lch(C2)%liste, nk)
      !call ajoutval(C1, -1.0d0/dy**2.0d0, lch(C1)%liste, nk)
   enddo


   !parcours des faces verticales de gauche et de droite
   do j=1,Ny
      C1 = 1+(j-1)*Nx
      C2 = Nx+(j-1)*Nx
      !call ajoutval(C1, -1.0d0/dx**2.0d0, lch(C1)%liste, nk)
      !call ajoutval(C2, -1.0d0/dx**2.0d0, lch(C2)%liste, nk)
      !call ajoutval(C1, -10.0d0**10.d0, lch(C1)%liste, nk)
      !call ajoutval(C2, -10.0d0**10.d0, lch(C2)%liste, nk)
      !call ajoutval(C1,  -1.0d0/dx**2, lch(C2)%liste, nk)
      !call ajoutval(C2,  -1.0d0/dx**2, lch(C1)%liste, nk)
   enddo



end subroutine projection

subroutine divergence(u, v, Nx, Ny, dx, dy, div)
   implicit none
   real (kind=8), dimension(:),intent(in), allocatable :: u, v
   real (kind=8), dimension(:), allocatable :: div
   real (kind=8), intent(in)  :: dx, dy
   integer                    :: C1, C2
   integer :: i, j, Nx, Ny, k

   div = 0.0d0

   k=1
   !parcours des faces verticales
   do j=1,Ny-1
      do i=1,Nx-1
         C1 = i+(j-1)*Nx
         C2 = C1+1
         div(k) = div(k) + (u(C2)-u(C1))/dx
         k = k+1
      enddo
    enddo

   !parcours des faces horizontales
   k=1
   do j=1,Ny-1
      do i=1,Nx-1
         C1 = i+(j-1)*(Nx-1)
         C2=C1+Nx-1
         div(k) = div(k) + (v(C2)-v(C1))/dy
         k = k+1
      enddo
   enddo

end subroutine divergence

subroutine divergence2(u, v, Nx, Ny, dx, dy, div)
   implicit none
   real (kind=8), dimension(:),intent(in), allocatable :: u, v
   real (kind=8), dimension(:), allocatable :: div
   real (kind=8), intent(in)  :: dx, dy
   integer                    :: C1, C2
   integer                    :: i, j, Nx, Ny, k

   !u : Nx*(Ny-1)
   !print*, Nx, Ny
   div = 0.0d0

   !print*, "Div horizontal"
   k=1
   do j=1,Ny-1
      do i=1,Nx-1
         C1 = i+(j-1)*Nx
         C2 = C1+1
         !print*, C1, C2, k
         div(k) = div(k) + (u(C2)-u(C1))/dx
         k=k+1
      enddo
   enddo


   !v : Ny*(Nx-1)
   !print*, "Div vertical"
   k=1
   do j=1,Ny-1
      do i=1,Nx-1
         C1 = i+(j-1)*(Nx-1)
         C2 = C1+Nx-1
         !print*, C1, C2, k
         div(k) = div(k) + (v(C2)-v(C1))/dy
         k=k+1
      enddo
      !print*, "---"
   enddo
end subroutine divergence2



!subroutine divergence2(u, v, Nx, Ny, dx, dy, div)
!   implicit none
!   real (kind=8), dimension(:),intent(in), allocatable :: u, v
!   real (kind=8), dimension(:), allocatable :: div
!   real (kind=8), intent(in)  :: dx, dy
!   integer                    :: C1, C2
!   integer                    :: i, j, Nx, Ny, k
!
!   !u : Nx*(Ny-1)
!   !print*, Nx, Ny
!   div = 0.0d0
!
!   k=1
!   do j=1,Ny-1
!      do i=1,Nx-1
!         C1 = i+(j-1)*Nx
!         C2 = C1+1
!         div(k) = div(k) + (u(C2)-u(C1))/dx
!         k=k+1
!      enddo
!   enddo
!
!
!   !v : Ny*(Nx-1)
!   k=1
!   do j=1,Nx-1
!      do i=1,Ny-1
!         C1 = i+(j-1)*(Nx-1)
!         C2 = C1+Nx-1
!         !print*, C1, C2, C1, k
!         div(k) = div(k) + (v(C2)-v(C1))/dy
!         k=k+1
!      enddo
!      !print*, "---"
!   enddo
!end subroutine divergence2

subroutine gradient(p, Nx, Ny, dx, dy, grad_x, grad_y)
   implicit none
   real (kind=8), dimension(:),intent(in), allocatable :: p
   real (kind=8), dimension(:), allocatable :: grad_x, grad_y
   real (kind=8), intent(in)  :: dx, dy
   integer :: i, j, Nx, Ny, C1, C2, C3, k
   
   !initialisation
   grad_x = 0.0d0
   grad_y = 0.0d0

   !gradient selon x
   do j=1, Ny
      do i=1,Nx-1
         C1 = i + (j-1)*Nx
         C2 = C1+1
         grad_x(C2+(j-1)) = (p(C2)-p(C1))/dx
      enddo
   enddo

   !gradient selon y
   do j=1,Ny-1
      do i=1, Nx

         C1 = i + (j-1)*Nx
         C2 = C1+Nx
         grad_y(C2) = (p(C2)-p(C1))/dy
      enddo
   enddo
end subroutine gradient

subroutine moyenne_u_v(u, v, Nx, Ny, Um, Vm)
   implicit none
   real (kind=8), dimension(:),intent(in), allocatable :: u, v
   real (kind=8), dimension(:), allocatable :: Um, Vm
   integer :: i, j, Nx, Ny, C1, C2, C3, k
   
   !initialisation
   Um = 0.0d0
   Vm = 0.0d0

   k=1
   do j=1,Ny-1
      do i=1,Nx-1
         C1 = i+(j-1)*Nx
         C2 = C1+1
         Um(k) = (u(C2)+u(C1))/2.0d0
         k=k+1
      enddo
   enddo
   
   k=1
   do j=1,Ny-1
      do i=1,Nx-1
         C1 = i+(j-1)*(Nx-1)
         C2 = C1+Nx-1
         Vm(k) = (v(C2)+v(C1))/2.0d0
         k=k+1
      enddo
   enddo
end subroutine moyenne_u_v

subroutine vorticite(u, v, Nx, Ny, dx, dy, omega)
   real(kind=8), dimension(:), allocatable, intent(in) :: u, v 
   real(kind=8), dimension(:), allocatable::  omega, omega_x, omega_y
   integer :: i, j, Nx, Ny, k, C1,C2
   real(kind=8), intent(in) :: dx, dy

   omega = 0.0d0

   allocate(omega_x((Nx-1)*Ny), omega_y((Ny-1)*Nx))

   omega_x = 0.0
   omega_y = 0.0

   !derivee de v % x
   k=1
   do j=1,Ny
      do i=1, Nx-1
         C1 = i + (j-1)*(Nx)
         C2 = C1 + 1
         omega_x(k) = (v(C2)-v(C1))/dx
         !omega(k) =  omega(k) + (v(C2)-v(C1))/dx
         k=k+1
      enddo
   enddo 

   !derivee de u % y
   k=1
   do j=1,Ny-1
      do i=1, Nx
         C1 = i + (j-1)*Nx
         C2 = C1 + Nx
         omega_y(k) = (u(C2)-u(C1))/dy
         !omega(k) =  omega(k) + (u(C2)-u(C1))/dy
         k=k+1
      enddo
   enddo 


   k = 1
   do j=1, Ny-1
      do i=1, Nx-1
         C1 = i+(j-1)*(Nx-1)
         C2 = C1 + Nx-1
         omega(k) = omega(k) + 0.5d0*(omega_x(C1) + omega_x(C2))
         k = k+1
      enddo
   enddo

   k = 1
   do j=1, Ny-1
      do i=1, Nx-1
         C1 = i+(j-1)*(Nx-1)
         C2 = C1 + 1
         print*, C1, C2, k
         omega(k) = omega(k) + 0.5d0*(omega_y(C1) + omega_y(C2))
         k = k+1
      enddo
      print*, "---"
   enddo

   print*, omega
   deallocate(omega_x, omega_y)

   
end subroutine vorticite



!subroutine vorticite(u, v, Nx, Ny, dx, dy, omega)
!   real(kind=8), dimension(:), allocatable, intent(in) :: u, v 
!   real(kind=8), dimension(:), allocatable::  omega, omega_x, omega_y
!   integer :: i, j, Nx, Ny, k, C1,C2
!   real(kind=8), intent(in) :: dx, dy
!
!   omega = 0.0d0
!
!   allocate(omega_x((Nx-2)*Ny), omega_y((Ny-2)*Nx))
!
!   omega_x = 0.0
!   omega_y = 0.0
!
!   !derivee de v % x
!   k=1
!   do j=1,Ny
!      do i=1, Nx-2
!         C1 = i + (j-1)*(Nx-1)
!         C2 = C1 + 1
!         omega_x(k) = (v(C2)-v(C1))/dx
!         omega(k) =  omega(k) + (v(C2)-v(C1))/dx
!         k=k+1
!      enddo
!   enddo 
!
!   !derivee de u % y
!   k=1
!   do j=1,Ny-2
!      do i=1, Nx
!         C1 = i + (j-1)*Nx
!         C2 = C1 + Nx
!         omega_y(k) = (u(C2)-u(C1))/dy
!         omega(k) =  omega(k) + (u(C2)-u(C1))/dy
!         k=k+1
!      enddo
!   enddo 
!
!   deallocate(omega_x, omega_y)
!
!   
!end subroutine vorticite




!=======================================================================
  subroutine ajoutval(j,val,ch,t)
    use deftype
    implicit none
    integer, intent(in)               :: j
    Type(chainej), Pointer            :: ch
    real(kind=kind(0.d0)), intent(in) :: val
    integer, intent(inout)            :: t
    logical                           :: droite
    Type(chainej), Pointer            :: temp, loclj 
    !
    !
!print*,'entree ajout val',t,j
       droite=.false.
       if (j>=ch%j) droite=.true.
       do while (associated(ch))
          temp=>ch
          if (j>ch%j) then
             if (droite) then
                ch=>ch%suiv
             else
                t=t+1
                call insertdroitelj(j,val,ch)
                return
             end if
          else
             if (j==ch%j) then 
                ch%val=ch%val+val
                !             print*,'ajout au point declare',t,j
                return
             else
                if (droite) then 
                   t=t+1
                   call insertgauchelj(j,val,ch)
                   return
                else
                   ch=>ch%prec
                end if
             end if
          end if
       end do
       !cas du bout de liste droite ou gauche
       allocate (loclj)
       t=t+1
       loclj%j=j
       loclj%val=val
       if (droite) then
!ajout a droite de temp
          loclj%prec=>temp
          nullify(loclj%suiv)
          temp%suiv=>loclj
       else
!ajout a gauche de temp
          loclj%suiv=>temp
          nullify(loclj%prec)
          temp%prec=>loclj
       end if
       ch=>temp
       return
  end subroutine ajoutval
!=======================================================================

!=======================================================================
 subroutine insertgauchelj(j,val,ch)
  use deftype
    implicit none
    integer, intent(in)               :: j
    Type(chainej), Pointer            :: ch
    real(kind=kind(0.d0)), intent(in) :: val
    Type(chainej), Pointer            :: temp, loclj
    !
!cas du bout de liste droite ou gauche
    temp=>ch%prec
       allocate (loclj)
       loclj%j=j
       loclj%val=val
       loclj%suiv=>ch
       ch%prec=>loclj
       loclj%prec=>temp
       temp%suiv=>loclj
       return
  end subroutine insertgauchelj
!=======================================================================

!=======================================================================
 subroutine insertdroitelj(j,val,ch)
  use deftype
    implicit none
    integer, intent(in)               :: j
    Type(chainej), Pointer            :: ch
    real(kind=kind(0.d0)), intent(in) :: val
    Type(chainej), Pointer            :: temp, loclj
    !
!cas du bout de liste droite ou gauche
    temp=>ch%suiv
       allocate (loclj)
       loclj%j=j
       loclj%val=val
       loclj%suiv=>temp
       temp%prec=>loclj
       loclj%prec=>ch
       ch%suiv=>loclj
       return
  end subroutine insertdroitelj
!=======================================================================


!=======================================================================
 subroutine test_symetrie(lch)
    use deftype
    implicit none
    Type(listechaine), dimension(:), intent(in)       :: lch

! 
    type (chainej), Pointer            :: chj,chjb
    integer :: i,k,tailleA
    real (kind=8) :: delta
    print*,'test symetrie'
  
    do k=1,size(lch)
       chj=>lch(k)%liste
       do while (associated(chj%prec))
          chj=>chj%prec
       end do
       do while (associated(chj))
          chjb=>lch(chj%j)%liste
          do while (associated(chjb%prec))
             chjb=>chjb%prec
          end do
          do while (associated(chjb))
             if (chjb%j==k) then
                delta=chjb%val-chj%val
!print*,'i,j,jb,val,valj',k,chj%j,chjb%j,chj%val,chjb%val
                if (abs(delta)>1.d-12) print*,'non symetrie'
             end if
             chjb=>chjb%suiv
          end do
         
          chj=>chj%suiv
       end do
    end do
  
  end subroutine test_symetrie
!=======================================================================
!=======================================================================
 subroutine produit(lch,u,v)
    use deftype
    implicit none
    Type(listechaine), dimension(:), intent(in)       :: lch
    real(kind=8),dimension(:), intent(in )  :: u
    real(kind=8),dimension(:), intent(out)  :: v
! 
    type (chainej), Pointer            :: chj
    integer :: i,k,tailleA
    
    v=0.d0
    do i=1,size(lch)
       chj=>lch(i)%liste
       do while (associated(chj%prec))
          chj=>chj%prec
       end do
       do while (associated(chj))
          v(i)=v(i)+chj%val*u(chj%j)
          chj=>chj%suiv
       end do
    end do
  end subroutine produit
!=======================================================================
!=======================================================================
 subroutine conversionchainemorsebis(lch,A)
    use deftype
    implicit none
    Type(listechaine), dimension(:), intent(in)       :: lch
    type(morse),dimension(:), intent(out)  :: A
! 
    type (chainej), Pointer            :: chj
    integer :: i,k,tailleA
    
    k=0
    do i=1,size(lch)
       chj=>lch(i)%liste
       do while (associated(chj%prec))
          chj=>chj%prec
       end do
       do while (associated(chj))
          k=k+1
   !       print*,'k',k,chj%j
          A(k)%i=i
          A(k)%j=chj%j
          A(k)%val=chj%val
          chj=>chj%suiv
       end do
    end do
tailleA=k
!print*,'verification de la taille', tailleA
end subroutine conversionchainemorsebis
!=======================================================================

!=======================================================================
 subroutine printChaine2(lch, N)
    use deftype
    implicit none
    Type(listechaine), dimension(:), intent(in)       :: lch
    real (kind=8), dimension(:,:), allocatable :: A
! 
    type (chainej), Pointer            :: chj
    integer :: i,j,k, N
    
    allocate(A(N, N))
    DO i=1,N
         DO j=1,N
            A(i, j) = 0.0d0
         ENDDO
    ENDDO


    k=0
    do i=1,size(lch)
       chj=>lch(i)%liste
       do while (associated(chj%prec))
          chj=>chj%prec
       end do
       do while (associated(chj))
          k=k+1
          !print*,i,chj%j, chj%val
          A(i, chj%j) = chj%val
          chj=>chj%suiv
       end do
    end do

    DO i=1,N
      print*, A(i, 1:N)
    ENDDO

end subroutine printChaine2
!=======================================================================
 subroutine desallouelistebis(ch)
    use deftype
    implicit none
    Type(listechaine), dimension(:), intent(inout)         :: ch
! 
    type (chainej), Pointer            :: chj,chjb,chjbb,chjbef
    integer :: i
    
! deallocation 
! deallocation partie droite de sous chaine
    do i=1,size(ch)
       chj=>ch(i)%liste

       chjb=>chj
       chjbb=>chj%prec
       do while (associated(chjb))
          chjbef=>chjb
          chjb=>chjb%suiv
          deallocate(chjbef)
       end do
       do while (associated(chjbb))
          chjbef=>chjbb
          chjbb=>chjbb%prec
          deallocate(chjbef)
       end do


    end do
!
  end subroutine desallouelistebis
!=======================================================================







end module operateur
