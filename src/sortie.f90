module sortie


contains

  !=================================================================
subroutine imprim2d (nt, Nx,Ny,phi,p,q,r,dx,dy)
!=================================================================
!

!
        implicit none
!
        real (kind = 8), dimension (1:), intent(in) :: phi,p,q,r
        integer, intent(in) :: nt, Nx,Ny
        real (kind=8), intent(in) :: dx,dy
        
        integer :: i, j, l, k, m, n
        character (len=49) :: str
!

        write (str(1:),1) 'N'
        k = len(trim(str(1:)))
        write (str(k+1:),2) nt
        m = len(trim(str(k+1:)))
        do i=1, m
         if (str(k+i:k+i)==' ') str(k+i:k+i)='0'
        end do

!
        open (unit = 11, file = './SOLU'//'TIO'//trim(str)//'.plt', status = 'unknown')
        WRITE (11, *) 'TITLE = "PHI"'
        WRITE (11, *) 'VARIABLES = "X", "Y", "p","u","v","r"'
        WRITE (11, *) 'ZONE T="FILET", I=',Nx,', J=',Ny,',F=POINT'
        WRITE (11,'(6E15.7)') (( (i-0.5)*dx, (j-0.5)*dy, phi(i+(j-1)*Nx),&
             p(i+(j-1)*Nx),q(i+(j-1)*Nx),r(i+(j-1)*Nx), i=1, Nx), j=1, Ny)

        close(11)
!
   1    format(a1)
   2    format(i5)
   3    format(a4)
!
!====================
end subroutine imprim2d
!====================


  !=================================================================
subroutine imprim2d2 (nt, Nx,Ny,phi,p,q,r,dx,dy)
!=================================================================
!

!
        implicit none
!
        real (kind = 8), dimension (1:), intent(in) :: phi,p,q,r
        integer, intent(in) :: nt, Nx,Ny
        real (kind=8), intent(in) :: dx,dy
        
        integer :: i, j, l, k, m, n
        character (len=49) :: str
!

        write (str(1:),1) 'N'
        k = len(trim(str(1:)))
        write (str(k+1:),2) nt
        m = len(trim(str(k+1:)))
        do i=1, m
         if (str(k+i:k+i)==' ') str(k+i:k+i)='0'
        end do

!
        open (unit = 11, file = './SOL'//'TIO'//trim(str)//'.plt', status = 'unknown')
        WRITE (11, *) 'TITLE = "PHI"'
        WRITE (11, *) 'VARIABLES = "X", "Y", "PHI","p","q","r"'
        WRITE (11, *) 'ZONE T="FILET", I=',Nx,', J=',Ny,',F=POINT'
        WRITE (11,'(6E15.7)') (( (i-0.5)*dx, (j-0.5)*dy, phi(i+(j-1)*Nx),&
             p(i+(j-1)*Nx),q(i+(j-1)*Nx),r(i+(j-1)*Nx), i=1, Nx), j=1, Ny)

        close(11)
!
   1    format(a1)
   2    format(i5)
   3    format(a4)
!
!====================
end subroutine imprim2d2
!====================
  !=================================================================
subroutine imprim3d (nt, Nx,Ny,Nz,phi,u,v,w,dx,dy,dz)
!=================================================================
!

!
        implicit none
!
        real (kind = 8), dimension (1:), intent(in) :: phi,u,v,w
        integer, intent(in) :: nt, Nx,Ny,Nz
        real (kind=8), intent(in) :: dx,dy,dz
        
        integer :: i, j, k, l, m, n
        character (len=49) :: str
!

        write (str(1:),1) 'N'
        k = len(trim(str(1:)))
        write (str(k+1:),2) nt
        m = len(trim(str(k+1:)))
        do i=1, m
         if (str(k+i:k+i)==' ') str(k+i:k+i)='0'
        end do

!
        open (unit = 11, file = './SOLU'//'TIO'//trim(str)//'.plt', status = 'unknown')
        WRITE (11, *) 'TITLE = "PHI"'
        WRITE (11, *) 'VARIABLES = "X", "Y", "Z", "PHI","u","v","w"'
        WRITE (11, *) 'ZONE T="FILET", I=',Nx,', J=',Ny,', K=',Nz,',F=POINT'
        WRITE (11,'(7E15.7)') ((( (i-0.5)*dx, (j-0.5)*dy,  (k-0.5)*dz,phi(i+(j-1)*Nx+(k-1)*Nx*Ny),&
             u(i+(j-1)*Nx+(k-1)*Nx*Ny),v(i+(j-1)*Nx+(k-1)*Nx*Ny),w(i+(j-1)*Nx+(k-1)*Nx*Ny), i=1, Nx)&
             , j=1, Ny), k=1, Nz)

        close(11)
!
   1    format(a1)
   2    format(i5)
   3    format(a4)
!
!====================
end subroutine imprim3d
!====================

subroutine printChaine(lch)
   use deftype
   implicit none
   Type(listechaine), dimension(:), intent(in)       :: lch
! 
   type (chainej), Pointer            :: chj
   integer :: i,k
   
   k=0
   do i=1,size(lch)
      chj=>lch(i)%liste
      do while (associated(chj%prec))
         chj=>chj%prec
      end do
      do while (associated(chj))
         print*,chj%val
         chj=>chj%suiv
      end do
   end do
!print*,'verification de la taille', tailleA
end subroutine printChaine

end module sortie


