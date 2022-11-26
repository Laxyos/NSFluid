program grad_valid 
    use deftype
    use operateur
    use solveur_openmp
    use sortie
    
!-------------------------------------------------------------------------
!                   VALIDATION DE LA SUBROUTINE gradient
!-------------------------------------------------------------------------

    implicit none
    real(kind=8), dimension(:), allocatable :: u, dux, duy, grad_x, grad_y
    real(kind=8)    :: Lx_domain, Ly_domain, dx, dy
    integer :: N, Nx, Ny, i, j

    N=40
    Nx = 40
    Ny = 80
    Lx_domain = 0.5d0
    Ly_domain = 1.0d0
    dx = Lx_domain/dble(Nx-1)
    dy = Ly_domain/dble(Ny-1)

    allocate(u((Nx-1)*(Ny-1)))
    allocate(dux(size(u)), duy(size(u)))
    allocate(grad_x(Nx*(Ny-1)), grad_y(Ny*(Nx-1)))


    u =   0.0d0
    dux = 0.0d0
    duy = 0.0d0

    do i=1, Nx-1
        do j=1, Ny-1
            u(i + (j-1)*(Nx-1)) = ((i-0.5d0)*dx)**2.0d0 * (j-0.5d0)*dy
            dux(i + (j-1)*(Nx-1)) = 2*(i-0.5d0)*dx*(j-0.5d0)*dy
            duy(i + (j-1)*(Nx-1)) = ((i-0.5d0)*dx)**2.0d0
        enddo
    enddo

    call gradient(u, Nx-1, Ny-1, dx, dy, grad_x, grad_y)

    call imprim2d (0, Nx-1,Ny-1,u,u,u,u,dx,dy)

    call imprim2d (1, Nx-1,Ny-1,dux,dux,dux,dux,dx,dy)
    call imprim2d (2, Nx,Ny-1,grad_x,grad_x,grad_x,grad_x,dx,dy)

    call imprim2d (3, Nx-1,Ny-1,duy,duy,duy,duy,dx,dy)
    call imprim2d (4, Nx-1,Ny,grad_y,grad_y,grad_y,grad_y,dx,dy)
    print*, grad_y
    

        
    
        
end program grad_valid
    