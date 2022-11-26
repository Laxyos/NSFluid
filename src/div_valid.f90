program div_valid 
    use deftype
    use operateur
    use solveur_openmp
    use sortie
    
!-------------------------------------------------------------------------
!                   VALIDATION DE LA SUBROUTINE divergence
!-------------------------------------------------------------------------

    implicit none
    real(kind=8), dimension(:), allocatable :: u, v, div_ex, div
    real(kind=8)    :: Lx_domain,Ly_domain, dx, dy
    integer :: N, Nx, Ny, i, j

    N=40
    Nx = 50
    Ny = 50
    Lx_domain = 1.0d0
    Ly_domain = 1.0d0
    dx = Lx_domain/dble(Nx-1)
    dy = Ly_domain/dble(Ny-1)
    

    allocate(u(Nx*(Ny-1)), v(Ny*(Nx-1)))
    allocate(div_ex((Nx-1)*(Ny-1)), div((Nx-1)*(Ny-1)))

    u = 0.0d0
    v = 0.0d0
    div_ex = 0.0d0
    div = 0.0d0

    do i=1,Nx
        do j=1, Ny-1
            print*, i+(j-1)*Nx
            u(i+(j-1)*Nx) = ((i-1.0d0)*dx)**2.0d0 * (j-0.5d0)*dy
        enddo
    enddo

    do i=1,Nx-1
        do j=1,Ny 
            v(i+(j-1)*(Nx-1)) = -2.0d0*(i-0.5d0)*dx * (j-1.0d0)*dy
        enddo
    enddo

    do i=1,Nx-1
        do j=1, Ny-1
            div_ex(i+(j-1)*(Nx-1)) = 2*((i-0.5d0)*dx)*((j-0.5d0)*dy-1.0d0)
            !div_ex(i+(j-1)*(Nx-1)) = 2.0d0*(i-0.5d0)*dx*(j-0.5d0)*dy
            !div_ex(i+(j-1)*(Nx-1)) = -2.0d0*(i-0.5d0)*dx
        enddo
    enddo

    print*, div_ex
    call divergence2(u, v, Nx, Ny, dx, dy, div)

    call imprim2d (0, Nx,Ny-1,u,u,u,u,dx,dy)
    call imprim2d (1, Nx-1,Ny, v,v,v,v,dx,dy)
    call imprim2d (2, Nx-1,Ny-1, div_ex,div_ex,div_ex,div_ex,dx,dy)
    call imprim2d (3, Nx-1,Ny-1, div,div,div,div,dx,dy)
    
end program div_valid
    