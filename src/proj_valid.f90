program proj_valid 
    use deftype
    use operateur
    use solveur_openmp
    use sortie
    
!-------------------------------------------------------------------------
!                   VALIDATION DE LA SUBROUTINE gradient
!-------------------------------------------------------------------------

    implicit none
    real(kind=8), dimension(:), allocatable :: u, uex, f 
    type (listechaine), dimension(:), allocatable :: mat_u
    real(kind=8)    :: L_domain, dx, dy, PI, eps, temp1, temp2, temp3, temp4, sigma
    integer :: N, Nx, Ny, i, j, nk, niter
    real(kind=8) :: x, y

    PI=4.D0*DATAN(1.D0)
    

    N=60
    Nx = N
    Ny = N
    L_domain = 1.0d0
    dx = L_domain/dble(Nx-1)
    dy = L_domain/dble(Ny-1)

    eps=1.d-12
    niter=3*max(Nx, Ny)

    allocate(mat_u((Nx-1)*(Ny-1)))
    allocate(u((Nx-1)*(Ny-1)), uex((Nx-1)*(Ny-1)), f((Nx-1)*(Ny-1)))

    u =   0.0d0
    uex =   0.0d0
    f = 0.0d0

    sigma = 0.1d0

    do i=1, Nx-1
        do j=1, Ny-1
            x = (i-0.5d0)*dx
            y = (j-0.5d0)*dy
            !uex(i+(j-1)*(Nx-1)) = dsin(PI*(i-0.5d0)*dx)*dsin(PI*(j-0.5d0)*dy)
            temp1 = dexp(-0.5d0*((x- 0.5d0)/sigma)**2.0d0)/(dsqrt(2.0d0*PI)*sigma)
            temp2 = dexp(-0.5d0*((y - 0.5d0)/sigma)**2.0d0)/(dsqrt(2.0d0*PI)*sigma)
            !temp1 = 1.0d0/(1.0d0+dexp(-(i-0.5d0)*dx))
            !temp2 = 1.0d0/(1.0d0+dexp(-(j-0.5d0)*dy))
            !uex(i+(j-1)*(Nx-1)) = temp1*temp2
            !f(i+(j-1)*(Nx-1)) = -2.0d0*PI*PI*dsin(PI*(i-0.5d0)*dx)*dsin(PI*(j-0.5d0)*dy)
            uex(i+(j-1)*(Nx-1)) = 2*dexp(-0.5d0*((x-0.5d0)**2.0d0+(y-0.5d0)**2.0d0)/sigma**2.0d0)

            temp3 = (-sigma**2.0d0+x**2.0d0-x+0.25d0)*uex(i+(j-1)*(Nx-1))/(sigma**4.0d0)
            temp4 = (-sigma**2.0d0+y**2.0d0-y+0.25d0)*uex(i+(j-1)*(Nx-1))/(sigma**4.0d0)
            !temp1 = temp3*temp1
            !temp1 = dexp(-0.5d0*((x-0.5d0)/sigma)**2.0d0)*temp3/(sigma**5.0d0*dsqrt(2.0d0*PI))
            !temp2 = dexp(-0.5d0*((y-0.5d0)/sigma)**2.0d0)*temp4/(sigma**5.0d0*dsqrt(2.0d0*PI))
            f(i+(j-1)*(Nx-1)) = temp3+temp4
        enddo
    enddo

    !f(1) = 0.0d0
    !f(9) = 0.0d0
    !f(2) = 0.0d0
    !f(3) = 0.0d0
    !f(4) = 0.0d0

    !do i=1, Nx-1
    !    f(i) = 0.0d0
    !    f(i+(Ny-2)*(Nx-1)) = 0.0d0
    !enddo

    !do i=1, Ny-1
    !    f(1+(i-1)*(Nx-1)) = 0.0d0
    !    f(Nx-1+(i-1)*(Nx-1)) = 0.0d0
    !    print*,  1+(i-1)*(Nx-1), Nx-1+(i-1)*(Nx-1)
    !enddo


    

    !do i=1,Nx-1
    !    f(i+(Ny-2)*(Nx-1)) = f(i+(Ny-2)*(Nx-1)) + 2.0d0/dy**2.0d0
    !enddo


    call projection(mat_u, dx, dy, Nx-1, Ny-1, nk)
    !call printChaine2(mat_u, (Nx-1)*(Ny-1))

    !call solve_par(mat_u,u,f,niter,eps,nk)
    call solve_par(mat_u,u,f,niter,eps,nk)

    call imprim2d (0, Nx-1,Ny-1,u,u,u,u,dx,dy)
    call imprim2d (1, Nx-1,Ny-1,uex,uex,uex,uex,dx,dy)
    call imprim2d (2, Nx-1,Ny-1,f,f,f,f,dx,dy)
    
    u = dabs(u-uex)
    print*, maxval(u)

        
end program proj_valid
    