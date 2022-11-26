!===========
program main 
!===========
      use deftype
      use operateur
      use solveur_openmp
      use sortie

    implicit none

    real (kind = 8), dimension (:), allocatable :: b_u, b_v, u, v, p, uex
    real (kind = 8), dimension (:), allocatable :: div, grad_x, grad_y
    type (listechaine), dimension(:), allocatable :: mat_ch
    type (listechaine), dimension(:), allocatable :: mat_u, mat_v, mat_p
    real (kind=8), dimension(:,:), allocatable :: Umat

    integer                     :: N,Nx,Ny,i, j, k,niter,nk, l
    real (kind=8) 				:: t, tmax
    real ( kind = 8)            :: epsdiv,eps,dx, dy, dt, L_domain, mu, PI, E


    PI=4.D0*DATAN(1.D0)

    k=0



    N=100
    Nx = N
    Ny = N
    L_domain = 1.0d0
    dx = L_domain/dble(Nx-1)
    dy = L_domain/dble(Ny-1)
    dt = 0.001d0
    tmax = 0.40d0
    t=0.0d0
    !mu = 0.1d0

    !VALIDATION
    mu = 1.0d0
    dt = 1.0d0



    allocate(mat_u(Nx*(Ny-1)))
    allocate(mat_v(Ny*(Nx-1)))
    allocate(mat_p((Nx-1)*(Ny-1)))
    allocate(u(Nx*(Ny-1)), v(Ny*(Nx-1)), uex(Nx*(Ny-1)))
    allocate(p(size(mat_p)),b_u(size(mat_u)), b_v(size(mat_v)))
    allocate(div((Nx-1)*(Ny-1)), grad_x(Nx*(Ny-1)), grad_y((Nx-1)*Ny))


    u = 0.0d0
    do i=1,Nx
        u(i+(Ny-2)*Nx) = 1.0d0
    enddo
    v = 0.0d0
    p = 0.0d0
    

    DO WHILE (t.le.tmax)

        call imprim2d (k, Nx,Ny-1,u,u,u,u,dx,dy)
        !call imprim2d (k, Nx-1,Ny,v,v,v,v,dx,dy)
        !call imprim2d (k, Nx-1,Ny-1,p,p,p,p,dx,dy)
        k=k+1

    !	Calcul de la prediction 
        eps=1.d-12
        niter=3*max(Nx, Ny)

        
      
        b_v =0.0d0
        b_u = 0.0d0

        !do i=1,Nx
        !    b_u(i+(Ny-2)*Nx) = 2.0d0/dy
        !enddo

        !VALIDATION PREDICTION
        do i=1, Nx
            do j=1, Ny-1
                b_u(i+(j-1)*Nx) = dsin(PI*(i-1)*dx)*dsin(PI*(j-0.5d0)*dy)*(1.0D0+2.0D0*PI*PI)
            enddo
        enddo



        call prediction_u(mat_u, dx, dy, Nx, Ny-1, mu, dt, nk)
        call solve_par(mat_u,u,b_u,niter,eps,nk)
        call imprim2d (k, Nx,Ny-1,u,u,u,u,dx,dy)

        !do j=1,Ny-1
        !    u(1+(j-1)*Nx) = 0.0d0
        !    u(Nx+(j-1)*Nx) = 0.0d0
        !enddo


        !SOLUTION EXACTE
        do i=1, Nx
            do j=1, Ny-1
                uex(i+(j-1)*Nx) = dsin(PI*(i-1)*dx)*dsin(PI*(j-0.5d0)*dy)
            enddo
        enddo

        E = 0.d0
        do i=1, Nx
            do j=1, Ny-1
                E = E + (u(i+(j-1)*Nx)-uex(i+(j-1)*Nx))**2.0d0
            enddo
        enddo

        E = dsqrt(E)/dble(N)
        u = dabs(u-uex)
        print*, maxval(u)


        !call imprim2d (k+1, Nx,Ny-1,u,u,u,u,dx,dy)
        call imprim2d (k+1, Nx,Ny-1,uex,uex,uex,uex,dx,dy)
        
        
        stop




        call prediction_v(mat_v, dx, dy, Nx-1, Ny, mu, dt, nk)
        call solve_par(mat_v,v,b_v,niter,eps,nk)



        !call printChaine2(mat_u, Nx*(Ny-1))
        !print*, "---------------------------------------------"
        !call printChaine2(mat_v, (Nx-1)*Ny)


        call projection(mat_p, dx, dy, Nx-1, Ny-1, nk)
        !call printChaine2(mat_p, (Nx-1)*(Ny-1))


        call divergence(u, v, Nx, Ny, dx, dy, div)
        call solve_par(mat_p,p,div/dt,niter,eps,nk)
        !print*, p
        call gradient(p, Nx, Ny, dx, dy, grad_x, grad_y)


        u = u - dt*grad_x 
        v = v - dt*grad_y

        !u = u - dt*grad_x
        !v = v - dt*grad_y
        
        
        !do i=1,Nx-1
        !    v(i + (Ny-1)*(Nx-1)) = 0.0d0
        !    v(i) = 0.0d0
        !enddo

        !do j=1,Ny-1
        !    u(1+(j-1)*Nx) = 0.0d0
        !    u(Nx+(j-1)*Nx) = 0.0d0
        !enddo

    !	print*, p

        print*, "temps:", t
        t = t + dt

    ENDDO
    
    !call imprim2d (0, Nx,Ny-1,u,u,u,u,dx,dy)



    deallocate(mat_u, mat_v, mat_p, p,b_u, b_v, grad_x, grad_y, div)
    stop

    
!
!===============
end program main
!===============
