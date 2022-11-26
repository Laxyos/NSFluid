!===========
program main 
!===========
      use deftype
      use operateur
      use solveur_openmp
      use sortie

    implicit none

    real (kind = 8), dimension (:), allocatable :: b_u, b_v, u, v, p, Um, Vm, omega
    real (kind = 8), dimension (:), allocatable :: div, grad_x, grad_y
    type (listechaine), dimension(:), allocatable :: mat_u, mat_v, mat_p
    real (kind=8), dimension(:,:), allocatable :: Umat

    integer                     :: N,Nx,Ny,i, j, k,niter,nk, l, it, idx
    real (kind=8) 				:: t, tmax
    real ( kind = 8)            :: epsdiv,eps,dx, dy, dt, Lx_domain,Ly_domain, mu

    k=0

    N=40
    Nx = N
    Ny = N
    Lx_domain = 1.0d0
    Ly_domain = 1.0d0
    dx = Lx_domain/dble(Nx-1)
    dy = Ly_domain/dble(Ny-1)
    dt = 0.01d0
    tmax = 10.0d0
    t=0.0d0
    mu = 1.0d0
    it = 0


    eps=1.d-12
    niter=3*max(Nx, Ny)

    allocate(mat_u(Nx*(Ny-1)))
    print*, size(mat_u)

    allocate(mat_v(Ny*(Nx-1)))
    allocate(mat_p((Nx-1)*(Ny-1)))
    allocate(Um(size(mat_p)), Vm(size(mat_p)), omega((Nx-2)*(Ny-2)))
    allocate(u(Nx*(Ny-1)), v(Ny*(Nx-1)))
    allocate(p(size(mat_p)),b_u(size(mat_u)), b_v(size(mat_v)))
    allocate(div((Nx-1)*(Ny-1)), grad_x(Nx*(Ny-1)), grad_y((Nx-1)*Ny))


    u = 0.0d0
    do i=1,Nx
        u(i+(Ny-2)*Nx) = 1.0d0
        !u(i) = -1.0d0
    enddo
    v = 0.0d0
    !do i=1, Ny
    !    v(1 + (i-1)*(Nx-1)) = 1.0d0
    !    v(Nx-1 + (i-1)*(Nx-1)) = -1.0d0
    !enddo
    p = 0.0d0

    open(32, file="min_u_vel")

    call moyenne_u_v(u, v, Nx, Ny, Um, Vm)
    !call vorticite(u, v, Nx-1, Ny-1, dx, dy, omega)
    call imprim2d (k, Nx-1,Ny-1,p,Um,Vm,dsqrt(Um**2.0d0 + Vm**2.0d0),dx,dy)
    !call imprim2d (k, Nx-2,Ny-2,omega, omega, omega, omega,dx,dy)

    DO WHILE (t.le.tmax)

        !call imprim2d (k, Nx-1,Ny,v,v,v,v,dx,dy)
        !call imprim2d (k, Nx-1,Ny-1,p,p,p,p,dx,dy)
        

    !	Calcul de la prediction 

        b_u = u/dt
        b_v = v/dt

        do i=1,Nx
            b_u(i+(Ny-2)*Nx) = b_u(i+(Ny-2)*Nx) + 2.0d0*mu*dt/dy**2.0d0
            !b_u(i) = b_u(i) - 2.0d0*mu*dt/dy**2.0d0
        enddo

        !do i=1, Nx-1
        !    b_v(1 + (i-1)*(Nx-1)) = b_v(1 + (i-1)*(Nx-1)) + 2.0d0*mu*dt/dx**2.0d0
        !    b_v(Nx-1 + (i-1)*(Nx-1)) = b_v(Nx-1 + (i-1)*(Nx-1)) - 2.0d0*mu*dt/dx**2.0d0
        !enddo


! ----------------- ETAPE DE PREDICTION -----------------------
        call prediction_u(mat_u, dx, dy, Nx, Ny-1, mu, dt, nk)
        call solve_par(mat_u,u,b_u,niter,eps,nk)

        call prediction_v(mat_v, dx, dy, Nx-1, Ny, mu, dt, nk)
        call solve_par(mat_v,v,b_v,niter,eps,nk)
!---------------------------------------------------------------

! ----------------- ETAPE DE PROJECTION -----------------------
        call projection(mat_p, dx, dy, Nx-1, Ny-1, nk)
        call divergence2(u, v, Nx, Ny, dx, dy, div)
        call solve_par(mat_p,p,div/dt,niter,eps,nk)

        call gradient(p, Nx-1, Ny-1, dx, dy, grad_x, grad_y)

        u = u - dt*grad_x 
        v = v - dt*grad_y
!---------------------------------------------------------------

        call moyenne_u_v(u, v, Nx, Ny, Um, Vm)
        !call vorticite(Um, Vm, Nx-1, Ny-1, dx, dy, omega)

        if(MOD(it, 10) == 0) then
            write(32, *) t, minval(u)
            k=k+1
            !call imprim2d (k, Nx-2,Ny-2,omega, omega, omega, omega,dx,dy)
            call imprim2d (k, Nx-1,Ny-1,p,Um,Vm,dsqrt(Um**2.0d0+Vm**2.0d0),dx,dy)
        endif
        print*, "temps:", t, "iteration:", it
        t = t + dt
        it = it + 1
    ENDDO
    

    open(33, file="vy_mid_prof")
    do i=1, Nx-1
        write(33, *) (i-0.5d0)*dx, v(i+(int(Ny/2)-1)*(Nx-1))
    enddo

    open(33, file="ux_mid_prof")
    do i=1, Ny-1
        write(33, *) (i-0.5d0)*dy, u(int(Nx/2) + (i-1)*Nx)
    enddo

    deallocate(mat_u, mat_v, mat_p, u, v, p,b_u, b_v, grad_x, grad_y, div, omega)
    stop

    
!
!===============
end program main
!===============
