!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use incomp_class,      only: incomp, bcond, dirichlet, clipped_neumann, slip
   use timetracker_class, only: timetracker
   use pgrid_class,       only: pgrid
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   use lpt_class,         only: lpt
   use partmesh_class,    only: partmesh
   implicit none
   private
   
   !> Solver and tracker objects
   type(incomp),      public :: fs
   type(timetracker), public :: time
   type(hypre_str),   public :: ps
   type(hypre_str),   public :: vs
   type(lpt),         public :: lp
   type(partmesh),    public :: pmesh
   
   !> Ensight and monitoring
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   type(monitor) :: mfile, cflfile, forcefile
   
   public :: simulation_init, simulation_run, simulation_final
   
   !> Work arrays
   real(WP), dimension(:,:,:), allocatable :: resU, resV, resW
   real(WP), dimension(:,:,:), allocatable :: Ui, Vi, Wi
   real(WP), dimension(:,:,:,:), allocatable :: SR
   real(WP), dimension(:,:,:), allocatable :: U, V, W, rho, visc
  
   !> Volume fraction arrays
   real(WP), dimension(:,:,:), allocatable :: particle_VF, fluid_VF
  
   !> Particle and flow properties
   real(WP) :: Uin, visc_val
   real(WP) :: mean_y, max_part_vel, dp
   type(bcond), pointer :: inflow_bc, outflow_bc
   type(event) :: ppevt
   integer :: step_count = 0
   integer, parameter :: collision_freq = 5
   integer, parameter :: injection_freq = 5

contains

   !> Boundary condition locators
   function left_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (i == pg%imin) isIn = .true.
   end function left_boundary

   function right_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (i == pg%imax+1) isIn = .true.
   end function right_boundary

   function ymin_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (j == pg%jmin) isIn = .true.
   end function ymin_boundary

   function ymax_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (j == pg%jmax+1) isIn = .true.
   end function ymax_boundary

   function zmin_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (k == pg%kmin) isIn = .true.
   end function zmin_boundary

   function zmax_boundary(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (k == pg%kmax+1) isIn = .true.
   end function zmax_boundary


   !> Calculate average particle position
   subroutine calc_barycenter
      use mpi_f08,  only: MPI_ALLREDUCE, MPI_SUM, MPI_INTEGER
      use parallel, only: MPI_REAL_WP
      implicit none
      integer :: i, ierr
      integer :: my_number_particles, number_particles
      real(WP) :: my_mean_y
      
      my_number_particles = 0
      my_mean_y = 0.0_WP
      do i = 1, lp%np_
         my_number_particles = my_number_particles + 1
         my_mean_y = my_mean_y + lp%p(i)%pos(2)
      end do
      
      call MPI_ALLREDUCE(my_number_particles, number_particles, 1, MPI_INTEGER, MPI_SUM, lp%cfg%comm, ierr)
      call MPI_ALLREDUCE(my_mean_y, mean_y, 1, MPI_REAL_WP, MPI_SUM, lp%cfg%comm, ierr)
      mean_y = mean_y / real(max(number_particles, 1), WP)
   end subroutine calc_barycenter
   
   !> Calculate maximum particle velocity
   subroutine calc_max_part_vel
      use mpi_f08,  only: MPI_ALLREDUCE, MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      integer :: i, ierr
      real(WP) :: my_max_vel, part_speed
      
      my_max_vel = 0.0_WP
      do i = 1, lp%np_
         part_speed = sqrt(lp%p(i)%vel(1)**2 + lp%p(i)%vel(2)**2 + lp%p(i)%vel(3)**2)
         if (part_speed > my_max_vel) my_max_vel = part_speed
      end do
      
      call MPI_ALLREDUCE(my_max_vel, max_part_vel, 1, MPI_REAL_WP, MPI_MAX, lp%cfg%comm, ierr)
   end subroutine calc_max_part_vel
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      integer :: i, j, k, n
      real(WP) :: ypos, ycenter, half_height
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resU(cfg%imino_:cfg%imaxo_, cfg%jmino_:cfg%jmaxo_, cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_, cfg%jmino_:cfg%jmaxo_, cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_, cfg%jmino_:cfg%jmaxo_, cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_, cfg%jmino_:cfg%jmaxo_, cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_, cfg%jmino_:cfg%jmaxo_, cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_, cfg%jmino_:cfg%jmaxo_, cfg%kmino_:cfg%kmaxo_))
         allocate(SR(6, cfg%imino_:cfg%imaxo_, cfg%jmino_:cfg%jmaxo_, cfg%kmino_:cfg%kmaxo_))
         allocate(U  (cfg%imino_:cfg%imaxo_, cfg%jmino_:cfg%jmaxo_, cfg%kmino_:cfg%kmaxo_))
         allocate(V  (cfg%imino_:cfg%imaxo_, cfg%jmino_:cfg%jmaxo_, cfg%kmino_:cfg%kmaxo_))
         allocate(W  (cfg%imino_:cfg%imaxo_, cfg%jmino_:cfg%jmaxo_, cfg%kmino_:cfg%kmaxo_))
         allocate(rho (cfg%imino_:cfg%imaxo_, cfg%jmino_:cfg%jmaxo_, cfg%kmino_:cfg%kmaxo_))
         allocate(visc(cfg%imino_:cfg%imaxo_, cfg%jmino_:cfg%jmaxo_, cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      ! Allocate volume fraction arrays
      allocate_volume_fraction: block
         allocate(particle_VF(cfg%imino_:cfg%imaxo_, cfg%jmino_:cfg%jmaxo_, cfg%kmino_:cfg%kmaxo_))
         allocate(fluid_VF(cfg%imino_:cfg%imaxo_, cfg%jmino_:cfg%jmaxo_, cfg%kmino_:cfg%kmaxo_))
         particle_VF = 0.0_WP
         fluid_VF = 1.0_WP  ! Start with full fluid
      end block allocate_volume_fraction
      
      ! Initialize time tracker
      initialize_timetracker: block
         time = timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size', time%dtmax)
         call param_read('Max cfl number', time%cflmax)
         time%dt = time%dtmax
         time%itmax = 2
      end block initialize_timetracker
      
      ! Create and initialize flow solver with non-periodic BCs
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg, gmres_pfmg
         
         real(WP) :: H
         ! Create flow solver
         fs = incomp(cfg=cfg, name='NS solver')
         
         ! Assign constant viscosity and density
         call param_read('Dynamic viscosity', visc_val)
         call param_read('Density', fs%rho)
         fs%visc = visc_val
         
         ! Initialize entire velocity field to zero first
         fs%U = 0.0_WP
         fs%V = 0.0_WP
         fs%W = 0.0_WP
         
         ! Add boundary conditions
         call fs%add_bcond(name='inflow', type=dirichlet, face='x', dir=-1, canCorrect=.false., locator=left_boundary)
         call fs%add_bcond(name='outflow', type=clipped_neumann, face='x', dir=+1, canCorrect=.true., locator=right_boundary)
         call fs%add_bcond(name='ymin', type=dirichlet, face='y', dir=-1, canCorrect=.false., locator=ymin_boundary)
         call fs%add_bcond(name='ymax', type=slip, face='y', dir=+1, canCorrect=.false., locator=ymax_boundary)
         call fs%add_bcond(name='zmin', type=dirichlet, face='z', dir=-1, canCorrect=.false., locator=zmin_boundary)
         call fs%add_bcond(name='zmax', type=dirichlet, face='z', dir=+1, canCorrect=.false., locator=zmax_boundary)
         
         ! Configure pressure solver
         ps = hypre_str(cfg=cfg, name='Pressure', method=pcg_pfmg, nst=7)
         ps%maxlevel = 10
         call param_read('Pressure iteration', ps%maxit)
         call param_read('Pressure tolerance', ps%rcvg)
         
         ! Configure implicit velocity solver
         vs = hypre_str(cfg=cfg, name='Velocity', method=gmres_pfmg, nst=7)
         call param_read('Implicit iteration', vs%maxit)
         call param_read('Implicit tolerance', vs%rcvg)
         
         ! Setup the solver
         call fs%setup(pressure_solver=ps, implicit_solver=vs)
         
         ! Set inflow profile for OPEN CHANNEL - MODIFIED
         call param_read('Inflow velocity', Uin)  ! Now represents maximum surface velocity
         call fs%get_bcond('inflow', inflow_bc)
         H = cfg%y(cfg%jmax+1) - cfg%y(cfg%jmin)  ! Total depth
         
         ! Initialize the entire domain with the parabolic profile
         do k = cfg%kmino_, cfg%kmaxo_
            do j = cfg%jmino_, cfg%jmaxo_
               do i = cfg%imino_, cfg%imaxo_
                  ypos = cfg%ym(j) - cfg%y(cfg%jmin)    ! Height from bottom
                  fs%U(i,j,k) = Uin * (2.0_WP * ypos/H - (ypos/H)**2)  ! Semi-parabolic profile
               end do
            end do
         end do
         
         ! Apply boundary conditions (will reset boundaries to correct values)
         call fs%apply_bcond(time%t, time%dt)
         
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui, Vi, Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver

      ! Initialize LPT solver
      initialize_lpt: block
         ! Create solver
         lp = lpt(cfg=cfg, name='LPT')
         
         ! Get particle properties
         call param_read('Particle density', lp%rho)
         call param_read('Particle diameter', dp)
         call param_read('Particle viscosity', lp%visc_p, default=visc_val)  ! ADDED
         call param_read('Drag model', lp%drag_model, default='Harper-Moore')
         ! Set filter scale
         lp%filter_width = 3.5_WP * cfg%min_meshsize
         
         ! Setup injection parameters
         setup_injection: block
            call param_read('Mass flow rate', lp%mfr)
            call param_read('Injection position', lp%inj_pos)
            call param_read('Injection velocity', lp%inj_vel)
            call param_read('Injection diameter', lp%inj_d)
            call param_read('Injection mean diameter', lp%inj_dmean, default=dp)
            call param_read('Injection diameter std', lp%inj_dsd, default=0.0_WP)
            call param_read('Injection min diameter', lp%inj_dmin, default=0.5_WP*dp)
            call param_read('Injection max diameter', lp%inj_dmax, default=1.5_WP*dp)
            call param_read('Injection diameter shift', lp%inj_dshift, default=0.0_WP)
         end block setup_injection
         
         ! Get initial particle volume fraction
         call lp%update_VF()
         
         ! Set collision parameters
         call param_read('Collision timescale', lp%tau_col, default=15.0_WP*time%dt)
         call param_read('Restitution coefficient', lp%e_n, default=0.7_WP)
         call param_read('Wall restitution', lp%e_w, default=lp%e_n)
         call param_read('Friction coefficient', lp%mu_f, default=0.0_WP)
         
         ! Set gravity
         call param_read('Gravity', lp%gravity)
      end block initialize_lpt

      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i
         pmesh = partmesh(nvar=1, nvec=2, name='lpt')
         pmesh%varname(1) = 'radius'
         pmesh%vecname(1) = 'velocity'
         pmesh%vecname(2) = 'Fcol'
         call lp%update_partmesh(pmesh)
         do i = 1, lp%np_
            pmesh%var(1,i) = 0.5_WP * lp%p(i)%d
            pmesh%vec(:,1,i) = lp%p(i)%vel
            pmesh%vec(:,2,i) = lp%p(i)%Acol
         end do
      end block create_pmesh   
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output
         ens_out = ensight(cfg=cfg, name='channelMIXXXX')
         
         ! Create event for Ensight output
         ens_evt = event(time=time, name='Ensight output')
         call param_read('Ensight output period', ens_evt%tper)
         
         ! Add variables to output
         call ens_out%add_vector('velocity', Ui, Vi, Wi)
         call ens_out%add_scalar('viscosity', fs%visc)
         call ens_out%add_scalar('particleVF', particle_VF)    ! Added particle volume fraction
         call ens_out%add_scalar('fluidVF', fluid_VF)          ! Added fluid volume fraction
         call ens_out%add_particle('particles', pmesh)
         
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      
      ! Create monitor files
      create_monitor: block
         ! Prepare flow field info
         call fs%get_cfl(time%dt, time%cfl)
         call fs%get_max()
         
         ! Create simulation monitor
         mfile = monitor(fs%cfg%amRoot, 'simulation')
         call mfile%add_column(time%n, 'Timestep number')
         call mfile%add_column(time%t, 'Time')
         call mfile%add_column(time%dt, 'Timestep size')
         call mfile%add_column(time%cfl, 'Maximum CFL')
         call mfile%add_column(fs%Umax, 'Umax')
         call mfile%add_column(fs%Vmax, 'Vmax')
         call mfile%add_column(fs%Wmax, 'Wmax')
         call mfile%add_column(fs%Pmax, 'Pmax')
         call mfile%add_column(fs%divmax, 'Maximum divergence')
         call mfile%add_column(fs%psolv%it, 'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr, 'Pressure error')
         call mfile%add_column(lp%np, 'Particle number')
         call mfile%add_column(mean_y, 'Mean y')
         call mfile%add_column(max_part_vel, 'Max particle vel')
         call mfile%add_column(lp%VFmean, 'Mean particle VF')  ! Added mean particle VF
         call mfile%add_column(lp%VFmax, 'Max particle VF')    ! Added max particle VF
         call mfile%write()
         
         ! Create CFL monitor
         cflfile = monitor(fs%cfg%amRoot, 'cfl')
         call cflfile%add_column(time%n, 'Timestep number')
         call cflfile%add_column(time%t, 'Time')
         call cflfile%add_column(fs%CFLc_x, 'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y, 'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z, 'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x, 'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y, 'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z, 'Viscous zCFL')
         call cflfile%add_column(fs%psolv%rerr, 'Pressure residual')
         call cflfile%write()
      end block create_monitor
      
      ! Initialize particle velocity statistics
      max_part_vel = 0.0_WP
      
   end subroutine simulation_init
   
   !> Time integration routine
   subroutine simulation_run
      implicit none
      integer :: i, j, k, n
      real(WP) :: top_boundary, buffer
      
      ! Perform time integration
      do while (.not.time%done())
         ! Increment time
         call fs%get_cfl(time%dt, time%cfl)
         call time%adjust_dt()
         
         ! Apply particle CFL condition
         if (max_part_vel > 0.0_WP) then
            time%dt = min(time%dt, 0.25_WP * cfg%min_meshsize / max_part_vel)
         end if
         
         call time%increment()
         step_count = step_count + 1

         ! Remember old velocity
         fs%Uold = fs%U
         fs%Vold = fs%V
         fs%Wold = fs%W
         
         ! Perform sub-iterations
         do while (time%it <= time%itmax)
            ! Build mid-time velocity
            fs%U = 0.5_WP*(fs%U + fs%Uold)
            fs%V = 0.5_WP*(fs%V + fs%Vold)
            fs%W = 0.5_WP*(fs%W + fs%Wold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU, resV, resW)
            
            ! Assemble explicit residual
            resU = -2.0_WP*(fs%rho*fs%U - fs%rho*fs%Uold) + time%dt*resU
            resV = -2.0_WP*(fs%rho*fs%V - fs%rho*fs%Vold) + time%dt*resV
            resW = -2.0_WP*(fs%rho*fs%W - fs%rho*fs%Wold) + time%dt*resW
            
            ! Form implicit residuals
            call fs%solve_implicit(time%dt, resU, resV, resW)
            
            ! Apply residuals
            fs%U = 2.0_WP*fs%U - fs%Uold + resU
            fs%V = 2.0_WP*fs%V - fs%Vold + resV
            fs%W = 2.0_WP*fs%W - fs%Wold + resW
            
            ! Apply boundary conditions
            call fs%apply_bcond(time%t, time%dt)
            
            ! Solve Poisson equation
            call fs%correct_mfr()
            call fs%get_div()
            fs%psolv%rhs = -fs%cfg%vol * fs%div * fs%rho / time%dt
            fs%psolv%sol = 0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol, resU, resV, resW)
            fs%P = fs%P + fs%psolv%sol
            fs%U = fs%U - time%dt * resU / fs%rho
            fs%V = fs%V - time%dt * resV / fs%rho
            fs%W = fs%W - time%dt * resW / fs%rho
            
            ! Increment sub-iteration counter
            time%it = time%it + 1
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui, Vi, Wi)
         call fs%get_div()

         ! Transfer fluid data to particle arrays
         do k = cfg%kmino_, cfg%kmaxo_
            do j = cfg%jmino_, cfg%jmaxo_
               do i = cfg%imino_, cfg%imaxo_
                  U(i,j,k)   = Ui(i,j,k)
                  V(i,j,k)   = Vi(i,j,k)
                  W(i,j,k)   = Wi(i,j,k)
                  rho(i,j,k) = fs%rho
                  visc(i,j,k) = fs%visc(i,j,k)
               end do
            end do
         end do

         ! Advance particles using LPT's built-in drag model
         call lp%advance(dt=time%dt, U=U, V=V, W=W, rho=rho, visc=visc)
         
         ! Remove particles that escape domain
         remove_particles: block
            top_boundary = cfg%ym(cfg%jmax+1)
            buffer = 1.0_WP * dp  ! Buffer of one particle diameter
            
            do n = 1, lp%np_
               ! Escape through top (bubble removal)
               if (lp%p(n)%pos(2) >= top_boundary - buffer) lp%p(n)%flag = 1
               ! Leave through outflow
               if (lp%p(n)%pos(1) > cfg%x(cfg%imax+1)) lp%p(n)%flag = 1
            end do
         end block remove_particles
         
         ! Perform expensive operations at reduced frequency
         if (mod(step_count, collision_freq) == 0) call lp%collide(dt=time%dt)
         if (mod(step_count, injection_freq) == 0) call lp%inject(dt=time%dt, avoid_overlap=.true.)
         
         call lp%sync()

         ! Update volume fractions for visualization
         update_volume_fractions: block
            ! Particle volume fraction comes directly from LPT solver
            particle_VF = lp%VF
            
            ! Fluid volume fraction = total cell volume - particle volume fraction
            ! Ensure physical range (0-1)
            do k = cfg%kmino_, cfg%kmaxo_
               do j = cfg%jmino_, cfg%jmaxo_
                  do i = cfg%imino_, cfg%imaxo_
                     fluid_VF(i,j,k) = max(0.0_WP, cfg%VF(i,j,k) - particle_VF(i,j,k))
                  end do
               end do
            end do
            
            ! Synchronize across processors
            call cfg%sync(particle_VF)
            call cfg%sync(fluid_VF)
         end block update_volume_fractions

         ! Output to ensight
         if (ens_evt%occurs()) then
            update_pmesh: block
               integer :: i
               call lp%update_partmesh(pmesh)
               do i = 1, lp%np_
                  pmesh%var(1,i)   = 0.5_WP * lp%p(i)%d
                  pmesh%vec(:,1,i) = lp%p(i)%vel
                  pmesh%vec(:,2,i) = lp%p(i)%Acol
               end do
            end block update_pmesh
            call ens_out%write_data(time%t)
         end if

         ! Monitoring
         call fs%get_max()
         call lp%get_max()
         call calc_barycenter()
         call calc_max_part_vel()
         call mfile%write()
         call cflfile%write()
      end do
      
      ! Output final particle state
      call lp%write(filename='part.file')
   end subroutine simulation_run
   
   !> Finalize the simulation
   subroutine simulation_final
      implicit none
      deallocate(resU, resV, resW, Ui, Vi, Wi, SR, U, V, W, rho, visc)
      deallocate(particle_VF, fluid_VF)
   end subroutine simulation_final
   
end module simulation