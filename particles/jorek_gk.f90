module domain
  real*8 :: psi_start, psi_end, psi_n_start, psi_n_end
end module

module mod_normalisations
real*8 :: zn_norm, rho_norm, t_norm, Tev_norm, v_norm, m_norm, E_norm
end module

program jorek_gk
use mod_particle_types
use mod_particle_sim
use mod_particle_allocation, only: allocate_particles_for_sim
use mod_initialise_particles
use particle_tracer, only: sim, events
use mod_particle_io
use mod_pusher_tools, only : particle_position_to_gc
use mod_particle_diagnostics
use mod_projection_functions
use mod_rhs_projections
use mpi
use mod_import_restart
use mod_atomic_elements
use mod_event
use mod_io_actions
use mod_fields_linear
use mod_project_particles
use mod_gc_variational
use nodes_elements
use mod_poisson_solver, only: poisson_solve_action
use mod_poisson_rhs, only: assemble_direct_poisson_rhs
use mod_uncoupled_projection, only: assemble_projection_rhs
use mod_simulation_data, only: type_MHD_SIM
use mod_sobseq_rng
use mod_pcg32_rng
use mod_random_seed
use mod_interp, only: mode_moivre, interp_RZ, interp_0, interp
use mod_basisfunctions
use basis_at_gaussian
use mod_normalisations
use phys_module, only: F0, tstep, nstep, nout, restart, rho_0, rho_1, rho_coef
use phys_module, only: CENTRAL_MASS, CENTRAL_DENSITY, xcase, xpoint, index_now, index_start
use phys_module, only: nstep_particles, nsubstep_particles, tstep_particles, nsubstep_electrons
use phys_module, only: part_group_configs, n_part_groups
use phys_module, only: filter_perp, filter_hyper, filter_par, filter_perp_n0, filter_hyper_n0, filter_par_n0
use phys_module, only: xtime, energies, mode, restart_particles, n_tht, n_leg
use phys_module, only: T_scale_factor, B_scale_factor, t_now
use phys_module, only: rst_format, ei_small_angle_scattering_Orb5, ei_small_angle_scattering_Lu
use constants,   only: MU_ZERO, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG
use mod_export_restart
use live_data

use mod_edge_domain
use mod_edge_elements
use data_structure
use equil_info
use mod_boundary, only: boundary_from_grid
use domain
use mod_fields

!$ use omp_lib

implicit none

!type(particle_sim)                                :: sim
!type(event), dimension(:), allocatable, target    :: events
type(event)                                       :: fieldreader, density_reader, partreader, partwriter
type(count_action)                                :: counter
type(projection), target                          :: jorek_feedback, project_profiles, project_density
type(poisson_solve_action)                        :: poisson
type(type_MHD_SIM), target                        :: poisson_mhd_sim
type(type_RHS)                                    :: deposition_rhs, poisson_rhs
type(type_edge_domain), allocatable, dimension(:) :: edge_domains
type(edge_elements)                               :: D_edge
type(write_particle_diagnostics)                  :: diag
type(type_node_list), target                      :: eq_node_list
type(type_element_list), target                   :: eq_element_list
class(fields_base), pointer                       :: eq_fields
!type(type_bnd_element_list) :: bnd_elm_list !< List of boundary elements
!type(type_bnd_node_list)    :: bnd_node_list !< List of boundary nodes.

logical :: update_electric_potential ! move to input?
integer :: ion_group, electron_group
integer :: n_particles_out
integer :: n_update_density
integer :: n_update_temperature
integer :: n_update_poisson
integer :: n_update_profiles
real*8  :: n_ions, n_electrons
integer :: n_particles_local, n_ions_local, n_electrons_local  
real*8  :: tstep_keep, particle_start_time
real*8  :: psi_axis, psi_bnd
real*8  :: total_particles, total_volume
real*8  :: W_mag(n_tor), W_kin(n_tor), W_tot, t0, t1, t2
real*4, external :: f_ions, f_electrons, f_density
real*8, external :: T_ions, T_electrons
integer   :: ifail, node_start, node_end
integer   :: i, j, k, l, m, in, jn, inode, i_elm, index_rhs, i_diagno(6)
integer   :: seed, i_rng, n_stream, ierr, n_particle_out, i_tor, i_tor_start
integer   :: n_ions_lost, n_electrons_lost, n_ions_lost_max, n_electrons_lost_max
integer, allocatable :: index_lost_ions(:), index_lost_electrons(:)
real*8    :: energy_local_lost_ions, energy_total_lost_ions, energy_local_lost_electrons, energy_total_lost_electrons, sum_energy
real*8    :: energy_local_ions, energy_local_electrons, energy_total_ions, energy_total_electrons
real*8    :: momentum_local_lost_ions, momentum_total_lost_ions, momentum_local_lost_electrons, momentum_total_lost_electrons, sum_momentum
real*8    :: momentum_local_ions, momentum_local_electrons, momentum_total_ions, momentum_total_electrons
real*8    :: potential_local_lost_ions, potential_total_lost_ions, potential_local_lost_electrons, potential_total_lost_electrons
real*8    :: potential_local_ions,      potential_total_ions,      potential_local_electrons,      potential_total_electrons
integer   :: total_ions_lost, total_electrons_lost, n_zeros
real*8    :: density_min, density_s, density_t, density_st, mass_ratio, psi_n, tstep_electrons, zero_fraction

character*18 :: fileout, filepart
real*8, allocatable  :: rhs_nodes(:,:), rhs_nodes_local(:,:)
real*8, allocatable  :: reservoir_E_ions(:), reservoir_E_electrons(:), reservoir_P_ions(:), reservoir_P_electrons(:), reservoir_delta_E(:), reservoir_delta_P(:)
integer              :: n_patches
integer, allocatable :: patch_list(:)

real*8, allocatable :: feedback_rhs_ions(:,:,:,:,:)


psi_n_start = 0.d0          ! limit domain to psi_start:psi_end (in normalised psi)
psi_n_end   = 1.5d0 

update_electric_potential = .true.

zero_fraction        = 0.1d0 ! create some space for new heating and fuelling particles

n_particles_out      = 500   ! in steps
n_update_density     = 999999
n_update_profiles    = 999999 !20
n_update_temperature = 999999 !20  ! must be a multiple of n_update_profiles
n_update_poisson     = 999999 !500

ion_group      = 1 ! to be derived from namelist input later
electron_group = 2

!call sim%initialize(skip_group_config =.true.) 
call sim%initialize()  

if (nsubstep_electrons < 1) then
  if (sim%my_id == 0) write(*,*) 'ERROR: nsubstep_electrons must be at least one'
  call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
endif

zn_norm  = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm = CENTRAL_MASS * ATOMIC_MASS_UNIT * zn_norm            ! rho_SI = rho_norm * rho
t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
E_norm   = 1.5d0 / MU_ZERO                                      ! E_SI   = E_norm * E_jorek
M_norm   = rho_norm * v_norm                                    ! momentum normalisation
Tev_norm = 1.d0 / (EL_CHG * MU_ZERO * zn_norm)                  ! T_ev [eV] = Tev_norm * T_jorek / 2.d0

write(*,*) ' part_num_groups : ',n_part_groups

if (restart_particles) restart = .true.
if (sim%my_id .eq. 0) write(*,*) 'RESTART = ',restart

n_ions            = part_group_configs(ion_group)%n_particles
n_electrons       = part_group_configs(electron_group)%n_particles
n_ions_local      = int(n_ions/sim%n_mpi) 
n_electrons_local = int(n_electrons/sim%n_mpi )

write(*,'(i3,A,2i9.2e12.4)') sim%my_id,' number of particles : ',n_ions_local, n_electrons_local, n_ions, n_electrons

open(113,file='energies.txt')

if (sim%my_id .eq. 0) call init_live_data()

!call import_restart(eq_node_list, eq_element_list, 'equil_base', rst_format, ierr, .true.)
!eq_fields%node_list    => eq_node_list
!eq_fields%element_list => eq_element_list

if (restart_particles) then
  deallocate(sim%groups)
  partreader = event(read_action(filename='restart_particles.h5'))
  call with(sim, partreader)
endif

fieldreader = event(read_jorek_fields_interp_linear(basename='equil_base', i=-1))
call with(sim, fieldreader)
if (sim%my_id .eq. 0) write(*,*) 'done reading equilibrium'

tstep     = tstep_keep
index_now = index_start

call det_modes()

write(*,*) 'SCALE_FACTORS : ',sim%my_id,B_scale_factor,T_scale_factor

if (.not. restart_particles) then
  do j=1, sim%fields%node_list%n_nodes
    sim%fields%node_list%node(j)%values(:,:,2) = 0.d0
    sim%fields%node_list%node(j)%values(:,:,1) = B_scale_factor * sim%fields%node_list%node(j)%values(:,:,1)
    sim%fields%node_list%node(j)%values(:,:,6) = T_scale_factor * sim%fields%node_list%node(j)%values(:,:,6)
  enddo
  F0 = B_scale_factor * F0 
endif

if (sim%my_id .eq. 0) then
  do i=1, index_start
    call write_live_data(i)
  enddo 
endif

zn_norm   = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm  = CENTRAL_MASS * ATOMIC_MASS_UNIT * zn_norm            ! rho_SI = rho_norm * rho
t_norm    = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm   * t_jorek
Tev_norm  = 1.d0 / (2.d0 * EL_CHG * MU_ZERO * zn_norm)           ! T_ev   = Tev_norm * T (factor 2 for electron temperature)
v_norm    = 1.d0 / t_norm                                        ! V_SI   = v_norm   * vpar_jorek * B
if (sim%my_id .eq. 0) write(*,*) ' zn_norm : ',zn_norm

if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)
call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)
call update_equil_state(sim%my_id, sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase)

psi_axis  = ES%psi_axis
psi_bnd   = ES%psi_bnd
psi_start = psi_axis + (psi_bnd - psi_axis) * psi_n_start
psi_end   = psi_axis + (psi_bnd - psi_axis) * psi_n_end

if (sim%my_id .eq. 0) write(*,'(A,6e14.6)') 'PSI_AXIS, PSI_BND : ',psi_axis, psi_bnd, psi_start, psi_end, psi_n_start, psi_n_end

!call grid_reduction(sim%fields%node_list,sim%fields%element_list,psi_n_start,psi_n_end,n_tht)

!allocate(patch_list(sim%fields%element_list%n_elements))
!if (xpoint) then
!  call define_patches_xpoint(sim%fields%node_list, sim%fields%element_list, patch_list, n_patches)
!else
!  call define_patches(sim%fields%node_list, sim%fields%element_list, patch_list, n_patches)
!endif

if (.not. restart_particles) then
  sim%groups(ion_group)%Z         = part_group_configs(ion_group)%Z
  sim%groups(ion_group)%mass      = part_group_configs(ion_group)%mass            ! atomic_weights(-2)  !< atomic mass units, -2 for deuterium
  sim%groups(ion_group)%n_particles = part_group_configs(ion_group)%n_particles  

  sim%groups(electron_group)%Z    = part_group_configs(electron_group)%Z          
  sim%groups(electron_group)%mass = part_group_configs(electron_group)%mass       ! atomic_weights(-2)/mass_ratio !< atomic mass units, -1 for electrons (here heavy electrons)
  sim%groups(electron_group)%n_particles = part_group_configs(electron_group)%n_particles  
  
 if (sim%my_id .eq. 0) write(*,'(A,2e12.4)') ' ion/electron mass : ',sim%groups(ion_group)%mass, sim%groups(electron_group)%mass

!  allocate(particle_gc_vpar::sim%groups(ion_group)%particles(n_ions_local))
!  allocate(particle_gc_vpar::sim%groups(electron_group)%particles(n_electrons_local))
  call allocate_particles_for_sim(sim) ! populate the particle arrays in the particle groups

  call cpu_time(t0)
  call initialise_particles_H_mu_psi(sim%groups(ion_group)%particles, sim%fields, sobseq_rng(), sim%groups(ion_group)%mass, &
                                     uniform_space=.true., uniform_space_rej_f=f_ions, &
!                                     uniform_space=.true., uniform_space_rej_f=f_density, &
!                                     uniform_space_rej_vars=[-2,1], charge = +1)
                                     uniform_space_rej_vars=[-2,1], charge = +1)
  call cpu_time(t1)
  write(*,*) ' cpu time initisalise ions : ', t1-t0

  call initialise_particles_H_mu_psi(sim%groups(electron_group)%particles, sim%fields, sobseq_rng(), sim%groups(electron_group)%mass, &
                                     uniform_space=.true., uniform_space_rej_f=f_electrons, &
!                                     uniform_space=.true., uniform_space_rej_f=f_density, &
!                                     uniform_space_rej_vars=[-2,1], charge = -1)
                                     uniform_space_rej_vars=[-2,1], charge = -1)
  call cpu_time(t2)
  write(*,*) ' cpu time initisalise electrons : ', t2-t1

  do i=1, size(sim%groups(ion_group)%particles)
    if (sim%groups(ion_group)%particles(i)%x(2) .lt. ES%z_xpoint(1) ) sim%groups(ion_group)%particles(i)%i_elm = 0
  enddo
  do i=1, size(sim%groups(electron_group)%particles)
    if (sim%groups(electron_group)%particles(i)%x(2) .lt. ES%z_xpoint(1) ) sim%groups(electron_group)%particles(i)%i_elm = 0
  enddo

  n_zeros = int(zero_fraction * size(sim%groups(ion_group)%particles)) 
  do i=size(sim%groups(ion_group)%particles) - n_zeros + 1, size(sim%groups(ion_group)%particles) 
    sim%groups(ion_group)%particles(i)%i_elm = 0
  enddo
  n_zeros = int(zero_fraction * size(sim%groups(electron_group)%particles)) 
  do i=size(sim%groups(electron_group)%particles) - n_zeros + 1, size(sim%groups(electron_group)%particles) 
    sim%groups(electron_group)%particles(i)%i_elm = 0
  enddo
  
  call density_integral(sim%fields%node_list,sim%fields%element_list,f_density,min(psi_axis,psi_bnd),max(psi_axis,psi_bnd),total_particles,total_volume) 

  total_particles = total_particles    * zn_norm
  if (sim%my_id .eq. 0) write(*,'(A,8e14.6)') ' total particles, volume : ',total_particles, total_volume

  call adjust_particle_weights(sim%groups(ion_group)%particles, total_particles)
  call adjust_particle_weights(sim%groups(electron_group)%particles, total_particles)

  if (sim%my_id .eq. 0) write(*,'(A,3e14.6)') ' Ion particle density was adjusted to      : ', total_particles, sim%groups(ion_group)%particles(1)%weight
  if (sim%my_id .eq. 0) write(*,'(A,3e14.6)') ' Electron particle density was adjusted to : ', total_particles, sim%groups(electron_group)%particles(1)%weight

  call with(sim, counter)        

  n_ions_local      = size(sim%groups(ion_group)%particles,1)
  n_electrons_local = size(sim%groups(electron_group)%particles,1)
  allocate(index_lost_ions(int(zero_fraction*n_ions_local)),     index_lost_electrons(int(zero_fraction*n_electrons_local)))

  n_ions_lost      = 0
  n_electrons_lost = 0
  n_ions_lost_max      = size(index_lost_ions)
  n_electrons_lost_max = size(index_lost_electrons)
endif  ! restart_particles

project_profiles = new_projection(sim%fields%node_list, sim%fields%element_list, &
                      filter    = filter_perp,    filter_hyper    = filter_hyper,    filter_parallel    = filter_par,    &
                      filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, filter_parallel_n0 = filter_par_n0, &
                      f=[proj_f(proj_one,      group = ion_group), proj_f(proj_one,      group = electron_group),        &
                         proj_f(proj_Pressure, group = ion_group), proj_f(proj_Pressure, group = electron_group),        &
                         proj_f(proj_vpar,     group = ion_group), proj_f(proj_vpar,     group = electron_group)],       &
!                      do_dirichlet_open_n0 = .false.,   do_dirichlet_corners_n0 = .false., do_neumann_n0 = .true.,      &
                      do_dirichlet =.false.,                                                                             &
                      fractional_digits = 9, calc_integrals=.false., to_vtk=.true., to_h5=.true., basename='profiles', nsub=5)

call with(sim, project_profiles)
call project_profiles%set_vtk_active(.false.)

project_density = new_projection(sim%fields%node_list, sim%fields%element_list, &
                      filter    = filter_perp,    filter_hyper    = filter_hyper,    filter_parallel    = filter_par,    &
                      filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, filter_parallel_n0 = filter_par_n0, &
                      f=[proj_f(proj_one, group = ion_group)],                                                           & 
!                      do_dirichlet_open_n0 = .false.,   do_dirichlet_corners_n0 = .false., do_neumann_n0 = .true. ,     &
                      do_dirichlet =.false.,                                                                             &
                      fractional_digits = 9, calc_integrals=.true., to_vtk=.true., to_h5=.true., basename='density', nsub=5)

!allocate(project_density%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 4))  !for ion and electron energy and momentum
!project_density%rhs = 0.d0

call with(sim, project_density)

call MPI_BARRIER(MPI_COMM_WORLD, ierr)

if (nstep .gt. 0 .and. update_electric_potential) then
  jorek_feedback = new_projection(sim%fields%node_list, sim%fields%element_list, &
                     filter    = filter_perp,    filter_hyper    = filter_hyper,    filter_parallel    = filter_par,    &
                     filter_n0 = filter_perp_n0, filter_hyper_n0 = filter_hyper_n0, filter_parallel_n0 = filter_par_n0, &
                     fractional_digits = 9,         &
!                     do_neumann_n0=.false., do_dirichlet_open_n0 = .true., do_dirichlet_corners_n0 = .true., &
                     do_dirichlet = .true.,                                                                                &
                     calc_integrals=.false., to_vtk=.true., to_h5 = .false., basename='projections')

  allocate(jorek_feedback%rhs(n_order+1, n_vertex_max, sim%fields%element_list%n_elements, n_tor, 1))

  jorek_feedback%rhs = 0.d0
  jorek_feedback%scaling_integral_weights = 0.00d0 ! subtract 1.d0 from the rhs of n=0 for adiabatic electrons

  allocate(feedback_rhs_ions,source=jorek_feedback%rhs)
  
  aux_node_list => jorek_feedback%node_list

  if (restart) then  
    fieldreader = event(read_jorek_fields_interp_linear(basename='jorek', i=-1))
    call with(sim, fieldreader)

    if (sim%my_id .eq. 0) call boundary_from_grid(sim%fields%node_list, sim%fields%element_list, bnd_node_list, bnd_elm_list, .false.)
    call broadcast_boundary(sim%my_id, bnd_elm_list, bnd_node_list)
    call update_equil_state(sim%my_id, sim%fields%node_list, sim%fields%element_list, bnd_elm_list, xpoint, xcase)

    tstep     = tstep_keep
    index_now = index_start

    do j=1, sim%fields%node_list%n_nodes
      sim%fields%node_list%node(j)%values(:,:,2) = 0.d0
      sim%fields%node_list%node(j)%values(:,:,1) = B_scale_factor * sim%fields%node_list%node(j)%values(:,:,1)
      sim%fields%node_list%node(j)%values(:,:,6) = T_scale_factor * sim%fields%node_list%node(j)%values(:,:,6)
    enddo
    F0 = B_scale_factor * F0 
  endif
  
  events = [new_event_ptr(jorek_feedback,  start = sim%time), event(stop_action(), start=1d12)  ]

! Call events at sim%time once to help event scheduler, before entering particle loop
! call with(sim, events, at=sim%time)

endif

allocate(rhs_nodes(4,sim%fields%node_list%n_nodes),rhs_nodes_local(4,sim%fields%node_list%n_nodes))

node_start = 1
node_end   = sim%fields%node_list%n_nodes

if (nstep .gt. 0 .and. update_electric_potential) then
  poisson_mhd_sim%my_id = sim%my_id
  poisson_mhd_sim%n_mpi = sim%n_mpi
  poisson_mhd_sim%n_tor = n_tor
  poisson_mhd_sim%freeboundary = .false.
  poisson_mhd_sim%restart = restart
  poisson_mhd_sim%sr_n_tor = 0
  poisson_mhd_sim%node_list => sim%fields%node_list
  poisson_mhd_sim%element_list => sim%fields%element_list
  poisson_mhd_sim%bnd_node_list => bnd_node_list
  poisson_mhd_sim%bnd_elm_list => bnd_elm_list
  poisson_mhd_sim%es => ES
  call poisson%setup(poisson_mhd_sim,MPI_COMM_WORLD)
  call poisson%construct_matrix()
endif

do i=1, nstep_particles

  particle_start_time = sim%time

  index_now = index_now + 1

  jorek_feedback%rhs = 0.d0

  jorek_feedback%rhs_gather_time = nsubstep_particles * tstep_particles
  
  reservoir_E_ions      = 0.d0
  reservoir_P_ions      = 0.d0
  reservoir_E_electrons = 0.d0
  reservoir_P_electrons = 0.d0
  reservoir_delta_E     = 0.d0
  reservoir_delta_P     = 0.d0
  project_density%rhs(:,:,:,:,1:4) = 0.d0

  call loop_particle_gc_local(sim, ion_group, 4, ion_group, electron_group,              &
                              jorek_feedback, project_profiles, project_density,         &
                              tstep_particles, nsubstep_particles, particle_start_time, &
                              index_lost_ions, n_ions_lost, n_ions_lost_max,             &
                              energy_local_ions,      energy_local_lost_ions,            &
                              momentum_local_ions,    momentum_local_lost_ions,          &
                              potential_local_ions,   potential_local_lost_ions,         &
                              .true.)   ! ions                           
                              
  feedback_rhs_ions = jorek_feedback%rhs               ! the ion contribution to the right hand side

!--- reservoir_E_ions      contains energy to be transferred to ions (from collisions with background ions)
!--- reservoir_E_electrons contains energy to be transferred to ions (from collisions with background electrons)
!--- i.e. the labels ion/electron refer to the background species
!  call MPI_ALLREDUCE(MPI_IN_PLACE, reservoir_E_ions,      size(reservoir_E_ions),      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!  call MPI_ALLREDUCE(MPI_IN_PLACE, reservoir_P_ions,      size(reservoir_P_ions),      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!  call MPI_ALLREDUCE(MPI_IN_PLACE, reservoir_E_electrons, size(reservoir_E_electrons), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!  call MPI_ALLREDUCE(MPI_IN_PLACE, reservoir_P_electrons, size(reservoir_P_electrons), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
 
!  if (sim%my_id .eq.0) write(*,'(A,4e18.10)') 'ions : reservoir_ions (E/P)      : ',sum(reservoir_E_ions),     sum(reservoir_P_ions)
!  if (sim%my_id .eq.0) write(*,'(A,4e18.10)') 'ions : reservoir_electrons (E/P) : ',sum(reservoir_E_electrons),sum(reservoir_P_electrons)
!  if (sim%my_id .eq.0) write(*,'(A,4e18.10)') 'ions : reservoir_delta (E/P)     : ',sum(reservoir_delta_E),    sum(reservoir_delta_P)
                            
  call with(sim, project_density)

 ! call empty_reservoir(sim, project_density, project_profiles, 1, 1, reservoir_E_ions,      reservoir_P_ions,      patch_list, n_patches)
 ! call empty_reservoir(sim, project_density, project_profiles, 1, 2, reservoir_E_electrons, reservoir_P_electrons, patch_list, n_patches)

  write(*,*) sim%my_id, ' nsubstep_electrons, t_step_electrons : ',nsubstep_electrons

  tstep_electrons = tstep_particles / nsubstep_electrons

  do k=1, nsubstep_electrons
                            
    reservoir_E_ions      = + reservoir_delta_E      ! use delta_E of ie collisions also for the ei collisions (to ensure conservation)
    reservoir_P_ions      = + reservoir_delta_P
    reservoir_E_electrons = 0.d0
    reservoir_P_electrons = 0.d0
    reservoir_delta_E     = 0.d0
    reservoir_delta_P     = 0.d0
    project_density%rhs(:,:,:,:,1:4) = 0.d0

    jorek_feedback%rhs = feedback_rhs_ions !???????????????? WHY THIS

    call loop_particle_gc_local(sim, electron_group, 0,  ion_group, electron_group,           &
                                jorek_feedback, project_profiles, project_density,            &
                                tstep_electrons, nsubstep_particles, particle_start_time,     &
                                index_lost_electrons, n_electrons_lost, n_electrons_lost_max, &
                                energy_local_electrons,    energy_local_lost_electrons,       &
                                momentum_local_electrons,  momentum_local_lost_electrons,     &
                                potential_local_electrons, potential_local_lost_electrons,    &
                                .true.)   ! electrons

    reservoir_E_ions = reservoir_E_ions + reservoir_delta_E      ! use delta_E of ie collisions also for the ei collisions (to ensure conservation)
    reservoir_P_ions = reservoir_P_ions + reservoir_delta_P

    do i_elm=1,sim%fields%element_list%n_elements    
      do in=1,n_vertex_max      
         inode = sim%fields%element_list%element(i_elm)%vertex(in)
         psi_n = (sim%fields%node_list%node(inode)%values(1,1,1) - ES%psi_axis) / (ES%psi_bnd - ES%psi_axis) 
         if ( (inode .lt. node_start) .or. (inode .gt. node_end) .or. (psi_n .gt. 1.1d0)) then
           jorek_feedback%rhs(:, in, i_elm, 1, 1)    = 0.d0
         endif
      enddo
    enddo

    if (update_electric_potential) then

      write(*,*) 'CHECK feedback : ',maxval(jorek_feedback%rhs),minval(jorek_feedback%rhs)

      call assemble_projection_rhs(sim%fields%node_list,sim%fields%element_list, &
           jorek_feedback%rhs(:,:,:,:,1),deposition_rhs,MPI_COMM_WORLD)
      call assemble_direct_poisson_rhs(deposition_rhs,poisson_rhs)
      call poisson%set_rhs(poisson_rhs)
      call poisson%solve()
      call poisson%gather()
      call poisson%store_phi()

      node_start = 1 
      node_end   = sim%fields%node_list%n_nodes
    else
      do j=1, node_end
        sim%fields%node_list%node(j)%values(:,1:4,var_u) = 0.d0
      enddo
    endif  

  enddo ! substep electrons

!--- reservoir_E_ions      contains energy to be transferred to electrons (from collisions with background ions)
!--- reservoir_E_electrons contains energy to be transferred to electrons (from collisions with background electrons)
!--- i.e. the labels ion/electron refer to the background species
!  call MPI_ALLREDUCE(MPI_IN_PLACE, reservoir_E_ions,      size(reservoir_E_ions),      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!  call MPI_ALLREDUCE(MPI_IN_PLACE, reservoir_P_ions,      size(reservoir_P_ions),      MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!  call MPI_ALLREDUCE(MPI_IN_PLACE, reservoir_E_electrons, size(reservoir_E_electrons), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!  call MPI_ALLREDUCE(MPI_IN_PLACE, reservoir_P_electrons, size(reservoir_P_electrons), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

!  if (sim%my_id .eq.0) write(*,'(A,4e18.10)') 'e    : reservoir_ions (E/P)      : ',sum(reservoir_E_ions),     sum(reservoir_P_ions)
!  if (sim%my_id .eq.0) write(*,'(A,4e18.10)') 'e    : reservoir_electrons (E/P) : ',sum(reservoir_E_electrons),sum(reservoir_P_electrons)
!  if (sim%my_id .eq.0) write(*,'(A,4e18.10)') 'e    : reservoir_deltas (E/P)    : ',sum(reservoir_delta_E),    sum(reservoir_delta_P)

  call with(sim, project_density)
  
  project_density%rhs = 0.d0

  sim%time = particle_start_time + tstep_particles * nsubstep_particles
  t_now    = sim%time

  if (sim%my_id == 0) xtime(index_now) = sim%time

  if (sim%my_id == 0) write(*,*) sim%my_id,'TIMESTEPS : ',index_now,xtime(index_now),sim%time,t_now

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(n_ions_lost,         total_ions_lost,       1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(n_electrons_lost,    total_electrons_lost,  1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(energy_local_ions,             energy_total_ions,             1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(energy_local_lost_ions,        energy_total_lost_ions,        1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(energy_local_electrons,        energy_total_electrons,        1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(energy_local_lost_electrons,   energy_total_lost_electrons,   1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(momentum_local_ions,           momentum_total_ions,           1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(momentum_local_lost_ions,      momentum_total_lost_ions,      1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(momentum_local_electrons,      momentum_total_electrons,      1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(momentum_local_lost_electrons, momentum_total_lost_electrons, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(potential_local_ions,           potential_total_ions,           1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(potential_local_lost_ions,      potential_total_lost_ions,      1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(potential_local_electrons,      potential_total_electrons,      1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_REDUCE(potential_local_lost_electrons, potential_total_lost_electrons, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
 
  if (sim%my_id .eq. 0) write(*,'(A,E14.6,2i10,20e18.10)') 'lost ions/electrons, energy : ',sim%time, total_ions_lost, total_electrons_lost, &
                              energy_total_lost_ions,   energy_total_lost_electrons,   energy_total_ions,      energy_total_electrons,       &
                              momentum_total_lost_ions, momentum_total_lost_electrons, momentum_total_ions,     momentum_total_electrons,    &
                              potential_total_lost_ions, potential_total_lost_electrons, potential_total_ions, potential_total_electrons
  if (sim%my_id .eq. 0) write(113,'(E14.6,2i10,20e18.10)') sim%time, total_ions_lost, total_electrons_lost,                                  &
                              energy_total_lost_ions,   energy_total_lost_electrons,   energy_total_ions,      energy_total_electrons,       &
                              momentum_total_lost_ions, momentum_total_lost_electrons, momentum_total_ions,    momentum_total_electrons,     &
                              potential_total_lost_ions, potential_total_lost_electrons, potential_total_ions, potential_total_electrons


  if (sim%my_id .eq. 0) then
    call itg_energy(node_list,element_list,min(psi_start,psi_end),max(psi_start,psi_end),W_kin,W_tot)
    energies(:,1,index_now) = W_kin(:)
    call write_live_data(index_now)
    write(*,'(A,32e14.6)') 'energies     : ',sim%time, W_kin
  endif

   
  if (mod(index_now,n_update_density).eq. 0) then                    ! update ion density in jorek background (i.e. sim%fields)
    if (sim%my_id .eq. 0) write(*,'(A,2e16.8)') 'updating density profile in sim%fields ',maxval(abs(project_density%rhs(1,:,:,1,1))),maxval(abs(project_density%rhs(1,:,:,1,2)))
    if (mod(index_now,nout).eq. 0) call project_density%set_vtk_active(.true.)
    call with(sim, project_density)
    if (sim%my_id .eq. 0) write(*,'(A,2e16.8)') 'after updating density profile ',maxval(abs(project_density%rhs(1,:,:,1,1))),maxval(abs(project_density%rhs(1,:,:,1,2)))
    call project_density%set_vtk_active(.false.)
    do j=1, node_end
      sim%fields%node_list%node(j)%values(1:n_tor,1:4,5) = project_density%node_list%node(j)%values(1:n_tor,1:4,5) / (CENTRAL_DENSITY * 1.d20)
    enddo
  endif

  if (mod(index_now,n_update_profiles).eq. 0) then
    if (sim%my_id .eq. 0) write(*,'(A)') 'updating profiles for background Maxwellians'

    call with(sim, project_profiles)

    if ((mod(index_now,n_update_temperature).eq. 0)) then             ! update temperature in jorek background (i.e. sim%fields)
      
      if (sim%my_id .eq. 0) write(*,'(A)') 'updating temperature profile in sim%fields'

      do j=1, sim%fields%node_list%n_nodes
  
        density_min = max(project_profiles%node_list%node(j)%values(1,1,2), 0.05*central_density * 1d20)  ! n=0 only for the density
        density_s   =     project_profiles%node_list%node(j)%values(1,2,2)
        density_t   =     project_profiles%node_list%node(j)%values(1,3,2)
        density_st  =     project_profiles%node_list%node(j)%values(1,4,2)
                    
        sim%fields%node_list%node(j)%values(:,1,6) =  project_profiles%node_list%node(j)%values(:,1,4) / density_min 
        sim%fields%node_list%node(j)%values(:,2,6) =  project_profiles%node_list%node(j)%values(:,2,4) / density_min    &
                                                   -  project_profiles%node_list%node(j)%values(:,1,4) * density_s / density_min**2
        sim%fields%node_list%node(j)%values(:,3,6) =  project_profiles%node_list%node(j)%values(:,3,4) / density_min    &
                                                   -  project_profiles%node_list%node(j)%values(:,1,4) * density_t / density_min**2
        sim%fields%node_list%node(j)%values(:,4,6) =  project_profiles%node_list%node(j)%values(:,4,4) / density_min    &
                                                   -  project_profiles%node_list%node(j)%values(:,2,4) * density_t  / density_min**2 &
                                                   -  project_profiles%node_list%node(j)%values(:,3,4) * density_s  / density_min**2 &
                                                   -  project_profiles%node_list%node(j)%values(:,1,4) * density_st / density_min**2 &
                                             + 2.d0 * project_profiles%node_list%node(j)%values(:,1,4) * density_s * density_t / density_min**3
                    
        sim%fields%node_list%node(j)%values(:,:,6) = sim%fields%node_list%node(j)%values(:,:,6) * zn_norm * MU_ZERO * 2.d0 ! factor 2, i.e. storing 2*Te (as in JOREK) but should change to Te,Ti                                             
              
        sim%fields%node_list%node(j)%values(1,1,6) = max(sim%fields%node_list%node(j)%values(1,1,6), 2.d0 * 10.d0 / Tev_norm) ! storing 2*Te in sim%fields    
      enddo
    endif
  endif

  ! Rebuild only at the explicit outer-step cadence, after this step's density
  ! and profile refreshes. The next physical solve factorizes the new matrix.
  if (update_electric_potential .and. mod(index_now,n_update_poisson).eq.0) &
    call poisson%construct_matrix()

  if (mod(index_now,nout).eq. 0) then
    call project_profiles%set_vtk_active(.true.)
    call with(sim, project_profiles)
    call project_profiles%set_vtk_active(.false.)
    if (sim%my_id .eq. 0) then
      write(fileout,'(A8,i5.5)') 'profiles',index_now  
      call export_restart(project_profiles%node_list, project_profiles%element_list, fileout)
    endif      
    if (sim%my_id .eq. 0) then
      write(fileout,'(A7,i5.5)') 'density',index_now  
      call export_restart(project_density%node_list, project_density%element_list, fileout)
    endif      
    if (sim%my_id .eq.0) then
      write(fileout,'(A5,i5.5)') 'jorek',index_now
      call export_restart(sim%fields%node_list, sim%fields%element_list, fileout)
    endif
  endif

  if (n_particle_out .gt. 0) then
    if (mod(index_now,n_particle_out) == 0) then
      write(filepart,'(A4,i5.5,A3)') 'part',index_now,'.h5'
      partwriter = event(write_action(filename=filepart))
      call with(sim, partwriter)
    endif
  endif

end do

if (nstep_particles .gt. 0) then
  
  close(111); close(112)

  if (sim%my_id .eq. 0) call export_restart(sim%fields%node_list, sim%fields%element_list, 'restart_jorek')

  partwriter = event(write_action(filename='restart_particles.h5'))
  call with(sim, partwriter)

  if ( sim%my_id .eq. 0 ) call finalize_live_data()

  if (update_electric_potential) then
    call poisson%finalize()
    if (associated(deposition_rhs%val)) deallocate(deposition_rhs%val)
    deposition_rhs%val => null()
    deposition_rhs%n = 0
    if (associated(poisson_rhs%val)) deallocate(poisson_rhs%val)
    poisson_rhs%val => null()
    poisson_rhs%n = 0
  endif
  
  call sim%finalize

endif

contains

subroutine electron_ion_collision_frequency(zni, Te, electron_mass, collision_frequency)
use constants
real*8                  :: zni
real*8                  :: Te
real*8                  :: electron_mass ! allowing for heavy electrons
real*8                  :: collision_frequency, Z, CLei

Z    = 1.d0
CLei = 14.d0

collision_frequency = 4.d0 / 3.d0 * sqrt(TWOPI)* Z**2 * EL_CHG**4 * zni * CLei &
                    / ((2.d0 * TWOPI * EPS_ZERO)**2 * sqrt(electron_mass) * Te**1.5)

return
end


subroutine electron_ion_collision_frequency_v(zni, Te, v, electron_mass, collision_frequency)
  use constants
  real*8                  :: zni
  real*8                  :: Te
  real*8                  :: v 
  real*8                  :: electron_mass ! allowing for heavy electrons
  real*8                  :: collision_frequency, Z, CLei, v_th
  
  Z    = 1.d0
  CLei = 14.d0
  v_th = sqrt(Te / electron_mass)

  collision_frequency = zni * Z**2* EL_CHG**4 * CLei &
                      / (TWOPI * EPS_ZERO**2 * electron_mass**2 * V_th**3 )

  collision_frequency = 0.25d0 * collision_frequency * (v_th/v)**3
      
return
end

subroutine electron_ion_small_angle_scattering_orb5(electron, electron_mass, zni, Te, delta_t, ran1, alfa)
  use mod_particle_types
  use constants 
  type(particle_gc_vpar)  :: electron
  real*8                  :: electron_mass 
  real*8                  :: zni        ! ion density [m^-3]
  real*8                  :: Te         ! electron temperature [eV]
  real*8                  :: delta_t
  real*8                  :: ran1
  real*8                  :: alfa
  real*8                  :: collision_frequency, v_perp, v_norm
  real*8                  :: pitch_angle, pitch_out, delta_pitch, delta_theta
  
  v_perp = sqrt(2.d0 * electron%B_norm * electron%mu)
  
  v_norm = sqrt(v_perp**2 + electron%vpar**2) 
  
  pitch_angle = electron%vpar / v_norm

  ! Kauffmann, EUTERPE, 2010
  !pitch_angle = cos(atan(v_perp / electron%vpar))

  call electron_ion_collision_frequency_v(zni, Te, v_norm, electron_mass, collision_frequency)

! T. Vernay, S. Brunner, L. Villard, B. F. McMillan, S. Jolliet, T. M. Tran, A. Bottino and J. P. Graves
! PHYSICS OF PLASMAS 17, 122301 2010

!  delta_pitch = -2.d0 * pitch_angle * collision_frequency * delta_t &  
!              + ran1 * sqrt(2.d0*(1.d0 - pitch_angle**2)*collision_frequency*delta_t)
  
!  pitch_out = pitch_angle + delta_pitch 
!  electron%vpar = min(1.d0,pitch_out) * v_norm

! Kauffmann, EUTERPE 2010
! Lanti 2019, ran1 = random distribution mean 0, variance 1, alfa random uniform 0-1
  delta_theta = 2.d0 * ran1 * sqrt(collision_frequency*delta_t)

  electron%vpar = v_norm * (pitch_angle * cos(delta_theta) &

                          - sin(delta_theta) * sin(TWOPI*alfa) * sqrt(1.d0 - pitch_angle**2))
  
  electron%mu   = (v_norm**2 - electron%vpar**2) / (2.d0 * electron%B_norm)
  
  !if (electron%mu .lt. 0.d0) stop 'negative energies are not acceptable'
  
 
  return
end
  

subroutine electron_ion_small_angle_scattering_lu(electron, collision_frequency, delta_t, ran1)
use mod_particle_types
type(particle_gc_vpar)  :: electron
real*8                  :: delta_t
real*8                  :: ran1
real*8                  :: collision_frequency, v_par, v_perp, v_norm
real*8                  :: pitch_angle, pitch_out

v_par  = electron%vpar
v_perp = sqrt(2.d0 * electron%B_norm * electron%mu)

v_norm = sqrt(v_perp**2 + electron%vpar**2) 

pitch_angle = electron%vpar / v_norm

pitch_out = pitch_angle * (1.d0 - collision_frequency * delta_t) &
 
          + (ran1 - 0.5d0) * sqrt(12.d0*(1.d0 - pitch_angle**2) * collision_frequency * delta_t)

electron%vpar = min(1.d0,pitch_out) * v_norm
electron%mu   = (v_norm**2 - electron%vpar**2) / (2.d0 * electron%B_norm)

v_par  = electron%vpar
v_perp = sqrt(2.d0 * electron%B_norm * electron%mu)

!electron%vpar = min(1.d0,pitch_out) / sqrt(max(0.01d0,1.d0 - pitch_out**2)) * v_perp

return

end

subroutine loop_particle_gc_local(sim, i_group, n_orbit, ion_group, electron_group ,                        &
                                 jorek_feedback, project_profiles, project_density,                         &
                                 timesteps, n_steps, particle_start_time,                                   &
                                 lost, n_lost, n_lost_max,                                                  &
                                 energy_local, energy_lost_local, momentum_local, momentum_lost_local,      &
                                 potential_energy_local, potential_energy_lost_local, update)
use mod_parameters
use mod_particle_sim
use mod_normalisations
use mod_project_particles
use mod_random_seed
use mod_sampling
use mod_interp, only: sincosperiod_moivre, interp
use mod_basisfunctions
use mod_collisions
use mod_particle_types, only: particle_gc_vpar, particle_kinetic_leapfrog, copy_particle_kinetic_leapfrog
use mod_event
use mod_pcg32_rng
use mod_interp, only: mode_moivre, interp_RZ, interp_0
use mod_gc_variational, only : convert_gc_vpar_to_kinetic, copy_particle_gc_vpar, push_gc_rk4
use constants,   only: MU_ZERO, ATOMIC_MASS_UNIT, K_BOLTZ, EL_CHG
use phys_module, only: F0, CENTRAL_MASS, CENTRAL_DENSITY, ei_small_angle_scattering_Orb5, ei_small_angle_scattering_Lu
use omp_lib
implicit none

class(particle_sim), target, intent(inout)    :: sim
integer, intent(in)                           :: i_group                    ! the group of particles to be advanced
integer, intent(in)                           :: n_orbit                    ! the number of points on a gyro-orbit for orbit averaging (n_orbit=0 implies guiding centre)
integer, intent(in)                           :: ion_group, electron_group  ! the index of the group representing ions and electrons
type(projection), target, intent(inout)       :: jorek_feedback             ! the projections using the physicsmodel (i.e. Poisson equation)
type(projection), target, intent(in)          :: project_profiles           ! projection of profiles
type(projection), target, intent(in)          :: project_density            ! projection of density
real*8, intent(in)                            :: timesteps, particle_start_time ! the particle time step and starting time [s]
integer, intent(in)                           :: n_steps                    ! the number of time steps
logical, intent(in)                           :: update                     ! switch to update particle position (.true.) or calculate only rhs contribution to Poisson (not used currenly)
integer                                       :: lost(*)                    ! stores the index of lost particles
integer                                       :: n_lost                     ! the current number of lost particles in lost array
integer                                       :: n_lost_max                 ! the maximum number of lost particles to store (size(lost))
real*8                                        :: energy_local,    energy_lost_local   ! the energy and energy loss in this timestep [J]
real*8                                        :: momentum_local,  momentum_lost_local ! momentum and momentum loss in this timeste [kg m/s] (?)
real*8                                        :: potential_energy_local, potential_energy_lost_local ! potential energy and potential energy loss [J]

type(count_action)                            :: counter
type(particle_gc_vpar)                        :: particle_tmp
type(particle_kinetic_leapfrog), allocatable  :: p_orbit(:) 

real*8, allocatable :: feedback_rhs(:,:,:,:,:)

real*8    :: t, E(3), B(3), psi, U, zn0, zne0, zni0, Te_eV, Te0_eV, Ti0_eV, Te0
real*8    :: vpar0, vpar_i0, vpar_e0, rz_old(2), st_old(2)
real*8    :: v_temp(3), T_eV, K_eV, collision_frequency
real*8    :: ran1(1), ran2(2), ran3(3), ran6(6), R3(6)
real*8    :: v_par, v_par_sum, v_perp_sum, mu_sum, Wk_sum
real*8    :: n_b, m_b, v_b(3), q_heat(3), q_heat_par, kT_b, coulomb_log, coulomb_log_ii
real*8    :: v, v_s, v_t, v_R, v_Z, v_E, v_P
real*8    :: P2(2), P3(3), P6(6), R_g, Z_g, R_s, R_t, Z_s, Z_t, xjac
real*8    :: HHZ(n_tor), HHZ_p(n_tor), HH(4,4), HH_s(4,4), HH_t(4,4)
real*8    :: w0, w1, mmm(3)

integer(kind=1)     :: q_b
integer   :: i, j, k, l, m, i_elm_old, i_elm, ifail, i_coll, n_coll_ei
integer   :: seed, i_rng, n_stream, ierr, nthreads
integer   :: i_tor, index_lm, i_elm_temp, n_phases, n_valid_phases
logical   :: already_lost
real*8    :: energy_before, energy_after
real*8    :: va_par, va_perp, va_par_out, va_perp_out, delta_E_ab, delta_P_ab, E_ab, P_ab
type(pcg32_rng), dimension(:), allocatable :: rng

w0 = omp_get_wtime()

if (sim%my_id .eq. 0)                                            write(*,*) sim%my_id,' particle loop, group : ',i_group
if (i_group .eq. electron_group) then
  if ((sim%my_id .eq. 0) .and. (ei_small_angle_scattering_Orb5)) write(*,*) 'using small-angle scattering Orb5 model'
  if ((sim%my_id .eq. 0) .and. (ei_small_angle_scattering_Lu))   write(*,*) 'using small-angle scattering Lu model (not recommendedd)'
  if (ei_small_angle_scattering_Orb5 .and. ei_small_angle_scattering_Lu) then
    write(*,*) ' select at most one small-angle scattering model '
    stop 
  endif
endif

zn_norm  = CENTRAL_DENSITY * 1.d20                              ! (number) density normalisation
rho_norm = CENTRAL_MASS * ATOMIC_MASS_UNIT * zn_norm            ! rho_SI = rho_norm * rho
t_norm   = sqrt((MU_ZERO * rho_norm))                           ! t_SI   = t_norm * t_jorek
v_norm   = 1.d0 / t_norm                                        ! V_SI   = v_norm * v_jorek
E_norm   = 1.5d0 / MU_ZERO                                      ! E_SI   = E_norm * E_jorek
M_norm   = rho_norm * v_norm                                    ! momentum normalisation
Tev_norm = 1.d0 / (EL_CHG * MU_ZERO * zn_norm)                  ! T_ev [eV] = Tev_norm * T_jorek / 2.d0

allocate(feedback_rhs,source=jorek_feedback%rhs)

feedback_rhs  = 0.d0

energy_local                = 0.d0
energy_lost_local           = 0.d0
momentum_local              = 0.d0
momentum_lost_local         = 0.d0
potential_energy_local      = 0.d0
potential_energy_lost_local = 0.d0

n_phases = max(n_orbit,1)
allocate(p_orbit(n_phases))

!call with(sim, counter)

seed = random_seed()
n_stream = 1
!$ n_stream = omp_get_max_threads()
allocate(rng(n_stream))
do i=1,n_stream
  call rng(i)%initialize(6, seed, n_stream, i)
end do

select type (particles => sim%groups(i_group)%particles)
type is (particle_gc_vpar)

if (sim%my_id .eq. 0) write(*,*) 'starting loop gc : ',size(particles,1), n_orbit, n_phases

#ifdef __GFORTRAN__
  !$omp parallel do default(shared) & 
#else
  !$omp parallel do default(none) &
#endif
  !$omp schedule(dynamic,10)                                                     &
  !$omp shared(sim, particles, n_steps, timesteps, particle_start_time, update,  &
  !$omp rho_norm, t_norm, v_norm, E_norm, M_norm, zn_norm, Tev_norm, n_phases,   &
  !$omp n_orbit, central_density, central_mass, F0, i_group, electron_group,     &
  !$omp rng, lost, n_lost, n_lost_max, ion_group,                                &
  !$omp ei_small_angle_scattering_Orb5, ei_small_angle_scattering_Lu,            &
  !$omp project_profiles, project_density)                                       &
  !$omp private(particle_tmp, i_rng, i,j,k,l,m, t, E, B, psi, U, rz_old, st_old, &
  !$omp i_elm_old, i_elm, zn0, zne0, zni0, Te_eV, Te0_eV, Ti0_eV,                &
  !$omp p_orbit, ran1, ran2,                                                     & 
  !$omp P2, P3, P6, R_g, R_s, R_t, Z_g, Z_s, Z_t, xjac, HH, HH_s, HH_t, index_lm,&
  !$omp R3, ran3, ran6, n_b, m_b, q_b, kT_b, v_b, q_heat, q_heat_par, n_valid_phases,&
  !$omp v_par, v_par_sum, mu_sum, v_perp_sum, Wk_sum,                            &
  !$omp vpar_i0, vpar_e0, vpar0, i_coll, n_coll_ei, collision_frequency,         &
  !$omp va_par, va_perp, va_par_out, va_perp_out, v_E, v_P,                      &
  !$omp ifail, v, v_s, v_t, v_R, v_Z, HHZ, HHZ_p, already_lost, energy_before, energy_after)    &
  !$omp reduction(+:feedback_rhs,                                                &
  !$omp energy_lost_local, energy_local, momentum_lost_local, momentum_local,    &
  !$omp potential_energy_local, potential_energy_lost_local)
  do j=1,size(particles,1)

    call copy_particle_gc_vpar(particles(j),particle_tmp)

    !i_rng = 1
    i_rng = omp_get_thread_num()+1

    already_lost = .false.
    if (particle_tmp%i_elm .le. 0) already_lost = .true.

    call push_gc_rk4(sim%fields, particle_tmp, sim%groups(i_group)%mass, timesteps, n_steps, n_orbit) ! add B, P-orbit, and U to output of push_gc_rk4

!    call sim%fields%calc_EBpsiU(sim%time, particles(j)%i_elm, particles(j)%st, particles(j)%x(3), E, B, psi, U)

    if  (particle_tmp%i_elm .gt. 0)  then
      call copy_particle_gc_vpar(particle_tmp, particles(j))
    else
      call copy_particle_gc_vpar(particle_tmp, particles(j))
      if (.not. already_lost) then
        energy_lost_local           = energy_lost_local           + (0.5d0 * particle_tmp%vpar**2 + particle_tmp%mu * particle_tmp%B_norm) * particle_tmp%weight
        momentum_lost_local         = momentum_lost_local         + particle_tmp%vpar * particle_tmp%weight
!        potential_energy_lost_local = potential_energy_lost_local + 0.5d0 * EL_CHG * particles(j)%q * F0 * U * particles(j)%weight
        !$omp critical
        if (n_lost .lt. n_lost_max) then
          if (i_group .eq. 2) write(*,*) ' lost : ',j, n_lost,n_lost_max
          n_lost       = n_lost + 1   ! check for size overrun
          lost(n_lost) = j
        endif
        !$omp end critical
      endif
    endif

    if (particle_tmp%i_elm .le. 0) cycle

    energy_local   = energy_local   + (0.5d0 * particle_tmp%vpar**2 + particle_tmp%mu * particle_tmp%B_norm) * particle_tmp%weight
    momentum_local = momentum_local + particle_tmp%vpar * particle_tmp%weight

    call sim%fields%calc_EBpsiU(sim%time, particle_tmp%i_elm, particle_tmp%st, particle_tmp%x(3), E, B, psi, U)

    potential_energy_local = potential_energy_local + 0.5d0 * EL_CHG * particle_tmp%q * F0 * U * particle_tmp%weight  ! should be gyro-averaged

    call convert_gc_vpar_to_kinetic(sim%fields%node_list, sim%fields%element_list, particle_tmp, B, sim%groups(i_group)%mass, n_phases, p_orbit, ifail)
        
    if (ifail .lt. 0) cycle
        
    n_valid_phases = 0
    v_par_sum  = 0.d0
    v_perp_sum = 0.d0
    mu_sum     = 0.d0
    Wk_sum     = 0.d0

    do i=1, n_phases

      if (p_orbit(i)%i_elm .le. 0) cycle

      i_elm   = p_orbit(i)%i_elm

      call basisfunctions(p_orbit(i)%st(1), p_orbit(i)%st(2), HH, HH_s, HH_t)
  
      call mode_moivre(p_orbit(i)%x(3), HHZ)

      do m=1,3
        call interp(sim%fields%node_list, sim%fields%element_list, p_orbit(i)%i_elm, &
             4+m, 1, p_orbit(i)%st(1), p_orbit(i)%st(2), P3(m))
      enddo
      zn0    = max(0.01d0,P3(1)) * zn_norm
      Te_eV  = max(P3(2) * Tev_norm / 2.d0, 10d0)   ! Temperature in sim%fields%node_list is 2*Te (should change to Te,Ti in model600)
      vpar0  = P3(3) * V_norm ! to be checked, Should be multiplied with B?

! should the background values be 3D?
      do m=1,6
        call interp(project_profiles%node_list, sim%fields%element_list, p_orbit(i)%i_elm, &
             m, 1, p_orbit(i)%st(1), p_orbit(i)%st(2), P6(m))
      enddo
      zni0   = max(0.01d0*zn_norm, P6(1))
      zne0   = max(0.01d0*zn_norm, P6(2))
      Ti0_eV = max(P6(3) / P6(1) / EL_CHG, 10d0)
      Te0_eV = max(P6(4) / P6(2) / EL_CHG, 10d0)    ! includes energy conservation correction
      vpar_i0 = P6(5) / P6(1)
      vpar_e0 = P6(6) / P6(2)

      do l=1,n_vertex_max
        do m=1,n_order+1
  
          index_lm = (l-1)*(n_order+1) + m
  
          v   = HH(l,m) * sim%fields%element_list%element(i_elm)%size(l,m) * particle_tmp%weight * particle_tmp%q * t_norm / F0

          do i_tor=1,n_tor
            feedback_rhs(m,l,i_elm,i_tor,1) = feedback_rhs(m,l,i_elm,i_tor,1) + HHZ(i_tor) * v 
          enddo

        enddo   !< order
      enddo     !< vertex
    enddo       !< phases

    if (i_group .eq. electron_group) then

      if (ei_small_angle_scattering_Lu) then

        call electron_ion_collision_frequency(zne0, Te0_ev*EL_CHG, sim%groups(electron_group)%mass * ATOMIC_MASS_UNIT, collision_frequency)
        n_coll_ei = 1 + int(collision_frequency * timesteps / 0.15d0)

        energy_before = 0.5d0 * particles(j)%vpar**2 + particles(j)%mu * particles(j)%B_norm

        do i_coll = 1, n_coll_ei
          call rng(i_rng)%next(ran1)
          call  electron_ion_small_angle_scattering_Lu(particles(j), collision_frequency, timesteps/real(n_coll_ei,8), ran1(1))
        enddo

        energy_after =  0.5d0 * particles(j)%vpar**2 + particles(j)%mu * particles(j)%B_norm

      elseif (ei_small_angle_scattering_Orb5) then

        n_coll_ei = 10

        do i_coll = 1, n_coll_ei

          call rng(i_rng)%next(ran3)
          ran2 = boxmueller_transform(ran3(1:2))
          call electron_ion_small_angle_scattering_orb5(particles(j), sim%groups(electron_group)%mass * ATOMIC_MASS_UNIT, &
                                                        zne0, Te0_ev*EL_CHG, timesteps/real(n_coll_ei,8), ran2(1), ran3(3))
        enddo
      
        energy_after =  0.5d0 * particles(j)%vpar**2 + particles(j)%mu * particles(j)%B_norm

      endif ! small angle scattering model

    endif   ! electron group

  enddo     ! particles
  !$omp end parallel do
    
end select

t = particle_start_time + n_steps*timesteps

jorek_feedback%rhs  = jorek_feedback%rhs  + feedback_rhs /real(n_phases,8)
 
deallocate(feedback_rhs)

energy_local                = energy_local        * sim%groups(i_group)%mass * ATOMIC_MASS_UNIT
energy_lost_local           = energy_lost_local   * sim%groups(i_group)%mass * ATOMIC_MASS_UNIT 
momentum_local              = momentum_local      * sim%groups(i_group)%mass * ATOMIC_MASS_UNIT
momentum_lost_local         = momentum_lost_local * sim%groups(i_group)%mass * ATOMIC_MASS_UNIT 

w1 = omp_get_wtime()
mmm = mpi_minmeanmax(w1-w0)
if (sim%my_id .eq. 0) write(*,"(A,3f9.4,A)") " Particle stepping complete in ", mmm, "s"
  
end subroutine

end program

function f_density(n, P, grad_P) result(f)
  use domain
  use phys_module
  use equil_info
  integer, intent(in) :: n
  real*8,  intent(in) :: P(n), grad_P(3,n)
  real*8 :: psi, psi_axis, psi_bnd, Z, Z_xpoint, psi_n
  real*8 :: zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz
  real*4 :: f
  logical :: xpoint_local
  integer :: xcase_local

  psi      = P(n)
  psi_axis = ES%psi_axis
  psi_bnd  = ES%psi_bnd
  Z        = 0.d0
  Z_xpoint = -99.d0
  xpoint_local   = .false.
  xcase_local    = 1

  call density(xpoint_local, xcase_local, Z, Z_xpoint, psi, psi_axis, psi_bnd,&
               zn,dn_dpsi,dn_dz,dn_dpsi2,dn_dz2,dn_dpsi_dz,dn_dpsi3,dn_dpsi_dz2, dn_dpsi2_dz)

  psi_n = (psi - psi_axis)/(psi_bnd - psi_axis)

!  if ((psi_n .ge. psi_n_start) .and. (psi_n .le. psi_n_end)) then
    f = zn
!  else
!    f = 0.d0
!  endif

  !write(*,'(A,8e16.6)') 'density : ',psi,psi_axis, psi_bnd, psi_n,f

end function f_density

function f_ions(n, P, grad_P) result(f)
  use domain
  use phys_module
  use equil_info
  integer, intent(in) :: n
  real*8,  intent(in) :: P(n), grad_P(3,n)
  real*8 :: psi_n
  real*4 :: f

  psi_n = (P(2) - ES%psi_axis)/(ES%psi_bnd - ES%psi_axis)

  f = rho_1 + (rho_0-rho_1)*(1.d0 + rho_coef(1) * psi_n + rho_coef(2) * psi_n**2 + rho_coef(3) * psi_n**3) &
            * (0.5d0 - 0.5d0*tanh((psi_n - 1.00d0)/0.05d0))

end function f_ions

function f_electrons(n, P, grad_P) result(f)
  use domain
  use phys_module
  use equil_info
  integer, intent(in) :: n
  real*8,  intent(in) :: P(n), grad_P(3,n)
  real*8 :: psi_n
  real*4 :: f

  psi_n = (P(2) - ES%psi_axis)/(ES%psi_bnd - ES%psi_axis)

  f = rho_1 + (rho_0-rho_1)*(1.d0 + rho_coef(1) * psi_n + rho_coef(2) * psi_n**2 + rho_coef(3) * psi_n**3) &
            * (0.5d0 - 0.5d0*tanh((psi_n - 1.00d0)/0.04d0))

end function f_electrons

function T_ions(psi) result(f)
  use domain
  use phys_module, only : xpoint, xcase
  use constants,   only : MU_ZERO, EL_CHG
  use equil_info
  use mod_normalisations, only: Tev_norm
  implicit none
  real*8, intent(in) :: psi
  real*8 :: psi_axis, psi_bnd, Z, Z_xpoint(2), psi_n
  real*8 :: zTi,dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz
  real*8 :: f
  logical :: xpoint_local
  integer :: xcase_local
    
  psi_axis = ES%psi_axis
  psi_bnd  = ES%psi_bnd
  Z        = 0.d0  
  Z_xpoint = ES%Z_xpoint     
  xpoint_local = xpoint
  xcase_local  = xcase
    
  psi_n = (psi - ES%psi_axis)/(ES%psi_bnd - ES%psi_axis)

  call temperature_i(xpoint_local, xcase_local, Z, Z_xpoint, psi, psi_axis, psi_bnd, &
                     zTi, dTi_dpsi,dTi_dz,dTi_dpsi2,dTi_dz2,dTi_dpsi_dz,dTi_dpsi3,dTi_dpsi_dz2, dTi_dpsi2_dz)             
    
!  zTi = Ti_1 + (Ti_0-Ti_1)*(1.d0 + Ti_coef(1) * psi_n + Ti_coef(2) * psi_n**2 + Ti_coef(3) * psi_n**3) &
!            * (0.5d0 - 0.5d0*tanh((psi_n - Ti_coef(5))/Ti_coef(4)))

  f = zTi * TeV_norm

end function T_ions

function T_electrons(psi) result(f)
  use domain
  use phys_module, only : xcase, xpoint
  use constants,   only : MU_ZERO, EL_CHG
  use equil_info
  use mod_normalisations, only : Tev_norm
  implicit none
  real*8, intent(in) :: psi
  real*8 :: psi_axis, psi_bnd, Z, Z_xpoint(2), psi_n
  real*8 :: zTe,dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz
  real*8 :: f
  logical :: xpoint_local
  integer :: xcase_local
  
  psi_axis = ES%psi_axis
  psi_bnd  = ES%psi_bnd
  Z        = 0.d0 
  Z_xpoint = ES%Z_xpoint
  xpoint_local   = xpoint
  xcase_local    = xcase

  psi_n = (psi - ES%psi_axis)/(ES%psi_bnd - ES%psi_axis)

  call temperature_e(xpoint_local, xcase_local, Z, Z_xpoint, psi, psi_axis, psi_bnd, &
                     zTe, dTe_dpsi,dTe_dz,dTe_dpsi2,dTe_dz2,dTe_dpsi_dz,dTe_dpsi3,dTe_dpsi_dz2, dTe_dpsi2_dz)             

!  zTe = Te_1 + (Te_0-Te_1)*(1.d0 + Te_coef(1) * psi_n + Te_coef(2) * psi_n**2 + Te_coef(3) * psi_n**3) &
!            * (0.5d0 - 0.5d0*tanh((psi_n - Te_coef(5))/Te_coef(4)))
  
  f = zTe * TeV_norm
  
end function T_electrons

subroutine itg_energy(node_list,element_list,psi_start,psi_end,W_kin,W_tot)
!---------------------------------------------------------------
!
!---------------------------------------------------------------
use data_structure
use gauss
use basis_at_gaussian
use phys_module

implicit none

type (type_node_list)    :: node_list
type (type_element_list) :: element_list
type (type_element)      :: element
type (type_node)         :: nodes(n_vertex_max)

real*8     :: x_g(n_gauss,n_gauss),  x_s(n_gauss,n_gauss),  x_t(n_gauss,n_gauss)
real*8     :: y_s(n_gauss,n_gauss),  y_t(n_gauss,n_gauss)
real*8     :: eq_s(n_gauss,n_gauss), eq_t(n_gauss,n_gauss), psi_eq(n_gauss,n_gauss)
integer    :: i, j, k, in, ms, mt, iv, inode, ife, n_elements, i_tor
real*8     :: W_kin(n_tor), W_tot, xjac, BigR, wst
real*8     :: u0_x, u0_y, psi_start, psi_end, psi_norm, psi_axis, psi_bnd
integer, parameter :: ivar_psi=1, ivar_u=2, ivar_rho=5

W_kin = 0.d0
W_tot = 0.d0

#ifdef __GFORTRAN__
   !$omp parallel do default(shared) & ! workaround for Error: �__vtab_mod_pcg32_rng_Pcg32_rng� not specified in enclosing �parallel�
#else
   !$omp parallel do default(none) &
#endif
   !$omp schedule(dynamic,10)                                                         &
   !$omp shared(element_list, node_list, psi_start, psi_end, H, H_s, H_t)             &
   !$omp private(ife, element, nodes, iv, inode, x_g, x_s, x_t, y_s, y_t, eq_s, eq_t, &
   !$omp         i, j, ms, mt, in, wst, bigR, xjac, u0_x, u0_y, psi_eq )              &
   !$omp reduction(+:W_kin,W_tot)
do ife =1,  element_list%n_elements

  element = element_list%element(ife)

  do iv = 1, n_vertex_max
    inode     = element%vertex(iv)
    nodes(iv) = node_list%node(inode)
  enddo

  x_g(:,:)  = 0.d0; x_s(:,:) = 0.d0; x_t(:,:) = 0.d0;
  y_s(:,:)  = 0.d0; y_t(:,:) = 0.d0;
  eq_s(:,:) = 0.d0; eq_t(:,:) = 0.d0; psi_eq(:,:) = 0.d0

  do i=1,n_vertex_max
    do j=1,n_order+1
      do ms=1, n_gauss
        do mt=1, n_gauss

          x_g(ms,mt) = x_g(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H(i,j,ms,mt)
          x_s(ms,mt) = x_s(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_s(i,j,ms,mt)
          x_t(ms,mt) = x_t(ms,mt) + nodes(i)%x(1,j,1) * element%size(i,j) * H_t(i,j,ms,mt)
          y_s(ms,mt) = y_s(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_s(i,j,ms,mt)
          y_t(ms,mt) = y_t(ms,mt) + nodes(i)%x(1,j,2) * element%size(i,j) * H_t(i,j,ms,mt)

          psi_eq(ms,mt) = psi_eq(ms,mt) + nodes(i)%values(1,j,ivar_psi) * element%size(i,j) * H(i,j,ms,mt)
          
        enddo
      enddo
    enddo
  enddo

  do in=1,n_tor

    eq_s(:,:) = 0.d0; eq_t(:,:) = 0.d0

    do ms=1, n_gauss
      do mt=1, n_gauss

        do i=1,n_vertex_max
          do j=1,n_order+1

            eq_s(ms,mt)  = eq_s(ms,mt)  + nodes(i)%values(in,j,ivar_u) * element%size(i,j) * H_s(i,j,ms,mt)
            eq_t(ms,mt)  = eq_t(ms,mt)  + nodes(i)%values(in,j,ivar_u) * element%size(i,j) * H_t(i,j,ms,mt)
        
          enddo
        enddo

      enddo
    enddo

    do ms=1, n_gauss
      do mt=1, n_gauss
          
        if (.not. ((psi_eq(ms,mt) .gt. psi_start) .and. (psi_eq(ms,mt) .lt. psi_end))) cycle

        wst = wgauss(ms)*wgauss(mt)

        xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        BigR = x_g(ms,mt)

        u0_x  = (   y_t(ms,mt) * eq_s(ms,mt) - y_s(ms,mt) * eq_t(ms,mt) ) / xjac
        u0_y  = ( - x_t(ms,mt) * eq_s(ms,mt) + x_s(ms,mt) * eq_t(ms,mt) ) / xjac

        W_kin(in) = W_kin(in) + (u0_x*u0_x   + u0_y*u0_y) * BigR**3 * xjac * wst
        W_tot     = W_tot     + (u0_x*u0_x   + u0_y*u0_y) * BigR    * xjac * wst        

      enddo
    enddo

  enddo  ! n_tor loop

enddo    ! elements loop
!$omp end parallel do

do in=1,n_tor
  if (mode(in) .ne. 0) then
    W_kin(in) = 0.5d0 * W_kin(in)
  endif
enddo

return
end

subroutine density_integral(node_list,element_list,f_density,psi_start,psi_end,total_particles,total_volume)
  !---------------------------------------------------------------
  !
  !---------------------------------------------------------------
  use data_structure
  use gauss
  use basis_at_gaussian
  use phys_module
  
  implicit none
  
  type (type_node_list)    :: node_list
  type (type_element_list) :: element_list
  real*8, intent(in)       :: psi_start, psi_end
  real*8, intent(out)      :: total_particles, total_volume
  real*4, external         :: f_density
  
  real*8     :: x_g(n_gauss,n_gauss),  x_s(n_gauss,n_gauss),  x_t(n_gauss,n_gauss)
  real*8     :: y_s(n_gauss,n_gauss),  y_t(n_gauss,n_gauss)
  real*8     :: eq_g(n_gauss,n_gauss), psi_eq(n_gauss,n_gauss)
  integer    :: i, j, ms, mt, iv, inode, ife, n_elements
  real*8     :: xjac, BigR, wst
  real*8     :: sum_particles, sum_volume
  integer, parameter :: ivar_psi=1, ivar_rho=5
  
 sum_particles = 0.d0
 sum_volume    = 0.d0
  
!$omp parallel do default(none) &
!$omp schedule(dynamic,10)                                                &
!$omp shared(element_list, node_list, psi_start, psi_end, H, H_s, H_t)    &
!$omp private(ife, iv, inode, x_g, x_s, x_t, y_s, y_t, eq_g,              &
!$omp         i, j, ms, mt, wst, bigR, xjac, psi_eq )                     &
!$omp reduction(+:sum_particles, sum_volume)
  do ife =1,  element_list%n_elements
    
    x_g(:,:)  = 0.d0;   x_s(:,:)    = 0.d0;   x_t(:,:) = 0.d0;
    y_s(:,:)  = 0.d0;   y_t(:,:)    = 0.d0;
    eq_g(:,:) = 0.d0;   psi_eq(:,:) = 0.d0
  
    do i=1,n_vertex_max

      inode = element_list%element(ife)%vertex(i)

      do j=1,n_order+1
        do ms=1, n_gauss
          do mt=1, n_gauss
  
            x_g(ms,mt) = x_g(ms,mt) + node_list%node(inode)%x(1,j,1) * element_list%element(ife)%size(i,j) * H(i,j,ms,mt)
            x_s(ms,mt) = x_s(ms,mt) + node_list%node(inode)%x(1,j,1) * element_list%element(ife)%size(i,j) * H_s(i,j,ms,mt)
            x_t(ms,mt) = x_t(ms,mt) + node_list%node(inode)%x(1,j,1) * element_list%element(ife)%size(i,j) * H_t(i,j,ms,mt)
            y_s(ms,mt) = y_s(ms,mt) + node_list%node(inode)%x(1,j,2) * element_list%element(ife)%size(i,j) * H_s(i,j,ms,mt)
            y_t(ms,mt) = y_t(ms,mt) + node_list%node(inode)%x(1,j,2) * element_list%element(ife)%size(i,j) * H_t(i,j,ms,mt)
  
            psi_eq(ms,mt) = psi_eq(ms,mt) + node_list%node(inode)%values(1,j,ivar_psi) * element_list%element(ife)%size(i,j) * H(i,j,ms,mt)
            eq_g(ms,mt)   = eq_g(ms,mt)   + node_list%node(inode)%values(1,j,ivar_rho) * element_list%element(ife)%size(i,j) * H(i,j,ms,mt)
            
          enddo
        enddo
      enddo
    enddo

    do ms=1, n_gauss
      do mt=1, n_gauss
            
!          if (.not. ((psi_eq(ms,mt) .gt. psi_start) .and. (psi_eq(ms,mt) .lt. psi_end))) cycle
  
        wst = wgauss(ms)*wgauss(mt)
  
        xjac = x_s(ms,mt)*y_t(ms,mt) - x_t(ms,mt)*y_s(ms,mt)
        BigR = x_g(ms,mt)

        sum_particles = sum_particles + BigR * xjac * wst * f_density(1, psi_eq(ms,mt), (/0.d0, 0.d0,0.d0/))    
        sum_volume    = sum_volume    + BigR * xjac * wst   
  
      enddo
    enddo

  enddo    ! elements loop
  !$omp end parallel do
  
  total_volume    = sum_volume    * TWOPI
  total_particles = sum_particles * TWOPI

  !write(*,'(A,4e14.6)') 'DENSITY INTEGRAL : ',total_particles, total_volume
  return
  end


subroutine grid_reduction(node_list,element_list,psi_n_start,psi_n_end,n_tht)
!----------------------------------------------------------------------
! select a part of the grid between psi_start and psi_end
!----------------------------------------------------------------------
use mod_settings
use data_structure, only: type_bnd_element_list, type_bnd_node_list 
use mod_neighbours
use equil_info

type(type_node_list)    :: node_list
type(type_element_list) :: element_list
  
real*8  :: psi_n_start, psi_n_end
integer :: n_tht
integer :: ivar_psi = 1
  
type(type_node_list)    :: newnode_list
type(type_element_list) :: newelement_list
integer                 :: n_remove_elements, n_remove_nodes, skip_index, inode, ielm, iv, k
integer                 :: remove_elements(n_nodes_max), remove_nodes(n_nodes_max), newnode_index(n_nodes_max)
integer                 :: i, node_start, node_end
real*8                  :: psi_n
logical                 :: add_node
  
node_start = 1
node_end   = node_list%n_nodes

do inode = 1, node_list%n_nodes
  psi_n = (node_list%node(inode)%values(1,1,1) - ES%psi_axis)/(ES%psi_bnd - ES%psi_axis)
  if (psi_n .lt. psi_n_start) node_start = inode + 1
  if (psi_n .lt. psi_n_end)   node_end   = inode
enddo
  
newnode_list    = node_list
newelement_list = element_list
  
n_remove_nodes  = node_start - 1 + node_list%n_nodes - node_end
  
do inode=1, node_start-1
  remove_nodes(inode) = inode
enddo
do inode=node_start, n_remove_nodes
  remove_nodes(inode) = node_end + inode - node_start + 1
enddo
  
write(*,*) ' number of nodes to be removed : ',n_remove_nodes, node_start, node_end
  
n_remove_elements = 0
  
do ielm=1, newelement_list%n_elements
  
  do iv=1, 4
    if ( any(remove_nodes(1:n_remove_nodes) == newelement_list%element(ielm)%vertex(iv))) then
      remove_elements(n_remove_elements+1) = ielm
      n_remove_elements = n_remove_elements + 1
      exit
    endif
  enddo
enddo
  
write(*,*) ' number of elements to be removed : ',n_remove_elements
  
!---------------------------- copy new grid into nodes/elements, optional remove some nodes/elements
node_list%n_nodes = 0
skip_index        = 0
  
do inode = 1, newnode_list%n_nodes
  
  psi_n = (newnode_list%node(inode)%values(1,1,1) - ES%psi_axis)/(ES%psi_bnd - ES%psi_axis)
  
  do i=1, n_remove_nodes
    add_node = .true.
    if (remove_nodes(i) == inode) then
      add_node = .false.
      exit
    endif
  enddo

  if (add_node) then
    newnode_index(inode)                      = node_list%n_nodes+1 
    node_list%node(node_list%n_nodes+1)       = newnode_list%node(inode) 
!      node_list%node(node_list%n_nodes+1)%index = newnode_list%node(inode)%index - skip_index
    node_list%n_nodes = node_list%n_nodes + 1
  else 
    skip_index = skip_index + newnode_list%node(inode)%index(4) - newnode_list%node(inode)%index(1) + 1
  endif 
enddo

write(*,*) 'new nodes : ',node_list%n_nodes,', expected nodes : ',newnode_list%n_nodes - n_remove_nodes

!--------------------------------------- index skipping above is not correct yet (due to axis indexes increasing by 3 instead of 4)
do inode = 1, node_list%n_nodes
  do k = 1, 4
    node_list%node(inode)%index(k)  = 4*(inode-1) + k
    node_list%node(inode)%axis_node = .false.
  enddo
enddo
  
element_list%n_elements = 0
do ielm = 1, newelement_list%n_elements
  if (.not. any(remove_elements(1:n_remove_elements) == ielm)) then
    element_list%element(element_list%n_elements+1) = newelement_list%element(ielm)    
    element_list%n_elements                         = element_list%n_elements + 1
  endif
enddo
  
do ielm = 1, element_list%n_elements
  do iv=1, n_vertex_max
    inode = element_list%element(ielm)%vertex(iv)
    element_list%element(ielm)%vertex(iv) = newnode_index(inode)  
  enddo
enddo
  
do inode = 1, n_tht
  node_list%node(inode)%boundary                     = 2
  node_list%node(node_list%n_nodes-inode+1)%boundary = 2
enddo

call update_neighbours(node_list,element_list, force_rtree_initialize=.true.)

write(*,'(A,4i8)') 'done grid_reduction :',node_list%n_nodes, element_list%n_elements,node_list%node(node_list%n_nodes)%index(4)
end
