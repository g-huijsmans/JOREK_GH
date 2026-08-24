program test_projection_diagnostics
  use constants, only: atomic_mass_unit
  use mod_particle_sim, only: particle_sim
  use mod_particle_types, only: particle_kinetic, particle_kinetic_leapfrog, &
                                particle_gc_vpar, particle_gc_Qin
  use mod_project_particles, only: projection
  use mod_rhs_projections, only: proj_Pressure, proj_vpar, kinetic_vpar_from_B
  implicit none

  real*8, parameter :: tol = 64.d0*epsilon(1.d0)
  type(particle_sim) :: sim
  type(particle_kinetic) :: kinetic
  type(particle_kinetic_leapfrog) :: leapfrog
  type(particle_gc_vpar) :: gc
  type(particle_gc_Qin) :: qin
  type(projection) :: projected_field
  real*8 :: expected, actual, max_abs_difference, max_rel_difference
  real*8 :: rhs_before(2), nodal_before

  max_abs_difference = 0.d0
  max_rel_difference = 0.d0
  allocate(sim%groups(1))
  sim%groups(1)%mass = 2.5d0

  ! Historical pressure formulas for centred and leapfrog full-orbit particles.
  kinetic%v = [3.d0, -4.d0, 12.d0]
  kinetic%weight = 7.d0
  expected = sim%groups(1)%mass*atomic_mass_unit*169.d0/3.d0
  call check_value('proj_Pressure particle_kinetic', &
                   proj_Pressure(sim,1,kinetic), expected)

  leapfrog%v = [-6.d0, 2.d0, -3.d0]
  leapfrog%weight = 11.d0
  expected = sim%groups(1)%mass*atomic_mass_unit*49.d0/3.d0
  call check_value('proj_Pressure particle_kinetic_leapfrog', &
                   proj_Pressure(sim,1,leapfrog), expected)

  ! mu has velocity-moment units here: 2*mu*B is v_perp squared.
  gc%vpar = -5.d0
  gc%mu = 0.75d0
  gc%B_norm = 4.d0
  gc%weight = 13.d0
  gc%i_elm = 1
  expected = sim%groups(1)%mass*atomic_mass_unit*31.d0/3.d0
  call check_value('proj_Pressure particle_gc_vpar', &
                   proj_Pressure(sim,1,gc), expected)

  qin%vpar = 8.d0
  qin%mu = 1.25d0
  qin%B_norm = 99.d0 ! deliberately stale: the Qin pusher does not update it
  qin%Bn_k = 2.d0
  qin%weight = 17.d0
  qin%i_elm = 1
  expected = sim%groups(1)%mass*atomic_mass_unit*69.d0/3.d0
  call check_value('proj_Pressure particle_gc_Qin', &
                   proj_Pressure(sim,1,qin), expected)

  ! The transformation is per physical particle; sample_rhs applies weight.
  expected = kinetic%weight*(sim%groups(1)%mass*atomic_mass_unit*169.d0/3.d0) + &
             gc%weight*(sim%groups(1)%mass*atomic_mass_unit*31.d0/3.d0)
  actual = kinetic%weight*proj_Pressure(sim,1,kinetic) + &
           gc%weight*proj_Pressure(sim,1,gc)
  call check_value('weighted pressure projection RHS contribution',actual,expected)

  ! Full-orbit v_parallel is v dot B/|B|, including sign and zero.
  call check_value('kinetic vpar positive', &
                   kinetic_vpar_from_B([1.d0,2.d0,3.d0],[0.d0,0.d0,2.d0]),3.d0)
  call check_value('kinetic vpar negative', &
                   kinetic_vpar_from_B([1.d0,2.d0,-3.d0],[0.d0,0.d0,2.d0]),-3.d0)
  call check_value('kinetic vpar zero', &
                   kinetic_vpar_from_B([1.d0,2.d0,0.d0],[0.d0,0.d0,2.d0]),0.d0)

  ! Both GC representations store physical v_parallel directly; Qin needs no
  ! conversion from Astar_k for this diagnostic.
  gc%vpar = 9.d0
  call check_value('gc_vpar positive',proj_vpar(sim,1,gc),9.d0)
  gc%vpar = -9.d0
  call check_value('gc_vpar negative',proj_vpar(sim,1,gc),-9.d0)
  gc%vpar = 0.d0
  call check_value('gc_vpar zero',proj_vpar(sim,1,gc),0.d0)

  qin%Astar_k = [101.d0,102.d0,103.d0]
  qin%vpar = 6.d0
  call check_value('Qin vpar positive',proj_vpar(sim,1,qin),6.d0)
  qin%vpar = -6.d0
  call check_value('Qin vpar negative',proj_vpar(sim,1,qin),-6.d0)
  qin%vpar = 0.d0
  call check_value('Qin vpar zero',proj_vpar(sim,1,qin),0.d0)

  ! set_vtk_active controls output eligibility only.  It leaves both RHS and
  ! already-projected nodal storage bit-for-bit unchanged.
  allocate(projected_field%rhs(1,1,1,1,2))
  projected_field%rhs(1,1,1,1,:) = [2.5d0,-7.5d0]
  rhs_before = projected_field%rhs(1,1,1,1,:)
  allocate(projected_field%node_list)
  projected_field%node_list%n_nodes = 1
  allocate(projected_field%node_list%node(1))
  allocate(projected_field%node_list%node(1)%values(1,1,1))
  projected_field%node_list%node(1)%values(1,1,1) = 42.d0
  nodal_before = projected_field%node_list%node(1)%values(1,1,1)

  if (.not. projected_field%vtk_active) error stop 'VTK output must be active by default'
  call projected_field%set_vtk_active(.false.)
  if (projected_field%vtk_active) error stop 'set_vtk_active(.false.) did not deactivate output'
  if (any(projected_field%rhs(1,1,1,1,:) /= rhs_before)) &
    error stop 'VTK deactivation changed the projection RHS'
  if (projected_field%node_list%node(1)%values(1,1,1) /= nodal_before) &
    error stop 'VTK deactivation changed projected nodal values'

  call projected_field%set_vtk_active(.true.)
  if (.not. projected_field%vtk_active) error stop 'set_vtk_active(.true.) did not activate output'
  if (any(projected_field%rhs(1,1,1,1,:) /= rhs_before)) &
    error stop 'VTK activation changed the projection RHS'
  if (projected_field%node_list%node(1)%values(1,1,1) /= nodal_before) &
    error stop 'VTK activation changed projected nodal values'

  write(*,'(A,ES12.4)') 'projection diagnostics max absolute difference: ',max_abs_difference
  write(*,'(A,ES12.4)') 'projection diagnostics max relative difference: ',max_rel_difference
  write(*,'(A)') 'projection diagnostics regression: PASS'

contains

  subroutine check_value(label, value, reference)
    character(len=*), intent(in) :: label
    real*8, intent(in) :: value, reference
    real*8 :: abs_difference, rel_difference, scale

    abs_difference = abs(value-reference)
    scale = max(abs(reference),tiny(1.d0))
    rel_difference = abs_difference/scale
    max_abs_difference = max(max_abs_difference,abs_difference)
    max_rel_difference = max(max_rel_difference,rel_difference)
    if (abs_difference .gt. tol*scale) then
      write(*,'(A)') trim(label)//': FAIL'
      write(*,'(A,ES24.16)') '  current   = ',value
      write(*,'(A,ES24.16)') '  reference = ',reference
      error stop 'projection diagnostic mismatch'
    end if
  end subroutine check_value

end program test_projection_diagnostics
