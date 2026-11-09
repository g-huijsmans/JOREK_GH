!> Bounded verification of Cartesian ballistic motion through a cubic R,Z mesh.
!> No particle arrays, RNG, MPI, slot pools or evolution state.
module mod_gk_neutral_trajectory
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use data_structure, only: type_node_list, type_element_list
use mod_interp, only: interp_RZ
use mod_boundary, only: get_st_on_bnd
implicit none
private
public :: verify_neutral_trajectory, side_correction, boundary_sides
public :: validate_neutral_start_coordinates
public :: WALL_CROSSING, INTERIOR_ENDPOINT, UNRESOLVED
real*8, parameter, public :: wall_intersection_tolerance=1.d-6
integer, parameter :: WALL_CROSSING=1, INTERIOR_ENDPOINT=2, UNRESOLVED=3
real*8, parameter :: st_tolerance=1.d-10, roundoff_tolerance=128.d0*epsilon(1.d0)
integer, parameter, public :: max_trace_elements=128
type :: trajectory_origin
  integer :: element=0
  real*8 :: x_old(3),v_old(3),st_old(2),dt
end type
type :: trajectory_point
  real*8 :: x(3),v(3),st(2)=0.d0
end type
! Geometry evidence; contains no particle, RNG or transaction state.
type :: side_correction
  integer :: status=UNRESOLVED, n_visits=0
  logical :: fallback_attempted=.false.,fallback_used=.false.
  integer :: fallback_faces=0,fallback_tangencies=0
  logical :: coalescible=.false.
  integer :: local_clusters=0
  real*8 :: cluster_tolerance=0.d0,cluster_extent=0.d0
  real*8 :: max_residual=0.d0
  ! Transient evidence: attempts, failures, total iterations, maximum iterations.
  integer :: refinement_counts(4)=0
  real*8 :: refinement_residuals(2)=0.d0
  real*8 :: trace_u(max_trace_elements)=0.d0
  character(len=160) :: reason=""
  logical :: corrected=.false., start_clamped=.false.
  integer :: original_side=0, side=0, original_element=0, element=0
  integer :: examined=0, nonexit=0, no_exit=0, outside=0, n_trace=0
  integer :: trace_elements(max_trace_elements)=0, trace_faces(max_trace_elements)=0
  real*8 :: trace_fractions(max_trace_elements)=0.d0
  real*8 :: original_st(2)=0.d0, st(2)=0.d0, wall(2)=0.d0, normal(3)=0.d0
  real*8 :: impact_x(3)=0.d0, impact_v(3)=0.d0
  real*8 :: original_vn=0.d0, vn=0.d0, lambda=0.d0
end type

contains
subroutine verify_neutral_trajectory(nodes,elements,start_element,start_st,x_old,v_old,dt,hit,endpoint_x,endpoint_v)
  type(type_node_list), intent(in) :: nodes
  type(type_element_list), intent(in) :: elements
  integer, intent(in) :: start_element
  real*8, intent(in) :: start_st(2),x_old(3),v_old(3),dt
  type(side_correction), intent(out) :: hit
  real*8, optional, intent(in) :: endpoint_x(3),endpoint_v(3)
  type(trajectory_origin) :: origin
  type(trajectory_point) :: endpoint
  integer :: error
  real*8 :: residual
  logical :: clamped,ok
  hit%reason='nonfinite state or nonpositive timestep'
  if(.not.all(ieee_is_finite([start_st,x_old,v_old,dt]))) return
  if(dt<=0.d0) return
  origin=trajectory_origin(start_element,x_old,v_old,start_st,dt)
  ! Only the local origin may be clipped, after independent FE validation.
  call validate_start_roundoff(nodes,elements,origin,clamped,residual,ok,hit%reason)
  if(.not.ok) return
  endpoint%x=x_old;endpoint%v=v_old
  call push_trajectory(endpoint,dt)
  if(present(endpoint_x)) endpoint%x=endpoint_x
  if(present(endpoint_v)) endpoint%v=endpoint_v
  call trace_trajectory(nodes,elements,origin,endpoint,dt,hit,error)
  hit%start_clamped=clamped
  hit%max_residual=max(hit%max_residual,residual)
  if(error/=0) then
    hit%status=UNRESOLVED
    select case(trim(hit%reason))
    case('outward crossing cannot be separated from roundoff/grazing', &
         'outward crossing does not progress beyond the element entry', &
         'earliest outward crossing is a corner or has an unresolved tie', &
         'unresolved tangency before any verified outward exit', &
         'coincident/grazing or unresolved polynomial intersection')
      call boundary_fallback(nodes,elements,origin,endpoint,dt,hit)
    end select
  endif
  if(hit%status==WALL_CROSSING)then
    ! Physical first-contact evidence, separate from the existing transport
    ! convention which discards unused flight and retains the endpoint angle.
    endpoint%x=origin%x_old;endpoint%v=origin%v_old
    call push_trajectory(endpoint,hit%lambda*dt)
    hit%impact_x=endpoint%x;hit%impact_x(:2)=hit%wall
    hit%impact_v=endpoint%v
  endif
end subroutine

! Excursions of at most 32 ulps of the unit reference cell are candidates only.
! Require both the original and clipped FE maps to reproduce the physical start
! to a precision-scaled tolerance, with positive Jacobians. No travel/root tolerance
! enters this check; the particle's position, velocity and saved provenance stay exact.
subroutine validate_start_roundoff(nodes,elements,origin,clamped,residual,ok,reason)
  type(type_node_list), intent(in) :: nodes
  type(type_element_list), intent(in) :: elements
  type(trajectory_origin), intent(inout) :: origin
  logical, intent(out) :: clamped,ok
  real*8, intent(out) :: residual
  character(len=*), intent(out) :: reason
  real*8, parameter :: coordinate_roundoff=32.d0*epsilon(1.d0)
  real*8 :: st(2),y(2),rs,rt,zs,zt,jac,tol
  integer :: k,e
  clamped=.false.;ok=.false.;residual=0.d0
  reason='starting reference coordinates outside roundoff interval'
  if(any(origin%st_old < -coordinate_roundoff).or.any(origin%st_old > 1.d0+coordinate_roundoff)) return
  if(all(origin%st_old>=0.d0).and.all(origin%st_old<=1.d0)) then
    ok=.true.
    return
  endif
  e=origin%element
  reason='invalid starting element/connectivity for coordinate roundoff check'
  if(e<1.or.e>elements%n_elements) return
  if(any(elements%element(e)%vertex<1).or.any(elements%element(e)%vertex>nodes%n_nodes)) return
  st=origin%st_old
  do k=1,2
    call interp_RZ(nodes,elements,e,st(1),st(2),y(1),rs,rt,y(2),zs,zt)
    jac=rs*zt-rt*zs
    reason='invalid starting FE geometry for coordinate roundoff check'
    if(.not.all(ieee_is_finite([y,rs,rt,zs,zt,jac])).or.jac<=0.d0.or.y(1)<=0.d0) return
    tol=128.d0*epsilon(1.d0)*max(1.d0,norm2(origin%x_old(:2)),norm2(y),norm2([rs,zs]),norm2([rt,zt]))
    residual=max(residual,norm2(y-origin%x_old(:2)))
    reason='starting position inconsistent with roundoff-clipped FE coordinates'
    if(residual>tol) return
    st=max(0.d0,min(1.d0,st))
  enddo
  origin%st_old=st;clamped=.true.;ok=.true.;reason=''
end subroutine

! Reuse the SAME start-roundoff validation for a staged reaction transaction.
! Only local reference coordinates are returned; no physical state is moved.
subroutine validate_neutral_start_coordinates(nodes,elements,element,x,st,checked_st,ok)
  type(type_node_list), intent(in) :: nodes
  type(type_element_list), intent(in) :: elements
  integer, intent(in) :: element
  real*8, intent(in) :: x(3),st(2)
  real*8, intent(out) :: checked_st(2)
  logical, intent(out) :: ok
  type(trajectory_origin) :: origin
  logical :: clamped
  real*8 :: residual
  character(len=160) :: reason
  checked_st=st;ok=.false.
  if(element<1.or.element>elements%n_elements) return
  if(.not.all(ieee_is_finite([x,st]))) return
  origin=trajectory_origin(element,x,[0.d0,0.d0,0.d0],st,0.d0)
  call validate_start_roundoff(nodes,elements,origin,clamped,residual,ok,reason)
  if(ok) checked_st=origin%st_old
end subroutine

! Direct inverse in the final traversed element, independent of the shared locator.
! Nine bounded Newton starts; require positive Jacobian, a unique interior solution,
! and reconstruction after any roundoff clipping. Failed inversion never recovers.
subroutine invert_endpoint(nodes,elements,e,x,travel,answer,residual,ok)
  type(type_node_list), intent(in) :: nodes
  type(type_element_list), intent(in) :: elements
  integer, intent(in) :: e
  real*8, intent(in) :: x(2),travel
  real*8, intent(out) :: answer(2),residual
  logical, intent(out) :: ok
  real*8 :: st(2),trial(2),y(2),yt(2),delta(2),rs,rt,zs,zt,jac,err,errt,alpha,tol
  integer :: a,b,it,ls
  logical :: found
  answer=0.d0;residual=huge(1.d0);ok=.false.;found=.false.
  tol=512.d0*epsilon(1.d0)*max(1.d0,norm2(x))+1.d-10*travel
  do a=1,3
    do b=1,3
      st=[real(a,8),real(b,8)]*0.25d0
      do it=1,64
        call interp_RZ(nodes,elements,e,st(1),st(2),y(1),rs,rt,y(2),zs,zt)
        jac=rs*zt-rt*zs
        if(.not.all(ieee_is_finite([y,rs,rt,zs,zt,jac])).or.jac<=0.d0) exit
        err=norm2(x-y)
        if(err<=8.d0*epsilon(1.d0)*max(1.d0,norm2(x))) then
          if(any(st < -st_tolerance).or.any(st>1.d0+st_tolerance)) exit
          st=max(0.d0,min(1.d0,st))
          call interp_RZ(nodes,elements,e,st(1),st(2),y(1),y(2))
          err=norm2(x-y)
          if(err>tol) exit
          if(found) then
            if(norm2(answer-st)>1.d-8) return
          endif
          if(err<residual) then
            answer=st;residual=err
          endif
          found=.true.
          exit
        endif
        delta=[zt*(x(1)-y(1))-rt*(x(2)-y(2)),rs*(x(2)-y(2))-zs*(x(1)-y(1))]/jac
        if(.not.all(ieee_is_finite(delta))) exit
        alpha=1.d0
        do ls=1,32
          trial=st+alpha*delta
          if(all(trial>=-st_tolerance).and.all(trial<=1.d0+st_tolerance)) then
            call interp_RZ(nodes,elements,e,trial(1),trial(2),yt(1),yt(2))
            errt=norm2(x-yt)
            if(ieee_is_finite(errt)) then
              if(errt<err) exit
            endif
          endif
          alpha=0.5d0*alpha
        enddo
        if(ls>32) exit
        st=trial
      enddo
    enddo
  enddo
  ok=found
end subroutine

pure subroutine push_trajectory(p, dt)
  type(trajectory_point), intent(inout) :: p
  real*8, intent(in) :: dt
  real*8 :: r_local, phi_local, radius, angle, c, s, vr, vphi

  r_local = p%x(1) + dt*p%v(1)
  phi_local = dt*p%v(3)
  radius = sqrt(r_local**2 + phi_local**2)
  angle = 0.d0
  ! At R=0 the basis is arbitrary; retaining the old basis avoids division by zero.
  if (radius > 0.d0) angle = atan2(phi_local, r_local)
  c = cos(angle)
  s = sin(angle)
  vr = p%v(1)
  vphi = p%v(3)
  p%x = [radius, p%x(2) + dt*p%v(2), p%x(3) + angle]
  p%v(1) = c*vr + s*vphi
  p%v(3) = -s*vr + c*vphi
end subroutine push_trajectory

! Adjacent boundary-node criterion shared in convention with boundary_from_grid.
function boundary_sides(nodes,elements,e) result(sides)
  type(type_node_list), intent(in) :: nodes
  type(type_element_list), intent(in) :: elements
  integer, intent(in) :: e
  logical :: sides(4)
  integer :: side,n1,n2
  do side=1,4
    n1=elements%element(e)%vertex(side)
    n2=elements%element(e)%vertex(mod(side,4)+1)
    sides(side)=nodes%node(n1)%boundary>0.and.nodes%node(n2)%boundary>0
  enddo
end function

!> Explicit-side normal: same positive-Jacobian convention as wall_normal_vector.
subroutine side_geometry(nodes,elements,e,side,u,st,wall,tangent,normal,ok)
  type(type_node_list), intent(in) :: nodes
  type(type_element_list), intent(in) :: elements
  integer, intent(in) :: e,side
  real*8, intent(in) :: u
  real*8, intent(out) :: st(2),wall(2),tangent(2),normal(3)
  logical, intent(out) :: ok
  real*8 :: rs,rt,zs,zt,jac,length
  logical :: s_const
  st=0.d0; s_const=.false.
  call get_st_on_bnd(u,side,st(1),st(2),s_const)
  call interp_RZ(nodes,elements,e,st(1),st(2),wall(1),rs,rt,wall(2),zs,zt)
  tangent=[rs,zs]; if(s_const) tangent=[rt,zt]
  select case(side)
  case(1); normal=[-zs,rs,0.d0]
  case(2); normal=[-zt,rt,0.d0]
  case(3); normal=[zs,-rs,0.d0]
  case(4); normal=[zt,-rt,0.d0]
  end select
  jac=rs*zt-rt*zs; length=norm2(normal)
  ok=all(ieee_is_finite([wall,tangent,normal,jac,length]))
  if(.not.ok) return
  ok=wall(1)>0.d0.and.jac>0.d0.and.length>0.d0
  if(ok) normal=normal/length
end subroutine

pure function polynomial(c,n,x) result(y)
  integer, intent(in) :: n
  real*8, intent(in) :: c(0:),x
  real*8 :: y
  integer :: i
  y=c(n)
  do i=n-1,0,-1
    y=y*x+c(i)
  enddo
end function

!> Isolate all roots on [0,1] using derivative roots to partition monotone intervals.
!> Degree at most six here. Bisection never leaves the bracket; near-stationary
!> candidates are retained for geometric validation, not silently discarded.
recursive subroutine bounded_roots(coeff,degree,roots,nroots,ok)
  integer, intent(in) :: degree
  real*8, intent(in) :: coeff(0:6)
  real*8, intent(out) :: roots(16)
  integer, intent(out) :: nroots
  logical, intent(out) :: ok
  real*8 :: c(0:6),derivative(0:6),critical(16),bounds(18),scale,tol,a,b,fa,fb,m,fm
  integer :: n,nc,i,it
  logical :: child_ok
  roots=0.d0; nroots=0; ok=.false.; n=degree
  scale=maxval(abs(coeff(:n)))
  if(.not.ieee_is_finite(scale).or.scale<=tiny(scale)) return
  ! Normalize the complete fixed-size buffer in one assignment.
  c=coeff/scale
  do while(n>0)
    if(c(n)/=0.d0) exit
    n=n-1
  enddo
  tol=64.d0*epsilon(1.d0)*sum(abs(c(:n)))
  if(n==0) then
    ok=abs(c(0))>tol
    return
  endif
  derivative=0.d0
  do i=1,n
    derivative(i-1)=i*c(i)
  enddo
  call bounded_roots(derivative,n-1,critical,nc,child_ok)
  if(.not.child_ok) return
  ok=.true.
  bounds(1)=0.d0; bounds(2:nc+1)=critical(:nc); bounds(nc+2)=1.d0
  do i=1,nc+2
    a=bounds(i); fa=polynomial(c,n,a)
    if(abs(fa)<=tol) call insert_root(a)
    if(i==nc+2) exit
    b=bounds(i+1); fb=polynomial(c,n,b)
    if((fa<0.d0.and.fb>0.d0).or.(fa>0.d0.and.fb<0.d0)) then
      do it=1,64
        m=0.5d0*(a+b); fm=polynomial(c,n,m)
        if(fm==0.d0.or.b-a<=8.d0*epsilon(1.d0)) exit
        if((fa<0.d0).eqv.(fm<0.d0)) then
          a=m; fa=fm
        else
          b=m
        endif
      enddo
      call insert_root(m)
    endif
  enddo
contains
  subroutine insert_root(x)
    real*8, intent(in) :: x
    integer :: j
    if(any(abs(roots(:nroots)-x)<=32.d0*epsilon(1.d0))) return
    if(nroots==size(roots)) then
      ok=.false.
      return
    endif
    j=nroots
    do while(j>0)
      if(roots(j)<=x) exit
      roots(j+1)=roots(j); j=j-1
    enddo
    roots(j+1)=x; nroots=nroots+1
  end subroutine
end subroutine

! Refine only an inaccurate candidate, using the original FE/ballistic equations.
! The neighbourhood uses candidate separation, NOT the final polynomial bracket:
! cancellation can put the true root outside that last floating-point bracket.
! Both coordinates remain local to their seed; unsuccessful/ambiguous motion is fatal.
subroutine refine_side_intersection(nodes,elements,e,side,event,dt,u,lam,u_radius,l_radius, &
                                    geom_tol,iterations,residual,ok,reason)
  type(type_node_list), intent(in) :: nodes
  type(type_element_list), intent(in) :: elements
  integer, intent(in) :: e,side
  type(trajectory_origin), intent(in) :: event
  real*8, intent(in) :: dt,u_radius,l_radius,geom_tol
  real*8, intent(inout) :: u,lam
  integer, intent(out) :: iterations
  real*8, intent(out) :: residual
  logical, intent(out) :: ok
  character(len=*), intent(out) :: reason
  integer, parameter :: max_iterations=8,max_backtracks=12
  type(trajectory_point) :: point
  real*8 :: seed_u,seed_l,st(2),wall(2),tangent(2),normal(3),f(2),d(2),tn,dn,det
  real*8 :: a(2),b(2),du,dl,alpha,trial_u,trial_l,trial_residual,old_residual
  integer :: it,back
  logical :: valid,accepted
  seed_u=u;seed_l=lam;iterations=0;residual=huge(1.d0);ok=.false.
  reason='intersection refinement has invalid input or isolation neighbourhood'
  if(.not.all(ieee_is_finite([u,lam,dt,u_radius,l_radius,geom_tol]))) return
  if(min(u_radius,l_radius,dt,geom_tol)<=0.d0) return
  if(min(u,lam)<0.d0.or.max(u,lam)>1.d0) return
  do it=0,max_iterations
    reason='intersection refinement has invalid FE or trajectory geometry'
    call side_geometry(nodes,elements,e,side,u,st,wall,tangent,normal,valid)
    if(.not.valid) return
    point%x=event%x_old;point%v=event%v_old
    call push_trajectory(point,lam*dt)
    if(.not.all(ieee_is_finite([point%x,point%v])).or.point%x(1)<=0.d0) return
    f=wall-point%x(:2);residual=norm2(f)
    if(.not.ieee_is_finite(residual)) return
    if(residual<=geom_tol) then
      ok=.true.;reason='';return
    endif
    reason='intersection refinement iteration limit exceeded'
    if(it==max_iterations) return
    d=dt*point%v(:2);tn=norm2(tangent);dn=norm2(d)
    reason='intersection refinement singular or poorly conditioned Jacobian'
    if(min(tn,dn)<=tiny(1.d0)) return
    ! Normalize columns to avoid squaring dimensional scales in the 2x2 solve.
    a=tangent/tn;b=d/dn;det=a(2)*b(1)-a(1)*b(2)
    if(.not.ieee_is_finite(det)) return
    if(abs(det)*min(tn,dn)/max(tn,dn)<=sqrt(epsilon(1.d0))) return
    du=(f(1)*b(2)-f(2)*b(1))/(det*tn)
    dl=(a(2)*f(1)-a(1)*f(2))/(det*dn)
    reason='intersection refinement nonfinite correction'
    if(.not.all(ieee_is_finite([du,dl]))) return
    ! A full Newton proposal outside the seed neighbourhood has uncertain root
    ! identity. Do not backtrack it into apparent convergence on another branch.
    reason='intersection refinement cannot preserve isolated root identity'
    if(abs(u+du-seed_u)>=u_radius.or.abs(lam+dl-seed_l)>=l_radius) return
    alpha=1.d0;accepted=.false.;old_residual=residual
    do back=0,max_backtracks
      trial_u=u+alpha*du;trial_l=lam+alpha*dl
      if(min(trial_u,trial_l)>=-st_tolerance.and.max(trial_u,trial_l)<=1.d0+st_tolerance) then
        trial_u=max(0.d0,min(1.d0,trial_u));trial_l=max(0.d0,min(1.d0,trial_l))
        call side_geometry(nodes,elements,e,side,trial_u,st,wall,tangent,normal,valid)
        if(valid) then
          point%x=event%x_old;point%v=event%v_old
          call push_trajectory(point,trial_l*dt)
          if(all(ieee_is_finite([point%x,point%v])).and.point%x(1)>0.d0) then
            trial_residual=norm2(wall-point%x(:2))
            if(ieee_is_finite(trial_residual).and.trial_residual<old_residual) then
              accepted=.true.;exit
            endif
          endif
        endif
      endif
      alpha=0.5d0*alpha
    enddo
    reason='intersection refinement backtracking failed to decrease residual'
    if(.not.accepted) return
    u=trial_u;lam=trial_l;iterations=iterations+1
  enddo
end subroutine

!> Earliest resolved outward face crossing after entry into one element.
!> Roots use the ORIGINAL substep trajectory and a monotonically increasing fraction.
!> Eliminate lambda through Z(lambda); for vZ=0 solve wall Z then radial motion.
subroutine first_element_exit(nodes,elements,e,event,endpoint,dt,entry_fraction,entered,hit,error)
  use mod_parameters, only: n_order
  type(type_node_list), intent(in) :: nodes
  type(type_element_list), intent(in) :: elements
  integer, intent(in) :: e
  type(trajectory_origin), intent(in) :: event
  type(trajectory_point), intent(in) :: endpoint
  real*8, intent(in) :: dt,entry_fraction
  logical, intent(in) :: entered
  type(side_correction), intent(out) :: hit
  integer, intent(out) :: error
  type(trajectory_point) :: trajectory,before,after
  real*8 :: ends(2,2),derivs(2,2),st(2),normal(3),tangent(2),wall(2),c(2,0:3),cr(0:3),cz(0:3)
  real*8 :: poly(0:6),roots(16),lambda(2),d(3),scale,geom_tol,r0,aa,bb,cc,disc,q,u,vn,speed
  real*8 :: lam,best,tie_tol,delta,old_wall(2),best_u,tangent_rate,tangent_fraction
  real*8 :: u_radius,l_radius,residual
  integer :: iterations
  integer :: side,i,k,m,nroots,nlambda,choices,best_side
  logical :: sides(4),ok,ambiguous
  character(len=100) :: reason
  hit%original_st=endpoint%st
  error=1; reason='missing or mismatched pre-substep event record'
  if(event%element<1.or.event%element>elements%n_elements) goto 800
  if(e<1.or.e>elements%n_elements) goto 800
  reason='invalid element vertex connectivity'
  if(any(elements%element(e)%vertex<1).or.any(elements%element(e)%vertex>nodes%n_nodes)) goto 800
  if(any(elements%element(event%element)%vertex<1).or. &
     any(elements%element(event%element)%vertex>nodes%n_nodes)) goto 800
  reason='invalid pre-substep state or timestep'
  if(event%dt/=dt.or..not.all(ieee_is_finite([event%x_old,event%v_old,event%st_old]))) goto 800
  if(event%x_old(1)<=0.d0.or.any(event%st_old<0.d0).or.any(event%st_old>1.d0)) goto 800
  reason='trajectory crossing requires the current cubic geometry representation'
  if(n_order/=3) goto 800
  trajectory=endpoint; trajectory%x=event%x_old; trajectory%v=event%v_old
  call push_trajectory(trajectory,dt)
  scale=max(norm2(dt*event%v_old),tiny(1.d0))
  geom_tol=512.d0*epsilon(1.d0)*max(1.d0,norm2(event%x_old(:2)),norm2(endpoint%x(:2)))+1.d-10*scale
  speed=norm2(endpoint%v)
  if(.not.all(ieee_is_finite([endpoint%x,endpoint%v,speed]))) goto 800
  reason='pre-substep trajectory does not reproduce the stored ballistic endpoint'
  if(.not.all(ieee_is_finite([trajectory%x,trajectory%v]))) goto 800
  if(norm2(trajectory%x-endpoint%x)>16.d0*geom_tol.or. &
     norm2(trajectory%v-endpoint%v)>roundoff_tolerance*max(speed,tiny(speed))) goto 800
  call interp_RZ(nodes,elements,event%element, &
                 event%st_old(1),event%st_old(2),old_wall(1),old_wall(2))
  reason='pre-substep position inconsistent with its retained element/st'
  if(norm2(old_wall-event%x_old(:2))>max(16.d0*geom_tol,1.d-8*max(1.d0,norm2(old_wall)))) goto 800
  hit%max_residual=norm2(old_wall-event%x_old(:2))
  sides=boundary_sides(nodes,elements,e); best=huge(1.d0); choices=0; best_side=0; best_u=0.d0; ambiguous=.false.
  tangent_fraction=huge(1.d0)
  do side=1,4
    hit%examined=hit%examined+1
    reason='boundary-node flags disagree with nonzero neighbour on candidate side'
    if(sides(side).and.elements%element(e)%neighbours(side)/=0) goto 800
    call side_geometry(nodes,elements,e,side,0.d0,st,ends(:,1),derivs(:,1),normal,ok)
    if(.not.ok) then
      reason='invalid candidate endpoint geometry'; goto 800
    endif
    call side_geometry(nodes,elements,e,side,1.d0,st,ends(:,2),derivs(:,2),normal,ok)
    if(.not.ok) then
      reason='invalid candidate endpoint geometry'; goto 800
    endif
    c(:,0)=ends(:,1); c(:,1)=derivs(:,1)
    c(:,2)=3.d0*(ends(:,2)-ends(:,1))-2.d0*derivs(:,1)-derivs(:,2)
    c(:,3)=2.d0*(ends(:,1)-ends(:,2))+derivs(:,1)+derivs(:,2)
    scale=max(norm2(dt*event%v_old),norm2(ends(:,2)-ends(:,1)),norm2(derivs(:,1)),norm2(derivs(:,2)))
    geom_tol=512.d0*epsilon(1.d0)*max(1.d0,maxval(abs(ends)))+1.d-10*scale
    cr=c(1,:)/scale; cz=c(2,:)/scale
    cr(0)=(c(1,0)-event%x_old(1))/scale; cz(0)=(c(2,0)-event%x_old(2))/scale
    d=dt*event%v_old/scale; r0=event%x_old(1)/scale
    poly=0.d0
    if(d(2)/=0.d0) then
      poly(:3)=2.d0*r0*(d(2)**2*cr-d(1)*d(2)*cz)
      do i=0,3
        do k=0,3
          poly(i+k)=poly(i+k)+d(2)**2*cr(i)*cr(k)-(d(1)**2+d(3)**2)*cz(i)*cz(k)
        enddo
      enddo
      call bounded_roots(poly,6,roots,nroots,ok)
    else
      poly(:3)=cz
      call bounded_roots(poly,3,roots,nroots,ok)
    endif
    reason='coincident/grazing or unresolved polynomial intersection'
    if(.not.ok) goto 800
    choices=0
    do i=1,nroots
      u=roots(i); nlambda=1
      if(d(2)/=0.d0) then
        lambda(1)=polynomial(cz,3,u)/d(2)
      else
        aa=d(1)**2+d(3)**2; bb=2.d0*r0*d(1)
        cc=-polynomial(cr,3,u)*(2.d0*r0+polynomial(cr,3,u))
        if(aa<=tiny(aa)) cycle
        disc=bb*bb-4.d0*aa*cc
        if(disc<0.d0) cycle
        q=-0.5d0*(bb+sign(sqrt(disc),bb))
        if(q==0.d0) then
          lambda(1)=-bb/(2.d0*aa)
        else
          lambda(1)=q/aa; lambda(2)=cc/q; nlambda=2
        endif
      endif
      do m=1,nlambda
        u=roots(i);lam=lambda(m)
        if(.not.ieee_is_finite(lam).or.lam<-st_tolerance.or.lam>1.d0+st_tolerance) then
          hit%outside=hit%outside+1
          cycle
        endif
        lam=max(0.d0,min(1.d0,lam))
        call side_geometry(nodes,elements,e,side,u,st,wall,tangent,normal,ok)
        reason='invalid geometry at candidate intersection'
        if(.not.ok) goto 800
        trajectory=endpoint; trajectory%x=event%x_old; trajectory%v=event%v_old
        call push_trajectory(trajectory,lam*dt)
        reason='candidate intersection residual exceeds geometry tolerance'
        hit%max_residual=max(hit%max_residual,norm2(trajectory%x(:2)-wall))
        if(norm2(trajectory%x(:2)-wall)>geom_tol) then
          u_radius=0.25d0;l_radius=0.25d0
          if(i>1) u_radius=min(u_radius,0.25d0*(roots(i)-roots(i-1)))
          if(i<nroots) u_radius=min(u_radius,0.25d0*(roots(i+1)-roots(i)))
          if(nlambda==2) l_radius=min(l_radius,0.25d0*abs(lambda(2)-lambda(1)))
          hit%refinement_counts(1)=hit%refinement_counts(1)+1
          hit%refinement_residuals(1)=max(hit%refinement_residuals(1),norm2(trajectory%x(:2)-wall))
          call refine_side_intersection(nodes,elements,e,side,event,dt,u,lam,u_radius,l_radius, &
                                        geom_tol,iterations,residual,ok,reason)
          hit%refinement_counts(3)=hit%refinement_counts(3)+iterations
          hit%refinement_counts(4)=max(hit%refinement_counts(4),iterations)
          hit%refinement_residuals(2)=max(hit%refinement_residuals(2),residual)
          if(.not.ok) then
            hit%refinement_counts(2)=hit%refinement_counts(2)+1
            goto 800
          endif
          call side_geometry(nodes,elements,e,side,u,st,wall,tangent,normal,ok)
          reason='invalid geometry after intersection refinement'
          if(.not.ok) goto 800
          trajectory=endpoint;trajectory%x=event%x_old;trajectory%v=event%v_old
          call push_trajectory(trajectory,lam*dt)
        endif
        reason='invalid refined candidate domain or residual'
        if(.not.all(ieee_is_finite([u,lam,trajectory%x,trajectory%v]))) goto 800
        if(trajectory%x(1)<=0.d0.or.min(u,lam)<0.d0.or.max(u,lam)>1.d0) goto 800
        if(norm2(trajectory%x(:2)-wall)>geom_tol) goto 800
        if(lam<entry_fraction-st_tolerance) cycle
        vn=dot_product(trajectory%v,normal)
        ! Entry and tangency are not outward exits. Do not relax normal tolerance.
        if(abs(vn)<=roundoff_tolerance*speed) tangent_fraction=min(tangent_fraction,lam)
        if(vn>=-roundoff_tolerance*speed) then
          hit%nonexit=hit%nonexit+1
          cycle
        endif
        tie_tol=max(st_tolerance,8.d0*geom_tol/max(-vn*dt,tiny(1.d0)))
        reason='outward crossing does not progress beyond the element entry'
        if(entered.and.lam<=entry_fraction+tie_tol)then
          call uncertain_exit();cycle
        endif
        if(min(u,1.d0-u)<=st_tolerance) then
          call uncertain_exit();cycle
        endif
        tangent_rate=abs(dt*dot_product(trajectory%v(:2),tangent))/sum(tangent**2)
        delta=min(1.d-4,0.1d0*min(u,1.d0-u)/max(tangent_rate,1.d0))
        if(lam>entry_fraction) delta=min(delta,0.5d0*(lam-entry_fraction))
        if(lam<1.d0) delta=min(delta,0.5d0*(1.d0-lam))
        reason='outward crossing cannot be separated from roundoff/grazing'
        if(-vn*dt*delta<=8.d0*geom_tol)then
          call uncertain_exit();cycle
        endif
        before=endpoint; before%x=event%x_old; before%v=event%v_old
        after=before
        call push_trajectory(before,(lam-delta)*dt)
        call push_trajectory(after,(lam+delta)*dt)
        if(dot_product(before%x(:2)-wall,normal(:2))<=0.d0.or. &
           dot_product(after%x(:2)-wall,normal(:2))>=0.d0)then
          call uncertain_exit();cycle
        endif
        choices=choices+1
        tie_tol=max(st_tolerance,8.d0*geom_tol/max(-vn*dt,tiny(1.d0)))
        if(abs(lam-best)<=tie_tol) then
          ! The same root may arise twice in a degenerate radial quadratic.
          if(side/=best_side.or.abs(u-best_u)>st_tolerance)then
            ambiguous=.true.;hit%coalescible=.true.
          endif
        else if(lam<best) then
          best=lam; best_side=side; best_u=u; ambiguous=.false.;hit%coalescible=.false.
          hit%side=side; hit%st=st; hit%wall=wall; hit%normal=normal
          hit%lambda=lam; hit%vn=vn
        endif
      enddo
    enddo
    if(choices==0) hit%no_exit=hit%no_exit+1
  enddo
  reason='unresolved tangency before any verified outward exit'
  if(tangent_fraction<huge(1.d0)) then
    if(tangent_fraction<=best+st_tolerance) goto 800
  endif
  reason='no verified outward face intersection during the remaining substep'
  if(best==huge(1.d0)) then
    error=2 ! All faces resolved, but there is no outward exit in this interval.
    return
  endif
  reason='earliest outward crossing is a corner or has an unresolved tie'
  if(ambiguous.or.best_side==0) goto 800
  hit%element=e; error=0
  return
800 continue
  hit%reason=reason
contains
  subroutine uncertain_exit()
    ! Keep chronological ordering even when one internal corner cannot be
    ! resolved at the much tighter ordinary intersection tolerance.
    if(lam<best-tie_tol)then
      best=lam;best_side=0;best_u=u;ambiguous=.true.;hit%coalescible=.true.
      hit%element=e;hit%side=side;hit%st=st;hit%wall=wall;hit%normal=normal
      hit%lambda=lam;hit%vn=vn
    else if(abs(lam-best)<=tie_tol)then
      ambiguous=.true.;hit%coalescible=.true.
    endif
  end subroutine
end subroutine

! Parameter-oriented vertices (sides 3/4 run opposite to polygon ordering).
function face_vertices(nodes,elements,e,side) result(vertices)
  type(type_node_list), intent(in) :: nodes
  type(type_element_list), intent(in) :: elements
  integer, intent(in) :: e,side
  integer :: vertices(2),v(4)
  v=elements%element(e)%vertex
  select case(side)
  case(1);vertices=v([1,2])
  case(2);vertices=v([2,3])
  case(3);vertices=v([4,3])
  case(4);vertices=v([1,4])
  end select
end function

! Duplicate node IDs occur at COMPASS mesh seams. A cubic face is uniquely
! determined by its two endpoints and endpoint tangents. Require EXACT agreement
! of those physical data; proximity at the crossing alone never repairs topology.
logical function identical_cubic_faces(nodes,elements,e,side,next,entry_side,reversed) result(same)
  type(type_node_list), intent(in) :: nodes
  type(type_element_list), intent(in) :: elements
  integer, intent(in) :: e,side,next,entry_side
  logical, intent(in) :: reversed
  real*8 :: u,v,orientation,st(2),a(2),b(2),da(2),db(2),normal(3)
  integer :: k
  logical :: ok
  same=.false.;orientation=1.d0;if(reversed) orientation=-1.d0
  if(any(elements%element(next)%vertex<1).or.any(elements%element(next)%vertex>nodes%n_nodes)) return
  do k=0,1
    u=real(k,8);v=u;if(reversed) v=1.d0-u
    call side_geometry(nodes,elements,e,side,u,st,a,da,normal,ok)
    if(.not.ok) return
    call side_geometry(nodes,elements,next,entry_side,v,st,b,db,normal,ok)
    if(.not.ok) return
    if(any(a/=b).or.any(da/=orientation*db)) return
  enddo
  same=.true.
end function

!> Bounded, local face traversal; never uses the locator's Newton path as a trajectory.
subroutine trace_trajectory(nodes,elements,event,endpoint,dt,hit,error)
  type(type_node_list), intent(in) :: nodes
  type(type_element_list), intent(in) :: elements
  type(trajectory_origin), intent(in) :: event
  type(trajectory_point), intent(in) :: endpoint
  real*8, intent(in) :: dt
  type(side_correction), intent(out) :: hit
  integer, intent(out) :: error
  type(side_correction) :: face
  integer :: refinement_counts(4)
  real*8 :: refinement_residuals(2),trace_u(max_trace_elements)
  type(trajectory_point) :: impact
  integer :: e,next,side,entry_side,h,n1(2),n2(2),visited(max_trace_elements),faces(max_trace_elements)
  real*8 :: fractions(max_trace_elements),cursor,u,st(2),wall(2),tangent(2),normal(3),tol,residual,max_residual
  integer :: clusters,cluster_next
  real*8 :: probe_fraction,cluster_tol,cluster_extent,max_cluster_extent,min_cluster_tol
  logical :: physical(4),ok,cluster_ok
  character(len=100) :: reason
  error=1;e=event%element;cursor=0.d0
  clusters=0;max_cluster_extent=0.d0;min_cluster_tol=huge(1.d0)
  visited=0;faces=0;fractions=0.d0;max_residual=0.d0
  refinement_counts=0;refinement_residuals=0.d0;trace_u=0.d0
  do h=1,max_trace_elements
    call first_element_exit(nodes,elements,e,event,endpoint,dt,cursor,h>1,face,error)
    refinement_counts(:3)=refinement_counts(:3)+face%refinement_counts(:3)
    refinement_counts(4)=max(refinement_counts(4),face%refinement_counts(4))
    refinement_residuals=max(refinement_residuals,face%refinement_residuals)
    hit%refinement_counts=refinement_counts;hit%refinement_residuals=refinement_residuals
    hit%trace_elements=visited;hit%trace_faces=faces;hit%trace_fractions=fractions;hit%trace_u=trace_u
    hit%trace_elements(h)=e;hit%n_trace=h-1
    max_residual=max(max_residual,face%max_residual)
    hit%n_visits=h; hit%max_residual=max_residual
    if(error==2) then
      call invert_endpoint(nodes,elements,e,endpoint%x(:2),dt*norm2(event%v_old),st,residual,ok)
      if(.not.ok) then
        hit%reason='no outward exit, but endpoint inversion is invalid or ambiguous'
        error=1
        return
      endif
      hit%reason='verified interior endpoint; no physical wall crossing'
      hit%status=INTERIOR_ENDPOINT;hit%element=e;hit%st=st;hit%wall=endpoint%x(:2)
      hit%n_trace=h-1;hit%lambda=1.d0;hit%max_residual=max(max_residual,residual)
      hit%trace_elements=visited;hit%trace_faces=faces;hit%trace_fractions=fractions
      error=0
      return
    else if(error/=0) then
      if(face%coalescible.and.trim(face%reason)=='earliest outward crossing is a corner or has an unresolved tie')then
        call resolve_internal_cluster(nodes,elements,e,event,face,cursor,cluster_next,probe_fraction, &
                                      cluster_tol,cluster_extent,cluster_ok)
        if(cluster_ok)then
          clusters=clusters+1;max_cluster_extent=max(max_cluster_extent,cluster_extent)
          min_cluster_tol=min(min_cluster_tol,cluster_tol)
          visited(h)=e;faces(h)=-face%side;fractions(h)=probe_fraction
          trace_u(h)=face%st(1);if(face%side==2.or.face%side==4)trace_u(h)=face%st(2)
          cursor=probe_fraction;e=cluster_next
          hit%local_clusters=clusters;hit%cluster_extent=max_cluster_extent;hit%cluster_tolerance=min_cluster_tol
          cycle
        endif
      endif
      hit%reason=face%reason
      return
    endif
    error=1;side=face%side;visited(h)=e;faces(h)=side;fractions(h)=face%lambda
    trace_u(h)=face%st(1);if(side==2.or.side==4) trace_u(h)=face%st(2)
    physical=boundary_sides(nodes,elements,e)
    next=elements%element(e)%neighbours(side)
    if(physical(side)) then
      ! first_element_exit has already required neighbour=0 for every physical side.
      hit=face;hit%corrected=.true.;hit%n_trace=h
      hit%local_clusters=clusters;hit%cluster_extent=max_cluster_extent
      if(clusters>0)hit%cluster_tolerance=min_cluster_tol
      hit%refinement_counts=refinement_counts;hit%refinement_residuals=refinement_residuals;hit%trace_u=trace_u
      hit%status=WALL_CROSSING;hit%n_visits=h;hit%max_residual=max_residual
      hit%trace_elements(:h)=visited(:h);hit%trace_faces(:h)=faces(:h);hit%trace_fractions(:h)=fractions(:h)
      error=0
      return
    endif
    reason='interior face has no valid neighbouring element'
    if(next<1.or.next>elements%n_elements) goto 800
    reason='interior face points back to the same element'
    if(next==e) goto 800
    ! A radial turn may legitimately revisit a cell. first_element_exit requires
    ! strictly later exits; the visit bound also terminates cyclic connectivity.
    reason='neighbour relation is not uniquely reciprocal'
    if(count(elements%element(next)%neighbours==e)/=1) goto 800
    entry_side=findloc(elements%element(next)%neighbours,e,dim=1)
    n1=face_vertices(nodes,elements,e,side);n2=face_vertices(nodes,elements,next,entry_side)
    u=face%st(1);if(side==2.or.side==4) u=face%st(2)
    reason='neighbour face has neither shared vertices nor identical cubic geometry'
    if(all(n1==n2(2:1:-1))) then
      u=1.d0-u
    else if(.not.all(n1==n2)) then
      if(identical_cubic_faces(nodes,elements,e,side,next,entry_side,.false.)) then
        ! Same physical parameterization despite duplicate node IDs.
      else if(identical_cubic_faces(nodes,elements,e,side,next,entry_side,.true.)) then
        u=1.d0-u
      else
        goto 800
      endif
    endif
    call side_geometry(nodes,elements,next,entry_side,u,st,wall,tangent,normal,ok)
    reason='invalid neighbour-entry geometry'
    if(.not.ok) goto 800
    tol=512.d0*epsilon(1.d0)*max(1.d0,norm2(wall))+ &
        1.d-10*max(norm2(dt*event%v_old),norm2(tangent))
    reason='neighbour geometry is discontinuous or normals are not opposite'
    if(norm2(wall-face%wall)>tol.or.norm2(normal+face%normal)>1.d-8) goto 800
    impact=endpoint;impact%x=event%x_old;impact%v=event%v_old
    call push_trajectory(impact,face%lambda*dt)
    reason='trajectory does not enter the neighbour through the shared face'
    if(dot_product(impact%v,normal)<=roundoff_tolerance*norm2(impact%v)) goto 800
    cursor=face%lambda;e=next
  enddo
  write(reason,'(A,I0,A)') 'element traversal limit exceeded (',max_trace_elements,' elements)'
800 continue
  hit%reason=reason;hit%element=e;hit%n_visits=min(h,max_trace_elements)
  hit%trace_elements=visited;hit%trace_faces=faces;hit%trace_fractions=fractions
end subroutine

! A tolerance is an error budget for intersections, never for particle location.
! The smallest local mapping scale caps the budget in thin/curved cells.
real*8 function local_intersection_tolerance(nodes,elements,e) result(tol)
  type(type_node_list),intent(in)::nodes
  type(type_element_list),intent(in)::elements
  integer,intent(in)::e
  real*8::s,t,r,z,rs,rt,zs,zt,jac,length
  integer::i
  tol=0.d0
  if(e<1.or.e>elements%n_elements)return
  if(any(elements%element(e)%vertex<1).or.any(elements%element(e)%vertex>nodes%n_nodes))return
  tol=wall_intersection_tolerance
  do i=1,5
    s=.5d0;t=.5d0
    if(i<=4)then
      s=0.d0;t=0.d0
      if(i==2.or.i==3)s=1.d0
      if(i>=3)t=1.d0
    endif
    call interp_RZ(nodes,elements,e,s,t,r,rs,rt,z,zs,zt)
    jac=rs*zt-rt*zs;length=norm2([rs,zs])+norm2([rt,zt])
    if(.not.all(ieee_is_finite([r,z,rs,rt,zs,zt,jac,length])).or.jac<=0.d0.or.length<=0.d0)then
      tol=0.d0;return
    endif
    ! det(J)/(||J_s||+||J_t||) is a lower bound on the smallest
    ! singular value of the local mapping, in metres per reference-cell unit.
    tol=min(tol,.01d0*jac/length)
  enddo
end function

! Resolve an indistinguishable INTERNAL vertex cluster in one local operation.
! The closed vertex star and a wall-free ball certify that advancing the SEARCH
! cursor cannot skip a physical contact. The physical trajectory is never moved.
subroutine resolve_internal_cluster(nodes,elements,e,event,face,cursor,next,probe_fraction,tol,extent,ok)
  type(type_node_list),intent(in)::nodes
  type(type_element_list),intent(in)::elements
  integer,intent(in)::e
  type(trajectory_origin),intent(in)::event
  type(side_correction),intent(in)::face
  real*8,intent(in)::cursor
  integer,intent(out)::next
  real*8,intent(out)::probe_fraction,tol,extent
  logical,intent(out)::ok
  integer,parameter::limit=16
  integer::patch(limit),corner(limit),count_cells,head,k,j,iv,side,nb,entry,v(4),found,ninside
  integer::incident(2,4)
  real*8::vertex(2),position(2,4),distance(4),radius,ctol,speed,st(2),wall(2),tangent(2),normal(3)
  real*8::ends(2,2),derivs(2,2),bez(2,4),gap(2),residual,chosen_st(2)
  logical::valid,physical(4),same
  type(trajectory_point)::probe,contact
  ok=.false.;next=0;probe_fraction=cursor;extent=0.d0
  tol=local_intersection_tolerance(nodes,elements,e)
  if(tol<=0.d0.or..not.face%coalescible)return
  v=elements%element(e)%vertex
  do j=1,4
    position(:,j)=nodes%node(v(j))%x(1,1,:)*elements%element(e)%size(j,1)
  enddo
  distance=[(norm2(position(:,j)-face%wall),j=1,4)]
  iv=minloc(distance,dim=1);vertex=position(:,iv)
  if(distance(iv)>tol/8.d0.or.nodes%node(v(iv))%boundary>0)return
  incident=reshape([1,4,1,2,2,3,3,4],[2,4])
  patch=0;corner=0;patch(1)=e;corner(1)=iv;count_cells=1;head=1
  do while(head<=count_cells)
    k=patch(head);v=elements%element(k)%vertex;iv=corner(head)
    ctol=local_intersection_tolerance(nodes,elements,k)
    if(ctol<=0.d0)return
    tol=min(tol,ctol)
    if(nodes%node(v(iv))%boundary>0)return
    do j=1,2
      side=incident(j,iv);nb=elements%element(k)%neighbours(side)
      if(nb<1.or.nb>elements%n_elements.or.nb==k)return
      if(count(elements%element(nb)%neighbours==k)/=1)return
      entry=findloc(elements%element(nb)%neighbours,k,dim=1)
      same=identical_cubic_faces(nodes,elements,k,side,nb,entry,.false.)
      if(.not.same)same=identical_cubic_faces(nodes,elements,k,side,nb,entry,.true.)
      if(.not.same)return ! A physical tolerance does not repair broken topology.
      if(any(patch(:count_cells)==nb))cycle
      if(count_cells==limit)return
      if(any(elements%element(nb)%vertex<1).or.any(elements%element(nb)%vertex>nodes%n_nodes))return
      do found=1,4
        position(:,found)=nodes%node(elements%element(nb)%vertex(found))%x(1,1,:)* &
                            elements%element(nb)%size(found,1)
      enddo
      distance=[(norm2(position(:,found)-vertex),found=1,4)]
      found=minloc(distance,dim=1)
      if(distance(found)>128.d0*epsilon(1.d0)*max(1.d0,norm2(vertex)))return
      count_cells=count_cells+1;patch(count_cells)=nb;corner(count_cells)=found
    enddo
    head=head+1
  enddo
  if(count_cells<3)return
  ! The full cubic Bernstein hull of every outer patch face must exclude the
  ! ball. This conservative test cannot mistake an unexamined wall for interior.
  radius=tol
  if(norm2(face%wall-vertex)>radius/8.d0)return
  do head=1,count_cells
    k=patch(head);iv=corner(head)
    do side=1,4
      if(any(incident(:,iv)==side))cycle
      call side_geometry(nodes,elements,k,side,0.d0,st,ends(:,1),derivs(:,1),normal,valid)
      if(.not.valid)return
      call side_geometry(nodes,elements,k,side,1.d0,st,ends(:,2),derivs(:,2),normal,valid)
      if(.not.valid)return
      bez(:,1)=ends(:,1);bez(:,2)=ends(:,1)+derivs(:,1)/3.d0
      bez(:,3)=ends(:,2)-derivs(:,2)/3.d0;bez(:,4)=ends(:,2)
      gap=max(0.d0,max(minval(bez,dim=2)-vertex,vertex-maxval(bez,dim=2)))
      if(norm2(gap)<=radius+128.d0*epsilon(1.d0)*max(1.d0,maxval(abs(bez))))return
    enddo
  enddo
  speed=norm2(event%v_old)
  if(speed<=0.d0)return
  probe_fraction=min(1.d0,max(cursor,face%lambda)+.5d0*tol/(event%dt*speed))
  if(probe_fraction<=cursor.or.probe_fraction<=face%lambda)return
  contact%x=event%x_old;contact%v=event%v_old
  call push_trajectory(contact,face%lambda*event%dt)
  extent=norm2(contact%x(:2)-vertex)+(probe_fraction-face%lambda)*event%dt*speed
  if(extent>=radius)return
  probe%x=event%x_old;probe%v=event%v_old
  call push_trajectory(probe,probe_fraction*event%dt)
  ninside=0
  ! Enumeration is bounded and the accepted interior cell must be unique;
  ! this is deterministic and cannot depend on a locator search order.
  do head=1,count_cells
    k=patch(head)
    call invert_endpoint(nodes,elements,k,probe%x(:2),event%dt*speed,st,residual,valid)
    if(.not.valid)cycle
    if(minval(st)<=st_tolerance.or.maxval(st)>=1.d0-st_tolerance)cycle
    ninside=ninside+1;next=k;chosen_st=st
  enddo
  if(ninside/=1)return
  ok=.true.
end subroutine

! Exceptional path only: the complete stored cubic boundary and exact ballistic
! trajectory determine the first physical exit independently of cell traversal.
! Higher precision resolves close internal-vertex crossings; no double-precision
! acceptance tolerance is enlarged. A nonmatching mesh corridor is never bypassed.
subroutine boundary_fallback(nodes,elements,event,endpoint,dt,hit)
  type(type_node_list), intent(in) :: nodes
  type(type_element_list), intent(in) :: elements
  type(trajectory_origin), intent(in) :: event
  type(trajectory_point), intent(in) :: endpoint
  real*8, intent(in) :: dt
  type(side_correction), intent(inout) :: hit
  integer, parameter :: qp=selected_real_kind(30,300)
  real(qp) :: c(2,0:3),other(2,0:3),poly(0:6),cr(0:3),cz(0:3),rr(16),ll(32),uu(32)
  real(qp) :: x0(2),d(3),r0,aa,bb,cc,disc,q,lam,u,turn,blo(2),bhi(2),bez(2,4)
  real(qp) :: best,bestu,sep,delta,tan(2),normal(2),point(2),vel(2),vn,speed,rate,sgn0,sgn1,tol
  real*8 :: position_tol,normal_speed
  real*8 :: st(2),wall(2),tangent(2),norm(3),othernorm(3),res,travel,controls(2,16),base(2),ds(2),dtan(2),mix(2)
  integer :: e,side,next,entry,nr,nc,i,j,k,beste,bests,n1(2),n2(2),ntouch
  logical :: sides(4),ok,matched,tie,positive_before,positive_after
  type(trajectory_point) :: impact
  hit%fallback_attempted=.true.;hit%fallback_faces=0;hit%fallback_tangencies=0
  hit%status=UNRESOLVED;hit%reason='boundary fallback: invalid or nonmatching mesh corridor'
  x0=real(event%x_old(:2),qp);d=real(dt,qp)*real(event%v_old,qp);r0=x0(1)
  aa=d(1)**2+d(3)**2;turn=0
  if(aa>0)turn=max(0._qp,min(1._qp,-r0*d(1)/aa))
  blo=[min(r0,sqrt((r0+d(1))**2+d(3)**2),sqrt((r0+turn*d(1))**2+(turn*d(3))**2)), &
       min(x0(2),x0(2)+d(2))]
  bhi=[max(r0,sqrt((r0+d(1))**2+d(3)**2)),max(x0(2),x0(2)+d(2))]
  ! Outward rounding of the broad-phase box only. It never accepts a root.
  tol=512._qp*real(epsilon(1.d0),qp)*max(1._qp,maxval(abs([blo,bhi])))
  blo=blo-tol;bhi=bhi+tol;best=huge(1._qp);beste=0;bests=0;bestu=0;tie=.false.;ntouch=0
  do e=1,elements%n_elements
    if(any(elements%element(e)%vertex<1).or.any(elements%element(e)%vertex>nodes%n_nodes))return
    sides=boundary_sides(nodes,elements,e)
    do side=1,4
      call cubic(e,side,c)
      bez(:,1)=c(:,0);bez(:,2)=c(:,0)+c(:,1)/3
      bez(:,4)=sum(c,dim=2);bez(:,3)=bez(:,4)-(c(:,1)+2*c(:,2)+3*c(:,3))/3
      if(any(maxval(bez,dim=2)<blo).or.any(minval(bez,dim=2)>bhi))cycle
      next=elements%element(e)%neighbours(side)
      call side_geometry(nodes,elements,e,side,0.d0,st,wall,tangent,norm,ok)
      if(.not.ok)return
      call side_geometry(nodes,elements,e,side,1.d0,st,wall,tangent,norm,ok)
      if(.not.ok)return
      if(.not.sides(side))then
        if(next<1.or.next>elements%n_elements.or.next==e)return
        if(count(elements%element(next)%neighbours==e)/=1)return
        if(any(elements%element(next)%vertex<1).or.any(elements%element(next)%vertex>nodes%n_nodes))return
        entry=findloc(elements%element(next)%neighbours,e,dim=1)
        call cubic(next,entry,other)
        ! Compare complete physical cubic data, including shared-ID faces.
        matched=all(c==other)
        if(.not.matched)then
          cr=[sum(other(1,:)),-other(1,1)-2*other(1,2)-3*other(1,3),other(1,2)+3*other(1,3),-other(1,3)]
          cz=[sum(other(2,:)),-other(2,1)-2*other(2,2)-3*other(2,3),other(2,2)+3*other(2,3),-other(2,3)]
          ! Algebraic reversal may round in extended precision: its bound is
          ! 64 extended ulps of coefficients, never a geometric double tolerance.
          matched=maxval(abs(c(1,:)-cr))<=64*epsilon(1._qp)*max(1._qp,maxval(abs(c))) .and. &
                  maxval(abs(c(2,:)-cz))<=64*epsilon(1._qp)*max(1._qp,maxval(abs(c)))
        endif
        if(.not.matched)return
        call side_geometry(nodes,elements,e,side,.5d0,st,wall,tangent,norm,ok)
        if(.not.ok)return
        call side_geometry(nodes,elements,next,entry,.5d0,st,wall,tangent,othernorm,ok)
        if(.not.ok.or.norm2(norm+othernorm)>1.d-8)return
        cycle
      endif
      if(next/=0)return
      hit%fallback_faces=hit%fallback_faces+1
      cr=c(1,:);cr(0)=cr(0)-r0;cz=c(2,:);cz(0)=cz(0)-x0(2);poly=0
      if(d(2)/=0)then
        poly(:3)=2*r0*(d(2)**2*cr-d(1)*d(2)*cz)
        do i=0,3
          do j=0,3
            poly(i+j)=poly(i+j)+d(2)**2*cr(i)*cr(j)-(d(1)**2+d(3)**2)*cz(i)*cz(j)
          enddo
        enddo
        call roots_q(poly,6,rr,nr,ok)
      else
        poly(:3)=cz
        call roots_q(poly,3,rr,nr,ok)
      endif
      hit%reason='boundary fallback: coincident or unresolved boundary polynomial'
      if(.not.ok)return
      nc=0
      do i=1,nr
        u=rr(i)
        if(d(2)/=0)then
          lam=eval_q(cz,3,u)/d(2);call add_candidate(u,lam)
        else
          bb=2*r0*d(1);cc=-eval_q(cr,3,u)*(2*r0+eval_q(cr,3,u))
          if(aa==0)cycle
          disc=bb*bb-4*aa*cc
          if(disc<0)cycle
          q=-(bb+sign(sqrt(disc),bb))/2
          if(q==0)then
            call add_candidate(u,-bb/(2*aa))
          else
            call add_candidate(u,q/aa);call add_candidate(u,cc/q)
          endif
        endif
      enddo
      do k=1,nc
        u=uu(k);lam=ll(k)
        call path(lam,point,vel)
        tol=4096*epsilon(1._qp)*max(1._qp,maxval(abs(c)),maxval(abs(point)))
        hit%reason='boundary fallback: high-precision candidate has a geometric residual'
        if(norm2(point-[eval_q(c(1,:),3,u),eval_q(c(2,:),3,u)])>tol)return
        tan=c(:,1)+2*u*c(:,2)+3*u*u*c(:,3)
        if(norm2(tan)==0)return
        normal=[-tan(2),tan(1)];if(side==3.or.side==4)normal=-normal
        normal=normal/norm2(normal);vn=dot_product(vel,normal);speed=norm2(vel)
        hit%reason='boundary fallback: boundary vertex contact has no unique face normal'
        if(min(u,1-u)<=real(st_tolerance,qp))return
        sep=1
        do j=1,nc
          if(j/=k.and.abs(ll(j)-lam)>256*epsilon(1._qp))sep=min(sep,abs(ll(j)-lam))
        enddo
        rate=abs(dot_product(d(:2),tan))/sum(tan**2)
        delta=min(1.e-4_qp,sep/8,.05_qp*min(u,1-u)/max(rate,1._qp))
        ! Positive time margins are used if available; start/end contacts use
        ! the analytic trajectory extension solely to classify the local root.
        if(lam>256*epsilon(1._qp))delta=min(delta,lam/4)
        if(1-lam>256*epsilon(1._qp))delta=min(delta,(1-lam)/4)
        hit%reason='boundary fallback: unresolved sidedness or tangent degeneracy'
        call signed_curve(lam-delta,u,c,side,sgn0,ok)
        if(.not.ok)return
        call signed_curve(lam+delta,u,c,side,sgn1,ok)
        if(.not.ok)return
        ! Extended arithmetic resolves the stored polynomial, but a physical
        ! classification must also survive uncertainty in its binary64 input.
        ! This is a separate sidedness floor, not a relaxed root tolerance.
        tol=max(tol,128*real(epsilon(1.d0),qp)*max(1._qp,maxval(abs(c)),maxval(abs(point))))
        if(min(abs(sgn0),abs(sgn1))<=tol)return
        positive_before=sgn0>0;positive_after=sgn1>0
        if(positive_before.eqv.positive_after)then
          ntouch=ntouch+1;cycle
        endif
        if(.not.positive_before)cycle ! genuine inward entry, not a loss
        if(vn>=0)return ! inconsistent orientation remains an error
        if(lam<best-4096*epsilon(1._qp))then
          best=lam;bestu=u;bests=side;beste=e;tie=.false.
        else if(abs(lam-best)<=4096*epsilon(1._qp))then
          if(e/=beste.or.side/=bests)tie=.true.
        endif
      enddo
      hit%reason='boundary fallback: invalid or nonmatching mesh corridor'
    enddo
  enddo
  hit%fallback_tangencies=ntouch
  hit%reason='boundary fallback: earliest physical contact has an unresolved tie'
  if(tie)return
  if(beste>0)then
    call side_geometry(nodes,elements,beste,bests,real(bestu,8),st,wall,tangent,norm,ok)
    if(.not.ok)return
    impact%x=event%x_old;impact%v=event%v_old;call push_trajectory(impact,real(best,8)*dt)
    res=norm2(impact%x(:2)-wall)
    position_tol=local_intersection_tolerance(nodes,elements,beste)
    normal_speed=abs(dot_product(impact%v,norm))
    if(position_tol<=0.d0.or.normal_speed<=roundoff_tolerance*norm2(impact%v))return
    if(res*norm2(impact%v)/normal_speed>position_tol/4.d0)return
    hit%status=WALL_CROSSING;hit%element=beste;hit%side=bests;hit%lambda=real(best,8)
    hit%st=st;hit%wall=wall;hit%normal=norm;hit%vn=dot_product(impact%v,norm)
    hit%max_residual=max(hit%max_residual,res);hit%corrected=.true.;hit%fallback_used=.true.
    hit%reason='verified first physical exit by exceptional cubic-boundary search'
  else
    ! Boundary-root enumeration certified no exit. A bounded FE inverse must
    ! still certify a unique physical endpoint; no negative locator is trusted.
    travel=dt*norm2(event%v_old)
    do e=1,elements%n_elements
      ! Convex hull of the 16 tensor-product Bezier control points. The
      ! stored Hermite derivative coefficients are offsets from each vertex.
      do i=1,4
        j=elements%element(e)%vertex(i)
        base=nodes%node(j)%x(1,1,:)*elements%element(e)%size(i,1)
        ds=nodes%node(j)%x(1,2,:)*elements%element(e)%size(i,2)
        dtan=nodes%node(j)%x(1,3,:)*elements%element(e)%size(i,3)
        mix=nodes%node(j)%x(1,4,:)*elements%element(e)%size(i,4)
        controls(:,4*i-3)=base;controls(:,4*i-2)=base+ds
        controls(:,4*i-1)=base+dtan;controls(:,4*i)=base+ds+dtan+mix
      enddo
      res=512*epsilon(1.d0)*max(1.d0,maxval(abs(controls)))
      if(any(endpoint%x(:2)<minval(controls,dim=2)-res).or. &
         any(endpoint%x(:2)>maxval(controls,dim=2)+res))cycle
      call invert_endpoint(nodes,elements,e,endpoint%x(:2),travel,st,res,ok)
      if(.not.ok)cycle
      hit%element=e;hit%st=st;hit%wall=endpoint%x(:2);hit%lambda=1.d0
      hit%status=INTERIOR_ENDPOINT;hit%fallback_used=.true.;hit%max_residual=max(hit%max_residual,res)
      hit%reason='verified no external exit and geometrically inverted interior endpoint';return
    enddo
    hit%reason='boundary fallback: no exit but endpoint cell remains unresolved'
  endif
contains
  subroutine cubic(e,s,c)
    integer,intent(in)::e,s
    real(qp),intent(out)::c(2,0:3)
    integer::a,b,k,va,vb
    real(qp)::p(2),q(2),da(2),db(2)
    select case(s)
    case(1);a=1;b=2;k=2
    case(2);a=2;b=3;k=3
    case(3);a=4;b=3;k=2
    case(4);a=1;b=4;k=3
    end select
    va=elements%element(e)%vertex(a);vb=elements%element(e)%vertex(b)
    p=real(nodes%node(va)%x(1,1,:),qp)*real(elements%element(e)%size(a,1),qp)
    q=real(nodes%node(vb)%x(1,1,:),qp)*real(elements%element(e)%size(b,1),qp)
    da=3*real(nodes%node(va)%x(1,k,:),qp)*real(elements%element(e)%size(a,k),qp)
    db=-3*real(nodes%node(vb)%x(1,k,:),qp)*real(elements%element(e)%size(b,k),qp)
    c(:,0)=p;c(:,1)=da;c(:,2)=3*(q-p)-2*da-db;c(:,3)=2*(p-q)+da+db
  end subroutine
  subroutine path(l,p,v)
    real(qp),intent(in)::l
    real(qp),intent(out)::p(2),v(2)
    p=[sqrt((r0+l*d(1))**2+(l*d(3))**2),x0(2)+l*d(2)]
    v=[((r0+l*d(1))*real(event%v_old(1),qp)+l*d(3)*real(event%v_old(3),qp))/p(1),real(event%v_old(2),qp)]
  end subroutine
  subroutine add_candidate(u,l)
    real(qp),intent(in)::u,l
    if(l<0.or.l>1)return
    if(nc>0)then
      if(any(abs(ll(:nc)-l)<=256*epsilon(1._qp).and.abs(uu(:nc)-u)<=256*epsilon(1._qp)))return
    endif
    nc=nc+1;ll(nc)=l;uu(nc)=u
  end subroutine
  subroutine signed_curve(l,seed,c,side,signed,ok)
    real(qp),intent(in)::l,seed,c(2,0:3)
    integer,intent(in)::side
    real(qp),intent(out)::signed
    logical,intent(out)::ok
    real(qp)::p(2),v(2),w(2),t(2),dd(2),n(2),z,dz,den,f,err,fe,de
    real(qp)::tn,wn,dn,third,r,tm,dm,wm,lower,du,dx,normal_error,distance_error,budget
    integer::it
    ok=.false.;signed=0;z=seed;call path(l,p,v)
    if(.not.all(ieee_is_finite([p,reshape(c,[8]),z])))return
    if(z<0.or.z>1)return
    budget=real(local_intersection_tolerance(nodes,elements,e),qp)/4
    if(budget<=0)return
    ! Absolute evaluation allowance in metres, including cancellation in p-w.
    ! Parameter increments alone cannot measure accuracy on a very short face.
    err=4096*epsilon(1._qp)*max(1._qp,maxval(abs(p)),sum(abs(c)))
    third=norm2(6*c(:,3))
    do it=1,24
      w=[eval_q(c(1,:),3,z),eval_q(c(2,:),3,z)]
      t=c(:,1)+2*z*c(:,2)+3*z*z*c(:,3);dd=2*c(:,2)+6*z*c(:,3)
      tn=norm2(t);wn=norm2(p-w);dn=norm2(dd)
      f=dot_product(p-w,t);den=dot_product(t,t)-dot_product(p-w,dd)
      fe=err*(tn+wn+err);de=err*(2*tn+dn+wn+2*err)
      if(.not.all(ieee_is_finite([f,den,fe,de,tn,wn,dn])))return
      if(den<=de.or.tn<=err)return
      dz=f/den
      ! Certify a unique local minimum: F=(w-p).t, F'=den and
      ! |F''| <= 3|t||w''|+|w-p||w'''| on [z-r,z+r].
      ! Taylor bounds and the residual bracket give |z-z_min| <= du.
      r=4*(abs(f)+fe)/(den-de)+128*epsilon(1._qp)
      if(r<min(z,1-z))then
        tm=tn+dn*r+third*r*r/2
        dm=dn+third*r
        wm=wn+tn*r+dn*r*r/2+third*r*r*r/6
        lower=den-de-r*(3*tm*dm+wm*third)
        if(lower>0.and.r*lower>abs(f)+fe)then
          du=(abs(f)+fe)/lower
          dx=tm*du+err
          ! Bound the change in the unit normal separately from position.
          ! The contact budget does NOT suffice to classify tiny signed probes.
          normal_error=dm*du+err
          if(normal_error<tn/2)then
            normal_error=2*normal_error/(tn-normal_error)
            distance_error=dx+(wn+err)*normal_error+err
            n=[-t(2),t(1)];if(side==3.or.side==4)n=-n
            signed=dot_product(p-w,n)/tn
            if(all(ieee_is_finite([dx,distance_error,signed])).and. &
               dx<=budget.and.distance_error<abs(signed)/16)then
              ok=.true.;return
            endif
          endif
        endif
      endif
      ! A stalled iterate is accepted only by the residual/conditioning/sign
      ! certificate above. Exhaustion or an uncertified stall remains unresolved.
      if(.not.ieee_is_finite(dz).or.z+dz==z)return
      z=z+dz
      if(z<0.or.z>1)return
    enddo
  end subroutine
end subroutine boundary_fallback

pure function eval_q(c,n,x)result(y)
  integer,parameter::qp=selected_real_kind(30,300)
  integer,intent(in)::n
  real(qp),intent(in)::c(0:),x
  real(qp)::y
  integer::i
  y=c(n)
  do i=n-1,0,-1
    y=y*x+c(i)
  enddo
end function

recursive subroutine roots_q(coeff,degree,roots,nroots,ok)
  integer,parameter::qp=selected_real_kind(30,300)
  integer,intent(in)::degree
  real(qp),intent(in)::coeff(0:6)
  real(qp),intent(out)::roots(16)
  integer,intent(out)::nroots
  logical,intent(out)::ok
  real(qp)::c(0:6),deriv(0:6),critical(16),bounds(18),scale,tol,a,b,fa,fb,m,fm
  integer::n,nc,i,it,j
  logical::child
  roots=0;nroots=0;ok=.false.;n=degree;scale=maxval(abs(coeff(:n)))
  if(scale==0)return
  c=coeff/scale
  do while(n>0)
    if(c(n)/=0)exit
    n=n-1
  enddo
  tol=128*epsilon(1._qp)*sum(abs(c(:n)))
  if(n==0)then
    ok=abs(c(0))>tol;return
  endif
  deriv=0
  do i=1,n
    deriv(i-1)=i*c(i)
  enddo
  call roots_q(deriv,n-1,critical,nc,child)
  if(.not.child)return
  bounds(1)=0;bounds(2:nc+1)=critical(:nc);bounds(nc+2)=1;ok=.true.
  do i=1,nc+2
    a=bounds(i);fa=eval_q(c,n,a)
    if(abs(fa)<=tol)call insert(a)
    if(i==nc+2)exit
    b=bounds(i+1);fb=eval_q(c,n,b)
    if((fa<0.and.fb>0).or.(fa>0.and.fb<0))then
      do it=1,128
        m=(a+b)/2;fm=eval_q(c,n,m)
        if(fm==0.or.b-a<=8*epsilon(1._qp))exit
        if((fa<0).eqv.(fm<0))then
          a=m;fa=fm
        else
          b=m
        endif
      enddo
      call insert(m)
    endif
  enddo
contains
  subroutine insert(x)
    real(qp),intent(in)::x
    if(any(abs(roots(:nroots)-x)<=256*epsilon(1._qp)))return
    if(nroots==16)then
      ok=.false.;return
    endif
    j=nroots
    do while(j>0)
      if(roots(j)<=x)exit
      roots(j+1)=roots(j);j=j-1
    enddo
    roots(j+1)=x;nroots=nroots+1
  end subroutine
end subroutine roots_q

end module mod_gk_neutral_trajectory
