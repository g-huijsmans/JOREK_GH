program verifier_geometry_test
use mpi
use data_structure
use mod_gk_neutral_trajectory
use mod_gk_neutrals,only:push_gk_neutral_ballistic
use mod_particle_types,only:particle_kinetic_leapfrog
use mod_interp,only:interp_RZ
implicit none
type(type_node_list)::nodes
type(type_element_list)::elements
type(side_correction)::hit
type(particle_kinetic_leapfrog)::p
integer::i,j,k,rank,np,ierr,bad,total,caseid,ntangent,fall
real*8::s,t,f,fp,sg,tt,x(3),v(3),st(2),dt,z,expected,rs,rt,zs,zt
logical::ok
call MPI_INIT(ierr);call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr);call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)
call init_threads();call init_node_list(nodes,8,32,n_var);elements%n_elements=1
! R=1+s; Z=t*(1+3*(s-.5)^2). Exact stored cubic face minimum Z=1.
elements%element(1)%vertex=[1,2,3,4];elements%element(1)%size=1;elements%element(1)%neighbours=0
 do i=1,4
  s=0;if(i==2.or.i==3)s=1
  t=0;if(i>=3)t=1
  sg=1-2*s;tt=1-2*t;f=1+3*(s-.5d0)**2;fp=6*(s-.5d0)
  nodes%node(i)%x=0;nodes%node(i)%boundary=1
  nodes%node(i)%x(1,1,:)=[1+s,t*f]
  nodes%node(i)%x(1,2,:)=[sg/3,sg*t*fp/3]
  nodes%node(i)%x(1,3,:)=[0.d0,tt*f/3]
  nodes%node(i)%x(1,4,:)=[0.d0,sg*tt*fp/9]
 enddo
bad=0;total=0;ntangent=0;fall=0
! Replicated read-only geometry, independent stack-local hits under OpenMP.
!$omp parallel do private(caseid,x,v,st,p,hit,z,expected,dt,ok) reduction(+:bad,total,ntangent,fall)
do k=1,160
 caseid=mod(k-1,10)+1;dt=1.d0;x=[1.25d0,.5d0,0.7d0];v=[.5d0,0.d0,0.d0];expected=-1
 select case(caseid)
 case(1) ! ordinary interior
 case(2) ! exactly tangent touch, no outward exit
  x(2)=1.d0
 case(3) ! curved crossing and re-entry, endpoint inside
  x(2)=1.03d0;expected=(.5d0-sqrt(.03d0/3)-.25d0)/.5d0
 case(4) ! very shallow pair of crossings
  x(2)=1.d0+2.d0**(-40);expected=(.5d0-sqrt(2.d0**(-40)/3)-.25d0)/.5d0
 case(5) ! inward start exactly on bottom wall
  x(2)=0;v=[0.d0,.2d0,0.d0]
 case(7) ! planar physical side
  v=[1.d0,0.d0,0.d0];expected=.75d0
 case(8) ! roots separated below binary64 geometric uncertainty
  x(2)=1.d0+epsilon(1.d0);expected=-2
 case(9) ! physical R trajectory has toroidal curvature, not a chord
  v=[0.d0,0.d0,1.d0];dt=2.d0;expected=sqrt(4.d0-x(1)**2)/2.d0
 case(10) ! remaining travel indistinguishable from a wall touch
  x=[2.d0,.5d0,.7d0];v=[epsilon(1.d0),0.d0,0.d0];expected=-2
 case(6) ! outward zero-time contact
  x(2)=0;v=[0.d0,-.2d0,0.d0];expected=0
 end select
 st=[x(1)-1,x(2)/(1+3*(x(1)-1.5d0)**2)]
 p%x=x;p%v=v;call push_gk_neutral_ballistic(p,dt)
 call verify_neutral_trajectory(nodes,elements,1,st,x,v,dt,hit,p%x,p%v)
 total=total+1;ntangent=ntangent+hit%fallback_tangencies
 if(hit%fallback_attempted)fall=fall+1
 ok=hit%status==INTERIOR_ENDPOINT
 if(expected==-2)ok=hit%status==UNRESOLVED
 if(expected>=0)ok=hit%status==WALL_CROSSING.and.abs(hit%lambda-expected)<1.d-8
 if(.not.ok)then
  bad=bad+1
  !$omp critical
  write(*,*)'FAILED_GEOMETRY',rank,caseid,hit%status,hit%lambda,expected,trim(hit%reason)
  !$omp end critical
 endif
enddo
!$omp end parallel do
! Duplicate physical faces with distinct node IDs do not define a unique
! earliest boundary element. They must not be accepted by arbitrary tie-breaking.
elements%n_elements=2;elements%element(2)=elements%element(1)
elements%element(2)%vertex=[5,6,7,8]
do i=1,4
 nodes%node(i+4)=nodes%node(i)
enddo
x=[1.25d0,1.d0+2.d0**(-40),.7d0];v=[.5d0,0.d0,0.d0];dt=1.d0
st=[.25d0,x(2)/(1+3*.25d0**2)]
call verify_neutral_trajectory(nodes,elements,1,st,x,v,dt,hit)
total=total+1
if(hit%status/=UNRESOLVED.or..not.hit%fallback_attempted)bad=bad+1
write(*,*)'COINCIDENT_PHYSICAL_FACES',rank,hit%status,trim(hit%reason)
write(*,*)'GEOMETRY_RESULT',rank,total,bad,fall,ntangent
call MPI_ALLREDUCE(MPI_IN_PLACE,bad,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
if(bad/=0)call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
call MPI_FINALIZE(ierr)
end program
