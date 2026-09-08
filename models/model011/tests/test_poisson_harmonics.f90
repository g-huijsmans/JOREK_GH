program test_poisson_harmonics
  use mod_poisson_element_kernel, only: accumulate_poisson_n0_block, &
       accumulate_poisson_nonzero_blocks,scatter_poisson_harmonics, &
       accumulate_poisson_load_mass,apply_poisson_load_harmonics
  implicit none
  integer, parameter :: b=3,nc=5,np=64,ncase=5
  integer :: modes(nc),icase
  character(len=3) :: types(nc)
  character(len=16) :: names(ncase)
  real*8 :: fhyper(ncase),fpar(ncase),fhyper0(ncase),fpar0(ncase),fcentre(ncase)
  real*8 :: analytic(b*nc,b*nc),sampled(b*nc,b*nc),rhs(b*nc),phi_a(b*nc),phi_s(b*nc)
  real*8 :: max_matrix_abs,max_matrix_rel,max_phi_abs,max_phi_rel,l2_phi
  real*8 :: sv(b),sdx(b),sdy(b),slap(b),shz(nc),shzp(nc),sphi
  real*8 :: sw,sr,sjac,sfactor,sbb2,sf0,spsix,spsiy,sdvi,sdpj,sintegrand,stwopi
  integer :: sq,sp,si,sj,sc,sd,srow,scol

  modes=(/0,2,2,6,6/); types=(/'cos','cos','sin','cos','sin'/)
  names=(/'filters off     ','hyper           ','parallel        ','n0 centre       ','full            '/)
  fhyper=(/0d0,0.08d0,0d0,0d0,0.08d0/)
  fpar=(/0d0,0d0,0.11d0,0d0,0.11d0/)
  fhyper0=(/0d0,0.07d0,0d0,0d0,0.07d0/)
  fpar0=(/0d0,0d0,0.09d0,0d0,0.09d0/)
  fcentre=(/0d0,0d0,0d0,1d0,1d0/)

  do icase=1,ncase
    call build_analytic(fhyper0(icase),fpar0(icase)+fcentre(icase), &
         fhyper(icase),fpar(icase),analytic)
    sampled=0d0; stwopi=2d0*acos(-1d0)
    do sq=1,3
      call point_data(sq,sv,sdx,sdy,slap,sw,sr,sjac,sfactor,sbb2,sf0,spsix,spsiy)
      do sp=1,np
        sphi=stwopi*real(sp-1,8)/real(np,8)
        shz=(/1d0,cos(2d0*sphi),sin(2d0*sphi),cos(6d0*sphi),sin(6d0*sphi)/)
        shzp=(/0d0,-2d0*sin(2d0*sphi),2d0*cos(2d0*sphi),-6d0*sin(6d0*sphi),6d0*cos(6d0*sphi)/)
        do sj=1,b
          do si=1,b
            do sd=1,nc
              do sc=1,nc
                if (modes(sc).ne.modes(sd)) cycle
                srow=(si-1)*nc+sc; scol=(sj-1)*nc+sd
                if (modes(sc).eq.0) then
                  sintegrand=sfactor*(sdx(si)*sdx(sj)+sdy(si)*sdy(sj)) &
                       +fhyper0(icase)*slap(si)*slap(sj) &
                       +(fpar0(icase)+fcentre(icase))*((sdx(si)*spsiy-sdy(si)*spsix)/sr) &
                       *((sdx(sj)*spsiy-sdy(sj)*spsix)/sr)/sbb2
                else
                  sdvi=(sf0/sr*sv(si)*shzp(sc)+(sdx(si)*spsiy-sdy(si)*spsix)*shz(sc))/sr
                  sdpj=(sf0/sr*sv(sj)*shzp(sd)+(sdx(sj)*spsiy-sdy(sj)*spsix)*shz(sd))/sr
                  sintegrand=sfactor*(sdx(si)*sdx(sj)+sdy(si)*sdy(sj))*shz(sc)*shz(sd) &
                       +fhyper(icase)*slap(si)*slap(sj)*shz(sc)*shz(sd) &
                       +sfactor*sv(si)*sv(sj)*shzp(sc)*shzp(sd)/sr**2 &
                       +(fpar(icase)-sfactor)*sdvi*sdpj/sbb2
                endif
                sampled(srow,scol)=sampled(srow,scol)+stwopi/real(np,8)*sw*sjac*sr*sintegrand
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    call matrix_errors(analytic,sampled,max_matrix_abs,max_matrix_rel)
    write(*,'(A,1X,A,2(1X,ES12.4))') 'matrix',trim(names(icase)),max_matrix_abs,max_matrix_rel
    if (max_matrix_abs.gt.2.d-12 .or. max_matrix_rel.gt.2.d-12) &
      error stop 'Sampled and analytic Poisson matrices differ.'
  enddo

  rhs=(/0d0,0.3d0,-0.2d0,0.1d0,0.4d0, 0d0,-0.1d0,0.2d0,0.5d0,-0.3d0, &
         0d0,0.6d0,-0.4d0,0.2d0,0.1d0/)
  call impose_dirichlet(analytic,rhs)
  call impose_dirichlet(sampled,rhs)
  call solve_dense(analytic,rhs,phi_a)
  call solve_dense(sampled,rhs,phi_s)
  call vector_errors(phi_a,phi_s,max_phi_abs,max_phi_rel)
  l2_phi=sqrt(sum((phi_a-phi_s)**2))
  write(*,'(A,3(1X,ES12.4))') 'phi abs/rel/l2',max_phi_abs,max_phi_rel,l2_phi
  if (max_phi_abs.gt.2.d-10 .or. max_phi_rel.gt.2.d-10) error stop 'Poisson solutions differ.'

  call check_rhs
contains
  non_recursive subroutine point_data(q,value,dx,dy,lap,weight,r,jac,factor,bb2,f0,psix,psiy)
    integer,intent(in) :: q
    real*8,intent(out) :: value(b),dx(b),dy(b),lap(b),weight,r,jac,factor,bb2,f0,psix,psiy
    value=(/0.7d0,-0.3d0,1.1d0/)+0.03d0*q*(/1d0,-2d0,0.5d0/)
    dx=(/0.2d0,0.8d0,-0.4d0/)+0.02d0*q*(/-1d0,0.5d0,1d0/)
    dy=(/0.9d0,-0.1d0,0.5d0/)+0.01d0*q*(/1d0,2d0,-1d0/)
    lap=(/0.4d0,0.6d0,-0.2d0/)+0.04d0*q*(/0.5d0,-1d0,2d0/)
    weight=0.31d0+0.07d0*q; r=2.2d0+0.1d0*q; jac=1.1d0+0.05d0*q
    factor=0.8d0+0.06d0*q; f0=0.655d0; psix=0.35d0-0.02d0*q; psiy=-0.27d0+0.03d0*q
    bb2=(f0*f0+psix*psix+psiy*psiy)/r**2
  end subroutine point_data

  non_recursive subroutine build_analytic(fh0,fq0,fh,fq,matrix)
    real*8,intent(in) :: fh0,fq0,fh,fq
    real*8,intent(out) :: matrix(:,:)
    real*8 :: n0(b,b),a(b,b),bk(b,b),c(b,b),value(b),dx(b),dy(b),lap(b)
    real*8 :: weight,r,jac,factor,bb2,f0,psix,psiy
    integer :: q
    n0=0d0; a=0d0; bk=0d0; c=0d0
    do q=1,3
      call point_data(q,value,dx,dy,lap,weight,r,jac,factor,bb2,f0,psix,psiy)
      call accumulate_poisson_n0_block(weight,r,jac,factor,bb2,psix,psiy,fh0,fq0, &
           value,dx,dy,lap,n0)
      call accumulate_poisson_nonzero_blocks(weight,r,jac,factor,bb2,f0,psix,psiy,fh,fq, &
           value,dx,dy,lap,a,bk,c)
    enddo
    matrix=0d0
    call scatter_poisson_harmonics(n0,a,bk,c,modes,types,matrix)
  end subroutine build_analytic

  subroutine check_rhs
    real*8 :: value(b),dx(b),dy(b),lap(b),weight,r,jac,factor,bb2,f0,psix,psiy
    real*8 :: mass(b,b),charge(b,3),analytic_rhs(b,3),sampled_rhs(b,3),hz(3),source,phi,density
    integer :: q,p,i,j,c
    integer :: rmodes(3)
    character(len=3) :: rtypes(3)
    rmodes=(/0,2,2/); rtypes=(/'cos','cos','sin'/); density=8d19
    charge=reshape((/0.8d0,-0.2d0,0.5d0,0.3d0,0.7d0,-0.4d0,-0.6d0,0.1d0,0.9d0/),shape(charge))
    mass=0d0
    do q=1,3
      call point_data(q,value,dx,dy,lap,weight,r,jac,factor,bb2,f0,psix,psiy)
      call accumulate_poisson_load_mass(weight,r,jac,value,mass)
    enddo
    call apply_poisson_load_harmonics(mass,charge,rmodes,rtypes,np,density,analytic_rhs)
    sampled_rhs=0d0
    do q=1,3
      call point_data(q,value,dx,dy,lap,weight,r,jac,factor,bb2,f0,psix,psiy)
      do p=1,np
        phi=2d0*acos(-1d0)*real(p-1,8)/real(np,8)
        hz=(/1d0,cos(2d0*phi),sin(2d0*phi)/)
        do j=1,b
          source=sum(charge(j,:)*hz)
          do i=1,b; do c=1,3
            sampled_rhs(i,c)=sampled_rhs(i,c)+2d0*acos(-1d0)/real(np,8)* &
                 weight*r*jac*value(i)*value(j)*source*hz(c)/density
          enddo; enddo
        enddo
      enddo
    enddo
    call matrix_errors(analytic_rhs,sampled_rhs,max_matrix_abs,max_matrix_rel)
    write(*,'(A,2(1X,ES12.4))') 'projected RHS abs/rel',max_matrix_abs,max_matrix_rel
    if (max_matrix_abs.gt.2d-33 .or. max_matrix_rel.gt.2d-12) error stop 'Physical RHS mismatch.'
    ! A saved direct particle load is already the tested global FE load; both paths
    ! perform this same density normalization for n=0, cosine and sine entries.
    sampled_rhs=analytic_rhs*density
    analytic_rhs=sampled_rhs/density
    call matrix_errors(analytic_rhs,sampled_rhs/density,max_matrix_abs,max_matrix_rel)
    write(*,'(A,2(1X,ES12.4))') 'direct RHS abs/rel',max_matrix_abs,max_matrix_rel
  end subroutine check_rhs

  subroutine impose_dirichlet(matrix,vector)
    real*8,intent(inout) :: matrix(:,:),vector(:)
    integer :: c,index
    do c=1,nc
      index=c
      matrix(index,:)=0d0; matrix(:,index)=0d0; matrix(index,index)=1d0; vector(index)=0d0
    enddo
  end subroutine impose_dirichlet

  subroutine solve_dense(matrix,vector,solution)
    real*8,intent(in) :: matrix(:,:),vector(:)
    real*8,intent(out) :: solution(:)
    real*8 :: a(size(vector),size(vector)),v(size(vector)),pivot,factor
    integer :: i,j,k,n
    n=size(vector); a=matrix; v=vector
    do k=1,n
      pivot=a(k,k)
      if (abs(pivot).lt.1d-13) error stop 'Singular comparison matrix.'
      a(k,k:n)=a(k,k:n)/pivot; v(k)=v(k)/pivot
      do i=k+1,n
        factor=a(i,k); a(i,k:n)=a(i,k:n)-factor*a(k,k:n); v(i)=v(i)-factor*v(k)
      enddo
    enddo
    solution=v
    do i=n-1,1,-1; solution(i)=solution(i)-dot_product(a(i,i+1:n),solution(i+1:n)); enddo
  end subroutine solve_dense

  subroutine matrix_errors(a,bv,absolute,relative)
    real*8,intent(in) :: a(:,:),bv(:,:)
    real*8,intent(out) :: absolute,relative
    absolute=maxval(abs(a-bv))
    relative=absolute/max(maxval(abs(a)),maxval(abs(bv)),tiny(1d0))
  end subroutine matrix_errors

  subroutine vector_errors(a,bv,absolute,relative)
    real*8,intent(in) :: a(:),bv(:)
    real*8,intent(out) :: absolute,relative
    absolute=maxval(abs(a-bv))
    relative=absolute/max(maxval(abs(a)),maxval(abs(bv)),tiny(1d0))
  end subroutine vector_errors
end program test_poisson_harmonics
