program test_poisson_harmonics
  use mod_poisson_element_kernel, only: accumulate_poisson_blocks, scatter_poisson_harmonics
  implicit none
  integer, parameter :: b=3, nc=5, np=64
  integer :: i, j, c, d, p
  integer :: modes(nc)
  character(len=3) :: types(nc)
  real*8 :: value(b), dx(b), dy(b), lap(b), a(b,b), kphi(b,b), expected_a(b,b), expected_b(b,b)
  real*8 :: matrix(b*nc,b*nc), sampled(b*nc,b*nc), hz(nc), hzp(nc), phi, err, cross, sampled_cross
  real*8 :: weight, r, jac, factor, bb2, psix, psiy, fh, fp, bi, bj

  modes = (/0,6,6,10,10/)
  types = (/'cos','cos','sin','cos','sin'/)
  value = (/0.7d0,-0.3d0,1.1d0/)
  dx    = (/0.2d0, 0.8d0,-0.4d0/)
  dy    = (/0.9d0,-0.1d0, 0.5d0/)
  lap   = (/0.4d0, 0.6d0,-0.2d0/)
  weight=0.73d0; r=2.4d0; jac=1.3d0; factor=1.7d0
  bb2=3.1d0; psix=0.35d0; psiy=-0.27d0; fh=0.08d0; fp=0.11d0

  a=0.d0; kphi=0.d0
  call accumulate_poisson_blocks(weight,r,jac,factor,bb2,psix,psiy,fh,fp,value,dx,dy,lap,a,kphi)

  do j=1,b
    bj=(dx(j)*psiy-dy(j)*psix)/r
    do i=1,b
      bi=(dx(i)*psiy-dy(i)*psix)/r
      expected_a(i,j)=weight*jac*r*(factor*(dx(i)*dx(j)+dy(i)*dy(j)) &
           +fh*lap(i)*lap(j)+fp*bi*bj/bb2)
      expected_b(i,j)=weight*jac/r*factor*value(i)*value(j)
    enddo
  enddo
  if (maxval(abs(a-expected_a)).gt.2.d-15) error stop 'Poisson A block mismatch.'
  if (maxval(abs(kphi-expected_b)).gt.2.d-15) error stop 'Poisson B block mismatch.'

  matrix=0.d0
  call scatter_poisson_harmonics(a,kphi,modes,types,np,matrix)

  ! Independently evaluate the old plane-sampled direct formulation.
  sampled=0.d0
  do p=1,np
    phi=2.d0*acos(-1.d0)*real(p-1,8)/real(np,8)
    hz=(/1.d0,cos(6.d0*phi),sin(6.d0*phi),cos(10.d0*phi),sin(10.d0*phi)/)
    hzp=(/0.d0,-6.d0*sin(6.d0*phi),6.d0*cos(6.d0*phi), &
                  -10.d0*sin(10.d0*phi),10.d0*cos(10.d0*phi)/)
    do j=1,b
      do i=1,b
        do d=1,nc
          do c=1,nc
            sampled((i-1)*nc+c,(j-1)*nc+d)=sampled((i-1)*nc+c,(j-1)*nc+d) &
              +a(i,j)*hz(c)*hz(d)+kphi(i,j)*hzp(c)*hzp(d)
          enddo
        enddo
      enddo
    enddo
  enddo

  err=maxval(abs(matrix-sampled))
  cross=0.d0
  sampled_cross=0.d0
  do j=1,b
    do i=1,b
      do d=1,nc
        do c=1,nc
          if (c.ne.d) then
            cross=max(cross,abs(matrix((i-1)*nc+c,(j-1)*nc+d)))
            sampled_cross=max(sampled_cross,abs(sampled((i-1)*nc+c,(j-1)*nc+d)))
          endif
        enddo
      enddo
    enddo
  enddo
  if (err.gt.2.d-12) error stop 'Analytic and sampled Poisson matrices differ.'
  if (cross.ne.0.d0) error stop 'Analytic Poisson matrix contains harmonic coupling.'
  if (sampled_cross.gt.2.d-12) error stop 'Sampled Poisson harmonic coupling exceeds roundoff.'
  write(*,'(A,ES12.4)') 'maximum analytic/sampled difference: ',err
  write(*,'(A,ES12.4)') 'maximum sampled cross-harmonic entry: ',sampled_cross
end program test_poisson_harmonics
