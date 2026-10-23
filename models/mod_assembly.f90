module mod_assembly
implicit none
!> This module provides routines to factorize code in boundary conditions from each models

contains
  !> Eliminate selected homogeneous boundary DOFs from a correction solve.
  !! values_new = values + delta, hence delta_c = -current_c. Remove both
  !! rows and columns and apply b_f <- b_f + A_fc*current_c before zeroing
  !! A_fc. Each rank visits only its owned rows, including columns owned by
  !! another rank. The unchanged block sparsity is shared by the product
  !! matrix and harmonic preconditioners.
  subroutine boundary_conditions_zero_dofs(selected, current, k, in, &
                                           index_min, index_max, rhs_loc, a_mat)
    use mod_parameters, only: n_var
    use mod_integer_types
    use data_structure, only: type_SP_MATRIX
    logical, intent(in) :: selected(:)
    real*8, intent(in) :: current(:)
    integer, intent(in) :: k, in, index_min, index_max
    real*8, intent(inout) :: rhs_loc(*)
    type(type_SP_MATRIX), intent(inout) :: a_mat
    integer :: nt, block, slot, i, local_i, j
    integer(kind=int_all) :: q, col, pos, entry, row

    if (k == 0 .or. in < a_mat%i_tor_min .or. in > a_mat%i_tor_max) return
    nt = a_mat%i_tor_max - a_mat%i_tor_min + 1
    block = n_var*nt
    slot = (k-1)*nt + in-a_mat%i_tor_min + 1
    do i=index_min,index_max
      local_i = i-index_min+1
      do q=1,a_mat%ijA_size(local_i)
        col = a_mat%irn_jcn(local_i,q)
        if (.not. selected(i) .and. .not. selected(col)) cycle
        pos = a_mat%ijA_index(local_i,q)
        if (selected(col)) then
          do j=1,block
            entry = pos-1 + (j-1)*block + slot
            row = int(block,int_all)*(i-1) + j
            if (.not. (selected(i) .and. j == slot)) &
              rhs_loc(row) = rhs_loc(row) + a_mat%val(entry)*current(col)
            a_mat%val(entry) = 0.d0
          enddo
        endif
        if (selected(i)) a_mat%val(pos+(slot-1)*block:pos+slot*block-1) = 0.d0
      enddo
      if (selected(i)) then
        call boundary_conditions_add_one_entry(i,k,in,i,k,in,1.d0,index_min,index_max,a_mat)
        call boundary_conditions_add_RHS(i,k,in,index_min,index_max,rhs_loc,-current(i), &
                                         a_mat%i_tor_min,a_mat%i_tor_max)
      endif
    enddo
  end subroutine boundary_conditions_zero_dofs

  !>
  !! Subroutine: boundary_conditions_add_one_entry
  !!
  !! Add one entry to the product and/or harminic matrix if index_node is local.
  !!
  !! @param index_node  row node index
  !! @param k           row var index
  !! @param in          row tor index
  !! @param index_node2 col node index
  !! @param k2          col var index
  !! @param in2         col tor index
  !! @param zbig        value
  !! @param use_murge   Use murge interface.
  !! @param use_murge_element Murge interface with elementary matrices.
  !! @param index_min   Minimal local element index
  !! @param index_max   Maximal local element index
  !!
  subroutine boundary_conditions_add_one_entry( &
       &   index_node,  k,  in,                 &
       &   index_node2, k2, in2,                &
       &   zbig, index_min, index_max, a_mat)
    use mod_parameters
    use mod_locate_irn_jcn
    use mod_integer_types
    use data_structure, only: type_SP_MATRIX

    integer,               intent(in)                 :: k,  in
    integer,               intent(in)                 :: k2, in2
    integer,               intent(in)                 :: index_node
    integer,               intent(in)                 :: index_node2
    real*8,                intent(in)                 :: zbig
    integer,               intent(in)                 :: index_min, index_max
    type(type_SP_MATRIX)                              :: a_mat
    
    logical                                           :: is_local
    integer(kind=int_all)                             :: ija_position, ilarge_vp
    integer                                           :: n_tor_local

    if ( (k==0) .or. (k2==0) ) return ! ignore calls for model family extensions not in use (variable number zero)

    n_tor_local = a_mat%i_tor_max - a_mat%i_tor_min +1
    if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then

       call locate_irn_jcn(index_node,index_node2,index_min,index_max,ijA_position,a_mat)
                             
       !-------- index dans A_mat
       ilarge_vp  = ijA_position  - 1 + ((k-1)*n_tor_local + in-a_mat%i_tor_min ) * n_var*n_tor_local + (k2-1)*n_tor_local + in2&
                    -a_mat%i_tor_min + 1 
                               
                             
       a_mat%irn(ilarge_vp) =  n_tor_local * n_var * (index_node -1) + (k -1)*n_tor_local + in - a_mat%i_tor_min + 1
       a_mat%jcn(ilarge_vp) =  n_tor_local * n_var * (index_node2-1) + (k2-1)*n_tor_local + in2 - a_mat%i_tor_min + 1
       a_mat%val(ilarge_vp) = ZBIG
    endif
  end subroutine boundary_conditions_add_one_entry

  !>
  !! Subroutine: boundary_conditions_add_RHS
  !!
  !! Add one entry to the product and/or harminic matrix if index_node is local.
  !!
  !! @param index_node  row node index
  !! @param k           row var index
  !! @param in          row tor index
  !! @param index_node2 col node index
  !! @param k2          col var index
  !! @param in2         col tor index
  !! @param zbig        value
  !! @param index_min   Minimal local element index
  !! @param index_max   Maximal local element index
  !!
  subroutine boundary_conditions_add_RHS( &
       &   index_node, k, in,             &
       &   index_min, index_max,          &
       &   rhs_loc,  val,                 &
       &   i_tor_min, i_tor_max)
    use mod_parameters
    use mod_integer_types
    integer,               intent(in)    :: k,  in
    integer,               intent(in)    :: index_node
    integer,               intent(in)    :: index_min, index_max
    integer,               intent(in)    :: i_tor_min, i_tor_max 
    real*8,                intent(in)    :: val
    real*8,                intent(inOUT) :: rhs_loc(*) 
    integer                              :: n_tor_local 
    logical                              :: is_local

    if ( (k==0) ) return ! ignore calls for model family extensions not in use (variable number zero)

    n_tor_local = i_tor_max - i_tor_min +1
    if ((index_node .ge. index_min) .and. (index_node .le. index_max)) then
       RHS_loc(n_tor_local*n_var * (index_node-1) + (k-1)*n_tor_local + in - i_tor_min + 1) = val
    endif
  end subroutine boundary_conditions_add_RHS
end module mod_assembly
