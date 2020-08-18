module data_mod
    use grid_function_mod
    use iso_c_binding, only : c_loc, c_ptr
    implicit none

    type :: grid_function_container_t
        class(grid_function_t), pointer :: obj
    end type grid_function_container_t

    type, public :: data_t
        integer(kind=4), public :: num_grid_functions, mpi_dtype
        type(grid_function_container_t), private :: grid_functions(20)
    contains
        procedure :: push           => push_data
        procedure :: deallocate     => deallocate_data
        procedure :: message_length => message_length_data
        procedure :: get_sub_array  => get_sub_array_data
        procedure :: put_sub_array  => put_sub_array_data
    end type data_t

    private

    contains

    subroutine push_data(this, grid_function)
        class(data_t), intent(inout) :: this
        class(grid_function_t), target, intent(in) :: grid_function

        this.num_grid_functions = this.num_grid_functions + 1
        allocate(this.grid_functions(this.num_grid_functions).obj, source=grid_function)
    end subroutine push_data

    subroutine deallocate_data(this)
        class(data_t), intent(inout) :: this
        
        integer(kind=4) :: i
        do i = 1, this.num_grid_functions
            deallocate(this.grid_functions(i).obj)
        end do
    end subroutine deallocate_data

    subroutine set_mpi_dtype_data(this, mpi_dtype)
        class(data_t), intent(inout) :: this
        integer(kind=4) :: mpi_dtype
        this.mpi_dtype = mpi_dtype
    end subroutine set_mpi_dtype_data

    function message_length_data(this) result(mlength)
        class(data_t), intent(in) :: this
        integer(kind=4) :: mlength
        integer(kind=4) :: i

        mlength = 0
        do i = 1, this.num_grid_functions
            mlength = mlength + this.grid_functions(i).obj.message_length()
        end do
    end function message_length_data

    subroutine get_sub_array_data(this, buf, idx, is, ie, js, je)
        class(data_t)  , intent(in) :: this
        integer(kind=4), intent(in)    :: is, ie, js, je
        integer(kind=8), intent(inout) :: buf(:)
        integer(kind=4), intent(inout) :: idx
        
        integer(kind=4) :: i, idx_max
        type(c_ptr) :: cptr

        idx_max = size(buf)
        cptr = c_loc(buf(1))

        do i = 1, this.num_grid_functions
            call this.grid_functions(i).obj.get_sub_array(cptr, idx, idx_max, is, ie, js, je)
        end do

    end subroutine get_sub_array_data

    subroutine put_sub_array_data(this, buf, idx, is, ie, js, je)
        class(data_t)  , intent(inout) :: this
        integer(kind=4), intent(in)    :: is, ie, js, je
        integer(kind=8), intent(in)    :: buf(:)
        integer(kind=4), intent(inout) :: idx
        
        integer(kind=4) :: i, idx_max
        type(c_ptr) :: cptr

        idx_max = size(buf)
        cptr = c_loc(buf(1))

        do i = 1, this.num_grid_functions
            call this.grid_functions(i).obj.put_sub_array(cptr, idx, idx_max, is, ie, js, je)
        end do

    end subroutine put_sub_array_data

end module