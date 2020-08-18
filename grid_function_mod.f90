module grid_function_mod
    use iso_c_binding, only : c_ptr, c_f_pointer
    implicit none

    type, public, abstract :: grid_function_t
    contains 
        procedure(message_length_s), deferred :: message_length
        procedure(get_sub_array_s) , deferred :: get_sub_array
        procedure(put_sub_array_s) , deferred :: put_sub_array
    end type grid_function_t

    abstract interface 
        function message_length_s(this) result(mlength)
            import :: grid_function_t

            class(grid_function_t), intent(in)    :: this
            integer(kind=4) :: mlength
        end function message_length_s
    end interface

    abstract interface 
        subroutine get_sub_array_s(this, buf, idx, idx_max, is, ie, js, je)
            import :: grid_function_t, c_ptr

            class(grid_function_t), intent(in)    :: this
            integer(kind=4)       , intent(in)    :: idx_max, is, ie, js, je
            type(c_ptr)           , intent(in)    :: buf
            integer(kind=4)       , intent(inout) :: idx
        end subroutine get_sub_array_s
    end interface

    abstract interface 
        subroutine put_sub_array_s(this, buf, idx, idx_max, is, ie, js, je)
            import :: grid_function_t, c_ptr

            class(grid_function_t), intent(inout) :: this
            integer(kind=4)       , intent(in)    :: idx_max, is, ie, js, je
            type(c_ptr)           , intent(in)    :: buf
            integer(kind=4)       , intent(inout) :: idx
        end subroutine put_sub_array_s
    end interface

    type, public, extends(grid_function_t) :: grid_function_2d_r8_t
        real(kind=8), pointer, contiguous :: array(:,:)
        contains
            procedure :: init           => init_2d_r8
            procedure :: message_length => message_length_2d_r8
            procedure :: get_sub_array  => get_sub_array_2d_r8
            procedure :: put_sub_array  => put_sub_array_2d_r8
    end type grid_function_2d_r8_t

    type, public, extends(grid_function_t) :: grid_function_2d_r4_t
        real(kind=4), pointer, contiguous :: array(:,:)
        contains
            procedure :: init           => init_2d_r4
            procedure :: message_length => message_length_2d_r4
            procedure :: get_sub_array  => get_sub_array_2d_r4
            procedure :: put_sub_array  => put_sub_array_2d_r4
    end type grid_function_2d_r4_t

    type, public, extends(grid_function_t) :: grid_function_4d_r8_t
        real(kind=8), pointer, contiguous :: array(:,:,:,:)
        contains
            procedure :: init           => init_4d_r8
            procedure :: message_length => message_length_4d_r8
            procedure :: get_sub_array  => get_sub_array_4d_r8
            procedure :: put_sub_array  => put_sub_array_4d_r8
    end type grid_function_4d_r8_t

    type, public, extends(grid_function_t) :: grid_function_4d_r4_t
        real(kind=4), pointer, contiguous :: array(:,:,:,:)
        contains
            procedure :: init           => init_4d_r4
            procedure :: message_length => message_length_4d_r4
            procedure :: get_sub_array  => get_sub_array_4d_r4
            procedure :: put_sub_array  => put_sub_array_4d_r4
    end type grid_function_4d_r4_t
    
    type, public, extends(grid_function_t) :: grid_function_3d_r8_t
        real(kind=8), pointer, contiguous :: array(:,:,:)
        contains
            procedure :: init           => init_3d_r8
            procedure :: message_length => message_length_3d_r8
            procedure :: get_sub_array  => get_sub_array_3d_r8
            procedure :: put_sub_array  => put_sub_array_3d_r8
    end type grid_function_3d_r8_t

    type, public, extends(grid_function_t) :: grid_function_3d_r4_t
        real(kind=4), pointer, contiguous :: array(:,:,:)
        contains
            procedure :: init           => init_3d_r4
            procedure :: message_length => message_length_3d_r4
            procedure :: get_sub_array  => get_sub_array_3d_r4
            procedure :: put_sub_array  => put_sub_array_3d_r4
    end type grid_function_3d_r4_t    
    
    type, public, extends(grid_function_t) :: grid_function_3d_r8_depth_t
        real(kind=8)   , pointer, contiguous :: array(:,:,:)
        integer(kind=4), pointer, contiguous :: km2_l(:,:)
        contains
            procedure :: init           => init_3d_r8_depth
            procedure :: message_length => message_length_3d_r8_depth
            procedure :: get_sub_array  => get_sub_array_3d_r8_depth
            procedure :: put_sub_array  => put_sub_array_3d_r8_depth
    end type grid_function_3d_r8_depth_t

    type, public, extends(grid_function_t) :: grid_function_3d_r4_depth_t
        real(kind=4)   , pointer, contiguous :: array(:,:,:)
        integer(kind=4), pointer, contiguous :: km2_l(:,:)
        contains
            procedure :: init           => init_3d_r4_depth
            procedure :: message_length => message_length_3d_r4_depth
            procedure :: get_sub_array  => get_sub_array_3d_r4_depth
            procedure :: put_sub_array  => put_sub_array_3d_r4_depth
    end type grid_function_3d_r4_depth_t
    
    type, public, extends(grid_function_t) :: grid_function_4d_r8_depth_t
        real(kind=8)   , pointer, contiguous :: array(:,:,:,:)
        integer(kind=4), pointer, contiguous :: km2_l(:,:)
        contains
            procedure :: init           => init_4d_r8_depth
            procedure :: message_length => message_length_4d_r8_depth
            procedure :: get_sub_array  => get_sub_array_4d_r8_depth
            procedure :: put_sub_array  => put_sub_array_4d_r8_depth
    end type grid_function_4d_r8_depth_t

    type, public, extends(grid_function_t) :: grid_function_4d_r4_depth_t
        real(kind=4)   , pointer, contiguous :: array(:,:,:,:)
        integer(kind=4), pointer, contiguous :: km2_l(:,:)
        contains
            procedure :: init           => init_4d_r4_depth
            procedure :: message_length => message_length_4d_r4_depth
            procedure :: get_sub_array  => get_sub_array_4d_r4_depth
            procedure :: put_sub_array  => put_sub_array_4d_r4_depth
    end type grid_function_4d_r4_depth_t

    private

    contains

    subroutine init_2d_r8(this, array, is, ie, js, je)
        class(grid_function_2d_r8_t), intent(inout) :: this
        integer(kind=4)             , intent(in)    :: is, ie, js, je
        real(kind=8), target        , intent(in)    :: array(is:ie,js:je)

        this.array => array       
    end subroutine init_2d_r8

    subroutine init_2d_r4(this, array, is, ie, js, je)
        class(grid_function_2d_r4_t), intent(inout) :: this
        integer(kind=4)             , intent(in)    :: is, ie, js, je
        real(kind=4), target        , intent(in)    :: array(is:ie,js:je)

        this.array => array       
    end subroutine init_2d_r4

    subroutine init_3d_r8(this, array, is, ie, js, je, ks, ke)
        class(grid_function_3d_r8_t), intent(inout) :: this
        integer(kind=4)             , intent(in)    :: is, ie, js, je, ks, ke
        real(kind=8), target        , intent(in)    :: array(is:ie,js:je,ks:ke)

        this.array => array        
    end subroutine init_3d_r8

    subroutine init_3d_r4(this, array, is, ie, js, je, ks, ke)
        class(grid_function_3d_r4_t), intent(inout) :: this
        integer(kind=4)             , intent(in)    :: is, ie, js, je, ks, ke
        real(kind=4), target        , intent(in)    :: array(is:ie,js:je,ks:ke)

        this.array => array        
    end subroutine init_3d_r4

    subroutine init_4d_r8(this, array, ls, le, is, ie, js, je, ks, ke)
        class(grid_function_4d_r8_t), intent(inout) :: this
        integer(kind=4)             , intent(in)    :: ls, le, is, ie, js, je, ks, ke
        real(kind=8), target        , intent(in)    :: array(ls:le,is:ie,js:je,ks:ke)

        this.array => array        
    end subroutine init_4d_r8

    subroutine init_4d_r4(this, array, ls, le, is, ie, js, je, ks, ke)
        class(grid_function_4d_r4_t), intent(inout) :: this
        integer(kind=4)             , intent(in)    :: ls, le, is, ie, js, je, ks, ke
        real(kind=4), target        , intent(in)    :: array(ls:le,is:ie,js:je,ks:ke)

        this.array => array        
    end subroutine init_4d_r4
    
    subroutine init_3d_r8_depth(this, array, km2_l, is, ie, js, je, ks, ke)
        class(grid_function_3d_r8_depth_t), intent(inout) :: this
        integer(kind=4)                   , intent(in)    :: is, ie, js, je, ks, ke
        real(kind=8), target              , intent(in)    :: array(is:ie,js:je,ks:ke)
        integer(kind=4), target           , intent(in)    :: km2_l(is:ie,js:je)

        this.array => array        
        this.km2_l => km2_l
    end subroutine init_3d_r8_depth

    subroutine init_3d_r4_depth(this, array, km2_l, is, ie, js, je, ks, ke)
        class(grid_function_3d_r4_depth_t), intent(inout) :: this
        integer(kind=4)                   , intent(in)    :: is, ie, js, je, ks, ke
        real(kind=4), target              , intent(in)    :: array(is:ie,js:je,ks:ke)
        integer(kind=4), target           , intent(in)    :: km2_l(is:ie,js:je)

        this.array => array
        this.km2_l => km2_l
    end subroutine init_3d_r4_depth
    
    subroutine init_4d_r8_depth(this, array, km2_l, ls, le, is, ie, js, je, ks, ke)
        class(grid_function_4d_r8_depth_t), intent(inout) :: this
        integer(kind=4)                   , intent(in)    :: ls, le, is, ie, js, je, ks, ke
        real(kind=8), target              , intent(in)    :: array(ls:le,is:ie,js:je,ks:ke)
        integer(kind=4), target           , intent(in)    :: km2_l(is:ie,js:je)

        this.array => array        
        this.km2_l => km2_l
    end subroutine init_4d_r8_depth

    subroutine init_4d_r4_depth(this, array, km2_l, ls, le, is, ie, js, je, ks, ke)
        class(grid_function_4d_r4_depth_t), intent(inout) :: this
        integer(kind=4)                   , intent(in)    :: ls, le, is, ie, js, je, ks, ke
        real(kind=4), target              , intent(in)    :: array(ls:le,is:ie,js:je,ks:ke)
        integer(kind=4), target           , intent(in)    :: km2_l(is:ie,js:je)

        this.array => array        
        this.km2_l => km2_l
    end subroutine init_4d_r4_depth

    function message_length_2d_r8(this) result(mlength)
        class(grid_function_2d_r8_t), intent(in) :: this
        integer(kind=4) :: mlength

        mlength = 1
    end function message_length_2d_r8

    function message_length_2d_r4(this) result(mlength)
        class(grid_function_2d_r4_t), intent(in) :: this
        integer(kind=4) :: mlength

        mlength = 1
    end function message_length_2d_r4

    function message_length_3d_r8(this) result(mlength)
        class(grid_function_3d_r8_t), intent(in) :: this
        integer(kind=4) :: mlength

        mlength = size(this.array,3)
    end function message_length_3d_r8

    function message_length_3d_r4(this) result(mlength)
        class(grid_function_3d_r4_t), intent(in) :: this
        integer(kind=4) :: mlength

        mlength = size(this.array,3)
    end function message_length_3d_r4

    function message_length_4d_r8(this) result(mlength)
        class(grid_function_4d_r8_t), intent(in) :: this
        integer(kind=4) :: mlength

        mlength = size(this.array,1) * size(this.array,4)
    end function message_length_4d_r8

    function message_length_4d_r4(this) result(mlength)
        class(grid_function_4d_r4_t), intent(in) :: this
        integer(kind=4) :: mlength

        mlength = size(this.array,1) * size(this.array,4)
    end function message_length_4d_r4
    
    function message_length_3d_r8_depth(this) result(mlength)
        class(grid_function_3d_r8_depth_t), intent(in) :: this
        integer(kind=4) :: mlength

        mlength = size(this.array,3)
    end function message_length_3d_r8_depth

    function message_length_3d_r4_depth(this) result(mlength)
        class(grid_function_3d_r4_depth_t), intent(in) :: this
        integer(kind=4) :: mlength

        mlength = size(this.array,3)
    end function message_length_3d_r4_depth

    function message_length_4d_r8_depth(this) result(mlength)
        class(grid_function_4d_r8_depth_t), intent(in) :: this
        integer(kind=4) :: mlength

        mlength = size(this.array,1) * size(this.array,4)
    end function message_length_4d_r8_depth

    function message_length_4d_r4_depth(this) result(mlength)
        class(grid_function_4d_r4_depth_t), intent(in) :: this
        integer(kind=4) :: mlength

        mlength = size(this.array,1) * size(this.array,4)
    end function message_length_4d_r4_depth
    
    subroutine get_sub_array_2d_r8(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_2d_r8_t), intent(in)    :: this
        integer(kind=4)             , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                 , intent(in)    :: buf
        integer(kind=4)             , intent(inout) :: idx

        integer(kind=4) :: i, j
        real(kind=8), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])

        do j = js, je
            do i = is, ie
                idx = idx + 1
                buffer(idx) = this.array(i,j)
            end do
        end do
    end subroutine get_sub_array_2d_r8

    subroutine get_sub_array_2d_r4(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_2d_r4_t), intent(in)    :: this
        integer(kind=4)             , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                 , intent(in)    :: buf
        integer(kind=4)             , intent(inout) :: idx

        integer(kind=4) :: i, j
        real(kind=4), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])

        do j = js, je
            do i = is, ie
                idx = idx + 1
                buffer(idx) = this.array(i,j)
            end do
        end do
    end subroutine get_sub_array_2d_r4

    subroutine get_sub_array_3d_r8(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_3d_r8_t), intent(in)    :: this
        integer(kind=4)             , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                 , intent(in)    :: buf
        integer(kind=4)             , intent(inout) :: idx

        integer(kind=4) :: i, j, k, ks, ke
        real(kind=8), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])
        ks = 1
        ke = size(this.array,3)

        do k = ks, ke
            do j = js, je
                do i = is, ie
                    idx = idx + 1
                    buffer(idx) = this.array(i,j,k)
                end do
            end do
        end do
    end subroutine get_sub_array_3d_r8

    subroutine get_sub_array_3d_r4(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_3d_r4_t), intent(in)    :: this
        integer(kind=4)             , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                 , intent(in)    :: buf
        integer(kind=4)             , intent(inout) :: idx

        integer(kind=4) :: i, j, k, ks, ke
        real(kind=4), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])
        ks = 1
        ke = size(this.array,3)

        do k = ks, ke
            do j = js, je
                do i = is, ie
                    idx = idx + 1
                    buffer(idx) = this.array(i,j,k)
                end do
            end do
        end do
    end subroutine get_sub_array_3d_r4

    subroutine get_sub_array_4d_r8(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_4d_r8_t), intent(in)    :: this
        integer(kind=4)             , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                 , intent(in)    :: buf
        integer(kind=4)             , intent(inout) :: idx

        integer(kind=4) :: l, i, j, k, ls, le, ks, ke
        real(kind=8), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])
        ls = 1
        le = size(this.array,1)
        ks = 1
        ke = size(this.array,4)

        do k = ks, ke
            do j = js, je
                do i = is, ie
                    do l = ls, le
                        idx = idx + 1
                        buffer(idx) = this.array(l,i,j,k)
                    end do
                end do
            end do
        end do
    end subroutine get_sub_array_4d_r8

    subroutine get_sub_array_4d_r4(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_4d_r4_t), intent(in)    :: this
        integer(kind=4)             , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                 , intent(in)    :: buf
        integer(kind=4)             , intent(inout) :: idx

        integer(kind=4) :: l, i, j, k, ls, le, ks, ke
        real(kind=4), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])
        ls = 1
        le = size(this.array,1)
        ks = 1
        ke = size(this.array,4)

        do k = ks, ke
            do j = js, je
                do i = is, ie
                    do l = ls, le
                        idx = idx + 1
                        buffer(idx) = this.array(l,i,j,k)
                    end do
                end do
            end do
        end do
    end subroutine get_sub_array_4d_r4
    
    subroutine get_sub_array_3d_r8_depth(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_3d_r8_depth_t), intent(in)    :: this
        integer(kind=4)                   , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                       , intent(in)    :: buf
        integer(kind=4)                   , intent(inout) :: idx

        integer(kind=4) :: i, j, k
        real(kind=8), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])

        do j = js, je
            do i = is, ie
                do k = 1, this.km2_l(i,j)
                    idx = idx + 1
                    buffer(idx) = this.array(i,j,k)
                end do
            end do
        end do
    end subroutine get_sub_array_3d_r8_depth

    subroutine get_sub_array_3d_r4_depth(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_3d_r4_depth_t), intent(in)    :: this
        integer(kind=4)                   , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                       , intent(in)    :: buf
        integer(kind=4)                   , intent(inout) :: idx

        integer(kind=4) :: i, j, k
        real(kind=4), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])

        do j = js, je
            do i = is, ie
                do k = 1, this.km2_l(i,j)
                    idx = idx + 1
                    buffer(idx) = this.array(i,j,k)
                end do
            end do
        end do
    end subroutine get_sub_array_3d_r4_depth
    
    subroutine get_sub_array_4d_r8_depth(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_4d_r8_depth_t), intent(in)    :: this
        integer(kind=4)                   , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                       , intent(in)    :: buf
        integer(kind=4)                   , intent(inout) :: idx

        integer(kind=4) :: l, i, j, k, ls, le
        real(kind=8), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])
        
        ls = 1
        le = size(this.array,1)

        do j = js, je
            do i = is, ie
                do k = 1, this.km2_l(i,j)
                    do l = ls, le
                        idx = idx + 1
                        buffer(idx) = this.array(l,i,j,k)
                    end do
                end do
            end do
        end do
    end subroutine get_sub_array_4d_r8_depth

    subroutine get_sub_array_4d_r4_depth(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_4d_r4_depth_t), intent(in)    :: this
        integer(kind=4)                   , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                       , intent(in)    :: buf
        integer(kind=4)                   , intent(inout) :: idx

        integer(kind=4) :: l, i, j, k, ls, le
        real(kind=4), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])
        
        ls = 1
        le = size(this.array,1)

        do j = js, je
            do i = is, ie
                do k = 1, this.km2_l(i,j)
                    do l = ls, le
                        idx = idx + 1
                        buffer(idx) = this.array(l,i,j,k)
                    end do
                end do
            end do
        end do
    end subroutine get_sub_array_4d_r4_depth

    subroutine put_sub_array_2d_r8(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_2d_r8_t), intent(inout) :: this
        integer(kind=4)             , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                 , intent(in)    :: buf
        integer(kind=4)             , intent(inout) :: idx

        integer(kind=4) :: i, j
        real(kind=8), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])

        do j = js, je
            do i = is, ie
                idx = idx + 1
                this.array(i,j) = buffer(idx)
            end do
        end do
    end subroutine put_sub_array_2d_r8

    subroutine put_sub_array_2d_r4(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_2d_r4_t), intent(inout) :: this
        integer(kind=4)             , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                 , intent(in)    :: buf
        integer(kind=4)             , intent(inout) :: idx

        integer(kind=4) :: i, j
        real(kind=4), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])

        do j = js, je
            do i = is, ie
                idx = idx + 1
                this.array(i,j) = buffer(idx)
            end do
        end do
    end subroutine put_sub_array_2d_r4

    subroutine put_sub_array_3d_r8(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_3d_r8_t), intent(inout) :: this
        integer(kind=4)             , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                 , intent(in)    :: buf
        integer(kind=4)             , intent(inout) :: idx

        integer(kind=4) :: i, j, k, ks, ke
        real(kind=8), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])
        ks = 1
        ke = size(this.array,3)

        do k = ks, ke
            do j = js, je
                do i = is, ie
                    idx = idx + 1
                    this.array(i,j,k) = buffer(idx)
                end do
            end do
        end do
    end subroutine put_sub_array_3d_r8

    subroutine put_sub_array_3d_r4(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_3d_r4_t), intent(inout) :: this
        integer(kind=4)             , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                 , intent(in)    :: buf
        integer(kind=4)             , intent(inout) :: idx

        integer(kind=4) :: i, j, k, ks, ke
        real(kind=4), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])
        ks = 1
        ke = size(this.array,3)

        do k = ks, ke
            do j = js, je
                do i = is, ie
                    idx = idx + 1
                    this.array(i,j,k) = buffer(idx)
                end do
            end do
        end do
    end subroutine put_sub_array_3d_r4

    subroutine put_sub_array_4d_r8(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_4d_r8_t), intent(inout) :: this
        integer(kind=4)             , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                 , intent(in)    :: buf
        integer(kind=4)             , intent(inout) :: idx

        integer(kind=4) :: l, i, j, k, ls, le, ks, ke
        real(kind=8), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])
        ls = 1
        le = size(this.array,1)
        ks = 1
        ke = size(this.array,4)

        do k = ks, ke
            do j = js, je
                do i = is, ie
                    do l = ls, le
                        idx = idx + 1
                        this.array(l,i,j,k) = buffer(idx)
                    end do
                end do
            end do
        end do
    end subroutine put_sub_array_4d_r8

    subroutine put_sub_array_4d_r4(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_4d_r4_t), intent(inout) :: this
        integer(kind=4)             , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                 , intent(in)    :: buf
        integer(kind=4)             , intent(inout) :: idx

        integer(kind=4) :: l, i, j, k, ls, le, ks, ke
        real(kind=4), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])
        ls = 1
        le = size(this.array,1)
        ks = 1
        ke = size(this.array,4)

        do k = ks, ke
            do j = js, je
                do i = is, ie
                    do l = ls, le
                        idx = idx + 1
                        this.array(l,i,j,k) = buffer(idx)
                    end do
                end do
            end do
        end do
    end subroutine put_sub_array_4d_r4
    
    subroutine put_sub_array_3d_r8_depth(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_3d_r8_depth_t), intent(inout) :: this
        integer(kind=4)                   , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                       , intent(in)    :: buf
        integer(kind=4)                   , intent(inout) :: idx

        integer(kind=4) :: i, j, k
        real(kind=8), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])

        do j = js, je
            do i = is, ie
                do k = 1, this.km2_l(i,j)
                    idx = idx + 1
                    this.array(i,j,k) = buffer(idx)
                end do
            end do
        end do
    end subroutine put_sub_array_3d_r8_depth

    subroutine put_sub_array_3d_r4_depth(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_3d_r4_depth_t), intent(inout) :: this
        integer(kind=4)                   , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                       , intent(in)    :: buf
        integer(kind=4)                   , intent(inout) :: idx

        integer(kind=4) :: i, j, k
        real(kind=4), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])
        
        do j = js, je
            do i = is, ie
                do k = 1, this.km2_l(i,j)
                    idx = idx + 1
                    this.array(i,j,k) = buffer(idx)
                end do
            end do
        end do
    end subroutine put_sub_array_3d_r4_depth
    
    subroutine put_sub_array_4d_r8_depth(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_4d_r8_depth_t), intent(inout) :: this
        integer(kind=4)                   , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                       , intent(in)    :: buf
        integer(kind=4)                   , intent(inout) :: idx

        integer(kind=4) :: l, i, j, k, ls, le
        real(kind=8), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])

        ls = 1
        le = size(this.array,1)

        do j = js, je
            do i = is, ie
                do k = 1, this.km2_l(i,j)
                    do l = ls, le
                        idx = idx + 1
                        this.array(l,i,j,k) = buffer(idx)
                    end do
                end do
            end do
        end do
    end subroutine put_sub_array_4d_r8_depth

    subroutine put_sub_array_4d_r4_depth(this, buf, idx, idx_max, is, ie, js, je)
        class(grid_function_4d_r4_depth_t), intent(inout) :: this
        integer(kind=4)                   , intent(in)    :: idx_max, is, ie, js, je
        type(c_ptr)                       , intent(in)    :: buf
        integer(kind=4)                   , intent(inout) :: idx

        integer(kind=4) :: l, i, j, k, ls, le
        real(kind=4), pointer :: buffer(:)

        call c_f_pointer(buf, buffer, [idx_max])
        
        ls = 1
        le = size(this.array,1)

        do j = js, je
            do i = is, ie
                do k = 1, this.km2_l(i,j)
                    do l = ls, le
                        idx = idx + 1
                        this.array(l,i,j,k) = buffer(idx)
                    end do
                end do
            end do
        end do
    end subroutine put_sub_array_4d_r4_depth

end module grid_function_mod
