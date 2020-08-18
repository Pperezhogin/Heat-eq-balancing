module halo_mod
    implicit none

    type, public :: halo_x_t
        
        ! second iterator
        integer(kind=4) :: num_ranks
        integer(kind=4), allocatable :: ranks(:)
        ! first iterator
        integer(kind=4), allocatable :: nstrips(:)

        integer(kind=4), allocatable :: i(:,:)
        integer(kind=4), allocatable :: js(:,:), je(:,:)
        logical(kind=1), allocatable :: remote_on_right(:,:)

        integer(kind=4), allocatable :: send_req(:)  , recv_req(:) 
        integer(kind=8), allocatable :: send_buf(:,:), recv_buf(:,:)

        integer(kind=4) :: max_message_size

        contains

        procedure, public :: init         => init_halo_x
        procedure, public :: cleanup      => cleanup_halo_x
        procedure, public :: allocate_buf => allocate_buf_halo_x        

    end type halo_x_t

    type, public :: halo_y_t
        
        ! second iterator
        integer(kind=4) :: num_ranks
        integer(kind=4), allocatable :: ranks(:)
        ! first iterator
        integer(kind=4), allocatable :: nstrips(:)

        integer(kind=4), allocatable :: j(:,:)
        integer(kind=4), allocatable :: is(:,:), ie(:,:)

        logical(kind=1), allocatable :: remote_on_top(:,:)

        logical(kind=1), allocatable :: right_border_send(:,:), left_border_send(:,:)
        logical(kind=1), allocatable :: right_border_recv(:,:), left_border_recv(:,:)

        integer(kind=4), allocatable :: send_req(:)  , recv_req(:) 
        integer(kind=8), allocatable :: send_buf(:,:), recv_buf(:,:)

        integer(kind=4) :: max_message_size

        contains

        procedure, public :: init         => init_halo_y
        procedure, public :: cleanup      => cleanup_halo_y
        procedure, public :: allocate_buf => allocate_buf_halo_y

    end type halo_y_t

    type, public :: halo_diag_t
        
        ! second iterator
        integer(kind=4) :: nranks_send, nranks_recv
        integer(kind=4), allocatable :: ranks_send(:), ranks_recv(:)

        ! first iterator
        integer(kind=4), allocatable :: nsend(:), nrecv(:)
        
        integer(kind=4), allocatable :: i_send(:,:), i_recv(:,:)
        integer(kind=4), allocatable :: j_send(:,:), j_recv(:,:)

        logical(kind=1), allocatable :: right_send(:,:), right_recv(:,:)
        logical(kind=1), allocatable :: top_send  (:,:), top_recv  (:,:)

        integer(kind=4), allocatable :: send_req(:)  , recv_req(:) 
        integer(kind=8), allocatable :: send_buf(:,:), recv_buf(:,:)

        integer(kind=4) :: max_message_size

        contains

        procedure, public :: init         => init_halo_diag
        procedure, public :: cleanup      => cleanup_halo_diag
        procedure, public :: allocate_buf => allocate_buf_halo_diag

    end type halo_diag_t

    private

    contains

    subroutine init_halo_x(this, num_ranks, ranks, nstrips, &
        i, js, je, remote_on_right, max_message_size)
        class(halo_x_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: num_ranks, max_message_size
        
        integer(kind=4), dimension(:  ), intent(in) :: ranks, nstrips
        integer(kind=4), dimension(:,:), intent(in) :: i, js, je
        logical(kind=1), dimension(:,:), intent(in) :: remote_on_right

        integer(kind=4) :: max_nstrips
        max_nstrips = maxval(nstrips(1:num_ranks))

        allocate(this.ranks  (num_ranks))
        allocate(this.nstrips(num_ranks))

        allocate(this.i (max_nstrips, num_ranks))
        allocate(this.js(max_nstrips, num_ranks))
        allocate(this.je(max_nstrips, num_ranks))

        allocate(this.remote_on_right(max_nstrips, num_ranks))

        allocate(this.send_req(num_ranks))
        allocate(this.recv_req(num_ranks))
        allocate(this.send_buf(1,1))
        allocate(this.recv_buf(1,1))

        this.num_ranks        = num_ranks
        this.max_message_size = max_message_size

        this.ranks   = ranks  (1:num_ranks)
        this.nstrips = nstrips(1:num_ranks)

        this.i  = i (1:max_nstrips,1:num_ranks)
        this.js = js(1:max_nstrips,1:num_ranks)
        this.je = je(1:max_nstrips,1:num_ranks)

        this.remote_on_right = remote_on_right(1:max_nstrips,1:num_ranks)

    end subroutine init_halo_x

    subroutine allocate_buf_halo_x(this, message_size)
        class(halo_x_t), intent(inout) :: this
        integer(kind=4), intent(in)    :: message_size

        if (size(this.send_buf,1) < message_size) then
            if (allocated(this.send_buf)) deallocate(this.send_buf)
            allocate(this.send_buf(message_size,this.num_ranks))

            if (allocated(this.recv_buf)) deallocate(this.recv_buf)
            allocate(this.recv_buf(message_size,this.num_ranks))
        end if

    end subroutine allocate_buf_halo_x

    subroutine cleanup_halo_x(this)
        class(halo_x_t), intent(inout) :: this

        if (allocated(this.ranks  )) deallocate(this.ranks)
        if (allocated(this.nstrips)) deallocate(this.nstrips)

        if (allocated(this.i )) deallocate(this.i )
        if (allocated(this.js)) deallocate(this.js)
        if (allocated(this.je)) deallocate(this.je)

        if (allocated(this.remote_on_right)) deallocate(this.remote_on_right)

        if (allocated(this.send_req)) deallocate(this.send_req)
        if (allocated(this.recv_req)) deallocate(this.recv_req)

        if (allocated(this.send_buf)) deallocate(this.send_buf)
        if (allocated(this.recv_buf)) deallocate(this.recv_buf)

    end subroutine cleanup_halo_x

    subroutine init_halo_y(this, num_ranks, ranks, nstrips, &
        j, is, ie, remote_on_top, right_border_send, left_border_send, &
        right_border_recv, left_border_recv, max_message_size)
        class(halo_y_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: num_ranks, max_message_size
        
        integer(kind=4), dimension(:  ), intent(in) :: ranks, nstrips
        integer(kind=4), dimension(:,:), intent(in) :: j, is, ie
        logical(kind=1), dimension(:,:), intent(in) :: remote_on_top
        logical(kind=1), dimension(:,:), intent(in) :: right_border_send, left_border_send
        logical(kind=1), dimension(:,:), intent(in) :: right_border_recv, left_border_recv

        integer(kind=4) :: max_nstrips
        max_nstrips = maxval(nstrips(1:num_ranks))

        allocate(this.ranks  (num_ranks))
        allocate(this.nstrips(num_ranks))

        allocate(this.j (max_nstrips, num_ranks))
        allocate(this.is(max_nstrips, num_ranks))
        allocate(this.ie(max_nstrips, num_ranks))

        allocate(this.remote_on_top    (max_nstrips, num_ranks))
        allocate(this.right_border_send(max_nstrips, num_ranks))
        allocate(this.right_border_recv(max_nstrips, num_ranks))
        allocate(this.left_border_send (max_nstrips, num_ranks))
        allocate(this.left_border_recv (max_nstrips, num_ranks))

        allocate(this.send_req(num_ranks))
        allocate(this.recv_req(num_ranks))
        allocate(this.send_buf(1,1))
        allocate(this.recv_buf(1,1))

        this.num_ranks        = num_ranks
        this.max_message_size = max_message_size

        this.ranks   = ranks  (1:num_ranks)
        this.nstrips = nstrips(1:num_ranks)

        this.j  = j (1:max_nstrips,1:num_ranks)
        this.is = is(1:max_nstrips,1:num_ranks)
        this.ie = ie(1:max_nstrips,1:num_ranks)

        this.remote_on_top     = remote_on_top    (1:max_nstrips,1:num_ranks)
        this.right_border_send = right_border_send(1:max_nstrips,1:num_ranks)
        this.right_border_recv = right_border_recv(1:max_nstrips,1:num_ranks)
        this.left_border_send  = left_border_send (1:max_nstrips,1:num_ranks)
        this.left_border_recv  = left_border_recv (1:max_nstrips,1:num_ranks)

    end subroutine init_halo_y

    subroutine allocate_buf_halo_y(this, message_size)
        class(halo_y_t), intent(inout) :: this
        integer(kind=4), intent(in)    :: message_size

        if (size(this.send_buf,1) < message_size) then
            if (allocated(this.send_buf)) deallocate(this.send_buf)
            allocate(this.send_buf(message_size,this.num_ranks))

            if (allocated(this.recv_buf)) deallocate(this.recv_buf)
            allocate(this.recv_buf(message_size,this.num_ranks))
        end if

    end subroutine allocate_buf_halo_y

    subroutine cleanup_halo_y(this)
        class(halo_y_t), intent(inout) :: this

        if (allocated(this.ranks  )) deallocate(this.ranks)
        if (allocated(this.nstrips)) deallocate(this.nstrips)

        if (allocated(this.j )) deallocate(this.j )
        if (allocated(this.is)) deallocate(this.is)
        if (allocated(this.ie)) deallocate(this.ie)

        if (allocated(this.remote_on_top    )) deallocate(this.remote_on_top    )
        if (allocated(this.right_border_send)) deallocate(this.right_border_send)
        if (allocated(this.right_border_recv)) deallocate(this.right_border_recv)
        if (allocated(this.left_border_send )) deallocate(this.left_border_send )
        if (allocated(this.left_border_recv )) deallocate(this.left_border_recv )

        if (allocated(this.send_req)) deallocate(this.send_req)
        if (allocated(this.recv_req)) deallocate(this.recv_req)

        if (allocated(this.send_buf)) deallocate(this.send_buf)
        if (allocated(this.recv_buf)) deallocate(this.recv_buf)

    end subroutine cleanup_halo_y

    subroutine init_halo_diag(this, nranks_send, nranks_recv, ranks_send, ranks_recv, &
        nsend, nrecv, i_send, i_recv, j_send, j_recv, right_send, right_recv, &
        top_send, top_recv, halo_size)

        class(halo_diag_t), intent(inout) :: this
        integer(kind=4)   , intent(in   ) :: nranks_send, nranks_recv, halo_size
        
        integer(kind=4), dimension(:  ), intent(in) :: ranks_send, nsend
        integer(kind=4), dimension(:  ), intent(in) :: ranks_recv, nrecv
        integer(kind=4), dimension(:,:), intent(in) :: i_send, i_recv
        integer(kind=4), dimension(:,:), intent(in) :: j_send, j_recv
        logical(kind=1), dimension(:,:), intent(in) :: right_send, top_send
        logical(kind=1), dimension(:,:), intent(in) :: right_recv, top_recv

        integer(kind=4) :: max_nsend, max_nrecv
        max_nsend = maxval(nsend(1:nranks_send))
        max_nrecv = maxval(nrecv(1:nranks_recv))
        this.max_message_size = halo_size * max(max_nsend, max_nrecv)

        allocate(this.ranks_send(nranks_send))
        allocate(this.ranks_recv(nranks_recv))

        allocate(this.nsend(nranks_send))
        allocate(this.nrecv(nranks_recv))

        allocate(this.i_send(max_nsend, nranks_send))
        allocate(this.j_send(max_nsend, nranks_send))
        allocate(this.i_recv(max_nrecv, nranks_recv))
        allocate(this.j_recv(max_nrecv, nranks_recv))

        allocate(this.right_send(max_nsend, nranks_send))
        allocate(this.right_recv(max_nrecv, nranks_recv))

        allocate(this.top_send  (max_nsend, nranks_send))
        allocate(this.top_recv  (max_nrecv, nranks_recv))

        allocate(this.send_req(nranks_send))
        allocate(this.recv_req(nranks_recv))
        allocate(this.send_buf(1,1))
        allocate(this.recv_buf(1,1))

        this.nranks_send = nranks_send
        this.nranks_recv = nranks_recv

        this.ranks_send = ranks_send(1:nranks_send)
        this.ranks_recv = ranks_recv(1:nranks_recv)

        this.nsend = nsend(1:nranks_send)
        this.nrecv = nrecv(1:nranks_recv)

        this.i_send = i_send(1:max_nsend,1:nranks_send)
        this.i_recv = i_recv(1:max_nrecv,1:nranks_recv)
        this.j_send = j_send(1:max_nsend,1:nranks_send)
        this.j_recv = j_recv(1:max_nrecv,1:nranks_recv)

        this.right_send = right_send(1:max_nsend,1:nranks_send)
        this.right_recv = right_recv(1:max_nsend,1:nranks_send)
        this.top_send   = top_send  (1:max_nsend,1:nranks_send)
        this.top_recv   = top_recv  (1:max_nsend,1:nranks_send)

    end subroutine init_halo_diag

    subroutine allocate_buf_halo_diag(this, message_size)
        class(halo_diag_t), intent(inout) :: this
        integer(kind=4), intent(in)    :: message_size

        if (size(this.send_buf,1) < message_size) then
            if (allocated(this.send_buf)) deallocate(this.send_buf)
            allocate(this.send_buf(message_size,this.nranks_send))

            if (allocated(this.recv_buf)) deallocate(this.recv_buf)
            allocate(this.recv_buf(message_size,this.nranks_recv))
        end if
    end subroutine allocate_buf_halo_diag

    subroutine cleanup_halo_diag(this)
        class(halo_diag_t), intent(inout) :: this
       
        if (allocated(this.ranks_send)) deallocate(this.ranks_send)
        if (allocated(this.ranks_recv)) deallocate(this.ranks_recv)

        if (allocated(this.nsend)) deallocate(this.nsend)
        if (allocated(this.nrecv)) deallocate(this.nrecv)

        if (allocated(this.i_send)) deallocate(this.i_send)
        if (allocated(this.j_send)) deallocate(this.j_send)
        if (allocated(this.i_recv)) deallocate(this.i_recv)
        if (allocated(this.j_recv)) deallocate(this.j_recv)

        if (allocated(this.right_send)) deallocate(this.right_send)
        if (allocated(this.right_recv)) deallocate(this.right_recv)

        if (allocated(this.top_send)) deallocate(this.top_send)
        if (allocated(this.top_recv)) deallocate(this.top_recv)

        if (allocated(this.send_req)) deallocate(this.send_req)
        if (allocated(this.recv_req)) deallocate(this.recv_req)

        if (allocated(this.send_buf)) deallocate(this.send_buf)
        if (allocated(this.recv_buf)) deallocate(this.recv_buf)

    end subroutine cleanup_halo_diag

end module halo_mod