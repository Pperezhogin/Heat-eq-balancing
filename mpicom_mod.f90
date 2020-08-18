module mpicom_mod
    use mpi
    use halo_mod
    use partition_mod
    use data_mod
    use grid_function_mod
    implicit none

    type, public :: mpicom_t
        integer(kind=4), public :: com, size, rank, ierr
        logical(kind=1), public :: use_mpi_barriers ! introduce barriers, which slows down model, but allows to measure exchange time

        private

        integer(kind=4) :: is  , ie  , js  , je
        integer(kind=4) :: is_g, ie_g, js_g, je_g

        integer(kind=4) :: max_gcx, max_gcy

        type(halo_x_t)    :: halo_x
        type(halo_y_t)    :: halo_y
        type(halo_diag_t) :: halo_diag

        ! gather-scatter data
        integer(kind=4), allocatable :: nblocks(:)
        integer(kind=4), allocatable :: is_gc0(:,:), ie_gc0(:,:), js_gc0(:,:), je_gc0(:,:)
        integer(kind=8), allocatable :: buf_g(:), buf_l(:)
        integer(kind=4) :: scatter_gcx

        type(data_t) :: data, glob_data

        real(kind=8) :: time_exch
        real(kind=8) :: time_gs, time_gs_copy, time_gs_send

        contains

        procedure, public  :: init                  => init_mpicom
        procedure, public  :: init_exchange         => init_exchange_mpicom
        procedure, public  :: cleanup               => cleanup_mpicom
        procedure, public  :: time_exchange         => time_exch_mpicom
        procedure, public  :: time_scatter          => time_scatter_mpicom
        procedure, public  :: print_diagnostics     => print_diagnostics_mpicom

        generic  , public  :: exchange_cross        => exchange_cross_2d_r8, exchange_cross_2d_r4, &
                                                       exchange_cross_3d_r8, exchange_cross_3d_r4, &
                                                       exchange_cross_4d_r8, exchange_cross_4d_r4
        generic  , public  :: exchange_halo         => exchange_halo_2d_r8,  exchange_halo_2d_r4,  &
                                                       exchange_halo_3d_r8,  exchange_halo_3d_r4,  &
                                                       exchange_halo_4d_r8,  exchange_halo_4d_r4,  &
                                                       exchange_halo_3d_r8_depth,  exchange_halo_3d_r4_depth, & ! array and km2_l should have the same gcx
                                                       exchange_halo_4d_r8_depth,  exchange_halo_4d_r4_depth    ! array and km2_l should have the same gcx
        generic  , public  :: scatter               => scatter_2d_r8,  scatter_2d_r4,  & ! WITH boundary = gcx
                                                       scatter_3d_r8,  scatter_3d_r4,  &
                                                       scatter_4d_r8,  scatter_4d_r4
        generic  , public  :: gather                => gather_2d_r8,  gather_2d_r4,  &   ! without boundary
                                                       gather_3d_r8,  gather_3d_r4,  &
                                                       gather_4d_r8,  gather_4d_r4

        procedure, private :: find_halo_x
        procedure, private :: find_halo_y
        procedure, private :: find_halo_diag
        procedure, private :: exchange_x
        procedure, private :: exchange_y
        procedure, private :: exchange_y_halo
        procedure, private :: scatter_mpicom
        procedure, private :: gather_mpicom
        procedure, private :: check_cross
        procedure, private :: check_halo
        procedure, private :: check_scatter

        procedure, private :: exchange_cross_2d_r8
        procedure, private :: exchange_cross_2d_r4
        procedure, private :: exchange_cross_3d_r8
        procedure, private :: exchange_cross_3d_r4
        procedure, private :: exchange_cross_4d_r8
        procedure, private :: exchange_cross_4d_r4

        procedure, private :: exchange_halo_2d_r8
        procedure, private :: exchange_halo_2d_r4
        procedure, private :: exchange_halo_3d_r8
        procedure, private :: exchange_halo_3d_r4
        procedure, private :: exchange_halo_4d_r8
        procedure, private :: exchange_halo_4d_r4
        procedure, private :: exchange_halo_3d_r8_depth
        procedure, private :: exchange_halo_3d_r4_depth
        procedure, private :: exchange_halo_4d_r8_depth
        procedure, private :: exchange_halo_4d_r4_depth

        procedure, private :: scatter_2d_r8
        procedure, private :: scatter_2d_r4
        procedure, private :: scatter_3d_r8
        procedure, private :: scatter_3d_r4
        procedure, private :: scatter_4d_r8
        procedure, private :: scatter_4d_r4

        procedure, private :: gather_2d_r8
        procedure, private :: gather_2d_r4
        procedure, private :: gather_3d_r8
        procedure, private :: gather_3d_r4
        procedure, private :: gather_4d_r8
        procedure, private :: gather_4d_r4

    end type mpicom_t

    private

    contains

    subroutine init_mpicom(this, mpi_com)
        class(mpicom_t), intent(inout) :: this 
        integer(kind=4), intent(in)    :: mpi_com

        this.com = mpi_com
        call mpi_comm_size(this.com, this.size, this.ierr)
        call mpi_comm_rank(this.com, this.rank, this.ierr)
        this.use_mpi_barriers = .true.

    end subroutine init_mpicom

    subroutine cleanup_mpicom(this)
        class(mpicom_t), intent(inout) :: this 

        call this.halo_x.cleanup()
        call this.halo_y.cleanup()
        call this.halo_diag.cleanup()

        if(allocated(this.nblocks)) deallocate(this.nblocks)
        if(allocated(this.buf_g))   deallocate(this.buf_g)
        if(allocated(this.buf_l))   deallocate(this.buf_l)
        if(allocated(this.is_gc0))    deallocate(this.is_gc0)
        if(allocated(this.ie_gc0))    deallocate(this.ie_gc0)
        if(allocated(this.js_gc0))    deallocate(this.js_gc0)
        if(allocated(this.je_gc0))    deallocate(this.je_gc0) 
        
        this.data.num_grid_functions = 0
        this.glob_data.num_grid_functions = 0

        this.time_exch = 0
        this.time_gs = 0
        this.time_gs_copy = 0
        this.time_gs_send = 0

    end subroutine cleanup_mpicom

    subroutine init_exchange_mpicom(this, partition, max_gcx, max_gcy)
        class(mpicom_t)   , intent(inout) :: this 
        type (partition_t), intent(in)    :: partition
        integer(kind=4)   , intent(in)    :: max_gcx, max_gcy

        integer(kind=4) :: max_nblocks

        call this.cleanup()

        this.max_gcx = min(partition.get_gcx(), max_gcx)
        this.max_gcy = min(partition.get_gcy(), max_gcy)

        call partition.global_idx(this.is_g,this.ie_g,this.js_g,this.je_g)
        call partition.local_idx (this.is,this.ie,this.js,this.je,this.rank)

        allocate(this.nblocks(0:this.size-1))
        call partition.get_blocks_number_on_rank(this.nblocks)
        max_nblocks = maxval(this.nblocks)

        allocate(this.is_gc0(max_nblocks,0:this.size-1),this.ie_gc0(max_nblocks,0:this.size-1))
        allocate(this.js_gc0(max_nblocks,0:this.size-1),this.je_gc0(max_nblocks,0:this.size-1))
        call partition.get_blocks_on_rank(this.is_gc0,this.ie_gc0,this.js_gc0,this.je_gc0,max_nblocks)

        call this.find_halo_x   (partition)
        call this.find_halo_y   (partition)
        call this.find_halo_diag(partition)

        this.data.num_grid_functions = 0
        this.glob_data.num_grid_functions = 0

        allocate(this.buf_g(1)) ! scatter_mpicom works only with allocated buffers
        allocate(this.buf_l(1))
        
        call this.check_cross  (partition)
        call this.check_halo   (partition)
        call this.check_scatter(partition)
        
        if(allocated(this.buf_g))   deallocate(this.buf_g)
        if(allocated(this.buf_l))   deallocate(this.buf_l)
        allocate(this.buf_g(1))
        allocate(this.buf_l(1))
        
        this.time_exch = 0
        this.time_gs = 0
        this.time_gs_copy = 0
        this.time_gs_send = 0

    end subroutine init_exchange_mpicom
    
    subroutine print_diagnostics_mpicom(this, partition, time_adv, time_full, message)
        class(mpicom_t), intent(in) :: this 
        type (partition_t), intent(in) :: partition
        real(kind=8), intent(in) :: time_adv, time_full
        character(*), intent(in) :: message
        
        integer(kind=4) :: i,j
        real   (kind=8), dimension(0:this.size-1) :: t_adv, t_full
        real   (kind=8), dimension(0:this.size-1) :: t_exch
        real   (kind=8), dimension(0:this.size-1) :: t_gs, t_gs_copy, t_gs_send
        integer(kind=4) :: rankij(this.is_g:this.ie_g, this.js_g:this.je_g)
        integer(kind=4) :: rank
        character(20)   :: numproc
        
        call partition.get_rankij(rankij,0,0)

        write(numproc,*) this.size
        numproc = adjustl(numproc)

        call MPI_Gather(time_adv   ,1,MPI_REAL8,t_adv  ,1,MPI_REAL8,0,this.com,this.ierr)
        call MPI_Gather(time_full  ,1,MPI_REAL8,t_full ,1,MPI_REAL8,0,this.com,this.ierr)
        
        call MPI_Gather(this.time_exch   ,1,MPI_REAL8,t_exch   ,1,MPI_REAL8,0,this.com,this.ierr)
        
        if (this.rank == 0) then
            open(20, file = 'diagnostics.out')
            
            write(20,*) message
            
            write(20,*) 'Full time', real(maxval(t_full))
            write(20,*) 'MPI exchange time', real(maxval(t_exch))
            
            write(20,*) ''
            
            write(20,*) 'Main loop time', real(maxval(t_adv))
            write(20,*) "Main loop load imbalance in percent:", & 
                (maxval(t_adv) * this.size - sum(t_adv)) / sum(t_adv) * 100.0_8
                
            write(20,*) ''
            
            write(20,*) ''

            do rank = 0, this.size-1
                write(20,'("rank ",I4," time_exch ",F8.4, &
                " main loop ",F8.4)') &
                rank, t_exch(rank), t_adv(rank)
            end do

            close(20)
        end if
    end subroutine print_diagnostics_mpicom

    function time_exch_mpicom(this) result(time)
        class(mpicom_t)   , intent(inout) :: this 
        real(kind=8) :: time

        call MPI_Allreduce(this.time_exch,time,1,MPI_DOUBLE_PRECISION,MPI_MAX,this.com,this.ierr)
    end function time_exch_mpicom

    function time_scatter_mpicom(this) result(time)
        class(mpicom_t)   , intent(inout) :: this 
        real(kind=8) :: time

        call MPI_Allreduce(this.time_gs,time,1,MPI_DOUBLE_PRECISION,MPI_MAX,this.com,this.ierr)
    end function time_scatter_mpicom

    subroutine find_halo_x(this, partition)
        class(mpicom_t)  , intent(inout) :: this
        type(partition_t), intent(in)    :: partition

        integer(kind=4), dimension(this.is_g-1:this.ie_g+1,this.js_g-1:this.je_g+1) :: rankij

        integer(kind=4), dimension(     this.size) :: ranks, nstrips
        integer(kind=4), dimension(1000,this.size) :: irb, js_exch, je_exch
        logical(kind=1), dimension(1000,this.size) :: remote_on_right

        integer(kind=4) :: i, j, irank, istrip
        integer(kind=4) :: rank_remote, rank_local
        integer(kind=4) :: message_size, max_message_size

        logical(kind=1) :: cond1, cond2
        logical(kind=1) :: remote_right  , remote_left
        logical(kind=1) :: on_strip_right, on_strip_left
        logical(kind=1) :: active_rank

        call partition.get_rankij(rankij,1,1)

        irank = 0
        rank_local = this.rank
        max_message_size = 0
        do rank_remote = 0, this.size-1
            if (rank_remote == rank_local) cycle
            
            istrip = 0
            message_size = 0
            active_rank = .false.
            
            do i = this.is-1, this.ie
                on_strip_left  = .false.
                on_strip_right = .false.
                do j = this.js, this.je
                    cond1 = (rankij(i  ,j) == rank_local )
                    cond2 = (rankij(i+1,j) == rank_remote)
                    remote_right = cond1 .and. cond2

                    cond1 = (rankij(i+1,j) == rank_local)
                    cond2 = (rankij(i  ,j) == rank_remote)
                    remote_left  = cond1 .and. cond2

                    if ((remote_right .or. remote_left) .and. not(active_rank)) then
                        active_rank = .true.
                        irank = irank + 1
                        ranks(irank) = rank_remote
                    end if

                    if (not(remote_right) .and. on_strip_right) on_strip_right = .false.
                    if (not(remote_left ) .and. on_strip_left ) on_strip_left  = .false.

                    if (remote_right) then
                        if (on_strip_right) then
                            je_exch(istrip, irank) = je_exch(istrip, irank) + 1
                        else
                            on_strip_right = .true.
                            
                            istrip = istrip + 1

                            js_exch (istrip, irank) = j
                            je_exch (istrip, irank) = j
                            irb     (istrip, irank) = i + 1

                            remote_on_right(istrip, irank) = .true.
                        end if

                        message_size = message_size + 1
                    end if

                    if (remote_left) then
                        if (on_strip_left) then
                            je_exch(istrip, irank) = je_exch(istrip, irank) + 1
                        else
                            on_strip_left  = .true.
                            
                            istrip = istrip + 1

                            js_exch (istrip, irank) = j
                            je_exch (istrip, irank) = j
                            irb     (istrip, irank) = i + 1
                            
                            remote_on_right(istrip, irank) = .false.
                        end if

                        message_size = message_size + 1
                    end if                    
                end do
            end do

            if (active_rank) nstrips(irank) = istrip
            if (message_size > max_message_size) max_message_size = message_size
        end do

        call this.halo_x.init(irank, ranks, nstrips, &
            irb, js_exch, je_exch, remote_on_right, max_message_size)

    end subroutine find_halo_x

    subroutine find_halo_y(this, partition)
        class(mpicom_t)  , intent(inout) :: this
        type(partition_t), intent(in)    :: partition

        ! rank depending on point (i,j), -1 land
        integer(kind=4), dimension(this.is_g-1:this.ie_g+1,this.js_g-1:this.je_g+1) :: rankij


        integer(kind=4), dimension(     this.size) :: ranks, nstrips
        integer(kind=4), dimension(1000,this.size) :: jtb, is_exch, ie_exch
        logical(kind=1), dimension(1000,this.size) :: remote_on_top
        logical(kind=1), dimension(1000,this.size) :: right_border_send, left_border_send
        logical(kind=1), dimension(1000,this.size) :: right_border_recv, left_border_recv

        integer(kind=4) :: i, j, irank , istrip
        integer(kind=4) :: rank_remote , rank_local
        integer(kind=4) :: message_size, max_message_size

        logical(kind=1) :: cond1, cond2, cond3
        logical(kind=1) :: remote_top  , remote_bottom
        logical(kind=1) :: on_strip_top, on_strip_bottom

        logical(kind=1) :: active_rank

        call partition.get_rankij(rankij,1,1)

        right_border_send = .false.
        right_border_recv = .false.
        left_border_send  = .false.
        left_border_recv  = .false.

        irank = 0
        rank_local = this.rank
        max_message_size = 0
        do rank_remote = 0, this.size-1
            if (rank_remote == rank_local) cycle
            
            istrip       = 0
            message_size = 0
            active_rank  = .false.
            
            do j = this.js-1, this.je
                on_strip_top    = .false.
                on_strip_bottom = .false.
                do i = this.is, this.ie+1
                    cond1 = (rankij(i,j  ) == rank_local )
                    cond2 = (rankij(i,j+1) == rank_remote)
                    remote_top    = cond1 .and. cond2

                    cond1 = (rankij(i,j+1) == rank_local )
                    cond2 = (rankij(i,j  ) == rank_remote)
                    remote_bottom = cond1 .and. cond2
                    
                    if ((remote_top .or. remote_bottom) .and. not(active_rank)) then
                        active_rank = .true.
                        irank = irank + 1
                        ranks(irank) = rank_remote
                    end if            

                    if (not(remote_top   ) .and. on_strip_top   ) then
                        on_strip_top    = .false.

                        cond1 = (rankij(i,j  ) /= rank_local )
                        cond2 = (rankij(i,j+1) /= rank_local )
                        cond3 = (rankij(i,j+1)  >         -1 )
                        right_border_recv(istrip, irank) = cond1 .and. cond2 .and. cond3

                        cond1 = (rankij(i,j  ) /= rank_remote)
                        cond2 = (rankij(i,j+1) /= rank_remote)
                        cond3 = (rankij(i,j  )  >         -1 )
                        right_border_send(istrip, irank) = cond1 .and. cond2 .and. cond3

                        if (right_border_recv(istrip, irank) .or. right_border_send(istrip, irank)) then
                                message_size = message_size + this.max_gcx
                        end if
                    end if
                
                    if (not(remote_bottom) .and. on_strip_bottom) then
                        on_strip_bottom = .false.

                        cond1 = (rankij(i,j  ) /= rank_local )
                        cond2 = (rankij(i,j+1) /= rank_local )
                        cond3 = (rankij(i,j  )  >         -1 )
                        right_border_recv(istrip, irank) = cond1 .and. cond2 .and. cond3

                        cond1 = (rankij(i,j  ) /= rank_remote)
                        cond2 = (rankij(i,j+1) /= rank_remote)
                        cond3 = (rankij(i,j+1)  >         -1 )
                        right_border_send(istrip, irank) = cond1 .and. cond2 .and. cond3

                        if (right_border_recv(istrip, irank) .or. right_border_send(istrip, irank)) then
                                message_size = message_size + this.max_gcx
                        end if
                    end if


                    if (remote_top) then
                        if (on_strip_top) then
                            ie_exch(istrip, irank) = ie_exch(istrip, irank) + 1
                        else
                            on_strip_top    = .true.

                            istrip = istrip + 1

                            is_exch (istrip, irank) = i
                            ie_exch (istrip, irank) = i
                            jtb     (istrip, irank) = j + 1

                            remote_on_top(istrip, irank) = .true.

                            cond1 = (rankij(i-1,j  ) /= rank_local )
                            cond2 = (rankij(i-1,j+1) /= rank_local )
                            cond3 = (rankij(i-1,j+1)  >         -1 )
                            left_border_recv(istrip, irank) = cond1 .and. cond2 .and. cond3

                            cond1 = (rankij(i-1,j  ) /= rank_remote)
                            cond2 = (rankij(i-1,j+1) /= rank_remote)
                            cond3 = (rankij(i-1,j  )  >         -1 )
                            left_border_send(istrip, irank) = cond1 .and. cond2 .and. cond3

                            if (left_border_recv(istrip, irank) .or. left_border_send(istrip, irank)) then
                                message_size = message_size + this.max_gcx
                            end if
                        end if

                        message_size = message_size + 1
                    end if

                    if (remote_bottom) then
                        if (on_strip_bottom) then
                            ie_exch(istrip, irank) = ie_exch(istrip, irank) + 1
                        else
                            on_strip_bottom = .true.

                            istrip = istrip + 1

                            is_exch (istrip, irank) = i
                            ie_exch (istrip, irank) = i
                            jtb     (istrip, irank) = j + 1
                            
                            remote_on_top(istrip, irank) = .false.

                            cond1 = (rankij(i-1,j  ) /= rank_local )
                            cond2 = (rankij(i-1,j+1) /= rank_local )
                            cond3 = (rankij(i-1,j  )  >         -1 )
                            left_border_recv(istrip, irank) = cond1 .and. cond2 .and. cond3

                            cond1 = (rankij(i-1,j  ) /= rank_remote)
                            cond2 = (rankij(i-1,j+1) /= rank_remote)
                            cond3 = (rankij(i-1,j+1)  >         -1 )
                            left_border_send(istrip, irank) = cond1 .and. cond2 .and. cond3

                            if (left_border_recv(istrip, irank) .or. left_border_send(istrip, irank)) then
                                message_size = message_size + this.max_gcx
                            end if
                        end if

                        message_size = message_size + 1
                    end if 

                end do
            end do

            if (active_rank) nstrips(irank) = istrip
            if (message_size > max_message_size) max_message_size = message_size
        end do

        call this.halo_y.init(irank, ranks, nstrips, &
            jtb, is_exch, ie_exch, remote_on_top, right_border_send, left_border_send, &
            right_border_recv, left_border_recv, max_message_size)

    end subroutine find_halo_y

    subroutine find_halo_diag(this, partition)
        class(mpicom_t)  , intent(inout) :: this
        type(partition_t), intent(in)    :: partition

        ! rank depending on point (i,j), -1 land
        integer(kind=4), dimension(this.is_g-1:this.ie_g+1,this.js_g-1:this.je_g+1) :: rankij

        integer(kind=4), dimension(this.size) :: ranks_send, nsend
        integer(kind=4), dimension(this.size) :: ranks_recv, nrecv

        integer(kind=4), dimension(1000,this.size) :: i_send, i_recv
        integer(kind=4), dimension(1000,this.size) :: j_send, j_recv
        logical(kind=1), dimension(1000,this.size) :: right_send, top_send
        logical(kind=1), dimension(1000,this.size) :: right_recv, top_recv

        integer(kind=4) :: i, j
        integer(kind=4) :: irank_send , isend
        integer(kind=4) :: irank_recv , irecv
        integer(kind=4) :: rank_remote, rank_local
        integer(kind=4) :: message_size, max_message_size

        logical(kind=1) :: cond1, cond2, cond3, cond4
        logical(kind=1) :: remote_right_bottom_send, remote_right_top_send
        logical(kind=1) :: remote_right_bottom_recv, remote_right_top_recv
        logical(kind=1) :: remote_left_bottom_send , remote_left_top_send
        logical(kind=1) :: remote_left_bottom_recv , remote_left_top_recv

        logical(kind=1) :: active_rank_send, send_flag
        logical(kind=1) :: active_rank_recv, recv_flag

        call partition.get_rankij(rankij,1,1)

        right_send = .false.
        right_recv = .false.
        top_send   = .false.
        top_recv   = .false. 

        irank_send = 0
        irank_recv = 0
        rank_local = this.rank
        max_message_size = 0
        do rank_remote = 0, this.size-1
            if (rank_remote == rank_local) cycle
            
            isend = 0
            irecv = 0
            message_size = 0

            active_rank_send = .false.
            active_rank_recv = .false.
            
            do j = this.js-1, this.je
                do i = this.is, this.ie
                    cond1 = (rankij(i-1,j  ) /= rank_local )
                    cond2 = (rankij(i  ,j  ) == rank_local )
                    cond3 = (rankij(i-1,j+1) == rank_remote)
                    cond4 = (rankij(i  ,j+1) ==         -1 )
                    remote_left_top_recv      = cond1 .and. cond2 .and. cond3 .and. cond4

                    cond1 = (rankij(i  ,j+1) /= rank_remote)
                    cond2 = (rankij(i  ,j  ) == rank_local )
                    cond3 = (rankij(i-1,j+1) == rank_remote)
                    cond4 = (rankij(i-1,j  ) ==         -1 )
                    remote_left_top_send      = cond1 .and. cond2 .and. cond3 .and. cond4


                    cond1 = (rankij(i+1,j  ) /= rank_local )
                    cond2 = (rankij(i  ,j  ) == rank_local )
                    cond3 = (rankij(i+1,j+1) == rank_remote)
                    cond4 = (rankij(i  ,j+1) ==         -1 )
                    remote_right_top_recv     = cond1 .and. cond2 .and. cond3 .and. cond4

                    cond1 = (rankij(i  ,j+1) /= rank_remote)
                    cond2 = (rankij(i  ,j  ) == rank_local )
                    cond3 = (rankij(i+1,j+1) == rank_remote)
                    cond4 = (rankij(i+1,j  ) ==         -1 )
                    remote_right_top_send     = cond1 .and. cond2 .and. cond3 .and. cond4


                    cond1 = (rankij(i-1,j+1) /= rank_local )
                    cond2 = (rankij(i  ,j+1) == rank_local )
                    cond3 = (rankij(i-1,j  ) == rank_remote)
                    cond4 = (rankij(i  ,j  ) ==         -1 )
                    remote_left_bottom_recv   = cond1 .and. cond2 .and. cond3 .and. cond4

                    cond1 = (rankij(i  ,j  ) /= rank_remote)
                    cond2 = (rankij(i  ,j+1) == rank_local )
                    cond3 = (rankij(i-1,j  ) == rank_remote)
                    cond4 = (rankij(i-1,j+1) ==         -1 )
                    remote_left_bottom_send   = cond1 .and. cond2 .and. cond3 .and. cond4


                    cond1 = (rankij(i+1,j+1) /= rank_local )
                    cond2 = (rankij(i  ,j+1) == rank_local )
                    cond3 = (rankij(i+1,j  ) == rank_remote)
                    cond4 = (rankij(i  ,j  ) ==         -1 )
                    remote_right_bottom_recv  = cond1 .and. cond2 .and. cond3 .and. cond4

                    cond1 = (rankij(i  ,j  ) /= rank_remote)
                    cond2 = (rankij(i  ,j+1) == rank_local )
                    cond3 = (rankij(i+1,j  ) == rank_remote)
                    cond4 = (rankij(i+1,j+1) ==         -1 )
                    remote_right_bottom_send  = cond1 .and. cond2 .and. cond3 .and. cond4

                    send_flag = remote_left_top_send .or. remote_right_top_send   & 
                        .or. remote_left_bottom_send .or. remote_right_bottom_send
                    recv_flag = remote_left_top_recv .or. remote_right_top_recv   & 
                        .or. remote_left_bottom_recv .or. remote_right_bottom_recv

                    if (send_flag .and. not(active_rank_send)) then
                        active_rank_send = .true.
                        irank_send = irank_send + 1
                        ranks_send(irank_send) = rank_remote
                    end if
                    
                    if (recv_flag .and. not(active_rank_recv)) then
                        active_rank_recv = .true.
                        irank_recv = irank_recv + 1
                        ranks_recv(irank_recv) = rank_remote
                    end if
                    

                    if (remote_left_top_recv) then
                        irecv = irecv + 1

                        i_recv(irecv, irank_recv) = i
                        j_recv(irecv, irank_recv) = j+1

                        right_recv(irecv, irank_recv) = .false.
                        top_recv  (irecv, irank_recv) = .true.

                        message_size = message_size + 1
                    end if

                    if (remote_left_top_send) then
                        isend = isend + 1

                        i_send(isend, irank_send) = i
                        j_send(isend, irank_send) = j+1

                        right_send(isend, irank_send) = .false.
                        top_send  (isend, irank_send) = .true.

                        message_size = message_size + 1
                    end if


                    if (remote_right_top_recv) then
                        irecv = irecv + 1

                        i_recv(irecv, irank_recv) = i+1
                        j_recv(irecv, irank_recv) = j+1

                        right_recv(irecv, irank_recv) = .true.
                        top_recv  (irecv, irank_recv) = .true.

                        message_size = message_size + 1
                    end if

                    if (remote_right_top_send) then
                        isend = isend + 1

                        i_send(isend, irank_send) = i+1
                        j_send(isend, irank_send) = j+1

                        right_send(isend, irank_send) = .true.
                        top_send  (isend, irank_send) = .true.

                        message_size = message_size + 1
                    end if


                    if (remote_left_bottom_recv) then
                        irecv = irecv + 1

                        i_recv(irecv, irank_recv) = i
                        j_recv(irecv, irank_recv) = j+1

                        right_recv(irecv, irank_recv) = .false.
                        top_recv  (irecv, irank_recv) = .false.

                        message_size = message_size + 1
                    end if

                    if (remote_left_bottom_send) then
                        isend = isend + 1

                        i_send(isend, irank_send) = i
                        j_send(isend, irank_send) = j+1

                        right_send(isend, irank_send) = .false.
                        top_send  (isend, irank_send) = .false.

                        message_size = message_size + 1
                    end if


                    if (remote_right_bottom_recv) then
                        irecv = irecv + 1

                        i_recv(irecv, irank_recv) = i+1
                        j_recv(irecv, irank_recv) = j+1

                        right_recv(irecv, irank_recv) = .true.
                        top_recv  (irecv, irank_recv) = .false.

                        message_size = message_size + 1
                    end if

                    if (remote_right_bottom_send) then
                        isend = isend + 1

                        i_send(isend, irank_send) = i+1
                        j_send(isend, irank_send) = j+1

                        right_send(isend, irank_send) = .true.
                        top_send  (isend, irank_send) = .false.

                        message_size = message_size + 1
                    end if                  

                end do
            end do

            if (active_rank_send) nsend(irank_send) = isend
            if (active_rank_recv) nrecv(irank_recv) = irecv
            if (message_size > max_message_size) max_message_size = message_size
        end do

        call this.halo_diag.init(irank_send, irank_recv, ranks_send, ranks_recv, &
        nsend, nrecv, i_send, i_recv, j_send, j_recv, right_send, right_recv, &
        top_send, top_recv, max_message_size)

    end subroutine find_halo_diag

    subroutine exchange_x(this, halo_x, data, exch_x)
        class(mpicom_t), intent(in)    :: this
        type (halo_x_t), intent(inout) :: halo_x
        integer(kind=4), intent(in   ) :: exch_x
        type(data_t)   , intent(inout) :: data

        integer(kind=4) :: irank, istrip
        integer(kind=4) :: message_size, max_message_size
        integer(kind=4) :: rank

        integer(kind=4) :: is_exch, ie_exch
        integer(kind=4) :: js_exch, je_exch

        if ((exch_x == 0) .or. (this.size == 1)) return

        max_message_size = halo_x.max_message_size * exch_x * data.message_length()
        call halo_x.allocate_buf(max_message_size)

        do irank = 1, halo_x.num_ranks
            rank = halo_x.ranks(irank)
            call MPI_Irecv(halo_x.recv_buf(1,irank), max_message_size, data.mpi_dtype, &
                rank, 0, this.com, halo_x.recv_req(irank), this.ierr)
        end do

        do irank = 1, halo_x.num_ranks
            rank = halo_x.ranks(irank)
            message_size = 0
            do istrip = 1, halo_x.nstrips(irank)
                if (halo_x.remote_on_right(istrip, irank)) then
                    is_exch = halo_x.i(istrip, irank) - exch_x
                    ie_exch = halo_x.i(istrip, irank) - 1                    
                else
                    is_exch = halo_x.i(istrip, irank)
                    ie_exch = halo_x.i(istrip, irank) + exch_x - 1
                end if
                js_exch = halo_x.js(istrip, irank)
                je_exch = halo_x.je(istrip, irank)

                call data.get_sub_array(halo_x.send_buf(:,irank), message_size, &
                    is_exch, ie_exch, js_exch, je_exch)
            end do

            call MPI_Isend(halo_x.send_buf(1,irank), message_size, data.mpi_dtype, &
                rank, 0, this.com, halo_x.send_req(irank), this.ierr)
        end do

        do
            call MPI_Waitany(halo_x.num_ranks, halo_x.recv_req, &
                irank, MPI_STATUSES_IGNORE, this.ierr)
            if (irank == MPI_UNDEFINED) exit

            message_size = 0
            do istrip = 1, halo_x.nstrips(irank)
                if (halo_x.remote_on_right(istrip, irank)) then
                    is_exch = halo_x.i(istrip, irank)
                    ie_exch = halo_x.i(istrip, irank) + exch_x - 1                    
                else
                    is_exch = halo_x.i(istrip, irank) - exch_x
                    ie_exch = halo_x.i(istrip, irank) - 1
                end if
                js_exch = halo_x.js(istrip, irank)
                je_exch = halo_x.je(istrip, irank)

                call data.put_sub_array(halo_x.recv_buf(:,irank), message_size, & 
                    is_exch, ie_exch, js_exch, je_exch)
            end do
        end do

        call MPI_Waitall(halo_x.num_ranks,halo_x.send_req,MPI_STATUSES_IGNORE,this.ierr)

    end subroutine exchange_x

    subroutine exchange_y(this, halo_y, data, exch_y)
        class(mpicom_t), intent(in)    :: this
        type (halo_y_t), intent(inout) :: halo_y
        integer(kind=4), intent(in   ) :: exch_y
        type(data_t)   , intent(inout) :: data

        integer(kind=4) :: irank, istrip
        integer(kind=4) :: message_size, max_message_size
        integer(kind=4) :: rank

        integer(kind=4) :: is_exch, ie_exch
        integer(kind=4) :: js_exch, je_exch

        if ((exch_y == 0) .or. (this.size == 1)) return

        max_message_size = halo_y.max_message_size * exch_y * data.message_length()
        call halo_y.allocate_buf(max_message_size)

        do irank = 1, halo_y.num_ranks
            rank = halo_y.ranks(irank)
            call MPI_Irecv(halo_y.recv_buf(1,irank), max_message_size, data.mpi_dtype, &
                rank, 0, this.com, halo_y.recv_req(irank), this.ierr) 
        end do

        do irank = 1, halo_y.num_ranks
            rank = halo_y.ranks(irank)
            message_size = 0
            do istrip = 1, halo_y.nstrips(irank)
                if (halo_y.remote_on_top(istrip, irank)) then
                    js_exch = halo_y.j (istrip, irank) - exch_y
                    je_exch = halo_y.j (istrip, irank) - 1                    
                else
                    js_exch = halo_y.j (istrip, irank)
                    je_exch = halo_y.j (istrip, irank) + exch_y - 1
                end if
                is_exch = halo_y.is(istrip, irank)
                ie_exch = halo_y.ie(istrip, irank)

                call data.get_sub_array(halo_y.send_buf(:,irank), message_size, &
                    is_exch, ie_exch, js_exch, je_exch)
            end do   

            call MPI_Isend(halo_y.send_buf(1,irank), message_size, data.mpi_dtype, &
                rank, 0, this.com, halo_y.send_req(irank), this.ierr)
        end do

        do 
            call MPI_Waitany(halo_y.num_ranks, halo_y.recv_req, &
                irank, MPI_STATUSES_IGNORE, this.ierr)
            if (irank == MPI_UNDEFINED) exit

            message_size = 0
            do istrip = 1, halo_y.nstrips(irank)
                if (halo_y.remote_on_top(istrip, irank)) then
                    js_exch = halo_y.j (istrip, irank)
                    je_exch = halo_y.j (istrip, irank) + exch_y - 1                    
                else
                    js_exch = halo_y.j (istrip, irank) - exch_y
                    je_exch = halo_y.j (istrip, irank) - 1
                end if
                is_exch = halo_y.is(istrip, irank)
                ie_exch = halo_y.ie(istrip, irank)

                call data.put_sub_array(halo_y.recv_buf(:,irank), message_size, & 
                    is_exch, ie_exch, js_exch, je_exch)
            end do
        end do

        call MPI_Waitall(halo_y.num_ranks,halo_y.send_req,MPI_STATUSES_IGNORE,this.ierr)

    end subroutine exchange_y

    subroutine exchange_y_halo(this, halo_y, halo_diag, data, exch_x, exch_y)
        class(mpicom_t   ), intent(in)    :: this
        type (halo_y_t   ), intent(inout) :: halo_y
        type (halo_diag_t), intent(inout) :: halo_diag
        
        type(data_t)   , intent(inout) :: data
        integer(kind=4), intent(in   ) :: exch_x, exch_y

        integer(kind=4) :: irank, istrip
        integer(kind=4) :: message_size, max_message_size, max_message_diag
        integer(kind=4) :: rank

        integer(kind=4) :: is_exch, ie_exch
        integer(kind=4) :: js_exch, je_exch

        if ((exch_y == 0) .or. (this.size == 1)) return

        max_message_size = halo_y.max_message_size    * exch_y          * data.message_length()
        max_message_diag = halo_diag.max_message_size * exch_y * exch_x * data.message_length()
        call halo_y.allocate_buf(max_message_size)
        call halo_diag.allocate_buf(max_message_diag)

        do irank = 1, halo_y.num_ranks
            rank = halo_y.ranks(irank)
            call MPI_Irecv(halo_y.recv_buf(1,irank), max_message_size, data.mpi_dtype, &
                rank, 0, this.com, halo_y.recv_req(irank), this.ierr) 
        end do

        do irank = 1, halo_diag.nranks_recv
            rank = halo_diag.ranks_recv(irank)
            call MPI_Irecv(halo_diag.recv_buf(1,irank), max_message_diag, data.mpi_dtype, &
                rank, 1, this.com, halo_diag.recv_req(irank), this.ierr) 
        end do

        do irank = 1, halo_diag.nranks_send
            rank = halo_diag.ranks_send(irank)
            message_size = 0
            do istrip = 1, halo_diag.nsend(irank)
                if (halo_diag.right_send(istrip,irank)) then
                    is_exch = halo_diag.i_send(istrip,irank) - exch_x
                    ie_exch = halo_diag.i_send(istrip,irank) - 1
                else
                    is_exch = halo_diag.i_send(istrip,irank)
                    ie_exch = halo_diag.i_send(istrip,irank) + exch_x - 1
                end if

                if (halo_diag.top_send  (istrip,irank)) then
                    js_exch = halo_diag.j_send(istrip,irank) - exch_y
                    je_exch = halo_diag.j_send(istrip,irank) - 1
                else
                    js_exch = halo_diag.j_send(istrip,irank)
                    je_exch = halo_diag.j_send(istrip,irank) + exch_y - 1
                end if

                call data.get_sub_array(halo_diag.send_buf(:,irank), message_size, &
                    is_exch, ie_exch, js_exch, je_exch)
            end do

            call MPI_Isend(halo_diag.send_buf(1,irank), message_size, data.mpi_dtype, &
                rank, 1, this.com, halo_diag.send_req(irank), this.ierr)
        end do

        do irank = 1, halo_y.num_ranks
            rank = halo_y.ranks(irank)
            message_size = 0
            do istrip = 1, halo_y.nstrips(irank)
                if (halo_y.remote_on_top(istrip, irank)) then
                    js_exch = halo_y.j (istrip, irank) - exch_y
                    je_exch = halo_y.j (istrip, irank) - 1                    
                else
                    js_exch = halo_y.j (istrip, irank)
                    je_exch = halo_y.j (istrip, irank) + exch_y - 1
                end if
                is_exch = halo_y.is(istrip, irank)
                ie_exch = halo_y.ie(istrip, irank)

                if (halo_y.right_border_send(istrip, irank)) ie_exch = ie_exch + exch_x
                if (halo_y.left_border_send (istrip, irank)) is_exch = is_exch - exch_x

                call data.get_sub_array(halo_y.send_buf(:,irank), message_size, &
                    is_exch, ie_exch, js_exch, je_exch)
            end do   

            call MPI_Isend(halo_y.send_buf(1,irank), message_size, data.mpi_dtype, &
                rank, 0, this.com, halo_y.send_req(irank), this.ierr)
        end do

        do
            call MPI_Waitany(halo_y.num_ranks, halo_y.recv_req, &
                irank, MPI_STATUSES_IGNORE, this.ierr)
            if (irank == MPI_UNDEFINED) exit

            message_size = 0
            do istrip = 1, halo_y.nstrips(irank)
                if (halo_y.remote_on_top(istrip, irank)) then
                    js_exch = halo_y.j (istrip, irank)
                    je_exch = halo_y.j (istrip, irank) + exch_y - 1                    
                else
                    js_exch = halo_y.j (istrip, irank) - exch_y
                    je_exch = halo_y.j (istrip, irank) - 1
                end if
                is_exch = halo_y.is(istrip, irank)
                ie_exch = halo_y.ie(istrip, irank)

                if (halo_y.right_border_recv(istrip, irank)) ie_exch = ie_exch + exch_x
                if (halo_y.left_border_recv (istrip, irank)) is_exch = is_exch - exch_x

                call data.put_sub_array(halo_y.recv_buf(:,irank), message_size, &
                    is_exch, ie_exch, js_exch, je_exch)
            end do
        end do

        do
            call MPI_Waitany(halo_diag.nranks_recv, halo_diag.recv_req, &
                irank, MPI_STATUSES_IGNORE, this.ierr)
            if (irank == MPI_UNDEFINED) exit

            message_size = 0
            do istrip = 1, halo_diag.nrecv(irank)
                if (halo_diag.right_recv(istrip,irank)) then
                    is_exch = halo_diag.i_recv(istrip,irank)
                    ie_exch = halo_diag.i_recv(istrip,irank) + exch_x - 1
                else
                    is_exch = halo_diag.i_recv(istrip,irank) - exch_x
                    ie_exch = halo_diag.i_recv(istrip,irank) - 1
                end if

                if (halo_diag.top_recv  (istrip,irank)) then
                    js_exch = halo_diag.j_recv(istrip,irank)
                    je_exch = halo_diag.j_recv(istrip,irank) + exch_y - 1
                else
                    js_exch = halo_diag.j_recv(istrip,irank) - exch_y
                    je_exch = halo_diag.j_recv(istrip,irank) - 1
                end if

                call data.put_sub_array(halo_diag.recv_buf(:,irank), message_size, &
                    is_exch, ie_exch, js_exch, je_exch)
            end do
        end do

        call MPI_Waitall(halo_y.num_ranks,halo_y.send_req,MPI_STATUSES_IGNORE,this.ierr)
        call MPI_Waitall(halo_diag.nranks_send,halo_diag.send_req,MPI_STATUSES_IGNORE,this.ierr)

    end subroutine exchange_y_halo

    subroutine scatter_mpicom(this, glob_data, loc_data)
        class(mpicom_t), intent(inout) :: this        
        class(data_t)  , intent(in)    :: glob_data
        class(data_t)  , intent(inout) :: loc_data

        integer(kind=4) :: displs(0:this.size), counts(0:this.size-1)
        integer(kind=4) :: message_size, mwidth, bufsize_g, bufsize_l
        integer(kind=4) :: irank, ib
        integer(kind=4) :: gcx
        real(kind=8) :: time

        time = MPI_Wtime()

        gcx = this.scatter_gcx
        
        displs(0) = 0
        do irank = 0, this.size - 1
            message_size = 0
            do ib = 1, this.nblocks(irank)
                message_size = message_size + (this.ie_gc0(ib,irank)-this.is_gc0(ib,irank)+1+2*gcx) * (this.je_gc0(ib,irank)-this.js_gc0(ib,irank)+1+2*gcx)
            end do
            displs(irank+1) = displs(irank) + message_size
            counts(irank) = message_size
        end do

        mwidth = glob_data.message_length()

        counts = counts * mwidth
        displs = displs * mwidth

        bufsize_g = sum(counts)
        bufsize_l = counts(this.rank)
        
        if ((size(this.buf_g) < bufsize_g) .and. (this.rank == 0)) then
            if (allocated(this.buf_g)) deallocate(this.buf_g)
            allocate(this.buf_g(bufsize_g)) 
        end if
        if (size(this.buf_l) < bufsize_l) then
            if (allocated(this.buf_l)) deallocate(this.buf_l)
            allocate(this.buf_l(bufsize_l)) 
        end if

        if (this.rank == 0) then
            message_size = 0
            do irank = 0, this.size - 1
                do ib = 1, this.nblocks(irank)
                    call glob_data.get_sub_array(this.buf_g(:), message_size, &
                        this.is_gc0(ib,irank)-gcx, this.ie_gc0(ib,irank)+gcx, this.js_gc0(ib,irank)-gcx, this.je_gc0(ib,irank)+gcx)
                end do
            end do
        end if

        this.time_gs_copy = this.time_gs_copy + (MPI_Wtime()-time)
        
        if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
        time = MPI_Wtime()
        
        call MPI_Scatterv(this.buf_g, counts, displs, glob_data.mpi_dtype, &
            this.buf_l, counts(this.rank), glob_data.mpi_dtype, 0, this.com, this.ierr)
            
        this.time_gs_send = this.time_gs_send + (MPI_Wtime() - time)

        time = MPI_Wtime()
        
        message_size = 0
        irank = this.rank
        do ib = 1, this.nblocks(irank)
            call loc_data.put_sub_array(this.buf_l(:), message_size, &
                this.is_gc0(ib,irank)-gcx, this.ie_gc0(ib,irank)+gcx, this.js_gc0(ib,irank)-gcx, this.je_gc0(ib,irank)+gcx)
        end do
        
        this.time_gs_copy = this.time_gs_copy + (MPI_Wtime()-time)
    
    end subroutine scatter_mpicom

    subroutine gather_mpicom(this, glob_data, loc_data)
        class(mpicom_t), intent(inout) :: this        
        class(data_t)  , intent(inout) :: glob_data
        class(data_t)  , intent(in)    :: loc_data

        integer(kind=4) :: displs(0:this.size), counts(0:this.size-1)
        integer(kind=4) :: message_size, mwidth, bufsize_g, bufsize_l
        integer(kind=4) :: irank, ib
        real(kind=8) :: time
        
        time = MPI_Wtime()

        displs(0) = 0
        do irank = 0, this.size - 1
            message_size = 0
            do ib = 1, this.nblocks(irank)
                message_size = message_size + (this.ie_gc0(ib,irank)-this.is_gc0(ib,irank)+1) * (this.je_gc0(ib,irank)-this.js_gc0(ib,irank)+1)
            end do
            displs(irank+1) = displs(irank) + message_size
            counts(irank) = message_size
        end do

        mwidth = glob_data.message_length()

        counts = counts * mwidth
        displs = displs * mwidth

        bufsize_g = sum(counts)
        bufsize_l = counts(this.rank)

        if ((size(this.buf_g) < bufsize_g) .and. (this.rank == 0)) then
            if (allocated(this.buf_g)) deallocate(this.buf_g)
            allocate(this.buf_g(bufsize_g)) 
        end if
        if (size(this.buf_l) < bufsize_l) then
            if (allocated(this.buf_l)) deallocate(this.buf_l)
            allocate(this.buf_l(bufsize_l)) 
        end if

        message_size = 0
        irank = this.rank
        do ib = 1, this.nblocks(irank)
            call loc_data.get_sub_array(this.buf_l(:), message_size, &
                this.is_gc0(ib,irank), this.ie_gc0(ib,irank), this.js_gc0(ib,irank), this.je_gc0(ib,irank))
        end do

        this.time_gs_copy = this.time_gs_copy + (MPI_Wtime()-time)
        
        if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
        time = MPI_Wtime()
        
        call MPI_Gatherv(this.buf_l, counts(this.rank), glob_data.mpi_dtype, this.buf_g, &
            counts, displs, glob_data.mpi_dtype, 0, this.com, this.ierr)
            
        this.time_gs_send = this.time_gs_send + (MPI_Wtime() - time)

        time = MPI_Wtime()
        
        if (this.rank == 0) then
            message_size = 0
            do irank = 0, this.size - 1
                do ib = 1, this.nblocks(irank)
                    call glob_data.put_sub_array(this.buf_g(:), message_size, &
                        this.is_gc0(ib,irank), this.ie_gc0(ib,irank), this.js_gc0(ib,irank), this.je_gc0(ib,irank))
                end do
            end do
        end if               

        this.time_gs_copy = this.time_gs_copy + (MPI_Wtime()-time)
        
    end subroutine gather_mpicom

    subroutine exchange_cross_2d_r8(this, array, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=8), intent(inout) :: array(:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_2d_r8_t) :: grid_function
        integer(kind=4) :: is, ie, js, je
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        call grid_function.init(array,is,ie,js,je)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL8
            
            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()
            
            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y(this.halo_y, this.data, exch_y)
            
            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_cross_2d_r8

    subroutine exchange_cross_2d_r4(this, array, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=4), intent(inout) :: array(:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_2d_r4_t) :: grid_function
        integer(kind=4) :: is, ie, js, je
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        call grid_function.init(array,is,ie,js,je)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL8
            
            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()
            
            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y(this.halo_y, this.data, exch_y)
            
            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_cross_2d_r4

    subroutine exchange_cross_3d_r8(this, array, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=8), intent(inout) :: array(:,:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_3d_r8_t) :: grid_function
        integer(kind=4) :: is, ie, js, je, ks, ke
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(array,3)
        call grid_function.init(array,is,ie,js,je,ks,ke)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL8
            
            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()
            
            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y(this.halo_y, this.data, exch_y)
            
            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_cross_3d_r8

    subroutine exchange_cross_3d_r4(this, array, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=4), intent(inout) :: array(:,:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_3d_r4_t) :: grid_function
        integer(kind=4) :: is, ie, js, je, ks, ke
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(array,3)
        call grid_function.init(array,is,ie,js,je,ks,ke)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL4
            
            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()
            
            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y(this.halo_y, this.data, exch_y)
            
            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_cross_3d_r4

    subroutine exchange_cross_4d_r8(this, array, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=8), intent(inout) :: array(:,:,:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_4d_r8_t) :: grid_function
        integer(kind=4) :: ls, le, is, ie, js, je, ks, ke
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(array,4)
        ls = 1
        le = size(array,1)
        call grid_function.init(array,ls,le,is,ie,js,je,ks,ke)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL8
            
            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()
            
            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y(this.halo_y, this.data, exch_y)
            
            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_cross_4d_r8

    subroutine exchange_cross_4d_r4(this, array, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=4), intent(inout) :: array(:,:,:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_4d_r4_t) :: grid_function
        integer(kind=4) :: ls, le, is, ie, js, je, ks, ke
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(array,4)
        ls = 1
        le = size(array,1)
        call grid_function.init(array,ls,le,is,ie,js,je,ks,ke)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL4
            
            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()
            
            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y(this.halo_y, this.data, exch_y)
            
            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_cross_4d_r4

    subroutine exchange_halo_2d_r8(this, array, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=8), intent(inout) :: array(:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_2d_r8_t) :: grid_function
        integer(kind=4) :: is, ie, js, je
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        call grid_function.init(array,is,ie,js,je)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL8

            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()

            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y_halo(this.halo_y, this.halo_diag, this.data, exch_x, exch_y)

            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_halo_2d_r8

    subroutine exchange_halo_2d_r4(this, array, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=4), intent(inout) :: array(:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_2d_r4_t) :: grid_function
        integer(kind=4) :: is, ie, js, je
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        call grid_function.init(array,is,ie,js,je)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL8

            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()

            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y_halo(this.halo_y, this.halo_diag, this.data, exch_x, exch_y)

            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_halo_2d_r4

    subroutine exchange_halo_3d_r8(this, array, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=8), intent(inout) :: array(:,:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_3d_r8_t) :: grid_function
        integer(kind=4) :: is, ie, js, je, ks, ke
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(array,3)
        call grid_function.init(array,is,ie,js,je,ks,ke)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL8

            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()

            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y_halo(this.halo_y, this.halo_diag, this.data, exch_x, exch_y)

            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_halo_3d_r8

    subroutine exchange_halo_3d_r4(this, array, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=4), intent(inout) :: array(:,:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_3d_r4_t) :: grid_function
        integer(kind=4) :: is, ie, js, je, ks, ke
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(array,3)
        call grid_function.init(array,is,ie,js,je,ks,ke)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL4

            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()

            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y_halo(this.halo_y, this.halo_diag, this.data, exch_x, exch_y)

            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_halo_3d_r4
    
    subroutine exchange_halo_3d_r8_depth(this, array, km2_l, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=8), intent(inout) :: array(:,:,:)
        integer(kind=4), intent(in)    :: km2_l(:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_3d_r8_depth_t) :: grid_function
        integer(kind=4) :: is, ie, js, je, ks, ke
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(array,3)
        call grid_function.init(array,km2_l,is,ie,js,je,ks,ke)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL8

            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()

            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y_halo(this.halo_y, this.halo_diag, this.data, exch_x, exch_y)

            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_halo_3d_r8_depth

    subroutine exchange_halo_3d_r4_depth(this, array, km2_l, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=4), intent(inout) :: array(:,:,:)
        integer(kind=4), intent(in)    :: km2_l(:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_3d_r4_depth_t) :: grid_function
        integer(kind=4) :: is, ie, js, je, ks, ke
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(array,3)
        call grid_function.init(array,km2_l,is,ie,js,je,ks,ke)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL4

            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()

            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y_halo(this.halo_y, this.halo_diag, this.data, exch_x, exch_y)

            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_halo_3d_r4_depth

    subroutine exchange_halo_4d_r8(this, array, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=8), intent(inout) :: array(:,:,:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_4d_r8_t) :: grid_function
        integer(kind=4) :: ls, le, is, ie, js, je, ks, ke
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(array,4)
        ls = 1
        le = size(array,1)
        call grid_function.init(array,ls,le,is,ie,js,je,ks,ke)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL8

            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()

            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y_halo(this.halo_y, this.halo_diag, this.data, exch_x, exch_y)

            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_halo_4d_r8

    subroutine exchange_halo_4d_r4(this, array, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=4), intent(inout) :: array(:,:,:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_4d_r4_t) :: grid_function
        integer(kind=4) :: ls, le, is, ie, js, je, ks, ke
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(array,4)
        ls = 1
        le = size(array,1)
        call grid_function.init(array,ls,le,is,ie,js,je,ks,ke)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL4

            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()

            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y_halo(this.halo_y, this.halo_diag, this.data, exch_x, exch_y)

            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_halo_4d_r4
    
    subroutine exchange_halo_4d_r8_depth(this, array, km2_l, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=8), intent(inout) :: array(:,:,:,:)
        integer(kind=4), intent(in)    :: km2_l(:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_4d_r8_depth_t) :: grid_function
        integer(kind=4) :: ls, le, is, ie, js, je, ks, ke
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(array,4)
        ls = 1
        le = size(array,1)
        call grid_function.init(array,km2_l,ls,le,is,ie,js,je,ks,ke)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL8

            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()

            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y_halo(this.halo_y, this.halo_diag, this.data, exch_x, exch_y)

            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_halo_4d_r8_depth

    subroutine exchange_halo_4d_r4_depth(this, array, km2_l, gcx, gcy, exch_x, exch_y, wait)
        class(mpicom_t), intent(inout) :: this
        integer(kind=4), intent(in   ) :: gcx, gcy, exch_x, exch_y
        real   (kind=4), intent(inout) :: array(:,:,:,:)
        integer(kind=4), intent(in)    :: km2_l(:,:)
        character(*)   , optional, intent(in) :: wait

        type(grid_function_4d_r4_depth_t) :: grid_function
        integer(kind=4) :: ls, le, is, ie, js, je, ks, ke
        real(kind=8) :: time

        if ((exch_x > this.max_gcx) .or. (exch_y > this.max_gcy)) then
            write(*,*) 'mpicom error: use smaller exch_x, exch_y'
            stop
        end if

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(array,4)
        ls = 1
        le = size(array,1)
        call grid_function.init(array,km2_l,ls,le,is,ie,js,je,ks,ke)
        call this.data.push(grid_function)

        if (not(present(wait))) then
            this.data.mpi_dtype = MPI_REAL4

            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()

            call this.exchange_x(this.halo_x, this.data, exch_x)
            call this.exchange_y_halo(this.halo_y, this.halo_diag, this.data, exch_x, exch_y)

            this.time_exch = this.time_exch + (MPI_Wtime()-time)

            call this.data.deallocate()
            this.data.num_grid_functions = 0
        end if

    end subroutine exchange_halo_4d_r4_depth

    subroutine scatter_2d_r8(this, global_array, gcx_g, gcy_g, local_array, gcx, gcy, wait)
        class(mpicom_t), intent(inout) :: this        
        real(kind=8)   , intent(in)    :: global_array(:,:)
        real(kind=8)   , intent(inout) :: local_array (:,:)
        integer(kind=4), intent(in)    :: gcx, gcy, gcx_g, gcy_g
        character(*)   , optional, intent(in) :: wait

        type(grid_function_2d_r8_t) :: glob_func, loc_func
        integer(kind=4) :: is, ie, js, je
        real(kind=8) :: time
        
        this.scatter_gcx = min(gcx_g,gcy_g,gcx,gcy)

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        call loc_func.init(local_array,is,ie,js,je)
        call this.data.push(loc_func)

        is = this.is_g - gcx_g
        ie = this.ie_g + gcx_g
        js = this.js_g - gcy_g
        je = this.je_g + gcy_g
        call glob_func.init(global_array,is,ie,js,je)
        call this.glob_data.push(glob_func)

        if (not(present(wait))) then
          this.data.mpi_dtype      = MPI_REAL8
          this.glob_data.mpi_dtype = MPI_REAL8

          if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
          time = MPI_Wtime()
          
          call this.scatter_mpicom(this.glob_data, this.data)
          
          this.time_gs = this.time_gs + (MPI_Wtime()-time)

          call this.data.deallocate()
          call this.glob_data.deallocate()
          this.data.num_grid_functions      = 0
          this.glob_data.num_grid_functions = 0
        end if

    end subroutine scatter_2d_r8

    subroutine scatter_2d_r4(this, global_array, gcx_g, gcy_g, local_array, gcx, gcy, wait)
        class(mpicom_t), intent(inout) :: this        
        real(kind=4)   , intent(in)    :: global_array(:,:)
        real(kind=4)   , intent(inout) :: local_array (:,:)
        integer(kind=4), intent(in)    :: gcx, gcy, gcx_g, gcy_g
        character(*)   , optional, intent(in) :: wait

        type(grid_function_2d_r4_t) :: glob_func, loc_func
        integer(kind=4) :: is, ie, js, je
        real(kind=8) :: time
        
        this.scatter_gcx = min(gcx_g,gcy_g,gcx,gcy)

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        call loc_func.init(local_array,is,ie,js,je)
        call this.data.push(loc_func)

        is = this.is_g - gcx_g
        ie = this.ie_g + gcx_g
        js = this.js_g - gcy_g
        je = this.je_g + gcy_g
        call glob_func.init(global_array,is,ie,js,je)
        call this.glob_data.push(glob_func)

        if (not(present(wait))) then
          this.data.mpi_dtype      = MPI_REAL4
          this.glob_data.mpi_dtype = MPI_REAL4

          if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
          time = MPI_Wtime()
          
          call this.scatter_mpicom(this.glob_data, this.data)
          
          this.time_gs = this.time_gs + (MPI_Wtime()-time)

          call this.data.deallocate()
          call this.glob_data.deallocate()
          this.data.num_grid_functions      = 0
          this.glob_data.num_grid_functions = 0
        end if

    end subroutine scatter_2d_r4

    subroutine scatter_3d_r8(this, global_array, gcx_g, gcy_g, local_array, gcx, gcy, wait)
        class(mpicom_t), intent(inout) :: this        
        real(kind=8)   , intent(in)    :: global_array(:,:,:)
        real(kind=8)   , intent(inout) :: local_array (:,:,:)
        integer(kind=4), intent(in)    :: gcx, gcy, gcx_g, gcy_g
        character(*)   , optional, intent(in) :: wait

        type(grid_function_3d_r8_t) :: glob_func, loc_func
        integer(kind=4) :: is, ie, js, je, ks, ke
        real(kind=8) :: time
        
        this.scatter_gcx = min(gcx_g,gcy_g,gcx,gcy)

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(local_array,3)
        call loc_func.init(local_array,is,ie,js,je,ks,ke)
        call this.data.push(loc_func)

        is = this.is_g - gcx_g
        ie = this.ie_g + gcx_g
        js = this.js_g - gcy_g
        je = this.je_g + gcy_g
        ks = 1
        ke = size(global_array,3)
        call glob_func.init(global_array,is,ie,js,je,ks,ke)
        call this.glob_data.push(glob_func)

        if (not(present(wait))) then
          this.data.mpi_dtype      = MPI_REAL8
          this.glob_data.mpi_dtype = MPI_REAL8
          
          if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
          time = MPI_Wtime()

          call this.scatter_mpicom(this.glob_data, this.data)
          
          this.time_gs = this.time_gs + (MPI_Wtime()-time)

          call this.data.deallocate()
          call this.glob_data.deallocate()
          this.data.num_grid_functions      = 0
          this.glob_data.num_grid_functions = 0
        end if

    end subroutine scatter_3d_r8

    subroutine scatter_3d_r4(this, global_array, gcx_g, gcy_g, local_array, gcx, gcy, wait)
        class(mpicom_t), intent(inout) :: this        
        real(kind=4)   , intent(in)    :: global_array(:,:,:)
        real(kind=4)   , intent(inout) :: local_array (:,:,:)
        integer(kind=4), intent(in)    :: gcx, gcy, gcx_g, gcy_g
        character(*)   , optional, intent(in) :: wait

        type(grid_function_3d_r4_t) :: glob_func, loc_func
        integer(kind=4) :: is, ie, js, je, ks, ke
        real(kind=8) :: time
        
        this.scatter_gcx = min(gcx_g,gcy_g,gcx,gcy)

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(local_array,3)
        call loc_func.init(local_array,is,ie,js,je,ks,ke)
        call this.data.push(loc_func)

        is = this.is_g - gcx_g
        ie = this.ie_g + gcx_g
        js = this.js_g - gcy_g
        je = this.je_g + gcy_g
        ks = 1
        ke = size(global_array,3)
        call glob_func.init(global_array,is,ie,js,je,ks,ke)
        call this.glob_data.push(glob_func)
        
        if (not(present(wait))) then
          this.data.mpi_dtype      = MPI_REAL4
          this.glob_data.mpi_dtype = MPI_REAL4
          
          if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
          time = MPI_Wtime()

          call this.scatter_mpicom(this.glob_data, this.data)
          
          this.time_gs = this.time_gs + (MPI_Wtime()-time)

          call this.data.deallocate()
          call this.glob_data.deallocate()
          this.data.num_grid_functions      = 0
          this.glob_data.num_grid_functions = 0
        end if

    end subroutine scatter_3d_r4

    subroutine scatter_4d_r8(this, global_array, gcx_g, gcy_g, local_array, gcx, gcy, wait)
        class(mpicom_t), intent(inout) :: this        
        real(kind=8)   , intent(in)    :: global_array(:,:,:,:)
        real(kind=8)   , intent(inout) :: local_array (:,:,:,:)
        integer(kind=4), intent(in)    :: gcx, gcy, gcx_g, gcy_g
        character(*)   , optional, intent(in) :: wait

        type(grid_function_4d_r8_t) :: glob_func, loc_func
        integer(kind=4) :: is, ie, js, je, ks, ke, ls, le
        real(kind=8) :: time
        
        this.scatter_gcx = min(gcx_g,gcy_g,gcx,gcy)

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(local_array,4)
        ls = 1
        le = size(local_array,1)
        call loc_func.init(local_array,ls,le,is,ie,js,je,ks,ke)
        call this.data.push(loc_func)

        is = this.is_g - gcx_g
        ie = this.ie_g + gcx_g
        js = this.js_g - gcy_g
        je = this.je_g + gcy_g
        ks = 1
        ke = size(global_array,4)
        ls = 1
        le = size(global_array,1)
        call glob_func.init(global_array,ls,le,is,ie,js,je,ks,ke)
        call this.glob_data.push(glob_func)

        if (not(present(wait))) then
          this.data.mpi_dtype      = MPI_REAL8
          this.glob_data.mpi_dtype = MPI_REAL8

          if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
          time = MPI_Wtime()
          
          call this.scatter_mpicom(this.glob_data, this.data)
          
          this.time_gs = this.time_gs + (MPI_Wtime()-time)

          call this.data.deallocate()
          call this.glob_data.deallocate()
          this.data.num_grid_functions      = 0
          this.glob_data.num_grid_functions = 0
        end if

    end subroutine scatter_4d_r8

    subroutine scatter_4d_r4(this, global_array, gcx_g, gcy_g, local_array, gcx, gcy, wait)
        class(mpicom_t), intent(inout) :: this        
        real(kind=4)   , intent(in)    :: global_array(:,:,:,:)
        real(kind=4)   , intent(inout) :: local_array (:,:,:,:)
        integer(kind=4), intent(in)    :: gcx, gcy, gcx_g, gcy_g
        character(*)   , optional, intent(in) :: wait

        type(grid_function_4d_r4_t) :: glob_func, loc_func
        integer(kind=4) :: is, ie, js, je, ks, ke, ls, le
        real(kind=8) :: time
        
        this.scatter_gcx = min(gcx_g,gcy_g,gcx,gcy)

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(local_array,4)
        ls = 1
        le = size(local_array,1)
        call loc_func.init(local_array,ls,le,is,ie,js,je,ks,ke)
        call this.data.push(loc_func)

        is = this.is_g - gcx_g
        ie = this.ie_g + gcx_g
        js = this.js_g - gcy_g
        je = this.je_g + gcy_g
        ks = 1
        ke = size(global_array,4)
        ls = 1
        le = size(global_array,1)
        call glob_func.init(global_array,ls,le,is,ie,js,je,ks,ke)
        call this.glob_data.push(glob_func)

        if (not(present(wait))) then
          this.data.mpi_dtype      = MPI_REAL4
          this.glob_data.mpi_dtype = MPI_REAL4
          
          if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
          time = MPI_Wtime()

          call this.scatter_mpicom(this.glob_data, this.data)
          
          this.time_gs = this.time_gs + (MPI_Wtime()-time)

          call this.data.deallocate()
          call this.glob_data.deallocate()
          this.data.num_grid_functions      = 0
          this.glob_data.num_grid_functions = 0
        end if

    end subroutine scatter_4d_r4

    subroutine gather_2d_r8(this, global_array, gcx_g, gcy_g, local_array, gcx, gcy, wait)
        class(mpicom_t), intent(inout) :: this        
        real(kind=8)   , intent(inout) :: global_array(:,:)
        real(kind=8)   , intent(in)    :: local_array (:,:)
        integer(kind=4), intent(in)    :: gcx, gcy, gcx_g, gcy_g
        character(*)   , optional, intent(in) :: wait

        type(grid_function_2d_r8_t) :: glob_func, loc_func
        integer(kind=4) :: is, ie, js, je
        real(kind=8) :: time

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        call loc_func.init(local_array,is,ie,js,je)
        call this.data.push(loc_func)

        is = this.is_g - gcx_g
        ie = this.ie_g + gcx_g
        js = this.js_g - gcy_g
        je = this.je_g + gcy_g
        call glob_func.init(global_array,is,ie,js,je)
        call this.glob_data.push(glob_func)

        if (not(present(wait))) then
            this.data.mpi_dtype      = MPI_REAL8
            this.glob_data.mpi_dtype = MPI_REAL8
            
            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()

            call this.gather_mpicom(this.glob_data, this.data)
            
            this.time_gs = this.time_gs + (MPI_Wtime()-time)

            call this.data.deallocate()
            call this.glob_data.deallocate()
            this.data.num_grid_functions      = 0
            this.glob_data.num_grid_functions = 0
        end if

    end subroutine gather_2d_r8

    subroutine gather_2d_r4(this, global_array, gcx_g, gcy_g, local_array, gcx, gcy, wait)
        class(mpicom_t), intent(inout) :: this        
        real(kind=4)   , intent(inout)  :: global_array(:,:)
        real(kind=4)   , intent(in)    :: local_array (:,:)
        integer(kind=4), intent(in)    :: gcx, gcy, gcx_g, gcy_g
        character(*)   , optional, intent(in) :: wait

        type(grid_function_2d_r4_t) :: glob_func, loc_func
        integer(kind=4) :: is, ie, js, je
        real(kind=8) :: time

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        call loc_func.init(local_array,is,ie,js,je)
        call this.data.push(loc_func)

        is = this.is_g - gcx_g
        ie = this.ie_g + gcx_g
        js = this.js_g - gcy_g
        je = this.je_g + gcy_g
        call glob_func.init(global_array,is,ie,js,je)
        call this.glob_data.push(glob_func)

        if (not(present(wait))) then
            this.data.mpi_dtype      = MPI_REAL4
            this.glob_data.mpi_dtype = MPI_REAL4
            
            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()

            call this.gather_mpicom(this.glob_data, this.data)
            
            this.time_gs = this.time_gs + (MPI_Wtime()-time)

            call this.data.deallocate()
            call this.glob_data.deallocate()
            this.data.num_grid_functions      = 0
            this.glob_data.num_grid_functions = 0
        end if

    end subroutine gather_2d_r4

    subroutine gather_3d_r8(this, global_array, gcx_g, gcy_g, local_array, gcx, gcy, wait)
        class(mpicom_t), intent(inout) :: this        
        real(kind=8)   , intent(inout) :: global_array(:,:,:)
        real(kind=8)   , intent(in)    :: local_array (:,:,:)
        integer(kind=4), intent(in)    :: gcx, gcy, gcx_g, gcy_g
        character(*)   , optional, intent(in) :: wait

        type(grid_function_3d_r8_t) :: glob_func, loc_func
        integer(kind=4) :: is, ie, js, je, ks, ke
        real(kind=8) :: time

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(local_array,3)
        call loc_func.init(local_array,is,ie,js,je,ks,ke)
        call this.data.push(loc_func)

        is = this.is_g - gcx_g
        ie = this.ie_g + gcx_g
        js = this.js_g - gcy_g
        je = this.je_g + gcy_g
        ks = 1
        ke = size(global_array,3)
        call glob_func.init(global_array,is,ie,js,je,ks,ke)
        call this.glob_data.push(glob_func)

        if (not(present(wait))) then
            this.data.mpi_dtype      = MPI_REAL8
            this.glob_data.mpi_dtype = MPI_REAL8
            
            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()
            
            call this.gather_mpicom(this.glob_data, this.data)
            this.time_gs = this.time_gs + (MPI_Wtime()-time)

            call this.data.deallocate()
            call this.glob_data.deallocate()
            this.data.num_grid_functions      = 0
            this.glob_data.num_grid_functions = 0
        end if

    end subroutine gather_3d_r8


    subroutine gather_3d_r4(this, global_array, gcx_g, gcy_g, local_array, gcx, gcy, wait)
        class(mpicom_t), intent(inout) :: this        
        real(kind=4)   , intent(inout) :: global_array(:,:,:)
        real(kind=4)   , intent(in)    :: local_array (:,:,:)
        integer(kind=4), intent(in)    :: gcx, gcy, gcx_g, gcy_g
        character(*)   , optional, intent(in) :: wait

        type(grid_function_3d_r4_t) :: glob_func, loc_func
        integer(kind=4) :: is, ie, js, je, ks, ke
        real(kind=8) :: time
        
        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(local_array,3)
        call loc_func.init(local_array,is,ie,js,je,ks,ke)
        call this.data.push(loc_func)

        is = this.is_g - gcx_g
        ie = this.ie_g + gcx_g
        js = this.js_g - gcy_g
        je = this.je_g + gcy_g
        ks = 1
        ke = size(global_array,3)
        call glob_func.init(global_array,is,ie,js,je,ks,ke)
        call this.glob_data.push(glob_func)

        if (not(present(wait))) then
            this.data.mpi_dtype      = MPI_REAL4
            this.glob_data.mpi_dtype = MPI_REAL4
            
            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()

            call this.gather_mpicom(this.glob_data, this.data)
            
            this.time_gs = this.time_gs + (MPI_Wtime()-time)

            call this.data.deallocate()
            call this.glob_data.deallocate()
            this.data.num_grid_functions      = 0
            this.glob_data.num_grid_functions = 0
        end if

    end subroutine gather_3d_r4

    subroutine gather_4d_r8(this, global_array, gcx_g, gcy_g, local_array, gcx, gcy, wait)
        class(mpicom_t), intent(inout) :: this        
        real(kind=8)   , intent(inout) :: global_array(:,:,:,:)
        real(kind=8)   , intent(in)    :: local_array (:,:,:,:)
        integer(kind=4), intent(in)    :: gcx, gcy, gcx_g, gcy_g
        character(*)   , optional, intent(in) :: wait

        type(grid_function_4d_r8_t) :: glob_func, loc_func
        integer(kind=4) :: is, ie, js, je, ks, ke, ls, le
        real(kind=8) :: time

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(local_array,4)
        ls = 1
        le = size(local_array,1)
        call loc_func.init(local_array,ls,le,is,ie,js,je,ks,ke)
        call this.data.push(loc_func)

        is = this.is_g - gcx_g
        ie = this.ie_g + gcx_g
        js = this.js_g - gcy_g
        je = this.je_g + gcy_g
        ks = 1
        ke = size(global_array,4)
        ls = 1
        le = size(global_array,1)
        call glob_func.init(global_array,ls,le,is,ie,js,je,ks,ke)
        call this.glob_data.push(glob_func)

        if (not(present(wait))) then
            this.data.mpi_dtype      = MPI_REAL8
            this.glob_data.mpi_dtype = MPI_REAL8
            
            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()

            call this.gather_mpicom(this.glob_data, this.data)
            
            this.time_gs = this.time_gs + (MPI_Wtime()-time)

            call this.data.deallocate()
            call this.glob_data.deallocate()
            this.data.num_grid_functions      = 0
            this.glob_data.num_grid_functions = 0
        end if

    end subroutine gather_4d_r8

    subroutine gather_4d_r4(this, global_array, gcx_g, gcy_g, local_array, gcx, gcy, wait)
        class(mpicom_t), intent(inout) :: this        
        real(kind=4)   , intent(inout) :: global_array(:,:,:,:)
        real(kind=4)   , intent(in)    :: local_array (:,:,:,:)
        integer(kind=4), intent(in)    :: gcx, gcy, gcx_g, gcy_g
        character(*)   , optional, intent(in) :: wait

        type(grid_function_4d_r4_t) :: glob_func, loc_func
        integer(kind=4) :: is, ie, js, je, ks, ke, ls, le
        real(kind=8) :: time

        is = this.is - gcx
        ie = this.ie + gcx
        js = this.js - gcy
        je = this.je + gcy
        ks = 1
        ke = size(local_array,4)
        ls = 1
        le = size(local_array,1)
        call loc_func.init(local_array,ls,le,is,ie,js,je,ks,ke)
        call this.data.push(loc_func)

        is = this.is_g - gcx_g
        ie = this.ie_g + gcx_g
        js = this.js_g - gcy_g
        je = this.je_g + gcy_g
        ks = 1
        ke = size(global_array,4)
        ls = 1
        le = size(global_array,1)
        call glob_func.init(global_array,ls,le,is,ie,js,je,ks,ke)
        call this.glob_data.push(glob_func)

        if (not(present(wait))) then
            this.data.mpi_dtype      = MPI_REAL4
            this.glob_data.mpi_dtype = MPI_REAL4
            
            if (this.use_mpi_barriers) call MPI_Barrier(this.com,this.ierr)
            time = MPI_Wtime()

            call this.gather_mpicom(this.glob_data, this.data)
            
            this.time_gs = this.time_gs + (MPI_Wtime()-time)

            call this.data.deallocate()
            call this.glob_data.deallocate()
            this.data.num_grid_functions      = 0
            this.glob_data.num_grid_functions = 0
        end if

    end subroutine gather_4d_r4

    subroutine check_cross(this, partition) 
        class(mpicom_t)  , intent(inout) :: this
        type(partition_t), intent(in   ) :: partition

        real(kind=8) :: array_2d_r8(this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy)
        real(kind=4) :: array_2d_r4(this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy)
        real(kind=8) :: array_3d_r8(this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,1:3)
        real(kind=4) :: array_3d_r4(this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,1:3)
        real(kind=8) :: array_4d_r8(1:5,this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,1:3)
        real(kind=4) :: array_4d_r4(1:5,this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,1:3)

        logical(kind=1) :: mask(this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy)
        integer(kind=4) :: rankij(this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy)

        logical(kind=1) :: active_x, active_y, direction
        integer(kind=4) :: i, j, k, l, exch_x, exch_y
        real   (kind=8) :: s, s_mpi
        real   (kind=8) :: value2d, value3d, value4d

        call partition.get_rankij(rankij,this.max_gcx,this.max_gcy)
        call partition.get_local_mask(mask,this.max_gcx,this.max_gcy,this.rank)

        s_mpi = 0.0_8
        exch_x = 0
        exch_y = 1
        active_x = .true.
        active_y = .true.
        direction = .true. !previous direction, true = x, false = y
        do while(active_x .or. active_y)
            if ((exch_x == this.max_gcx) .and. (s_mpi < 1.0d-10)) then
                active_x  = .false.
                direction = .false.
            end if
            if ((exch_y == this.max_gcy) .and. (s_mpi < 1.0d-10)) then 
                active_y = .false.
                direction = .true.
            end if

            if (s_mpi > 1.0d-10) then
                if (direction) then
                    active_x = .false.
                    exch_x = exch_x - 1
                    direction = not(direction)
                else
                    active_y = .false.
                    exch_y = exch_y - 1
                    direction = not(direction)
                end if
            end if
            if (    direction  .and. active_x) exch_x = exch_x + 1
            if (not(direction) .and. active_y) exch_y = exch_y + 1

            array_2d_r8 = -1.0d+10
            array_2d_r4 = -1.0e+10
            array_3d_r8 = -1.0d+10
            array_3d_r4 = -1.0e+10
            array_4d_r8 = -1.0d+10
            array_4d_r4 = -1.0e+10

            do k = 1, 3
                do j = this.js-this.max_gcy, this.je+this.max_gcy
                    do i = this.is-this.max_gcx, this.ie+this.max_gcx
                        if ((rankij(i,j) == this.rank) .or. (rankij(i,j) == -1)) then
                            do l = 1, 5
                                array_2d_r8(i,j)     = real(i+j,8)
                                array_2d_r4(i,j)     = real(i+j,4)
                                array_3d_r8(i,j,k)   = real(i+j+k,8)
                                array_3d_r4(i,j,k)   = real(i+j+k,4)
                                array_4d_r8(l,i,j,k) = real(i+j+k+l,8)
                                array_4d_r4(l,i,j,k) = real(i+j+k+l,4)
                            end do
                        end if
                    end do
                end do
            end do

            call this.exchange_cross(array_4d_r4,this.max_gcx,this.max_gcy,exch_x,exch_y,'wait')
            call this.exchange_cross(array_2d_r4,this.max_gcx,this.max_gcy,exch_x,exch_y,'wait')
            call this.exchange_cross(array_3d_r4,this.max_gcx,this.max_gcy,exch_x,exch_y)

            call this.exchange_cross(array_4d_r8,this.max_gcx,this.max_gcy,exch_x,exch_y,'wait')
            call this.exchange_cross(array_2d_r8,this.max_gcx,this.max_gcy,exch_x,exch_y,'wait')
            call this.exchange_cross(array_3d_r8,this.max_gcx,this.max_gcy,exch_x,exch_y)

            s = 0.0_8
                do k = 1, 3
                    do j = this.js-this.max_gcy, this.je+this.max_gcy
                        do i = this.is-this.max_gcx, this.ie+this.max_gcx
                            if (mask(i,j)) then
                                do l = 1, 5
                                if (exch_x < 2 .or. exch_y < 2) then
                                    value4d = 0.125_8 *                                    &
                                        (array_4d_r8(l,i,j+1,k) + array_4d_r8(l,i,j-1,k) + & 
                                         array_4d_r8(l,i+1,j,k) + array_4d_r8(l,i-1,j,k) + &
                                         4.0_8 * array_4d_r8(l,i,j,k))                     &
                                          + 0.125_4 *                                      &
                                        (array_4d_r4(l,i,j+1,k) + array_4d_r4(l,i,j-1,k) + & 
                                         array_4d_r4(l,i+1,j,k) + array_4d_r4(l,i-1,j,k) + &
                                         4.0_4 * array_4d_r4(l,i,j,k))
                                    value3d = 0.125_8 *                                &
                                        (array_3d_r8(i,j+1,k) + array_3d_r8(i,j-1,k) + & 
                                         array_3d_r8(i+1,j,k) + array_3d_r8(i-1,j,k) + &
                                         4.0_8 * array_3d_r8(i,j,k))                   &
                                          + 0.125_4 *                                  &
                                        (array_3d_r4(i,j+1,k) + array_3d_r4(i,j-1,k) + & 
                                         array_3d_r4(i+1,j,k) + array_3d_r4(i-1,j,k) + &
                                         4.0_4 * array_3d_r4(i,j,k))
                                    value2d = 0.125_8 *                            &
                                        (array_2d_r8(i,j+1) + array_2d_r8(i,j-1) + & 
                                         array_2d_r8(i+1,j) + array_2d_r8(i-1,j) + &
                                         4.0_8 * array_2d_r8(i,j))                 &
                                          + 0.125_4 *                              &
                                        (array_2d_r4(i,j+1) + array_2d_r4(i,j-1) + & 
                                         array_2d_r4(i+1,j) + array_2d_r4(i-1,j) + &
                                         4.0_4 * array_2d_r4(i,j))
                                else
                                    value4d = 0.03125_8 *                                                  &
                                        (                                                                  &    
                                                 array_4d_r8(l,i,j+2,k) +         array_4d_r8(l,i,j-2,k) + &
                                                 array_4d_r8(l,i+2,j,k) +         array_4d_r8(l,i-2,j,k) + &
                                        4.0_8  * array_4d_r8(l,i,j+1,k) + 4.0_8 * array_4d_r8(l,i,j-1,k) + &
                                        4.0_8  * array_4d_r8(l,i+1,j,k) + 4.0_8 * array_4d_r8(l,i-1,j,k) + &
                                        12.0_8 * array_4d_r8(l,i,  j,k)                                    &
                                        )                                                                  &
                                          + 0.03125_4 *                                                    &
                                        (                                                                  &    
                                                 array_4d_r4(l,i,j+2,k) +         array_4d_r4(l,i,j-2,k) + &
                                                 array_4d_r4(l,i+2,j,k) +         array_4d_r4(l,i-2,j,k) + &
                                        4.0_4  * array_4d_r4(l,i,j+1,k) + 4.0_4 * array_4d_r4(l,i,j-1,k) + &
                                        4.0_4  * array_4d_r4(l,i+1,j,k) + 4.0_4 * array_4d_r4(l,i-1,j,k) + &
                                        12.0_4 * array_4d_r4(l,i,  j,k)                                    &
                                        )
                                    value3d = 0.03125_8 *                                              &
                                        (                                                              &    
                                                 array_3d_r8(i,j+2,k) +         array_3d_r8(i,j-2,k) + &
                                                 array_3d_r8(i+2,j,k) +         array_3d_r8(i-2,j,k) + &
                                        4.0_8  * array_3d_r8(i,j+1,k) + 4.0_8 * array_3d_r8(i,j-1,k) + &
                                        4.0_8  * array_3d_r8(i+1,j,k) + 4.0_8 * array_3d_r8(i-1,j,k) + &
                                        12.0_8 * array_3d_r8(i,  j,k)                                  &
                                        )                                                              &
                                          + 0.03125_4 *                                                &
                                        (                                                              &    
                                                 array_3d_r4(i,j+2,k) +         array_3d_r4(i,j-2,k) + &
                                                 array_3d_r4(i+2,j,k) +         array_3d_r4(i-2,j,k) + &
                                        4.0_4  * array_3d_r4(i,j+1,k) + 4.0_4 * array_3d_r4(i,j-1,k) + &
                                        4.0_4  * array_3d_r4(i+1,j,k) + 4.0_4 * array_3d_r4(i-1,j,k) + &
                                        12.0_4 * array_3d_r4(i,  j,k)                                  &
                                        )
                                    value2d = 0.03125_8 *                                          &
                                        (                                                          &    
                                                 array_2d_r8(i,j+2) +         array_2d_r8(i,j-2) + &
                                                 array_2d_r8(i+2,j) +         array_2d_r8(i-2,j) + &
                                        4.0_8  * array_2d_r8(i,j+1) + 4.0_8 * array_2d_r8(i,j-1) + &
                                        4.0_8  * array_2d_r8(i+1,j) + 4.0_8 * array_2d_r8(i-1,j) + &
                                        12.0_8 * array_2d_r8(i,  j)                                &
                                        )                                                          &
                                          + 0.03125_4 *                                            &
                                        (                                                          &    
                                                 array_2d_r4(i,j+2) +         array_2d_r4(i,j-2) + &
                                                 array_2d_r4(i+2,j) +         array_2d_r4(i-2,j) + &
                                        4.0_4  * array_2d_r4(i,j+1) + 4.0_4 * array_2d_r4(i,j-1) + &
                                        4.0_4  * array_2d_r4(i+1,j) + 4.0_4 * array_2d_r4(i-1,j) + &
                                        12.0_4 * array_2d_r4(i,  j)                                &
                                        )
                                end if

                                s = s + abs(2.0_8 * real(i+j+k,8) - value3d) + abs(2.0_8 * real(i+j,8) - value2d) + abs(2.0_8 * real(i+j+k+l,8) - value4d)
                                end do
                            end if
                        end do
                    end do
                end do

            call MPI_Allreduce(s, s_mpi, 1, MPI_REAL8, MPI_SUM, this.com, this.ierr)

            if ((s_mpi < 1.0d-10) .and. active_x .and. active_y) then
                direction = not(direction)
            end if
        end do

        ! update maximum avaliable gcx
        this.max_gcx = exch_x
        this.max_gcy = exch_y

        if (s_mpi < 1.0d-10) then
            if (this.rank == 0) write(*,*) 'mpicom exchange_cross: ok, maximum gcx, gcy:', exch_x, exch_y
        else
            if (this.rank == 0) write(*,*) 'mpicom exchange_cross: error', s_mpi
            stop
        end if

    end subroutine check_cross

    subroutine check_halo(this, partition) 
        class(mpicom_t)  , intent(inout) :: this
        type(partition_t), intent(in   ) :: partition

        real(kind=8) :: array_2d_r8(this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy)
        real(kind=4) :: array_2d_r4(this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy)
        real(kind=8) :: array_3d_r8(this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,1:3)
        real(kind=4) :: array_3d_r4(this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,1:3)
        real(kind=8) :: array_4d_r8(1:5,this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,1:3)
        real(kind=4) :: array_4d_r4(1:5,this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,1:3)
        real(kind=8) :: array_3d_r8d(this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,1:3)
        real(kind=4) :: array_3d_r4d(this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,1:3)
        real(kind=8) :: array_4d_r8d(1:5,this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,1:3)
        real(kind=4) :: array_4d_r4d(1:5,this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,1:3)

        logical(kind=1) :: mask(this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy)
        integer(kind=4) :: rankij(this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy)
        integer(kind=4) :: km2_l(this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy)

        logical(kind=1) :: active_x, active_y, direction
        integer(kind=4) :: i, j, k, l, exch_x, exch_y
        real   (kind=8) :: s, s_mpi
        real   (kind=8) :: value2d, value3d, value4d, value3dd, value4dd

        call partition.get_rankij(rankij,this.max_gcx,this.max_gcy)
        call partition.get_local_mask(mask,this.max_gcx,this.max_gcy,this.rank)

        s_mpi = 0.0_8
        exch_x = 0
        exch_y = 1
        active_x = .true.
        active_y = .true.
        direction = .true. !previous direction, true = x, false = y
        do while(active_x .or. active_y)
            if ((exch_x == this.max_gcx) .and. (s_mpi < 1.0d-10)) then
                active_x  = .false.
                direction = .false.
            end if
            if ((exch_y == this.max_gcy) .and. (s_mpi < 1.0d-10)) then 
                active_y = .false.
                direction = .true.
            end if

            if (s_mpi > 1.0d-10) then
                if (direction) then
                    active_x = .false.
                    exch_x = exch_x - 1
                    direction = not(direction)
                else
                    active_y = .false.
                    exch_y = exch_y - 1
                    direction = not(direction)
                end if
            end if
            if (    direction  .and. active_x) exch_x = exch_x + 1
            if (not(direction) .and. active_y) exch_y = exch_y + 1

            array_2d_r8  = -1.0d+10
            array_2d_r4  = -1.0e+10
            array_3d_r8  = -1.0d+10
            array_3d_r4  = -1.0e+10
            array_4d_r8  = -1.0d+10
            array_4d_r4  = -1.0e+10
            array_3d_r8d = -1.0d+10
            array_3d_r4d = -1.0e+10
            array_4d_r8d = -1.0d+10
            array_4d_r4d = -1.0e+10
            
            do k = 1, 3
                do j = this.js-this.max_gcy, this.je+this.max_gcy
                    do i = this.is-this.max_gcx, this.ie+this.max_gcx
                        if ((rankij(i,j) == this.rank) .or. (rankij(i,j) == -1)) then
                            do l = 1, 5
                                array_2d_r8(i,j)     = real(i+j,8)
                                array_2d_r4(i,j)     = real(i+j,4)
                                array_3d_r8(i,j,k)   = real(i+j+k,8)
                                array_3d_r4(i,j,k)   = real(i+j+k,4)
                                array_4d_r8(l,i,j,k) = real(i+j+k+l,8)
                                array_4d_r4(l,i,j,k) = real(i+j+k+l,4)
                                
                                array_3d_r8d(i,j,k)   = real(i+j+k,8)
                                array_3d_r4d(i,j,k)   = real(i+j+k,4)
                                array_4d_r8d(l,i,j,k) = real(i+j+k+l,8)
                                array_4d_r4d(l,i,j,k) = real(i+j+k+l,4)
                            end do
                        end if
                    end do
                end do
            end do
            
            km2_l = 3

            call this.exchange_halo(array_4d_r4,this.max_gcx,this.max_gcy,exch_x,exch_y,'wait')
            call this.exchange_halo(array_2d_r4,this.max_gcx,this.max_gcy,exch_x,exch_y,'wait')
            call this.exchange_halo(array_3d_r4,this.max_gcx,this.max_gcy,exch_x,exch_y,'wait')
            call this.exchange_halo(array_3d_r4d,km2_l,this.max_gcx,this.max_gcy,exch_x,exch_y,'wait')
            call this.exchange_halo(array_4d_r4d,km2_l,this.max_gcx,this.max_gcy,exch_x,exch_y)

            call this.exchange_halo(array_4d_r8,this.max_gcx,this.max_gcy,exch_x,exch_y,'wait')
            call this.exchange_halo(array_2d_r8,this.max_gcx,this.max_gcy,exch_x,exch_y,'wait')
            call this.exchange_halo(array_3d_r8,this.max_gcx,this.max_gcy,exch_x,exch_y,'wait')
            call this.exchange_halo(array_3d_r8d,km2_l,this.max_gcx,this.max_gcy,exch_x,exch_y,'wait')
            call this.exchange_halo(array_4d_r8d,km2_l,this.max_gcx,this.max_gcy,exch_x,exch_y)

            s = 0.0_8
            do k = 1, 3
                do j = this.js-this.max_gcy, this.je+this.max_gcy
                    do i = this.is-this.max_gcx, this.ie+this.max_gcx
                        if (mask(i,j)) then
                            do l = 1, 5 
                            if (exch_x < 2 .or. exch_y < 2) then
                                value4d = 0.0625_8 *                                                                        &
                                   (array_4d_r8(l,i+1,j+1,k) + 2.0_8 * array_4d_r8(l,i,j+1,k) + array_4d_r8(l,i-1,j+1,k)  + & 
                           2.0_8 * (array_4d_r8(l,i+1,j  ,k) + 2.0_8 * array_4d_r8(l,i,j  ,k) + array_4d_r8(l,i-1,j  ,k)) + & 
                                    array_4d_r8(l,i+1,j-1,k) + 2.0_8 * array_4d_r8(l,i,j-1,k) + array_4d_r8(l,i-1,j-1,k)) + &
                                          0.0625_4 *                                                                  &
                                   (array_4d_r4(l,i+1,j+1,k) + 2.0_8 * array_4d_r4(l,i,j+1,k) + array_4d_r4(l,i-1,j+1,k)  + & 
                           2.0_8 * (array_4d_r4(l,i+1,j  ,k) + 2.0_8 * array_4d_r4(l,i,j  ,k) + array_4d_r4(l,i-1,j  ,k)) + & 
                                    array_4d_r4(l,i+1,j-1,k) + 2.0_8 * array_4d_r4(l,i,j-1,k) + array_4d_r4(l,i-1,j-1,k))
                                value4dd = 0.0625_8 *                                                                       &
                                   (array_4d_r8d(l,i+1,j+1,k) + 2.0_8 * array_4d_r8d(l,i,j+1,k) + array_4d_r8d(l,i-1,j+1,k)  + & 
                           2.0_8 * (array_4d_r8d(l,i+1,j  ,k) + 2.0_8 * array_4d_r8d(l,i,j  ,k) + array_4d_r8d(l,i-1,j  ,k)) + & 
                                    array_4d_r8d(l,i+1,j-1,k) + 2.0_8 * array_4d_r8d(l,i,j-1,k) + array_4d_r8d(l,i-1,j-1,k)) + &
                                          0.0625_4 *                                                                  &
                                   (array_4d_r4d(l,i+1,j+1,k) + 2.0_8 * array_4d_r4d(l,i,j+1,k) + array_4d_r4d(l,i-1,j+1,k)  + & 
                           2.0_8 * (array_4d_r4d(l,i+1,j  ,k) + 2.0_8 * array_4d_r4d(l,i,j  ,k) + array_4d_r4d(l,i-1,j  ,k)) + & 
                                    array_4d_r4d(l,i+1,j-1,k) + 2.0_8 * array_4d_r4d(l,i,j-1,k) + array_4d_r4d(l,i-1,j-1,k))
                                value3d = 0.0625_8 *                                                                  &
                                   (array_3d_r8(i+1,j+1,k) + 2.0_8 * array_3d_r8(i,j+1,k) + array_3d_r8(i-1,j+1,k)  + & 
                           2.0_8 * (array_3d_r8(i+1,j  ,k) + 2.0_8 * array_3d_r8(i,j  ,k) + array_3d_r8(i-1,j  ,k)) + & 
                                    array_3d_r8(i+1,j-1,k) + 2.0_8 * array_3d_r8(i,j-1,k) + array_3d_r8(i-1,j-1,k)) + &
                                          0.0625_4 *                                                                  &
                                   (array_3d_r4(i+1,j+1,k) + 2.0_8 * array_3d_r4(i,j+1,k) + array_3d_r4(i-1,j+1,k)  + & 
                           2.0_8 * (array_3d_r4(i+1,j  ,k) + 2.0_8 * array_3d_r4(i,j  ,k) + array_3d_r4(i-1,j  ,k)) + & 
                                    array_3d_r4(i+1,j-1,k) + 2.0_8 * array_3d_r4(i,j-1,k) + array_3d_r4(i-1,j-1,k))
                                value3dd = 0.0625_8 *                                                                    &
                                   (array_3d_r8d(i+1,j+1,k) + 2.0_8 * array_3d_r8d(i,j+1,k) + array_3d_r8d(i-1,j+1,k)  + & 
                           2.0_8 * (array_3d_r8d(i+1,j  ,k) + 2.0_8 * array_3d_r8d(i,j  ,k) + array_3d_r8d(i-1,j  ,k)) + & 
                                    array_3d_r8d(i+1,j-1,k) + 2.0_8 * array_3d_r8d(i,j-1,k) + array_3d_r8d(i-1,j-1,k)) + &
                                          0.0625_4 *                                                                     &
                                   (array_3d_r4d(i+1,j+1,k) + 2.0_8 * array_3d_r4d(i,j+1,k) + array_3d_r4d(i-1,j+1,k)  + & 
                           2.0_8 * (array_3d_r4d(i+1,j  ,k) + 2.0_8 * array_3d_r4d(i,j  ,k) + array_3d_r4d(i-1,j  ,k)) + & 
                                    array_3d_r4d(i+1,j-1,k) + 2.0_8 * array_3d_r4d(i,j-1,k) + array_3d_r4d(i-1,j-1,k))
                                value2d = 0.0625_8 *                                                            &
                                   (array_2d_r8(i+1,j+1) + 2.0_8 * array_2d_r8(i,j+1) + array_2d_r8(i-1,j+1)  + & 
                           2.0_8 * (array_2d_r8(i+1,j  ) + 2.0_8 * array_2d_r8(i,j  ) + array_2d_r8(i-1,j  )) + & 
                                    array_2d_r8(i+1,j-1) + 2.0_8 * array_2d_r8(i,j-1) + array_2d_r8(i-1,j-1)) + & 
                                          0.0625_4 *                                                            &
                                   (array_2d_r4(i+1,j+1) + 2.0_8 * array_2d_r4(i,j+1) + array_2d_r4(i-1,j+1)  + &  
                           2.0_8 * (array_2d_r4(i+1,j  ) + 2.0_8 * array_2d_r4(i,j  ) + array_2d_r4(i-1,j  )) + &  
                                    array_2d_r4(i+1,j-1) + 2.0_8 * array_2d_r4(i,j-1) + array_2d_r4(i-1,j-1))
                            else
                                value4d = 0.00390625_8 *                                                                                                               &
        (array_4d_r8(l,i+2,j+2,k) + 4.0_8 * array_4d_r8(l,i+1,j+2,k) + 6.0_8 * array_4d_r8(l,i,j+2,k) + 4.0_8 * array_4d_r8(l,i-1,j+2,k) + array_4d_r8(l,i-2,j+2,k)  + & 
4.0_8 * (array_4d_r8(l,i+2,j+1,k) + 4.0_8 * array_4d_r8(l,i+1,j+1,k) + 6.0_8 * array_4d_r8(l,i,j+1,k) + 4.0_8 * array_4d_r8(l,i-1,j+1,k) + array_4d_r8(l,i-2,j+1,k)) + & 
6.0_8 * (array_4d_r8(l,i+2,j  ,k) + 4.0_8 * array_4d_r8(l,i+1,j  ,k) + 6.0_8 * array_4d_r8(l,i,j  ,k) + 4.0_8 * array_4d_r8(l,i-1,j  ,k) + array_4d_r8(l,i-2,j  ,k)) + & 
4.0_8 * (array_4d_r8(l,i+2,j-1,k) + 4.0_8 * array_4d_r8(l,i+1,j-1,k) + 6.0_8 * array_4d_r8(l,i,j-1,k) + 4.0_8 * array_4d_r8(l,i-1,j-1,k) + array_4d_r8(l,i-2,j-1,k)) + & 
         array_4d_r8(l,i+2,j-2,k) + 4.0_8 * array_4d_r8(l,i+1,j-2,k) + 6.0_8 * array_4d_r8(l,i,j-2,k) + 4.0_8 * array_4d_r8(l,i-1,j-2,k) + array_4d_r8(l,i-2,j-2,k)) + &
       0.00390625_4 *                                                                                                                                                  &
        (array_4d_r4(l,i+2,j+2,k) + 4.0_4 * array_4d_r4(l,i+1,j+2,k) + 6.0_4 * array_4d_r4(l,i,j+2,k) + 4.0_4 * array_4d_r4(l,i-1,j+2,k) + array_4d_r4(l,i-2,j+2,k)  + & 
4.0_4 * (array_4d_r4(l,i+2,j+1,k) + 4.0_4 * array_4d_r4(l,i+1,j+1,k) + 6.0_4 * array_4d_r4(l,i,j+1,k) + 4.0_4 * array_4d_r4(l,i-1,j+1,k) + array_4d_r4(l,i-2,j+1,k)) + & 
6.0_4 * (array_4d_r4(l,i+2,j  ,k) + 4.0_4 * array_4d_r4(l,i+1,j  ,k) + 6.0_4 * array_4d_r4(l,i,j  ,k) + 4.0_4 * array_4d_r4(l,i-1,j  ,k) + array_4d_r4(l,i-2,j  ,k)) + & 
4.0_4 * (array_4d_r4(l,i+2,j-1,k) + 4.0_4 * array_4d_r4(l,i+1,j-1,k) + 6.0_4 * array_4d_r4(l,i,j-1,k) + 4.0_4 * array_4d_r4(l,i-1,j-1,k) + array_4d_r4(l,i-2,j-1,k)) + & 
         array_4d_r4(l,i+2,j-2,k) + 4.0_4 * array_4d_r4(l,i+1,j-2,k) + 6.0_4 * array_4d_r4(l,i,j-2,k) + 4.0_4 * array_4d_r4(l,i-1,j-2,k) + array_4d_r4(l,i-2,j-2,k))     
                                value4dd = 0.00390625_8 *                                                                                                              &
        (array_4d_r8d(l,i+2,j+2,k) + 4.0_8 * array_4d_r8d(l,i+1,j+2,k) + 6.0_8 * array_4d_r8d(l,i,j+2,k) + 4.0_8 * array_4d_r8d(l,i-1,j+2,k) + array_4d_r8d(l,i-2,j+2,k)  + & 
4.0_8 * (array_4d_r8d(l,i+2,j+1,k) + 4.0_8 * array_4d_r8d(l,i+1,j+1,k) + 6.0_8 * array_4d_r8d(l,i,j+1,k) + 4.0_8 * array_4d_r8d(l,i-1,j+1,k) + array_4d_r8d(l,i-2,j+1,k)) + & 
6.0_8 * (array_4d_r8d(l,i+2,j  ,k) + 4.0_8 * array_4d_r8d(l,i+1,j  ,k) + 6.0_8 * array_4d_r8d(l,i,j  ,k) + 4.0_8 * array_4d_r8d(l,i-1,j  ,k) + array_4d_r8d(l,i-2,j  ,k)) + & 
4.0_8 * (array_4d_r8d(l,i+2,j-1,k) + 4.0_8 * array_4d_r8d(l,i+1,j-1,k) + 6.0_8 * array_4d_r8d(l,i,j-1,k) + 4.0_8 * array_4d_r8d(l,i-1,j-1,k) + array_4d_r8d(l,i-2,j-1,k)) + & 
         array_4d_r8d(l,i+2,j-2,k) + 4.0_8 * array_4d_r8d(l,i+1,j-2,k) + 6.0_8 * array_4d_r8d(l,i,j-2,k) + 4.0_8 * array_4d_r8d(l,i-1,j-2,k) + array_4d_r8d(l,i-2,j-2,k)) + &
       0.00390625_4 *                                                                                                                                                  &
        (array_4d_r4d(l,i+2,j+2,k) + 4.0_4 * array_4d_r4d(l,i+1,j+2,k) + 6.0_4 * array_4d_r4d(l,i,j+2,k) + 4.0_4 * array_4d_r4d(l,i-1,j+2,k) + array_4d_r4d(l,i-2,j+2,k)  + & 
4.0_4 * (array_4d_r4d(l,i+2,j+1,k) + 4.0_4 * array_4d_r4d(l,i+1,j+1,k) + 6.0_4 * array_4d_r4d(l,i,j+1,k) + 4.0_4 * array_4d_r4d(l,i-1,j+1,k) + array_4d_r4d(l,i-2,j+1,k)) + & 
6.0_4 * (array_4d_r4d(l,i+2,j  ,k) + 4.0_4 * array_4d_r4d(l,i+1,j  ,k) + 6.0_4 * array_4d_r4d(l,i,j  ,k) + 4.0_4 * array_4d_r4d(l,i-1,j  ,k) + array_4d_r4d(l,i-2,j  ,k)) + & 
4.0_4 * (array_4d_r4d(l,i+2,j-1,k) + 4.0_4 * array_4d_r4d(l,i+1,j-1,k) + 6.0_4 * array_4d_r4d(l,i,j-1,k) + 4.0_4 * array_4d_r4d(l,i-1,j-1,k) + array_4d_r4d(l,i-2,j-1,k)) + & 
         array_4d_r4d(l,i+2,j-2,k) + 4.0_4 * array_4d_r4d(l,i+1,j-2,k) + 6.0_4 * array_4d_r4d(l,i,j-2,k) + 4.0_4 * array_4d_r4d(l,i-1,j-2,k) + array_4d_r4d(l,i-2,j-2,k))                
                                value3d = 0.00390625_8 *                                                                                                     &
        (array_3d_r8(i+2,j+2,k) + 4.0_8 * array_3d_r8(i+1,j+2,k) + 6.0_8 * array_3d_r8(i,j+2,k) + 4.0_8 * array_3d_r8(i-1,j+2,k) + array_3d_r8(i-2,j+2,k)  + & 
4.0_8 * (array_3d_r8(i+2,j+1,k) + 4.0_8 * array_3d_r8(i+1,j+1,k) + 6.0_8 * array_3d_r8(i,j+1,k) + 4.0_8 * array_3d_r8(i-1,j+1,k) + array_3d_r8(i-2,j+1,k)) + & 
6.0_8 * (array_3d_r8(i+2,j  ,k) + 4.0_8 * array_3d_r8(i+1,j  ,k) + 6.0_8 * array_3d_r8(i,j  ,k) + 4.0_8 * array_3d_r8(i-1,j  ,k) + array_3d_r8(i-2,j  ,k)) + & 
4.0_8 * (array_3d_r8(i+2,j-1,k) + 4.0_8 * array_3d_r8(i+1,j-1,k) + 6.0_8 * array_3d_r8(i,j-1,k) + 4.0_8 * array_3d_r8(i-1,j-1,k) + array_3d_r8(i-2,j-1,k)) + & 
         array_3d_r8(i+2,j-2,k) + 4.0_8 * array_3d_r8(i+1,j-2,k) + 6.0_8 * array_3d_r8(i,j-2,k) + 4.0_8 * array_3d_r8(i-1,j-2,k) + array_3d_r8(i-2,j-2,k)) + &
       0.00390625_4 *                                                                                                                                        &
        (array_3d_r4(i+2,j+2,k) + 4.0_4 * array_3d_r4(i+1,j+2,k) + 6.0_4 * array_3d_r4(i,j+2,k) + 4.0_4 * array_3d_r4(i-1,j+2,k) + array_3d_r4(i-2,j+2,k)  + & 
4.0_4 * (array_3d_r4(i+2,j+1,k) + 4.0_4 * array_3d_r4(i+1,j+1,k) + 6.0_4 * array_3d_r4(i,j+1,k) + 4.0_4 * array_3d_r4(i-1,j+1,k) + array_3d_r4(i-2,j+1,k)) + & 
6.0_4 * (array_3d_r4(i+2,j  ,k) + 4.0_4 * array_3d_r4(i+1,j  ,k) + 6.0_4 * array_3d_r4(i,j  ,k) + 4.0_4 * array_3d_r4(i-1,j  ,k) + array_3d_r4(i-2,j  ,k)) + & 
4.0_4 * (array_3d_r4(i+2,j-1,k) + 4.0_4 * array_3d_r4(i+1,j-1,k) + 6.0_4 * array_3d_r4(i,j-1,k) + 4.0_4 * array_3d_r4(i-1,j-1,k) + array_3d_r4(i-2,j-1,k)) + & 
         array_3d_r4(i+2,j-2,k) + 4.0_4 * array_3d_r4(i+1,j-2,k) + 6.0_4 * array_3d_r4(i,j-2,k) + 4.0_4 * array_3d_r4(i-1,j-2,k) + array_3d_r4(i-2,j-2,k)) 
                                value3dd = 0.00390625_8 *                                                                                                    &
        (array_3d_r8d(i+2,j+2,k) + 4.0_8 * array_3d_r8d(i+1,j+2,k) + 6.0_8 * array_3d_r8d(i,j+2,k) + 4.0_8 * array_3d_r8d(i-1,j+2,k) + array_3d_r8d(i-2,j+2,k)  + & 
4.0_8 * (array_3d_r8d(i+2,j+1,k) + 4.0_8 * array_3d_r8d(i+1,j+1,k) + 6.0_8 * array_3d_r8d(i,j+1,k) + 4.0_8 * array_3d_r8d(i-1,j+1,k) + array_3d_r8d(i-2,j+1,k)) + & 
6.0_8 * (array_3d_r8d(i+2,j  ,k) + 4.0_8 * array_3d_r8d(i+1,j  ,k) + 6.0_8 * array_3d_r8d(i,j  ,k) + 4.0_8 * array_3d_r8d(i-1,j  ,k) + array_3d_r8d(i-2,j  ,k)) + & 
4.0_8 * (array_3d_r8d(i+2,j-1,k) + 4.0_8 * array_3d_r8d(i+1,j-1,k) + 6.0_8 * array_3d_r8d(i,j-1,k) + 4.0_8 * array_3d_r8d(i-1,j-1,k) + array_3d_r8d(i-2,j-1,k)) + & 
         array_3d_r8d(i+2,j-2,k) + 4.0_8 * array_3d_r8d(i+1,j-2,k) + 6.0_8 * array_3d_r8d(i,j-2,k) + 4.0_8 * array_3d_r8d(i-1,j-2,k) + array_3d_r8d(i-2,j-2,k)) + &
       0.00390625_4 *                                                                                                                                        &
        (array_3d_r4d(i+2,j+2,k) + 4.0_4 * array_3d_r4d(i+1,j+2,k) + 6.0_4 * array_3d_r4d(i,j+2,k) + 4.0_4 * array_3d_r4d(i-1,j+2,k) + array_3d_r4d(i-2,j+2,k)  + & 
4.0_4 * (array_3d_r4d(i+2,j+1,k) + 4.0_4 * array_3d_r4d(i+1,j+1,k) + 6.0_4 * array_3d_r4d(i,j+1,k) + 4.0_4 * array_3d_r4d(i-1,j+1,k) + array_3d_r4d(i-2,j+1,k)) + & 
6.0_4 * (array_3d_r4d(i+2,j  ,k) + 4.0_4 * array_3d_r4d(i+1,j  ,k) + 6.0_4 * array_3d_r4d(i,j  ,k) + 4.0_4 * array_3d_r4d(i-1,j  ,k) + array_3d_r4d(i-2,j  ,k)) + & 
4.0_4 * (array_3d_r4d(i+2,j-1,k) + 4.0_4 * array_3d_r4d(i+1,j-1,k) + 6.0_4 * array_3d_r4d(i,j-1,k) + 4.0_4 * array_3d_r4d(i-1,j-1,k) + array_3d_r4d(i-2,j-1,k)) + & 
         array_3d_r4d(i+2,j-2,k) + 4.0_4 * array_3d_r4d(i+1,j-2,k) + 6.0_4 * array_3d_r4d(i,j-2,k) + 4.0_4 * array_3d_r4d(i-1,j-2,k) + array_3d_r4d(i-2,j-2,k))
                                value2d = 0.00390625_8 *                                                                                           &
        (array_2d_r8(i+2,j+2) + 4.0_8 * array_2d_r8(i+1,j+2) + 6.0_8 * array_2d_r8(i,j+2) + 4.0_8 * array_2d_r8(i-1,j+2) + array_2d_r8(i-2,j+2)  + & 
4.0_8 * (array_2d_r8(i+2,j+1) + 4.0_8 * array_2d_r8(i+1,j+1) + 6.0_8 * array_2d_r8(i,j+1) + 4.0_8 * array_2d_r8(i-1,j+1) + array_2d_r8(i-2,j+1)) + & 
6.0_8 * (array_2d_r8(i+2,j  ) + 4.0_8 * array_2d_r8(i+1,j  ) + 6.0_8 * array_2d_r8(i,j  ) + 4.0_8 * array_2d_r8(i-1,j  ) + array_2d_r8(i-2,j  )) + & 
4.0_8 * (array_2d_r8(i+2,j-1) + 4.0_8 * array_2d_r8(i+1,j-1) + 6.0_8 * array_2d_r8(i,j-1) + 4.0_8 * array_2d_r8(i-1,j-1) + array_2d_r8(i-2,j-1)) + & 
         array_2d_r8(i+2,j-2) + 4.0_8 * array_2d_r8(i+1,j-2) + 6.0_8 * array_2d_r8(i,j-2) + 4.0_8 * array_2d_r8(i-1,j-2) + array_2d_r8(i-2,j-2)) + &
       0.00390625_4 *                                                                                                                              &
        (array_2d_r4(i+2,j+2) + 4.0_4 * array_2d_r4(i+1,j+2) + 6.0_4 * array_2d_r4(i,j+2) + 4.0_4 * array_2d_r4(i-1,j+2) + array_2d_r4(i-2,j+2)  + & 
4.0_4 * (array_2d_r4(i+2,j+1) + 4.0_4 * array_2d_r4(i+1,j+1) + 6.0_4 * array_2d_r4(i,j+1) + 4.0_4 * array_2d_r4(i-1,j+1) + array_2d_r4(i-2,j+1)) + & 
6.0_4 * (array_2d_r4(i+2,j  ) + 4.0_4 * array_2d_r4(i+1,j  ) + 6.0_4 * array_2d_r4(i,j  ) + 4.0_4 * array_2d_r4(i-1,j  ) + array_2d_r4(i-2,j  )) + & 
4.0_4 * (array_2d_r4(i+2,j-1) + 4.0_4 * array_2d_r4(i+1,j-1) + 6.0_4 * array_2d_r4(i,j-1) + 4.0_4 * array_2d_r4(i-1,j-1) + array_2d_r4(i-2,j-1)) + & 
         array_2d_r4(i+2,j-2) + 4.0_4 * array_2d_r4(i+1,j-2) + 6.0_4 * array_2d_r4(i,j-2) + 4.0_4 * array_2d_r4(i-1,j-2) + array_2d_r4(i-2,j-2)) 
                            end if

                            s = s + abs(2.0_8 * real(i+j+k,8) - value3dd) + abs(2.0_8 * real(i+j+k,8) - value3d) + abs(2.0_8 * real(i+j,8) - value2d) + abs(2.0_8 * real(i+j+k+l,8) - value4d) + abs(2.0_8 * real(i+j+k+l,8) - value4dd)
                            end do
                        end if
                    end do
                end do
            end do

            call MPI_Allreduce(s, s_mpi, 1, MPI_REAL8, MPI_SUM, this.com, this.ierr)

            if ((s_mpi < 1.0d-10) .and. active_x .and. active_y) then
                direction = not(direction)
            end if
        end do

        ! update maximum avaliable gcx
        this.max_gcx = exch_x
        this.max_gcy = exch_y

        if (s_mpi < 1.0d-10) then
            if (this.rank == 0) write(*,*) 'mpicom exchange_halo : ok, maximum gcx, gcy:', exch_x, exch_y
        else
            if (this.rank == 0) write(*,*) 'mpicom exchange_halo : error', s_mpi
            stop
        end if

    end subroutine check_halo
    
    subroutine check_scatter(this, partition) 
        class(mpicom_t)  , intent(inout) :: this
        type(partition_t), intent(in   ) :: partition
        
        real(kind=8) :: loc_2d_r8  (this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy)
        real(kind=8) :: glob_2d_r8 (this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy)
        real(kind=8) :: glob2_2d_r8(this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy)
        real(kind=8) :: loc_2d_r4  (this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy)
        real(kind=8) :: glob_2d_r4 (this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy)
        real(kind=8) :: glob2_2d_r4(this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy)

        real(kind=8) :: loc_3d_r8  (this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,3)
        real(kind=8) :: glob_3d_r8 (this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy,3)
        real(kind=8) :: glob2_3d_r8(this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy,3)
        real(kind=8) :: loc_3d_r4  (this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,3)
        real(kind=8) :: glob_3d_r4 (this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy,3)
        real(kind=8) :: glob2_3d_r4(this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy,3)
        
        real(kind=8) :: loc_4d_r8  (5,this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,3)
        real(kind=8) :: glob_4d_r8 (5,this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy,3)
        real(kind=8) :: glob2_4d_r8(5,this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy,3)
        real(kind=8) :: loc_4d_r4  (5,this.is-this.max_gcx:this.ie+this.max_gcx,this.js-this.max_gcy:this.je+this.max_gcy,3)
        real(kind=8) :: glob_4d_r4 (5,this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy,3)
        real(kind=8) :: glob2_4d_r4(5,this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy,3)
        integer(kind=4) :: rankij(this.is_g-this.max_gcx:this.ie_g+this.max_gcx,this.js_g-this.max_gcy:this.je_g+this.max_gcy)

        integer(kind=4) :: i, j, k, l
        real   (kind=8) :: s

        call partition.get_rankij(rankij,this.max_gcx, this.max_gcy)

        if (this.rank == 0) then
            do k = 1, 3
                do j = this.js_g-this.max_gcy, this.je_g+this.max_gcy
                    do i = this.is_g-this.max_gcx, this.ie_g+this.max_gcx
                        do l = 1, 5
                            glob_2d_r8(i,j)     = real(i+j,8)
                            glob_3d_r8(i,j,k)   = real(i+j+k,8)
                            glob_4d_r8(l,i,j,k) = real(i+j+k+l,8)

                            glob_2d_r4(i,j)     = real(i+j,4)
                            glob_3d_r4(i,j,k)   = real(i+j+k,4)
                            glob_4d_r4(l,i,j,k) = real(i+j+k+l,4)
                        end do
                    end do
                end do
            end do
        end if

        loc_2d_r8   = -1.0_8
        glob2_2d_r8 = -1.0_8
        loc_3d_r8   = -1.0_8
        glob2_3d_r8 = -1.0_8
        loc_4d_r8   = -1.0_8
        glob2_4d_r8 = -1.0_8

        loc_2d_r4   = -1.0_4
        glob2_2d_r4 = -1.0_4
        loc_3d_r4   = -1.0_4
        glob2_3d_r4 = -1.0_4
        loc_4d_r4   = -1.0_4
        glob2_4d_r4 = -1.0_4

        call this.scatter(glob_2d_r8, this.max_gcx, this.max_gcy, loc_2d_r8, this.max_gcx, this.max_gcy, 'wait')
        call this.scatter(glob_3d_r8, this.max_gcx, this.max_gcy, loc_3d_r8, this.max_gcx, this.max_gcy, 'wait')
        call this.scatter(glob_4d_r8, this.max_gcx, this.max_gcy, loc_4d_r8, this.max_gcx, this.max_gcy)
        
        call this.scatter(glob_2d_r4, this.max_gcx, this.max_gcy, loc_2d_r4, this.max_gcx, this.max_gcy, 'wait')
        call this.scatter(glob_3d_r4, this.max_gcx, this.max_gcy, loc_3d_r4, this.max_gcx, this.max_gcy, 'wait')
        call this.scatter(glob_4d_r4, this.max_gcx, this.max_gcy, loc_4d_r4, this.max_gcx, this.max_gcy)
        
        call this.gather(glob2_2d_r8, this.max_gcx, this.max_gcy, loc_2d_r8, this.max_gcx, this.max_gcy, 'wait')
        call this.gather(glob2_3d_r8, this.max_gcx, this.max_gcy, loc_3d_r8, this.max_gcx, this.max_gcy, 'wait')
        call this.gather(glob2_4d_r8, this.max_gcx, this.max_gcy, loc_4d_r8, this.max_gcx, this.max_gcy)

        call this.gather(glob2_2d_r4, this.max_gcx, this.max_gcy, loc_2d_r4, this.max_gcx, this.max_gcy, 'wait')
        call this.gather(glob2_3d_r4, this.max_gcx, this.max_gcy, loc_3d_r4, this.max_gcx, this.max_gcy, 'wait')
        call this.gather(glob2_4d_r4, this.max_gcx, this.max_gcy, loc_4d_r4, this.max_gcx, this.max_gcy)

        if (this.rank == 0) then
            s = 0.0_8
            do k = 1, 3
                do j = this.js_g, this.je_g
                    do i = this.is_g, this.ie_g
                        if (rankij(i,j) > - 1) then
                            do l = 1, 5
                                s = s + abs(glob2_2d_r8(i,j) - glob_2d_r8(i,j)) + abs(glob2_3d_r8(i,j,k) - glob_3d_r8(i,j,k)) + abs(glob2_4d_r8(l,i,j,k) - glob_4d_r8(l,i,j,k)) + &
                                        abs(glob2_2d_r4(i,j) - glob_2d_r4(i,j)) + abs(glob2_3d_r4(i,j,k) - glob_3d_r4(i,j,k)) + abs(glob2_4d_r4(l,i,j,k) - glob_4d_r4(l,i,j,k))
                            end do
                        end if
                    end do
                end do
            end do
        end if

        call mpi_Bcast(s, 1, MPI_REAL8, 0, this.com, this.ierr)

        if (s < 1.0d-10) then
            if (this.rank == 0) write(*,*) 'mpicom scatter-gather: ok'
        else
            if (this.rank == 0) write(*,*) 'mpicom scatter-gather: error', s
            stop
        end if

    end subroutine check_scatter

end module mpicom_mod
