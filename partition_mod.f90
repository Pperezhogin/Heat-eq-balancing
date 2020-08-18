module partition_mod
    implicit none

    type, public :: partition_t
        private
        
        logical(kind=1) :: test_partition
        logical(kind=1) :: partition_hard
        integer(kind=4) :: iter_partition

        integer(kind=4) :: is_g, ie_g, js_g, je_g
        integer(kind=4) :: max_gcx, max_gcy

        integer(kind=4), allocatable :: is_l(:), ie_l(:), js_l(:), je_l(:) ! boundary index for each rank

        logical(kind=1), allocatable :: mask(:,:)

        integer(kind=4) :: nbx, nby, nprocs
        integer(kind=4), allocatable, dimension(:)   :: ibx_hilbert, iby_hilbert
        integer(kind=4), allocatable, dimension(:,:) :: rank_block
        logical(kind=1), allocatable, dimension(:,:) :: mask_block
        real   (kind=8), allocatable, dimension(:,:) :: work_block
        real   (kind=8), allocatable, dimension(:,:) :: work2d_block
        real   (kind=8), allocatable, dimension(:,:) :: work3d_block
        integer(kind=4), allocatable, dimension(:,:) :: is_block, ie_block
        integer(kind=4), allocatable, dimension(:,:) :: js_block, je_block

        contains

        procedure, public :: init_1block_per_core
        procedure, public :: init_hilbert                     
        procedure, public :: cleanup                   
        procedure, public :: get_gcx                   
        procedure, public :: get_gcy                   
        procedure, public :: local_idx                 
        procedure, public :: global_idx                
        procedure, public :: get_rankij                
        procedure, public :: get_local_mask            
        procedure, public :: get_blocks_number_on_rank 
        procedure, public :: get_blocks_on_rank        
        procedure, public :: write_blocks_info

    end type partition_t

    private

    contains 
    
    subroutine init_1block_per_core(this, is_g, ie_g, js_g, je_g, max_gcx, max_gcy, mask, &
        weight, weight2d, weight3d, nprocs, myrank, test_partition, iter_partition, partition_hard)
        class(partition_t), intent(inout) :: this
        integer(kind=4)   , intent(in)    :: is_g, ie_g, js_g, je_g
        integer(kind=4)   , intent(in)    :: max_gcx, max_gcy
        logical(kind=1)   , intent(in)    :: mask(is_g:ie_g,js_g:je_g)
        real   (kind=8)   , intent(in)    :: weight(is_g:ie_g,js_g:je_g)
        real   (kind=8)   , intent(in)    :: weight2d(is_g:ie_g,js_g:je_g)
        real   (kind=8)   , intent(in)    :: weight3d(is_g:ie_g,js_g:je_g)
        integer(kind=4)   , intent(in)    :: nprocs
        integer(kind=4)   , intent(in)    :: myrank
        logical(kind=1)   , intent(in)    :: test_partition
        logical(kind=1)   , intent(in)    :: partition_hard
        integer(kind=4)   , intent(in)    :: iter_partition
        
        integer(kind=4) :: nx, ny
        integer(kind=4) :: nblocks_x, nblocks_y
        integer(kind=4) :: rest_x, rest_y

        integer(kind=4) :: i, j
        integer(kind=4) :: istart, iend, jstart, jend
        integer(kind=4) :: ibx, iby, ib, irank
        integer(kind=4) :: counter

        integer(kind=4), allocatable :: block_x(:), block_y(:)
        integer(kind=4), allocatable :: nblocks(:)
        
        this.test_partition = test_partition
        this.partition_hard = partition_hard
        this.iter_partition = iter_partition
        
        this.nprocs = nprocs

        this.is_g = is_g
        this.ie_g = ie_g
        this.js_g = js_g
        this.je_g = je_g

        nx = ie_g - is_g + 1
        ny = je_g - js_g + 1

        allocate(this.mask(is_g:ie_g,js_g:je_g))
        this.mask = mask
        allocate(nblocks(min(nx,ny)/max(max_gcx,max_gcy)))
        nblocks = 0
    
        allocate(block_x(nx),block_y(ny))
        do nblocks_x = 1, min(nx,ny)/max(max_gcx,max_gcy)
            nblocks_y = nblocks_x
      
            block_x(:) = nx / nblocks_x
            block_y(:) = ny / nblocks_y
            rest_x = nx - block_x(1) * nblocks_x
            rest_y = ny - block_y(1) * nblocks_y
            if (rest_x > 0) block_x(1:rest_x) = block_x(1:rest_x) + 1
            if (rest_y > 0) block_y(1:rest_y) = block_y(1:rest_y) + 1
            
            ib = 0
            do iby = 1, nblocks_y
                    do ibx = 1, nblocks_x
                        istart = sum(block_x(1:ibx-1)) + is_g
                        iend   = sum(block_x(1:ibx  )) + is_g - 1

                        jstart = sum(block_y(1:iby-1)) + js_g
                        jend   = sum(block_y(1:iby  )) + js_g - 1

                        if (count(mask(istart:iend,jstart:jend))>0) ib = ib + 1
                    end do
            end do
            
            nblocks(nblocks_x) = ib
            
            if (ib == nprocs) exit
        end do
        
        if (count(nblocks == nprocs) == 0) then
            if (myrank == 0) then
                write(*,*) 'partition error 1: incorrect number of cores, use:'
                write(*,*) nblocks
            end if
            stop
        end if
        
        this.nbx = nblocks_x
        this.nby = nblocks_y
        this.max_gcx = block_x(nblocks_x)
        this.max_gcy = block_y(nblocks_y)
        
        if (this.max_gcx < max_gcx .or. this.max_gcy < max_gcy) then
            if (myrank == 0) write(*,*) 'partition error 2: call pperezhogin'
            stop
        end if

        allocate(this.rank_block(nblocks_x,nblocks_y))
        allocate(this.mask_block(nblocks_x,nblocks_y))
        allocate(this.work_block(nblocks_x,nblocks_y))
        allocate(this.work2d_block(nblocks_x,nblocks_y))
        allocate(this.work3d_block(nblocks_x,nblocks_y))
        allocate(this.is_block  (nblocks_x,nblocks_y))
        allocate(this.ie_block  (nblocks_x,nblocks_y))
        allocate(this.js_block  (nblocks_x,nblocks_y))
        allocate(this.je_block  (nblocks_x,nblocks_y))

        this.rank_block = -1
        this.mask_block = .false.
        
        this.work_block = 0.
        this.work2d_block = 0.
        this.work3d_block = 0.
        
        irank = 0
        do iby = 1, nblocks_y
            do ibx = 1, nblocks_x
                istart = sum(block_x(1:ibx-1)) + is_g
                iend   = sum(block_x(1:ibx  )) + is_g - 1

                jstart = sum(block_y(1:iby-1)) + js_g
                jend   = sum(block_y(1:iby  )) + js_g - 1

                this.is_block(ibx,iby) = istart
                this.ie_block(ibx,iby) = iend
                this.js_block(ibx,iby) = jstart
                this.je_block(ibx,iby) = jend

                if (count(mask(istart:iend,jstart:jend))>0) then
                    this.mask_block(ibx,iby) = .true.
                    this.work_block(ibx,iby) = sum(weight(istart:iend,jstart:jend))
                    this.work2d_block(ibx,iby) = sum(weight2d(istart:iend,jstart:jend))
                    this.work3d_block(ibx,iby) = sum(weight3d(istart:iend,jstart:jend))
                    this.rank_block(ibx,iby) = irank
                    irank = irank + 1
                end if
            end do
        end do
        
        allocate(this.is_l(0:nprocs-1),this.ie_l(0:nprocs-1),this.js_l(0:nprocs-1),this.je_l(0:nprocs-1))
        do irank = 0, nprocs-1
            call this.local_idx(this.is_l(irank),this.ie_l(irank),this.js_l(irank),this.je_l(irank),irank)
        end do       
        
        deallocate(block_x, block_y)
        deallocate(nblocks)
        
    end subroutine init_1block_per_core

    subroutine init_hilbert(this, is_g, ie_g, js_g, je_g, max_gcx, max_gcy, nblocks_x, nblocks_y, mask, &
        weight, weight2d, weight3d, nprocs, myrank, test_partition, iter_partition, partition_hard)
        class(partition_t), intent(inout) :: this
        integer(kind=4)   , intent(in)    :: is_g, ie_g, js_g, je_g
        integer(kind=4)   , intent(in)    :: max_gcx, max_gcy
        integer(kind=4)   , intent(inout) :: nblocks_x, nblocks_y
        logical(kind=1)   , intent(in)    :: mask(is_g:ie_g,js_g:je_g)
        real   (kind=8)   , intent(in)    :: weight(is_g:ie_g,js_g:je_g)
        real   (kind=8)   , intent(in)    :: weight2d(is_g:ie_g,js_g:je_g)
        real   (kind=8)   , intent(in)    :: weight3d(is_g:ie_g,js_g:je_g)
        integer(kind=4)   , intent(in)    :: nprocs
        integer(kind=4)   , intent(in)    :: myrank
        logical(kind=1)   , intent(in)    :: test_partition
        logical(kind=1)   , intent(in)    :: partition_hard
        integer(kind=4)   , intent(in)    :: iter_partition
        
        integer(kind=4) :: nx, ny
        integer(kind=4) :: hilbert_level
        integer(kind=4) :: rest_x, rest_y

        integer(kind=4) :: i, j
        integer(kind=4) :: istart, iend, jstart, jend
        integer(kind=4) :: ibx, iby, ib
        integer(kind=4) :: counter

        integer(kind=4), allocatable :: block_x(:), block_y(:)
        
        if (nprocs == 1) then
            call init_1block_per_core(this, is_g, ie_g, js_g, je_g, max_gcx, max_gcy, mask, &
                weight, weight2d, weight3d, nprocs, myrank, test_partition, iter_partition, partition_hard)
                nblocks_x = 1
                nblocks_y = 1
            return
        end if
        
        this.test_partition = test_partition
        this.partition_hard = partition_hard
        this.iter_partition = iter_partition
        
        this.nprocs = nprocs

        this.is_g = is_g
        this.ie_g = ie_g
        this.js_g = js_g
        this.je_g = je_g

        nx = ie_g - is_g + 1
        ny = je_g - js_g + 1
        
        allocate(this.mask(is_g:ie_g,js_g:je_g))
        this.mask = mask

        ! correct the number of blocks if it is not a power of 2
        nblocks_x = min(max(nblocks_x,nblocks_y), min(nx,ny)/max(max_gcx,max_gcy))
        hilbert_level = floor(log(real(nblocks_x,8))/log(2.d0) + 1.d-10)
        nblocks_x = 2 ** hilbert_level
        nblocks_y = nblocks_x
        this.nbx = nblocks_x
        this.nby = nblocks_y
        call hilbert(this.ibx_hilbert,this.iby_hilbert,hilbert_level)
        
        allocate(block_x(nblocks_x),block_y(nblocks_y))
        block_x(:) = nx / nblocks_x
        block_y(:) = ny / nblocks_y
        rest_x = nx - block_x(1) * nblocks_x
        rest_y = ny - block_y(1) * nblocks_y
        if (rest_x > 0) block_x(1:rest_x) = block_x(1:rest_x) + 1
        if (rest_y > 0) block_y(1:rest_y) = block_y(1:rest_y) + 1

        this.max_gcx = block_x(nblocks_x)
        this.max_gcy = block_y(nblocks_y)
        
        if (this.max_gcx < max_gcx .or. this.max_gcy < max_gcy) then
            if (myrank == 0) write(*,*) 'partition error 4: use less blocks'
            stop
        end if

        allocate(this.rank_block(nblocks_x,nblocks_y))
        allocate(this.mask_block(nblocks_x,nblocks_y))
        allocate(this.work_block(nblocks_x,nblocks_y))
        allocate(this.work2d_block(nblocks_x,nblocks_y))
        allocate(this.work3d_block(nblocks_x,nblocks_y))
        allocate(this.is_block  (nblocks_x,nblocks_y))
        allocate(this.ie_block  (nblocks_x,nblocks_y))
        allocate(this.js_block  (nblocks_x,nblocks_y))
        allocate(this.je_block  (nblocks_x,nblocks_y))

        this.rank_block = -1
        this.mask_block = .false.
        
        this.work_block = 0.
        this.work2d_block = 0.
        this.work3d_block = 0.
        do iby = 1, nblocks_y
            do ibx = 1, nblocks_x
                istart = sum(block_x(1:ibx-1)) + is_g
                iend   = sum(block_x(1:ibx  )) + is_g - 1

                jstart = sum(block_y(1:iby-1)) + js_g
                jend   = sum(block_y(1:iby  )) + js_g - 1

                this.is_block(ibx,iby) = istart
                this.ie_block(ibx,iby) = iend
                this.js_block(ibx,iby) = jstart
                this.je_block(ibx,iby) = jend

                if (count(mask(istart:iend,jstart:jend))>0) then
                    this.mask_block(ibx,iby) = .true.
                    this.work_block(ibx,iby) = sum(weight(istart:iend,jstart:jend))
                    this.work2d_block(ibx,iby) = sum(weight2d(istart:iend,jstart:jend))
                    this.work3d_block(ibx,iby) = sum(weight3d(istart:iend,jstart:jend))
                end if
            end do
        end do

        deallocate(block_x, block_y)
        
        call balance_hilbert(this,nprocs,myrank)

    end subroutine init_hilbert

    subroutine cleanup(this)
        class(partition_t), intent(inout) :: this
        
        if (allocated(this.mask))          deallocate(this.mask)
        if (allocated(this.ibx_hilbert))   deallocate(this.ibx_hilbert)
        if (allocated(this.iby_hilbert))   deallocate(this.iby_hilbert)
        if (allocated(this.rank_block))    deallocate(this.rank_block)
        if (allocated(this.mask_block))    deallocate(this.mask_block)
        if (allocated(this.work_block))    deallocate(this.work_block)
        if (allocated(this.work2d_block))  deallocate(this.work2d_block)
        if (allocated(this.work3d_block))  deallocate(this.work3d_block)
        if (allocated(this.is_block))      deallocate(this.is_block)
        if (allocated(this.ie_block))      deallocate(this.ie_block)
        if (allocated(this.js_block))      deallocate(this.js_block)
        if (allocated(this.je_block))      deallocate(this.je_block)

    end subroutine cleanup

    subroutine balance_hilbert(this, nprocs, myrank)
        class(partition_t), intent(inout) :: this
        integer(kind=4), intent(in) :: nprocs
        integer(kind=4), intent(in) :: myrank

        real(kind=8) :: work_mean, work_cumulate
        real(kind=8) :: work_rank(0:nprocs-1)
        integer(kind=4) :: ib, ibx, iby, irank, nb, balance_iter
        integer(kind=4) :: i, j, sum1, sum2
        logical(kind=1) :: warning_message

        work_mean = sum(this.work_block) / nprocs

        this.rank_block = -1
        work_cumulate = 0.
        irank = 0
        nb = 0
        do ib = 1,this.nbx*this.nby
            ibx = this.ibx_hilbert(ib)
            iby = this.iby_hilbert(ib)

            if (this.mask_block(ibx,iby)) then
                nb = nb + 1
                if (work_cumulate + this.work_block(ibx,iby)/2 > work_mean * (irank+1)) then
                    irank = irank + 1
                end if
                work_cumulate = work_cumulate + this.work_block(ibx,iby)
                this.rank_block(ibx,iby) = irank
            end if
        end do

        if (nb < nprocs) then
            if (myrank == 0) write(*,*) 'partition error 5: use no more than', nb, 'cores, or use more blocks'
            stop
        end if
        
        
        if (this.test_partition) write(*,*) '% imbalance before optimization of partition', &
            real(100*imbalance(this.work_block, this.rank_block, this.nbx, this.nby, this.nprocs))
        if (nprocs > 1) then
            call remove_not_connected_subdomains(this.rank_block, this.work_block, this.nbx, this.nby, nprocs)
            if (this.test_partition) write(*,*) '% imbalance after removing subdomains', &
                real(100*imbalance(this.work_block, this.rank_block, this.nbx, this.nby, this.nprocs))
            do balance_iter = 1,this.iter_partition
                if (this.test_partition) write(*,*) 'starting iteration of partition optimization', balance_iter
                if (this.partition_hard) then
                    call balance_all_ranks(this.rank_block, this.work_block, this.nbx, this.nby, nprocs)
                    if (this.test_partition)  write(*,*) '% imbalance after balancing_all_ranks', & 
                        real(100*imbalance(this.work_block, this.rank_block, this.nbx, this.nby, this.nprocs))
                else
                    call balance_min_load_rank(this.rank_block, this.work_block, this.nbx, this.nby, nprocs)
                    if (this.test_partition)  write(*,*) '% imbalance after balancing_min_load', & 
                        real(100*imbalance(this.work_block, this.rank_block, this.nbx, this.nby, this.nprocs))
                    call balance_max_load_rank(this.rank_block, this.work_block, this.nbx, this.nby, nprocs)
                    if (this.test_partition)  write(*,*) '% imbalance after balancing_min_load', & 
                        real(100*imbalance(this.work_block, this.rank_block, this.nbx, this.nby, this.nprocs))
                end if
                
                call remove_not_connected_subdomains(this.rank_block, this.work_block, this.nbx, this.nby, nprocs)
                if (this.test_partition) write(*,*) '% imbalance after removing subdomains', &
                real(100*imbalance(this.work_block, this.rank_block, this.nbx, this.nby, this.nprocs))
            end do
        end if
        
        if (this.test_partition) write(*,*) 'end of partition optimization'

        allocate(this.is_l(0:nprocs-1),this.ie_l(0:nprocs-1),this.js_l(0:nprocs-1),this.je_l(0:nprocs-1))
        do irank = 0, nprocs-1
            call this.local_idx(this.is_l(irank),this.ie_l(irank),this.js_l(irank),this.je_l(irank),irank)
        end do
        
        warning_message = .false.
        do irank = 0, nprocs-1
            if (count(this.rank_block == irank) == 0) then
                if (myrank == 0) write(*,*) 'partition error 6: use more blocks'
                stop
            end if
            
            if (count(this.rank_block == irank) < 5 .and. nprocs>1) warning_message = .true.
        end do
        
        if (myrank == 0 .and. warning_message) write(*,*) 'partition WARNING 1: try more blocks. If impossible, try partition_hard = .true. in parallel.nml'
        
    end subroutine balance_hilbert

    function get_gcx(this) result(gcx)
        class(partition_t), intent(in)  :: this
        integer(kind=4) :: gcx

        gcx = this.max_gcx

    end function get_gcx

    function get_gcy(this) result(gcy)
        class(partition_t), intent(in)  :: this
        integer(kind=4) :: gcy

        gcy = this.max_gcy

    end function get_gcy

    subroutine local_idx(this, is, ie, js, je, rank)
        class(partition_t), intent(in)  :: this
        integer(kind=4)   , intent(in)  :: rank
        integer(kind=4)   , intent(out) :: is, ie, js, je

        integer(kind=4) :: ibx, iby

        is =   huge(0)
        js =   huge(0)
        ie = - huge(0)
        je = - huge(0)

        do iby = 1, this.nby
            do ibx = 1, this.nbx
                if (this.rank_block(ibx,iby) .eq. rank) then
                    if (this.is_block(ibx,iby) < is) is = this.is_block(ibx,iby)
                    if (this.ie_block(ibx,iby) > ie) ie = this.ie_block(ibx,iby)
                    if (this.js_block(ibx,iby) < js) js = this.js_block(ibx,iby)
                    if (this.je_block(ibx,iby) > je) je = this.je_block(ibx,iby)
                end if
            end do
        end do

    end subroutine local_idx

    subroutine global_idx(this, is, ie, js, je)
        class(partition_t), intent(in)  :: this
        integer(kind=4)   , intent(out) :: is, ie, js, je

        is = this.is_g
        ie = this.ie_g
        js = this.js_g
        je = this.je_g

    end subroutine global_idx

    ! rank depending on point (i,j), -1 land
    subroutine get_rankij(this, rankij, gcx, gcy)
        class(partition_t), intent(in)  :: this
        integer(kind=4)   , intent(in)  :: gcx, gcy
        integer(kind=4)   , intent(out) :: rankij(this.is_g-gcx:this.ie_g+gcx,this.js_g-gcy:this.je_g+gcy)

        integer(kind=4) :: ibx, iby, i, j

        rankij = -1

        do iby = 1, this.nby
            do ibx = 1, this.nbx
                do j = this.js_block(ibx,iby),this.je_block(ibx,iby)
                    do i = this.is_block(ibx,iby),this.ie_block(ibx,iby)
                        rankij(i,j) = this.rank_block(ibx,iby)
                    end do
                end do
            end do
        end do

    end subroutine get_rankij

    ! rank depending on point (i,j), -1 land
    subroutine get_local_mask(this, mask, gcx, gcy, rank)
        class(partition_t), intent(in)  :: this
        integer(kind=4)   , intent(in)  :: gcx, gcy, rank
        logical(kind=1)   , intent(out) :: mask(this.is_l(rank)-gcx:this.ie_l(rank)+gcx,this.js_l(rank)-gcy:this.je_l(rank)+gcy)

        integer(kind=4) :: i, j
        integer(kind=4) :: rankij(this.is_g-gcx:this.ie_g+gcx,this.js_g-gcy:this.je_g+gcy)

        call this.get_rankij(rankij,gcx,gcy)

        mask = .false.
        do j = this.js_l(rank), this.je_l(rank)
            do i = this.is_l(rank), this.ie_l(rank)
                if (rankij(i,j) == rank .and. this.mask(i,j))  mask(i,j) = .true.
            end do
        end do
    end subroutine get_local_mask

    subroutine get_blocks_number_on_rank(this, nblocks)
        class(partition_t), intent(in)  :: this
        integer(kind=4), intent(inout)  :: nblocks(0:this.nprocs-1)
    
        integer(kind=4) :: ibx, iby, rank

        nblocks = 0
        do iby = 1, this.nby
            do ibx = 1, this.nbx
                rank = this.rank_block(ibx,iby)
                if (rank > -1) nblocks(rank) = nblocks(rank) + 1
            end do
        end do

    end subroutine get_blocks_number_on_rank

    subroutine get_blocks_on_rank(this, is, ie, js, je, max_block_number)
        class(partition_t), intent(in)  :: this
        integer(kind=4), dimension(max_block_number,0:this.nprocs-1), intent(inout) :: is, ie, js, je
        integer(kind=4), intent(in) :: max_block_number

        integer(kind=4) :: ibx, iby, rank
        integer(kind=4) :: nblocks(0:this.nprocs-1)

        nblocks = 0
        do iby = 1, this.nby
            do ibx = 1, this.nbx
                rank = this.rank_block(ibx,iby)
                if (rank > -1) then
                    nblocks(rank) = nblocks(rank) + 1
                    is(nblocks(rank),rank) = this.is_block(ibx,iby)
                    ie(nblocks(rank),rank) = this.ie_block(ibx,iby)
                    js(nblocks(rank),rank) = this.js_block(ibx,iby)
                    je(nblocks(rank),rank) = this.je_block(ibx,iby)
                end if
            end do
        end do

    end subroutine get_blocks_on_rank

    subroutine write_blocks_info(this)
        class(partition_t), intent(in)  :: this

        integer(kind=4) :: irank, i, j, ibx, iby
        integer(kind=4) :: wet_points(0:this.nprocs-1), all_points(0:this.nprocs-1)
        integer(kind=4) :: rankij(this.is_g:this.ie_g,this.js_g:this.je_g)
        integer(kind=4) :: blocks_on_rank(0:this.nprocs-1)
        real   (kind=8) :: work_procij (this.is_g:this.ie_g,this.js_g:this.je_g)
        real   (kind=8) :: work_blockij(this.is_g:this.ie_g,this.js_g:this.je_g)
        real   (kind=8) :: work_proc(0:this.nprocs-1)
        character(20)   :: numproc

        write(numproc,*) this.nprocs
        numproc = adjustl(numproc)

        call this.get_rankij(rankij,0,0)
        
        do irank = 0, this.nprocs-1
            wet_points(irank) = count(rankij == irank .and. this.mask)
            all_points(irank) = (this.ie_l(irank) - this.is_l(irank) + 1) * (this.je_l(irank) - this.js_l(irank) + 1)
        end do
        
        open(20,file = 'partition_info.out')

        write(20,*) 'partition on', this.nbx, 'x', this.nby, 'blocks'
        write(20,*) 'the number of cores', this.nprocs    
        write(20,*) 'the number of nonzero blocks', count(this.mask_block)
        write(20,*) 'minimum ratio of wet points / all points, %', int(100*minval(real(wet_points) / real(all_points))), 'on proc', minloc(real(wet_points) / real(all_points))-1
        write(20,*) 'maximum ratio of wet points / all points, %', int(100*maxval(real(wet_points) / real(all_points))), 'on proc', maxloc(real(wet_points) / real(all_points))-1
        write(20,*) 'mean    ratio of wet points / all points, %', int(SUM(real(wet_points) / real(all_points)) / this.nprocs * 100)
        write(20,*) 'work imbalance  :', int(100*imbalance(this.work_block, this.rank_block, this.nbx, this.nby, this.nprocs)),'%'
        write(20,*) 'work2d imbalance:', int(100*imbalance(this.work2d_block, this.rank_block, this.nbx, this.nby, this.nprocs)),'%'
        write(20,*) 'work3d imbalance:', int(100*imbalance(this.work3d_block, this.rank_block, this.nbx, this.nby, this.nprocs)),'%'
        
        do irank = 0, this.nprocs-1
            blocks_on_rank(irank) = count(this.rank_block == irank)
        end do
        
        write(20,*) 'min blocks:', minval(blocks_on_rank)
        write(20,*) 'max blocks:', maxval(blocks_on_rank)

        write(20,*) ''

        write(20,*) 'number of blocks on each rank'

        do irank = 0, this.nprocs-1
            write(20,*) 'rank:', irank, 'blocks:', blocks_on_rank(irank)
        end do

        close(20)

        open(30, file = 'rank_gnuplot.out')

        do j = this.js_g,this.je_g
            do i = this.is_g,this.ie_g
                write(30,*) i, j, rankij(i,j)
            end do
            write(30,*) ''
        end do
        close(30)
        
    end subroutine write_blocks_info

    subroutine hilbert(i, j, level)
        integer(kind=4), intent(in ) :: level
        integer(kind=4), allocatable, intent(inout) :: i(:), j(:) 

        integer(kind=4) :: ilevel, length
        integer(kind=4) :: ls1, le1, ls2, le2, ls3, le3, ls4, le4, idx
        integer(kind=4), dimension(4**level,2) :: A, B, C, D, AA, BB, CC, DD
        integer(kind=4), dimension(1      , 2) :: north, east, south, west

        north(1,1) =  0
        north(1,2) =  1
        east (1,1) =  1
        east (1,2) =  0
        south(1,1) =  0
        south(1,2) = -1
        west (1,1) = -1
        west (1,2) =  0

        A = 0
        B = 0
        C = 0
        D = 0
        length = 0

        do ilevel = 1, level
            ls1 = 1
            le1 = length

            ls2 =     length + 2
            le2 = 2 * length + 1

            ls3 = 2 * length + 3
            le3 = 3 * length + 2

            ls4 = 3 * length + 4
            le4 = 4 * length + 3

            length = le4

            AA(ls1:le1,:) = B(ls1:le1,:)
            AA(ls2:le2,:) = A(ls1:le1,:)
            AA(ls3:le3,:) = A(ls1:le1,:)
            AA(ls4:le4,:) = C(ls1:le1,:)
            
            AA(le1+1,:) = north(1,:)
            AA(le2+1,:) = east (1,:)
            AA(le3+1,:) = south(1,:)


            BB(ls1:le1,:) = A(ls1:le1,:)
            BB(ls2:le2,:) = B(ls1:le1,:)
            BB(ls3:le3,:) = B(ls1:le1,:)
            BB(ls4:le4,:) = D(ls1:le1,:)
            
            BB(le1+1,:) = east (1,:)
            BB(le2+1,:) = north(1,:)
            BB(le3+1,:) = west (1,:)            

            
            CC(ls1:le1,:) = D(ls1:le1,:)
            CC(ls2:le2,:) = C(ls1:le1,:)
            CC(ls3:le3,:) = C(ls1:le1,:)
            CC(ls4:le4,:) = A(ls1:le1,:)
            
            CC(le1+1,:) = west (1,:)
            CC(le2+1,:) = south(1,:)
            CC(le3+1,:) = east (1,:)


            DD(ls1:le1,:) = C(ls1:le1,:)
            DD(ls2:le2,:) = D(ls1:le1,:)
            DD(ls3:le3,:) = D(ls1:le1,:)
            DD(ls4:le4,:) = B(ls1:le1,:)
            
            DD(le1+1,:) = south(1,:)
            DD(le2+1,:) = west (1,:)
            DD(le3+1,:) = north(1,:)

            A = AA
            B = BB
            C = CC
            D = DD
        end do

        length = length + 1
        allocate(i(length), j(length))

        i(1) = 1
        j(1) = 1
        do idx = 1, length-1
            i(idx+1) = i(idx) + A(idx,1)
            j(idx+1) = j(idx) + A(idx,2)
        end do
    end subroutine hilbert
    
    subroutine balance_max_load_rank(rank_block, work_block, nbx, nby, nprocs)
        ! find max_load_rank. then find min-load neighbour and try brute force
        integer(kind=4), intent(inout) :: rank_block(nbx,nby)
        real(kind=8)   , intent(in)    :: work_block(nbx,nby)
        integer(kind=4), intent(in)    :: nbx, nby, nprocs

        real(kind=8)    :: work_rank(0:nprocs-1), w_block(nbx*nby)
        integer(kind=4) :: ibx_arr(nbx*nby), iby_arr(nbx*nby)
        integer(kind=4) :: rank_of_max_work, rank_of_min_work
        integer(kind=4) :: nblocks

        call compute_work_distribution(work_rank,work_block,rank_block,nbx,nby,nprocs)
        rank_of_max_work = max_location(work_rank)-1

        call find_min_load_neighbour_rank(rank_of_min_work,nblocks,ibx_arr,iby_arr, &
            rank_block,work_rank,rank_of_max_work,nbx,nby,nprocs)
            
        nblocks = min(nblocks,20)
        
        call try_brute_force(rank_of_max_work,rank_of_min_work,nblocks,ibx_arr,iby_arr,rank_block, &
            work_rank,work_block,nbx,nby,nprocs)

    end subroutine balance_max_load_rank

    subroutine balance_min_load_rank(rank_block, work_block, nbx, nby, nprocs)
        ! find min_load_rank. then find max-load neighbour and try brute force
        integer(kind=4), intent(inout) :: rank_block(nbx,nby)
        real(kind=8)   , intent(in)    :: work_block(nbx,nby)
        integer(kind=4), intent(in)    :: nbx, nby, nprocs

        real(kind=8)    :: work_rank(0:nprocs-1), w_block(nbx*nby)
        integer(kind=4) :: ibx_arr(nbx*nby), iby_arr(nbx*nby)
        integer(kind=4) :: rank_of_max_work, rank_of_min_work
        integer(kind=4) :: nblocks

        call compute_work_distribution(work_rank,work_block,rank_block,nbx,nby,nprocs)
        rank_of_min_work = min_location(work_rank)-1

        call find_max_load_neighbour_rank(rank_of_max_work,nblocks,ibx_arr,iby_arr, &
            rank_block,work_rank,rank_of_min_work,nbx,nby,nprocs)
            
        nblocks = min(nblocks,20)
        
        call try_brute_force(rank_of_max_work,rank_of_min_work,nblocks,ibx_arr,iby_arr,rank_block, &
            work_rank,work_block,nbx,nby,nprocs)

    end subroutine balance_min_load_rank
    
    subroutine balance_all_ranks(rank_block, work_block, nbx, nby, nprocs)
        ! find min_load_rank. then find max-load neighbour and try brute force
        integer(kind=4), intent(inout) :: rank_block(nbx,nby)
        real(kind=8)   , intent(in)    :: work_block(nbx,nby)
        integer(kind=4), intent(in)    :: nbx, nby, nprocs

        real(kind=8)    :: work_rank(0:nprocs-1), w_block(nbx*nby)
        real(kind=8)    :: mean_work
        integer(kind=4) :: ibx_arr(nbx*nby), iby_arr(nbx*nby)
        integer(kind=4) :: rank_of_interest, rank_of_max_work, rank_of_min_work
        integer(kind=4) :: nblocks
        integer(kind=4) :: lol

        do rank_of_interest = 0, nprocs-1
            call compute_work_distribution(work_rank,work_block,rank_block,nbx,nby,nprocs)
            mean_work = SUM(work_rank) / nprocs
            
            if (work_rank(rank_of_interest) < mean_work) then
                call find_max_load_neighbour_rank(rank_of_max_work,nblocks,ibx_arr,iby_arr, &
                   rank_block,work_rank,rank_of_interest,nbx,nby,nprocs)
                rank_of_min_work = rank_of_interest
            else
                call find_min_load_neighbour_rank(rank_of_min_work,nblocks,ibx_arr,iby_arr, &
                    rank_block,work_rank,rank_of_interest,nbx,nby,nprocs)
                rank_of_max_work = rank_of_interest
            end if
            
            nblocks = min(nblocks,20)
            
            call try_brute_force(rank_of_max_work,rank_of_min_work,nblocks,ibx_arr,iby_arr,rank_block, &
                work_rank,work_block,nbx,nby,nprocs)           
        end do
            
    end subroutine balance_all_ranks
    
    subroutine try_brute_force(rank_of_max_work,rank_of_min_work,nblocks,ibx_arr,iby_arr,&
        rank_block,work_rank,work_block,nbx,nby,nprocs)
        ! blocks initially located on rank_of_max_work
        integer(kind=4), intent(in)    :: rank_of_max_work, rank_of_min_work
        integer(kind=4), intent(in)    :: nblocks
        integer(kind=4), intent(in)    :: ibx_arr(:), iby_arr(:)
        integer(kind=4), intent(inout) :: rank_block(nbx,nby)
        real(kind=8)   , intent(in)    :: work_rank(0:nprocs-1)
        real(kind=8)   , intent(in)    :: work_block(nbx,nby)
        integer(kind=4), intent(in)    :: nbx, nby, nprocs
        
        real(kind=8)    :: w_block(nbx*nby)
        integer(kind=4) :: best_irank(nbx*nby), irank(nbx*nby)
        real(kind=8)    :: w0(2), w0_no_blocks(2), w1(2)
        
        integer(kind=4) :: iblock, num_cases, icase, icase_1
        real(kind=8)    :: max_work
        real(kind=8)    :: imbalance_before, imbalance_after
        
        if (nblocks == 0) return
        
        imbalance_before = imbalance(work_block, rank_block, nbx, nby, nprocs)
        
        w0(1) = work_rank(rank_of_max_work); w0(2) = work_rank(rank_of_min_work)
        w0_no_blocks = w0
        
        do iblock = 1,nblocks
           w_block(iblock) = work_block(ibx_arr(iblock),iby_arr(iblock))
           w0_no_blocks(1) = w0_no_blocks(1) - w_block(iblock) ! blocks initially located on rank_of_max_work
        end do
        
        num_cases = 2**nblocks;
        max_work = maxval(w0);
        if(w0(2) > w0(1)) write(*,*) 'partition error 7: call pperezhogin' 
        best_irank = 0 ! 0 "rank" corresponds to rank_of_max_work, 1 to rank_of_min_work
        ! i.e. we set initial configureation
        
        do icase = 0,num_cases-1
            icase_1 = icase
            do iblock = 1,nblocks
               irank(iblock) = mod(icase_1,2)
               icase_1 = icase_1/2
            end do
            
            w1 = w0_no_blocks;
            do iblock = 1,nblocks
               if (irank(iblock) == 0) w1(1) = w1(1) + w_block(iblock)
               if (irank(iblock) == 1) w1(2) = w1(2) + w_block(iblock)
            end do
            
            if (maxval(w1) < max_work) then
               best_irank = irank; 
               max_work = maxval(w1);
            end if
        end do
        
        do iblock = 1,nblocks
           if (best_irank(iblock) == 0) rank_block(ibx_arr(iblock),iby_arr(iblock)) =  rank_of_max_work;
           if (best_irank(iblock) == 1) rank_block(ibx_arr(iblock),iby_arr(iblock)) =  rank_of_min_work;
        end do
        
        imbalance_after = imbalance(work_block, rank_block, nbx, nby, nprocs)
        
        if (imbalance_after > imbalance_before + 1e-10) write(*,*) 'wrong balance_min_load_rank', imbalance_before,'->',imbalance_after
        
    end subroutine try_brute_force

    subroutine remove_not_connected_subdomains(rank_block, work_block, nbx, nby, nprocs)
        integer(kind=4), intent(inout) :: rank_block(nbx,nby)
        real(kind=8)   , intent(in)    :: work_block(nbx,nby)
        integer(kind=4), intent(in)    :: nbx, nby, nprocs

        integer(kind=4) :: irank

        do irank = 0,nprocs-1
            call remove_subdomain(rank_block,work_block,nbx,nby,irank,nprocs)
        end do

    end subroutine remove_not_connected_subdomains
    
    subroutine remove_subdomain(rank_block, work_block, nbx, nby, rank_of_interest, nprocs)
        ! algorithm of clustering: The shortest non-closed path algorithm
        integer(kind=4), intent(inout) :: rank_block(nbx,nby)
        real(kind=8)   , intent(in)    :: work_block(nbx,nby)
        integer(kind=4), intent(in)    :: nbx, nby, rank_of_interest, nprocs

        integer(kind=4) :: ibx_on_rank(nbx*nby), iby_on_rank(nbx*nby)
        integer(kind=4) :: graph(nbx*nby,4), connections(nbx*nby), distance(nbx*nby,4)
        integer(kind=4) :: color(nbx*nby), color_2d(nbx,nby)
        integer(kind=4) :: rank_new(nbx,nby), ranks(4)
        real   (kind=8) :: work_rank(0:nprocs-1)
        logical(kind=1) :: mark(nbx*nby), mask_of_remove(nbx,nby)
        real   (kind=8), allocatable :: work_color(:)

        integer(kind=4) :: ibx, iby, nb, ib, ib_marked, ib_non_marked, ib_new, ib_old
        integer(kind=4) :: dist, dist_best, icolor, color_of_max_work, nranks
        logical(kind=1) :: new_connection

        ! form 1d array of blocks belonging to rank_of_interest
        nb = 0;
        do iby = 1, nby
            do ibx = 1, nbx
                if (rank_block(ibx,iby) == rank_of_interest) then
                    nb = nb + 1
                    ibx_on_rank(nb) = ibx
                    iby_on_rank(nb) = iby
                end if
            end do
        end do

        ! find graph of closest points
        graph = 0
        connections = 0
        distance = 0
        mark = .false.
        mark(1) = .true.

        new_connection = .true.
        do while (new_connection)
            dist_best = nbx+nby
            new_connection = .false.
            do ib_marked = 1,nb
                do ib_non_marked = 1,nb
                    if (.not. mark(ib_non_marked) .and. mark(ib_marked)) then
                        dist = abs(ibx_on_rank(ib_non_marked) - ibx_on_rank(ib_marked)) + abs(iby_on_rank(ib_non_marked) - iby_on_rank(ib_marked))
                        if (dist < dist_best) then
                            ib_new = ib_non_marked
                            ib_old = ib_marked
                            new_connection = .true.
                            dist_best = dist
                        end if
                    end if
                end do
            end do

            if (new_connection) then
                mark(ib_new) = .true.
                ! remove distance bigger than 1
                if (dist_best == 1) then
                    connections(ib_new) = connections(ib_new) + 1
                    connections(ib_old) = connections(ib_old) + 1

                    graph(ib_new,connections(ib_new)) = ib_old
                    graph(ib_old,connections(ib_old)) = ib_new
                    distance(ib_new,connections(ib_new)) = dist_best
                    distance(ib_old,connections(ib_old)) = dist_best
                end if
            end if
        end do

        ! paint all subgrahs
        color  = -1
        icolor =  0
        do ib = 1,nb
            if (color(ib) == -1) then
                color(ib) = icolor
                call extend_color(graph,connections,color,ib)
                icolor = icolor + 1
            end if
        end do

        ! transform color to 2d array
        color_2d = -1
        nb = 0;
        do iby = 1,nby
           do ibx = 1,nbx
              if (rank_block(ibx,iby) == rank_of_interest) then
                nb = nb+1;
                color_2d(ibx,iby) = color(nb);
              end if
           end do
        end do

        allocate(work_color(0:icolor-1))

        call compute_work_distribution(work_color,work_block,color_2d,nbx,nby,icolor)
        color_of_max_work = max_location(work_color)-1
        deallocate(work_color)

        ! remove subdomains which have less work than subdomain of maximum work
        mask_of_remove = .false.
        where(color_2d > -1 .and. color_2d .ne. color_of_max_work)
            mask_of_remove = .true.
        end where

        rank_new = rank_block
        do while (count(mask_of_remove) > 0)
            rank_block = rank_new
            call compute_work_distribution(work_rank,work_block,rank_block,nbx,nby,nprocs)
            do iby = 1, nby
                do ibx = 1, nbx
                    if (.not. mask_of_remove(ibx,iby)) cycle

                    call find_neighbour_ranks(nranks,ranks,rank_block,nbx,nby,ibx,iby)
                    if (nranks > 0) then
                        rank_new(ibx,iby) = ranks(min_location(work_rank(ranks(1:nranks))))
                        mask_of_remove(ibx,iby) = .false.
                    end if
                end do
            end do
        end do      
        rank_block = rank_new

    end subroutine remove_subdomain

    recursive subroutine extend_color(graph, connections, color, ib)
        integer(kind=4), intent(in)    :: graph(:,:), connections(:)
        integer(kind=4), intent(in)    :: ib
        integer(kind=4), intent(inout) :: color(:)

        integer(kind=4) :: ic, ib_new

        do ic = 1,connections(ib)
            ib_new = graph(ib,ic)
            if (color(ib_new) == -1) then
               color(ib_new) = color(ib)
               call extend_color(graph,connections,color,ib_new)
            end if
        end do

    end subroutine extend_color

    subroutine find_min_load_neighbour_rank(min_load_rank,nblocks,ibx_arr,iby_arr, &
        rank_block,work_rank,rank_of_interest,nbx,nby,nprocs)
        ! for rank_of_interest find neighbouring rank with minimum work (min_load_rank)
        ! returns list of blocks located on rank_of_interest adjoining with min_load_rank
        integer(kind=4), intent(out):: min_load_rank, nblocks
        integer(kind=4), intent(out):: ibx_arr(:), iby_arr(:)
        integer(kind=4), intent(in) :: rank_block(nbx,nby)
        real(kind=8)   , intent(in) :: work_rank(0:nprocs-1)
        integer(kind=4), intent(in) :: nbx, nby, rank_of_interest, nprocs

        real(kind=8)    :: min_load
        integer(kind=4) :: ibx, iby, nranks, ranks(4)

        min_load = work_rank(rank_of_interest) ! neighbouring rank should have smaller work
        min_load_rank = -1;
        
        do iby = 1,nby
            do ibx = 1,nbx
                if (rank_block(ibx,iby) .ne. rank_of_interest) cycle
                
                call find_neighbour_ranks(nranks,ranks,rank_block,nbx,nby,ibx,iby)
                if (nranks>0) then
                    if (minval(work_rank(ranks(1:nranks))) < min_load) then
                        min_load_rank = ranks(min_location(work_rank(ranks(1:nranks))))
                        min_load      =             minval(work_rank(ranks(1:nranks)))
                    end if
                end if
            end do
        end do
        
        nblocks = 0
        if (min_load_rank < 0) return
        do iby = 1,nby
            do ibx = 1,nbx
                if (rank_block(ibx,iby) .ne. rank_of_interest) cycle
                
                call find_neighbour_ranks(nranks,ranks,rank_block,nbx,nby,ibx,iby)
                if (nranks>0) then
                    if (count(ranks(1:nranks) == min_load_rank) > 0) then
                        if (count_my_ranks_in_neighbourhood(rank_of_interest, rank_block, nbx, nby, ibx, iby) <= 6) then
                            nblocks = nblocks + 1
                            ibx_arr(nblocks) = ibx
                            iby_arr(nblocks) = iby
                        end if
                    end if
                end if
            end do
        end do
    end subroutine find_min_load_neighbour_rank

    subroutine find_max_load_neighbour_rank(max_load_rank,nblocks,ibx_arr,iby_arr, &
        rank_block,work_rank,rank_of_interest,nbx,nby,nprocs)
        ! for rank_of_interest find neighbouring rank with maximum work (max_load_rank)
        ! returns list of blocks located on max_load_rank adjoining with rank_of_interest
        integer(kind=4), intent(out):: max_load_rank, nblocks
        integer(kind=4), intent(out):: ibx_arr(:), iby_arr(:)
        integer(kind=4), intent(in) :: rank_block(nbx,nby)
        real(kind=8)   , intent(in) :: work_rank(0:nprocs-1)
        integer(kind=4), intent(in) :: nbx, nby, rank_of_interest, nprocs

        real(kind=8)    :: max_load
        integer(kind=4) :: ibx, iby, nranks, ranks(4)

        max_load = work_rank(rank_of_interest) ! neighbouring rank should have bigger work
        max_load_rank = -1

        do iby = 1,nby
            do ibx = 1,nbx
                if (rank_block(ibx,iby) .ne. rank_of_interest) cycle
                
                call find_neighbour_ranks(nranks,ranks,rank_block,nbx,nby,ibx,iby)
                if (nranks>0) then
                    if (maxval(work_rank(ranks(1:nranks))) > max_load) then
                        max_load_rank = ranks(max_location(work_rank(ranks(1:nranks))))
                        max_load      =             maxval(work_rank(ranks(1:nranks)))
                    end if
                end if
            end do
        end do

        nblocks = 0
        if (max_load_rank < 0) return        
        do iby = 1,nby
            do ibx = 1,nbx
                if (rank_block(ibx,iby) .ne. max_load_rank) cycle
                
                call find_neighbour_ranks(nranks,ranks,rank_block,nbx,nby,ibx,iby)
                if (nranks>0) then                
                    if (count(ranks(1:nranks) == rank_of_interest) > 0) then
                        if (count_my_ranks_in_neighbourhood(max_load_rank, rank_block, nbx, nby, ibx, iby) <= 6) then
                            nblocks = nblocks + 1
                            ibx_arr(nblocks) = ibx
                            iby_arr(nblocks) = iby
                        end if
                    end if
                end if
            end do
        end do
    end subroutine find_max_load_neighbour_rank

    subroutine find_neighbour_ranks(nranks, ranks, rank_block, nbx, nby, ibx, iby)
        integer(kind=4), intent(out) :: nranks
        integer(kind=4), intent(out) :: ranks(4)
        integer(kind=4), intent(in)  :: rank_block(nbx,nby)
        integer(kind=4), intent(in)  :: nbx, nby, ibx, iby

        nranks = 0
        ranks = -1
        
        if (ibx<nbx) then
            if ((rank_block(ibx+1,iby) .ne. rank_block(ibx,iby)) .and. (rank_block(ibx+1,iby) > -1)) then
                nranks = nranks + 1
                ranks(nranks) = rank_block(ibx+1,iby) 
            end if
        end if

        if (ibx>1) then
            if ((rank_block(ibx-1,iby) .ne. rank_block(ibx,iby)) .and. (rank_block(ibx-1,iby) > -1)) then
                nranks = nranks + 1
                ranks(nranks) = rank_block(ibx-1,iby) 
            end if
        end if

        if (iby<nby) then
            if ((rank_block(ibx,iby+1) .ne. rank_block(ibx,iby)) .and. (rank_block(ibx,iby+1) > -1)) then
                nranks = nranks + 1
                ranks(nranks) = rank_block(ibx,iby+1) 
            end if
        end if

        if (iby>1) then
            if ((rank_block(ibx,iby-1) .ne. rank_block(ibx,iby)) .and. (rank_block(ibx,iby-1) > -1)) then
                nranks = nranks + 1
                ranks(nranks) = rank_block(ibx,iby-1) 
            end if
        end if

    end subroutine find_neighbour_ranks
    
    integer(kind=4) function count_my_ranks_in_neighbourhood(my_rank, rank_block, nbx, nby, ibx, iby)
        integer(kind=4), intent(in)  :: my_rank
        integer(kind=4), intent(in)  :: rank_block(nbx,nby)
        integer(kind=4), intent(in)  :: nbx, nby, ibx, iby
        
        integer(kind=4) :: count_rank

        count_rank = 0
        
        if (rank_block(ibx,iby) == my_rank) count_rank = count_rank + 1
        
        if (ibx<nbx) then
            if (rank_block(ibx+1,iby) == my_rank) count_rank = count_rank + 1
        end if

        if (ibx>1) then
            if (rank_block(ibx-1,iby) == my_rank) count_rank = count_rank + 1
        end if
        
        if (iby<nby) then
            if (rank_block(ibx,iby+1) == my_rank) count_rank = count_rank + 1
        end if
        
        if (iby>1) then
            if (rank_block(ibx,iby-1) == my_rank) count_rank = count_rank + 1
        end if
        
        if (ibx<nbx .and. iby<nby) then
            if (rank_block(ibx+1,iby+1) == my_rank) count_rank = count_rank + 1
        end if
        
        if (ibx<nbx .and. iby>1) then
            if (rank_block(ibx+1,iby-1) == my_rank) count_rank = count_rank + 1
        end if
        
        if (ibx>1 .and. iby<nby) then
            if (rank_block(ibx-1,iby+1) == my_rank) count_rank = count_rank + 1
        end if
        
        if (ibx>1 .and. iby>1) then
            if (rank_block(ibx-1,iby-1) == my_rank) count_rank = count_rank + 1
        end if

        count_my_ranks_in_neighbourhood = count_rank
    end function count_my_ranks_in_neighbourhood

    subroutine compute_work_distribution(work_rank, work_block, rank_block, nbx, nby, nprocs)
        real(kind=8)   , intent(out) :: work_rank(0:nprocs-1)
        real(kind=8)   , intent(in)  :: work_block(1:nbx,1:nby)
        integer(kind=4), intent(in)  :: rank_block(1:nbx,1:nby)
        integer(kind=4), intent(in)  :: nbx, nby, nprocs
        
        integer(kind=4) :: ibx, iby, irank

        work_rank = 0.
        do iby = 1,nby
            do ibx = 1,nbx
                if (rank_block(ibx,iby) > -1) then
                    irank = rank_block(ibx,iby)
                    work_rank(irank) = work_rank(irank) + work_block(ibx,iby)
                end if
            end do
        end do
    end subroutine compute_work_distribution

    real(kind=8) function imbalance(work_block, rank_block, nbx, nby, nprocs)
        real(kind=8)   , intent(in)  :: work_block(1:nbx,1:nby)
        integer(kind=4), intent(in)  :: rank_block(1:nbx,1:nby)
        integer(kind=4), intent(in)  :: nbx, nby, nprocs
        
        real(kind=8)    :: work_rank(0:nprocs-1)
        integer(kind=4) :: ibx, iby, irank

        work_rank = 0.
        do iby = 1,nby
            do ibx = 1,nbx
                if (rank_block(ibx,iby) > -1) then
                    irank = rank_block(ibx,iby)
                    work_rank(irank) = work_rank(irank) + work_block(ibx,iby)
                end if
            end do
        end do

        imbalance = (maxval(work_rank) - sum(work_rank)/size(work_rank)) / (sum(work_rank)/size(work_rank))
    end function imbalance

    integer(kind=4) function min_location(array)
        real(kind=8) :: array(:)
        
        integer(kind=4) :: i
        real(kind=8) :: minimum

        minimum = huge(1.0_8)
        do i = 1, size(array)
            if (array(i) < minimum) then
                minimum = array(i)
                min_location = i
            end if
        end do
    end function min_location

    integer(kind=4) function max_location(array)
        real(kind=8) :: array(:)
        
        integer(kind=4) :: i
        real(kind=8) :: maximum

        maximum = - huge(1.0_8)
        do i = 1, size(array)
            if (array(i) > maximum) then
                maximum= array(i)
                max_location = i
            end if
        end do
    end function max_location

end module partition_mod
