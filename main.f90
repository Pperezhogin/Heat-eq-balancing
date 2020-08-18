program heat_eq_parallel
    use MPI
    use partition_mod
    use mpicom_mod
    
    integer, parameter :: nx = 500, ny = 500, nz = 39
    
    ! exchanges data structures
    type(partition_t) :: partition
    type(mpicom_t)    :: mpicom
    
    ! global fields
    integer(kind=4) :: depth_g(0:nx+1,0:ny+1)
    real   (kind=8) :: T_g(0:nx+1,0:ny+1,1:nz)
    
    ! local fields
    logical(kind=1), allocatable, dimension(:,:)   :: mask
    integer(kind=4), allocatable, dimension(:,:)   :: depth
    integer(kind=4), allocatable, dimension(:,:)   :: exch ! depth to exchange
    real   (kind=8), allocatable, dimension(:,:,:) :: T, Tn
    
    integer(kind=4) :: i, j, k
    integer(kind=4) :: is, ie, js, je, ke

    call read_depth()
    call init_parallel() ! see params in parallel.nml
    
    call set_initial_condition()    
    call time_loop()
    call get_results()
        
    call MPI_Finalize(mpicom.ierr)
    contains
    
    subroutine time_loop
        integer(kind=4) :: nt
        real   (kind=8) :: time_1, time_2, time_3D
        
        time_1 = MPI_Wtime()
        time_3D = 0._8
        do nt = 1, 100000
            time_2 = MPI_Wtime()
            do j = js, je
                do i = is, ie
                    do k = 1,depth(i,j) ! only wet points corresponding to given CPU are nonzeros 
                        Tn(i,j,k) = 0.25_8 * (T(i+1,j,k) + T(i-1,j,k) + T(i,j+1,k) + T(i,j-1,k))
                    end do
                end do
            end do
            
            !T(:,:,1:ke) = Tn(:,:,1:ke) ! may be expensive (includes land points)    
            do j = js, je
                do i = is, ie
                    do k = 1,depth(i,j)
                        T(i,j,k) = Tn(i,j,k)
                    end do
                end do
            end do
            
            time_3D = time_3D + (MPI_Wtime() - time_2)
            call mpicom.exchange_halo(T,exch,1,1,1,1)
            
            ! solution and runtime diagnostics
            if (mod(nt,1000) == 0) call write_error(nt, time_1)
            if (nt == 50000) then
                call mpicom.print_diagnostics(partition, time_3D, MPI_Wtime()-time_1, "50000 time steps chronometry")
                mpicom.use_mpi_barriers = .false.
                time_1 = MPI_Wtime()
                if (mpicom.rank == 0) write(*,*) ''
            end if
        end do
    end subroutine
    
    subroutine set_initial_condition()
        ! Temperature = 1 on boundary, 0 inside the domain. Exact solution = 1.
        if (mpicom.rank == 0) then
            T_g = 1._8
            do j = 1, ny
                do i = 1, nx
                    do k = 1,depth_g(i,j)
                        T_g(i,j,k) = 0._8
                    end do
                end do
            end do
        end if
        
        allocate(T (is-1:ie+1,js-1:je+1,nz))
        allocate(Tn(is-1:ie+1,js-1:je+1,nz))
        
        call mpicom.scatter(T_g,1,1,T,1,1) ! scatter with b.c.
        Tn = T        
    end subroutine
    
    subroutine get_results()
        call mpicom.gather(T_g,1,1,T,1,1)
        
        if (mpicom.rank == 0) then
            write(*,*) 'Error on root:', maxval(abs(T_g-1._8))
            
            open(1, file = 'solution_gnuplot.out')
            do j = 1,ny
                do i = 1,nx
                    write(1,*) i, j, T_g(i,j,1)
                end do
                write(1,*) ''
            end do
            close(1)
        end if
    end subroutine
    
    subroutine write_error(nt, time_1)
        integer(kind=4) :: nt
        real(kind=8) :: time_1
        
        real(kind=8) :: err_l, err
        
        err_l = 0._8
        do j = js, je
            do i = is, ie
                do k = 1,depth(i,j)
                    err_l = max(abs(T(i,j,k)-1._8),err_l)
                end do
            end do
        end do
        
        call MPI_Allreduce(err_l,err,1,MPI_REAL8,MPI_MAX,mpicom.com,mpicom.ierr)
        
        if (mpicom.rank == 0) write(*,*) 'Iter, error, time', nt, err, real(MPI_Wtime() - time_1,4)
    end subroutine
    
    subroutine read_depth()
        open(1, file = "depth.txt")
        read(1,*) depth_g(1:nx,1:ny)
        close(1)
    end subroutine
    
    subroutine init_parallel()        
        integer(kind=4) :: mpierr, comm, ncpus
        logical(kind=1) :: test_partition, hilbert, partition_hard
        integer(kind=4) :: test_numprocs, iter_partition, nblocks_x, nblocks_y, gcx, gcy
        real   (kind=8) :: weight_3d
        NAMELIST /PARALLEL_PARAMS/ test_partition, test_numprocs, iter_partition, partition_hard, hilbert, nblocks_x, nblocks_y, gcx, gcy, weight_3d
        
        logical(kind=1) :: mask_g(nx,ny)
        real   (kind=8) :: work(nx,ny), work_2d(nx,ny), work_3d(nx,ny)
        
        open(1,file="parallel.nml")
        read(1,NML=PARALLEL_PARAMS)
        close(1)
        
        comm = MPI_COMM_WORLD
        call mpi_init(mpierr)
        call mpi_comm_size(comm,ncpus,mpierr)
        call mpi_comm_rank(comm,myid,mpierr)
        
        ! init weights
        mask_g = depth_g(1:nx,1:ny)>0
        work_2d = real(abs(mask_g),8)
        work_3d = real(depth_g(1:nx,1:ny),8)
        work_3d = work_3d / (sum(work_3d)/count(mask_g))
        if (weight_3d > 100.0) then
            work = work_3d
        else
            work = weight_3d * work_3d + work_2d
        end if

        ! try partition
        if (test_partition) then
            if (hilbert) then
                call partition.init_hilbert(1,nx,1,ny,gcx,gcy,nblocks_x,nblocks_y,mask_g, &
                work,work_2d,work_3d,test_numprocs,myid,test_partition,iter_partition,partition_hard)       
            else
                call partition.init_1block_per_core(1,nx,1,ny,gcx,gcy,mask_g, &
                work,work_2d,work_3d,test_numprocs,myid,test_partition,iter_partition,partition_hard)
            end if
            if(myid == 0) call partition.write_blocks_info()
            if(myid == 0) write(*,*) 'if no convergence, try increasing iter_partition in parallel.nml'
            if(myid == 0) write(*,*) 'if no error/warning occur, partition is fine'
            stop
        end if
        
        ! partition init
        if (hilbert) then
            call partition.init_hilbert(1,nx,1,ny,gcx,gcy,nblocks_x,nblocks_y,mask_g, &
            work,work_2d,work_3d,ncpus,myid,test_partition,iter_partition,partition_hard)       
        else
            call partition.init_1block_per_core(1,nx,1,ny,gcx,gcy,mask_g, &
            work,work_2d,work_3d,ncpus,myid,test_partition,iter_partition,partition_hard)
        end if
        if(myid == 0) call partition.write_blocks_info()
        
        ! init exchanges on partition
        call mpicom.init(comm)
        call mpicom.init_exchange(partition,gcx,gcy) 

        ! get local idx
        call partition.local_idx(is,ie,js,je,mpicom.rank)
        
        allocate(mask(is:ie,js:je))
        allocate(depth(is:ie,js:je))
        
        ! get local mask and depth fields
        call partition.get_local_mask(mask,0,0,mpicom.rank)
        depth = depth_g(is:ie,js:je)*abs(mask)
        
        ke = maxval(depth) ! maximum integer depth on local CPU
        
        allocate(exch(is-1:ie+1,js-1:je+1))
        exch = depth_g(is-1:ie+1,js-1:je+1)
    end subroutine
    
end program
