#ifndef NXSIZE
#define NXSIZE 101
#endif

#ifndef NYSIZE
#define NYSIZE 641
#endif

#ifndef NUM_PROCS
#define NUM_PROCS 64
#endif

#ifndef NZSLICE 
#define NZSLICE 10
#endif

program main
  !< Main program
  !-------------------------------------------------
  use solver     ,  only: iterate_one_more_time_step
  use solver     ,  only: control
  use convergence,  only: converged
  use solver, only:  start_run
  use solver, only: finish_run
  
#ifdef _OPENACC
  integer :: devnum
  integer :: setDevice
#endif
  
#ifdef _OPENACC
  devnum = setDevice(nb_procs,rank)
#endif
  

!--------Start---------!
  call start_run()

  do while ((control%current_iter <= control%max_iters) .and. (.not. converged(control)) .and. (.not. control%Halt))
     call iterate_one_more_time_step()
  end do

  call finish_run()
!--------Stop---------!

end program main


  
  !!$acc if defined(_OPENACC)
#ifdef _OPENACC

  function setDevice(nprocs,myrank) 

    use iso_c_binding
    use openacc
    implicit none
    include 'mpif.h'

    interface
      function gethostid() BIND(C)
        use iso_c_binding
        integer (C_INT) :: gethostid
      end function gethostid
    end interface

    integer, intent(in) :: nprocs, myrank
    integer, allocatable :: hostids(:), localprocs(:)
    integer :: hostid, ierr, numdev, mydev, i, numlocal
    integer :: setDevice

    allocate(hostids(nprocs))
    allocate(localprocs(nprocs))
  
  

  ! get the hostids so we can determine what other processes are on this node
    hostid = gethostid()
    CALL mpi_allgather(hostid,1,MPI_INTEGER,hostids,1,MPI_INTEGER, &
                     MPI_COMM_WORLD,ierr)

  ! determine which processors are on this node
    numlocal=0
    localprocs=0
    do i=1,nprocs
      if (hostid .eq. hostids(i)) then
        localprocs(i)=numlocal
        numlocal = numlocal+1
      endif
    enddo

  ! get the number of devices on this node
    numdev = acc_get_num_devices(ACC_DEVICE_NVIDIA)

    if (numdev .lt. 1) then
      print *, 'ERROR: There are no devices available on this host.  &
                ABORTING.', myrank
      stop
    endif

  ! print a warning if the number of devices is less then the number
  ! of processes on this node.  Having multiple processes share devices is not
  ! recommended.
    if (numdev .lt. numlocal) then
     if (localprocs(myrank+1).eq.1) then
       ! print the message only once per node
     print *, 'WARNING: The number of process is greater then the number  &
               of GPUs.', myrank
     endif
     mydev = mod(localprocs(myrank+1),numdev)
    else
     mydev = localprocs(myrank+1)
    endif

    call acc_set_device_num(mydev,ACC_DEVICE_NVIDIA)
    call acc_init(ACC_DEVICE_NVIDIA)
    setDevice = mydev
 
    deallocate(hostids)
    deallocate(localprocs)

    end function setDevice
    !!$acc end if
#endif



