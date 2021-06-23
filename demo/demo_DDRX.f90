! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2021

program demo

    use specfab    
    use netcdf
    
    implicit none

    ! Numerics
    real, parameter    :: dt = 0.1 ! 1/10 year
    integer, parameter :: Nt = 1000 ! Number of time steps

    ! Constants and argv strings    
    integer, parameter :: dp = 8
    integer :: ii,tt ! loop vars
    character(len=5) :: arg_exp 

    ! Fabric state and evolution
    integer, parameter :: Lcap = 4 ! Expansion series truncation
    complex(kind=dp), allocatable :: nlm(:), nlmiso(:), dndt(:,:) ! Series expansion coefs and evolution matrix
    real(kind=dp) :: a2(3,3), a4(3,3,3,3), da2dt(3,3), da4dt(3,3,3,3)
    real(kind=dp) :: ugrad(3,3), tau(3,3) ! Large-scale deformation tensors

    ! DRX
    real(kind=dp), parameter :: Gamma0 = 1.0 ! Sets DRX time scale

    ! For dumping state to netCDF
    complex(kind=dp), allocatable   :: nlm_save(:,:)
    real(kind=dp)                   :: a2_true_save(3,3,Nt), a2_save(3,3,Nt), a4_save(3,3,3,3,Nt), a4_true_save(3,3,3,3,Nt)
    character(len=30) :: fname_sol
    integer :: ncid, c_did, time_did, dim_did, pair_did ! Dimension IDs
    integer :: id_cre,id_cim,id_lm, id_a2, id_a2_true, id_a4, id_a4_true ! Var IDs

    if (command_argument_count() .ne. 1) then
        print *,'usage: ./demo ugrad'
        call exit(0)
    end if

    !-------------------------------------------------------------------
    ! Velocity gradient tensor
    !-------------------------------------------------------------------

    call get_command_argument(1, arg_exp)
    select case (arg_exp)

        ! RECALL COLUMN FIRST IN FORTRAN
        
        ! Uniaxial compression (uc) and uniaxial extension (ue)
        case ('uc_xx', 'ue_xx')
            ugrad = reshape([-1.0,0.,0., 0.,0.5,0., 0.,0.,0.5], [3, 3])
        case ('uc_yy', 'ue_yy')
            ugrad = reshape([0.5,0.,0., 0.,-1.0,0., 0.,0.,0.5], [3, 3])
        case ('uc_zz', 'ue_zz')
            ugrad = reshape([0.5,0.,0., 0.,0.5,0., 0.,0.,-1.], [3, 3])

        ! Confined compression (cc)
        case ('cc_zx')
            ugrad = reshape([1.,0.,0., 0.,0.,0., 0.,0.,-1.], [3, 3]) ! confined in y
        case ('cc_zy')
            ugrad = reshape([0.,0.,0., 0.,1.,0., 0.,0.,-1.], [3, 3]) ! confined in x
        case ('cc_yx')
            ugrad = reshape([1.,0.,0., 0.,-1.,0., 0.,0.,0.], [3, 3]) ! confined in z
                    
        ! Simple shear (ss)
        case ('ss_xz')
            ugrad = reshape([0.,0.,0., 0.,0.,0., 1.,0.,0.], [3, 3]) 
        case ('ss_xy')
            ugrad = reshape([0.,0.,0., 1.,0.,0., 0.,0.,0.], [3, 3]) 
        case ('ss_yz')
            ugrad = reshape([0.,0.,0., 0.,0.,0., 0.,1.,0.], [3, 3]) 

        case default
            print *,'argv error: valid "ugrad" are "uc_ij", "ue_ij", "cc_ij", "ss_ij" where i,j are any of x,y,z'
            call exit(0)
    end select

    select case (arg_exp)
        case ('ue_xx','ue_yy','ue_zz')
            ugrad = -1*ugrad
    end select

    ! Of course tau != eps, but we use eps to define the stress tensors *as if* stress and strain-rate were coaxial (Glen's law)
    tau = (ugrad+transpose(ugrad))/2 
            
    !-------------------------------------------------------------------
    ! Initialize
    !-------------------------------------------------------------------

    call initspecfab(Lcap) ! nlm_len is now defined (number of expansion coeffcients, i.e. #DOFs)

    allocate(nlm_save(nlm_len,Nt))
    allocate(dndt(nlm_len,nlm_len))
    
    allocate(nlm(nlm_len))
    allocate(nlmiso(nlm_len))
    nlm    = [(0,ii=1,nlm_len)] ! Expansion coefs "n_l^m" are saved in the 1D array "nlm". Corresponding (l,m) values for the i'th coef (i.e. nlm(i)) are (l,m) = (lm(1,i),lm(2,i))
    nlmiso = [(0,ii=1,nlm_len)] 
    
    ! Normalize such that integral over ODF is 1
    nlmiso(1) = (1,0)
    nlmiso(1) = nlmiso(1)/f_ev_c0(nlmiso(1)) 

    ! Init with isotropy
    nlm = nlmiso      ! Spectral expansion
    a2 = a2_ij(nlm)   ! Tensorial expansion
    a4 = a4_ijkl(nlm) ! ...
    
    ! Save initial state
    nlm_save(:,1)       = nlm
    a2_true_save(:,:,1) = a2 
    a2_save(:,:,1)      = a2 
    a4_true_save(:,:,:,:,1) = a4
    a4_save(:,:,:,:,1) = a4
    
    !-------------------------------------------------------------------
    ! Integrate
    !-------------------------------------------------------------------
    
    ! Model fabric evolution by representing fabric both spectrally and tensorially 

    write(*,"(A13,I4,A5,F12.10,A4,I2,A10,I3,A1)") 'Numerics: Nt=', Nt, ', dt=', dt, ', L=', Lcap, ' (nlm_len=',nlm_len,')'

    do tt = 2, Nt

        !print *, '*** a^(2), spec.: ', a2_true_save(:,:,tt-1)
        !print *, '*** a^(2), tens.: ', a2_save(:,:,tt-1)

        !write(*,"(A9,I3)") '*** Step ', tt

        !print *,'*** nlm consistency check ***'
        !print *, nlm
        !print *, a2_to_nlm(a2_ij(nlm))
        !print *, a4_to_nlm(a2_ij(nlm), a4_ijkl(nlm))
        
        ! Spectral expansion
        dndt = Gamma0 * dndt_ij_DDRX(nlm, tau) ! Tau is assumed constant, but since DRX rate depends on the state itself (i.e. DRX is nonlinear), it must be called for each time step.
        nlm_save(:,tt) = nlm + dt * matmul(dndt, nlm) ! Spectral coefs evolve by a linear transformation
        a2_true_save(:,:,tt) = a2_ij(nlm)
        a4_true_save(:,:,:,:,tt) = a4_ijkl(nlm)

        ! Tensorial expansion
        call daidt_DDRX(tau, Gamma0, a2, a4, da2dt, da4dt) ! Sets rate-of-change matrices da2dt and da4dt
        a2_save(:,:,tt)     = a2 + dt*da2dt 
        a4_save(:,:,:,:,tt) = a4 + dt*da4dt
        
        ! Set previous state = present state for next loop entry
        nlm = nlm_save(:,tt) 
        a2  = a2_save(:,:,tt)
        a4  = a4_save(:,:,:,:,tt)  
        
    end do
    
    !-------------------------------------------------------------------
    ! Dump solution to netCDF
    !-------------------------------------------------------------------
    
    write (fname_sol,"('solutions/DDRX_',A5,'.nc')") arg_exp
    call check( nf90_create(fname_sol, NF90_CLOBBER, ncid) )
    
    call check(nf90_put_att(ncid,NF90_GLOBAL, "tsteps", Nt))
    call check(nf90_put_att(ncid,NF90_GLOBAL, "dt",     dt))
    call check(nf90_put_att(ncid,NF90_GLOBAL, "L",      Lcap))
    call check(nf90_put_att(ncid,NF90_GLOBAL, "ugrad",  reshape(ugrad, [size(ugrad)]) ))
    
    call check( nf90_def_dim(ncid, "DOF",    nlm_len,   c_did) )
    call check( nf90_def_dim(ncid, "tstep",  Nt,        time_did) )
    call check( nf90_def_dim(ncid, "dim",    3,         dim_did) )
    call check( nf90_def_dim(ncid, "pair",   2,         pair_did) )

    call check( nf90_def_var(ncid, "lm",    NF90_INT,    [pair_did, c_did], id_lm) )
    
    call check( nf90_def_var(ncid, "c_re",  NF90_DOUBLE, [c_did,   time_did], id_cre) )
    call check( nf90_def_var(ncid, "c_im",  NF90_DOUBLE, [c_did,   time_did], id_cim) )
    
    call check( nf90_def_var(ncid, "a2",      NF90_DOUBLE, [dim_did,dim_did, time_did], id_a2) )
    call check( nf90_def_var(ncid, "a2_true", NF90_DOUBLE, [dim_did,dim_did, time_did], id_a2_true) )    
    call check( nf90_def_var(ncid, "a4",      NF90_DOUBLE, [dim_did,dim_did,dim_did,dim_did, time_did], id_a4) )
    call check( nf90_def_var(ncid, "a4_true", NF90_DOUBLE, [dim_did,dim_did,dim_did,dim_did, time_did], id_a4_true) )    
    
    call check( nf90_enddef(ncid) )

    call check( nf90_put_var(ncid, id_lm,    lm(:,1:nlm_len)) )
        
    call check( nf90_put_var(ncid, id_cre,   real(nlm_save)) )
    call check( nf90_put_var(ncid, id_cim,   aimag(nlm_save)) )

    call check( nf90_put_var(ncid, id_a2,      a2_save) )
    call check( nf90_put_var(ncid, id_a2_true, a2_true_save) )
    call check( nf90_put_var(ncid, id_a4,      a4_save) )
    call check( nf90_put_var(ncid, id_a4_true, a4_true_save) )
    
    call check( nf90_close(ncid) )

    print *, 'Solution dumped in ', fname_sol
    print *, "Plot result:"
    write(*,"(A26,A5)") "python3 plot_demo_DDRX.py ", arg_exp

contains

  
    subroutine check(status)
        implicit none
        integer, intent (in) :: status
        if(status /= nf90_noerr) then 
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check  
    
end program

