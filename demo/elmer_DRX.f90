! N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <dlilien90@gmail.com>, 2021
! IBOF taken from Elmer/ICE, credit to Fabien Gillet-Chaulet, Olivier Gagliardini, and others
! See Gillet-Chaulet et al., 2005 and 2006 for this particular implementation

program demo

    use specfab    
    use tensorproducts
    use netcdf
    
    implicit none

    ! Numerics
    real, parameter    :: dt = 0.01 ! 1/10 year
    integer, parameter :: Nt = 5 ! Number of time steps

    ! Constants and argv strings    
    integer, parameter :: dp = 8
    integer :: ii,tt ! loop vars
    character(len=5) :: arg_exp 

    ! Fabric state and evolution
    integer, parameter :: Lcap = 4 ! Expansion series truncation
    complex(kind=dp), allocatable :: nlm(:), nlmiso(:), dndt(:,:) ! Series expansion coefs and evolution matrix
    real(kind=dp) :: a2(3,3), a4(3,3,3,3)
    real(kind=dp) :: a2check(3,3), a4check(3,3,3,3), da2(3,3), da2check(3,3)
    real(kind=dp) :: ai2(6), ae2(5), ae4(9), f
    real(kind=dp) :: ugrad(3,3), tau(3,3) ! Large-scale deformation tensors

    ! DRX
    real(kind=dp), parameter :: Gamma0 = 0.1 ! Sets DRX time scale

    ! For dumping state to netCDF
    complex(kind=dp), allocatable   :: nlm_save(:,:)
    real(kind=dp)                   :: a2_save(3,3,Nt), a2_true_save(3,3,Nt)
    character(len=30) :: fname_sol
    integer :: ncid, c_did, time_did, dim_did, pair_did ! Dimension IDs
    integer :: id_cre,id_cim,id_lm, id_a2, id_a2_true ! Var IDs
    integer :: i,j,k,l

    if (command_argument_count() .ne. 1) then
        print *,'usage: ./elmer_DRX ugrad'
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
    
    nlmiso(1) = (1,0)
    nlmiso(1) = nlmiso(1)/f_ev_c0(nlmiso(1)) ! Normalize such that integral over ODF is 1

    ! Init with isotropy
    nlm = nlmiso
    a2 = a2_ij(nlm) ! Init corresponding tensorial formulation
    
    nlm_save(:,1)  = nlm
    a2_true_save(:,:,1) = a2 
    a2_save(:,:,1)      = a2 
    
    !-------------------------------------------------------------------
    ! Integrate
    !-------------------------------------------------------------------
    
    ! Model fabric evolution by representing fabric both spectrally and tensorially 
    ! ...in order to comapre the effect of modelling d/dt a^(2) only using the python plotting script

    write(*,"(A13,I4,A5,F12.10,A4,I2,A10,I3,A1)") 'Numerics: Nt=', Nt, ', dt=', dt, ', L=', Lcap, ' (nlm_len=',nlm_len,')'

    do tt = 2, Nt

        ! print *, '*** a^(2),  spectral: ', a2_true_save(:,:,tt-1)
        ! print *, '*** a^(2), tensorial: ', a2_save(:,:,tt-1)

        write(*,*) ''
        write(*,"(A9,I3)") '*** Step ', tt

        ! print *,'*** nlm consistency check ***'
        ! print *, nlm
        ! print *, a4_to_nlm(a2_ij(nlm), a4_ijkl(nlm))
        
        ! Spectral representation
        dndt = Gamma0 * dndt_ij_DRX(nlm, tau) ! Tau is assumed constant, but since DRX rate depends on the state itself (i.e. DRX is nonlinear), it must be called for each time step.
        nlm_save(:,tt) = nlm + dt * matmul(dndt, nlm) ! Spectral coefs evolve by a linear transformation
        a2_true_save(:,:,tt) = a2_ij(nlm)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! These variables are not used, this is only to check the IBOF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        a2check = a2_ij(nlm_save(:, tt-1))
        a4check = a4_ijkl(nlm_save(:, tt-1))
        ai2 = (/ a2check(1, 1), a2check(2,2), a2check(3,3),&
                 a2check(1,2), a2check(2,3), a2check(1,3) /)
        ae4 = 0.0
        CALL IBOF(ai2, ae4)
        ae2 = a2_to_ae2(a2check)
        a4 = ae4_to_a4(ae2, ae4)
        da2 = da2dt_DRX(tau, a2check, a4)
        da2check = da2dt_DRX(tau, a2check, a4check)
        WRITE(*,*) 'Difference in a4 caused by IBOF'

        ! First check that the resulting a4 preserves the trace rules
        do i=1,3
          do j=1,3
            if (abs((a2check(i,j) - (a4(i,j,1,1) + a4(i,j,2,2) + a4(i,j,3,3))) / a2check(i,j)).GE.1.0e-6) THEN
                write(*,*) 'Error in',i,j,' component of a2 compared to a4'
                write(*,*) a2check(i,j), (a4(i,j,1,1) + a4(i,j,2,2) + a4(i,j,3,3))
            END IF
          END DO
        END DO

        ! And then see what this does to a4
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                IF ((abs(a4(i,j,k,l)).GT.1e-8).OR.(abs(a4check(i,j,k,l)).GT.1e-8)) THEN
                  IF (abs((a4(i,j,k,l)-a4check(i,j,k,l)) / a4(i,j,k,l)).GT.1.0e-4) THEN
                    ! dont print too much because symmetry
                    IF ((i.LE.j).AND.(j.LE.k).AND.(k.LE.l)) THEN
                      write(*,'(I1,x,I1,x,I1,x,I1,x,A,ES12.4,x,ES12.4,x,F10.5)') i,j,k,l,&
                          'ibof, true, % diff: ',a4(i,j,k,l),&
                        a4check(i,j,k,l),abs((a4(i,j,k,l)-a4check(i,j,k,l))/a4check(i,j,k,l) * 100)
                    END IF
                  END IF
                END IF
              end do
            end do
          end do
        end do
        WRITE(*,*) 'Difference in da2dt caused by IBOF'
        do i=1,3
            do j=1,3
                IF ((abs(da2(i,j)).GT.1e-8).OR.(abs(da2check(i,j)).GT.1e-8)) THEN
                    IF (abs((da2(i,j) - da2check(i,j)) / da2(i,j)).GT.1.0e-4) THEN
                      write(*,'(I1,x,I1,x,A,ES12.4,x,ES12.4,x,F10.5)') i,j,&
                          'ibof, true, % diff: ', da2(i,j), da2check(i,j),&
                          abs((da2(i,j) -da2check(i,j)) / da2check(i,j) * 100)
                    END IF
                END IF
            end do
        end do


        ai2 = (/ a2(1, 1), a2(2,2), a2(3,3), a2(1,2), a2(2,3), a2(1,3) /)
        ae4 = 0.0
        CALL IBOF(ai2, ae4)
        write(*,*) '<D>'
        write(*,*) (doubleinner22(MATMUL(tau,tau),a2)-doubleinner22(tau,doubleinner42(a4, tau))) / doubleinner22(tau,tau)
        write(*,*) doubleinner22(MATMUL(tau,tau),a2), doubleinner22(tau,doubleinner42(a4, tau))    
        ae2 = a2_to_ae2(a2)
        a4 = ae4_to_a4(ae2, ae4)
        da2 = da2dt_DRX(tau, a2, a4)
        WRITE(*,*) 'Difference in da2dt overall (IBOF and difference in tens/spec a2)'
        do i=1,3
            do j=1,3
                IF ((abs(da2(i,j)).GT.1e-8).OR.(abs(da2check(i,j)).GT.1e-8)) THEN
                    IF (abs((da2(i,j) - da2check(i,j)) / da2(i,j)).GT.1.0e-4) THEN
                      write(*,'(I1,x,I1,x,A,ES12.4,x,ES12.4,x,F10.5)') i,j,&
                          'ibof, true, % diff: ', da2(i,j), da2check(i,j),&
                          abs((da2(i,j) -da2check(i,j)) / da2check(i,j) * 100)
                    END IF
                END IF
            end do
        end do
        a2_save(:,:,tt) = a2 + dt * Gamma0 * da2dt_DRX(tau, a2, a4)


        ! Set previous state = present state for next loop entry
        nlm = nlm_save(:,tt) 
        a2  = a2_save(:,:,tt) 
        
    end do
    
    !-------------------------------------------------------------------
    ! Dump solution to netCDF
    !-------------------------------------------------------------------
    
    write (fname_sol,"('solutions/ELMER_DRX_',A5,'.nc')") arg_exp
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
    
    call check( nf90_enddef(ncid) )

    call check( nf90_put_var(ncid, id_lm,    lm(:,1:nlm_len)) )
        
    call check( nf90_put_var(ncid, id_cre,   real(nlm_save)) )
    call check( nf90_put_var(ncid, id_cim,   aimag(nlm_save)) )

    call check( nf90_put_var(ncid, id_a2,      a2_save) )
    call check( nf90_put_var(ncid, id_a2_true, a2_true_save) )
    
    call check( nf90_close(ncid) )

    print *, 'Solution dumped in ', fname_sol
    print *, "Plot result:"
    write(*,"(A25,A5)") "python3 plot_elmer_DRX.py  ", arg_exp

contains

  
    subroutine check(status)
        implicit none
        integer, intent (in) :: status
        if(status /= nf90_noerr) then 
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check  

    ! This is directly copied from elmerfem/elmerice/Solvers/GolfLaw.F90
    subroutine IBOF(a2,a4)
       implicit none
       Real(dp),dimension(6),intent(in):: a2  
       Real(dp),dimension(9),intent(out):: a4  
       Real(dp):: a_11,a_22,a_33,a_12,a_13,a_23
       Real(dp):: b_11,b_22,b_12,b_13,b_23
       Real(dp):: aPlusa

       Real(dp),dimension(21) :: vec
       Real(dp),dimension(3,21) :: Mat
       Real(dp),dimension(6) :: beta
       Real(dp) :: Inv2,Inv3
       integer :: i,j      

       
       a_11=a2(1)
       a_22=a2(2)
       a_33=a2(3)
       a_12=a2(4)
       a_23=a2(5)
       a_13=a2(6)


!       write(*,*) 'a2'
!       write(*,*) a2

       !Coefficients 

      Mat(1,1)=0.217774509809788e+02_dp
      Mat(1,2)=-.297570854171128e+03_dp
      Mat(1,3)=0.188686077307885e+04_dp
      Mat(1,4)=-.272941724578513e+03_dp
      Mat(1,5)=0.417148493642195e+03_dp
      Mat(1,6)=0.152038182241196e+04_dp
      Mat(1,7)=-.137643852992708e+04_dp
      Mat(1,8)=-.628895857556395e+03_dp
      Mat(1,9)=-.526081007711996e+04_dp
      Mat(1,10)=-.266096234984017e+03_dp
      Mat(1,11)=-.196278098216953e+04_dp
      Mat(1,12)=-.505266963449819e+03_dp
      Mat(1,13)=-.110483041928547e+03_dp
      Mat(1,14)=0.430488193758786e+04_dp
      Mat(1,15)=-.139197970442470e+02_dp
      Mat(1,16)=-.144351781922013e+04_dp
      Mat(1,17)=-.265701301773249e+03_dp
      Mat(1,18)=-.428821699139210e+02_dp
      Mat(1,19)=-.443236656693991e+01_dp
      Mat(1,20)=0.309742340203200e+04_dp
      Mat(1,21)=0.386473912295113e+00_dp
      Mat(2,1)=-.514850598717222e+00_dp
      Mat(2,2)=0.213316362570669e+02_dp
      Mat(2,3)=-.302865564916568e+03_dp
      Mat(2,4)=-.198569416607029e+02_dp
      Mat(2,5)=-.460306750911640e+02_dp
      Mat(2,6)=0.270825710321281e+01_dp
      Mat(2,7)=0.184510695601404e+03_dp
      Mat(2,8)=0.156537424620061e+03_dp
      Mat(2,9)=0.190613131168980e+04_dp
      Mat(2,10)=0.277006550460850e+03_dp
      Mat(2,11)=-.568117055198608e+02_dp
      Mat(2,12)=0.428921546783467e+03_dp
      Mat(2,13)=0.142494945404341e+03_dp
      Mat(2,14)=-.541945228489881e+04_dp
      Mat(2,15)=0.233351898912768e+02_dp
      Mat(2,16)=0.104183218654671e+04_dp
      Mat(2,17)=0.331489412844667e+03_dp
      Mat(2,18)=0.660002154209991e+02_dp
      Mat(2,19)=0.997500770521877e+01_dp
      Mat(2,20)=0.560508628472486e+04_dp
      Mat(2,21)=0.209909225990756e+01_dp
      Mat(3,1)=0.203814051719994e+02_dp
      Mat(3,2)=-.283958093739548e+03_dp
      Mat(3,3)=0.173908241235198e+04_dp
      Mat(3,4)=-.195566197110461e+03_dp
      Mat(3,5)=-.138012943339611e+03_dp
      Mat(3,6)=0.523629892715050e+03_dp
      Mat(3,7)=0.859266451736379e+03_dp
      Mat(3,8)=-.805606471979730e+02_dp
      Mat(3,9)=-.468711180560599e+04_dp
      Mat(3,10)=0.889580760829066e+01_dp
      Mat(3,11)=-.782994158054881e+02_dp
      Mat(3,12)=-.437214580089117e+02_dp
      Mat(3,13)=0.112996386047623e+01_dp
      Mat(3,14)=0.401746416262936e+04_dp
      Mat(3,15)=0.104927789918320e+01_dp
      Mat(3,16)=-.139340154288711e+03_dp
      Mat(3,17)=-.170995948015951e+02_dp
      Mat(3,18)=0.545784716783902e+00_dp
      Mat(3,19)=0.971126767581517e+00_dp
      Mat(3,20)=0.141909512967882e+04_dp
      Mat(3,21)=0.994142892628410e+00_dp

       
       ! calcul des invariants
       Inv2=0.5_dp*(1._dp-(a_11*a_11+a_22*a_22+a_33*a_33+ &
            2._dp*(a_12*a_12+a_13*a_13+a_23*a_23)))
            
       Inv3=a_11*(a_22*a_33-a_23*a_23)+a_12*(a_23*a_13-a_12*a_33)+ &
             a_13*(a_12*a_23-a_22*a_13)
       
     ! polynome complet de degre 5 des 2 invariants.
         vec(1)=1._dp
         vec(2)=Inv2
         vec(3)=vec(2)*vec(2)
         vec(4)=Inv3
         vec(5)=vec(4)*vec(4)
         vec(6)=vec(2)*vec(4)
         vec(7)=vec(3)*vec(4)
         vec(8)=vec(2)*vec(5)
         vec(9)=vec(2)*vec(3)
         vec(10)=vec(5)*vec(4)
         vec(11)=vec(9)*vec(4)
         vec(12)=vec(3)*vec(5)
         vec(13)=vec(2)*vec(10)
         vec(14)=vec(3)*vec(3)
         vec(15)=vec(5)*vec(5)
         vec(16)=vec(14)*vec(4)
         vec(17)=vec(12)*vec(2)
         vec(18)=vec(12)*vec(4)
         vec(19)=vec(2)*vec(15)
         vec(20)=vec(14)*vec(2)
         vec(21)=vec(15)*vec(4)

       ! calcul des beta_bar (cf annexe C Chung)
       ! attention beta(1)=beta_bar_3 (Chung); beta(2)=beta_bar_4; beta(3)=beta_bar_6
       !           beta(4)=beta_bar_1        ; beta(5)=beta_bar_2; beta(6)=beta_bar_5

       ! calcul des trois beta en fonction du polynome
         beta(:)=0._dp
         Do i=1,3
          Do j=1,21
            beta(i)=beta(i)+Mat(i,j)*vec(j)
          End do
         End do
          
       ! calcul des 3 autres pour avoir la normalisation
         beta(4)=3._dp*(-1._dp/7._dp+beta(1)*(1._dp/7._dp+4._dp*Inv2/7._dp+8._dp*Inv3/3._dp)/5._dp- &
                  beta(2)*(0.2_dp-8._dp*Inv2/15._dp-14._dp*Inv3/15._dp)- &
                  beta(3)*(1._dp/35._dp-24._dp*Inv3/105._dp-4._dp*Inv2/35._dp+ &
                  16._dp*Inv2*Inv3/15._dp+8._dp*Inv2*Inv2/35._dp))/5._dp

         beta(5)=6._dp*(1._dp-0.2_dp*beta(1)*(1._dp+4._dp*Inv2)+ &
                  7._dp*beta(2)*(1._dp/6._dp-Inv2)/5._dp- &
                  beta(3)*(-0.2_dp+2._dp*Inv3/3._dp+4._dp*Inv2/5._dp- &
                  8._dp*Inv2*Inv2/5._dp))/7._dp

         beta(6)=-4._dp*beta(1)/5._dp-7._dp*beta(2)/5._dp- &
                   6._dp*beta(3)*(1._dp-4._dp*Inv2/3._dp)/5._dp

        ! pour avoir les beta_bar
        Do i=1,6
         beta(i)=beta(i)/3._dp
        End do
         beta(2)=beta(2)/2._dp
         beta(5)=beta(5)/2._dp
         beta(6)=beta(6)/2._dp

        !! calcul des 5 b=a.a
        b_11=a_11*a_11+a_12*a_12+a_13*a_13
        b_22=a_22*a_22+a_12*a_12+a_23*a_23
        b_12=a_11*a_12+a_12*a_22+a_13*a_23
        b_13=a_11*a_13+a_12*a_23+a_13*a_33
        b_23=a_12*a_13+a_22*a_23+a_23*a_33

        !Calcul des 9 termes de a4

        a4(1)=3._dp*beta(4)+6._dp*beta(5)*a_11+3._dp*beta(1)*a_11*a_11+&
         6._dp*beta(2)*b_11+6._dp*beta(6)*a_11*b_11+3._dp*beta(3)*b_11*b_11
        a4(2)=3._dp*beta(4)+6._dp*beta(5)*a_22+3._dp*beta(1)*a_22*a_22+&
         6._dp*beta(2)*b_22+6._dp*beta(6)*a_22*b_22+3._dp*beta(3)*b_22*b_22

        a4(3)=beta(4)+beta(5)*(a_22+a_11)+beta(1)*(a_11*a_22+2._dp*a_12*a_12)+&
         beta(2)*(b_22+b_11)+beta(6)*(a_11*b_22+a_22*b_11+4._dp*a_12*b_12)+&
         beta(3)*(b_11*b_22+2._dp*b_12*b_12)


         a4(4)=beta(5)*a_23+beta(1)*(a_11*a_23+2._dp*a_12*a_13)+beta(2)*b_23+&
          beta(6)*(a_11*b_23+a_23*b_11+2._dp*(a_12*b_13+a_13*b_12))+beta(3)*&
          (b_11*b_23+2._dp*b_12*b_13)
         a4(5)=beta(5)*a_13+beta(1)*(a_22*a_13+2._dp*a_12*a_23)+beta(2)*b_13+&
          beta(6)*(a_22*b_13+a_13*b_22+2._dp*(a_12*b_23+a_23*b_12))+beta(3)*&
          (b_22*b_13+2._dp*b_12*b_23)


         a4(6)=3._dp*beta(5)*a_13+3._dp*beta(1)*a_11*a_13+3._dp*beta(2)*b_13+&
          3._dp*beta(6)*(a_11*b_13+a_13*b_11)+3._dp*beta(3)*b_11*b_13
         a4(7)=3._dp*beta(5)*a_12+3._dp*beta(1)*a_11*a_12+3._dp*beta(2)*b_12+&
          3._dp*beta(6)*(a_11*b_12+a_12*b_11)+3._dp*beta(3)*b_11*b_12
         a4(8)=3._dp*beta(5)*a_23+3._dp*beta(1)*a_22*a_23+3._dp*beta(2)*b_23+&
          3._dp*beta(6)*(a_22*b_23+a_23*b_22)+3._dp*beta(3)*b_22*b_23
         a4(9)=3._dp*beta(5)*a_12+3._dp*beta(1)*a_22*a_12+3._dp*beta(2)*b_12+&
          3._dp*beta(6)*(a_22*b_12+a_12*b_22)+3._dp*beta(3)*b_22*b_12


         End 

    
end program
