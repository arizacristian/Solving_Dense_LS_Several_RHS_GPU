
program main
!!

!!                     compiling whit makefile
!		        make PLS2019_magma_d_v5


!! 2019 v5_3, work space for LVF


USE SBHSNE2019_magma_d_v5
use Sub_CPU_Interfaces
USE mod_IO_S 

    implicit none
    integer :: Npar
    character(len=100) :: arg1, arg2,arg3,arg4
    integer  :: n, n1,m,  Nc,k,L, Nb, Lmax,nrhs,maxNRHS
    !integer  ::  i, j, desL,desC
    integer :: opc,opcA,opcB
    integer :: LDA, LDB
    DOUBLE PRECISION, allocatable :: X(:,:),A(:,:), A0(:,:)
    DOUBLE PRECISION, allocatable :: B(:,:),B0(:,:)
    DOUBLE PRECISION, allocatable :: H(:,:),D(:,:), hteo(:,:)

    DOUBLE PRECISION, allocatable :: EEFF(:,:)
    DOUBLE PRECISION, allocatable :: FFF(:,:)
    DOUBLE PRECISION, allocatable :: workLVF(:)
    integer, allocatable          :: Iseed(:)
    integer :: memstat   !! allocacao
    !! matrix gem
    DOUBLE PRECISION, allocatable :: diag(:)
   INTEGER,  ALLOCATABLE          :: a_seed(:)
   INTEGER :: i_seed
    !DOUBLE PRECISION :: EEr, r
    DOUBLE PRECISION, allocatable :: work(:),work2n(:)

    integer :: lwork,info
    !integer :: ista2, iend2,jsta2, jend2, ibs,ibe,iAs,iAe,jAs,jAe,deEf
    !REal ( kind = 8 ) :: sMR, sMRp
    REal ( kind = 8 ) ::  start , T1, T2, finish ,  tfac,tslv,t4,t5 !,t3, ttime,T2omp,finishomp
    REal ( kind = 8 ) :: tmm,tmmsyrk	!,ttimeRef,tfacRef,tslvRef
    !integer*8 :: tF, tH 
    integer*8 :: Ofac, Osol,OMM
    DOUBLE PRECISION :: dOfac, dOsol,dOMM,dOMMsyrk
    REal ( kind = 8 ) :: GflpSfac, GflpSslv, GflpST,GflpMM,GflpMMsyrk
    REal ( kind = 8 ) :: MenLoc, MenLoc0, MenLocSET !,MenTot,maxMenLoc  !! calc max memory 
    REal ( kind = 8 ) :: B2GB  !! calc max memory
    DOUBLE PRECISION ::  dmen_trf, dmen_trs,dMen_trf_work,dMen_trs_work 
    integer ::  sumL1,sumk !,rowHH, colHH, rowFF, colFF, sumL
    integer :: jnrhs,rank	!je1,js1, js2,je2, cont,
    logical :: ejecutar
    

    character(len=1) :: normtyp
    character(len=16)             :: okay 

    !!   openMP
   !integer :: myid, nthreads,maxTh,nth
   !integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM,omp_get_max_threads
   REal ( kind = 8 ) :: omp_get_wtime
!     .. Local Scalars ..
   real ( kind = 8 ) ::  ANORM, AXBNORM, EPS, RESID, XNORM 
   double precision              :: dBnorm, derror, tol !,dAnorm, dRnorm, dXnorm
!     .. External Functions ..
   real  ( kind = 8 ) ::    DLANGE, dlamch,DLANSY
!! report
   logical :: exist
   logical :: writefile
   character(len=100) :: filereport
!! ------ control
      INTEGER            NOUT, NIN, bytes,numarg
      PARAMETER          ( NIN = 1,     bytes=8)
      integer :: nmin,nmax,nstep,nn,nk,itermin,itermax
      integer :: kdiv,kmin,kmax,kstep,nnRHS,nkrhs
      integer   :: ceil
      DOUBLE PRECISION, ALLOCATABLE :: TIME_ITER_LVF(:),TIME_ITER_LVS(:),TIME_ITER_mm(:),TIME_ITER_MMsyrk(:)
      DOUBLE PRECISION, ALLOCATABLE :: TIME_ITER_LVS_RHS(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: RESID_RHS(:,:),derror_RHS(:,:),SM_abs_AX_B(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: MenLoc_rhs(:)
      integer, ALLOCATABLE :: max_rhs(:)
      DOUBLE PRECISION ::   THRESH, SM_abs
      integer  :: ITERATIONS, LNCH
      LOGICAL          ::  CHECK
      INTEGER          :: PASSED, TESTS_FAILED, TESTS_PASSED!, TESTS_WOCHK !,TESTS_SKIPPED
      !integer  :: rank, sumk
        
!     /* names */
      CHARACTER*128     INFILE, OUTFILE, PROGNAM
      INTEGER*4 hostnm, status
      character*80 name
!-----|---------------------------------------------------------------|
     ! get hostname
     status = hostnm( name )
    
!-----|---------------------------------------------------------------|



 B2GB=1024.0d0**3

!! +++++++++++++++++++++++++++++++++++++     
            !!    /* Get program name and in/out file names  */
      CALL GETARG(0, PROGNAM)
    numarg=2
	Npar=IARGC()
	!print *, "numero de argumentos ingresados =", Npar
	if (Npar < numarg ) then
  		print *, "Erro no uso"
  		GO TO 111  ! exit
    else 
!!    /* Get  in/out file names  */
      CALL GETARG(1, INFILE)
      CALL GETARG(2, OUTFILE)
	endif
	
!!      /*     Check whether in/out file names are valid
      IF( INFILE.EQ.'' .OR. INFILE.EQ.' ' .OR.&
        OUTFILE.EQ.'' .OR. OUTFILE.EQ.' ' )then
         write(*,*)" in OR out file names are INvalid"
        goto 111        !! EXIT
      ENDIF
!      /*     Open input file, skip header and read NOUT

         OPEN( NIN, FILE = INFILE, STATUS = 'OLD', ACTION = 'READ' )
         READ( NIN, FMT = * )
         READ( NIN, FMT = * ) NOUT

      
!      /*     Open output file
      IF( NOUT.NE.0 .AND. NOUT.NE.6  )THEN
        OPEN( NOUT, FILE = OUTFILE, STATUS = 'UNKNOWN', ACTION = 'WRITE' )
      ENDIF
      !     /*     ==== Read  data ====================================

         READ( NIN, FMT = * ) itermin,itermax
!!        Read task dimensions array
         READ( NIN, FMT = * ) nmin,nmax,nstep
         READ( NIN, FMT = * ) ceil
         READ( NIN, FMT = * ) THRESH
         READ( NIN, FMT = * ) kdiv,kmin,kmax,kstep
         !READ( NIN, FMT = * ) Lmin,Lmax,Lstep
         
        CLOSE( NIN )
      
      
        CHECK = ( THRESH.GT.0.0E+0 )
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     /*     Write table header to output

         WRITE( NOUT,   *          ) '# PLS LVF & PLS LVS performance results, pogram=',trim(PROGNAM)
         WRITE( NOUT,   *          ) '# Par in:'
         WRITE( NOUT, '(a,I8)'     ) '# NOUT, device out                   =', NOUT
         WRITE( NOUT, '(a,2I8)'    ) '# itermin,itermax                    =', itermin,itermax
         WRITE( NOUT, '(a,3I8)'    ) '# Loop N, n_mim, n_max, n_step       =', nmin,nmax,nstep
         WRITE( NOUT, '(a,4I8)'    ) '# Loop RHS:: kdiv,kmin,kmax,kstep    =', kdiv,kmin,kmax,kstep
         WRITE( NOUT, '(a,3I8)'    ) '# ceiling LDA yes=1 ceil             =', ceil
         WRITE( NOUT, '(a,F10.3,A,L1)') '# THRESH, check(THRESH>0.0)          =', THRESH," ; CHECK= ",CHECK
         !WRITE( NOUT, '(a,3I8)') 'Loop L, L_min, L_max, L_step         =', Lmin,Lmax,Lstep
         WRITE( NOUT, '(A,A,a)'    ) '# run in host name                   =', "  ",trim(name)
         WRITE( NOUT, '(A)'    ) ''
         WRITE( NOUT, '(A)'    ) ''




    !! Initialize MAGMA
    call magmaf_init()


         write(NOUT,514) "#  n","jNrhs","LDA","Nc","maxRHS","L","it","#OK","Erro","RESID","S_AB(AX-B)","tfac","tslv","tsolT", &
                            "GFPS_F","GFPS_S","GFPS_SV","GFPS_Mrk","GFPS_SMM","M_f(MB)","M_s(MB)","TM (GB)"


!     /*     ====  get L and Nc ====================================

    


    !call cpu_time ( start )
    !start=MPI_WTIME()
    start=omp_get_wtime()

    

    


ejecutar=.true.
if (ejecutar .eqv. .true.) then


        nk=0 
        NN=(nmax-nmin)/nstep+1
        ITERATIONS=itermax
        !!  NRHS
        nnRHS=(kmax-kmin)/kstep+1


  do N = nmin , nmax , nstep
    nk = nk+1

     call get_L_Nc(n,L,nc)   !! cal Nc and L
     
        if(ceil==1)then
                LDA= ceiling((real(N))/16)*16       !Align arrays on 128-byte boundaries (16 elem3nt)
                if(mod(LDA,256)==0)LDA=LDA+16
                !print *,"N=",N," LDA",LDA
        else
                LDA=N
        endif
        LDB=LDA 
        !print *, "LDA= ", LDA
        
     
   
    
    !print *, "call get_L_Nc n= ", n, "L =",L,"Nc =",Nc
    TESTS_PASSED=0
    TESTS_FAILED=0
    ITERATIONS=itermax-int(real(int((10*(real(nk)/real(nn)))))/10*real(itermax))
    if(ITERATIONS < itermin)ITERATIONS=itermin
    ALLOCATE( TIME_ITER_lvf( ITERATIONS ) )
    ALLOCATE( TIME_ITER_lvs( ITERATIONS ) )
    ALLOCATE( TIME_ITER_mm( ITERATIONS ) )
    ALLOCATE( TIME_ITER_mmsyrk( ITERATIONS ) )
    ALLOCATE( TIME_ITER_lvs_RHS( ITERATIONS,NNRHS ) )
    ALLOCATE( RESID_RHS( ITERATIONS,NNRHS ) )
    ALLOCATE( derror_RHS( ITERATIONS,NNRHS ) )
    ALLOCATE( SM_abs_AX_B( ITERATIONS,NNRHS ) )
    ALLOCATE( MenLoc_rhs(NNRHS))
    ALLOCATE( max_rhs(NNRHS))
    allocate (work2n(2*n))
    
    DO  LNCH = 1, ITERATIONS
      MenLoc=0.0    ; MenLoc0=0.0 
        !print *, "Gerando X(m,n)"

	!!! control genetaiing problem
      opcA=1   !!!  get A
      opcB=3   !!  for get HS


      if(opcA==1)then        !!! 
        !print *, " n= ",n
        m=n
	    allocate (X(M,n))
		    MenLoc0=((real(m)*real(n))*bytes)/B2GB ; MenLoc=MenLoc+MenLoc0

        call random_seed ( SIZE = i_seed )
                !print *, " i_seed ",i_seed
        allocate ( a_seed(1:i_seed) )
   	    call random_seed ( GET =a_seed )
                !print *, " a_seed ",a_seed
  		!a_seed=(/ -7298322282102941213,  4048769279036287019,  4561622708110880578,  1265253876178542821,  8732455838606679523, &
  		!           7809519287067318445,  -604845786819106709, -7644001736620384689,  7652339853076689456,  8412705608369769683, &
  		!           2611745309122114329,  2627111933639368995, -7774130911972220177,  2473139245134731268, -8803254973330007476, &
  		!          -3242008929124045442,                    0 /)
                  !!!  senai
                 !a_seed=(/    287027030,  -719361131,   574274270,   292048305,   185733336, -1598963619,   572469522, &
	!		  1446716853,   437591706,  1398099429,   570932571, -1177695979 /)
                  !print *, " a_seed ",a_seed
		   !!atribuir o novo vetor-semente (SEED) para cada porcesador
  	    CALL RANDOM_SEED(put=a_seed)                
 	    CALL random_number(X)  !! entre 0 e 1
            !!call mat_output (m,n,X,"X_m>n")  
  	    X=20.0d0*X-10.0d0
 	        !call mat_output (m,n,X,"X_m>n")
 	      !do j=1,min(m,n)
 	      !     X(j,j)=10.0*X(j,j)
 	      !enddo
 	      !call mat_output (m,n,X,"X_m>n")

	  !print *, "Gerando BLOCK hermitian cieficient matrix  (n,n)"
	  allocate (A(LDA,n))
		MenLoc0=((real(LDA)*real(n))*bytes)/B2GB ; MenLoc=MenLoc+MenLoc0
	  !A=0.0
         !print *, " Star X^T*X "
         T1 =omp_get_wtime()  
         !print *, "X (1:10,1)", X(1:10,1)    
         !CALL magmaf_dGEMM_cpu('T','N',N,N,M,1.0d0,X,M,X,M,0.0d0, A,N )  !  X^T*X
         !CALL magmaf_BIG_dgemm_cpu('T','N',N,N,M,1.0d0,X,M,X,M,0.0d0, A,N)
      call magmaf_BIG_dsyrk_cpu("L","T",N,M,1.0d0,X,M,0.0d0,A,LDA)  !! L(yes), U(?)
         !print *, "A (1:10,1)", A(1:10,1)                      
           !call mat_output (n,n,A,"A")
          !! symetric
          !call make_dsym("L",n,A)
          !call test_dsym(n,A,info)
         			!call random_number(A)
         !call mat_output (n,n,A,"A")
         T2 =omp_get_wtime()
         !print *, "time set A X^T*X= ",T2-T1
         TIME_ITER_MMsyrk(LNCH)= T2-T1
          !!  O sgemm= 2*m*k*n
         GflpMM  =(2.0d0*dble(n)*dble(n)*dble(n))/((T2-T1)*1000000000D0)
         
         !print *, "Gflops set A X^T*X= ",GflpMM

	    if(opcB==2 .or. opcB==3)then
        	deallocate (X)
        	MenLoc0=-((real(m)*real(n))*bytes)/B2GB ; MenLoc=MenLoc+MenLoc0
        endif

    elseif(opcA==2 .or. opcA==3)then  !! direct gem sym matrix
       m=n
      
            !print *, "Gerando sym matrix using dlagsy (n,n)"
       allocate (A(LDA,n))
		    MenLoc0=((real(LDA)*real(n))*bytes)/B2GB ; MenLoc=MenLoc+MenLoc0
	   A=0.0
       allocate (diag(n),iseed(4))
       call random_number(diag)
       !diag=diag*1000.0d0
       iseed(1)=200 ; iseed(2)=500; iseed(3)=2000;iseed(4)=4555
       call dlagsy( N, n-1, diag, A, N, iseed, work2n, INFO )
            !print *, "info =",info
            !call gem_blkTri(nc,L,A)
            !call mat_output (n,n,A,"A")

    end if

    IF( CHECK .AND. L==1 ) THEN
        !print *, " start copy A to A0, for verification"
       allocate (A0(LDA,n))
       MenLoc0=((real(LDA)*real(n))*bytes)/B2GB ; MenLoc=MenLoc+MenLoc0
       !A0=A       !!  for verification
       call dlacpy("T",n,n,A(1,1),LDA,A0(1,1),LDA)
    ENDIF

    !print *, "solver lapack dsytrf+dsytrs2"
    !!   solve con lapack, referecia 
    !LWORK =64*n
    !allocate (WORK( LWORK ))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!print *, ""
!print *, ""
!print *, ""


    !print *, "solver art 2010 serial version 2019  F+S,  nrhs=",nrhs, " 2 subrutines kernel DsyTRF DsyTRS"
    !allocate (IPIVnc(nc))
   
    allocate (EEFF(nc,nc*L), stat=memstat) !! work space  EEFF
    if ( memstat /= 0 ) then
     stop "Error in memory allocation for EEFF"
    endif
		MenLoc0=((real(nc)*real(nc*L))*bytes)/B2GB ; MenLoc=MenLoc+MenLoc0
    sumL1=((L-1)*((L-1)+1))/2
    !!allocate (FFF(nc,nc*sumL1), stat=memstat) !! work space  FFF
    allocate (FFF(nc,nc*(sumL1+(L-1))), stat=memstat) !! work space  FFF
    if ( memstat /= 0 ) then
     stop "Error in memory allocation for FFF"
    endif
		MenLoc0=((real(nc)*real(nc*(sumL1+(L-1))))*bytes)/B2GB ; MenLoc=MenLoc+MenLoc0

            dMen_trf=(dble(LDA)*dble(n)*dble(bytes))/(1024.0**2)
            dMen_trf=dMen_trf+((dble(nc)*dble(nc*L))*dble(bytes))/(1024.0**2)
            dMen_trf=dMen_trf+((dble(nc)*dble(nc*(sumL1+(L-1))))*dble(bytes))/(1024.0**2)



    call get_lwork_magma_dpoLVF(L,nc,lwork)
    allocate (workLVF(lwork), stat=memstat) !! work space  FFF
    if ( memstat /= 0 ) then
     stop "Error in memory allocation for workLVF"
    endif
        MenLoc0=((real(lwork))*bytes)/B2GB; ; MenLoc=MenLoc+MenLoc0
        dMen_trf_work=(dble(lwork)*dble(bytes))/(1024.0**2)
        
        
    MenLocSET=MenLoc

            !print *, " Mem Set A0,workspace  = ", MenLocSET, " GB "
            !print *, ""
	        !print *, "++++ Star   magmaf_Outcore_PoTRF_d_V5_4 "
         !T1 = MPI_WTIME()
         T1 =omp_get_wtime()

    call magmaf_Outcore_PoTRF_d_V5_4(n,Nc,L,A,LDA,EEFF,FFF,workLVF)


        
    
         !T2 = MPI_WTIME()
	    T2 =omp_get_wtime()
         !print *, "++++ finish magmaf_Outcore_PoTRF_d_V5_4 "
        !print *, ""
        !print *, " time magmaPoTRF ", T2-T1," s"
         tfac=t2-t1
         TIME_ITER_lvf( LNCH ) = tfac
         !call mat_output (n,n,A,"A  output")
         do k=1,L
             rank=nc
             if(K==L)rank=n-(L-1)*nc
             !call mat_output (nc,nc,EEFF(1:nc,(k-1)*nc+1:k*nc),"loop LL^T EEFF")
         enddo
         do k=1,L-1
           sumk=((k-1)*(k))/2*nc
           !print *," 1FFj=",k, " sumk=",sumk
           !call mat_output (nc,k*nc,FFF(1:nc,sumk+1:sumk+k*nc),"FF")
         enddo


!!!!!!!!!!!   N RHS


    nkRHS=0
            !!      --------------------          RHS               ------------------------------------
    do k = kmin , kmax , kstep   !!  jnrhs  !! loop 
      
      MenLoc=0.0
      nkRHS=nkRHS+1
      jnrhs=(N/kdiv)*k
      !jnrhs=N/10
      !write(*,'(a,I7,a,I5,a,I5,a,I5)') " N= ",N," L =",L," jnrhs=",jnrhs
             !print *, "**********   JNRHS=",jnrhs, "  *******************   K=",k
      opcB=3   !!  for get HS
	
	!!!!  calc allocated memory and exe if less lim
		if(opcB==1)then  !!! RHS
			MenLoc0=((real(LDB)*real(jnrhs))*bytes*3.0d0)/B2GB ; MenLoc=MenLocSET+MenLoc0
                	MenLoc0=((real(m)*real(jnrhs))*bytes*1.0d0)/B2GB ; MenLoc=MenLoc+MenLocSET+MenLoc0

		elseif(opcB==2)then  !! direct gem sym matrix
			MenLoc0=((real(LDB)*real(jnrhs))*bytes*3.0d0)/B2GB ; MenLoc=MenLoc+MenLocSET+MenLoc0

		elseif(opcB==3)then  !! b random
			MenLoc0=((real(LDB)*real(jnrhs))*bytes*1.0d0)/B2GB ; MenLoc=MenLoc+MenLocSET+MenLoc0
		end if
		IF( CHECK )then
			MenLoc0=((real(LDB)*real(jnrhs))*bytes*1.0d0)/B2GB ; MenLoc=MenLoc+MenLocSET+MenLoc0
		ENDIF

		print *, " Total memory whit B, B0 = ", MenLoc, " GB "



      ejecutar=.true.
      if(MenLoc<450.0d0 .and. ejecutar.eqv. .true. )then

      if(opcB==1)then        !!! RHS
      
	    allocate (h(LDB,jnrhs),hteo(LDB,jnrhs),b(LDB,jnrhs),D(m,jnrhs))

 	    !hteo=0.0
        call gem_dmodel(n,jnrhs,hteo(1:n,:))
       	    !print *, "  X*H=d "
	        T1 = omp_get_wtime()
	        !! dados observados
	        !d=matmul(X,hteo)
	    !CALL dGEMM('N','N',M,Jnrhs,N,1.0d0,X,M,hteo,N,0.0d0, d,M )
        !CALL magmaf_dGEMM_cpu('N','N',M,Jnrhs,N,1.0d0,X,M,hteo,N,0.0d0, d,M )
        CALL magmaf_BIG_dgemm_cpu('N','N',M,Jnrhs,N,1.0d0,X,M,hteo,LDB,0.0d0, d,M)
	            !call vec_output (m,d,"dados obs")

	        !print *, "  A^T D=B "  
 	    !b(1:n,:)=matmul(transpose(X),d)
	    !CALL dGEMM('T','N',N,Jnrhs,M,1.0d0,X,M,D,M,0.0d0, B,N )
        !CALL magmaf_dGEMM_cpu('T','N',N,Jnrhs,M,1.0d0,X,M,D,M,0.0d0, B,N )
        CALL magmaf_BIG_dgemm_cpu('T','N',N,Jnrhs,M,1.0d0,X,M,D,M,0.0d0, B,LDB)
	    !call vec_output (n,b,"X^T d")
	    T2 = omp_get_wtime()

	    !print *, " setup take ", T2-T1," s"
	    
      elseif(opcB==2)then  !! direct gem sym matrix
       
         allocate (h(LDB,jnrhs),hteo(LDB,jnrhs),b(LDB,jnrhs))
             
 	     hteo=0.0
         call gem_dmodel(n,jnrhs,hteo(1:n,:))
              !call mat_output (n,jnrhs,hteo,"hteo")
	     !print *, "  B=A*H " 
        !b=0.0d0
        !b(1:n,:)=matmul(A0(1:n,1:n),hteo)
        !CALL magmaf_dGEMM_cpu('N','N',N,Jnrhs,N,1.0d0,A0,N,hteo,N,0.0d0, B,N )
        CALL magmaf_BIG_dgemm_cpu('N','N',N,Jnrhs,N,1.0d0,A0,LDA,hteo,LDB,0.0d0, B,LDB)

      elseif(opcB==3)then  !! b random
     
        allocate (B(LDB,jnrhs) ) !!     
             !print *, "  random_number(B) "       
        call random_number(B)
           
      end if

      
      !!!  anterior pos de ejecutar e condicional size

        IF( CHECK )then
            allocate (B0(LDB,jnrhs) )
            !print *, " start copy B to B0, for verification"
           !B0=B
           call dlacpy("f",n,jnrhs,B,LDB,B0,LDB	)
           !print *, " end   copy B to B0, for verification"
        ENDIF  !!  copy B in B0
    
        dMen_trs_work=0
        MenLoc_rhs(nkRHS)=MenLoc

        !h=0.0d0 
            !print *, ""
            print *, "+++++ Star   magmaf_Outcore_PoTRS_d_v5 nrhs=",jnrhs
            !T4=MPI_WTIME()
            T4 =omp_get_wtime()
        !call solve_BHSNEs2019_magmaPoTRS_d_v5(n,Nc,L,A,b,jnrhs,EEFF,FFF)
        call magmaf_Outcore_PoTRS_d_v5(n,Nc,L,jnrhs,A,LDA,b,LDB,EEFF,FFF,maxNRHs)
        
            max_rhs(nkRHS)=maxNRHS
            
            !T5=MPI_WTIME() 
            T5 =omp_get_wtime()
            !print *, "+++++ finish magmaf_Outcore_PoTRS_d_v5"
            !print *, ""
            !print *, " time magmaPoTRS  ", T5-T4," s"
            tslv=T5-t4
            TIME_ITER_lvS_RHS(LNCH,nkRHS ) = tslv


!     System   A*x=b
!     Compute residual ||A * X  - B|| / ( ||X|| * ||A|| * eps * N )
        IF( CHECK )then
            EPS = dLAMCH('E' )
      		    !print *, "EPS= ",EPS
            normtyp='I'  !!  1, I
                !print *, "normtyp= '",normtyp,"'"
            !ANORM = dLANGE( normtyp, n, n, A0, LDA, WORK2n )
            if ( L == 1) then 
               ANORM = dLANSY( normtyp, 'L', n, A0, LDA, WORK2n )
            else  !! L>1
               ANORM = dLANSY( normtyp,'L', n, A, LDA, WORK2n )
            endif

      		    !print *, "ANORM= ",ANORM
            dBnorm = dLANGE( normtyp, n, jnrhs, B0, LDB, WORK2n )
                !print *, 'Bnorm',dBnorm
            XNORM = dLANGE( normtyp, n, jnrhs, B, LDB, WORK2n )
      		    !print *, "XNORM= ",XNORM
                !T1=MPI_WTIME()
                T1 =omp_get_wtime()
                !call mat_output (n,n,A0,"A0")
            !CALL magmaf_BIG_dgemm_cpu('N','N',N,jNRHS,N,1.0d0,A0,LDA,B,LDB,-1.0d0, B0,LDB)
            if ( L == 1) then
               CALL magmaf_BIG_dsymm_cpu('L','L',N,jNRHS,1.0d0,A0,LDA,B,LDB,-1.0d0, B0,LDB)
            else  !! L>1
               CALL magmaf_BIG_dsymm_cpu('L','L',N,jNRHS,1.0d0,A,LDA,B,LDB,-1.0d0, B0,LDB)
            endif
            
                !call mat_output (n,jNRHS,B0,"B0")
                !T2=MPI_WTIME()
                T2 =omp_get_wtime()
                tmm=t2-t1
                !print *, "time dsymm= ",tmm
                TIME_ITER_mm( LNCH ) = tmm
            AXBNORM = dLANGE( normtyp, n, jnrhs, B0, LDB, WORK2n )
      		    !print *, "AXBNORM= ",AXBNORM
            derror = AXBNORM / ((ANORM*XNORM + dBnorm) * dble(n))
            derror_RHS(LNCH,nkRHS)=derror
            RESID = AXBNORM / ( ANORM*XNORM*EPS*dble( n ) )
            RESID_RHS(LNCH,nkRHS)=RESID  
            SM_abs_AX_B(LNCH,nkRHS)=sum(abs(B0(1:n,:)))
                !WRITE(*,'(A,20ES11.3)') "(A * X  - B)(1:10 ,1)", B0(1:10,1)
                !WRITE(*,'(A,20ES11.3)') "(A * X  - B)(N-9:N,1)", B0(N-9:N,1)
       		    !print *, "RESID= ",RESID
                !print *, "sum(abs(A * X  - B))=",sum(abs(B0(1:N,:))), "sum(A * X  - B)=",sum(B0(1:N,:))




            tol = 60 * eps
            if (derror < tol .or. RESID <  THRESH ) then
                okay = "ok"
            else
                okay = "FAILED"
            endif
            PASSED = 1
            IF( RESID.GT.THRESH )THEN
                        PASSED = 0
            ENDIF
            TESTS_PASSED = TESTS_PASSED + PASSED
            TESTS_FAILED = TESTS_FAILED + ( 1 - PASSED )
        ELSE
            RESID = 0.0E+0
            !TESTS_WOCHK = TESTS_WOCHK + 1
            PASSED = 0
        ENDIF       !!CHECK
      endif   !!   lim mem RAM

        if(opcB==1 .or. opcB==2)deallocate(B,Hteo,H)
        if(opcB==1)deallocate (D)
        if(opcB==3)then
	    if(allocated(B)) deallocate(B)
        endif
	IF( CHECK )then
	   if(allocated(B0)) deallocate(B0)
        ENDIF

    enddo !! jrhs
    
   if(allocated(workLVF))deallocate(workLVF)
   if(allocated(X))deallocate(X)
   if(allocated(A))deallocate(A)
   if(allocated(A0)) deallocate(A0)
   if(allocated(EEFF))deallocate(EEFF)
   if(allocated(FFF))deallocate(FFF)
   deallocate ( a_seed )

   ENDDO !!    LUNCH
                  !!  fac 
        CALL CALC_AVERAGE_TIME( TIME_ITER_lvF, ITERATIONS, tfac )
        
        CALL CALC_AVERAGE_TIME( TIME_ITER_mm, ITERATIONS, tmm )

	CALL CALC_AVERAGE_TIME( TIME_ITER_MMsyrk, ITERATIONS, tmmsyrk )
        
	!! fact
         call get_Ofac(L,nc,Ofac)
         dOfac=dble(Ofac)
         GflpSfac=dOfac/(tfac*1.0D9)
        
	!! MM syrk
        dOMMsyrk=2*dble(n)*dble(n)*dble(n)
        GflpMMsyrk=dOMMsyrk/(tmmsyrk*1.0D9)
        
        !!      --------------------          RESULT              ------------------------------------
         nkRHS=0
         do k=kmin,kmax,kstep   !!  jnrhs  !! loop 
            nkRHS=nkRHS+1
            jnrhs=(N/kdiv)*k
            
            CALL CALC_AVERAGE_TIME( TIME_ITER_lvS_rhs(:,nkrhs), ITERATIONS, tslv )
            
           
	        Osol=2*(INT8(L)**2)*(INT8(Nc)**2)-INT8(L)*INT8(Nc) 
            Osol=Osol*INT8(jnrhs)
            OMM=2*int8(n)*int8(n)*int8(jnrhs)

            dMen_trs=dMen_trf+((dble(LDA)*dble(jnrhs))*dble(bytes))/(1024.0**2)
            
            dOsol=2*(dble(L)**2)*(dble(Nc)**2)-dble(L)*dble(Nc) 
            dOsol=dOsol*dble(jnrhs)    !! spotrs= nrhs(2n**2)
            dOMM=2*dble(n)*dble(n)*dble(jnrhs)                !! sgemm= 2*m*k*n

            if( CHECK ) then
               RESID=SUM(RESID_RHS(1:ITERATIONS,nkRHS))/dble(ITERATIONS)
               derror=SUM(derror_RHS(1:ITERATIONS,nkRHS))/dble(ITERATIONS)
               SM_abs=SUM(SM_abs_AX_B(1:ITERATIONS,nkRHS))/dble(ITERATIONS)
            endif
            
	        GflpSslv=dOsol/(tslv*1.0D9)
	        GflpST  =(dOfac+dOsol)/((tfac+tslv)*1.0D9)
            GflpMM  =dOMM/(tmm*1.0D9)


             !print *, "write info"
            !write(*,512) "n","jNrhs","LDA","it","PASED","RESID","tfac","tslv","tsolT","GflpSfac","GflpSslv","GflpST","GflpMM"
            write(NOUT,414) n,jNrhs,LDA,Nc,max_rhs(nkRHS),L,iterations,TESTS_PASSED,derror,RESID,SM_abs,tfac,tslv,tfac+tslv,&
                  GflpSfac,GflpSslv,GflpST,GflpMMsyrk,GflpMM,dmen_trf+dMen_trf_work,dmen_trs+dMen_trs_work,MenLoc_rhs(nkRHS)
                            
                            
              !write(*,513)"m","n","nc","jNrhs","RESID","tfac","tslv","tsolT","GflpSfac","GflpSslv","GflpST","GflpMM","Ofac","Osol"
              !write(*,413) m,n,nc,jNrhs,RESID,tfac,tslv,tfac+tslv,GflpSfac,GflpSslv,GflpST,GflpMM,Ofac,Osol

    !! Print header
        !print '(a5, "  ", a5, "  ", a12, "  ", a12, "  ", a12, "  ", a12, "  ", a7)', &
        !  "n", "nrhs",  "error", "RESID", "trf=Gflop/s","trs=Gflop/s", "status"

        !print '(i5, "  ", i5, "  ", "  ", es13.4, "  ", es13.4, "  ", f12.4, "  ", f12.4, "  ", a)', &
        !      n, nrhs, derror, RESID, GflpSfac,GflpSslv, okay
              
              
              
         enddo   !! nrhs result
        WRITE( NOUT, '(A)'    ) ''                    
            DEALLOCATE( TIME_ITER_lvf )
            DEALLOCATE( TIME_ITER_lvs )
            DEALLOCATE( TIME_ITER_lvs_rhs )
            DEALLOCATE( TIME_ITER_mm )
	    DEALLOCATE( TIME_ITER_mmsyrk )
            DEALLOCATE( RESID_RHS)
            DEALLOCATE( derror_RHS)
            DEALLOCATE(SM_abs_AX_B)
            DEALLOCATE( MenLoc_rhs)
            DEALLOCATE( max_rhs)
             deallocate (work2n)
  

  enddo !!  N, loop for N


 
413   format((4I7,8ES13.4,2I18))  !! 
513   format((4(a7),8(a13),2(a18)))  !! 
414   format((5I8,3I4,6ES11.3,8F9.2))  !! 
514   format((5(a8),3(a4),6(a11),8(a9)))  !! 




endif  !!ejecutar

!print *,""
!print '(4I7,3ES13.4)',m,n,Nc,L,sMR,sMR*1024.0,sMR*1024.0**2

    finish=omp_get_wtime()
    
111  CONTINUE
    print *, "TOTAL time ",finish-start

    call magmaf_finalize()

end program main

subroutine get_lwork_magma_dpoLVF(L,nc,lwork)


implicit none
integer*8   :: L,nc,lwork

    if(L>3)then
        lwork=nc*nc*(3*L-2)
    else
        lwork=nc
    endif

end subroutine get_lwork_magma_dpoLVF




subroutine gem_dmodel(n,nrhs,h)

   implicit none
   integer :: n, nrhs,k,i
   real (kind=8)   :: h(n,nrhs)
   real (kind=8)   :: r

   h=0.0d0
   h(1:n,1)=1.0d0  !! first model, one vector
   if(nrhs>1)then
      do k=2,nrhs
         if(k==2)then
            h(n/2,2)=1000.0d0  !! delta
         elseif(k==3)then
            CALL RANDOM_NUMBER(h(:,3))  !! RANDOM
         else
              do i=1,n
                    CALL RANDOM_NUMBER(r)
                    h(i,k)=dble(i)-dble(k)/10.d0+5.d0+r
             enddo
        endif
      enddo
   endif

end subroutine gem_dmodel



subroutine gem_blkTri(nc,L,BHC)

   implicit none
   integer :: nc, L,k!,i
   real (kind=8)   :: BHC(nc*L,nc*L)

  !! low 
   do k=1,L-1
       BHC((k+1)*nc+1:L*nc,(k-1)*nc+1:k*nc)=0.0d0
   enddo

  !! uper
   do k=1,L-1
       BHC((k-1)*nc+1:k*nc,(k+1)*nc+1:L*nc)=0.0d0
   enddo

end subroutine gem_blkTri






    subroutine get_Ofac(L,nc,Ofac)


    implicit none
    integer   :: L,nc
    integer*8 :: Ofac, Ocho,Ofb,sm,Ommpm


    !Ocho=(INT8(nc)*(INT8(nc)+1)*(2*INT8(nc)+1))/6   !!! igual
    Ocho=( 2*INT8(nc)*INT8(nc)*INT8(nc) + 3*INT8(nc)*INT8(nc) + INT8(nc) )/6 !! lawn41.ps

    !Ofb=( 2*INT8(nc)**2-INT8(nc) ) * INT8(nc)
    Ofb=2*INT8(nc)**3  !!  !!law41

    !Ommpm=INT8(nc)**3+INT8(nc)**2   !! only Low part
    Ommpm=2*INT8(nc)**3            ! full

    sm=Ocho+Ofb
    if(L==1)then
	Ofac=Ocho

    elseif(L==2)then
	Ofac=sm+ Ommpm  + 2*Ocho

	!elseif(L==3)then

    else
        ! row 10 delH
	Ofac= (( INT8(L)*(INT8(L)-1)*(INT8(L)-2)*2*INT8(nc)**3))/6

	!!   iHjj and iFFjj
	! row 13 iHjj
	!Ofac=Ofac+ ( (INT8(L)-1)*(INT8(L)-2)*sm)/2
	! row 16 iFjj
	!Ofac=Ofac+ ( (INT8(L))*(INT8(L)-1)*sm)/2
	Ofac=Ofac + INT8(L)*Ocho   !!  fac of all Ajj
        Ofac=Ofac + INT8(L-1)*Ofb    !! slv of 1FFjj
        Ofac=Ofac + INT8(L-2)*Ofb    !! slv of 1HHjj



	! row 14  iEHj
	Ofac=Ofac+  (    (INT8(L)-1)  *  (INT8(L)-2)  *  Ommpm  )/2

	! row 17 iEFj
	Ofac=Ofac+ ( INT8(L)*(INT8(L)-1)*(Ommpm))/2

	! row 18  recLev FFj
	Ofac=Ofac+ ( INT8(L)*(INT8(L)-1)*(INT8(L)-2)*INT8(nc)**3)/3
	! row 20  recLev HHj
	Ofac=Ofac+ ( (INT8(L)-1)*(INT8(L)-2)*(INT8(L)-3)*2*INT8(nc)**3)/6

	! row 25  LL^T (1EFj)
	Ofac=Ofac+ ( INT8(L-1)*Ocho )   !! exepc A11

    endif

end subroutine

