
module Sub_CPU_Interfaces
use magma
implicit none
contains


subroutine magmaf_dgemm_cpu(TRANSA,TRANSB,M,N,K,ALPHA,hA,LDA,hB,LDB,BETA,hC,LDC)  !!   work

!     .. Scalar Arguments ..
       DOUBLE PRECISION ALPHA,BETA
       INTEGER LDA,LDB,LDC,M,N,K
       CHARACTER TRANSA,TRANSB
!     .. Array Arguments ..CPU
       DOUBLE PRECISION hA(LDA,*),hB(LDB,*),hC(LDC,*)
!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dA, dB, dC
       magma_devptr_t          :: queue
       INTEGER LDDA,LDDB,LDDC,dev
!  =====================================================================
!     .. Local Scalars ..
       INTEGER INFO,NCOLA,NROWA,NROWB,NCOLB
       LOGICAL NOTA,NOTB
!     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame
!     .. Parameters ..
       double precision                    ::  zero,one
       parameter                       ( zero=0.0d0,one=1.0d+0)

       print *, "Start magmaf_dgemm_cpu"

!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
       nota = lsame(transa,'N')
       notb = lsame(transb,'N')
       IF (nota) THEN
           nrowa = m
           ncola = k
       ELSE
           nrowa = k
           ncola = m
       END IF
       IF (notb) THEN
           nrowb = k
           ncolb = n
       ELSE
           nrowb = n
           ncolb = k
       END IF


!    Quick return if possible.
!
       IF ((m.EQ.0) .OR. (n.EQ.0) .OR.(((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one)))then
           print *, "Quick return magmaf_dgemm_cpu"
           RETURN
       ENDIF

    lddc = ceiling(real(m)/32)*32
    ldda = ceiling(real(nrowa)/32)*32
    lddb = ceiling(real(nrowb)/32)*32
    info = 0

    dev=0
    call magmaf_queue_create( dev, queue )
      !!  C := alpha*op( A )*op( B ) + beta*C,
      !! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*ncola )
        if (info .ne. 0) then
            print *, "dgemm_cpu Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dB, lddb*ncolb )
        if (info .ne. 0) then
            print *, "dgemm_cpu Error: magmaf_dmalloc( dB  ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dC, lddc*n )
        if (info .ne. 0) then
            print *, "dgemm_cpu Error: magmaf_dmalloc( dC  ) failed, info = ", info
            goto 1000
        endif
           !! C := alpha*op( A )*op( B ) + beta*C,
    	!! dA = hA
        call magmaf_dsetmatrix( nrowa, ncola, hA, lda, dA, ldda, queue )
	!! dB = hB
        call magmaf_dsetmatrix( nrowB, ncolb, hB, ldB, dB, lddB, queue )
        if(BETA/=zero)then
		!! dC = hC
        	call magmaf_dsetmatrix( m, n, hC, ldC, dC, lddC, queue )
         endif
	!call magmablasf_dgemm(TRANSA,TRANSB,m,n,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )
	call magmaf_dgemm(TRANSA,TRANSB,m,n,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddC,queue )

 	!! hC = dC
        call magmaf_dgetmatrix( m, n, dC, lddC, hC, ldC, queue )

   
!! cleanup:
1000 continue
  !! Free GPU memory
        info = magmaf_free( dA )
        if (info .ne. 0) then
            print *, 'magmaf_dgemm_cpu Error: magmaf_free( dA ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dB )
        if (info .ne. 0) then
            print *, 'magmaf_dgemm_cpu Error: magmaf_free( dB ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dC )
        if (info .ne. 0) then
            print *, 'magmaf_dgemm_cpu Error: magmaf_free( dC ) failed, info = ', info
            stop
        endif
    call magmaf_queue_destroy( queue )
    print *, "END magmaf_dgemm_cpu"
end subroutine


subroutine magmablasf_dgemm_cpu(TRANSA,TRANSB,M,N,K,ALPHA,hA,LDA,hB,LDB,BETA,hC,LDC)

!     .. Scalar Arguments ..
       DOUBLE PRECISION ALPHA,BETA
       INTEGER LDA,LDB,LDC,M,N,K
       CHARACTER TRANSA,TRANSB
!     .. Array Arguments ..CPU
       DOUBLE PRECISION hA(LDA,*),hB(LDB,*),hC(LDC,*)
!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dA, dB, dC
       magma_devptr_t          :: queue
       INTEGER LDDA,LDDB,LDDC,dev
!  =====================================================================
!     .. Local Scalars ..
       INTEGER INFO,NCOLA,NROWA,NROWB,NCOLB
       LOGICAL NOTA,NOTB
!     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame
!     .. Parameters ..
       double precision                    ::  zero,one
       parameter                       ( zero=0.0d0,one=1.0d+0)

       print *, "Start magmaf_dgemm_cpu"

!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
       nota = lsame(transa,'N')
       notb = lsame(transb,'N')
       IF (nota) THEN
           nrowa = m
           ncola = k
       ELSE
           nrowa = k
           ncola = m
       END IF
       IF (notb) THEN
           nrowb = k
           ncolb = n
       ELSE
           nrowb = n
           ncolb = k
       END IF


!    Quick return if possible.
!
       IF ((m.EQ.0) .OR. (n.EQ.0) .OR.(((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one)))then
           print *, "Quick return magmaf_dgemm_cpu"
           RETURN
       ENDIF

    lddc = ceiling(real(m)/32)*32
    ldda = ceiling(real(nrowa)/32)*32
    lddb = ceiling(real(nrowb)/32)*32
    info = 0

    dev=0
    call magmaf_queue_create( dev, queue )
      !!  C := alpha*op( A )*op( B ) + beta*C,
      !! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*ncola )
        if (info .ne. 0) then
            print *, "dgemm_cpu Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dB, lddb*ncolb )
        if (info .ne. 0) then
            print *, "dgemm_cpu Error: magmaf_dmalloc( dB  ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dC, lddc*n )
        if (info .ne. 0) then
            print *, "dgemm_cpu Error: magmaf_dmalloc( dC  ) failed, info = ", info
            goto 1000
        endif
           !! C := alpha*op( A )*op( B ) + beta*C,
    	!! dA = hA
        call magmaf_dsetmatrix( nrowa, ncola, hA, lda, dA, ldda, queue )
	!! dB = hB
        call magmaf_dsetmatrix( nrowB, ncolb, hB, ldB, dB, lddB, queue )
        if(BETA/=zero)then
		!! dC = hC
        	call magmaf_dsetmatrix( m, n, hC, ldC, dC, lddC, queue )
         endif
	!call magmablasf_dgemm(TRANSA,TRANSB,m,n,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )
	call magmablasf_dgemm(TRANSA,TRANSB,m,n,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddC,queue )

 	!! hC = dC
        call magmaf_dgetmatrix( m, n, dC, lddC, hC, ldC, queue )

   
!! cleanup:
1000 continue
  !! Free GPU memory
        info = magmaf_free( dA )
        if (info .ne. 0) then
            print *, 'Error: magmaf_free( dA ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dB )
        if (info .ne. 0) then
            print *, 'Error: magmaf_free( dB ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dC )
        if (info .ne. 0) then
            print *, 'Error: magmaf_free( dX ) failed, info = ', info
            stop
        endif
    call magmaf_queue_destroy( queue );

end subroutine


subroutine magmaf_BIG_dgemm_cpu(TRANSA,TRANSB,M,N,K,ALPHA,hA,LDA,hB,LDB,BETA,hC,LDC)
!        C(m,n) := alpha*op( A(m,k) )*op( B(k,n) ) + beta*C(m,n)
!     .. Scalar Arguments ..
       DOUBLE PRECISION ALPHA,BETA
       INTEGER LDA,LDB,LDC,M,N,K
       CHARACTER TRANSA,TRANSB
!     .. Array Arguments ..CPU
       DOUBLE PRECISION hA(LDA,*),hB(LDB,*),hC(LDC,*)
!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dA, dB, dC
       magma_devptr_t          :: queue
       INTEGER LDDA,LDDB,LDDC,dev
!  =====================================================================
!     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame
!     .. Local Scalars ..
       LOGICAL NOTA,NOTB
       INTEGER INFO
       INTEGER*8 :: mem_size
       !INTEGER :: h_groups, v_groups, height, width
       double precision                    ::  zero,ONE
       parameter                       ( zero=0.0d0,ONE=1.0D0)
       double precision                    :: LimMem,totM, MaxMem
       integer                             :: j,i,jb,ib,NCOLA,NROWA,NROWB,ncolb
       integer                             :: nb,mb !,maxnb
!  =====================================================================

     !print *, "+ Star magmaf_BIG_dgemm_cpu"
    dev=0
    call magmaf_queue_create( dev, queue )

    nota = lsame(transa,'N')
    notb = lsame(transb,'N')

!     Quick return if possible.
       IF ((m.EQ.0) .OR. (n.EQ.0) .OR. (((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one)))then
           !print *, "Quick return magmaf_BIG_dgemm_cpu"
           RETURN
       ENDIF

    mem_size=magmaf_mem_size( queue )    !!    mem_size=11 821 973 504  integer
    !mem_size=int8(0.5d0*(1024**3))

    LimMem=dble(mem_size)*0.90d0 !   GB  90% of mex men
    MaxMem=( dble(m)*dble(n) + dble(m)*dble(k) +dble(k)*dble(n) )*8.0d0 
    !print *, "memory full problem= ",MaxMem, "limit memory device=",LimMem

    !print *, "mem_size=",mem_size
     
    !call matrixPartitioner(M,  N,  K, mem_size, h_groups, v_groups,height, width)

    ! print *, "h_groups=",h_groups, "v_groups=",v_groups,"height=",height, "width=",width

     !print *,"Partitioner Mem in dev=", ((dble(width)*dble(height)+dble(height)*dble(k)+dble(width)*dble(k))*8.0)/1e9 ," GB"

if(maxmem<LimMem)then
   !! full problem
      !print *, "Full problem"

       IF (nota) THEN
           nrowa = m
           ncola = k
       ELSE
           nrowa = k
           ncola = m
       END IF
       IF (notb) THEN
           nrowb = k
           ncolb = n
       ELSE
           nrowb = n
           ncolb = k
       END IF


    lddc = ceiling(real(m)/32)*32
    ldda = ceiling(real(nrowa)/32)*32
    lddb = ceiling(real(nrowb)/32)*32
    info = 0



      !! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*ncola )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dB, lddb*ncolb )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dC, lddc*n )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dC  ) failed, info = ", info
            goto 1000
        endif

    	!! dA = hA
        !call magmaf_dsetmatrix( nrowa, ncola, hA, lda, dA, ldda, queue )
        call magmaf_dsetmatrix_async( nrowa, ncola, hA, lda, dA, ldda, queue )
	!! dB = hB
        !call magmaf_dsetmatrix( nrowb, ncolb, hB, ldB, dB, lddB, queue )
        call magmaf_dsetmatrix_async( nrowb, ncolb, hB, ldB, dB, lddB, queue )
        if(BETA/=zero)then
		!! dC = hC
        	!call magmaf_dsetmatrix( m, n, hC, ldC, dC, lddC, queue )
                call magmaf_dsetmatrix_async( m, n, hC, ldC, dC, lddC, queue )
         endif
        call magmaf_queue_sync(queue)
	!call magmablasf_dgemm(TRANSA,TRANSB,m,n,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )
        call magmaf_dgemm(TRANSA,TRANSB,m,n,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )

 	!! hC = dC
        call magmaf_dgetmatrix( m, n, dC, lddC, hC, ldC, queue )
        !call magmaf_dgetmatrix_async( m, n, dC, lddC, hC, ldC, queue )

ELSE   !!  	by block

    if( (dble(m)*dble(512) + dble(m)*dble(k) +dble(k)*dble(512)) *8.0d0 <LimMem)then   !!  by colums block
        !!   full A in device and nb=>512,  =>> C(m,nb)=C(m,nb)+a(m,k)*b(k,nb)

       IF (notb) THEN    
	     !!  C := alpha*A*B + beta*C.  or C := alpha*A**T*B + beta*C.
          nb=(LimMem/8-m*k)/(m+k)-127
          nb = ceiling(real(nb)/128)*128
          !print *, "full A in device and nb=",nb
          lddc = ceiling(real(m)/32)*32
          ldda = ceiling(real(m)/32)*32
          lddb = ceiling(real(k)/32)*32
       
          !print *,"Mem blk=", ((dble(m)*dble(nb)+dble(m)*dble(k)+dble(k)*dble(nb))*8.0)/1e9 ," GB"       
          !print *,"Mem in device=", ((dble(lddc)*dble(nb)+dble(ldda)*dble(k)+dble(lddb)*dble(nb))*8.0)/1e9 ," GB"
             !! Allocate GPU memory
          info = magmaf_dmalloc( dA, ldda*k )
          if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
          endif

          info = magmaf_dmalloc( dB, lddb*nb )
          if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            goto 1000
          endif

          info = magmaf_dmalloc( dC, lddc*nb )
          if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dC  ) failed, info = ", info
            goto 1000
          endif
          !print *, "end Allocate GPU memory"
	  !!  num blk

	  !! dA = hA   !! full A
          call magmaf_dsetmatrix( m, k, hA, lda, dA, ldda, queue )


          DO j=1,n,nb
                jb = min( nb, n-j+1 )
		!! dB = hB
        	!call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                call magmaf_dsetmatrix_async( k, jb, hB(1,j), ldB, dB, lddB, queue )

                		!print *, "j=",j,"jb=",jb
        	if(BETA/=zero)then
			!! dC = hC
        		!call magmaf_dsetmatrix( m, jb, hC(1,j), ldC, dC, lddC, queue )
                        call magmaf_dsetmatrix_async( m, jb, hC(1,j), ldC, dC, lddC, queue )
                 endif
                call magmaf_queue_sync(queue) 	
		!call magmablasf_dgemm(TRANSA,TRANSB,m,n,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )
		call magmaf_dgemm(TRANSA,TRANSB,m,jb,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddC,queue )
              		!print *, "C(1:m,j:j+jb)=C(1:m",j,":",j+jb-1,")*A(1:m,1:k)*B(1:k,",j,":",j+jb-1,")"
 		!! hC = dC
        	!call magmaf_dgetmatrix( m, jb, dC, lddC, hC(1,j), ldC, queue )
                call magmaf_dgetmatrix_async( m, jb, dC, lddC, hC(1,j), ldC, queue )


          ENDDO
       ELSE
          !!  C := alpha*A*B**T + beta*C.  or C := alpha*A**T*B**T + beta*C.
             PRINT *, " NOT IMPLEMENT YET, full A in device"
             PRINT *, " C := alpha*A*B**T + beta*C.  or C := alpha*A**T*B**T + beta*C."
       ENDIF  
    
    ELSE   !!  BY BLK(mb,nb)
        !!   cal C by C(mb,nb)=C(mb,nb)+op(a())*op(b())
       IF (nota) THEN
           nrowa = m
           ncola = k
       ELSE
           nrowa = k
           ncola = m
       END IF
       IF (notb) THEN
           nrowb = k
       ELSE
           nrowb = n
       END IF

       call magmaf_get_dgemm_mb_nb(limMem,m,n,k,mb,nb)
       !call matrixPartitioner(M,  N,  K, int8(limMem), h_groups, v_groups,mb, nb)
       !print *, "BIG_dgemm; mb", mb," nb",nb

      lddc = ceiling(real(mb)/32)*32
      IF (nota) THEN
          ldda = ceiling(real(mb)/32)*32
       ELSE
          ldda = ceiling(real(k)/32)*32
       END IF
       IF (notb) THEN
           lddb = ceiling(real(k)/32)*32
       ELSE
           lddb = ceiling(real(nb)/32)*32
       END IF
      

      !! Allocate GPU memory  
      info = magmaf_dmalloc( dC, lddc*nb )
      if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dC  ) failed, info = ", info
            goto 1000
      endif



      IF (notb) THEN


           IF (nota) THEN
		          ! C(mb,nb) := alpha*A(mb,k)*B(k,nb) + beta*C(mb,nb). 

       		!ldda = ceiling(real(mb)/32)*32
       		!lddb = ceiling(real(k)/32)*32
       		info = 0
       		!print *," Mem in dev=", ((dble(mb)*dble(nb)+dble(mb)*dble(k)+dble(nb)*dble(k))*8.0)/1e9 ," GB"
       		totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(k)+dble(lddb)*dble(nb))*8.0D0)/1e9
       		!print *, "BIG_dgemm; Memory in device=", totM," GB"

		!!  num blk
      		info = magmaf_dmalloc( dA, ldda*k )
      		if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            		goto 1000
      		endif
      		info = magmaf_dmalloc( dB, lddb*nb )
      		if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            		goto 1000
      		endif
      		!print *, "end Allocate GPU memory"

        	DO j=1,n,nb
                	jb = min( nb, n-j+1 )
			!! dB = hB
        		!call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                	call magmaf_dsetmatrix_async( k, jb, hB(1,j), ldB, dB, lddB, queue )
           		DO i=1,m,mb
    				ib = min( mb, m-i+1 )
				!! dA = hA
        			!call magmaf_dsetmatrix( ib, k, hA(i,1), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( ib, k, hA(i,1), lda, dA, ldda, queue )
                				!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib
        			if(BETA/=zero)then
					!! dC = hC
        				call magmaf_dsetmatrix( ib, jb, hC(i,j), ldC, dC, lddC, queue )
                 		endif
                		call magmaf_queue_sync(queue) 	
				!call magmablasf_dgemm(TRANSA,TRANSB,m,n,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )
				call magmaf_dgemm(TRANSA,TRANSB,ib,jb,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"
 				!! hC = dC
        			!call magmaf_dgetmatrix( ib, jb, dC, lddC, hC(i,j), ldC, queue )
                		call magmaf_dgetmatrix_async( ib, jb, dC, lddC, hC(i,j), ldC, queue )

             		ENDDO

          	ENDDO
           ELSE  !   nota
	          !  C := alpha*A**T*B + beta*C
                   !print *,"C := alpha*A**T*B + beta*C"
       		!ldda = ceiling(real(k)/32)*32
       		!lddb = ceiling(real(k)/32)*32
       		info = 0
       		!print *," Mem in dev=", ((dble(mb)*dble(nb)+dble(mb)*dble(k)+dble(nb)*dble(k))*8.0)/1e9 ," GB"
       		totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(mb)+dble(lddb)*dble(nb))*8.0D0)/1e9
       		!print *, "BIG_dgemm; Memory in device=", totM," GB"


      		!! Allocate GPU memory
        	info = magmaf_dmalloc( dA, ldda*mb )
        	if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            		goto 1000
        	endif

        	info = magmaf_dmalloc( dB, lddb*nb )
        	if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            		goto 1000
        	endif

        	!print *, "end Allocate GPU memory"
		!!  num blk


        	DO j=1,n,nb
                	jb = min( nb, n-j+1 )
			!! dB = hB
        		!call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                	call magmaf_dsetmatrix_async( k, jb, hB(1,j), ldB, dB, lddB, queue )
           		DO i=1,m,mb
    				ib = min( mb, m-i+1 )
				!! dA = hA
        			!call magmaf_dsetmatrix( k,ib, hA(1,i), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( k,ib, hA(1,i), lda, dA, ldda, queue )
                		!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib
        			if(BETA/=zero)then
					!! dC = hC
        				call magmaf_dsetmatrix( ib, jb, hC(i,j), ldC, dC, lddC, queue )
                 		endif
                		call magmaf_queue_sync(queue) 	
				!call magmablasf_dgemm(TRANSA,TRANSB,m,n,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )
				call magmaf_dgemm(TRANSA,TRANSB,ib,jb,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"
 				!! hC = dC
        			!call magmaf_dgetmatrix( ib, jb, dC, lddC, hC(i,j), ldC, queue )
                		call magmaf_dgetmatrix_async( ib, jb, dC, lddC, hC(i,j), ldC, queue )

             		ENDDO  ! i

          	ENDDO	! j

           ENDIF

	ELSE  !notB
           IF (nota) THEN
          	! C := alpha*A*B**T + beta*C

       		!ldda = ceiling(real(mb)/32)*32
       		!lddb = ceiling(real(nb)/32)*32
       		info = 0
       		!print *," Mem in dev=", ((dble(mb)*dble(nb)+dble(mb)*dble(k)+dble(nb)*dble(k))*8.0)/1e9 ," GB"
       		totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(k)+dble(lddb)*dble(k))*8.0)/1e9
       		!print *, "BIG_dgemm; Memory in device=", totM," GB"


      		!! Allocate GPU memory
        	info = magmaf_dmalloc( dA, ldda*k )
        	if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            		goto 1000
        	endif

        	info = magmaf_dmalloc( dB, lddb*k )
        	if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            		goto 1000
        	endif

        	!print *, "end Allocate GPU memory"
		!!  num blk


        	DO j=1,n,nb
                	jb = min( nb, n-j+1 )
			!! dB = hB
        		!call magmaf_dsetmatrix( jb, k, hB(j,1), ldB, dB, lddB, queue )
                	call magmaf_dsetmatrix_async( jb, k, hB(j,1), ldB, dB, lddB, queue )
           		DO i=1,m,mb
    				ib = min( mb, m-i+1 )
				!! dA = hA
        			!call magmaf_dsetmatrix( ib, k, hA(i,1), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( ib, k, hA(i,1), lda, dA, ldda, queue )
                				!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib
        			if(BETA/=zero)then
					!! dC = hC
        				call magmaf_dsetmatrix( ib, jb, hC(i,j), ldC, dC, lddC, queue )
                 		endif
                		call magmaf_queue_sync(queue) 	
				!call magmablasf_dgemm(TRANSA,TRANSB,m,n,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )
				call magmaf_dgemm(TRANSA,TRANSB,ib,jb,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"
 				!! hC = dC
        			!call magmaf_dgetmatrix( ib, jb, dC, lddC, hC(i,j), ldC, queue )
                		call magmaf_dgetmatrix_async( ib, jb, dC, lddC, hC(i,j), ldC, queue )

             		ENDDO

          	ENDDO

           ELSE  !  nota
		!!C := alpha*A**T*B**T + beta*C

       		!ldda = ceiling(real(k)/32)*32
       		!lddb = ceiling(real(nb)/32)*32
       		info = 0
       		!print *," Mem in dev=", ((dble(mb)*dble(nb)+dble(mb)*dble(k)+dble(nb)*dble(k))*8.0)/1e9 ," GB"
       		totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(mb)+dble(lddb)*dble(k))*8.0)/1e9
       		!print *, "BIG_dgemm; Memory in device=", totM," GB"


      		!! Allocate GPU memory
        	info = magmaf_dmalloc( dA, ldda*mb )
        	if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            		goto 1000
        	endif

        	info = magmaf_dmalloc( dB, lddb*k )
        	if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            		goto 1000
        	endif
        	!print *, "end Allocate GPU memory"
		!!  num blk


        	DO j=1,n,nb
                	jb = min( nb, n-j+1 )
			!! dB = hB
        		!call magmaf_dsetmatrix( jb, k, hB(j,1), ldB, dB, lddB, queue )
                	call magmaf_dsetmatrix_async( jb, k, hB(j,1), ldB, dB, lddB, queue )
           		DO i=1,m,mb
    				ib = min( mb, m-i+1 )
				!! dA = hA
        			!call magmaf_dsetmatrix( k, ib, hA(1,i), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( k, ib, hA(1,i), lda, dA, ldda, queue )
                				!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib
        			if(BETA/=zero)then
					!! dC = hC
        				call magmaf_dsetmatrix( ib, jb, hC(i,j), ldC, dC, lddC, queue )
                 		endif
                		call magmaf_queue_sync(queue) 	
				!call magmablasf_dgemm(TRANSA,TRANSB,m,n,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )
				call magmaf_dgemm(TRANSA,TRANSB,ib,jb,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"
 				!! hC = dC
        			!call magmaf_dgetmatrix( ib, jb, dC, lddC, hC(i,j), ldC, queue )
                		call magmaf_dgetmatrix_async( ib, jb, dC, lddC, hC(i,j), ldC, queue )

             		ENDDO

          	ENDDO

           END IF !  nota
       END IF  !notB

    ENDIF  !! full A in device

ENDIF  !!  maxmem<LimMem





 call magmaf_queue_sync(queue)
   
!! cleanup:
1000 continue
  !! Free GPU memory
        info = magmaf_free( dA )
        if (info .ne. 0) then
            print *, 'magmaf_BIG_dgemm_cpu Error: magmaf_free( dA ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dB )
        if (info .ne. 0) then
            print *, 'magmaf_BIG_dgemm_cpu Error: magmaf_free( dB ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dC )
        if (info .ne. 0) then
            print *, 'magmaf_BIG_dgemm_cpu Error: magmaf_free( dC ) failed, info = ', info
            stop
        endif
    call magmaf_queue_destroy( queue )

    !print *, "+ END magmaf_BIG_dgemm_cpu"

end subroutine


!  =====================================================================
       SUBROUTINE magmaf_BIG_dpotrf_cpu( UPLO, N, A, LDA, INFO )

!  
 !     .. Scalar Arguments ..
       CHARACTER          UPLO
       INTEGER            INFO, LDA, N
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION   A( LDA,*   )
 !     ..
!  
 !  =====================================================================

 !     .. Parameters ..
       DOUBLE PRECISION   ONE, neg_one
       parameter( one = 1.0d0, neg_one=-1.0d0 )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            UPPER
       INTEGER            J, JB, NB
 !     ..
 !     .. External Functions ..
       LOGICAL            LSAME
       INTEGER            ILAENV
       EXTERNAL           lsame, ilaenv
 !     ..
 !     .. External Subroutines ..
       !EXTERNAL           magmaf_BIG_dgemm_cpu, magmaf_dsyrk, magmaf_dtrsm, xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max, min
 !     ..
 !     .. Executable Statements ..
!  
 !     Test the input parameters.
!  
      print *, "* USING magmaf_BIG_dpotrf_cpu"

       info = 0
       upper = lsame( uplo, 'U' )
       IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( lda.LT.max( 1, n ) ) THEN
          info = -4
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'DPOTRF', -info )
          RETURN
       END IF
!  
 !     Quick return if possible
!  
       IF( n.EQ.0 ) &
         RETURN
!  
 !     Determine the block size for this environment.
!  
       nb = magmaf_get_dpotrf_nb( n );
       print *, "NB=",nb
       IF( nb.LE.1 .OR. nb.GE.n ) THEN
!  
 !        Use lapack code.
!  
          CALL dpotrf( uplo, n, a, lda, info )
       ELSE

!  
 !        /* Use  blocked code. */
!  
          IF( upper ) THEN
!  
 !           Compute the Cholesky factorization A = U**T*U.
!  
             DO 10 j = 1, n, nb
!  
 !              Update and factorize the current diagonal block and test
 !              for non-positive-definiteness.
!  
                jb = min( nb, n-j+1 )
                CALL magmaf_big_dsyrk_cpu( 'Upper', 'Transpose', jb, j-1, -one, a( 1, j ), lda, one, a( j, j ), lda )
                !CALL magmaf_dsyrk_cpu( 'Upper', 'Transpose', jb, j-1, -one, a( 1, j ), lda, one, a( j, j ), lda )
                CALL dpotrf( 'Upper', jb, a( j, j ), lda, info )
                IF( info.NE.0 ) &
                  GO TO 30
                IF( j+jb.LE.n ) THEN
!  
 !                 Compute the current block row.
!  
                   CALL magmaf_BIG_dgemm_cpu( 'Transpose', 'No transpose', jb, n-j-jb+1, &
                              j-1, -one, a( 1, j ), lda, a( 1, j+jb ), 			 &
                              lda, one, a( j, j+jb ), lda )
                   CALL magmaf_dtrsm_CPU( 'Left', 'Upper', 'Transpose', 'Non-unit',		&
                              jb, n-j-jb+1, one, a( j, j ), lda,			&
                              a( j, j+jb ), lda )
                END IF
    10       CONTINUE
!  
          ELSE
             print *, "Compute the Cholesky factorization A = L*L**T."
!  
 !           Compute the Cholesky factorization A = L*L**T.
!  
             DO 20 j = 1, n, nb
!  
 !              Update and factorize the current diagonal block and test
 !              for non-positive-definiteness.
!  
                jb = min( nb, n-j+1 )
                CALL magmaf_big_dsyrk_cpu( 'Lower', 'No transpose', jb, j-1, -one,a( j, 1 ), lda, one, a( j, j ), lda )
                !CALL magmaf_dsyrk_cpu( 'L', 'N', jb, j-1, -one,a( j, 1 ), lda, one, a( j, j ), lda )
                !CALL dsyrk( 'L', 'N', jb, j-1, -one,a( j, 1 ), lda, one, a( j, j ), lda )
                CALL dpotrf( 'Lower', jb, a( j, j ), lda, info )
                IF( info.NE.0 )then
                  PRINT *, "INFO=",INFO, " EXIT magmaf_BIG_dpotrf_cpu"
                  GO TO 30
                ENDIF
                IF( j+jb.LE.n ) THEN
!  
 !                 Compute the current block column.
!                        C := alpha* A * A^T + beta*C,
!  					     !'No transpose','Transpose'
                   CALL magmaf_BIG_dgemm_cpu('N','T',n-j-jb+1,jb,j-1,-one,a(j+jb,1),lda,a(j,1),lda,one,a(j+jb,j),lda) ! funciona 
                   !CALL magmaf_dgemm_cpu('N','T',n-j-jb+1,jb,j-1,-one,a(j+jb,1),lda,a(j,1),lda,one,a(j+jb,j),lda)  ! funciona
                   !CALL dgemm('N','T',n-j-jb+1,jb,j-1,-one,a(j+jb,1),lda,a(j,1),lda,one,a(j+jb,j),lda)  !! funciona
                   CALL magmaf_dtrsm_CPU('Right','Lower','Transpose','Non-unit',n-j-jb+1, jb, one, a( j, j ),lda,a( j+jb, j ),lda )
                   !CALL dtrsm( 'Right', 'Lower', 'Transpose', 'Non-unit',n-j-jb+1, jb, one, a( j, j ), lda,a( j+jb, j ), lda )
                END IF
    20       CONTINUE
          END IF
       END IF
       GO TO 40
!  
    30 CONTINUE
       info = info + j - 1
!  
    40 CONTINUE
       RETURN

      print *, "* END magmaf_BIG_dpotrf_cpu"

!  
 !     End of DPOTRF
!  
       END subroutine magmaf_BIG_dpotrf_cpu




!  =====================================================================
       SUBROUTINE magmaf_BIG_dpotrf_cpu2( UPLO, N, A, LDA, INFO )

!  
 !     .. Scalar Arguments ..
       CHARACTER          UPLO
       INTEGER            INFO, LDA, N
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION   A( LDA,*   )
 !     ..
!  
!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dT,dT2
       magma_devptr_t          :: queue
       INTEGER 		       :: LDDnb,LDDn,dev
!  =====================================================================

 !     .. Parameters ..
       DOUBLE PRECISION   ONE, neg_one
       parameter( one = 1.0d0, neg_one=-1.0d0 )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            UPPER
       INTEGER            J, JB, NB
 !     ..
 !     .. External Functions ..
       LOGICAL            LSAME
       INTEGER            ILAENV
       EXTERNAL           lsame, ilaenv
 !     ..
 !     .. External Subroutines ..
       !EXTERNAL           magmaf_BIG_dgemm_cpu, magmaf_dsyrk, magmaf_dtrsm, xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max, min
 !     ..
 !     .. Executable Statements ..
! 

    dev=0
    call magmaf_queue_create( dev, queue )

 
 !     Test the input parameters.
!  

      print *, "***** USING magmaf_BIG_dpotrf_cpu2"

       info = 0
       upper = lsame( uplo, 'U' )
       IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( lda.LT.max( 1, n ) ) THEN
          info = -4
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'DPOTRF', -info )
          RETURN
       END IF
!  
 !     Quick return if possible
!  
       IF( n.EQ.0 ) &
         RETURN
!  
 !     Determine the block size for this environment.
!  
       nb = magmaf_get_dpotrf_nb( n );
       print *, "NB=",nb
       IF( nb.LE.1 .OR. nb.GE.n ) THEN
!  
 !        Use lapack code.
!  
          CALL dpotrf( uplo, n, a, lda, info )
       ELSE

!  
 !        /* Use  blocked code. */
!  
          !! Allocate GPU memory
        lddnb = ceiling(real(nb)/32)*32
        lddn = ceiling(real(n)/32)*32
        info = 0
        
    
          !! Allocate GPU memory
        info = magmaf_dmalloc( dT, lddnb*nb )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dT ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dT2, lddn*lddnb )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dT2  ) failed, info = ", info
            goto 1000
        endif


          IF( upper ) THEN
!  
 !           Compute the Cholesky factorization A = U**T*U.
!  
             DO 10 j = 1, n, nb
!  
 !              Update and factorize the current diagonal block and test
 !              for non-positive-definiteness.
!  
                jb = min( nb, n-j+1 )
                CALL magmaf_big_dsyrk_cpu( 'Upper', 'Transpose', jb, j-1, -one, a( 1, j ), lda, one, a( j, j ), lda )
                !CALL magmaf_dsyrk_cpu( 'Upper', 'Transpose', jb, j-1, -one, a( 1, j ), lda, one, a( j, j ), lda )
                CALL dpotrf( 'Upper', jb, a( j, j ), lda, info )
                IF( info.NE.0 ) &
                  GO TO 30
                IF( j+jb.LE.n ) THEN
!  
 !                 Compute the current block row.
!  
                   CALL magmaf_BIG_dgemm_cpu( 'Transpose', 'No transpose', jb, n-j-jb+1, &
                              j-1, -one, a( 1, j ), lda, a( 1, j+jb ), 			 &
                              lda, one, a( j, j+jb ), lda )
                   CALL magmaf_dtrsm_CPU( 'Left', 'Upper', 'Transpose', 'Non-unit',		&
                              jb, n-j-jb+1, one, a( j, j ), lda,			&
                              a( j, j+jb ), lda )
                END IF
    10       CONTINUE
!  
          ELSE
             print *, "Compute the Cholesky factorization A = L*L**T."
!  
 !           Compute the Cholesky factorization A = L*L**T.
!  
             DO 20 j = 1, n, nb
!  
 !              Update and factorize the current diagonal block and test
 !              for non-positive-definiteness.
!  
                jb = min( nb, n-j+1 )

                !CALL magmaf_big_dsyrk_cpu( 'Lower', 'No transpose', jb, j-1, -one,a( j, 1 ), lda, one, a( j, j ), lda )
                if(j-1>1)then
                  print *, " big_dsyrk_cpu j,jb,j-1=",j,jb,j-1
		  call magmaf_dsetmatrix( jb, jb, A(j,j), ldA, dT, lddnb, queue )
                  call magmaf_dsetmatrix( jb, j-1, A(j,1), ldA, dT2, lddnb, queue )
		  call magmaf_dsyrk('Lower', 'No transpose',jb, j-1, -one,dT2,lddnb,one,dT,lddnb,queue )
                  call magmaf_dgetmatrix( jb, jb, dT, lddnb, A(j,j),ldA, queue )
                endif

		IF( j+jb.LE.n ) THEN
!  
 !                 Compute the current block column.  big-panel
!                        C := alpha* A * A^T + beta*C,
!  					     !'No transpose','Transpose'
                   print *, " BIG panel BIG_dgemm_cpu m,n,k", n-j-jb+1,jb,j-1
                   CALL magmaf_BIG_dgemm_cpu('N','T',n-j-jb+1,jb,j-1,-one,a(j+jb,1),lda,a(j,1),lda,one,a(j+jb,j),lda) ! funciona 
		 END IF


		print *, " dpotrf j,jb=",j,jb
                CALL dpotrf( 'Lower', jb, a( j, j ), lda, info )
                IF( info.NE.0 )then
                  PRINT *, "INFO=",INFO, " EXIT magmaf_BIG_dpotrf_cpu2"
                  GO TO 30
                ENDIF


                IF( j+jb.LE.n ) THEN
                   print *, " magmaf_dtrsm n-j-jb+1, jb,j+jb =>",n-j-jb+1, jb,j+jb
                   !CALL magmaf_dtrsm_CPU('Right','Lower','Transpose','Non-unit',n-j-jb+1, jb, one, a( j, j ),lda,a( j+jb, j ),lda )
                    call magmaf_dsetmatrix( jb, jb, A(j,j), ldA, dT, lddnb, queue )
                    call magmaf_dsetmatrix( n-j-jb+1, jb, a( j+jb, j ), ldA, dT2, lddn, queue )
		    call magmaf_dtrsm('Right','Lower','Transpose','Non-unit',n-j-jb+1, jb, one,dT,LDDnb,dT2,LDDn,queue)
		    call magmaf_dgetmatrix( n-j-jb+1, jb, dT2, lddn, a( j+jb, j ),ldA, queue )


                END IF
    20       CONTINUE
          END IF
       END IF
       GO TO 40
!  
    30 CONTINUE
       info = info + j - 1
!  
    40 CONTINUE



!! cleanup:
1000 continue
  !! Free GPU memory
        info = magmaf_free( dT )
        if (info .ne. 0) then
            print *, 'TRS Error: magmaf_free( dA ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dT2 )
        if (info .ne. 0) then
            print *, 'TRS Error: magmaf_free( dB ) failed, info = ', info
            stop
        endif
      print *, "***** END magmaf_BIG_dpotrf_cpu2"

    call magmaf_queue_destroy( queue );

       RETURN
!  
 !     End of DPOTRF
!  
       END subroutine magmaf_BIG_dpotrf_cpu2


!  =====================================================================
       SUBROUTINE magmaf_BIG_dpotrf_cpu3( UPLO, N, A, LDA, INFO )

!  
 !     .. Scalar Arguments ..
       CHARACTER          UPLO
       INTEGER            INFO, LDA, N
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION   A( LDA,*   )
 !     ..
!  
!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dT,dT2
       magma_devptr_t          :: queue,queue2
       INTEGER 		       :: LDDnb,LDDn,dev
!  =====================================================================

 !     .. Parameters ..
       DOUBLE PRECISION   ONE, neg_one
       parameter( one = 1.0d0, neg_one=-1.0d0 )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            UPPER
       INTEGER            J, JB, NB
 !     ..
 !     .. External Functions ..
       LOGICAL            LSAME
       INTEGER            ILAENV
       EXTERNAL           lsame, ilaenv
 !     ..
 !     .. External Subroutines ..
       !EXTERNAL           magmaf_BIG_dgemm_cpu, magmaf_dsyrk, magmaf_dtrsm, xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max, min
 !     ..
 !     .. Executable Statements ..
! 

    dev=0
    call magmaf_queue_create( dev, queue )
    call magmaf_queue_create( dev, queue2 )

 
 !     Test the input parameters.
!  

      print *, "***** USING magmaf_BIG_dpotrf_cpu3"

       info = 0
       upper = lsame( uplo, 'U' )
       IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( lda.LT.max( 1, n ) ) THEN
          info = -4
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'DPOTRF', -info )
          RETURN
       END IF
!  
 !     Quick return if possible
!  
       IF( n.EQ.0 ) &
         RETURN
!  
 !     Determine the block size for this environment.
!  
       nb = magmaf_get_dpotrf_nb( n );
       print *, "NB=",nb
       IF( nb.LE.1 .OR. nb.GE.n ) THEN
!  
 !        Use lapack code.
!  
          CALL dpotrf( uplo, n, a, lda, info )
       ELSE

!  
 !        /* Use  blocked code. */
!  
          !! Allocate GPU memory
        lddnb = ceiling(real(nb)/32)*32
        lddn = ceiling(real(n)/32)*32
        info = 0
        
    
          !! Allocate GPU memory
        info = magmaf_dmalloc( dT, lddnb*nb )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dT ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dT2, lddn*lddnb )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dT2  ) failed, info = ", info
            goto 1000
        endif


          IF( upper ) THEN
!  
 !           Compute the Cholesky factorization A = U**T*U.
!  
             DO 10 j = 1, n, nb
!  
 !              Update and factorize the current diagonal block and test
 !              for non-positive-definiteness.
!  
                jb = min( nb, n-j+1 )
                CALL magmaf_big_dsyrk_cpu( 'Upper', 'Transpose', jb, j-1, -one, a( 1, j ), lda, one, a( j, j ), lda )
                !CALL magmaf_dsyrk_cpu( 'Upper', 'Transpose', jb, j-1, -one, a( 1, j ), lda, one, a( j, j ), lda )
                CALL dpotrf( 'Upper', jb, a( j, j ), lda, info )
                IF( info.NE.0 ) &
                  GO TO 30
                IF( j+jb.LE.n ) THEN
!  
 !                 Compute the current block row.
!  
                   CALL magmaf_BIG_dgemm_cpu( 'Transpose', 'No transpose', jb, n-j-jb+1, &
                              j-1, -one, a( 1, j ), lda, a( 1, j+jb ), 			 &
                              lda, one, a( j, j+jb ), lda )
                   CALL magmaf_dtrsm_CPU( 'Left', 'Upper', 'Transpose', 'Non-unit',		&
                              jb, n-j-jb+1, one, a( j, j ), lda,			&
                              a( j, j+jb ), lda )
                END IF
    10       CONTINUE
!  
          ELSE
             print *, "Compute the Cholesky factorization A = L*L**T."
!  
 !           Compute the Cholesky factorization A = L*L**T.
!  
             DO 20 j = 1, n, nb
!  
 !              Update and factorize the current diagonal block and test
 !              for non-positive-definiteness.
!  
                jb = min( nb, n-j+1 )

                !CALL magmaf_big_dsyrk_cpu( 'Lower', 'No transpose', jb, j-1, -one,a( j, 1 ), lda, one, a( j, j ), lda )
                if(j-1>1)then
                  print *, " big_dsyrk_cpu j,jb,j-1=",j,jb,j-1
		  call magmaf_dsetmatrix_async( jb, jb, A(j,j), ldA, dT, lddnb, queue )
                  call magmaf_dsetmatrix_async( jb, j-1, A(j,1), ldA, dT2, lddnb, queue )
		  call magmaf_dsyrk('Lower', 'No transpose',jb, j-1, -one,dT2,lddnb,one,dT,lddnb,queue )
                  call magmaf_queue_sync(queue)
                  call magmaf_dgetmatrix_async( jb, jb, dT, lddnb, A(j,j),ldA, queue2 )
                endif

		IF( j+jb.LE.n ) THEN
                   !  
 !                 Compute the current block column.  big-panel
!                        C := alpha* A * A^T + beta*C,
!  					     !'No transpose','Transpose'
                   print *, " BIG panel BIG_dgemm_cpu m,n,k", n-j-jb+1,jb,j-1
                   CALL magmaf_BIG_dgemm_cpu('N','T',n-j-jb+1,jb,j-1,-one,a(j+jb,1),lda,a(j,1),lda,one,a(j+jb,j),lda) ! funciona 
                   !call magmaf_queue_sync(queue2)

                END IF


                call magmaf_queue_sync(queue2)



		print *, " dpotrf j,jb=",j,jb
                CALL dpotrf( 'Lower', jb, a( j, j ), lda, info )
                IF( info.NE.0 )then
                  PRINT *, "INFO=",INFO, " EXIT magmaf_BIG_dpotrf_cpu2"
                  GO TO 30
                ENDIF
		call magmaf_dsetmatrix_async( jb, jb, A(j,j), ldA, dT, lddnb, queue2 )

                IF( j+jb.LE.n ) THEN


                   print *, " magmaf_dtrsm n-j-jb+1, jb,j+jb =>",n-j-jb+1, jb,j+jb
                   !CALL magmaf_dtrsm_CPU('Right','Lower','Transpose','Non-unit',n-j-jb+1, jb, one, a( j, j ),lda,a( j+jb, j ),lda )
                    !call magmaf_dsetmatrix_async( jb, jb, A(j,j), ldA, dT, lddnb, queue2 )
                    call magmaf_dsetmatrix_async( n-j-jb+1, jb, a( j+jb, j ), ldA, dT2, lddn, queue )
		    call magmaf_dtrsm('Right','Lower','Transpose','Non-unit',n-j-jb+1, jb, one,dT,LDDnb,dT2,LDDn,queue)
		    call magmaf_dgetmatrix_async( n-j-jb+1, jb, dT2, lddn, a( j+jb, j ),ldA, queue )


                END IF
    20       CONTINUE
          END IF
       END IF
       GO TO 40
!  
    30 CONTINUE
       info = info + j - 1
!  
    40 CONTINUE



!! cleanup:
1000 continue
  !! Free GPU memory
        info = magmaf_free( dT )
        if (info .ne. 0) then
            print *, 'TRS Error: magmaf_free( dA ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dT2 )
        if (info .ne. 0) then
            print *, 'TRS Error: magmaf_free( dB ) failed, info = ', info
            stop
        endif
      print *, "* END magmaf_BIG_dpotrf_cpu3"

    call magmaf_queue_destroy( queue );
    call magmaf_queue_destroy( queue2 );

       RETURN
!  
 !     End of DPOTRF
!  
       END subroutine magmaf_BIG_dpotrf_cpu3





!!!!!!!!!!!!!!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=

	!subroutine magmaf_OoC_DPoTRS(UPLO, N, NRHS, A, LDA, B, LDB, INFO)
	subroutine magmaf_OoC_DPoTRS(UPLO, N, NRHS, A, LDA, B, LDB, INFO,L,NB,maxNRHS)


!!  DPOTRS solves a system of linear equations A*X = B with a symmetric
!   positive definite matrix A using the Cholesky factorization
!   A = U**T*U or A = L*L**T computed by DPOTRF.
!
!! 
!!===========================================================


	implicit none
!     .. Scalar Arguments ..
       CHARACTER          UPLO
       INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )

!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dA,dB,dC
       magma_devptr_t          :: queue
       INTEGER 		       :: LDDA,LDDB,LDDC,dev
!  =====================================================================
!     .. Parameters ..
       DOUBLE PRECISION   ONE, neg_one
       parameter( one = 1.0d+0 )
       parameter( neg_one = -1.0d+0 )
!     ..
!     .. Local Scalars ..
       LOGICAL            UPPER
!     ..
!     .. External Functions ..
       LOGICAL            LSAME
       EXTERNAL           lsame
!!  aux
        INTEGER*8 		:: mem_size
        double precision   	:: LimMem,totM
        INTEGER                 :: NB, maxNRHS, nbmax,L,j,i,rank,jb,jr,k, nbmin
	double precision        ::  MaxMem, MemA
!     .. Executable Statements ..

    dev=0
    call magmaf_queue_create( dev, queue )

  
    mem_size=magmaf_mem_size( queue )    !!    mem_size=11 821 973 504  integer
    !!!mem_size=int8(1.0d0*(1000**3))
	!mem_size=221515520_8                 !!  simulando <1GB
        !mem_size=2147483648_8                 !!  simulando 2GB
        !print *, " men size=",mem_size
        ! mem_size = 11821973504_8   !!  full k40c
         !mem_size =  5021973504_8  !! max size 5000*5000
         !mem_size =  1321973504_8  !! max size 5000*5000
	 !mem_size =   303063630_8  !! max size 
         !mem_size =   223696213_8  !! max size 2000*2000
         !mem_size =   120197350_8  !!  max size 1500*1500

    LimMem=dble(mem_size)*0.90d0 !   GB  90% of mex men

    ldda = ceiling(real(n)/32)*32
    !print *, "ldda=",ldda
    MaxMem=( dble(ldda)*dble(n) + dble(ldda)*dble(nrhs) )*8.0d0 
    !print *, "MaxMem=",MaxMem
    MemA=( dble(ldda)*dble(n) + dble(ldda)*dble(512) )*8.0d0 
    !print *, "MemA=",MemA
    
    
    call get_L_Nb(n,L,nbmin)
    
         !! full problem
    IF(maxmem<LimMem)then
        !print *, "full magmaf_OoC_DPoTRS"
        ldda = ceiling(real(n)/32)*32
        lddb = ldda
        lddc = ldda
        info = 0
        maxNRHS=NRHS
        L=1
        NB=N
    
          !! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*n )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dB, lddb*nrhs )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            goto 1000
        endif


        !! dA = hA
        call magmaf_dsetmatrix( n, n, A, lda, dA, ldda, queue )
        !! dB = hB
        call magmaf_dsetmatrix( n, nrhs, B, ldb, dB, lddB, queue )


	call magmaf_dpotrs_gpu( uplo, n, nrhs, dA, ldda, dB, ldda, info )


         !! hB = dB
        call magmaf_dgetmatrix( n, nrhs, dB, lddb, B, ldb, queue )
        


    ELSEIF(MemA<LimMem)then !!  OUT OF CORE, full A in device
        !print *, "full A in device, magmaf_OoC_DPoTRS, blk NRHS"
        ldda = ceiling(real(n)/32)*32
        lddb = ldda
        lddc = ldda
        info = 0
        
        L=1
        NB=N
        
        maxNRHS=int( (LimMem/8.0d0-dble(n)**2) / (2.0d0*dble(n)) )
        !print *, "maxNRHS=",maxNRHS
    
          !! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*n )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dB, lddb*MIN(maxNRHS,nrhs) )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            goto 1000
        endif

        !! dA = A
        call magmaf_dsetmatrix( n, n, A, lda, dA, ldda, queue )
     
	DO  jr=1,nrhs,maxNRHS
            jb = min( maxNRHS, nrhs-jr+1 )
            !print *, "jb=",jb, "jr=",jr
            !! dB = hB
            call magmaf_dsetmatrix( n, Jb, B(1,jr), ldb, dB, lddB, queue )


	    call magmaf_dpotrs_gpu( uplo, n, jb, dA, ldda, db, lddb, info )

            !! hB = B
            call magmaf_dgetmatrix( n, Jb, dB, lddb, B(1,jr), ldb, queue )
	ENDDO


    ELSE  !    OUT OF CORE, BLK
        !print *, "blk magmaf_OoC_DPoTRS"
        upper = lsame( uplo, 'U' )
	!!    get nb and maxNRHS
	!!   1, all nrhs and nb x nb, nb?   need 3 matrix  1 nb x nb, and 2 nb x nrhs
	nbmax=int(sqrt(dble(nrhs)**2+limMem/8.0d0)-dble(nrhs) )
        !print *, "nbmax=",nbmax
	!    test   nb too small
	k=2
	!do while(nbmax<1024)
	do while(nbmax<nbmin)
             nbmax=int(sqrt(dble(nrhs/k)**2+limMem/8.0d0)-dble(nrhs/k) )
             k=k+1
        enddo
        !print *, "nbmax=",nbmax
        !!  nb multiplo 128 < nbmax
        nb=nbmax-127
        !print *, "nb=",nb
        nb = ceiling(real(nb)/128)*128
        !print *, "nb=",nb
        maxNRHS=int( (LimMem/8.0d0-dble(nb)**2) / (2.0d0*dble(nb)) )
        if(mod(n,nb)==0)then
           L=n/nb
        else
	   L=n/nb+1
        endif
        !print *, "nb=",nb,"L=",L,"maxNRHS=",maxNRHS
        
     

        ldda = ceiling(real(nb)/32)*32
        lddb = ldda
        lddc = ldda
        info = 0

        totM=((dble(ldda)*dble(nb)+2.0*dble(lddb)*dble(min(maxNRHS,nrhs)))*8.0)/1e9
        !print *, "n=",n," ldda = ", ldda, " nrhs", nrhs," Memory in device=", totM," GB"

          !! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*nb )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dB, lddb*MIN(maxNRHS,nrhs) )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dC, lddc*MIN(maxNRHS,nrhs) )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dC  ) failed, info = ", info
            goto 1000
        endif




 	IF( upper ) THEN		
!           Solve A*X = B where A = U**T *U.
!
!           Solve U**T *X = B, overwriting B with X.

!           Solve U*X = B, overwriting B with X.

        ELSE 
!           Solve A*X = B where A = L*L**T.
!
!            print *, "Solve L*X = B, overwriting B with X."
!



            DO jr=1,nrhs,maxNRHS
               jb = min( maxNRHS, nrhs-jr+1 ) 
               			!print *, "jr loop    jb=",jb, "jr=",jr
               DO j=1,L
                   rank=nb
   		   if(j==L)rank=n-(L-1)*nb
                   		!print *, "j  loop    j=",j,"rank=",rank
		   ! set Bj
                   		!print *, "set Bj=  B(",(j-1)*nb+1,",",jr,")"
                   !! dB = hB
                   call magmaf_dsetmatrix( rank, Jb, B((j-1)*nb+1,jr), ldb, dB, lddB, queue )
                   
                   !!  Bj-Lj(i=1,j-1)Y(i=1,j-1)
    		   DO i=j-1,1,-1
                         		!print *, ""
                         		!print *, "i  loop    j=",j, "i=",i,"rank=",rank
		         ! set Lji
                         		!print *, "set Lji=  A(",(j-1)*nb+1,",",(i-1)*nb+1,")"
                         call magmaf_dsetmatrix( rank, nb, A((j-1)*nb+1,(i-1)*nb+1), ldb, dA, lddA, queue )
                         		!print *, "set Yi=  B(",(i-1)*nb+1,",",jr,")"
                         ! set Yi
                         if(i < j-1)then
                             call magmaf_dsetmatrix( nb, jb, B((i-1)*nb+1,jr), ldb, dC, lddC, queue )
                         endif
                         ! Bj=Bj-LjiYi  
                         call magmaf_dgemm(  'N', 'N', rank, jb,nb, neg_one, dA, ldda, dC, lddC,one,dB, lddB, queue)

                    ENDDO
                    !! set Ljj
                    		!print *, "set Ljj      j=",j,"rank=",rank
                    		!print *, "set Ljj=  A(",(j-1)*nb+1,",",(j-1)*nb+1,")"
                    call magmaf_dsetmatrix( rank, rank, A((j-1)*nb+1,(j-1)*nb+1), ldA, dA, lddA, queue )
                   		!magmablasf_dtrsm( side, uplo, transA, diag, m, n, alpha, dA, ldda, dB, lddb,  queue )
		   		!magmaf_dtrsm    ( side, uplo, trans, diag, m, n, alpha, dA, ldda, dB, lddb, queue  )
    		    call magmaf_dtrsm ( 'Left', 'Lower', 'N', 'N', rank, jb, one, dA, ldda, dB, lddb, queue)
		    if(j<L)call magmablasf_dlacpy( "F", rank, jb, dB, lddb, dC, lddb, queue )

		    !!  get  Yj 
                    		!print *, "get Bj=  B(",(j-1)*nb+1,",",jr,")"
                    call magmaf_dgetmatrix( rank, Jb, dB, lddB, B((j-1)*nb+1,jr), ldB, queue )
			

	       ENDDO  !  J=1,L
  

!
!       	    print *, "Solve L**T *X = B, overwriting B with X."
!

	       DO j=L,1,-1
                   rank=nb
   		   if(j==L)rank=n-(L-1)*nb
                   		!print *, "j  loop    j=",j,"rank=",rank
                   		!print *, "jr loop    jb=",jb, "jr=",jr
		   ! set Bj
                   		!print *, "set Bj=  B(",(j-1)*nb+1,",",jr,")"
                   !! dB = hB
                   if(j<L)call magmaf_dsetmatrix( rank, Jb, B((j-1)*nb+1,jr), ldb, dB, lddB, queue )

                   !!  Bj-L(i=j+1,L)j^T Xi
    		   DO i=j+1,L
                        rank=nb
                        if(i==L)rank=n-(L-1)*nb
                        		!print *, "i  loop    j=",j, "i=",i,"rank=",rank
			! set Lij^T
                        		!print *, "set (Lij)^T=  A(",(i-1)*nb+1,",",(j-1)*nb+1,")"
                        call magmaf_dsetmatrix( rank, nb, A((i-1)*nb+1,(j-1)*nb+1), ldb, dA, lddA, queue )
                        ! set Xi
					!print *, "set Xi=  B(",(i-1)*nb+1,",",jr,")"
                        if(i > j+1)then
                           call magmaf_dsetmatrix( rank, jb, B((i-1)*nb+1,jr), ldb, dC, lddC, queue )
                        endif

                        !!  Bj=Bj-Lij^T Xi    
                        call magmaf_dgemm(  'T', 'N', nb, jb,rank, neg_one, dA, ldda, dC, lddC,one,dB, lddB, queue)

                   ENDDO  !!i=L,j+1,-1
                   rank=nb
                   if(j==L)rank=n-(L-1)*nb
                   		! print *, "set Ljj^T      j=",j,"rank=",rank
                    !! set Ljj
                   		!print *, "set Ljj=  A(",(j-1)*nb+1,",",(j-1)*nb+1,")"
                   call magmaf_dsetmatrix( rank, rank, A((j-1)*nb+1,(j-1)*nb+1), ldA, dA, lddA, queue )
                   		!magmablasf_dtrsm( side, uplo, transA, diag, m, n, alpha, dA, ldda, dB, lddb,  queue )
		   		!magmaf_dtrsm    ( side, uplo, transA, diag, m, n, alpha, dA, ldda, dB, lddb, queue  )
    		   call magmaf_dtrsm ( 'Left', 'Lower', 'T', 'N', rank, jb, one, dA, ldda, dB, lddb, queue)
		   if(j>1)call magmablasf_dlacpy( "F", rank, jb, dB, lddb, dC, lddc, queue )
		!!  get  Xj 
                    		!print *, "get Bj=  B(",(j-1)*nb+1,",",jr,")"
                    call magmaf_dgetmatrix( rank, Jb, dB, lddb, B((j-1)*nb+1,jr), ldb, queue )
                ENDDO  ! j=L,1,-1



	    ENDDO  !! JR



       ENDIF  ! ! uper
        info = magmaf_free( dC )
        if (info .ne. 0) then
            print *, 'TRS Error: magmaf_free( dC ) failed, info = ', info
            stop
        endif


    ENDIF  !! size: full, full A, blk(A)
	

!! cleanup:
1000 continue
  !! Free GPU memory
        info = magmaf_free( dA )
        if (info .ne. 0) then
            print *, 'TRS Error: magmaf_free( dA ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dB )
        if (info .ne. 0) then
            print *, 'TRS Error: magmaf_free( dB ) failed, info = ', info
            stop
        endif


    call magmaf_queue_destroy( queue );


409   format((9ES13.4))
509   format((9(a13)))
!finish=omp_get_wtime() 
!write(*,509) "t copy =","t dM","taM","tpo","total"
!write(*,409) tcpy,tdM,taM,tpo,finish-start

	return

	end subroutine magmaf_OoC_DPoTRS




   subroutine get_L_Nb(n,L,nc)
        !!  for use GPU out of core
        !! igual a get l nc
	implicit none
        !integer :: n,L,nc
        integer, intent(in) :: n
        integer, intent(out) :: L,nc
        integer         :: nc2
!  =====================================================================
!     .. Device Array Arguments ..GPU
        magma_devptr_t          :: queue
        INTEGER                 :: dev
!!  aux
        INTEGER*8 		:: mem_size
        double precision   	:: LimMem	!,totM
	double precision        ::  MemA, MenL2	!,MaxMem

    	dev=0
    	call magmaf_queue_create( dev, queue )
  
    	mem_size=magmaf_mem_size( queue )    !!    mem_size=11 821 973 504  integer
         !mem_size=221515520_8                 !!  simulando <1GB
         !mem_size=2147483648_8                 !!  simulando 2GB
        !print *, " men size=",mem_size
        ! mem_size = 11821973504_8   !!  full k40c
         !mem_size =  5021973504_8  !! max size 5000*5000
         !mem_size =  1321973504_8  !! max size 5000*5000
	 !mem_size =   303063630_8  !! max size   2368*2368
         !mem_size =   221973504_8  !! max size 2000*2000
         !mem_size =    120197350_8  !! max size 1500*1500 L>3, (L=1, 3500)

    	LimMem=dble(mem_size)*0.90d0 !   GB  90% of max men
        !print *, "LimMem=",LimMem
    	memA=( dble(n)*dble(n) )*8.0d0 
        !print *, "memA=",memA
        menL2=(( 3.0d0*dble(n/2)*dble(n/2) )*8.0d0 )
         !print *, "menL2=",menL2

    	if(memA<LimMem)then
            !print *, "L=1"
            L=1
            nc=n
        elseif(menL2<LimMem*0.7d0)then   !!   3 matrix nc*nc, 30% for nrhs
            !print *, "L=2"
	        L=2
            IF ( mod (n,L) .ne. 0 ) then
               Nc=n/L+1
            else
                Nc=n/L
            endif
            !print *, "(( 3.0d0*dble(n/2)*dble(n/2) )*8.0d0 =",(( 3.0d0*dble(nc)*dble(nc) )*8.0d0)

        else   !   6 matrix nc*nc
             !print *, "L>2"
             nc=int(sqrt(LimMem/(6.0d0*8.0d0)))
             IF ( mod (n,nc) .ne. 0 ) then
               L=n/nc+1
             else
               L=n/nc
             endif
             nc2=n/L
             !print *, "nc=",nc,"nc2=",nc2
             nc=min(nc,nc2)
             do while(nc*L<n)
                 nc=nc+1
             enddo
         endif

         NC = ceiling(real(nc)/32)*32

! ajust, new L (if necesary)
    do while((L-1)*Nc>=N)
       !print*," L=L-1"
       L=L-1
    enddo

    !Verfic
    if(L*Nc<N .or. (L-1)*Nc>=N)then
       print*," erro cal Nc"
    endif


         call magmaf_queue_destroy( queue )
         return
    end subroutine get_L_Nb




!!!!!!!!!!!!!!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=

	subroutine magmaf_OoC_DPoTRS2(UPLO, N, NRHS, A, LDA, B, LDB, INFO)   !!   


!!  DPOTRS solves a system of linear equations A*X = B with a symmetric
!   positive definite matrix A using the Cholesky factorization
!   A = U**T*U or A = L*L**T computed by DPOTRF.
!
!! 
!!===========================================================


	implicit none
!     .. Scalar Arguments ..
       CHARACTER          UPLO
       INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )

!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dA,dB,dC
       magma_devptr_t          :: queue
       INTEGER 		       :: LDDA,LDDB,LDDC,dev
!  =====================================================================
!     .. Parameters ..
       DOUBLE PRECISION   ONE, neg_one
       parameter( one = 1.0d+0 )
       parameter( neg_one = -1.0d+0 )
!     ..
!     .. Local Scalars ..
       LOGICAL            UPPER
!     ..
!     .. External Functions ..
       LOGICAL            LSAME
       EXTERNAL           lsame
!!  aux
        INTEGER*8 		:: mem_size
        double precision   	:: LimMem,totM
        INTEGER                 :: NB, maxNRHS, nbmax,L,j,i,rank,jb,jr,k
	double precision        ::  MaxMem, MemA
!     .. Executable Statements ..

    dev=0
    call magmaf_queue_create( dev, queue )
  
    mem_size=magmaf_mem_size( queue )    !!    mem_size=11 821 973 504  integer
    !mem_size=int8(0.450d0*(1000**3)
!mem_size=221515520_8                 !!  simulando <1GB
        !print *, " men size=",mem_size
        ! mem_size = 11821973504_8   !!  full k40c
         !mem_size =  5021973504_8  !! max size 5000*5000
         !mem_size =  1321973504_8  !! max size 5000*5000
         !mem_size =   221973504_8  !! max size 2000*2000
         ! mem_size =   120197350_8  !!  max size 1500*1500


    LimMem=dble(mem_size)*1.00d0 !   GB  90% of mex men

    ldda = ceiling(real(n)/32)*32
    !print *, "ldda=",ldda
    MaxMem=( dble(ldda)*dble(n) + dble(ldda)*dble(nrhs) )*8.0d0 
    !print *, "MaxMem=",MaxMem
    MemA=( dble(ldda)*dble(n) + dble(ldda)*dble(512) )*8.0d0 
    !print *, "MemA=",MemA
         !! full problem
    IF(maxmem<LimMem)then
        !print *, "full magmaf_OoC_DPoTRS"
        ldda = ceiling(real(n)/32)*32
        lddb = ldda
        lddc = ldda
        info = 0
        
    
          !! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*n )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dB, lddb*nrhs )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            goto 1000
        endif


        !! dA = hA
        call magmaf_dsetmatrix( n, n, A, lda, dA, ldda, queue )
        !! dB = hB
        call magmaf_dsetmatrix( n, nrhs, B, ldb, dB, lddB, queue )


	call magmaf_dpotrs_gpu( uplo, n, nrhs, dA, ldda, dB, ldda, info )


         !! hB = dB
        call magmaf_dgetmatrix( n, nrhs, dB, lddb, B, ldb, queue )
        


    ELSEIF(MemA<LimMem)then !!  OUT OF CORE, full A in device
        !print *, "full A in device, magmaf_OoC_DPoTRS"
        ldda = ceiling(real(n)/32)*32
        lddb = ldda
        lddc = ldda
        info = 0
        
        maxNRHS=int( (LimMem/8.0d0-dble(n)**2) / (2.0d0*dble(n)) )
        !print *, "maxNRHS=",maxNRHS
    
          !! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*n )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dB, lddb*MIN(maxNRHS,nrhs) )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            goto 1000
        endif

        !! dA = A
        call magmaf_dsetmatrix( n, n, A, lda, dA, ldda, queue )
     
	DO  jr=1,nrhs,maxNRHS
            jb = min( maxNRHS, nrhs-jr+1 )
            print *, "jb=",jb, "jr=",jr
            !! dB = hB
            call magmaf_dsetmatrix( n, Jb, B(1,jr), ldb, dB, lddB, queue )


	    call magmaf_dpotrs_gpu( uplo, n, jb, dA, ldda, db, lddb, info )

            !! hB = B
            call magmaf_dgetmatrix( n, Jb, dB, lddb, B(1,jr), ldb, queue )
	ENDDO


    ELSE  !    OUT OF CORE, BLK
        !print *, "blk magmaf_OoC_DPoTRS"
        upper = lsame( uplo, 'U' )
	!!    get nb and maxNRHS
	!!   1, all nrhs and nb x nb, nb?   need 3 matrix  1 nb x nb, and 2 nb x nrhs
	nbmax=int(sqrt(dble(nrhs)**2+limMem/8.0d0)-dble(nrhs) )
        !print *, "nbmax=",nbmax
	!    test   nb too small
	k=2
	do while(nbmax<512)
             nbmax=int(sqrt(dble(nrhs/k)**2+limMem/8.0d0)-dble(nrhs/k) )
             k=k+1
        enddo
            !print *, "nbmax=",nbmax
        !!  nb multiplo 128 < nbmax
        nb=nbmax-127
            !print *, "nb=",nb
        nb = ceiling(real(nb)/128)*128
        !nb=4000
            !print *, "nb=",nb
        maxNRHS=int( (LimMem/8.0d0-dble(nb)**2) / (2.0d0*dble(nb)) )
            !print *, "maxNRHS=",maxNRHS
        if(mod(n,nb)==0)then
           L=n/nb
        else
	   L=n/nb+1
        endif
            !print *, "nb=",nb,"L=",L
        
     

        ldda = ceiling(real(nb)/32)*32
        lddb = ldda
        lddc = ldda
        info = 0

        totM=((dble(ldda)*dble(nb)+2.0*dble(lddb)*dble(min(maxNRHS,nrhs)))*8.0)/1e9
            !print *, "n=",n," ldda = ", ldda, " nrhs", nrhs," Memory in device=", totM," GB"

          !! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*nb )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dB, lddb*MIN(maxNRHS,nrhs) )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dC, lddc*MIN(maxNRHS,nrhs) )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dC  ) failed, info = ", info
            goto 1000
        endif


 	IF( upper ) THEN		
!           Solve A*X = B where A = U**T *U.
!
!           Solve U**T *X = B, overwriting B with X.

!           Solve U*X = B, overwriting B with X.

        ELSE 

!           Solve A*X = B where A = L*L**T.
!
            !print *, "Solve L*X = B, overwriting B with X."
!
            DO j=1,L
                !print *, ""
                !print *, ""
                rank=nb
   		if(j==L)rank=n-(L-1)*nb
                    !print *, "j  loop    j=",j,"rank=",rank
                DO jr=1,nrhs,maxNRHS
           	   jb = min( maxNRHS, nrhs-jr+1 ) 
                   !print *, "jr loop    jb=",jb, "jr=",jr
		   ! set Bj
                        !print *, "set Bj=  B(",(j-1)*nb+1,",",jr,")"
                   !! dB = hB
                   call magmaf_dsetmatrix( rank, Jb, B((j-1)*nb+1,jr), ldb, dB, lddB, queue )
                   
                   !!  Bj-Lj(i=1,j-1)Y(i=1,j-1)
    		   DO i=j-1,1,-1
                        !print *, ""
                        !print *, "i  loop    j=",j, "i=",i,"rank=",rank
			! set Lji
                            !print *, "set Lji=  A(",(j-1)*nb+1,",",(i-1)*nb+1,")"
                        call magmaf_dsetmatrix( rank, nb, A((j-1)*nb+1,(i-1)*nb+1), ldb, dA, lddA, queue )
                            !print *, "set Yi=  B(",(i-1)*nb+1,",",jr,")"
                        ! set Yi
                        call magmaf_dsetmatrix( nb, jb, B((i-1)*nb+1,jr), ldb, dC, lddC, queue )

                        ! Bj=Bj-LjiYi  
                        call magmaf_dgemm(  'N', 'N', rank, jb,nb, neg_one, dA, ldda, dC, lddC,one,dB, lddB, queue)

                   ENDDO
                    !! set Ljj
                        !print *, "set Ljj      j=",j,"rank=",rank
                        !print *, "set Ljj=  A(",(j-1)*nb+1,",",(j-1)*nb+1,")"
                   call magmaf_dsetmatrix( rank, rank, A((j-1)*nb+1,(j-1)*nb+1), ldA, dA, lddA, queue )
                   		!magmablasf_dtrsm( side, uplo, transA, diag, m, n, alpha, dA, ldda, dB, lddb,  queue )
		   		!magmaf_dtrsm    ( side, uplo, trans, diag, m, n, alpha, dA, ldda, dB, lddb, queue  )
    		   call magmaf_dtrsm ( 'Left', 'Lower', 'N', 'N', rank, jb, one, dA, ldda, dB, lddb, queue)

		!!  get  Yj 
                        !print *, "get Bj=  B(",(j-1)*nb+1,",",jr,")"
                    call magmaf_dgetmatrix( rank, Jb, dB, lddB, B((j-1)*nb+1,jr), ldB, queue )
			

		ENDDO  !  JR
   	    ENDDO  !  J=1,L
	print *, ""
	print *, ""
	print *, ""

!
       	    !print *, "Solve L**T *X = B, overwriting B with X."
!

	    DO j=L,1,-1
                !print *, ""
                !print *, ""
                rank=nb
   		if(j==L)rank=n-(L-1)*nb
                !print *, "j  loop    j=",j,"rank=",rank
		DO jr=1,nrhs,maxNRHS
           	   jb = min( maxNRHS, nrhs-jr+1 ) 
                    !print *, "jr loop    jb=",jb, "jr=",jr
		   ! set Bj
                    ! print *, "set Bj=  B(",(j-1)*nb+1,",",jr,")"
                   !! dB = hB
                   call magmaf_dsetmatrix( rank, Jb, B((j-1)*nb+1,jr), ldb, dB, lddB, queue )

                   !!  Bj-L(i=j+1,L)j^T Xi
    		   DO i=j+1,L
                        rank=nb
                        if(i==L)rank=n-(L-1)*nb
                            !print *, "i  loop    j=",j, "i=",i,"rank=",rank
			! set Lij^T
                            !print *, "set (Lij)^T=  A(",(i-1)*nb+1,",",(j-1)*nb+1,")"
                        call magmaf_dsetmatrix( rank, nb, A((i-1)*nb+1,(j-1)*nb+1), ldb, dA, lddA, queue )
                        ! set Xi
			                !print *, "set Xi=  B(",(i-1)*nb+1,",",jr,")"
                        call magmaf_dsetmatrix( rank, jb, B((i-1)*nb+1,jr), ldb, dC, lddC, queue )

                        !!  Bj=Bj-Lij^T Xi    
                        call magmaf_dgemm(  'T', 'N', nb, jb,rank, neg_one, dA, ldda, dC, lddC,one,dB, lddB, queue)

                   ENDDO  !!i=L,j+1,-1
                   rank=nb
                   if(j==L)rank=n-(L-1)*nb
                        !print *, "set Ljj^T      j=",j,"rank=",rank
                    !! set Ljj
                        !print *, "set Ljj=  A(",(j-1)*nb+1,",",(j-1)*nb+1,")"
                   call magmaf_dsetmatrix( rank, rank, A((j-1)*nb+1,(j-1)*nb+1), ldA, dA, lddA, queue )
                   		!magmablasf_dtrsm( side, uplo, transA, diag, m, n, alpha, dA, ldda, dB, lddb,  queue )
		   		!magmaf_dtrsm    ( side, uplo, transA, diag, m, n, alpha, dA, ldda, dB, lddb, queue  )
    		   call magmaf_dtrsm ( 'Left', 'Lower', 'T', 'N', rank, jb, one, dA, ldda, dB, lddb, queue)

		!!  get  Xj 
                        !print *, "get Bj=  B(",(j-1)*nb+1,",",jr,")"
                    call magmaf_dgetmatrix( rank, Jb, dB, lddb, B((j-1)*nb+1,jr), ldb, queue )
                ENDDO  ! JR
	    ENDDO  !!j=L,1,-1



       ENDIF  ! ! uper
        info = magmaf_free( dC )
        if (info .ne. 0) then
            print *, 'TRS Error: magmaf_free( dC ) failed, info = ', info
            stop
        endif

    ENDIF  !! size: full, full A, blk(A)
	

!! cleanup:
1000 continue
  !! Free GPU memory
        info = magmaf_free( dA )
        if (info .ne. 0) then
            print *, 'TRS Error: magmaf_free( dA ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dB )
        if (info .ne. 0) then
            print *, 'TRS Error: magmaf_free( dB ) failed, info = ', info
            stop
        endif


    call magmaf_queue_destroy( queue );


409   format((9ES13.4))
509   format((9(a13)))
!finish=omp_get_wtime() 
!write(*,509) "t copy =","t dM","taM","tpo","total"
!write(*,409) tcpy,tdM,taM,tpo,finish-start

	return

	end subroutine magmaf_OoC_DPoTRS2





subroutine magmaf_get_dgemm_mb_nb(limmem,m,n,k,mb,nb)

implicit none
!  =====================================================================
!     .. Scalar Arguments ..
   !CHARACTER        :: CRwise
   DOUBLE PRECISION :: limMem
   integer          :: m,n,k
   integer          :: mb,nb, nbmax !,L1,L2
   logical          :: okay
   DOUBLE PRECISION :: Mem,tem

okay=.false.


do while(okay .eqv. .false.)


!! c(mb,nb)=c(mb,nb)+A(mb,k)*B(k,nb)  in GPU  mb=nb
    nbmax=int(sqrt(dble(k)**2+limMem/8.0d0)-dble(k) )
    nb=nbmax; mb=nb
    !print *, "nbmax=",nbmax
    if(nbmax>m/2 .and. nbmax<m)then
                mb=m/2
    elseif(nbmax>m)then
		mb=m
    endif
    if(nbmax>n/2 .and. nbmax<n)then
           nb=n/2
    elseif(nbmax>n)then
	nb=n
    endif
    call get_mem_dgemm(mb,nb,k,Mem)
    if(mem<limmem)then
         !print *, "c(mb,nb)=c(mb,nb)+A(mb,k)*B(k,nb)  in GPU  mb=nb"
         !print *, "mb=",mb,"nb=",nb
          !print *, " men dgemm in dev",Mem
         okay=.true.
         goto 1000
    endif

1000 continue

end do

    nb = ceiling(real(nb)/128)*128
    mb = ceiling(real(mb)/128)*128

     !print *, "ceiling mb=",mb,"nb=",nb
    call get_mem_dgemm(mb,nb,k,tem)
      !print *, " men dgemm in dev",tem

end subroutine magmaf_get_dgemm_mb_nb


subroutine get_mem_dgemm(m,n,k,Mem)
implicit none
!  =====================================================================
!     .. Scalar Arguments ..
DOUBLE PRECISION :: mem
integer          :: m,n,k

Mem = ( dble(m)*dble(n) + dble(m)*dble(k) +dble(k)*dble(n) )*8.0d0 

end subroutine get_mem_dgemm


	subroutine matrixPartitioner(M,  N,  K, mem_size, h_groups, v_groups,height, width) 
!!  https://csgitlab.ucd.ie/manumachu/libhclooc/blob/master/src/hclgemmOOC.cpp
	implicit none
!  =====================================================================
!     .. Scalar Arguments ..
	INTEGER :: M,  N,  K, h_groups, v_groups, height, width
	INTEGER*8 :: mem_size
!     .. Scalar work ..
        INTEGER   :: temp_h, temp_v,criterion
        INTEGER*8 :: t1, te, t2, t3, t4
!  =====================================================================
        criterion=0
	h_groups=2
        mem_size=mem_size/8

        t1=M*N
        te=N*K
        t2=te*h_groups
        t3 = mem_size * h_groups
        t4 = M * K

	v_groups = ceiling( (2 * t1 + t2) / ( real(t3 - 2 * t4)) )
	criterion = h_groups * v_groups

        temp_h = 3
        temp_v=huge(temp_v)

        !print*,"h_groups",h_groups, "v_groups", v_groups, "criterion", criterion
        !print*, "te", te,"t2", t2,"t3", t3
        
	do while(temp_h <= M .and. (temp_v <= 0 .or. temp_v > 1) )
	    t2 = te * temp_h
            t3 = mem_size * temp_h
            !print *, "te =", te," t2 = ",t2," t3 = ",t3," h = ",temp_h," v = ",temp_v," criterion", criterion

            temp_v = ceiling((2 * t1 + t2) / (real (t3 - 2 * t4)))
	    if (temp_v > 0)then
                 !print *,"temp_h",temp_h," temp_v", temp_v
                 if ( v_groups <= 0 .or. (temp_v < v_groups .and. criterion >= temp_h * temp_v))then
                    ! keeping temp_h * temp_v minimum prevent matrix A from excessive partitioning
                    ! when matrix A is split into too many small slices, communication costs of in-core
                    ! computation outweighs computation costs and it degrades total performance.

                    h_groups = temp_h
                    v_groups = temp_v
                    criterion = temp_h * temp_v
                 endif
            endif
             temp_h=temp_h+1
        enddo


        if ( v_groups <= 0 .or. v_groups > N)then
            print *, "ERROR: No workload distribution to fit into the accelerator memory." 
            stop
        endif
         !print *, "h_groups= ",h_groups ," v_groups=",v_groups


         !height = (int *) malloc(sizeof(int) * (*h_groups));
         !width = (int *) malloc(sizeof(int) * (*v_groups));
           height=M/h_groups-127
           width=N/v_groups-127

            height = ceiling(real(height)/128)*128
             width = ceiling(real(width)/128)*128
        !print *, "height= ",height ," width=",width

        !*height = (int *) malloc(sizeof(int) * (*h_groups));
        !*width = (int *) malloc(sizeof(int) * (*v_groups));

        !int i;
        !for (i = 0; i < *h_groups; i++) {
        !    (*height)[i] = M / (*h_groups);
        !}

        !for (i = 0; i < M % (*h_groups); i++) {
        !    (*height)[i]++;
        !}

        !for (i = 0; i < *v_groups; i++) {
        !    (*width)[i] = N / (*v_groups);
        !}

        !for (i = 0; i < N % (*v_groups); i++) {
        !    (*width)[i]++;
        !}


	end subroutine matrixPartitioner

subroutine magmaf_BIG_dsyrk_cpu(UPLO,TRANS,N,K,ALPHA,hA,LDA,BETA,hC,LDC)
!DSYRK  performs one of the symmetric rank k operations
!
!    C := alpha*A*A**T + beta*C,    TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.
!
! or
!
!    C := alpha*A**T*A + beta*C,    TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.
!
! where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
! and  A  is an  n by k  matrix in the first case and a  k by n  matrix
! in the second case.

!     .. Scalar Arguments ..
       DOUBLE PRECISION ALPHA,BETA
       INTEGER LDA,LDC,N,K
       CHARACTER TRANS, UPLO
!     .. Array Arguments ..CPU
       DOUBLE PRECISION hA(LDA,*),hC(LDC,*)
!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dA, dB, dC
       magma_devptr_t          :: queue
       INTEGER LDDA,LDDB,LDDC,dev
!  =====================================================================
!     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame
!     .. Local Scalars ..
       LOGICAL            UPPER
       LOGICAL NOTA
       INTEGER INFO
       INTEGER*8 :: mem_size
       !INTEGER :: h_groups, v_groups, height, width
       double precision                    ::  zero,one
       parameter                       ( zero=0.0d0,one=1.0d0)
       double precision                    :: LimMem,totM, MaxMem
       integer                             :: j,i,jb,ib,NCOLA,NROWA !,NROWB
       integer                             :: nb	!,mb ,maxnb
!  =====================================================================

			!print *, "start  magmaf_BIG_dsyrk_cpu n=",n,"k=",k
    dev=0
    call magmaf_queue_create( dev, queue )

    nota  = lsame(trans,'N')
    upper = lsame( uplo, 'U' )

!     Quick return if possible.
!
       IF ((n.EQ.0) .OR. (((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one)))then
          !print *, "Quick return"
          RETURN
       endif

    mem_size=magmaf_mem_size( queue )    !!    mem_size=11 821 973 504  integer
    !mem_size=int8(0.10d0*(1024**3))
    !mem_size=2500
    !print *, "mem_size=",mem_size,"GB=",mem_size/1e9
    LimMem=dble(mem_size)*0.90d0 !   GB  90% of mex men
    MaxMem=( dble(n)*dble(n) + dble(n)*dble(k)  )*8.0d0 
    !print *, "memory full problem= ",MaxMem, "limit memory device=",LimMem

    
     
    !call matrixPartitioner(M,  N,  K, mem_size, h_groups, v_groups,height, width)

    ! print *, "h_groups=",h_groups, "v_groups=",v_groups,"height=",height, "width=",width

     !print *,"Partitioner Mem in dev=", ((dble(width)*dble(height)+dble(height)*dble(k)+dble(width)*dble(k))*8.0)/1e9 ," GB"

if(maxmem<LimMem)then
   !! full problem
      !print *, "Full problem  magmaf_BIG_dsyrk_cpu"
    lddc = ceiling(real(n)/32)*32
    ldda = ceiling(real(n)/32)*32
    info = 0



      !print *, " Allocate GPU memory"
        info = magmaf_dmalloc( dA, ldda*k )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif


        info = magmaf_dmalloc( dC, lddc*n )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dC  ) failed, info = ", info
            goto 1000
        endif

        !print *, " dset dA"
    	!! dA = hA
        call magmaf_dsetmatrix( n, k, hA, lda, dA, ldda, queue )
        !call magmaf_dsetmatrix_async( n, k, hA, lda, dA, ldda, queue )
	!print *, " finish dset dA"
        if(BETA/=zero)then
                !print *, " dset dC"
		!! dC = hC
        	!call magmaf_dsetmatrix( n, n, hC, ldC, dC, lddC, queue )
                call magmaf_dsetmatrix_async( n, n, hC, ldC, dC, lddC, queue )
         endif
        !call magmaf_queue_sync(queue)
        !print *, " magmaf_dsyrk"
        call magmaf_dsyrk(UPLO,TRANS,n,k,ALPHA,dA,ldda,BETA,dC,lddc,queue )

 	!! hC = dC
        call magmaf_dgetmatrix( n, n, dC, lddC, hC, ldC, queue )
        !call magmaf_dgetmatrix_async( n, n, dC, lddC, hC, ldC, queue )

ELSE   !!  	by block

    if( (dble(n)*dble(512) + dble(n)*dble(k) +dble(k)*dble(512)) *8.0d0 <LimMem/1000)then   !!  by colums block
        !!   full A in device and nb=>512,  =>> C(n,nb)=beta*C(n,nb)+alpha*A*A**T 

       IF (notA) THEN    
	     !!  C := alpha*A*A**T + beta*C.
          
       ELSE
          !!  C := alpha*A*B**T + beta*C.  or C := alpha*A**T*B**T + beta*C.
             PRINT *, " NOT IMPLEMENT YET, full A in device"
             PRINT *, " C := alpha*A*B**T + beta*C.  or C := alpha*A**T*B**T + beta*C."

       ENDIF  
    
    ELSE   !!  BY BLK(mb,nb)
        !!   cal C by C(nb,nb)=beta*C(nb,nb)+alpha*A(nb,k)**T*A(k,nb)  or beta*C(nb,nb)+alpha*A(nb,n)*A(**T
       IF (nota) THEN
           nrowa = n
           ncola = k
       ELSE
           nrowa = k
           ncola = n
       END IF


       call magmaf_get_dgemm_mb_nb(limMem,n,n,k,nb,nb)
       !call matrixPartitioner(M,  N,  K, int8(limMem), h_groups, v_groups,mb, nb)
       !print *, "magmaf_BIG_dsyrk_cpu; nb", nb," nb",nb

      lddc = ceiling(real(nb)/32)*32
      IF (nota) THEN
          ldda = ceiling(real(nb)/32)*32
          lddb = ceiling(real(k)/32)*32
       ELSE
          ldda = ceiling(real(k)/32)*32
          lddb = ceiling(real(nb)/32)*32
       END IF

      

 
            !! Allocate GPU memory
      info = magmaf_dmalloc( dC, lddc*nb )
      if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dC  ) failed, info = ", info
            goto 1000
      endif



      IF (upper) THEN
           !print *, "cal UP "

           IF (nota) THEN
		          ! C(nb,nb) := alpha*A*A**T + beta*C(nb,nb). 
                 print *, "not tested yet "
                ldda = ceiling(real(nb)/32)*32
       		lddb = ceiling(real(nb)/32)*32
       		info = 0
       		print *," Mem in dev=", ((dble(nb)*dble(nb)+dble(nb)*dble(k)+dble(nb)*dble(k))*8.0)/1e9 ," GB"
       		totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(nb)+dble(lddb)*dble(nb))*8.0)/1e9
       		print *, "magmaf_BIG_dsyrk_cpu; Memory in device=", totM," GB"


      		!! Allocate GPU memory
        	info = magmaf_dmalloc( dA, ldda*k )
        	if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            		goto 1000
        	endif

        	info = magmaf_dmalloc( dB, lddb*k )
        	if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            		goto 1000
        	endif
        	print *, "end Allocate GPU memory"
		!!  num blk


        	DO j=1,n,nb
                	jb = min( nb, n-j+1 )
			!! dA = hA(j,1)
        		!call magmaf_dsetmatrix( k, jb, hA(1,j), ldA, dA, lddA, queue )
                	call magmaf_dsetmatrix_async( jb,k, hA(j,1), ldA, dA, lddA, queue )
           		DO i=1,j,nb
                            ib = min( nb, n-i+1 )
                             		print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib
                            if(BETA/=zero)then
				!! dC = hC
        			call magmaf_dsetmatrix( ib, jb, hC(i,j), ldC, dC, lddC, queue )
                  	    endif
                            if(j==i)then
                               print *, "magmaf_dsyrk"
                               call magmaf_dsyrk( "U", "N", nb, k, alpha, dA, lddA, beta, dC, lddc, queue )

                            else
                                print *, "set B"
				!! dB = hA(1,i)
        			!call magmaf_dsetmatrix( k, jb, hA(1,i), ldA, dA, lddA, queue )
                		call magmaf_dsetmatrix_async( ib,k, hA(i,1), ldA, dB, lddB, queue )
                				
                                print *, "magmaf_dgemm"
                		call magmaf_queue_sync(queue) 	
				!call magmablasf_dgemm(TRANSA,TRANSB,m,n,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )
				call magmaf_dgemm("N","T",ib,jb,k,ALPHA,dB,lddB,dA,lddA,BETA,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"
 				
                            endif
			    !! hC = dC
        		    !call magmaf_dgetmatrix( ib, jb, dC, lddC, hC(i,j), ldC, queue )
                	    call magmaf_dgetmatrix_async( ib, jb, dC, lddC, hC(i,j), ldC, queue )
             		ENDDO  !  i

          	ENDDO  !   j

           ELSE  !   nota
	          !  C := alpha*A**T*A + beta*C

                ldda = ceiling(real(k)/32)*32
       		lddb = ceiling(real(k)/32)*32
       		info = 0
       		!print *," Mem in dev=", ((dble(nb)*dble(nb)+dble(nb)*dble(k)+dble(nb)*dble(k))*8.0)/1e9 ," GB"
       		totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(nb)+dble(lddb)*dble(nb))*8.0)/1e9
       		!print *, "magmaf_BIG_dsyrk_cpu; Memory in device=", totM," GB"


      		!! Allocate GPU memory
        	info = magmaf_dmalloc( dA, ldda*nb )
        	if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            		goto 1000
        	endif

        	info = magmaf_dmalloc( dB, lddb*nb )
        	if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            		goto 1000
        	endif
        	!print *, "end Allocate GPU memory"
		!!  num blk


        	DO j=1,n,nb
                	jb = min( nb, n-j+1 )
			!! dA = hA(j,1)
        		!call magmaf_dsetmatrix( k, jb, hA(1,j), ldA, dA, lddA, queue )
                	call magmaf_dsetmatrix_async( k, jb, hA(1,j), ldA, dA, lddA, queue )
           		DO i=1,j,nb
                            ib = min( nb, n-i+1 )
                             		print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib
                            if(BETA/=zero)then
				!! dC = hC
        			call magmaf_dsetmatrix( ib, jb, hC(i,j), ldC, dC, lddC, queue )
                  	    endif
                            if(j==i)then
                               print *, "magmaf_dsyrk"
                               call magmaf_dsyrk( "U", "T", nb, k, alpha, dA, lddA, beta, dC, lddc, queue )

                            else
                                print *, "set B"
				!! dB = hA(1,i)
        			!call magmaf_dsetmatrix( k, jb, hA(1,i), ldA, dA, lddA, queue )
                		call magmaf_dsetmatrix_async( k, ib, hA(1,i), ldA, dB, lddB, queue )
                				
                                print *, "magmaf_dgemm"
                		call magmaf_queue_sync(queue) 	
				!call magmablasf_dgemm(TRANSA,TRANSB,m,n,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )
				call magmaf_dgemm("T","N",ib,jb,k,ALPHA,dB,lddB,dA,lddA,BETA,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"
 				
                            endif
			    !! hC = dC
        		    !call magmaf_dgetmatrix( ib, jb, dC, lddC, hC(i,j), ldC, queue )
                	    call magmaf_dgetmatrix_async( ib, jb, dC, lddC, hC(i,j), ldC, queue )
             		ENDDO  !  i

          	ENDDO  !   j

           ENDIF

	ELSE  ! lower
            !print *, "cal lower "
           IF (nota) THEN

          	! C(nb,nb) := alpha*A*A**T + beta*C(nb,nb). ; A(n,k)
		 print *, "not TESTED "
		ldda = ceiling(real(NB)/32)*32
       		lddb = ceiling(real(nb)/32)*32
       		info = 0
       		print *," Mem in dev=", ((dble(nb)*dble(nb)+dble(nb)*dble(k)+dble(nb)*dble(k))*8.0)/1e9 ," GB"
       		totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(nb)+dble(lddb)*dble(k))*8.0)/1e9
       		print *, "magmaf_BIG_dsyrk_cpu; Memory in device=", totM," GB"


      		!! Allocate GPU memory
        	info = magmaf_dmalloc( dA, ldda*k )
        	if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            		goto 1000
        	endif

        	info = magmaf_dmalloc( dB, lddb*k )
        	if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            		goto 1000
        	endif
        	!print *, "end Allocate GPU memory"
		!!  num blk


        	DO j=1,n,nb
                	jb = min( nb, n-j+1 )
			!! dA = hA(j,1)
        		!call magmaf_dsetmatrix( k, jb, hA(1,j), ldA, dA, lddA, queue )
                	call magmaf_dsetmatrix_async( jb, k, hA(j,1), ldA, dA, lddA, queue )
           		DO i=j,n,nb
                            ib = min( nb, n-i+1 )
                             		print *, "i=",i,"ib=",ib,"j=",j,"jb=",jb
                            if(BETA/=zero)then
				!! dC = hC
        			call magmaf_dsetmatrix( ib, jb, hC(i,j), ldC, dC, lddC, queue )
                  	    endif
                            if(j==i)then
                               print *, "magmaf_dsyrk"
                               call magmaf_dsyrk( "L", "N", nb, k, alpha, dA, lddA, beta, dC, lddc, queue )

                            else
                                print *, "set B"
				!! dB = hA(1,i)
        			!call magmaf_dsetmatrix( k, jb, hA(1,i), ldA, dA, lddA, queue )
                		call magmaf_dsetmatrix_async( ib, k, hA(i,1), ldA, dB, lddB, queue )
                				
                                print *, "magmaf_dgemm"
                		call magmaf_queue_sync(queue) 	
				!call magmablasf_dgemm(TRANSA,TRANSB,m,n,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )
				call magmaf_dgemm("N","T",ib,jb,k,ALPHA,dB,lddB,dA,lddA,BETA,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"
 				
                            endif
			    !! hC = dC
        		    !call magmaf_dgetmatrix( ib, jb, dC, lddC, hC(i,j), ldC, queue )
                	    call magmaf_dgetmatrix_async( ib, jb, dC, lddC, hC(i,j), ldC, queue )
             		ENDDO  !  i

          	ENDDO  !   j

       		
           ELSE  !  nota
		!print *, "C := alpha*A**T*A + beta*C"

       		ldda = ceiling(real(k)/32)*32
       		lddb = ceiling(real(k)/32)*32
       		info = 0
       		!print *," Mem in dev=", ((dble(nb)*dble(nb)+dble(nb)*dble(k)+dble(nb)*dble(k))*8.0)/1e9 ," GB"
       		totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(nb)+dble(lddb)*dble(nb))*8.0)/1e9
       		!print *, "1 magmaf_BIG_dsyrk_cpu; Memory in device=", totM," GB"


      		!! Allocate GPU memory
        	info = magmaf_dmalloc( dA, ldda*nb )
        	if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            		goto 1000
        	endif

        	info = magmaf_dmalloc( dB, lddb*nb )
        	if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            		goto 1000
        	endif
        	!print *, "end Allocate GPU memory"
		!!  num blk


        	DO j=1,n,nb
                	jb = min( nb, n-j+1 )
			!! dA = hA(j,1)
        		!call magmaf_dsetmatrix( k, jb, hA(1,j), ldA, dA, lddA, queue )
                	call magmaf_dsetmatrix_async( k, jb, hA(1,j), ldA, dA, lddA, queue )
           		DO i=j,n,nb
                            ib = min( nb, n-i+1 )
                             		!print *, "i=",i,"ib=",ib,"j=",j,"jb=",jb
                            if(BETA/=zero)then
				!! dC = hC
        			call magmaf_dsetmatrix( ib, jb, hC(i,j), ldC, dC, lddC, queue )
                  	    endif
                            if(j==i)then
                               !print *, "magmaf_dsyrk"
                               call magmaf_dsyrk( "L", "T", nb, k, alpha, dA, lddA, beta, dC, lddc, queue )

                            else
                                !print *, "set B"
				!! dB = hA(1,i)
        			!call magmaf_dsetmatrix( k, jb, hA(1,i), ldA, dA, lddA, queue )
                		call magmaf_dsetmatrix_async( k, ib, hA(1,i), ldA, dB, lddB, queue )
                				
                                !print *, "magmaf_dgemm"
                		call magmaf_queue_sync(queue) 	
				!call magmablasf_dgemm(TRANSA,TRANSB,m,n,k,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )
				call magmaf_dgemm("T","N",ib,jb,k,ALPHA,dB,lddB,dA,lddA,BETA,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"
 				
                            endif
			    !! hC = dC
        		    !call magmaf_dgetmatrix( ib, jb, dC, lddC, hC(i,j), ldC, queue )
                	    call magmaf_dgetmatrix_async( ib, jb, dC, lddC, hC(i,j), ldC, queue )
             		ENDDO  !  i

          	ENDDO  !   j

           END IF !  nota
       END IF  !upper

    ENDIF  !! full A in device

ENDIF  !!  maxmem<LimMem





 call magmaf_queue_sync(queue)
   
!! cleanup:
1000 continue
  !! Free GPU memory

     
        info = magmaf_free( dA )
        if (info .ne. 0) then
            print *, 'magmaf_BIG_dsyrk_cpu Error: magmaf_free( dA ) failed, info = ', info
            stop
        endif
        info = magmaf_free( dC )
        if (info .ne. 0) then
            print *, 'magmaf_BIG_dsyrk_cpu Error: magmaf_free( dC ) failed, info = ', info
            stop
        endif
      if(maxmem>LimMem)then
        info = magmaf_free( dB )
        if (info .ne. 0) then
            print *, 'magmaf_BIG_dsyrk_cpu Error: magmaf_free( dB ) failed, info = ', info
            stop
        endif
      endif

    call magmaf_queue_destroy( queue )

end subroutine


subroutine magmaf_dtrsm_cpu(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,hA,LDA,hB,LDB)  
!Purpose:
!
!     DTRSM  solves one of the matrix equations
!
!        op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!     where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!     non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!        op( A ) = A   or   op( A ) = A**T.
!
!     The matrix X is overwritten on B.
!


!	SIDE is CHARACTER*1
!            On entry, SIDE specifies whether op( A ) appears on the left
!            or right of X as follows:
!  
!               SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!  
!               SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!     UPLO is CHARACTER*1
!            On entry, UPLO specifies whether the matrix A is an upper or
!            lower triangular matrix as follows:
!  
!               UPLO = 'U' or 'u'   A is an upper triangular matrix.
!  
!               UPLO = 'L' or 'l'   A is a lower triangular matrix.
!     TRANSA is CHARACTER*1
!            On entry, TRANSA specifies the form of op( A ) to be used in
!            the matrix multiplication as follows:
!  
!               TRANSA = 'N' or 'n'   op( A ) = A.
!  
!               TRANSA = 'T' or 't'   op( A ) = A**T.
!  
!               TRANSA = 'C' or 'c'   op( A ) = A**T.
!    DIAG is CHARACTER*1
!            On entry, DIAG specifies whether or not A is unit triangular
!            as follows:
!  
!               DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!  
!               DIAG = 'N' or 'n'   A is not assumed to be unit
!                                   triangular.
!     M is INTEGER
!            On entry, M specifies the number of rows of B. M must be at
!            least zero.
!      N is INTEGER
!            On entry, N specifies the number of columns of B.  N must be
!            at least zero.
!      ALPHA is DOUBLE PRECISION.
!            On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!            zero then  A is not referenced and  B need not be set before
!            entry.
!     A is DOUBLE PRECISION array, dimension ( LDA, k ),
!            where k is m when SIDE = 'L' or 'l'
!              and k is n when SIDE = 'R' or 'r'.
!            Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!            upper triangular part of the array  A must contain the upper
!            triangular matrix  and the strictly lower triangular part of
!            A is not referenced.
!            Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!            lower triangular part of the array  A must contain the lower
!            triangular matrix  and the strictly upper triangular part of
!            A is not referenced.
!            Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!            A  are not referenced either,  but are assumed to be  unity.
!      LDA is INTEGER
!            On entry, LDA specifies the first dimension of A as declared
!            in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!            LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!            then LDA must be at least max( 1, n ).
!      B is DOUBLE PRECISION array, dimension ( LDB, N )
!            Before entry,  the leading  m by n part of the array  B must
!            contain  the  right-hand  side  matrix  B,  and  on exit  is
!            overwritten by the solution matrix  X.
!      LDB is INTEGER
!            On entry, LDB specifies the first dimension of B as declared
!            in  the  calling  (sub)  program.   LDB  must  be  at  least
!            max( 1, m ).

!  =====================================================================
!     .. Scalar Arguments ..
       DOUBLE PRECISION ALPHA
       INTEGER LDA,LDB,N,M
       CHARACTER DIAG,SIDE,TRANSA,UPLO
!     .. Array Arguments ..CPU
       DOUBLE PRECISION hA(LDA,*),hB(LDB,*)
!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dA,  dB
       magma_devptr_t          :: queue
       INTEGER LDDA,LDDB,dev
!  =====================================================================
!     .. Local Scalars ..
       INTEGER INFO
    LOGICAL LSIDE
    double precision                    :: one, neg_one, zero
    parameter                       (one = 1.0d0, neg_one = -1.0d0, zero=0.0d0)
!     .. External Functions ..
       LOGICAL LSAME

       EXTERNAL lsame
!    .. Local Scalars ..
       INTEGER NROWA


       lside = lsame(side,'L')
       IF (lside) THEN
           nrowa = m
       ELSE
           nrowa = n
       END IF

    lddB = ceiling(real(M)/32)*32
    ldda = ceiling(real(NROWA)/32)*32

    info = 0

    dev=0
    call magmaf_queue_create( dev, queue )

      !! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*nrowa )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif
        info = magmaf_dmalloc( dB, lddb*n )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dC  ) failed, info = ", info
            goto 1000
        endif

    	!! dA = hA
        call magmaf_dsetmatrix( nrowa, nrowa, hA, lda, dA, ldda, queue )

        if(alpha/=zero)then
	    !! dB = hB
            call magmaf_dsetmatrix( m, n, hB, ldB, dB, lddB, queue )

            call magmaf_dtrsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,dA,LDDA,dB,LDDB,queue)
        else
            call magmaf_dtrsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,dA,LDDA,dB,LDDB,queue)
        endif

 	!! hb = db
        call magmaf_dgetmatrix( m, n, dB, lddB, hB, ldB, queue )

   
!! cleanup:
1000 continue
  !! Free GPU memory
        info = magmaf_free( dA )
        if (info .ne. 0) then
            print *, 'Error: magmaf_free( dA ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dB )
        if (info .ne. 0) then
            print *, 'Error: magmaf_free( dX ) failed, info = ', info
            stop
        endif
    call magmaf_queue_destroy( queue );

end subroutine




subroutine magmaf_dsyrk_cpu(UPLO,TRANS,N,K,ALPHA,hA,LDA,BETA,hC,LDC)  
![in]	TRANS	
!          TRANS is CHARACTER*1
!           On entry,  TRANS  specifies the operation to be performed as
!           follows:
!              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.
!              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.
!              TRANS = 'C' or 'c'   C := alpha*A**T*A + beta*C.
![in]	K	
!          K is INTEGER
!           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
!           of  columns   of  the   matrix   A,   and  on   entry   with
!           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
!           of rows of the matrix  A.  K must be at least zero.
!  =====================================================================
!     .. Scalar Arguments ..
       DOUBLE PRECISION ALPHA,BETA
       INTEGER LDA,LDC,N,K
       CHARACTER UPLO,TRANS
!     .. Array Arguments ..CPU
       DOUBLE PRECISION hA(LDA,*),hC(LDC,*)
!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dA,  dC
       magma_devptr_t          :: queue
       INTEGER LDDA,LDDC,dev
!  =====================================================================
!     .. Local Scalars ..
       INTEGER INFO
    double precision                    :: one, neg_one, zero
    parameter                       (one = 1.0d0, neg_one = -1.0d0, zero=0.0d0)

      PRINT *, "Start magmaf_dsyrk_cpu"
!     Quick return if possible.
!
       IF ((n.EQ.0) .OR. (((alpha.EQ.zero).OR.(k.EQ.0)).AND. (beta.EQ.one)))THEN
           PRINT *, "Quick return magmaf_dsyrk_cpu"
           RETURN
       ENDIF

    lddc = ceiling(real(ldc)/32)*32
    ldda = ceiling(real(lda)/32)*32

    info = 0

    dev=0
    call magmaf_queue_create( dev, queue )

      !! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*k )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif
        info = magmaf_dmalloc( dC, lddc*n )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dC  ) failed, info = ", info
            goto 1000
        endif

    	!! dA = hA
        call magmaf_dsetmatrix( n, k, hA, lda, dA, ldda, queue )

        if(BETA/=zero)then
	    !! dC = hC
            call magmaf_dsetmatrix( n, n, hC, ldC, dC, lddC, queue )

            call magmaf_dsyrk(UPLO,TRANS,N,k,ALPHA,dA,lddA,BETA,dC,lddc,queue )
        else
            call magmaf_dsyrk(UPLO,TRANS,N,k,ALPHA,dA,lddA,zero,dC,lddc,queue )
        endif

 	!! hC = dC
        call magmaf_dgetmatrix( n, n, dC, lddC, hC, ldC, queue )

   
!! cleanup:
1000 continue
  !! Free GPU memory
        info = magmaf_free( dA )
        if (info .ne. 0) then
            print *, 'Error: magmaf_free( dA ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dC )
        if (info .ne. 0) then
            print *, 'Error: magmaf_free( dX ) failed, info = ', info
            stop
        endif
    call magmaf_queue_destroy( queue );

   PRINT *, "END magmaf_dsyrk_cpu"
end subroutine


     SUBROUTINE magmaf_dpotrs_cpu( UPLO, N, NRHS, hA, LDA, hB, LDB, INFO )

!     .. Scalar Arguments ..
       CHARACTER          UPLO
       INTEGER            INFO, LDA, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
       DOUBLE PRECISION   hA( LDA, * ), hB( LDB, * )
!     ..
!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dA,  dB
       magma_devptr_t          :: queue
       INTEGER LDDA,LDDB,dev
!  =====================================================================


    lddb = ceiling(real(n)/32)*32
    ldda = ceiling(real(n)/32)*32

    info = 0

    dev=0
    call magmaf_queue_create( dev, queue )


!! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*n )
        if (info .ne. 0) then
            print *, "Error: magmaf_dpotrs_cpu magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif
        info = magmaf_dmalloc( dB, lddb*nrhs )
        if (info .ne. 0) then
            print *, "Error: magmaf_dpotrs_cpu magmaf_dmalloc( dB  ) failed, info = ", info
            goto 1000
        endif

    	!! dA = hA
        call magmaf_dsetmatrix( N, N, hA, lda, dA, ldda, queue )

    	!! dB = hB
        call magmaf_dsetmatrix( N, NRHS, hB, ldb, dB, lddb, queue )

    !! Call magma solve (gpu)
     call magmaf_dpotrs_gpu(UPLO, N, NRHS, dA, LDDA, dB, LDDB, INFO)
        	if (info .ne. 0) then
            		print *, "++++++++ Error +++++++++: magmaf_dpotrs_gpu failed, info = ", info
            		stop
        	endif

      !! hB = dB
        call magmaf_dgetmatrix( N, NRHS, db, lddb, hb, ldb, queue )

	!! cleanup:
	1000 continue
  	!! Free GPU memory
        info = magmaf_free( dA )
        if (info .ne. 0) then
            print *, 'Error: dpotrs_cpu magmaf_free( dA ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dB )
        if (info .ne. 0) then
            print *, 'Error: dpotrs_cpu magmaf_free( dB ) failed, info = ', info
            stop
        endif
     call magmaf_queue_destroy( queue );
      
      END SUBROUTINE magmaf_dpotrs_cpu


     SUBROUTINE magmaf_dpotrf_cpu( UPLO, N, hA, LDA, INFO )

!     .. Scalar Arguments ..
       CHARACTER          UPLO
       INTEGER            INFO, LDA, N
!     ..
!     .. Array Arguments ..
       DOUBLE PRECISION   hA( LDA, * )
!     ..
!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dA
       magma_devptr_t          :: queue
       INTEGER LDDA,dev
!  =====================================================================


    
    ldda = ceiling(real(n)/32)*32

    info = 0

    dev=0
    call magmaf_queue_create( dev, queue )


!! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*n )
        if (info .ne. 0) then
            print *, "Error: dpotrf_cpu  magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif


    	!! dA = hA
        call magmaf_dsetmatrix( N, N, hA, lda, dA, ldda, queue )



    !! Call magma solve (gpu)
     call magmaf_dpotrf_gpu(UPLO, N, dA, LDDA,  INFO)
        	if (info .ne. 0) then
            		print *, "++++++++ Error +++++++++: magmaf_dpotrf_gpu failed, info = ", info
            		stop
        	endif

      !! hB = dB
        call magmaf_dgetmatrix( N, N, dA, lddA, hA, lda, queue )

	!! cleanup:
	1000 continue
  	!! Free GPU memory
        info = magmaf_free( dA )
        if (info .ne. 0) then
            print *, 'Error: dpotrf_cpu magmaf_free( dA ) failed, info = ', info
            stop
        endif

     call magmaf_queue_destroy( queue );
      
      END SUBROUTINE magmaf_dpotrf_cpu





      SUBROUTINE magmaf_dsymm_cpu(SIDE,UPLO,M,N,ALPHA,hA,LDA,hB,LDB,BETA,hC,LDC)  !!!  working
 ! DSYMM  performs one of the matrix-matrix operations, 
 !    SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
 ! or
 !    SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
 !             where alpha and beta are scalars,  A is a symmetric matrix and  
 !                                                B and C are  m by n matrices.
!     .. Scalar Arguments ..
       DOUBLE PRECISION ALPHA,BETA
       INTEGER LDA,LDB,LDC,M,N
       CHARACTER SIDE,UPLO
!     ..
!     .. Array Arguments ..
       DOUBLE PRECISION hA(LDA,*),hB(LDB,*),hC(LDC,*)
!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dA,dB,dC
       magma_devptr_t          :: queue
       INTEGER LDDA,LDDB,LDDC,dev
!  =====================================================================
!     .. Local Scalars ..
       INTEGER INFO,NROWA
    double precision                    :: one,  zero
    parameter                       (one = 1.0d0,  zero=0.0d0)
!  =====================================================================
!    .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame

       IF (lsame(side,'L')) THEN
           nrowa = m
       ELSE
           nrowa = n
       END IF


    lddc = ceiling(real(m)/32)*32
    ldda = ceiling(real(nrowa)/32)*32
    lddb = ceiling(real(ldb)/32)*32  !???
    info = 0

    dev=0
    call magmaf_queue_create( dev, queue )

          !! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*nrowa )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dB, lddb*n )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dB  ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dC, lddc*n )
        if (info .ne. 0) then
            print *, "Error: magmaf_dmalloc( dC  ) failed, info = ", info
            goto 1000
        endif

    	!! dA = hA
        call magmaf_dsetmatrix( nrowa, nrowa, hA, lda, dA, ldda, queue )
	!! dB = hB
        call magmaf_dsetmatrix( M, N, hB, ldB, dB, lddB, queue )

        if(BETA/=zero)then
	    !! dC = hC
            call magmaf_dsetmatrix( M, n, hC, ldC, dC, lddC, queue )

            call magmaf_dsymm(SIDE,UPLO,M,N,ALPHA,dA,LDDA,dB,LDDB,BETA,dC,LDDC,queue)
        else
            call magmaf_dsymm(SIDE,UPLO,M,N,ALPHA,dA,LDDA,dB,LDDB,zero,dC,LDDC,queue)
        endif


 	!! hC = dC
        call magmaf_dgetmatrix( m, n, dC, lddC, hC, ldC, queue )

   
!! cleanup:
1000 continue
  !! Free GPU memory
        info = magmaf_free( dA )
        if (info .ne. 0) then
            print *, 'Error: magmaf_free( dA ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dB )
        if (info .ne. 0) then
            print *, 'Error: magmaf_free( dB ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dC )
        if (info .ne. 0) then
            print *, 'Error: magmaf_free( dX ) failed, info = ', info
            stop
        endif
    call magmaf_queue_destroy( queue );


     END SUBROUTINE magmaf_dsymm_cpu



subroutine magmaf_BIG_dsymm_cpu(SIDE,UPLO,M,N,ALPHA,hA,LDA,hB,LDB,BETA,hC,LDC)  !!! working
!	DSYMM  performs one of the matrix-matrix operations
!
!	    C := alpha*A*B + beta*C,
!
!	 or
!
!	    C := alpha*B*A + beta*C,
!
!	 where alpha and beta are scalars,  A is a symmetric matrix and  B and
!	 C are  m by n matrices.

!     .. Scalar Arguments ..
       DOUBLE PRECISION ALPHA,BETA
       INTEGER LDA,LDB,LDC,M,N
       CHARACTER SIDE,UPLO
!     .. Array Arguments ..CPU
       DOUBLE PRECISION hA(LDA,*),hB(LDB,*),hC(LDC,*)
        !! ARRAY WORK
        DOUBLE PRECISION, allocatable, dimension(:,:) :: TEMP
!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dA, dB, dC
       magma_devptr_t          :: queue
       INTEGER LDDA,LDDB,LDDC,dev
!  =====================================================================
!     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame
!     .. Local Scalars ..
       LOGICAL UPPER
       INTEGER INFO,lddmb
       INTEGER freeA
       INTEGER*8 :: mem_size
       !INTEGER :: h_groups, v_groups, height, width
       double precision                    ::  zero, ONE
       parameter                       ( zero=0.0d0 , ONE=1.0D0)
       double precision                    :: LimMem,totM, MaxMem
       integer                             :: j,i,jb,ib,NROWA !,NCOLA,NROWB
       integer                             :: nb,mb	!,maxnb
!  =====================================================================

    !print *, "magmaf_BIG_dsymm_cpu "
    dev=0
    call magmaf_queue_create( dev, queue )

    UPPER = lsame(UPLO,'U')
       IF (lsame(side,'L')) THEN
           nrowa = m
       ELSE
           nrowa = n
       END IF

      freeA=0
    mem_size=magmaf_mem_size( queue )    !!    mem_size=11 821 973 504  integer
      !mem_size=1921515520_8  ! 2GB   max size 6000*6000
      !mem_size=120197350_8  !!  max size 1500*1500
    !print *, "		mem_size=",mem_size
    !print*, "		huge(mem_size)=",huge(mem_size)
    !print*, "		huge(ldda)=",huge(ldda)
    !mem_size=int8(0.3d0*(1024**3))
    !mem_size=int8(3000)
    LimMem=dble(mem_size)*0.90d0 !   GB  90% of mex men
    MaxMem=( dble(m)*dble(n) + dble(nrowa)*dble(nrowa) +dble(m)*dble(n) )*8.0d0 
    !print *, "		memory full problem= ",MaxMem, "limit mem device=",LimMem

    !print *, "		mem_size=",mem_size
     
    !call matrixPartitioner(M,  N,  K, mem_size, h_groups, v_groups,height, width)

    ! print *, "h_groups=",h_groups, "v_groups=",v_groups,"height=",height, "width=",width

     !print *,"Partitioner Mem in dev=", ((dble(width)*dble(height)+dble(height)*dble(k)+dble(width)*dble(k))*8.0)/1e9 ," GB"

    if(maxmem<LimMem)then
       !! full problem
       !print *, "		Full problem"
        lddc = ceiling(real(m)/32)*32
        ldda = ceiling(real(m)/32)*32
        lddb = ceiling(real(M)/32)*32
        info = 0



      !! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*M )
        if (info .ne. 0) then
            print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dB, lddb*n )
        if (info .ne. 0) then
            print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dB  ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dC, lddc*n )
        if (info .ne. 0) then
            print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dC  ) failed, info = ", info
            goto 1000
        endif
		!print *, "end dmalloc dA,dB,dC"
        if(ALPHA/=zero)then
    	    !! dA = hA
            !call magmaf_dsetmatrix( m, k, hA, lda, dA, ldda, queue )
            call magmaf_dsetmatrix_async( m, M, hA, lda, dA, ldda, queue )
	    !! dB = hB
            !call magmaf_dsetmatrix( k, n, hB, ldB, dB, lddB, queue )
            call magmaf_dsetmatrix_async( M, n, hB, ldB, dB, lddB, queue )
        ENDIF
        if(BETA/=zero)then
		!! dC = hC
        	!call magmaf_dsetmatrix( m, n, hC, ldC, dC, lddC, queue )
                call magmaf_dsetmatrix_async( m, n, hC, ldC, dC, lddC, queue )
         endif
        call magmaf_queue_sync(queue)
	!call magmablasf_dSYmm(SIDE,UPLO,m,n,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue  )
        call magmaf_dSYmm(SIDE,UPLO,m,n,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )

 	!! hC = dC
        call magmaf_dgetmatrix( m, n, dC, lddC, hC, ldC, queue )
        !call magmaf_dgetmatrix_async( m, n, dC, lddC, hC, ldC, queue )

    ELSE   !!  	by block,

        if( (dble(m)*dble(512) + dble(nrowa)*dble(nrowa) +dble(M)*dble(512)) *8.0d0 <LimMem)then   !!  by colums block
        !!   full A in device and nb=>512,  =>> C(m,nb)=C(m,nb)+a(m,M)*b(M,nb)
           !print *, "BIG_dsymm_cpu full A in device and nb=>512"
           IF (lsame(side,'L')) THEN    
	     !!  C := alpha*A*B + beta*C.    A SYM
               nb=(LimMem/8-nrowa*nrowa)/(nrowa+nrowa)-127
               nb = ceiling(real(nb)/128)*128
               print *, "full A in device and nb=",nb
               lddc = ceiling(real(m)/32)*32
               ldda = ceiling(real(m)/32)*32
               lddb = ceiling(real(m)/32)*32
       
           !print *,"BIG_dsymm_cpu Mem blk=",((dble(m)*dble(nb)+dble(nrowa)*dble(nrowa)+dble(nrowa)*dble(nb))*8.0)/1e9," GB"       
           !print *,"BIG_dsymm_cpu Mem in device=",((dble(lddc)*dble(nb)+dble(ldda)*dble(nrowa)+dble(lddb)*dble(nb))*8.0)/1e9," GB"
                  !! Allocate GPU memory
               info = magmaf_dmalloc( dA, ldda*m )
               if (info .ne. 0) then
                  print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA ) failed, info = ", info
                  goto 1000
               endif
               print *, "dA=",ldda*m*8/1e9,"GB"
              info = magmaf_dmalloc( dB, lddb*nb )
              if (info .ne. 0) then
                  print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dB  ) failed, info = ", info
                  goto 1000
              endif

              info = magmaf_dmalloc( dC, lddc*nb )
              if (info .ne. 0) then
                 print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dC  ) failed, info = ", info
                 goto 1000
              endif
              !print *, "end Allocate GPU memory"
	      !!  num blk

	      !! dA = hA   !! full A
              call magmaf_dsetmatrix( m, M, hA, lda, dA, ldda, queue )


              DO j=1,n,nb
                jb = min( nb, n-j+1 )
		!! dB = hB
        	!call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                call magmaf_dsetmatrix_async( M, jb, hB(1,j), ldB, dB, lddB, queue )

                		!print *, "j=",j,"jb=",jb
        	if(BETA/=zero)then
			!! dC = hC
        		!call magmaf_dsetmatrix( m, jb, hC(1,j), ldC, dC, lddC, queue )
                        call magmaf_dsetmatrix_async( m, jb, hC(1,j), ldC, dC, lddC, queue )
                 endif
                call magmaf_queue_sync(queue) 	
		!call magmablasf_dSYmm("L",UPLO,m,jb,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddC,queue )
		call magmaf_dsymm("L",UPLO,m,jb,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddC,queue )
              		!print *, "C(1:m,j:j+jb)=C(1:m",j,":",j+jb-1,")*A(1:m,1:k)*B(1:k,",j,":",j+jb-1,")"
 		!! hC = dC
        	!call magmaf_dgetmatrix( m, jb, dC, lddC, hC(1,j), ldC, queue )
                call magmaf_dgetmatrix_async( m, jb, dC, lddC, hC(1,j), ldC, queue )


              ENDDO
           ELSE  !  side R
          !!  C := alpha*B*A + beta*C.  A is SYM
             PRINT *, " NOT IMPLEMENT YET, full A in device"
             PRINT *, "  C := alpha*B*A + beta*C.  A is SYM"
           ENDIF  
    
        ELSE   !!  BY BLK(mb,nb)
             
            !!   cal C by C(mb,nb)=C(mb,nb)+A(mb,m)*B(m,nb) or  C(mb,nb)=C(mb,nb)+B(mb,n)*A(n,nb)
            print *, "BIG_dsymm_cpu	BY BLK(mb,nb)"

           !call magmaf_get_dgemm_mb_nb(limMem,m,n,nrowa,mb,nb)
	   call magmaf_get_BIG_dsymm_mb_nb(limmem,SIDE,m,n,mb,nb)
             !mb=2000; nb=mb
            !mb=min(mb,nb)
            !nb=mb
            ALLOCATE(TEMP(MB,NB))
           !call matrixPartitioner(M,  N,  K, int8(limMem), h_groups, v_groups,mb, nb)
           print *, "		BIG_dsymm; mb", mb," nb",nb

           lddc = ceiling(real(mb)/32)*32
           !ldda = ceiling(real(mb)/32)*32
           lddB = ceiling(real(m)/32)*32

           !! Allocate GPU memory  
           info = magmaf_dmalloc( dC, lddc*nb )
           if (info .ne. 0) then
             print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dC  ) failed, info = ", info
             goto 1000
           endif
           !print *, "alloco dC; lddc*nb", lddc,nb,lddc*nb," GB=",(dble(lddc)*dble(nb)*8.0d0)/1e9


           IF (lsame(side,'L')) THEN
               print *, "	LEFT"
		          ! C := alpha*A*B + beta*C,   A sym
               IF (upper) THEN
                   print *, "UPPER"
       		   !ldda = ceiling(real(m)/32)*32
       		   lddb = ceiling(real(m)/32)*32
       		   info = 0
       		   print *," Mem in dev=", ((dble(mb)*dble(nb)+dble(mb)*dble(nrowa)+dble(nb)*dble(nrowa))*8.0)/1e9 ," GB"
       		   totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(mb)+dble(lddb)*dble(mb))*8.0D0)/1e9
       		   print *, "BIG_dsymm; Memory in device=", totM," GB"

		   !!  num blk

      		   info = magmaf_dmalloc( dB, lddb*mb )
      		   if (info .ne. 0) then
            		print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dB  ) failed, info = ", info
            		goto 1000
      		   endif
      		   !print *, "end Allocate GPU memory"
                   !print *, "alloco dB; lddb*mb", lddb,mb,lddb*mb," GB=",(dble(lddb)*dble(mb)*8.0d0)/1e9
        	   DO j=1,n,nb
                	jb = min( nb, n-j+1 )
           		DO i=1,m,mb
    			    ib = min( mb, m-i+1 )
                            !if(BETA/=zero)then
			    !    !! dC = hC
        	    	    !    call magmaf_dsetmatrix( ib, jb, hC(i,j), ldC, dC, lddC, queue )
                 	    !endif
                            call magmaf_dscal( LDDC*NB, 0.0d0, dC, 1, queue )

                            if(i-1>1)then
                                 !print *, "rec^T  (i-1,mb)        =>  j=",j,"jb=",jb,"i=",i,"ib=",ib,"i-1",i-1
                                 ldda = ceiling(real(i-1)/32)*32
                                 !print *, "ldda = ceiling(real(i-1)/32)*32=",ldda
      		   		info = magmaf_dmalloc( dA, ldda*mb )
      		   		if (info .ne. 0) then
            				print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA  ) failed, info = ", info
            				goto 1000
      		   		endif
                                !print *, "alloco dA; ldda*mb", ldda,mb,ldda*mb," GB=",(dble(ldda)*dble(mb)*8.0d0)/1e9
                                !print *,"    rec^T hA(",1,":",I,",",i,":",I+IB-1,")"
				!! dA = hA
        			!call magmaf_dsetmatrix( ib, k, hA(i,1), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( i-1,ib, hA(1,i), lda, dA, ldda, queue )
                			!print *,"    rec^T hB(",1,":",I,",",j,":",J+JB-1,")"	
                                !! dB = hB
        		        !call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                	        call magmaf_dsetmatrix_async( i-1, jb, hB(1,j), ldB, dB, lddB, queue )

                		call magmaf_queue_sync(queue) 	
				call magmaf_dgemm('T','N',ib,jb,i-1,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"

                                info = magmaf_free( dA )
        			if (info .ne. 0) then
            				print *, 'magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            				stop
        			endif
                                !print *, "free     dA; ldda*mb", ldda,mb,ldda*mb," GB=",(dble(ldda)*dble(mb)*8.0d0)/1e9
                                freeA=1
			    endif

			    ldda = ceiling(real(mb)/32)*32
      		   	    info = magmaf_dmalloc( dA, ldda*mb )
      		   	    if (info .ne. 0) then
            			print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA  ) failed, info = ", info
            			goto 1000
      	       		    endif
                            !print *, "alloco dA; ldda*mb", ldda,mb,ldda*mb," GB=",(dble(ldda)*dble(mb)*8.0d0)/1e9
                            !print *, "square (mb,mb)         =>  j=",j,"jb=",jb,"i=",i,"ib=",ib
				                                ! square
                            print *, "  sqrt hA(",i,":",I+IB-1,",",i,":",I+IB-1,")"
			   !! dA = hA
        		   !call magmaf_dsetmatrix( ib, mb, hA(i,1), lda, dA, ldda, queue )
               		   call magmaf_dsetmatrix_async( ib,ib, hA(i,i), lda, dA, ldda, queue )
                               			!! dB = hB
                           !print *, "    sqrt hB(",i,":",I+IB-1,",",j,":",J+JB-1,")"
        		   !call magmaf_dsetmatrix( ib, jb, hB(1,j), ldB, dB, lddB, queue )
                	   call magmaf_dsetmatrix_async( ib, jb, hB(i,j), ldB, dB, lddB, queue )
                				!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib

                	   call magmaf_queue_sync(queue) 	
			   call magmaf_dsymm('L','U',ib,jb,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"

                           info = magmaf_free( dA )
        		   if (info .ne. 0) then
            			print *, 'magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            			stop
        		   endif
                           
                           freeA=1
				!print *, " m-mb=",m-mb
			    if(i<m-mb)then
                                !print *, "rec    (mb,m-(i+mb)+1) =>  j=",j,"jb=",jb,"i=",i,"ib=",ib,"m-(i+mb)+1",m-(i+mb)+1
                                 ldda = ceiling(real(mb)/32)*32
      		   		info = magmaf_dmalloc( dA, ldda*(m-(i+mb)+1) )
      		   		if (info .ne. 0) then
            				print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA ) failed, info = ", info
            				goto 1000
      		   		endif
                                    !print *,"  rec hA(",i,":",I+IB-1,",",i+mb,":",i+mb+m-(i+mb)+1-1,")"
 				! rect
                                !! dA = hA
        			!call magmaf_dsetmatrix( ib, mb, hA(i,1), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( ib,m-(i+mb)+1, hA(i,i+mb), lda, dA, ldda, queue )
                                    !print *,"  rec hB(",i+mb,":",I+MB+m-(i+mb)+1-1,",",j,":",J+JB-1,")"
                                			!! dB = hB
        		        !call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                	        call magmaf_dsetmatrix_async( m-(i+mb)+1, jb, hB(i+mb,j), ldB, dB, lddB, queue )
                				!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib

                		call magmaf_queue_sync(queue) 	
				call magmaf_dgemm('N','N',ib,jb,m-(i+mb)+1,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )

                                info = magmaf_free( dA )
        			if (info .ne. 0) then
            				print *, 'magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            				stop
        			endif
                                freeA=1
                            endif
				!print *, "++get C(",i,",",j,")"
                             				!! hC = dC
        		    !call magmaf_dgetmatrix( ib, jb, dC, lddC, TEMP(1,1), MB, queue )
                	    call magmaf_dgetmatrix_async( ib, jb, dC, lddC, TEMP(1,1), MB, queue )
                            hC(i:i+ib-1,j:j+jb-1)=temp(1:ib,1:jb)+beta*hC(i:i+ib-1,j:j+jb-1)
                                !print *,""
             		ENDDO
			!print *,""
			!print *,""
                        !print *,""
          	   ENDDO
               ELSE  !! LEFT , LOWER
                   print *, "		LOWER"
       		   !ldda = ceiling(real(m)/32)*32
       		   lddb = ceiling(real(m)/32)*32  !! max size  m*nb
		   lddmb = ceiling(real(mb)/32)*32
       		   info = 0
       		   print *," 		Mem in dev=", ((dble(mb)*dble(nb)+dble(mb)*dble(nrowa)+dble(nb)*dble(nrowa))*8.0)/1e9 ," GB"
       		   totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(nrowa)+dble(lddb)*dble(mb))*8.0D0)/1e9
       		   print *, "		BIG_dsymm; Memory in device=??", totM," GB"

		   !!  num blk                         LEFT , LOWER

      		   info = magmaf_dmalloc( dB, lddb*nb ) !! max size B m*nb
      		   if (info .ne. 0) then
            		print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dB  ) failed, info = ", info
            		goto 1000
      		   endif
                   !print *, "alloco dB; lddb*nb", lddb,nb,lddb*nb," GB=",(dble(lddb)*dble(nb)*8.0d0)/1e9

 		   info = magmaf_dmalloc( dA, lddb*lddmb)  ! max size of A mb*b or m*mb
      		   if (info .ne. 0) then
            		print *, "Error: magmaf_dmalloc( dA  ) failed, info = ", info
            		goto 1000
      		   endif
		   !print *, "alloco dA; lddb*mb", lddb,mb,lddb*mb," GB=",(dble(lddb)*dble(mb)*8.0d0)/1e9
				totM=((dble(lddc)*dble(nb)+dble(lddb)*dble(mb)+dble(lddb)*dble(nb))*8.0D0)/1e9
       		   		!print *, "1 BIG_dsymm; Memory in device=??", totM," GB"

      		   !print *, "end Allocate GPU memory"

        	   DO j=1,n,nb			!!!!!!!!!!!   LEFT , LOWER
                	jb = min( nb, n-j+1 )
           		DO i=1,m,mb
    			    ib = min( mb, m-i+1 )
				!print *, ""
                                !print *, ""
                                !print *, "j=",j," jb=",jb," i=",i," ib=",ib
                            !if(BETA/=zero)then
			    !    !! dC = hC
        	    	    !    call magmaf_dsetmatrix( ib, jb, hC(i,j), ldC, dC, lddC, queue )
                 	    !endif
                            call magmaf_dscal( LDDC*NB, 0.0d0, dC, 1, queue )

                            if(i-1>1)then   !! dgemm 
                                 !print *, "rec    (i-1,mb)        =>  j=",j,"jb=",jb,"i=",i,"ib=",ib,"i-1",i-1
                                 ldda = ceiling(real(mb)/32)*32
				!print *, "1 ldda*i=",ldda,i,ldda*i
      		   		
                                	!print *,"    rec   hA(",1,":",I,",",i,":",I+IB-1,")"
				!! dA = hA
        			!call magmaf_dsetmatrix( ib,i-1, hA(i,1), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( ib,i-1, hA(i,1), lda, dA, ldda, queue )
                			!print *,"    rec   hB(",1,":",I,",",j,":",J+JB-1,")"	
                                !! dB = hB
        		        !call magmaf_dsetmatrix( i-1, jb, hB(1,j), ldB, dB, lddB, queue )
                	        call magmaf_dsetmatrix_async( i-1, jb, hB(1,j), ldB, dB, lddB, queue )

                		call magmaf_queue_sync(queue) 	
				call magmaf_dgemm('N','N',ib,jb,i-1,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"

                               
			    endif
										!!! LEFT , LOWER
			    ldda = ceiling(real(mb)/32)*32
				!print *, "2 ldda*mb=",ldda,mb,ldda*mb
      		   	    
                            !print *, "square (mb,mb)         =>  j=",j,"jb=",jb,"i=",i,"ib=",ib
				                                ! square
                            !print *, "  sqrt hA(",i,":",I+IB-1,",",i,":",I+IB-1,")"
			   !! dA = hA
        		   !call magmaf_dsetmatrix( ib,ib, hA(i,i), lda, dA, ldda, queue )
               		   call magmaf_dsetmatrix_async( ib,ib, hA(i,i), lda, dA, ldda, queue )
                               			!! dB = hB
                           !print *, "     hB(",i,":",I+IB-1,",",j,":",J+JB-1,")"
        		   !call magmaf_dsetmatrix( ib, jb, hB(i,j), ldB, dB, lddB, queue )
                	   call magmaf_dsetmatrix_async( ib, jb, hB(i,j), ldB, dB, lddB, queue )
                				!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib

                	   call magmaf_queue_sync(queue) 	
			   call magmaf_dsymm('L','L',ib,jb,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )
              !print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"

                          
				!print *, " m-mb=",m-mb
			    if(i<m-mb)then
                                !print *, "rec    (mb,m-(i+mb)+1) =>  j=",j,"jb=",jb,"i=",i,"ib=",ib,"m-(i+mb)+1",m-(i+mb)+1
                                 ldda = ceiling(real((m-(i+mb)+1))/32)*32
                                 !print *, "3 ldda = ceiling(real((m-(i+mb)+1))/32)*32=",ldda
				!print *, "3 ldda*mb=",ldda,mb,ldda*mb
      		   		
                                    !print *,"  rec hA(",i,":",I+IB-1,",",i+mb,":",i+mb+m-(i+mb)+1-1,")"
 				! rect
                                !! dA = hA
        			!call magmaf_dsetmatrix( m-(i+mb)+1,ib, hA(i+mb,i), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( m-(i+mb)+1,ib, hA(i+mb,i), lda, dA, ldda, queue )
                                    !print *,"  rec hB(",i+mb,":",I+MB+m-(i+mb)+1-1,",",j,":",J+JB-1,")"
                                			!! dB = hB
        		        !call magmaf_dsetmatrix( m-(i+mb)+1, jb, hB(i+mb,j), ldB, dB, lddB, queue )
                	        call magmaf_dsetmatrix_async( m-(i+mb)+1, jb, hB(i+mb,j), ldB, dB, lddB, queue )
                				!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib

                		call magmaf_queue_sync(queue) 	
				call magmaf_dgemm('T','N',ib,jb,m-(i+mb)+1,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )

                                
                            endif
				!print *, "++get C(",i,":",i+ib-1,",",j,":",j+jb-1,")"
                             				!! hC = dC
        		    !call magmaf_dgetmatrix( ib, jb, dC, lddC, TEMP(1,1), MB, queue )
                	    call magmaf_dgetmatrix_async( ib, jb, dC, lddC, TEMP(1,1), MB, queue )
                            hC(i:i+ib-1,j:j+jb-1)=temp(1:ib,1:jb)+beta*hC(i:i+ib-1,j:j+jb-1)
                                !print *,""
             		ENDDO   !i
			!print *,""
			!print *,""
                        !print *,""
          	   ENDDO  !j 
                                   !!  end  LEFT , LOWER
               ENDIF
	   ELSE  !side  R   not yet
          
          	! C := alpha*B*A + beta*C

      

           
           END IF  !side

        ENDIF  !! full A in device

    ENDIF  !!  maxmem<LimMem





    call magmaf_queue_sync(queue)
   
!! cleanup:
1000 continue
  !! Free GPU memory

        ! print *, "freeA=",freeA
        info = magmaf_free( dA )
        if (info .ne. 0) then
            print *, 'Magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dB )
        if (info .ne. 0) then
            print *, 'Magmaf_BIG_dsymm_cpu Error: magmaf_free( dB ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dC )
        if (info .ne. 0) then
            print *, 'Magmaf_BIG_dsymm_cpu Error: magmaf_free( dC ) failed, info = ', info
            stop
        endif
    call magmaf_queue_destroy( queue )

end subroutine magmaf_BIG_dsymm_cpu


subroutine magmaf_get_BIG_dsymm_mb_nb(limmem,SIDE,m,n,mb,nb)

implicit none
!  =====================================================================
!     .. Scalar Arguments ..
       CHARACTER SIDE
   DOUBLE PRECISION :: limMem
   integer          :: m,n,nrowa
   integer          :: mb,nb, nbmax,Lddm,lddr
   logical          :: okay
   DOUBLE PRECISION :: Mem,tem

!     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame


       IF (lsame(side,'L')) THEN
           nrowa = m
       ELSE
           nrowa = n
       END IF

	!! test full problem size
        Lddm=ceiling(real(m)/128)*128
        Lddr=ceiling(real(nrowa)/128)*128
        Mem = ( dble(Lddm)*dble(n) + dble(Lddm)*dble(nrowa) +dble(Lddr)*dble(n) )*8.0d0 
	if(mem<limmem)then
		print *, "non necesary OutOfCore sub"
		return
        endif


		!! first, asuming mb=nb
   !! c(mb,nb)=c(mb,nb)+A(mb,nrowa)*B(nrowa,nb)  in GPU  mb=nb
    !! c(mb,nb)=c(mb,nb)+B(nb,nrowa)*A(nrowa,mb)  in GPU  mb=nb

    nbmax=int(sqrt(dble(nrowa)**2+limMem/8.0d0)-dble(nrowa) )
    print *, "nbmax=",nbmax
    nb=nbmax; mb=nb

    if(nbmax>n)then
	nb=n
        !!increse mb
	call get_mem_dgemm(mb,nb,nrowa,Mem)
	do while(mem < limmem)
		mb=mb+1
		call get_mem_dgemm(mb+128,nb+128,nrowa,Mem)
        enddo
    endif
    if(nbmax>m)then
	mb=m
	!!increse nb
	call get_mem_dgemm(mb,nb,nrowa,Mem)
	do while(mem < limmem)
		nb=nb+1
		call get_mem_dgemm(mb+128,nb+128,nrowa,Mem)
        enddo
    endif



okay=.false.
    nb = ceiling(real(nb)/128)*128
    mb = ceiling(real(mb)/128)*128

do while(okay .eqv. .false.)
	


    call get_mem_dgemm(mb,nb,nrowa,Mem)
    if(mem<limmem)then
         !print *, "c(mb,nb)=c(mb,nb)+A(mb,k)*B(k,nb)  in GPU  mb=nb"
         !print *, "mb=",mb,"nb=",nb
         ! print *, " men dgemm in dev",Mem
         okay=.true.
         goto 1000
    else
	nb=nb-128
	mb=mb-128

    endif

1000 continue

end do



     !print *, "magmaf_get_BIG_dsymm_mb_nb  mb=",mb,"nb=",nb
    call get_mem_dgemm(mb,nb,nrowa,tem)
      !print *, " men dgemm in dev",tem



   return
end subroutine magmaf_get_BIG_dsymm_mb_nb






    subroutine make_dsym(UPLO,n,A)
    
    implicit none
    integer :: n
    character :: UPLO
    DOUBLE PRECISION :: A(n,n)
   !  =====================================================================
!     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame
!     .. Local Scalars ..
       LOGICAL            UPPER
       integer :: i,j


    upper = lsame( uplo, 'U' )


       if(upper)then

	   do j=1,n
              do i=j,n
                 A(i,j)=A(j,i)
              enddo
           enddo

       else

	   do j=1,n
              do i=1,j
                 A(i,j)=A(j,i)
              enddo
           enddo

       endif


    end subroutine 




    subroutine test_dsym(n,A,info)


    implicit none
    integer :: n
    DOUBLE PRECISION :: A(n,n)
   !  =====================================================================
       integer :: i,j,info

     info=0
	   do j=1,n
              do i=1,j-1
                 if(A(i,j)/=A(j,i))then
		     info=info+1
		 endif
              enddo
           enddo

        if(info/=0)then
	   print *," Error, matrix not symetric, info=",info
        else
           print *," Matrix are symetric, info=",info," sum(A)=",sum(A)
        end if
    end subroutine


subroutine magmaf_BIG_dsymm_cpu_old(SIDE,UPLO,M,N,ALPHA,hA,LDA,hB,LDB,BETA,hC,LDC)  !!! working
!	DSYMM  performs one of the matrix-matrix operations
!
!	    C := alpha*A*B + beta*C,
!
!	 or
!
!	    C := alpha*B*A + beta*C,
!
!	 where alpha and beta are scalars,  A is a symmetric matrix and  B and
!	 C are  m by n matrices.

!     .. Scalar Arguments ..
       DOUBLE PRECISION ALPHA,BETA
       INTEGER LDA,LDB,LDC,M,N
       CHARACTER SIDE,UPLO
!     .. Array Arguments ..CPU
       DOUBLE PRECISION hA(LDA,*),hB(LDB,*),hC(LDC,*)
        !! ARRAY WORK
        DOUBLE PRECISION, allocatable, dimension(:,:) :: TEMP
!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dA, dB, dC
       magma_devptr_t          :: queue
       INTEGER LDDA,LDDB,LDDC,dev
!  =====================================================================
!     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame
!     .. Local Scalars ..
       LOGICAL UPPER
       INTEGER INFO
       INTEGER freeA
       INTEGER*8 :: mem_size
       !INTEGER :: h_groups, v_groups, height, width
       double precision                    ::  zero, ONE
       parameter                       ( zero=0.0d0 , ONE=1.0D0)
       double precision                    :: LimMem,totM, MaxMem
       integer                             :: j,i,jb,ib,NROWA !,NCOLA,NROWB
       integer                             :: nb,mb !,maxnb
!  =====================================================================

    print *, "magmaf_BIG_dsymm_cpu "
    dev=0
    call magmaf_queue_create( dev, queue )

    UPPER = lsame(UPLO,'U')
       IF (lsame(side,'L')) THEN
           nrowa = m
       ELSE
           nrowa = n
       END IF

      freeA=0
    mem_size=magmaf_mem_size( queue )    !!    mem_size=11 821 973 504  integer
    !  mem_size=1921515520_8  ! 2GB   max size 6000*6000
    ! mem_size=120197350_8  !!  max size 1500*1500


    print *, "mem_size=",mem_size
    print*, "huge(mem_size)=",huge(mem_size)
    print*, "huge(ldda)=",huge(ldda)
    !mem_size=int8(0.3d0*(1024**3))
    !mem_size=int8(3000)
    LimMem=dble(mem_size)*0.90d0 !   GB  90% of mex men
    MaxMem=( dble(m)*dble(n) + dble(nrowa)*dble(nrowa) +dble(m)*dble(n) )*8.0d0 
    print *, "memory full problem= ",MaxMem, "limit memory device=",LimMem

    print *, "mem_size=",mem_size
     
    !call matrixPartitioner(M,  N,  K, mem_size, h_groups, v_groups,height, width)

    ! print *, "h_groups=",h_groups, "v_groups=",v_groups,"height=",height, "width=",width

     !print *,"Partitioner Mem in dev=", ((dble(width)*dble(height)+dble(height)*dble(k)+dble(width)*dble(k))*8.0)/1e9 ," GB"

    if(maxmem<LimMem)then
       !! full problem
       !print *, "Full problem"
        lddc = ceiling(real(m)/32)*32
        ldda = ceiling(real(m)/32)*32
        lddb = ceiling(real(M)/32)*32
        info = 0



      !! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*M )
        if (info .ne. 0) then
            print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dB, lddb*n )
        if (info .ne. 0) then
            print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dB  ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dC, lddc*n )
        if (info .ne. 0) then
            print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dC  ) failed, info = ", info
            goto 1000
        endif
		!print *, "end dmalloc dA,dB,dC"
        if(ALPHA/=zero)then
    	    !! dA = hA
            !call magmaf_dsetmatrix( m, k, hA, lda, dA, ldda, queue )
            call magmaf_dsetmatrix_async( m, M, hA, lda, dA, ldda, queue )
	    !! dB = hB
            !call magmaf_dsetmatrix( k, n, hB, ldB, dB, lddB, queue )
            call magmaf_dsetmatrix_async( M, n, hB, ldB, dB, lddB, queue )
        ENDIF
        if(BETA/=zero)then
		!! dC = hC
        	!call magmaf_dsetmatrix( m, n, hC, ldC, dC, lddC, queue )
                call magmaf_dsetmatrix_async( m, n, hC, ldC, dC, lddC, queue )
         endif
        call magmaf_queue_sync(queue)
	!call magmablasf_dSYmm(SIDE,UPLO,m,n,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue  )
        call magmaf_dSYmm(SIDE,UPLO,m,n,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )

 	!! hC = dC
        call magmaf_dgetmatrix( m, n, dC, lddC, hC, ldC, queue )
        !call magmaf_dgetmatrix_async( m, n, dC, lddC, hC, ldC, queue )

    ELSE   !!  	by block,

        if( (dble(m)*dble(512) + dble(nrowa)*dble(nrowa) +dble(M)*dble(512)) *8.0d0 <LimMem)then   !!  by colums block
        !!   full A in device and nb=>512,  =>> C(m,nb)=C(m,nb)+a(m,M)*b(M,nb)
           print *, "full A in device and nb=>512"
           IF (lsame(side,'L')) THEN    
	     !!  C := alpha*A*B + beta*C.    A SYM
               nb=(LimMem/8-nrowa*nrowa)/(nrowa+nrowa)-127
               nb = ceiling(real(nb)/128)*128
               !print *, "full A in device and nb=",nb
               lddc = ceiling(real(m)/32)*32
               ldda = ceiling(real(m)/32)*32
               lddb = ceiling(real(m)/32)*32
       
               print *,"Mem blk=", ((dble(m)*dble(nb)+dble(nrowa)*dble(nrowa)+dble(nrowa)*dble(nb))*8.0)/1e9 ," GB"       
               print *,"Mem in device=", ((dble(lddc)*dble(nb)+dble(ldda)*dble(nrowa)+dble(lddb)*dble(nb))*8.0)/1e9 ," GB"
                  !! Allocate GPU memory
               info = magmaf_dmalloc( dA, ldda*m )
               if (info .ne. 0) then
                  print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA ) failed, info = ", info
                  goto 1000
               endif
               !print *, "dA=",ldda*m*8/1e9,"GB"
              info = magmaf_dmalloc( dB, lddb*nb )
              if (info .ne. 0) then
                  print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dB  ) failed, info = ", info
                  goto 1000
              endif

              info = magmaf_dmalloc( dC, lddc*nb )
              if (info .ne. 0) then
                 print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dC  ) failed, info = ", info
                 goto 1000
              endif
              !print *, "end Allocate GPU memory"
	      !!  num blk

	      !! dA = hA   !! full A
              call magmaf_dsetmatrix( m, M, hA, lda, dA, ldda, queue )


              DO j=1,n,nb
                jb = min( nb, n-j+1 )
		!! dB = hB
        	!call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                call magmaf_dsetmatrix_async( M, jb, hB(1,j), ldB, dB, lddB, queue )

                		!print *, "j=",j,"jb=",jb
        	if(BETA/=zero)then
			!! dC = hC
        		!call magmaf_dsetmatrix( m, jb, hC(1,j), ldC, dC, lddC, queue )
                        call magmaf_dsetmatrix_async( m, jb, hC(1,j), ldC, dC, lddC, queue )
                 endif
                call magmaf_queue_sync(queue) 	
		!call magmablasf_dSYmm("L",UPLO,m,jb,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddC,queue )
		call magmaf_dsymm("L",UPLO,m,jb,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddC,queue )
              		!print *, "C(1:m,j:j+jb)=C(1:m",j,":",j+jb-1,")*A(1:m,1:k)*B(1:k,",j,":",j+jb-1,")"
 		!! hC = dC
        	!call magmaf_dgetmatrix( m, jb, dC, lddC, hC(1,j), ldC, queue )
                call magmaf_dgetmatrix_async( m, jb, dC, lddC, hC(1,j), ldC, queue )


              ENDDO
           ELSE  !  side R
          !!  C := alpha*B*A + beta*C.  A is SYM
             PRINT *, " NOT IMPLEMENT YET, full A in device"
             PRINT *, "  C := alpha*B*A + beta*C.  A is SYM"
           ENDIF  
    
        ELSE   !!  BY BLK(mb,nb)
             
            !!   cal C by C(mb,nb)=C(mb,nb)+op(a())*op(b())
            print *, "BY BLK(mb,nb)"

           call magmaf_get_dgemm_mb_nb(limMem,m,n,nrowa,mb,nb)
             !mb=2000; nb=mb
            !mb=min(mb,nb)
            !nb=mb
            ALLOCATE(TEMP(MB,NB))
           !call matrixPartitioner(M,  N,  K, int8(limMem), h_groups, v_groups,mb, nb)
           print *, "BIG_dsymm; mb", mb," nb",nb

           lddc = ceiling(real(mb)/32)*32
           !ldda = ceiling(real(mb)/32)*32
           lddB = ceiling(real(m)/32)*32

           !! Allocate GPU memory  
           info = magmaf_dmalloc( dC, lddc*nb )
           if (info .ne. 0) then
             print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dC  ) failed, info = ", info
             goto 1000
           endif
           print *, "alloco dC; lddc*nb", lddc,nb,lddc*nb," GB=",(dble(lddc)*dble(nb)*8.0d0)/1e9


           IF (lsame(side,'L')) THEN
               print *, "LEFT"
		          ! C := alpha*A*B + beta*C,   A sym
               IF (upper) THEN   !!   "LEFT"    "UPPER"
                   print *, "UPPER"
       		   !ldda = ceiling(real(m)/32)*32
       		   lddb = ceiling(real(m)/32)*32
       		   info = 0
       		   print *," Mem in dev=", ((dble(mb)*dble(nb)+dble(mb)*dble(nrowa)+dble(nb)*dble(nrowa))*8.0)/1e9 ," GB"
       		   totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(mb)+dble(lddb)*dble(mb))*8.0D0)/1e9
       		   print *, "BIG_dsymm; Memory in device=", totM," GB"

		   !!  num blk

      		   info = magmaf_dmalloc( dB, lddb*mb )
      		   if (info .ne. 0) then
            		print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dB  ) failed, info = ", info
            		goto 1000
      		   endif
      		   !print *, "end Allocate GPU memory"
                   print *, "alloco dB; lddb*mb", lddb,mb,lddb*mb," GB=",(dble(lddb)*dble(mb)*8.0d0)/1e9
        	   DO j=1,n,nb
                	jb = min( nb, n-j+1 )
           		DO i=1,m,mb
    			    ib = min( mb, m-i+1 )
                            !if(BETA/=zero)then
			    !    !! dC = hC
        	    	    !    call magmaf_dsetmatrix( ib, jb, hC(i,j), ldC, dC, lddC, queue )
                 	    !endif
                            call magmaf_dscal( LDDC*MB, 0.0d0, dC, 1, queue )

                            if(i-1>1)then
                                 print *, "rec^T  (i-1,mb)        =>  j=",j,"jb=",jb,"i=",i,"ib=",ib,"i-1",i-1
                                 ldda = ceiling(real(i-1)/32)*32
                                 print *, "ldda = ceiling(real(i-1)/32)*32=",ldda
      		   		info = magmaf_dmalloc( dA, ldda*mb )
      		   		if (info .ne. 0) then
            				print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA  ) failed, info = ", info
            				goto 1000
      		   		endif
                                print *, "alloco dA; ldda*mb", ldda,mb,ldda*mb," GB=",(dble(ldda)*dble(mb)*8.0d0)/1e9
                                print *,"    rec^T hA(",1,":",I,",",i,":",I+IB-1,")"
				!! dA = hA
        			!call magmaf_dsetmatrix( ib, k, hA(i,1), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( i-1,ib, hA(1,i), lda, dA, ldda, queue )
                			print *,"    rec^T hB(",1,":",I,",",j,":",J+JB-1,")"	
                                !! dB = hB
        		        !call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                	        call magmaf_dsetmatrix_async( i-1, jb, hB(1,j), ldB, dB, lddB, queue )

                		call magmaf_queue_sync(queue) 	
				call magmaf_dgemm('T','N',ib,jb,i-1,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"

                                info = magmaf_free( dA )
        			if (info .ne. 0) then
            				print *, 'magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            				stop
        			endif
                                print *, "free     dA; ldda*mb", ldda,mb,ldda*mb," GB=",(dble(ldda)*dble(mb)*8.0d0)/1e9
                                freeA=1
			    endif

			    ldda = ceiling(real(mb)/32)*32
      		   	    info = magmaf_dmalloc( dA, ldda*mb )
      		   	    if (info .ne. 0) then
            			print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA  ) failed, info = ", info
            			goto 1000
      	       		    endif
                            print *, "alloco dA; ldda*mb", ldda,mb,ldda*mb," GB=",(dble(ldda)*dble(mb)*8.0d0)/1e9
                            print *, "square (mb,mb)         =>  j=",j,"jb=",jb,"i=",i,"ib=",ib
				                                ! square
                            print *, "  sqrt hA(",i,":",I+IB-1,",",i,":",I+IB-1,")"
			   !! dA = hA
        		   !call magmaf_dsetmatrix( ib, mb, hA(i,1), lda, dA, ldda, queue )
               		   call magmaf_dsetmatrix_async( ib,ib, hA(i,i), lda, dA, ldda, queue )
                               			!! dB = hB
                           print *, "    sqrt hB(",i,":",I+IB-1,",",j,":",J+JB-1,")"
        		   !call magmaf_dsetmatrix( ib, jb, hB(1,j), ldB, dB, lddB, queue )
                	   call magmaf_dsetmatrix_async( ib, jb, hB(i,j), ldB, dB, lddB, queue )
                				!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib

                	   call magmaf_queue_sync(queue) 	
			   call magmaf_dsymm('L','U',ib,jb,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"

                           info = magmaf_free( dA )
        		   if (info .ne. 0) then
            			print *, 'magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            			stop
        		   endif
                           
                           freeA=1
				!print *, " m-mb=",m-mb
			    if(i<m-mb)then
                                print *, "rec    (mb,m-(i+mb)+1) =>  j=",j,"jb=",jb,"i=",i,"ib=",ib,"m-(i+mb)+1",m-(i+mb)+1
                                 ldda = ceiling(real(mb)/32)*32
      		   		info = magmaf_dmalloc( dA, ldda*(m-(i+mb)+1) )
      		   		if (info .ne. 0) then
            				print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA ) failed, info = ", info
            				goto 1000
      		   		endif
                                    print *,"  rec hA(",i,":",I+IB-1,",",i+mb,":",i+mb+m-(i+mb)+1-1,")"
 				! rect
                                !! dA = hA
        			!call magmaf_dsetmatrix( ib, mb, hA(i,1), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( ib,m-(i+mb)+1, hA(i,i+mb), lda, dA, ldda, queue )
                                    print *,"  rec hB(",i+mb,":",I+MB+m-(i+mb)+1-1,",",j,":",J+JB-1,")"
                                			!! dB = hB
        		        !call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                	        call magmaf_dsetmatrix_async( m-(i+mb)+1, jb, hB(i+mb,j), ldB, dB, lddB, queue )
                				!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib

                		call magmaf_queue_sync(queue) 	
				call magmaf_dgemm('N','N',ib,jb,m-(i+mb)+1,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )

                                info = magmaf_free( dA )
        			if (info .ne. 0) then
            				print *, 'magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            				stop
        			endif
                                freeA=1
                            endif
				print *, "++get C(",i,",",j,")"
                             				!! hC = dC
        		    !call magmaf_dgetmatrix( ib, jb, dC, lddC, hC(i,j), ldC, queue )
                	    call magmaf_dgetmatrix_async( ib, jb, dC, lddC, TEMP(1,1), MB, queue )
                            hC(i:i+ib-1,j:j+jb-1)=temp(1:ib,1:jb)+beta*hC(i:i+ib-1,j:j+jb-1)
                                !print *,""
             		ENDDO
			!print *,""
			!print *,""
                        !print *,""
          	   ENDDO
               ELSE  !!  "LEFT"  LOWER
                   print *, "LOWER"
       		   !ldda = ceiling(real(m)/32)*32
       		   lddb = ceiling(real(m)/32)*32  !! max size  m*nb
       		   info = 0
       		   print *," Mem in dev=", ((dble(mb)*dble(nb)+dble(mb)*dble(nrowa)+dble(nb)*dble(nrowa))*8.0)/1e9 ," GB"
       		   totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(nrowa)+dble(lddb)*dble(mb))*8.0D0)/1e9
       		   print *, "BIG_dsymm; Memory in device=??", totM," GB"

		   !!  num blk

      		   info = magmaf_dmalloc( dB, lddb*nb )
      		   if (info .ne. 0) then
            		print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dB  ) failed, info = ", info
            		goto 1000
      		   endif
                   print *, "alloco dB; lddb*mb", lddb,nb,lddb*nb," GB=",(dble(lddb)*dble(nb)*8.0d0)/1e9
      		   !print *, "end Allocate GPU memory"

        	   DO j=1,n,nb
                	jb = min( nb, n-j+1 )
           		DO i=1,m,mb
    			    ib = min( mb, m-i+1 )
				print *, ""
                                print *, ""
                                print *, "j=",j," jb=",jb," i=",i," ib=",ib
                            !if(BETA/=zero)then
			    !    !! dC = hC
        	    	    !    call magmaf_dsetmatrix( ib, jb, hC(i,j), ldC, dC, lddC, queue )
                 	    !endif
                            call magmaf_dscal( LDDC*nB, 0.0d0, dC, 1, queue )

                            if(i-1>1)then   !! dgemm 
                                 print *, "rec    (i-1,mb)        =>  j=",j,"jb=",jb,"i=",i,"ib=",ib,"i-1",i-1
                                 ldda = ceiling(real(mb)/32)*32
				print *, "1 ldda*i=",ldda,i,ldda*i
      		   		info = magmaf_dmalloc( dA, ldda*i)
      		   		if (info .ne. 0) then
            				print *, "Error: magmaf_dmalloc( dA  ) failed, info = ", info
            				goto 1000
      		   		endif
                                print *, "alloco dA; ldda*i", ldda,i,lddb*i," GB=",(dble(ldda)*dble(i)*8.0d0)/1e9
				totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(i)+dble(lddb)*dble(mb))*8.0D0)/1e9
       		   		print *, "1 BIG_dsymm; Memory in device=??", totM," GB"
                                	print *,"    rec   hA(",1,":",I,",",i,":",I+IB-1,")"
				!! dA = hA
        			!call magmaf_dsetmatrix( ib, k, hA(i,1), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( ib,i-1, hA(i,1), lda, dA, ldda, queue )
                			print *,"    rec   hB(",1,":",I,",",j,":",J+JB-1,")"	
                                !! dB = hB
        		        !call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                	        call magmaf_dsetmatrix_async( i-1, jb, hB(1,j), ldB, dB, lddB, queue )

                		call magmaf_queue_sync(queue) 	
				call magmaf_dgemm('N','N',ib,jb,i-1,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )
              			print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"

                                info = magmaf_free( dA )
        			if (info .ne. 0) then
            				print *, '1 magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            				stop
        			endif
                                print *, "free dA"
                                freeA=1
			    endif
				
			    ldda = ceiling(real(mb)/32)*32
				print *, "2 ldda*mb=",ldda,mb,ldda*mb
      		   	    info = magmaf_dmalloc( dA, ldda*mb )
      		   	    if (info .ne. 0) then
            			print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA  ) failed, info = ", info
            			goto 1000
      	       		    endif
				print *, "alloco dA; ldda*mb", ldda,mb,ldda*mb," GB=",(dble(ldda)*dble(nb)*8.0d0)/1e9
                            	totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(mb)+dble(lddb)*dble(nb))*8.0D0)/1e9
       		   		print *, "2 BIG_dsymm; Memory in device=??", totM," GB"
                            !print *, "square (mb,mb)         =>  j=",j,"jb=",jb,"i=",i,"ib=",ib
				                                ! square
                            print *, "  sqrt hA(",i,":",I+IB-1,",",i,":",I+IB-1,")"
			   !! dA = hA
        		   !call magmaf_dsetmatrix( ib, mb, hA(i,1), lda, dA, ldda, queue )
               		   call magmaf_dsetmatrix_async( ib,ib, hA(i,i), lda, dA, ldda, queue )
                               			!! dB = hB
                           print *, "     hB(",i,":",I+IB-1,",",j,":",J+JB-1,")"
        		   !call magmaf_dsetmatrix( ib, jb, hB(1,j), ldB, dB, lddB, queue )
                	   call magmaf_dsetmatrix_async( ib, jb, hB(i,j), ldB, dB, lddB, queue )
                				!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib

                	   call magmaf_queue_sync(queue) 	
			   call magmaf_dsymm('L','L',ib,jb,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"

                           info = magmaf_free( dA )
        		   if (info .ne. 0) then
            			print *, '2 magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            			stop
        		   endif
                           print *, "free dA"
                           freeA=1
				!print *, " m-mb=",m-mb
			    if(i<m-mb)then
                                !print *, "rec    (mb,m-(i+mb)+1) =>  j=",j,"jb=",jb,"i=",i,"ib=",ib,"m-(i+mb)+1",m-(i+mb)+1
                                 ldda = ceiling(real((m-(i+mb)+1))/32)*32
                                 print *, "3 ldda = ceiling(real((m-(i+mb)+1))/32)*32=",ldda
				print *, "3 ldda*mb=",ldda,mb,ldda*mb
      		   		info = magmaf_dmalloc( dA, ldda*mb )
      		   		if (info .ne. 0) then
            				print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            				goto 1000
      		   		endif
					print *, "alloco dA; ldda*mb", ldda,mb,ldda*mb," GB=",(dble(ldda)*dble(mb)*8.0d0)/1e9
                                       totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(mb)+dble(lddb)*dble(mb))*8.0D0)/1e9
       		   		print *, "3 BIG_dsymm; Memory in device=??", totM," GB"
                                    print *,"  rec hA(",i,":",I+IB-1,",",i+mb,":",i+mb+m-(i+mb)+1-1,")"
 				! rect
                                !! dA = hA
        			!call magmaf_dsetmatrix( ib, mb, hA(i,1), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( m-(i+mb)+1,ib, hA(i+mb,i), lda, dA, ldda, queue )
                                    !print *,"  rec hB(",i+mb,":",I+MB+m-(i+mb)+1-1,",",j,":",J+JB-1,")"
                                			!! dB = hB
        		        !call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                	        call magmaf_dsetmatrix_async( m-(i+mb)+1, jb, hB(i+mb,j), ldB, dB, lddB, queue )
                				!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib

                		call magmaf_queue_sync(queue) 	
				call magmaf_dgemm('T','N',ib,jb,m-(i+mb)+1,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )

                                info = magmaf_free( dA )
        			if (info .ne. 0) then
            				print *, '3 magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            				stop
        			endif
      				print *, "free dA"
                                freeA=1
                            endif
				print *, "++get C(",i,":",i+ib-1,",",j,":",j+jb-1,")"
                             				!! hC = dC
        		    !call magmaf_dgetmatrix( ib, jb, dC, lddC, hC(i,j), ldC, queue )
                	    call magmaf_dgetmatrix_async( ib, jb, dC, lddC, TEMP(1,1), MB, queue )
                            hC(i:i+ib-1,j:j+jb-1)=temp(1:ib,1:jb)+beta*hC(i:i+ib-1,j:j+jb-1)
                                print *,""
             		ENDDO   !i
			print *,""
			print *,""
                        print *,""
          	   ENDDO  !j 

               ENDIF
	   ELSE  !side  R
          
		PRINT *, " NOT IMPLEMENT YET, "

           END IF  !side

        ENDIF  !! full A in device

    ENDIF  !!  maxmem<LimMem





    call magmaf_queue_sync(queue)
   
!! cleanup:
1000 continue
  !! Free GPU memory
     if(freeA<1)then
                print *, "freeA=",freeA
        info = magmaf_free( dA )
        if (info .ne. 0) then
            print *, 'Magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            stop
        endif
     endif
        info = magmaf_free( dB )
        if (info .ne. 0) then
            print *, 'Magmaf_BIG_dsymm_cpu Error: magmaf_free( dB ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dC )
        if (info .ne. 0) then
            print *, 'Magmaf_BIG_dsymm_cpu Error: magmaf_free( dC ) failed, info = ', info
            stop
        endif
    call magmaf_queue_destroy( queue )

     print *, "END magmaf_BIG_dsymm_cpu_old"
end subroutine magmaf_BIG_dsymm_cpu_old


subroutine magmaf_BIG_dsymm_cpu_1(SIDE,UPLO,M,N,ALPHA,hA,LDA,hB,LDB,BETA,hC,LDC)  !!! working
!	DSYMM  performs one of the matrix-matrix operations
!
!	    C := alpha*A*B + beta*C,
!
!	 or
!
!	    C := alpha*B*A + beta*C,
!
!	 where alpha and beta are scalars,  A is a symmetric matrix and  B and
!	 C are  m by n matrices.

!     .. Scalar Arguments ..
       DOUBLE PRECISION ALPHA,BETA
       INTEGER LDA,LDB,LDC,M,N
       CHARACTER SIDE,UPLO
!     .. Array Arguments ..CPU
       DOUBLE PRECISION hA(LDA,*),hB(LDB,*),hC(LDC,*)
        !! ARRAY WORK
        DOUBLE PRECISION, allocatable, dimension(:,:) :: TEMP
!  =====================================================================
!     .. Device Array Arguments ..GPU
       magma_devptr_t          :: dA, dB, dC
       magma_devptr_t          :: queue
       INTEGER LDDA,LDDB,LDDC,dev
!  =====================================================================
!     .. External Functions ..
       LOGICAL LSAME
       EXTERNAL lsame
!     .. Local Scalars ..
       LOGICAL UPPER
       INTEGER INFO
       INTEGER freeA
       INTEGER*8 :: mem_size
       !INTEGER ::  h_groups, v_groups, height, width
       double precision                    ::  zero, ONE
       parameter                       ( zero=0.0d0 , ONE=1.0D0)
       double precision                    :: LimMem,totM, MaxMem
       integer                             :: j,i,jb,ib,NROWA !,NCOLA,NROWB
       integer                             :: nb,mb  !,maxnb
!  =====================================================================

    print *, "magmaf_BIG_dsymm_cpu "
    dev=0
    call magmaf_queue_create( dev, queue )

    UPPER = lsame(UPLO,'U')
       IF (lsame(side,'L')) THEN
           nrowa = m
       ELSE
           nrowa = n
       END IF

      freeA=0
    mem_size=magmaf_mem_size( queue )    !!    mem_size=11 821 973 504  integer
    !  mem_size=1921515520_8  ! 2GB   max size 6000*6000
    !mem_size=120197350_8  !!  max size 1500*1500
    print *, "mem_size=",mem_size
    print*, "huge(mem_size)=",huge(mem_size)
    print*, "huge(ldda)=",huge(ldda)
    !mem_size=int8(0.3d0*(1024**3))
    !mem_size=int8(3000)
    LimMem=dble(mem_size)*0.90d0 !   GB  90% of mex men
    MaxMem=( dble(m)*dble(n) + dble(nrowa)*dble(nrowa) +dble(m)*dble(n) )*8.0d0 
    print *, "memory full problem= ",MaxMem, "limit memory device=",LimMem

    print *, "mem_size=",mem_size
     
    !call matrixPartitioner(M,  N,  K, mem_size, h_groups, v_groups,height, width)

    ! print *, "h_groups=",h_groups, "v_groups=",v_groups,"height=",height, "width=",width

     !print *,"Partitioner Mem in dev=", ((dble(width)*dble(height)+dble(height)*dble(k)+dble(width)*dble(k))*8.0)/1e9 ," GB"

    if(maxmem<LimMem)then
       !! full problem
       !print *, "Full problem"
        lddc = ceiling(real(m)/32)*32
        ldda = ceiling(real(m)/32)*32
        lddb = ceiling(real(M)/32)*32
        info = 0



      !! Allocate GPU memory
        info = magmaf_dmalloc( dA, ldda*M )
        if (info .ne. 0) then
            print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dB, lddb*n )
        if (info .ne. 0) then
            print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dB  ) failed, info = ", info
            goto 1000
        endif

        info = magmaf_dmalloc( dC, lddc*n )
        if (info .ne. 0) then
            print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dC  ) failed, info = ", info
            goto 1000
        endif
		!print *, "end dmalloc dA,dB,dC"
        if(ALPHA/=zero)then
    	    !! dA = hA
            !call magmaf_dsetmatrix( m, k, hA, lda, dA, ldda, queue )
            call magmaf_dsetmatrix_async( m, M, hA, lda, dA, ldda, queue )
	    !! dB = hB
            !call magmaf_dsetmatrix( k, n, hB, ldB, dB, lddB, queue )
            call magmaf_dsetmatrix_async( M, n, hB, ldB, dB, lddB, queue )
        ENDIF
        if(BETA/=zero)then
		!! dC = hC
        	!call magmaf_dsetmatrix( m, n, hC, ldC, dC, lddC, queue )
                call magmaf_dsetmatrix_async( m, n, hC, ldC, dC, lddC, queue )
         endif
        call magmaf_queue_sync(queue)
	!call magmablasf_dSYmm(SIDE,UPLO,m,n,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue  )
        call magmaf_dSYmm(SIDE,UPLO,m,n,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddc,queue )

 	!! hC = dC
        call magmaf_dgetmatrix( m, n, dC, lddC, hC, ldC, queue )
        !call magmaf_dgetmatrix_async( m, n, dC, lddC, hC, ldC, queue )

    ELSE   !!  	by block,

        if( (dble(m)*dble(512) + dble(nrowa)*dble(nrowa) +dble(M)*dble(512)) *8.0d0 <LimMem)then   !!  by colums block
        !!   full A in device and nb=>512,  =>> C(m,nb)=C(m,nb)+a(m,M)*b(M,nb)
           print *, "full A in device and nb=>512"
           IF (lsame(side,'L')) THEN    
	     !!  C := alpha*A*B + beta*C.    A SYM
               nb=(LimMem/8-nrowa*nrowa)/(nrowa+nrowa)-127
               nb = ceiling(real(nb)/128)*128
               print *, "full A in device and nb=",nb
               lddc = ceiling(real(m)/32)*32
               ldda = ceiling(real(m)/32)*32
               lddb = ceiling(real(m)/32)*32
       
               print *,"Mem blk=", ((dble(m)*dble(nb)+dble(nrowa)*dble(nrowa)+dble(nrowa)*dble(nb))*8.0)/1e9 ," GB"       
               print *,"Mem in device=", ((dble(lddc)*dble(nb)+dble(ldda)*dble(nrowa)+dble(lddb)*dble(nb))*8.0)/1e9 ," GB"
                  !! Allocate GPU memory
               info = magmaf_dmalloc( dA, ldda*m )
               if (info .ne. 0) then
                  print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA ) failed, info = ", info
                  goto 1000
               endif
               print *, "dA=",ldda*m*8/1e9,"GB"
              info = magmaf_dmalloc( dB, lddb*nb )
              if (info .ne. 0) then
                  print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dB  ) failed, info = ", info
                  goto 1000
              endif

              info = magmaf_dmalloc( dC, lddc*nb )
              if (info .ne. 0) then
                 print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dC  ) failed, info = ", info
                 goto 1000
              endif
              !print *, "end Allocate GPU memory"
	      !!  num blk

	      !! dA = hA   !! full A
              call magmaf_dsetmatrix( m, M, hA, lda, dA, ldda, queue )


              DO j=1,n,nb
                jb = min( nb, n-j+1 )
		!! dB = hB
        	!call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                call magmaf_dsetmatrix_async( M, jb, hB(1,j), ldB, dB, lddB, queue )

                		!print *, "j=",j,"jb=",jb
        	if(BETA/=zero)then
			!! dC = hC
        		!call magmaf_dsetmatrix( m, jb, hC(1,j), ldC, dC, lddC, queue )
                        call magmaf_dsetmatrix_async( m, jb, hC(1,j), ldC, dC, lddC, queue )
                 endif
                call magmaf_queue_sync(queue) 	
		!call magmablasf_dSYmm("L",UPLO,m,jb,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddC,queue )
		call magmaf_dsymm("L",UPLO,m,jb,ALPHA,dA,ldda,dB,lddb,BETA,dC,lddC,queue )
              		!print *, "C(1:m,j:j+jb)=C(1:m",j,":",j+jb-1,")*A(1:m,1:k)*B(1:k,",j,":",j+jb-1,")"
 		!! hC = dC
        	!call magmaf_dgetmatrix( m, jb, dC, lddC, hC(1,j), ldC, queue )
                call magmaf_dgetmatrix_async( m, jb, dC, lddC, hC(1,j), ldC, queue )


              ENDDO
           ELSE  !  side R
          !!  C := alpha*B*A + beta*C.  A is SYM
             PRINT *, " NOT IMPLEMENT YET, full A in device"
             PRINT *, "  C := alpha*B*A + beta*C.  A is SYM"
           ENDIF  
    
        ELSE   !!  BY BLK(mb,nb)
             
            !!   cal C by C(mb,nb)=C(mb,nb)+op(a())*op(b())
            print *, "BY BLK(mb,nb)"

           call magmaf_get_dgemm_mb_nb(limMem,m,n,nrowa,mb,nb)
             !mb=2000; nb=mb
            !mb=min(mb,nb)
            !nb=mb
            ALLOCATE(TEMP(MB,NB))
           !call matrixPartitioner(M,  N,  K, int8(limMem), h_groups, v_groups,mb, nb)
           print *, "BIG_dsymm; mb", mb," nb",nb

           lddc = ceiling(real(mb)/32)*32
           !ldda = ceiling(real(mb)/32)*32
           lddB = ceiling(real(m)/32)*32

           !! Allocate GPU memory  
           info = magmaf_dmalloc( dC, lddc*nb )
           if (info .ne. 0) then
             print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dC  ) failed, info = ", info
             goto 1000
           endif
           print *, "alloco dC; lddc*nb", lddc,nb,lddc*nb," GB=",(dble(lddc)*dble(nb)*8.0d0)/1e9


           IF (lsame(side,'L')) THEN
               print *, "LEFT"
		          ! C := alpha*A*B + beta*C,   A sym
               IF (upper) THEN   !!   "LEFT"    "UPPER"
                   print *, "UPPER"
       		   !ldda = ceiling(real(m)/32)*32
       		   lddb = ceiling(real(m)/32)*32
       		   info = 0
       		   print *," Mem in dev=", ((dble(mb)*dble(nb)+dble(mb)*dble(nrowa)+dble(nb)*dble(nrowa))*8.0)/1e9 ," GB"
       		   totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(mb)+dble(lddb)*dble(mb))*8.0D0)/1e9
       		   print *, "BIG_dsymm; Memory in device=", totM," GB"

		   !!  num blk

      		   info = magmaf_dmalloc( dB, lddb*mb )
      		   if (info .ne. 0) then
            		print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dB  ) failed, info = ", info
            		goto 1000
      		   endif
      		   !print *, "end Allocate GPU memory"
                   print *, "alloco dB; lddb*mb", lddb,mb,lddb*mb," GB=",(dble(lddb)*dble(mb)*8.0d0)/1e9
        	   DO j=1,n,nb
                	jb = min( nb, n-j+1 )
           		DO i=1,m,mb
    			    ib = min( mb, m-i+1 )
                            !if(BETA/=zero)then
			    !    !! dC = hC
        	    	    !    call magmaf_dsetmatrix( ib, jb, hC(i,j), ldC, dC, lddC, queue )
                 	    !endif
                            call magmaf_dscal( LDDC*MB, 0.0d0, dC, 1, queue )

                            if(i-1>1)then
                                 print *, "rec^T  (i-1,mb)        =>  j=",j,"jb=",jb,"i=",i,"ib=",ib,"i-1",i-1
                                 ldda = ceiling(real(i-1)/32)*32
                                 print *, "ldda = ceiling(real(i-1)/32)*32=",ldda
      		   		info = magmaf_dmalloc( dA, ldda*mb )
      		   		if (info .ne. 0) then
            				print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA  ) failed, info = ", info
            				goto 1000
      		   		endif
                                print *, "alloco dA; ldda*mb", ldda,mb,ldda*mb," GB=",(dble(ldda)*dble(mb)*8.0d0)/1e9
                                print *,"    rec^T hA(",1,":",I,",",i,":",I+IB-1,")"
				!! dA = hA
        			!call magmaf_dsetmatrix( ib, k, hA(i,1), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( i-1,ib, hA(1,i), lda, dA, ldda, queue )
                			print *,"    rec^T hB(",1,":",I,",",j,":",J+JB-1,")"	
                                !! dB = hB
        		        !call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                	        call magmaf_dsetmatrix_async( i-1, jb, hB(1,j), ldB, dB, lddB, queue )

                		call magmaf_queue_sync(queue) 	
				call magmaf_dgemm('T','N',ib,jb,i-1,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"

                                info = magmaf_free( dA )
        			if (info .ne. 0) then
            				print *, 'magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            				stop
        			endif
                                print *, "free     dA; ldda*mb", ldda,mb,ldda*mb," GB=",(dble(ldda)*dble(mb)*8.0d0)/1e9
                                freeA=1
			    endif

			    ldda = ceiling(real(mb)/32)*32
      		   	    info = magmaf_dmalloc( dA, ldda*mb )
      		   	    if (info .ne. 0) then
            			print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA  ) failed, info = ", info
            			goto 1000
      	       		    endif
                            print *, "alloco dA; ldda*mb", ldda,mb,ldda*mb," GB=",(dble(ldda)*dble(mb)*8.0d0)/1e9
                            print *, "square (mb,mb)         =>  j=",j,"jb=",jb,"i=",i,"ib=",ib
				                                ! square
                            print *, "  sqrt hA(",i,":",I+IB-1,",",i,":",I+IB-1,")"
			   !! dA = hA
        		   !call magmaf_dsetmatrix( ib, mb, hA(i,1), lda, dA, ldda, queue )
               		   call magmaf_dsetmatrix_async( ib,ib, hA(i,i), lda, dA, ldda, queue )
                               			!! dB = hB
                           print *, "    sqrt hB(",i,":",I+IB-1,",",j,":",J+JB-1,")"
        		   !call magmaf_dsetmatrix( ib, jb, hB(1,j), ldB, dB, lddB, queue )
                	   call magmaf_dsetmatrix_async( ib, jb, hB(i,j), ldB, dB, lddB, queue )
                				!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib

                	   call magmaf_queue_sync(queue) 	
			   call magmaf_dsymm('L','U',ib,jb,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"

                           info = magmaf_free( dA )
        		   if (info .ne. 0) then
            			print *, 'magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            			stop
        		   endif
                           
                           freeA=1
				!print *, " m-mb=",m-mb
			    if(i<m-mb)then
                                print *, "rec    (mb,m-(i+mb)+1) =>  j=",j,"jb=",jb,"i=",i,"ib=",ib,"m-(i+mb)+1",m-(i+mb)+1
                                 ldda = ceiling(real(mb)/32)*32
      		   		info = magmaf_dmalloc( dA, ldda*(m-(i+mb)+1) )
      		   		if (info .ne. 0) then
            				print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA ) failed, info = ", info
            				goto 1000
      		   		endif
                                    print *,"  rec hA(",i,":",I+IB-1,",",i+mb,":",i+mb+m-(i+mb)+1-1,")"
 				! rect
                                !! dA = hA
        			!call magmaf_dsetmatrix( ib, mb, hA(i,1), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( ib,m-(i+mb)+1, hA(i,i+mb), lda, dA, ldda, queue )
                                    print *,"  rec hB(",i+mb,":",I+MB+m-(i+mb)+1-1,",",j,":",J+JB-1,")"
                                			!! dB = hB
        		        !call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                	        call magmaf_dsetmatrix_async( m-(i+mb)+1, jb, hB(i+mb,j), ldB, dB, lddB, queue )
                				!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib

                		call magmaf_queue_sync(queue) 	
				call magmaf_dgemm('N','N',ib,jb,m-(i+mb)+1,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )

                                info = magmaf_free( dA )
        			if (info .ne. 0) then
            				print *, 'magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            				stop
        			endif
                                freeA=1
                            endif
				print *, "++get C(",i,",",j,")"
                             				!! hC = dC
        		    !call magmaf_dgetmatrix( ib, jb, dC, lddC, hC(i,j), ldC, queue )
                	    call magmaf_dgetmatrix_async( ib, jb, dC, lddC, TEMP(1,1), MB, queue )
                            hC(i:i+ib-1,j:j+jb-1)=temp(1:ib,1:jb)+beta*hC(i:i+ib-1,j:j+jb-1)
                                print *,""
             		ENDDO
			print *,""
			print *,""
                        print *,""
          	   ENDDO
               ELSE  !!  "LEFT"  LOWER
                   print *, "LOWER"
       		   !ldda = ceiling(real(m)/32)*32
       		   lddb = ceiling(real(m)/32)*32  !! max size  m*nb
       		   info = 0
       		   print *," Mem in dev=", ((dble(mb)*dble(nb)+dble(mb)*dble(nrowa)+dble(nb)*dble(nrowa))*8.0)/1e9 ," GB"
       		   totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(nrowa)+dble(lddb)*dble(mb))*8.0D0)/1e9
       		   print *, "BIG_dsymm; Memory in device=??", totM," GB"

		   !!  num blk

      		   info = magmaf_dmalloc( dB, lddb*nb )
      		   if (info .ne. 0) then
            		print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dB  ) failed, info = ", info
            		goto 1000
      		   endif
                   print *, "alloco dB; lddb*nb", lddb,nb,lddb*nb," GB=",(dble(lddb)*dble(nb)*8.0d0)/1e9
                   info = magmaf_dmalloc( dA, lddb*mb)
      		   		if (info .ne. 0) then
            				print *, "Error: magmaf_dmalloc( dA  ) failed, info = ", info
            				goto 1000
      		   		endif
                                print *, "alloco dA; lddb*mb", lddb,mb,lddb*mb," GB=",(dble(lddb)*dble(mb)*8.0d0)/1e9
				totM=((dble(lddc)*dble(nb)+dble(lddb)*dble(mb)+dble(lddb)*dble(nb))*8.0D0)/1e9
       		   		print *, "1 BIG_dsymm; Memory in device=??", totM," GB"

      		   !print *, "end Allocate GPU memory"

        	   DO j=1,n,nb
                	jb = min( nb, n-j+1 )
           		DO i=1,m,mb
    			    ib = min( mb, m-i+1 )
				print *, ""
                                print *, ""
                                print *, "j=",j," jb=",jb," i=",i," ib=",ib
                            !if(BETA/=zero)then
			    !    !! dC = hC
        	    	    !    call magmaf_dsetmatrix( ib, jb, hC(i,j), ldC, dC, lddC, queue )
                 	    !endif
                            call magmaf_dscal( LDDC*nB, 0.0d0, dC, 1, queue )

                            if(i-1>1)then   !! dgemm 
                                 print *, "rec    (i-1,mb)        =>  j=",j,"jb=",jb,"i=",i,"ib=",ib,"i-1",i-1
                                 ldda = ceiling(real(mb)/32)*32
				print *, "1 ldda*i=",ldda,i,ldda*i
      		   		!info = magmaf_dmalloc( dA, ldda*i)
      		   		!if (info .ne. 0) then
            		!		print *, "Error: magmaf_dmalloc( dA  ) failed, info = ", info
            		!		goto 1000
      		   	!	endif
                         !       print *, "alloco dA; ldda*i", ldda,i,lddb*i," GB=",(dble(ldda)*dble(i)*8.0d0)/1e9
			!	totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(i)+dble(lddb)*dble(mb))*8.0D0)/1e9
       		   	!	print *, "1 BIG_dsymm; Memory in device=??", totM," GB"
                                	print *,"    rec   hA(",1,":",I,",",i,":",I+IB-1,")"
				!! dA = hA
        			!call magmaf_dsetmatrix( ib, k, hA(i,1), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( ib,i-1, hA(i,1), lda, dA, ldda, queue )
                			print *,"    rec   hB(",1,":",I,",",j,":",J+JB-1,")"	
                                !! dB = hB
        		        !call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                	        call magmaf_dsetmatrix_async( i-1, jb, hB(1,j), ldB, dB, lddB, queue )

                		call magmaf_queue_sync(queue) 	
				call magmaf_dgemm('N','N',ib,jb,i-1,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )
              			print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"

                              !  info = magmaf_free( dA )
        		!	if (info .ne. 0) then
            		!		print *, '1 magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            		!		stop
        		!	endif
                         !       print *, "free dA"
                          !      freeA=1
			    endif
				
			    ldda = ceiling(real(mb)/32)*32
				print *, "2 ldda*mb=",ldda,mb,ldda*mb
      		   	    !info = magmaf_dmalloc( dA, ldda*mb )
      		   	    !if (info .ne. 0) then
            		!	print *, "BIG_dsymm_cpu Error: magmaf_dmalloc( dA  ) failed, info = ", info
            		!	goto 1000
      	       		 !   endif
			!	print *, "alloco dA; ldda*mb", ldda,mb,ldda*mb," GB=",(dble(ldda)*dble(nb)*8.0d0)/1e9
                         !   	totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(mb)+dble(lddb)*dble(nb))*8.0D0)/1e9
       		   	!	print *, "2 BIG_dsymm; Memory in device=??", totM," GB"
                            !print *, "square (mb,mb)         =>  j=",j,"jb=",jb,"i=",i,"ib=",ib
				                                ! square
                            print *, "  sqrt hA(",i,":",I+IB-1,",",i,":",I+IB-1,")"
			   !! dA = hA
        		   !call magmaf_dsetmatrix( ib, mb, hA(i,1), lda, dA, ldda, queue )
               		   call magmaf_dsetmatrix_async( ib,ib, hA(i,i), lda, dA, ldda, queue )
                               			!! dB = hB
                           print *, "     hB(",i,":",I+IB-1,",",j,":",J+JB-1,")"
        		   !call magmaf_dsetmatrix( ib, jb, hB(1,j), ldB, dB, lddB, queue )
                	   call magmaf_dsetmatrix_async( ib, jb, hB(i,j), ldB, dB, lddB, queue )
                				!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib

                	   call magmaf_queue_sync(queue) 	
			   call magmaf_dsymm('L','L',ib,jb,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )
              			!print *, "C(i:i+ib,j:j+jb)=C(",i,":",i+ib-1,",",j,":",j+jb-1,")*A(",i,":",i+ib-1,",1:k)*B(1:k,",j,":",j+jb-1,")"

                           !info = magmaf_free( dA )
        		   !if (info .ne. 0) then
            		!	print *, '2 magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            		!	stop
        		 !  endif
                          ! print *, "free dA"
                           !freeA=1
				!print *, " m-mb=",m-mb
			    if(i<m-mb)then
                                !print *, "rec    (mb,m-(i+mb)+1) =>  j=",j,"jb=",jb,"i=",i,"ib=",ib,"m-(i+mb)+1",m-(i+mb)+1
                                 ldda = ceiling(real((m-(i+mb)+1))/32)*32
                                 print *, "3 ldda = ceiling(real((m-(i+mb)+1))/32)*32=",ldda
				print *, "3 ldda*mb=",ldda,mb,ldda*mb
      		   		!info = magmaf_dmalloc( dA, ldda*mb )
      		   		!if (info .ne. 0) then
            		!		print *, "Error: magmaf_dmalloc( dA ) failed, info = ", info
            		!		goto 1000
      		   	!	endif
			!		print *, "alloco dA; ldda*mb", ldda,mb,ldda*mb," GB=",(dble(ldda)*dble(mb)*8.0d0)/1e9
                         !              totM=((dble(lddc)*dble(nb)+dble(ldda)*dble(mb)+dble(lddb)*dble(mb))*8.0D0)/1e9
       		   		print *, "3 BIG_dsymm; Memory in device=??", totM," GB"
                                    print *,"  rec hA(",i,":",I+IB-1,",",i+mb,":",i+mb+m-(i+mb)+1-1,")"
 				! rect
                                !! dA = hA
        			!call magmaf_dsetmatrix( ib, mb, hA(i,1), lda, dA, ldda, queue )
               			 call magmaf_dsetmatrix_async( m-(i+mb)+1,ib, hA(i+mb,i), lda, dA, ldda, queue )
                                    !print *,"  rec hB(",i+mb,":",I+MB+m-(i+mb)+1-1,",",j,":",J+JB-1,")"
                                			!! dB = hB
        		        !call magmaf_dsetmatrix( k, jb, hB(1,j), ldB, dB, lddB, queue )
                	        call magmaf_dsetmatrix_async( m-(i+mb)+1, jb, hB(i+mb,j), ldB, dB, lddB, queue )
                				!print *, "j=",j,"jb=",jb,"i=",i,"ib=",ib

                		call magmaf_queue_sync(queue) 	
				call magmaf_dgemm('T','N',ib,jb,m-(i+mb)+1,ALPHA,dA,ldda,dB,lddb,one,dC,lddC,queue )

                                !info = magmaf_free( dA )
        			!if (info .ne. 0) then
            		!		print *, '3 magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            		!		stop
        		!	endif
      			!	print *, "free dA"
                         !       freeA=1
                            endif
				print *, "++get C(",i,":",i+ib-1,",",j,":",j+jb-1,")"
                             				!! hC = dC
        		    call magmaf_dgetmatrix( ib, jb, dC, lddC, TEMP(1,1), MB, queue )
                	    !call magmaf_dgetmatrix_async( ib, jb, dC, lddC, TEMP(1,1), MB, queue )
                            hC(i:i+ib-1,j:j+jb-1)=temp(1:ib,1:jb)+beta*hC(i:i+ib-1,j:j+jb-1)
                                print *,""
             		ENDDO   !i
			print *,""
			print *,""
                        print *,""
          	   ENDDO  !j 

               ENDIF
	   ELSE  !side  R
          
		PRINT *, " NOT IMPLEMENT YET, "

           END IF  !side

        ENDIF  !! full A in device

    ENDIF  !!  maxmem<LimMem





    call magmaf_queue_sync(queue)
   
!! cleanup:
1000 continue
  !! Free GPU memory
     !if(freeA<1)then
                print *, "freeA=",freeA
        info = magmaf_free( dA )
        if (info .ne. 0) then
            print *, 'Magmaf_BIG_dsymm_cpu Error: magmaf_free( dA ) failed, info = ', info
            stop
        endif
    ! endif
        info = magmaf_free( dB )
        if (info .ne. 0) then
            print *, 'Magmaf_BIG_dsymm_cpu Error: magmaf_free( dB ) failed, info = ', info
            stop
        endif

        info = magmaf_free( dC )
        if (info .ne. 0) then
            print *, 'Magmaf_BIG_dsymm_cpu Error: magmaf_free( dC ) failed, info = ', info
            stop
        endif
    call magmaf_queue_destroy( queue )

     print *, "END magmaf_BIG_dsymm_cpu_1"
end subroutine magmaf_BIG_dsymm_cpu_1


end module

