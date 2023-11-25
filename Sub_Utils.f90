
!===============================================================================
! Copyright 2015-2018 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

      SUBROUTINE CALC_AVERAGE_TIME( TIME_ITER, ITERATIONS, TIME )
!     ==== Subroutine arguments =======================================
      INTEGER, intent(in)    ::     ITERATIONS
      DOUBLE PRECISION TIME_ITER( ITERATIONS )
      DOUBLE PRECISION TIME
!
!     This subroutine computes average value of time array, excluding 
!     quarter of best times and quarter of worst times
!
!     ==== Local scalars ===============================================
      DOUBLE PRECISION FULLTIME
      INTEGER          I, LOWBND
!
!     ==== Executable statements =======================================
!
!     Sort array by increase index
      !!CALL SORT_ARRAY( TIME_ITER, ITERATIONS, 'U' )
      CALL DSort(TIME_ITER, ITERATIONS)
!
!     Calculate average time, excluding upper and lower 1/5 of array
      LOWBND = ITERATIONS / 5
      FULLTIME = 0.0D+0
      DO 10 I = 1 + LOWBND, ITERATIONS - LOWBND
         FULLTIME = FULLTIME + TIME_ITER( I )
   10 CONTINUE
      TIME = FULLTIME / DBLE( ITERATIONS - 2 * LOWBND )
!
!     ==== End of CALC_AVERAGE_TIME ====================================
      END
!===============================================================================
!
!=======================================================================
!
      SUBROUTINE SORT_ARRAY( ARRAY, DIM, ORDER )
!     ==== Subroutine parameters =======================================
      INTEGER DIM
      DOUBLE PRECISION ARRAY( DIM )
      CHARACTER ORDER
!
!     This is subroutine for array sorting
!
!     ==== Local scalars ===============================================
      INTEGER I, J
      DOUBLE PRECISION TMP
      print *,"in sort array"
      print *,"dim", dim
      print *," ARRAY",  ARRAY
!
!     ==== Executable statements =======================================
      DO 20 I = 2, DIM
         TMP = ARRAY( I )
         J = I - 1
         print *," J",  J
            DO WHILE ( J.GE.1 )!.AND. ARRAY( J ).GT.TMP )     ##  bug
                print *,"start do while J",  J, "ARRAY( J )=",ARRAY( J ), "TMP=",TMP
                IF( ARRAY( J ).GT.TMP)THEN
                ARRAY( J + 1 ) = ARRAY( J )
                ENDIF
                J = J - 1
                print *,"do while J",  J
            ENDDO
         
         ARRAY( J + 1 ) = TMP
!
   20 CONTINUE
!
!     ==== End of SORT_ARRAY ===========================================
      END


SUBROUTINE DSORT (A, N)
    implicit none
    INTEGER N,  K, L
    DOUBLE PRECISION :: A(N)
    DOUBLE PRECISION :: TEMP
    DO  K = 1, N -1
        DO  L = K+1, N
            !print *, "K=",k," L=",L
            IF(A(K).GT.A(L)) THEN
                TEMP = A(K)
                A(K) = A(L)
                A(L) = TEMP
            ENDIF
        ENDDO
    ENDDO

   RETURN
END

