      Module fx_EigJacob
      
      implicit none
     
      
      contains
      
      
C =============================================================
C MS FORTRAN Program, file EIGJACOB.FOR
C
C Subroutine Jacobi, finding eigenvalues and eigenvectors
C for symmetric matrix by iteration known as the Jacobi method.
C -------------------------------------------------------------

$undefine Debug_Print

C -------------------------------------------------------------
C Litt. ref:
C   Lee W.Johnson & R.Dean Riess: NUMERICAL ANALYSIS
C            Addison Wesley 1982, ISBN 0.201-10392-3
C   Atkinson and Harley : An Introduction to numerical Methods
C            with Pascal, Addison-Wesley 1983.
C
C The code below is based on a subroutine found in
C Turbo Pascal Numerical Toolbox,
C but slightly modified regarding export of Eigenvalues and
C ordering of Eigenvector matrix.
C
C                     NIVA        B. Bjerkeng        16.10.1991
C Printed: %DATE
C =============================================================

$if defined  Debug_print
C ----------------- subroutine for printout of test results:
      subroutine matprint(Unit, Name, RDim, Rows, Columns, Matrix )

      integer Unit
      character*(*) Name
      integer RDim, Rows, Columns
      real*8  Matrix (RDim, Columns)
      
      
      integer Row, Col

      write(Unit,*)' ------- ', Name,':'
      do Row = 1, Rows
          write (Unit, '(1x,6G12.5)' ) (Matrix(Row,Col),Col=1,Columns)
      enddo
      end
$endif


C ====================================================================
      subroutine Jacobi( Eigvec_initiated, Dimen, MDim, MaxIter,
     &         Tolerance, Mat, Eigenvectors, W, Iter, Error, Accuracy )
      

C In:
      logical Eigvec_initiated  ! Signals if Eigenvectors carries
C                                 guess of values as input.
      integer Dimen             ! Dimensioning size of matrix
      integer MDim              ! Active (used) size of matrix <=Dimen
      integer MaxIter           ! Max. number of iterations
      real*8  Tolerance         ! Error tolerance

C           The error tolerance is acceptable absolute value of largest
C           off_diagonal element relative to largest diagonal value.
C           In any case the iteration stops when the accuracy stops
C           improving during 3 subsequent steps.

C In/Out:
      real*8  Mat (Dimen,Dimen)
C           In:   Symmetrix matrix
C           Out:  Matrix transformed by rotations to diagonal matrix,
C                 with Diagonal values containing Eigenvalue estimates.

C Out:
      real*8  W(Dimen,Dimen) !Work matrix, only used if Eigvec_initiated
      real*8  Eigenvectors(Dimen,Dimen)
C                  (I,N),I=1,Dimen is Eigenvector for Eigenvalue nr. N
      integer Iter
      integer Error         ! Error status
      real*8  Accuracy      ! Achieved accuracy > Tolerance if error=5

C -------------------------------------------------------------
C The subroutine uses the Jacobi eigenvalue iteration method
C to find Eigenvalues(k), k=1..Dimen
C with corresponding ( Eigenvectors(i,k),k=1..Dimen ), defined
C by the equation:
C   Matrix(i,j)*Eigenvector(j,k) = Eigenvector(i,k)*Eigenvalue(k)
C -------------------------------------------------------------


      integer Row, Col, Diag, Index, I, K
      real*8  X,Y
      real*8  SinTheta, CosTheta
      real*8  LargestDiag, ZeroDeviation, Prev_accuracy(2)
      logical Done, Rotated
      integer M
C --------------------------------------------------------------

$if defined  PRINT_IN_SUBROUTINE
      call matprint(999, 'Matrix Entered:'  , Dimen, MDim, MDim,  Mat )
$endif

      Error = 0
      M = Min(Dimen, MDim)

C -------- Checks symmetry:
      do Row = 1, M, - 1
         do Col = Row+1, M
             X =  Mat(Row,Col) - Mat(Col,Row)
             Y =  ( Mat(Row,Col) + Mat(Col,Row) )/2.0
             if ( X/10.+Y .ne. Y ) then
               Error = 1  ! Matrix not symmetric within number precision
             else
               Mat(Row,Col) = Y
               Mat(Col,Row) = Y
             endif
         enddo
      enddo


      Iter = 0
C --------- Initiate rotation:
      if ( .not. Eigvec_initiated ) then
C     ....... Eigenvectors to the identity matrix:
         do Diag = 1, M
            Eigenvectors(Diag,Diag) = 1
            do Index = 1, Diag - 1
               Eigenvectors(Diag,Index) = 0
               Eigenvectors(Index,Diag) = 0
            enddo
         enddo
      else
C     ....... Rotate input matrix by Eigvec, and continue from there:
C             This is done by the formula:
C         RotatedMatrix(i,k) = Eigvec~(i,r)*InputMatrix(r,c)*Eigvec(c,k)
C             where Transposed Eigvec~(i,r) = Eigvec(r,i)
C             summed over s and s:
C      1. ...... collect new terms in work matrix:
         do I=1, M
            do K=1, I
               X = 0.0
               do Row = 1, M
                   do Col = 1, M
                      X = X + Eigenvectors(Row,I)
     &                        *Mat(Row,Col)*Eigenvectors(Col,K)
                   enddo
               enddo
               W(I, K) = X
            enddo
         enddo
C      2. ...... move to original matrix:
         do I = 1, M
            do K = 1, I-1
               Mat(I, K) = W(I, K)
               Mat(K, I) = W(I, K)
            enddo
            Mat(I, I) = W(I, I)
         enddo
      endif

C ------- Rotate matrix successively until solution is obtained:
      Done = .FALSE.
      do while ( .NOT. Done .and. Iter .lt. Maxiter )
         Iter = Iter+1
         LargestDiag = 0
         do Diag = 1, M
            LargestDiag = Max ( LargestDiag, ABS(Mat(Diag,Diag) ) )
         enddo

         Prev_accuracy(1) = Prev_accuracy(2)
         Prev_accuracy(2) = Accuracy
         Accuracy  =  0.0

C     ..... Check upper offdiagonal elements and rotate where necessary:
         Rotated = .false.
         do Row = 1, M - 1
            do Col = Row + 1, M
               ZeroDeviation = ABS( Mat(Row,Col)/LargestDiag  )
               Accuracy  =  Max( Accuracy, ZeroDeviation )
               if (  ZeroDeviation .gt. ABS(Tolerance) ) THEN
                  call CalcRotation( Mat(Row, Row), Mat(Row, Col),
     &                         Mat(Col, Col), SinTheta, CosTheta )
                  call RotateMatrix( Dimen, M, SinTheta, CosTheta,
     &                          Row, Col, Mat )
                  call RotateEigenvectors( Dimen, M, SinTheta, CosTheta,
     &                                Row, Col, Eigenvectors )
                  Rotated = .true.
               endif
            enddo
         enddo
         if ( Iter .gt.2. and. Accuracy .ge. Prev_accuracy(1)) THEN
              Done = .true.
         else
              Done = .not. Rotated
         endif
      enddo

$if defined  PRINT_IN_SUBROUTINE
      call matprint(999, 'Matrix Rotated:'  , Dimen, MDim, MDim,  Mat )
      call matprint(999, 'Eigenvectors in columns:'   ,
     &                            Dimen, MDim, MDim,  Eigenvectors)
$endif

C Eigenvalues are exported as diagonal values of rotated matrix Mat

      if ( Iter .gt. MaxIter ) then
          Error = 2
      endif
      end subroutine



C ==============================================================
      Subroutine CalcRotation ( RowRow, RowCol, ColCol,
     &                               SinTheta, CosTheta )
      

C In:
      real*8 RowRow, RowCol, ColCol

C Out:
      real*8 SinTheta, CosTheta

C This procedure calculates the sine and cosine of the
C angle Theta through which to rotate the matrix Mat.
C Given the tangent of 2*Theta, the tangent of Theta can
C be calculated with the quadratic formula.  The cosine
C and sine are easily calculable from the tangent. The
C rotation must be such that the Row, Column element is
C zero. RowRow is the Row,Row element RowCol is the
C Row,Column element ColCol is the Column,Column element
C of Mat.
C---------------------------------------------------------
      real*8 DiagSum
      real*8 TangentTwoTheta, TangentTheta, Dummy

C---------------------------------------------------------

      DiagSum = ABS(RowRow ) + ABS(ColCol)
      if ( ABS (RowRow - ColCol)/10. + Diagsum .ne. 0.0 )  then
C     ...... significantly different diagonal elements:
         TangentTwoTheta = (RowRow - ColCol) / (2 * RowCol)
         Dummy = Sqrt( TangentTwoTheta**2 + 1.0D0)
         if (TangentTwoTheta .lt. 0. ) then ! Choose root nearer to zero
            TangentTheta = -TangentTwoTheta - Dummy
         else
            TangentTheta = -TangentTwoTheta + Dummy
         endif
         CosTheta = 1.D0 / Sqrt(1.D0 + (TangentTheta)**2)
         SinTheta = CosTheta * TangentTheta
      else
C     ...... almost equal diagonal elements, assume them identical:
         CosTheta = Sqrt(0.5D0)
         if ( RowCol .lt. 0 ) then
            SinTheta = - CosTheta
         else
            SinTheta =   CosTheta
         endif
      endif
      end subroutine


C ============================================================
      Subroutine RotateMatrix( Dimen, M, SinTheta, CosTheta,
     &                         Row, Col, Mat )


C In:
      integer Dimen, M
      real*8  SinTheta
      real*8  CosTheta
      integer Row
      integer Col

C Out:
      real*8 Mat ( Dimen, Dimen)

C------------------------------------------------------------
C This procedure rotates the matrix Mat through an angle
C Theta.  The rotation matrix is the identity matrix execept
C for the Row,Row Row,Col Col,Col and Col,Row elements.
C The rotation will make the Row,Col element of Mat
C to be zero.
C------------------------------------------------------------

      real*8 CosSqr, SinSqr, SinCos
      real*8 MatRowRow, MatColCol, MatRowCol, MatRowIndex, MatColIndex
      Integer Index

C------------------------------------------------------------

      CosSqr = CosTheta**2
      SinSqr = SinTheta**2
      SinCos = SinTheta * CosTheta

      MatRowRow =   Mat(Row, Row) * CosSqr + 2 * Mat(Row, Col) * SinCos
     &            + Mat(Col, Col) * SinSqr
      MatColCol =   Mat(Row, Row) * SinSqr - 2 * Mat(Row, Col) * SinCos
     &            + Mat(Col, Col) * CosSqr
      MatRowCol = ( Mat(Col, Col) - Mat(Row, Row) ) * SinCos
     &            + Mat(Row, Col) * (CosSqr - SinSqr)

      do Index = 1, M
         if ( Index.ne.Row .and. Index.ne.Col ) then
            MatRowIndex =   Mat(Row, Index) * CosTheta
     &                    + Mat(Col, Index) * SinTheta
            MatColIndex = - Mat(Row, Index) * SinTheta
     &                    + Mat(Col, Index) * CosTheta
            Mat(Row, Index) = MatRowIndex
            Mat(Index, Row) = MatRowIndex
            Mat(Col, Index) = MatColIndex
            Mat(Index, Col) = MatColIndex
        end if
      end do
      Mat(Row, Row) = MatRowRow
      Mat(Col, Col) = MatColCol
      Mat(Row, Col) = MatRowCol
      Mat(Col, Row) = MatRowCol

      end subroutine


C ============================================================
      subroutine RotateEigenvectors( Dimen, M, SinTheta, CosTheta,
     &                               Row, Col, Eigenvectors )


C In:
      integer  Dimen, M
      real*8   SinTheta
      real*8   CosTheta
      integer  Row
      integer  Col
C Out:
      real*8   Eigenvectors ( Dimen, Dimen )

C------------------------------------------------------------
C This procedure rotates the Eigenvectors matrix through an
C angle Theta.  The rotation matrix is the identity matrix
C except for the Row,Row Row,Col Col,Col and Col,Row
C elements.  The Eigenvectors matrix will be the product of
C all the rotation matrices which operate on Mat.
C------------------------------------------------------------

      real*8 EigenvectorsIndexRow
      real*8 EigenvectorsIndexCol
      integer Index

C -----------------------------------------------------------
      do Index = 1, M

C In the PASCAL Original:
c        EigenvectorsRowIndex =   CosTheta * Eigenvectors(Row, Index)
c    &                          + SinTheta * Eigenvectors(Col, Index)
c        EigenvectorsColIndex = - SinTheta * Eigenvectors(Row, Index)
c    &                          + CosTheta * Eigenvectors(Col, Index)
c         Eigenvectors(Row, Index) = EigenvectorsRowIndex
c         Eigenvectors(Col, Index) = EigenvectorsColIndex

C Now changed to:
         EigenvectorsIndexRow =   CosTheta * Eigenvectors(Index, Row)
     &                          + SinTheta * Eigenvectors(Index, Col)
         EigenvectorsIndexCol = - SinTheta * Eigenvectors(Index, Row)
     &                          + CosTheta * Eigenvectors(Index, Col)
         Eigenvectors(Index, Row) = EigenvectorsIndexRow
         Eigenvectors(Index, Col) = EigenvectorsIndexCol
C In order to get eigenvectors by column ( given by second index)
      end do
      end subroutine
      
      end Module
