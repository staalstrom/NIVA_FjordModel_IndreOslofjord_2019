      Module fx_Rand_Vec
      
      implicit none      

      ! F95 module in fixed format


C ******** Modified by Birger Bjerkeng, NIVA:
C 1. moved position of statement ' data TEST ...' to get 
C    compilation without errors.
C 2. Rewritten subroutines to allow 
C    parallell propagation of
C    several independent sequences
C 3. prepared separate test progrm
C Test runs OK in WATCOM FORTRAN
C ********


C This random number generator originally appeared in "Toward a Universal 
C Random Number Generator" by George Marsaglia and Arif Zaman. 
C Florida State University Report: FSU-SCRI-87-50 (1987)
C 
C It was later modified by F. James and published in "A Review of Pseudo-
C random Number Generators" 
C 
C THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
C       (However, a newly discovered technique can yield 
C         a period of 10^600. But that is still in the development stage.)
C
C It passes ALL of the tests for random number generators and has a period 
C   of 2^144, is completely portable (gives bit identical results on all 
C   machines with at least 24-bit mantissas in the floating point 
C   representation). 
C 
C The algorithm is a combination of a Fibonacci sequence (with lags of 97
C   and 33, and operation "subtraction plus one, modulo one") and an 
C   "arithmetic sequence" (using subtraction).
C
C On a Vax 11/780, this random number generator can produce a number in 
C    13 microseconds. 
C======================================================================== 

      CONTAINS
      

C ***************************************************************************      
      subroutine RMARIN( IJ, UC, NSeries)
      integer NSeries
      integer IJ(2,NSeries)
      real*8 UC(98,NSeries)

C This is the initialization routine for the random number generator RANMAR()
C This version of the subroutine intializes NSeries independent sequences,
C from specified pair of seeds in IJ.

C The initiation resets arrays IJ and UC, these arrays saves
C the current state of each sequence, and are to be used in later calls
C to RANVAL (below) for generating the random sequences.

C Note that Input values of IJ is not preserved.

C NOTE: The seed variables for each sequence consists of two values
C IJ(i,..),i=1,2, which can have values between:    
C    0 <= IJ(1,...) <= 31328
C    0 <= IJ(2,...) <= 30081

C The random number sequences created by these two seeds are of sufficient 
C length to complete an entire calculation with. For example, if several 
C different groups are working on different parts of the same calculation,
C each group could be assigned its own IJ seed. This would leave each group
C with 30000 choices for the second seed. That is to say, this random 
C number generator can create 900 million different subsequences -- with 
C each subsequence having a length of approximately 10^30.
C 
C Use IJ(1,...) = 1802 & IJ(2,...) = 9373 to test 
C the random number generator. 

C The subroutine RANMAR should be used to generate 20000 random numbers.
C Then display the next six random numbers generated multiplied by 4096*4096
C If the random number generator is working properly, the random numbers
C should be:
C           6533892.0  14220222.0  7275067.0
C           6172232.0  8354498.0   10633180.0


      logical TEST
      common /raset1/ TEST
      data TEST /.FALSE./
      integer idim, IJ1, IJ2, m, i, j, k, l, II, JJ
      real*8 S,T

      do 10 idim = 1,NSeries
         IJ1=IJ(1,idim)
         IJ2=IJ(2,idim)
         if( IJ1 .lt. 0  .or.  IJ1 .gt. 31328  .or.
     *       IJ2 .lt. 0  .or.  IJ2 .gt. 30081 ) then
             print '(A)', ' The first random number seed must have a va
     *lue between 0 and 31328'
             print '(A)',' The second seed must have a value between 0 
     *and 30081'
             stop
         endif

         i = mod(IJ1/177, 177) + 2
         j = mod(IJ1    , 177) + 2
         k = mod(IJ2/169, 178) + 1
         l = mod(IJ2,     169) 

         do 2 ii = 1, 97
            s = 0.0
            t = 0.5
            do 3 jj = 1, 24
               m = mod(mod(i*j, 179)*k, 179)
               i = j
               j = k
               k = m
               l = mod(53*l+1, 169)
               if (mod(l*m, 64) .ge. 32) then
                  s = s + t
               endif
               t = 0.5 * t
3           continue
            UC(ii,idim) = s
2        continue

         UC (98,idim) = 362436.0 / 16777216.0
         IJ(1,idim) = 97
         IJ(2,idim) = 33

10    continue

      TEST = .TRUE.
      return
      end subroutine


C ***************************************************************************      
      subroutine ranmar (IJ, UC, NSeries, RVEC, RDim, LEN)
      integer NSeries
      real*8 UC(98,NSeries)
      integer IJ(2,NSeries)
      integer RDim, Len
      real*8 RVEC(RDim,NSeries)

C This is the random number generator proposed by George Marsaglia in 
C Florida State University Report: FSU-SCRI-87-50
C It was slightly modified by F. James to produce an array of pseudorandom
C numbers.

C This version has been rearranged by B.Bjerkeng, Norwegian Institute 
C for Water Research, to allow independent and/or parallell 
C progagation of several sequences at the same time.

C The subroutine will generate LEN new values in each series, updating 
C the arrays UC and IJ so that the next call will create the 
C continuation sequence, independently of the size of LEN. 
C i.e. the call generates value nr. X(i)+1 to X(i)+LEN in each of the
C sequences i=1, NSeries, where X(i) is the sum of LEN values in 
C previous calls involving each series since initiation. 

C The array RVEC has been set up as 2-dimensional array with first 
C dimension independent of LEN, so that the series can be propagated 
C independently of each other.  For instance, the call
C      call ranmar (IJ(1,i), UC(1,i), 2, RVEC(1,i), RDim, LEN)
C will create the next LEN values of series the 2 series nr. i and i+1
C keeping the other sequences at the same state.

C It is vital that Arrays UC and IJ is kept unchanged between calls 
C to RANMAR. 
      
      
      real*8 C
      
      real*8 CD, CM
      parameter ( CD = 7654321.0 / 16777216.0,
     &            CM = 16777213.0 /16777216.0  )
      
      integer I97, J97
      logical TEST
      common /raset1/ TEST
 
      integer idim, ivec
      real*8 uni

      if( .NOT. TEST ) then
         print '(A)',' Call the init routine (RMARIN) before calling RAN        
     *MAR'  
         stop
      endif
      
      do 200 idim = 1, NSeries
         C = UC(98,idim)
         I97 = IJ(1,idim)
         J97 = IJ(2,idim)
         do 100 ivec = 1, LEN
            uni = UC(I97,idim) - UC(J97,idim)
            if( uni .lt. 0.0 ) uni = uni + 1.0
            UC(I97,idim) = uni
            I97 = I97 - 1
            if(I97 .eq. 0) I97 = 97
            J97 = J97 - 1
            if(J97 .eq. 0) J97 = 97
            C = C - CD
            if( C .lt. 0.0 ) C = C + CM
            uni = uni - C
            if( uni .lt. 0.0 ) uni = uni + 1.0
            RVEC(ivec,idim) = uni
100      continue
         UC(98,idim) = C
         IJ(1,idim) = I97 
         IJ(2,idim) = J97 
200   continue
      return
      end subroutine

      end module
