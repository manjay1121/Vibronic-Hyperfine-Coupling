program main
        implicit none
        integer,parameter::natoms=64
        integer,parameter::ndim=natoms*3
        real*8::Re(natoms,3),Rg(natoms,3),mass(natoms), coord_ES(ndim),coord_GS(ndim),mass3N(ndim)
        real*8::cart_ES2GS(natoms,3)
        integer,parameter::nline=mod(ndim,5)
        integer::iread,j,i,k
        integer::mode(ndim),dummy
        real*8:: vec_gs(ndim,ndim), vec_es(ndim,ndim),delta_q(ndim),displace_massWeight(ndim)
        real*8 :: freq(ndim)
        character*10::atom_name
        real*8,parameter::wavenumber2au=4.55633e-6
        real*8::HRF(ndim)
        real*8::temp
        integer,parameter::Ncols=ndim-6
        real*8 :: tempA(ndim,Ncols), tempQ(ndim,Ncols), tempM(ndim)
        integer,parameter::m=ndim
        integer,parameter::n=ndim-6
        real*8::q(m,n),r(n,n),v(m,n),u(m,n),a(m,n),lambda(Ncols)
        real*8::inverse_lamda(ncols,ncols),temp_mat(n,m),mass_mat(m,m)
        real*8::E_reorg(Ncols),T_rot(3,3)
        real*8,parameter::a02SI=5.291772d-11
        real*8,parameter:: autocmwave = 2.19474d5 ! Energy (cm^-1/au)
        real*8,parameter::aumass2SI=9.1093837d-31
        real*8,parameter::c_light=2.9998d10 !cm/s
        real*8,parameter::au2Joule=4.359744722207d-18 !Joule
        real*8,parameter::hbar=1.054571d-34  !reduced plank constant
        real*8::womfreq(ndim),d2v_ref(ndim,ndim)
        real*8::hess(ndim,ndim)
        real*8::q0(ndim)
        real*8,parameter:: m_au = 1822.88839d0    ! atomic unit mass
        


        print *,mod(ndim,5)
        print *,ndim/5
        print *,'natoms= ', natoms
        print *, nline

        open(16,file="freq.dat")
        do i=1,ndim
        read(16,*) freq(i)
        enddo
        close(16)

        !.....................................................
        !rest the excited state normal modes
!        open(12,file="ES_eigenvector.dat")
!        do iread=1,ndim/5
!        read(12,*) mode(5*iread-4:5*iread)
!        do j=1,ndim
!        read(12,*) dummy, vec_es(5*iread-4:5*iread,j)
!        enddo
!        enddo
!
!        !read rest modes
!        read(12,*) mode(ndim-mod(ndim,5)+1:ndim)
!        do j=1,ndim
!        read(12,*) dummy, vec_es(ndim-mod(ndim,5)+1:ndim,j)
!        enddo
!
!        close(12)
        !.....................................................
        !read the groud state normal modes

        open(13,file="GS_eigenvector.dat")
        do iread=1,ndim/5
        read(13,*) mode(5*iread-4:5*iread)
        print *, mode(5*iread-4:5*iread)
        do j=1,ndim
        read(13,*) dummy, vec_gs(5*iread-4:5*iread,j)
        enddo
        enddo

        if( mod(ndim,5) .ne.0 ) then
          !read rest modes
          read(13,*) mode(ndim-mod(ndim,5)+1:ndim)
          do j=1,ndim
          read(13,*) dummy, vec_gs(ndim-mod(ndim,5)+1:ndim,j)
          enddo
        endif
        close(13)


        open(12,file="GS_geom.xyz")
        do i=1,natoms
        read(12,*) atom_name, mass(i),Rg(i,:) !coord_GS(3*i-2:3*i)
        enddo
        close(12)

        open(12,file="ES_geom.xyz")
        do i=1,natoms
        read(12,*) atom_name, mass(i),Re(i,:) !coord_ES(3*i-2:3*i)
        enddo
        close(12)


        do i=1,natoms
        mass3N(3*i-2:3*i)=mass(i)*m_au
        enddo

        do i=1,natoms
        coord_GS(3*i-2:3*i)=Rg(i,:)
        coord_ES(3*i-2:3*i)=Re(i,:)
        enddo


        !do i=1,ndim
        !print *,sqrt(d2v_ref(i,i))*autocmwave
        !enddo


        call align_Kabsch_sep(natoms,mass3N,Re,Rg,cart_ES2GS,T_rot)


        do i=1,natoms
        coord_GS(3*i-2:3*i)=Rg(i,:)
        coord_ES(3*i-2:3*i)=cart_ES2GS(i,:)
        enddo

       ! print *, coord_GS
       ! print *,'================================'
       ! print *,coord_ES



      !  !Gram-Schimit orthonormalization matrix
      !  !============================================

     !   ! Initialize Q and V
     !   do i=1,m
     !   do j=1,n
     !   a(i,j)=vec_gs(j+6,i)
     !   enddo
     !   enddo
     !   q = 0.0
     !   v = a

     !   ! Compute R and Q using Gram-Schmidt orthonormalization
     !   do k = 1, n
     !   r(k,k) = norm2(v(:,k))
     !   q(:,k) = v(:,k) / r(k,k)
     !   do j = k+1, n
     !   r(k,j) = dot_product(q(:,k), v(:,j))
     !   v(:,j) = v(:,j) - r(k,j) * q(:,k)
     !   end do
     !   end do


     !   do i=1,m
     !   do j=1,n
     !   vec_gs(j+6,i)=q(i,j)
     !   enddo
     !   enddo



        !=============================================

        !calculate Huang-Rhys factor
        delta_q=0.0

        !mass weighted coordinate
        do i=1,ndim
        displace_massWeight(i)=(coord_ES(i)-coord_GS(i))*sqrt(mass3N(i))
        enddo

        do j=1,ndim

        write(96,*) sum(vec_gs(j,:)**2)

        !temp=0.d0
        !do i=1,ndim
        !temp=temp+(displace_massWeight(i)*vec_gs(j,i))**2
        !enddo
        !temp=sqrt(temp)

        temp=dot_product(vec_gs(j,:), displace_massWeight(:))
        delta_q(j)=temp
        !print *,delta_q(j)

        enddo! j=1,ndim


        open(24,file="HR_factor.dat")
        !do i=7,ndim
        do i=1,ndim
        HRF(i)=delta_q(i)**2*(freq(i)/autocmwave)/2
        write(24,*) freq(i)*0.00012,freq(i),delta_q(i),delta_q(i)*sqrt(freq(i)/autocmwave),HRF(i)
        !if (HRF(i)>0.00001) then
                !write(24,*) freq(i)*0.00012,freq(i),HRF(i)
        !endif
        enddo
        close(24)

        print *,'sum HR factor = ',sum(HRF(7:ndim))


        !q0=matmul(transpose(vec_gs),displace_massWeight)
        q0=matmul(vec_gs,displace_massWeight)

        q0(:)=q0(:)*(freq(:)/autocmwave)**0.5

        print *,q0**2/2.d0

        print *,'overall = ',sum(q0**2/2.)
        





contains

        ! Compute the L2-norm of a vector
        function norm2(x) result(norm)
                real*8, intent(in) :: x(:)
                real*8 :: norm
                norm = sqrt(sum(x**2))
        end function norm2

        ! Compute the dot product of two vectors
        function dot_product(x, y) result(dot)
                real*8, intent(in) :: x(:), y(:)
                real*8:: dot
                dot = sum(x * y)
        end function dot_product

endprogram main



!
subroutine align_Kabsch_sep( nps,fmass,cart_A, cart_B, cart_A2B ,T_rot)
   !     ------------------------------------------------------------
   !       cart_A * T_rot   compared with Cart_B
   !       The approximate rotating matrix to map cart_A to Cart_B
   !     -----------------------------------------------------------
   implicit none
   integer,intent(in )::nps
   real*8,intent(in):: cart_A(nps, 3), cart_B(nps, 3)
   real*8,intent(out):: cart_A2B(nps, 3)
   real*8:: T_rot(3,3)
   integer:: i, j
   real*8:: A_cen(3), B_cen(3), mass_tot
   real*8:: mod_cart_A(nps,3), mod_cart_B(nps,3)
   real*8:: cart_temp(nps,3)
   real*8:: Fmat(3,3)
   character*1:: jobUT, jobV
   integer:: lwork, info
   real*8:: Vmat(3,3), Umat_trans(3,3), sigma(3), work(5*3+3)
   real*8:: sign_mat(3,3), egn(3), mid_mat(3,3)
   real*8:: mul_egn
   integer:: jsign
   real*8:: T_rot_cp(3,3), cart_A2B_cp(nps,3), distance1, distance2
   real(8) ::fmass(nps*3)        
   !           
   A_cen = 0.0d0
   B_cen = 0.0d0
   mass_tot = 0.0d0
   do i = 1, nps
      A_cen = A_cen + cart_A(i,:)*fmass(3*i)
      B_cen = B_cen + cart_B(i,:)*fmass(3*i)
      mass_tot = mass_tot + fmass(3*i)
   end do
   A_cen = A_cen / mass_tot
   B_cen = B_cen / mass_tot
   !
   !        A_cen = cart_A(2,:)
   !        B_cen = cart_B(2,:)
   !
   do i = 1, nps
      mod_cart_A(i,:) = cart_A(i,:) - A_cen
      mod_cart_B(i,:) = cart_B(i,:) - B_cen
   end do
   !
   cart_temp = mod_cart_A
   do i = 1, nps
      cart_temp(i,:) = cart_temp(i,:) * fmass(3*i)
   end do
   cart_temp = cart_temp / mass_tot
   !
   Fmat = matmul(transpose( cart_temp ), mod_cart_B)
   !
   jobUT = 'A'
   jobV = 'A'
   lwork = 5*3+3
   !
   ! --------------------------------------------------------------------
   ! Notice in the following, we do the SVD
   ! as  Fmat = Vmat * Sigma * Umat_trans
   ! --------------------------------------------------------------------
   !
   call DGESVD( jobV, jobUT, 3, 3, Fmat, 3, Sigma, Vmat, 3, &
      !                  --------------------------------------------------- 
   Umat_trans, 3, work, lwork, info)
   !
   ! ----------------------------------------------------------------------------
   !      SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
   !     $                   WORK, LWORK, INFO )
   !*
   !*  -- LAPACK driver routine (version 3.3.1) --
   !*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
   !*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
   !*  -- April 2011                                                      --
   !*
   !*     .. Scalar Arguments ..
   !      CHARACTER          JOBU, JOBVT
   !      INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
   !*     ..
   !*     .. Array Arguments ..
   !      DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),
   !     $                   VT( LDVT, * ), WORK( * )
   !*     ..
   !*
   !*  Purpose
   !*  =======
   !*
   !*  DGESVD computes the singular value decomposition (SVD) of a real
   !*  M-by-N matrix A, optionally computing the left and/or right singular
   !*  vectors. The SVD is written
   !*
   !*       A = U * SIGMA * transpose(V)
   !*
   !*  where SIGMA is an M-by-N matrix which is zero except for its
   !*  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
   !*  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
   !*  are the singular values of A; they are real and non-negative, and
   !*  are returned in descending order.  The first min(m,n) columns of
   !*  U and V are the left and right singular vectors of A.
   !*
   !*  Note that the routine returns V**T, not V.
   !*
   !*  Arguments
   !*  =========
   !*
   !*  JOBU    (input) CHARACTER*1
   !*          Specifies options for computing all or part of the matrix U:
   !*          = 'A':  all M columns of U are returned in array U:
   !*          = 'S':  the first min(m,n) columns of U (the left singular
   !*                  vectors) are returned in the array U;
   !*          = 'O':  the first min(m,n) columns of U (the left singular
   !*                  vectors) are overwritten on the array A;
   !*          = 'N':  no columns of U (no left singular vectors) are
   !*                  computed.
   !*
   !*  JOBVT   (input) CHARACTER*1
   !*          Specifies options for computing all or part of the matrix
   !*          V**T:
   !*          = 'A':  all N rows of V**T are returned in the array VT;
   !*          = 'S':  the first min(m,n) rows of V**T (the right singular
   !*                  vectors) are returned in the array VT;
   !*          = 'O':  the first min(m,n) rows of V**T (the right singular
   !*                  vectors) are overwritten on the array A;
   !*          = 'N':  no rows of V**T (no right singular vectors) are
   !*                  computed.
   !*
   !*          JOBVT and JOBU cannot both be 'O'.
   !*
   !*  M       (input) INTEGER
   !*          The number of rows of the input matrix A.  M >= 0.
   !*
   !*  N       (input) INTEGER
   !*          The number of columns of the input matrix A.  N >= 0.
   !*
   !*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   !*          On entry, the M-by-N matrix A.
   !*          On exit,
   !*          if JOBU = 'O',  A is overwritten with the first min(m,n)
   !*                          columns of U (the left singular vectors,
   !*                          stored columnwise);
   !*          if JOBVT = 'O', A is overwritten with the first min(m,n)
   !*                          rows of V**T (the right singular vectors,
   !*                          stored rowwise);
   !*          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
   !*                          are destroyed.
   !*
   !*  LDA     (input) INTEGER
   !*          The leading dimension of the array A.  LDA >= max(1,M).
   !*
   !*  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
   !*          The singular values of A, sorted so that S(i) >= S(i+1).
   !*
   !*  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
   !*          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
   !*          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
   !*          if JOBU = 'S', U contains the first min(m,n) columns of U
   !*          (the left singular vectors, stored columnwise);
   !*          if JOBU = 'N' or 'O', U is not referenced.
   !*
   !*  LDU     (input) INTEGER
   !*          The leading dimension of the array U.  LDU >= 1; if
   !*          JOBU = 'S' or 'A', LDU >= M.
   !*
   !*  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
   !*          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
   !*          V**T;
   !*          if JOBVT = 'S', VT contains the first min(m,n) rows of
   !*          V**T (the right singular vectors, stored rowwise);
   !*          if JOBVT = 'N' or 'O', VT is not referenced.
   !*
   !*  LDVT    (input) INTEGER
   !*          The leading dimension of the array VT.  LDVT >= 1; if
   !*          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
   !*
   !*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   !*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
   !*          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
   !*          superdiagonal elements of an upper bidiagonal matrix B
   !*          whose diagonal is in S (not necessarily sorted). B
   !*          satisfies A = U * B * VT, so it has the same singular values
   !*          as A, and singular vectors related by U and VT.
   !*
   !*  LWORK   (input) INTEGER
   !*          The dimension of the array WORK.
   !*          LWORK >= MAX(1,5*MIN(M,N)) for the paths (see comments inside code):
   !*             - PATH 1  (M much larger than N, JOBU='N') 
   !*             - PATH 1t (N much larger than M, JOBVT='N')
   !*          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)) for the other paths
   !*          For good performance, LWORK should generally be larger.
   !*
   !*          If LWORK = -1, then a workspace query is assumed; the routine
   !*          only calculates the optimal size of the WORK array, returns
   !*          this value as the first entry of the WORK array, and no error
   !*          message related to LWORK is issued by XERBLA.
   !*
   !*  INFO    (output) INTEGER
   !*          = 0:  successful exit.
   !*          < 0:  if INFO = -i, the i-th argument had an illegal value.
   !*          > 0:  if DBDSQR did not converge, INFO specifies how many
   !*                superdiagonals of an intermediate bidiagonal form B
   !*                did not converge to zero. See the description of WORK
   !*                above for details.
   !*
   !*  -------------------------------------------------------------------
   !
   sign_mat = matmul( transpose(Umat_trans), transpose(Vmat) )
   call mat_diagonal( sign_mat, egn, 1)
   mul_egn = 1.0d0
   do i = 1, 3
      mul_egn = mul_egn * egn(i)
   end do
   if( mul_egn .lt. 0.0d0 )then
      jsign = 1
   else
      jsign = -1
   endif
   mid_mat = 0.0d0
   mid_mat(1,1) = 1.0d0
   mid_mat(2,2) = 1.0d0
   mid_mat(3,3) = dble(jsign)
   !
   T_rot = matmul(transpose( Umat_trans ), mid_mat )
   T_rot = matmul(T_rot, transpose(Vmat) )
   !
   cart_A2B = matmul( mod_cart_A, transpose(T_rot) )
   distance1 = 0.0d0
   do i = 1, nps
      distance1 = distance1 + ( cart_A2B(i,1) - mod_cart_B(i,1) )**2  &
         +  ( cart_A2B(i,2) - mod_cart_B(i,2) )**2  &
         +  ( cart_A2B(i,3) - mod_cart_B(i,3) )**2
   enddo
   !
   cart_A2B_cp = cart_A2B
   T_rot_cp = T_rot
   !
   !       -----------------------------------------------------------
   !
   if( mul_egn .lt. 0.0d0 )then
      jsign = -1
   else
      jsign = 1
   endif
   mid_mat = 0.0d0
   mid_mat(1,1) = 1.0d0
   mid_mat(2,2) = 1.0d0
   mid_mat(3,3) = dble(jsign)
   !
   T_rot = matmul(transpose( Umat_trans ), mid_mat )
   T_rot = matmul(T_rot, transpose(Vmat) )
   !
   cart_A2B = matmul( mod_cart_A, transpose(T_rot) )
   distance2 = 0.0d0
   do i = 1, nps
      distance2 = distance2 + ( cart_A2B(i,1) - mod_cart_B(i,1) )**2  &
         +  ( cart_A2B(i,2) - mod_cart_B(i,2) )**2  &
         +  ( cart_A2B(i,3) - mod_cart_B(i,3) )**2
   enddo
   !
   !       ------------------------------------------------------------
   !
   !        write(*,*) distance1, distance2
   !        write(*,*) cart_A2B
   !        write(*,*) cart_A2B_cp
   !        write(*,*) '-----------------------'
   !
   if (distance2 .gt. distance1)then
      cart_A2B = cart_A2B_cp
      T_rot = T_rot_cp
   end if
   !print *,'distance 2',distance2
   !
   do i = 1, nps
      cart_A2B(i,:) = cart_A2B(i,:) + B_cen
   end do
   !
   return
end subroutine align_Kabsch_sep
!


!....................................................
!@use lapack to diagnal matrix
!with lapack subrotine dsyev
!@So make sure you have already install lapack in 
!your computer, or you can add -mkl when you compile 
!your code when you use ifort compiler
!...................................................
subroutine mat_diagonal(d2v,womfreq,ndim)
   implicit none
   integer,intent(in):: ndim
   real*8,intent(inout):: d2v(3*ndim,3*ndim)
   real*8,intent(out):: womfreq(3*ndim)
   character*1::jobz, uplo
   real*8:: work(9*ndim-1)
   integer:: lwork, info

   jobz = 'V'
   uplo = 'U'
   lwork = 9*ndim-1
   call dsyev(jobz,uplo,3*ndim,d2v,3*ndim,womfreq,work,lwork,info)

   return
end subroutine mat_diagonal

