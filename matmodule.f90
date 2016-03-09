!MatModule.f90: A Fortran 2003 module to make Fortran matrix manipulation easier.
!Copyright (C) 2016 Amir Hossein Pasdar

!This program is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.

!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.

!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE matmodule
    INTEGER, PARAMETER :: dp=SELECTED_REAL_KIND(14)
    TYPE :: realptr
        REAL (KIND=dp), POINTER :: rptr
    END TYPE realptr
    ! -----------------------Eigen type------------------------------
    TYPE, PUBLIC :: eigen
    END TYPE eigen
    TYPE, PUBLIC, EXTENDS(eigen) :: eigen_r_s
        REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: eigvec
        REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: eigval
    CONTAINS
        PROCEDURE, PUBLIC :: geteigval => geteigval_r_s
        PROCEDURE, PUBLIC :: geteigvec => geteigvec_r_s
    END TYPE eigen_r_s
    TYPE, PUBLIC, EXTENDS(eigen) :: eigen_c_h
        COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: eigvec
        REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: eigval
    CONTAINS
        PROCEDURE, PUBLIC :: geteigval => geteigval_h_c
        PROCEDURE, PUBLIC :: geteigvec => geteigvec_h_c
    END TYPE eigen_c_h
    TYPE, PUBLIC, EXTENDS(eigen) :: eigen_c
        COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: leigvec, reigvec
        COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: eigval
    CONTAINS
        PROCEDURE, PUBLIC :: geteigval => geteigval_c
        PROCEDURE, PUBLIC :: getleigvec => getleigvec_c
        PROCEDURE, PUBLIC :: getreigvec => getreigvec_c
    END TYPE eigen_c
    TYPE, PUBLIC, EXTENDS(eigen) :: eigen_r
        REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: leigvec, reigvec
        REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: eigval
    CONTAINS
        PROCEDURE, PUBLIC :: geteigval => geteigval_r
        PROCEDURE, PUBLIC :: getleigvec => getleigvec_r
        PROCEDURE, PUBLIC :: getreigvec => getreigvec_r
    END TYPE eigen_r
    !---------------------------Finish eigen type---------------------------
    INTERFACE matsolve
        MODULE PROCEDURE solve_r
        MODULE PROCEDURE solve_c
    END INTERFACE matsolve
    INTERFACE matsolvex
        MODULE PROCEDURE solve_rx
        MODULE PROCEDURE solve_cx
    END INTERFACE matsolvex
    INTERFACE matwrite
        MODULE PROCEDURE matwrite_r
        MODULE PROCEDURE matwrite_c
        MODULE PROCEDURE matwrite_vc
        MODULE PROCEDURE matwrite_vr
    END INTERFACE matwrite
    INTERFACE OPERATOR (.meq.)
        MODULE PROCEDURE realmateq
        MODULE PROCEDURE complexmateq
        MODULE PROCEDURE realveceq
        MODULE PROCEDURE complexveceq
    END INTERFACE
    INTERFACE OPERATOR (.d.)
        MODULE PROCEDURE matdividerealreal
        MODULE PROCEDURE matdividecomplexcomplex
        MODULE PROCEDURE matdividerealcomplex
        MODULE PROCEDURE matdividecomplexreal
    END INTERFACE
    INTERFACE OPERATOR (.dx.)
        MODULE PROCEDURE xmatdividerealreal
        MODULE PROCEDURE xmatdividecomplexcomplex
        MODULE PROCEDURE xmatdividerealcomplex
        MODULE PROCEDURE xmatdividecomplexreal
    END INTERFACE
    INTERFACE OPERATOR (.m.)
        MODULE PROCEDURE matmultiplyrealreal
        MODULE PROCEDURE matmultiplycomplexcomplex
        MODULE PROCEDURE matmultiplyrealcomplex
        MODULE PROCEDURE matmultiplycomplexreal
        MODULE PROCEDURE matmultiplyintint
    END INTERFACE
    INTERFACE eig
        MODULE PROCEDURE eigc
        MODULE PROCEDURE eigr
    END INTERFACE eig
    INTERFACE lu
        MODULE PROCEDURE lur
        MODULE PROCEDURE luc
    END INTERFACE lu
    INTERFACE det
        MODULE PROCEDURE detr
        MODULE PROCEDURE detc
    END INTERFACE det
    INTERFACE inv
        MODULE PROCEDURE invr
        MODULE PROCEDURE invc
    END INTERFACE inv
    PRIVATE :: solve_r, solve_c, matwrite_c, matwrite_r, matwrite_vc, matwrite_vr, matdividerealreal,&
        matdividecomplexcomplex, matdividerealcomplex, matdividecomplexreal, eigc, eigr, detc, detr, invc,&
        invr, luc, lur, matmultiplycomplexcomplex, matmultiplycomplexreal, matmultiplyrealcomplex, matmultiplyrealreal,&
        matmultiplyintint, xmatmultiplycomplexcomplex, xmatmultiplycomplexreal, xmatmultiplyrealcomplex, xmatmultiplyrealreal,&
        solve_rx, solve_cx
CONTAINS
    !####################################################################################################
    SUBROUTINE matwrite_r(a,outfile)
        IMPLICIT NONE
        REAL(KIND=dp), DIMENSION (:,:), INTENT (IN) :: a
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: outfile
        INTEGER::dim1,dim2
        INTEGER :: i,j,m,STATUS
        dim1=UBOUND(a,1)
        dim2=UBOUND(a,2)
        IF (PRESENT(outfile)) THEN
            OPEN (UNIT=9,FILE=outfile,ACTION='write',IOSTAT=STATUS)
        ELSE
            OPEN (UNIT=9,FILE='MatrixOutput',ACTION='write',IOSTAT=STATUS)
        ENDIF
        ! building the outer line of the table
        WRITE(9,300,ADVANCE='no') ' '
        DO m=1,dim2
            WRITE(9,200,ADVANCE='no')
        ENDDO
        WRITE(9,300)'|'
        ! building outer line finished
        DO i=1,dim1
            ! writing the indices in a nice way
            WRITE(9,300,ADVANCE='no') ' '
            DO m=1,dim2
                WRITE(9,400,ADVANCE='no')i,m
            ENDDO
            WRITE(9,300)'|'
            ! writing the indices in a nice way finished
            ! writing the line after indices
            WRITE(9,300,ADVANCE='no') ' '
            DO m=1,dim2
                WRITE(9,500,ADVANCE='no')
            ENDDO
            WRITE(9,300)'|'
            ! writing the line after indices
            ! writing the matrix values
            DO j=1,dim2
                WRITE (9,100,ADVANCE='no') a(i,j)
            ENDDO
            WRITE(9,'(a2)')'|'
            ! writing the matrix values finished
            ! writing the outer lower line
            WRITE(9,300,ADVANCE='no') ' '
            DO m=1,dim2
                WRITE(9,200,ADVANCE='no')
            ENDDO
            WRITE(9,300)'|'
            ! writing the outer lower line finished
        ENDDO
        100 FORMAT (1x,'|',1x,E13.5E3,1x)
        200 FORMAT ('|================')
        300 FORMAT (a1)
        !400 format ('|',1x,'i=',i3,2x,'j=',i3,1x) !suitable for only up to 3 digits
        400 FORMAT ('| ','i=',i4,2x,'j=',i4,1x)  !suitable for up to 4 digits
        500 FORMAT ('|~~~~~~~~~~~~~~~~')
    END SUBROUTINE matwrite_r
    !####################################################################################################
    SUBROUTINE matwrite_c(a,outfile)
        IMPLICIT NONE
        COMPLEX(KIND=dp), DIMENSION (:,:), INTENT (IN) :: a
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: outfile
        INTEGER::dim1,dim2
        INTEGER :: i,j,m,STATUS
        dim1=UBOUND(a,1)
        dim2=UBOUND(a,2)
        IF (PRESENT(outfile)) THEN
            OPEN (UNIT=9,FILE=outfile,ACTION='write',IOSTAT=STATUS)
        ELSE
            OPEN (UNIT=9,FILE='MatrixOutput',ACTION='write',IOSTAT=STATUS)
        ENDIF
        ! building the outer line of the table
        WRITE(9,300,ADVANCE='no') ' '
        DO m=1,dim2
            WRITE(9,200,ADVANCE='no')
        ENDDO
        WRITE(9,300)'|'
        ! building outer line finished
        DO i=1,dim1
            ! writing the indices in a nice way
            WRITE(9,300,ADVANCE='no') ' '
            DO m=1,dim2
                WRITE(9,400,ADVANCE='no')i,m
            ENDDO
            WRITE(9,300)'|'
            ! writing the indices in a nice way finished
            ! writing the line after indices
            WRITE(9,300,ADVANCE='no') ' '
            DO m=1,dim2
                WRITE(9,500,ADVANCE='no')
            ENDDO
            WRITE(9,300)'|'
            ! writing the line after indices
            ! writing the matrix values
            DO j=1,dim2
                WRITE (9,100,ADVANCE='no') a(i,j)
            ENDDO
            WRITE(9,'(a2)')'|'
            ! writing the matrix values finished
            ! writing the outer lower line
            WRITE(9,300,ADVANCE='no') ' '
            DO m=1,dim2
                WRITE(9,200,ADVANCE='no')
            ENDDO
            WRITE(9,300)'|'
            ! writing the outer lower line finished
        ENDDO
        100 FORMAT (1x,'|',ss,e13.5e3,1x,sp,e13.5e3,1x,'j')
        200 FORMAT ('|==============================')
        300 FORMAT (a1)
        !400 format ('|',1x,'i=',i3,2x,'j=',i3,1x) !suitable for only up to 3 digits
        400 FORMAT ('|',6x,'i=',i4,6x,'j=',i4,6x)  !suitable for up to 4 digits
        500 FORMAT ('|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    END SUBROUTINE matwrite_c
    !####################################################################################################
    SUBROUTINE matwrite_vr(a,outfile)
        IMPLICIT NONE
        REAL(KIND=dp), DIMENSION (:), INTENT (IN) :: a
        REAL(KIND=dp), DIMENSION (1,UBOUND(a,1)) :: aa
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: outfile
        aa(1,:)=a
        IF (PRESENT(outfile)) THEN
            CALL matwrite(aa,outfile)
        ELSE
            CALL matwrite(aa)
        ENDIF
    END SUBROUTINE matwrite_vr
    !####################################################################################################
    SUBROUTINE matwrite_vc(a,outfile)
        IMPLICIT NONE
        COMPLEX(KIND=dp), DIMENSION (:), INTENT (IN) :: a
        COMPLEX(KIND=dp), DIMENSION (1,UBOUND(a,1)) :: aa
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: outfile
        aa(1,:)=a
        IF (PRESENT(outfile)) THEN
            CALL matwrite(aa,outfile)
        ELSE
            CALL matwrite(aa)
        ENDIF
    END SUBROUTINE matwrite_vc
    !####################################################################################################
    SUBROUTINE solve_r(a,b,x)
        IMPLICIT NONE
        REAL(KIND=dp), INTENT(IN), DIMENSION(:,:) :: a,b
        REAL(KIND=dp), INTENT(OUT), DIMENSION(UBOUND(a,2),UBOUND(b,2)) :: x
        INTEGER, DIMENSION (:), ALLOCATABLE :: ipiv
        REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: aa,bb
        INTEGER :: info, n, nrhs, lda, ldb
        ! The order: N, NRHS, A, LDA, IPIV, B, LDB, INFO
        lda=UBOUND(a,1)
        n=UBOUND(a,2)
        IF (n .ne. lda) PRINT*,'.D. : The matrix of coefficients is NOT square! Something will go wrong ...'
        ldb=UBOUND(b,1)
        nrhs=UBOUND(b,2)
        IF (n .ne. ldb) PRINT*,'.D. : The matrix of b does not have the same leading dimension as a!&
            & Something will definitely GO wrong ...'
        ALLOCATE (aa(lda,n),bb(ldb,nrhs),ipiv(n))
        aa=a
        bb=b
        CALL dgesv (n,nrhs,aa,lda,ipiv,bb,ldb,info)
        IF (info.ne.0) PRINT*, 'The linear system of equations was not solved successfully.'
        x=bb
        DEALLOCATE (aa,bb,ipiv)
    END SUBROUTINE solve_r
    !####################################################################################################
    SUBROUTINE solve_c(a,b,x)
        IMPLICIT NONE
        COMPLEX(KIND=dp), INTENT(IN), DIMENSION(:,:) :: a,b
        COMPLEX(KIND=dp), INTENT(OUT), DIMENSION(UBOUND(a,2),UBOUND(b,2)) :: x
        INTEGER, DIMENSION (:), ALLOCATABLE :: ipiv
        COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: aa,bb
        INTEGER :: info, n, nrhs, lda, ldb
        ! The order: N, NRHS, A, LDA, IPIV, B, LDB, INFO
        lda=UBOUND(a,1)
        n=UBOUND(a,2)
        IF (n .ne. lda) PRINT*,'.D. : The matrix of coefficients is NOT square! Something will go wrong ...'
        ldb=UBOUND(b,1)
        nrhs=UBOUND(b,2)
        IF (n .ne. ldb) PRINT*,'.D. : The matrix of b does not have the same leading dimension as a!&
            & Something will definitely GO wrong ...'
        ALLOCATE (aa(lda,n),bb(ldb,nrhs),ipiv(n))
        aa=a
        bb=b
        CALL zgesv (n,nrhs,aa,lda,ipiv,bb,ldb,info)
        IF (info.ne.0) PRINT*, 'The linear system of equations was not solved successfully.'
        x=bb
        DEALLOCATE (aa,bb,ipiv)
    END SUBROUTINE solve_c
    !####################################################################################################
    SUBROUTINE fullpvt (a,b,ap,bp,indx)
        IMPLICIT NONE
        REAL (KIND=dp), DIMENSION (:,:), INTENT (IN), TARGET :: a,b
        REAL (KIND=dp) :: r_pvt, c_pvt
        INTEGER :: m,n,l,k,i,j,rp(1),cp(1)
        TYPE (realptr), DIMENSION(UBOUND(a,1),UBOUND(a,2)), INTENT (OUT) :: ap
        TYPE (realptr), DIMENSION(UBOUND(b,1),UBOUND(b,2)), INTENT (OUT) :: bp
        TYPE (realptr), DIMENSION(UBOUND(a,1)) :: ptra_r
        TYPE (realptr), DIMENSION(UBOUND(a,2)) :: ptra_c
        TYPE (realptr), DIMENSION(UBOUND(b,2)) :: ptrb
        INTEGER, DIMENSION(UBOUND(a,1)), INTENT(OUT) :: indx
        m=UBOUND(a,1)
        n=UBOUND(a,2)
        l=UBOUND(b,1)
        k=UBOUND(b,2)
        IF (m.ne.n) PRINT*,'FullPVT: The matrix of coefficients is not square! Something will go wrong...'
        IF (l.ne.m) PRINT*,'FullPVT: a & b don''t have the same number of rows! Check it!'
        FORALL (i=1:m,j=1:n)
            ap(i,j)%rptr=>a(i,j)
        END FORALL
        FORALL (i=1:l,j=1:k)
            bp(i,j)%rptr=>b(i,j)
        END FORALL
        !ap=>a
        !bp=>b
        FORALL (i=1:m)
            indx(i)=i
        END FORALL
        DO i=1,m
            cp=MAXLOC(ABS(a(i:,i)))
            rp=MAXLOC(ABS(a(i,i:)))
            r_pvt=a(i,rp(1))
            c_pvt=a(cp(1),i)
            IF (ABS(r_pvt)>ABS(c_pvt)) THEN
                !allocate (ptra(m),stat=istat); if (istat.ne.0) print*,'Allocation of the temporary pointer failed!'
                DO j=1,m
                    ptra_r(j)%rptr=>ap(j,i)%rptr
                    ap(j,i)%rptr=>ap(j,rp(1))%rptr
                    ap(j,rp(1))%rptr=>ptra_r(j)%rptr
                ENDDO
                indx(i)=rp(1)
                indx(rp(1))=i
            ELSEIF(ABS(r_pvt)<= ABS(c_pvt) .and. i.ne.cp(1)) THEN
                DO j=1,n
                    ptra_c(j)%rptr=>ap(i,j)%rptr
                    ap(i,j)%rptr=>ap(cp(1),j)%rptr
                    ap(cp(1),j)%rptr=>ptra_c(j)%rptr
                ENDDO
                DO j=1,k
                    ptrb(j)%rptr=>bp(i,j)%rptr
                    bp(i,j)%rptr=>bp(cp(1),j)%rptr
                    bp(cp(1),j)%rptr=>ptrb(j)%rptr
                ENDDO
            ENDIF
        ENDDO
    END SUBROUTINE fullpvt
    !####################################################################################################
    FUNCTION matdividerealreal(b,a)
        IMPLICIT NONE
        REAL (KIND=dp), INTENT (IN), DIMENSION(:,:) :: a,b
        REAL (KIND=dp), DIMENSION(UBOUND(a,2),UBOUND(b,2)) :: matdividerealreal
        CALL matsolve(a,b,matdividerealreal)
    END FUNCTION matdividerealreal
    !####################################################################################################
    FUNCTION matdividecomplexcomplex(b,a)
        IMPLICIT NONE
        COMPLEX (KIND=dp), INTENT (IN), DIMENSION(:,:) :: a,b
        COMPLEX (KIND=dp), DIMENSION(UBOUND(a,2),UBOUND(b,2)) :: matdividecomplexcomplex
        CALL matsolve(a,b,matdividecomplexcomplex)
    END FUNCTION matdividecomplexcomplex
    !####################################################################################################
    FUNCTION matdividecomplexreal(b,a)
        IMPLICIT NONE
        COMPLEX (KIND=dp), INTENT (IN), DIMENSION(:,:) :: a
        REAL (KIND=dp), INTENT (IN), DIMENSION(:,:) :: b
        COMPLEX (KIND=dp), DIMENSION(UBOUND(a,2),UBOUND(b,2)) :: matdividecomplexreal
        COMPLEX (KIND=dp), DIMENSION(UBOUND(b,1),UBOUND(b,2)) :: bb
        bb=b
        CALL matsolve(a,bb,matdividecomplexreal)
    END FUNCTION matdividecomplexreal
    !####################################################################################################
    FUNCTION matdividerealcomplex(b,a)
        IMPLICIT NONE
        REAL (KIND=dp), INTENT (IN), DIMENSION(:,:) :: a
        COMPLEX (KIND=dp), INTENT (IN), DIMENSION(:,:) :: b
        COMPLEX (KIND=dp), DIMENSION(UBOUND(a,2),UBOUND(b,2)) :: matdividerealcomplex
        COMPLEX (KIND=dp), DIMENSION(UBOUND(a,1),UBOUND(a,2)) :: aa
        aa=a
        CALL matsolve(aa,b,matdividerealcomplex)
    END FUNCTION matdividerealcomplex
    !####################################################################################################
    FUNCTION xmatdividerealreal(b,a)
        IMPLICIT NONE
        REAL (KIND=dp), INTENT (IN), DIMENSION(:,:) :: a,b
        REAL (KIND=dp), DIMENSION(UBOUND(a,2),UBOUND(b,2)) :: xmatdividerealreal
        CALL matsolvex(a,b,xmatdividerealreal)
    END FUNCTION xmatdividerealreal
    !####################################################################################################
    FUNCTION xmatdividecomplexcomplex(b,a)
        IMPLICIT NONE
        COMPLEX (KIND=dp), INTENT (IN), DIMENSION(:,:) :: a,b
        COMPLEX (KIND=dp), DIMENSION(UBOUND(a,2),UBOUND(b,2)) :: xmatdividecomplexcomplex
        CALL matsolvex(a,b,xmatdividecomplexcomplex)
    END FUNCTION xmatdividecomplexcomplex
    !####################################################################################################
    FUNCTION xmatdividecomplexreal(b,a)
        IMPLICIT NONE
        COMPLEX (KIND=dp), INTENT (IN), DIMENSION(:,:) :: a
        REAL (KIND=dp), INTENT (IN), DIMENSION(:,:) :: b
        COMPLEX (KIND=dp), DIMENSION(UBOUND(a,2),UBOUND(b,2)) :: xmatdividecomplexreal
        COMPLEX (KIND=dp), DIMENSION(UBOUND(b,1),UBOUND(b,2)) :: bb
        bb=b
        CALL matsolvex(a,bb,xmatdividecomplexreal)
    END FUNCTION xmatdividecomplexreal
    !####################################################################################################
    FUNCTION xmatdividerealcomplex(b,a)
        IMPLICIT NONE
        REAL (KIND=dp), INTENT (IN), DIMENSION(:,:) :: a
        COMPLEX (KIND=dp), INTENT (IN), DIMENSION(:,:) :: b
        COMPLEX (KIND=dp), DIMENSION(UBOUND(a,2),UBOUND(b,2)) :: xmatdividerealcomplex
        COMPLEX (KIND=dp), DIMENSION(UBOUND(a,1),UBOUND(a,2)) :: aa
        aa=a
        CALL matsolvex(aa,b,xmatdividerealcomplex)
    END FUNCTION xmatdividerealcomplex
    !####################################################################################################
    FUNCTION eigc(a)
        IMPLICIT NONE
        COMPLEX (KIND=dp), INTENT(IN), DIMENSION(:,:) :: a
        CLASS (eigen), POINTER :: eigc
        CLASS (eigen_c_h), POINTER :: eigc1
        CLASS (eigen_c), POINTER :: eigc2
        CHARACTER, PARAMETER :: jobvl='V', jobvr='V'
        INTEGER :: info, lda, ldvl, ldvr, lwork, n, i, j
        REAL (KIND=dp) , DIMENSION(2*UBOUND(a,2)) :: rwork
        REAL (KIND=dp) , DIMENSION(3*UBOUND(a,2)-2) :: rworks
        COMPLEX (KIND=dp), DIMENSION(UBOUND(a,2)) :: w
        REAL (KIND=dp), DIMENSION(UBOUND(a,2)) :: wr
        COMPLEX (KIND=dp), DIMENSION(UBOUND(a,1),UBOUND(a,2)) :: vl, vr, aa
        COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: work
        LOGICAL :: sym=.true.
        lda=UBOUND(a,1)
        ldvl=lda
        ldvr=lda
        n=UBOUND(a,2)
        IF (n.ne.lda) PRINT*, 'EIG: The matrix is not square. Check your work.'
        aa=a
        DO i=1,lda
            DO j=1,n
                IF (a(i,j).ne.CONJG(a(j,i))) THEN
                    sym=.false.
                    EXIT
                ENDIF
            ENDDO
        ENDDO
        IF (sym) THEN
            !SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
            lwork=-1
            ALLOCATE(work(n))
            CALL ZHEEV( 'V', 'L', n, aa, lda, wr, work, lwork, rworks, info)
            lwork=work(1)
            DEALLOCATE(work)
            ALLOCATE(work(lwork))
            CALL ZHEEV( 'V', 'L', n, aa, lda, wr, work, lwork, rworks, info)
            IF (info.ne.0) THEN
                PRINT*,'The alghorithm of EIG did not end successfully.'
                STOP
            ENDIF
            ALLOCATE(eigc1)
            ALLOCATE (eigc1%eigval(n),eigc1%eigvec(n,n))
            eigc1%eigval=wr
            eigc1%eigvec=aa
            eigc=>eigc1
        ELSE
            !SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
            lwork=-1
            ALLOCATE(work(n))
            CALL ZGEEV (jobvl, jobvr, n, aa, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
            lwork=work(1)
            DEALLOCATE(work)
            ALLOCATE(work(lwork))
            CALL ZGEEV (jobvl, jobvr, n, aa, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
            IF (info.ne.0) THEN
                PRINT*,'The alghorithm of EIG did not end successfully.'
                STOP
            ENDIF
            ALLOCATE(eigc2)
            ALLOCATE (eigc2%eigval(n),eigc2%reigvec(n,n),eigc2%leigvec(n,n))
            eigc2%eigval=w
            eigc2%leigvec=vl
            eigc2%reigvec=vr
            eigc=>eigc2
        ENDIF
    END FUNCTION eigc
    !###################################################################################################
    FUNCTION eigr(a)
        IMPLICIT NONE
        REAL (KIND=dp), INTENT(IN), DIMENSION(:,:) :: a
        COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: ac
        COMPLEX (KIND=dp), PARAMETER :: jj=(0._dp,1._dp)
        REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: ar
        CLASS (eigen), POINTER :: eigr
        CLASS (eigen_r_s), POINTER :: eigr1
        CLASS (eigen_c), POINTER :: eigr2
        CLASS (eigen_r), POINTER :: eigr3
        CHARACTER, PARAMETER :: jobs='V', uplo='L', jobvl='V', jobvr='V'
        INTEGER :: info, lda, lwork, n, i, j, ldvl, ldvr
        REAL (KIND=dp), DIMENSION(UBOUND(a,2)) :: w, wi
        REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: work
        REAL (KIND=dp), DIMENSION(UBOUND(a,1),UBOUND(a,2)) :: vl, vr
        LOGICAL :: sym=.true.
        lda=UBOUND(a,1)
        n=UBOUND(a,2)
        ldvl=lda
        ldvr=lda
        IF (n.ne.lda) PRINT*, 'EIG: The matrix is not square. Check your work.'
        DO i=1,lda
            DO j=1,n
                IF (a(i,j).ne.a(j,i)) THEN
                    sym=.false.
                    EXIT
                ENDIF
            ENDDO
        ENDDO
        ALLOCATE (ar(lda,n))
        ar=a
        IF (sym) THEN
            !SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
            lwork=-1
            ALLOCATE(work(n))
            CALL dsyev( jobs, uplo, n, ar, lda, w, work, lwork, info)
            lwork=work(1)
            DEALLOCATE(work)
            ALLOCATE(work(lwork))
            CALL dsyev( jobs, uplo, n, ar, lda, w, work, lwork, info)
            IF (info.ne.0) THEN
                PRINT*,'The alghorithm of EIG did not end successfully.'
                STOP
            ENDIF
            ALLOCATE(eigr1)
            ALLOCATE (eigr1%eigval(n),eigr1%eigvec(n,n))
            eigr1%eigval=w
            eigr1%eigvec=ar
            eigr => eigr1
        ELSE
            lwork=-1
            ALLOCATE(work(n))
            CALL DGEEV( JOBVL, JOBVR, N, ar, LDA, W, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
            lwork=work(1)
            DEALLOCATE(work)
            ALLOCATE(work(lwork))
            CALL DGEEV( JOBVL, JOBVR, N, ar, LDA, W, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
            IF (info.ne.0) THEN
                PRINT*,'The alghorithm of EIG did not end successfully.'
                STOP
            ENDIF
            IF (MAXVAL(ABS(wi)).eq.0._dp) THEN
                ALLOCATE(eigr3)
                ALLOCATE (eigr3%eigval(n),eigr3%reigvec(n,n),eigr3%leigvec(n,n))
                eigr3%eigval=w
                eigr3%leigvec=vl
                eigr3%reigvec=vr
                eigr => eigr3
            ELSE
                ALLOCATE(eigr2)
                ALLOCATE (eigr2%eigval(n),eigr2%reigvec(n,n),eigr2%leigvec(n,n))
                eigr2%eigval=w+wi*jj
                j=1
                DO WHILE (j.le.n)
                    IF (wi(j).eq.0._dp) THEN
                        eigr2%leigvec(:,j)=vl(:,j)
                        eigr2%reigvec(:,j)=vr(:,j)
                        j=j+1
                    ELSE
                        eigr2%leigvec(:,j)=vl(:,j)+jj*vl(:,j+1)
                        eigr2%reigvec(:,j)=vr(:,j)+jj*vr(:,j+1)
                        eigr2%leigvec(:,j+1)=vl(:,j)-jj*vl(:,j+1)
                        eigr2%reigvec(:,j+1)=vr(:,j)-jj*vr(:,j+1)
                        j=j+2
                    ENDIF
                ENDDO
                eigr => eigr2
            ENDIF
        ENDIF
    END FUNCTION eigr
    !###################################################################################################
    FUNCTION eye(n)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: n
        REAL(KIND=dp),DIMENSION(n,n) :: eye
        INTEGER :: i
        eye=0._dp
        FORALL (i=1:n)
            eye(i,i)=1._dp
        END FORALL
    END FUNCTION eye
    !###################################################################################################
    FUNCTION lur(a)
        IMPLICIT NONE
        REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: a
        REAL(KIND=dp), DIMENSION(UBOUND(a,1),UBOUND(a,2)) :: lur
        INTEGER :: info, m, n, lda
        INTEGER, DIMENSION(MIN(UBOUND(a,1),UBOUND(a,2))) :: ipiv
        m=UBOUND(a,1)
        n=UBOUND(a,2)
        lda=m
        !SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
        lur=a
        CALL dgetrf(m,n,lur,lda,ipiv,info)
        IF (info.ne.0) PRINT*, 'The algorithm of LU did not finish successfully.'
    END FUNCTION lur
    !###################################################################################################
    FUNCTION luc(a)
        IMPLICIT NONE
        COMPLEX(KIND=dp), DIMENSION(:,:), INTENT(IN) :: a
        COMPLEX(KIND=dp), DIMENSION(UBOUND(a,1),UBOUND(a,2)) :: luc
        INTEGER :: info, m, n, lda
        INTEGER, DIMENSION(MIN(UBOUND(a,1),UBOUND(a,2))) :: ipiv
        m=UBOUND(a,1)
        n=UBOUND(a,2)
        lda=m
        !SUBROUTINE zGETRF( M, N, A, LDA, IPIV, INFO )
        luc=a
        CALL zgetrf(m,n,luc,lda,ipiv,info)
        IF (info.ne.0) PRINT*, 'The algorithm  of LU did not finish successfully.'
    END FUNCTION luc
    !###################################################################################################
    FUNCTION detr(a)
        IMPLICIT NONE
        REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: a
        REAL(KIND=dp), DIMENSION(UBOUND(a,1),UBOUND(a,2)) :: lur
        REAL(KIND=dp) :: detr
        INTEGER :: i, m, n
        m=UBOUND(a,1)
        n=UBOUND(a,2)
        IF (m.ne.n) PRINT*,'For computing the determinent, the matrix must be square.'
        lur=lu(a)
        detr=1._dp
        DO i=1,m
            detr=detr*lur(i,i)
        ENDDO
    END FUNCTION detr
    !###################################################################################################
    FUNCTION detc(a)
        IMPLICIT NONE
        COMPLEX(KIND=dp), DIMENSION(:,:), INTENT(IN) :: a
        COMPLEX(KIND=dp), DIMENSION(UBOUND(a,1),UBOUND(a,2)) :: luc
        COMPLEX(KIND=dp) :: detc
        INTEGER :: i, m, n
        m=UBOUND(a,1)
        n=UBOUND(a,2)
        IF (m.ne.n) PRINT*,'For computing the determinent, the matrix must be square.'
        luc=lu(a)
        detc=1._dp
        DO i=1,m
            detc=detc*luc(i,i)
        ENDDO
    END FUNCTION detc
    !###################################################################################################
    FUNCTION invr(a)
        IMPLICIT NONE
        REAL(KIND=dp), DIMENSION(:,:), INTENT(IN) :: a
        REAL(KIND=dp), DIMENSION(UBOUND(a,1),UBOUND(a,2)) :: invr
        REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: work
        INTEGER :: info, m, n, lda, lwork
        INTEGER, DIMENSION(UBOUND(a,1)) :: ipiv
        m=UBOUND(a,1)
        n=UBOUND(a,2)
        IF (m.ne.n) PRINT*,'For computing the inverse, the matrix must be square.'
        lda=m
        invr=a
        CALL dgetrf(m,n,invr,lda,ipiv,info)
        IF (info.ne.0) PRINT*, 'The first phase of inversion did not end successfully.'
        !SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
        lwork=n !For each cpu lwork must be calculated.
        ALLOCATE(work(lwork))
        lwork=-1
        CALL dgetri(n,invr,lda,ipiv,work,lwork,info)
        lwork=work(1)
        DEALLOCATE(work)
        ALLOCATE(work(lwork))
        CALL dgetri(n,invr,lda,ipiv,work,lwork,info)
        IF (info.ne.0) PRINT*, 'The algorithm of inversion did not end successfully.'
    END FUNCTION invr
    !###################################################################################################
    FUNCTION invc(a)
        IMPLICIT NONE
        COMPLEX(KIND=dp), DIMENSION(:,:), INTENT(IN) :: a
        COMPLEX(KIND=dp), DIMENSION(UBOUND(a,1),UBOUND(a,2)) :: invc
        COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE :: work
        INTEGER :: info, m, n, lda, lwork, istat
        INTEGER, DIMENSION(UBOUND(a,1)) :: ipiv
        m=UBOUND(a,1)
        n=UBOUND(a,2)
        IF (m.ne.n) PRINT*,'For computing the inverse, the matrix must be square.'
        lda=m
        invc=a
        CALL zgetrf(m,n,invc,lda,ipiv,info)
        IF (info.ne.0) PRINT*, 'The first phase of inversion did not end successfully.'
        !SUBROUTINE zGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
        lwork=n !For each cpu lwork must be calculated.
        ALLOCATE(work(lwork))
        lwork=-1
        CALL zgetri(n,invc,lda,ipiv,work,lwork,info)
        lwork=work(1)
        DEALLOCATE(work)
        ALLOCATE(work(lwork))
        CALL zgetri(n,invc,lda,ipiv,work,lwork,info)
        IF (info.ne.0) PRINT*, 'The algorithm of inversion did not end successfully.'
    END FUNCTION invc
    !###################################################################################################
    FUNCTION geteigval_c(this)
        CLASS (eigen_c), INTENT(IN) :: this
        COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: geteigval_c
        INTEGER :: n
        n=UBOUND(this%eigval,1)
        ALLOCATE(geteigval_c(n))
        geteigval_c=this%eigval
    END FUNCTION geteigval_c
    !###################################################################################################
    FUNCTION getleigvec_c(this)
        CLASS (eigen_c), INTENT(IN) :: this
        COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: getleigvec_c
        INTEGER :: n, m
        m=UBOUND(this%leigvec,1)
        n=UBOUND(this%leigvec,2)
        ALLOCATE(getleigvec_c(m,n))
        getleigvec_c=this%leigvec
    END FUNCTION getleigvec_c
    !###################################################################################################
    FUNCTION getreigvec_c(this)
        CLASS (eigen_c), INTENT(IN) :: this
        COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: getreigvec_c
        INTEGER :: n, m
        m=UBOUND(this%reigvec,1)
        n=UBOUND(this%reigvec,2)
        ALLOCATE(getreigvec_c(m,n))
        getreigvec_c=this%reigvec
    END FUNCTION getreigvec_c
    !###################################################################################################
    FUNCTION geteigval_r(this)
        CLASS (eigen_r), INTENT(IN) :: this
        REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: geteigval_r
        INTEGER :: n
        n=UBOUND(this%eigval,1)
        ALLOCATE(geteigval_r(n))
        geteigval_r=this%eigval
    END FUNCTION geteigval_r
    !###################################################################################################
    FUNCTION getleigvec_r(this)
        CLASS (eigen_r), INTENT(IN) :: this
        REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: getleigvec_r
        INTEGER :: n, m
        m=UBOUND(this%leigvec,1)
        n=UBOUND(this%leigvec,2)
        ALLOCATE(getleigvec_r(m,n))
        getleigvec_r=this%leigvec
    END FUNCTION getleigvec_r
    !###################################################################################################
    FUNCTION getreigvec_r(this)
        CLASS (eigen_r), INTENT(IN) :: this
        COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: getreigvec_r
        INTEGER :: n, m
        m=UBOUND(this%reigvec,1)
        n=UBOUND(this%reigvec,2)
        ALLOCATE(getreigvec_r(m,n))
        getreigvec_r=this%reigvec
    END FUNCTION getreigvec_r
    !###################################################################################################
    FUNCTION geteigval_r_s(this)
        CLASS (eigen_r_s), INTENT(IN) :: this
        REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: geteigval_r_s
        INTEGER :: n
        n=UBOUND(this%eigval,1)
        ALLOCATE(geteigval_r_s(n))
        geteigval_r_s=this%eigval
    END FUNCTION geteigval_r_s
    !###################################################################################################
    FUNCTION geteigvec_r_s(this)
        CLASS (eigen_r_s), INTENT(IN) :: this
        REAL (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: geteigvec_r_s
        INTEGER :: n, m
        m=UBOUND(this%eigvec,1)
        n=UBOUND(this%eigvec,2)
        ALLOCATE(geteigvec_r_s(m,n))
        geteigvec_r_s=this%eigvec
    END FUNCTION geteigvec_r_s
    !###################################################################################################
    FUNCTION geteigval_h_c(this)
        CLASS (eigen_c_h), INTENT(IN) :: this
        REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: geteigval_h_c
        INTEGER :: n
        n=UBOUND(this%eigval,1)
        ALLOCATE(geteigval_h_c(n))
        geteigval_h_c=this%eigval
    END FUNCTION geteigval_h_c
    !###################################################################################################
    FUNCTION geteigvec_h_c(this)
        CLASS (eigen_c_h), INTENT(IN) :: this
        COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: geteigvec_h_c
        INTEGER :: n, m
        m=UBOUND(this%eigvec,1)
        n=UBOUND(this%eigvec,2)
        ALLOCATE(geteigvec_h_c(m,n))
        geteigvec_h_c=this%eigvec
    END FUNCTION geteigvec_h_c
    !####################################################################################################
    FUNCTION matmultiplyrealreal(b,a)
        IMPLICIT NONE
        REAL (KIND=dp), INTENT (IN), DIMENSION(:,:) :: a,b
        REAL (KIND=dp), DIMENSION(UBOUND(b,1),UBOUND(a,2)) :: matmultiplyrealreal
        matmultiplyrealreal=MATMUL(b,a)
    END FUNCTION matmultiplyrealreal
    !####################################################################################################
    FUNCTION matmultiplycomplexcomplex(b,a)
        IMPLICIT NONE
        COMPLEX (KIND=dp), INTENT (IN), DIMENSION(:,:) :: a,b
        COMPLEX (KIND=dp), DIMENSION(UBOUND(b,1),UBOUND(a,2)) :: matmultiplycomplexcomplex
        matmultiplycomplexcomplex=MATMUL(b,a)
    END FUNCTION matmultiplycomplexcomplex
    !####################################################################################################
    FUNCTION matmultiplycomplexreal(b,a)
        IMPLICIT NONE
        COMPLEX (KIND=dp), INTENT (IN), DIMENSION(:,:) :: a
        REAL (KIND=dp), INTENT (IN), DIMENSION(:,:) :: b
        COMPLEX (KIND=dp), DIMENSION(UBOUND(b,1),UBOUND(a,2)) :: matmultiplycomplexreal
        matmultiplycomplexreal=MATMUL(b,a)
    END FUNCTION matmultiplycomplexreal
    !####################################################################################################
    FUNCTION matmultiplyrealcomplex(b,a)
        IMPLICIT NONE
        REAL (KIND=dp), INTENT (IN), DIMENSION(:,:) :: a
        COMPLEX (KIND=dp), INTENT (IN), DIMENSION(:,:) :: b
        COMPLEX (KIND=dp), DIMENSION(UBOUND(b,1),UBOUND(a,2)) :: matmultiplyrealcomplex
        matmultiplyrealcomplex=MATMUL(b,a)
    END FUNCTION matmultiplyrealcomplex
    !####################################################################################################
    FUNCTION matmultiplyintint(b,a)
        IMPLICIT NONE
        INTEGER (KIND=SELECTED_INT_KIND(15)), INTENT (IN), DIMENSION(:,:) :: b,a
        INTEGER (KIND=SELECTED_INT_KIND(15)), DIMENSION(UBOUND(b,1),UBOUND(a,2)) :: matmultiplyintint
        matmultiplyintint=MATMUL(b,a)
    END FUNCTION matmultiplyintint
    !###################################################################################################
    SUBROUTINE solve_rx(a,b,x)
        IMPLICIT NONE
        REAL(KIND=dp), INTENT(IN), DIMENSION(:,:) :: a,b
        REAL(KIND=dp), INTENT(OUT), DIMENSION(UBOUND(a,2),UBOUND(b,2)) :: x
        INTEGER, DIMENSION (:), ALLOCATABLE :: ipiv, iwork
        REAL(KIND=dp), DIMENSION (:), ALLOCATABLE :: r, c, ferr, berr, work
        REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: bb, aa, af
        INTEGER :: info, lda, ldaf, ldb, ldx, n, nrhs
        CHARACTER  ::   equed='N', fact='N', trans='N'
        REAL(KIND=dp) :: rcond
        lda=UBOUND(a,1)
        n=UBOUND(a,2)
        ldaf=n
        IF (n .ne. lda) PRINT*,'.DX. : The matrix of coefficients is NOT square! Something will go wrong ...'
        ldb=UBOUND(b,1)
        nrhs=UBOUND(b,2)
        ldx=n
        IF (n .ne. ldb) PRINT*,'.DX. : The matrix of b does not have the same leading dimension as a!&
            & Something will definitely GO wrong ...'
        ALLOCATE (aa(lda,n), af(lda,n), bb(ldb,nrhs), ipiv(n), iwork(n), r(n), c(n),&
            ferr(nrhs), berr(nrhs), work(4*n))
        aa=a
        bb=b
        !SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R, C, B, LDB, X,&
            !& LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )
        CALL dgesvx (fact, trans, n, nrhs, aa, lda, af, ldaf, ipiv, equed, r, c, bb, ldb, x,&
            ldx, rcond, ferr, berr, work, iwork, info)
        PRINT*,'Condition number of this linear system is',rcond
    END SUBROUTINE solve_rx
    !###################################################################################################
    SUBROUTINE solve_cx(a,b,x)
        IMPLICIT NONE
        COMPLEX (KIND=dp), INTENT(IN), DIMENSION(:,:) :: a,b
        COMPLEX (KIND=dp), INTENT(OUT), DIMENSION(UBOUND(a,2),UBOUND(b,2)) :: x
        INTEGER, DIMENSION (:), ALLOCATABLE :: ipiv
        REAL(KIND=dp), DIMENSION (:), ALLOCATABLE :: r, c, ferr, berr, rwork
        COMPLEX (KIND=dp), DIMENSION(:,:), ALLOCATABLE :: bb, aa, af
        COMPLEX (KIND=dp), DIMENSION (:), ALLOCATABLE :: work
        INTEGER :: info, lda, ldaf, ldb, ldx, n, nrhs
        CHARACTER  ::   equed='N', fact='N', trans='N'
        REAL(KIND=dp) :: rcond
        lda=UBOUND(a,1)
        n=UBOUND(a,2)
        ldaf=n
        IF (n .ne. lda) PRINT*,'.DX. : The matrix of coefficients is NOT square! Something will go wrong ...'
        ldb=UBOUND(b,1)
        nrhs=UBOUND(b,2)
        ldx=n
        IF (n .ne. ldb) PRINT*,'.DX. : The matrix of b does not have the same leading dimension as a!&
            & Something will definitely GO wrong ...'
        ALLOCATE (aa(lda,n), af(lda,n), bb(ldb,nrhs), ipiv(n), rwork(2*n), r(n), c(n),&
            ferr(nrhs), berr(nrhs), work(2*n))
        aa=a
        bb=b
        !SUBROUTINE ZGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R, C, B, LDB, X,&
            !& LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
        CALL zgesvx (fact, trans, n, nrhs, aa, lda, af, ldaf, ipiv, equed, r, c, bb, ldb, x,&
            ldx, rcond, ferr, berr, work, rwork, info)
        PRINT*,'Condition number of this linear system is',rcond
    END SUBROUTINE solve_cx
    !###################################################################################################
    LOGICAL FUNCTION realmateq(a,b)
        IMPLICIT NONE
        REAL(KIND=dp), INTENT(IN), DIMENSION(:,:) :: a,b
        INTEGER :: i, j, ma, na, mb, nb
        realmateq=.true.
        ma=UBOUND(a,1)
        na=UBOUND(a,2)
        mb=UBOUND(b,1)
        nb=UBOUND(b,2)
        i=1
        IF (ma /= mb .or. na /= nb) THEN
            realmateq=.false.
        ELSE
            DO WHILE (i <= ma .and. realmateq)
                DO j=1,na
                    IF (a(i,j) /= b(i,j)) THEN
                        realmateq=.false.
                        EXIT
                    ENDIF
                ENDDO
                i=i+1
            ENDDO
        ENDIF
    END FUNCTION realmateq
    !###################################################################################################
    LOGICAL FUNCTION complexmateq(a,b)
        IMPLICIT NONE
        COMPLEX(KIND=dp), INTENT(IN), DIMENSION(:,:) :: a,b
        INTEGER :: i, j, ma, na, mb, nb
        complexmateq=.true.
        ma=UBOUND(a,1)
        na=UBOUND(a,2)
        mb=UBOUND(b,1)
        nb=UBOUND(b,2)
        i=1
        IF (ma /= mb .or. na /= nb) THEN
            complexmateq=.false.
        ELSE
            DO WHILE (i <= ma .and. complexmateq)
                DO j=1,na
                    IF (a(i,j) /= b(i,j)) THEN
                        complexmateq=.false.
                        EXIT
                    ENDIF
                ENDDO
                i=i+1
            ENDDO
        ENDIF
    END FUNCTION complexmateq
    !###################################################################################################
    LOGICAL FUNCTION realveceq(a,b)
        IMPLICIT NONE
        REAL(KIND=dp), INTENT(IN), DIMENSION(:) :: a,b
        INTEGER :: i, ma, mb
        realveceq=.true.
        ma=UBOUND(a,1)
        mb=UBOUND(b,1)
        i=1
        IF (ma /= mb) THEN
            realveceq=.false.
        ELSE
            DO i=1,ma
                IF (a(i) /= b(i)) THEN
                    realveceq=.false.
                    EXIT
                ENDIF
            ENDDO
        ENDIF
    END FUNCTION realveceq
    !###################################################################################################
    LOGICAL FUNCTION complexveceq(a,b)
        IMPLICIT NONE
        COMPLEX (KIND=dp), INTENT(IN), DIMENSION(:) :: a,b
        INTEGER :: i, ma, mb
        complexveceq=.true.
        ma=UBOUND(a,1)
        mb=UBOUND(b,1)
        i=1
        IF (ma /= mb) THEN
            complexveceq=.false.
        ELSE
            DO i=1,ma
                IF (a(i) /= b(i)) THEN
                    complexveceq=.false.
                    EXIT
                ENDIF
            ENDDO
        ENDIF
    END FUNCTION complexveceq
    !###################################################################################################
END MODULE matmodule
