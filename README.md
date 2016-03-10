
A Fortran 2003 module to make Fortran matrix manipulation easier.

A) What is it?
Have you ever wished Fortran would have even an easier interface? Have you ever missed the almighty backslash “\” operator which is available in MatLab? Have you ever wished more abstractions in matrix operations in Fortran? This module is a humble effort to fulfill your wishes!

B) What can it do? 

For all of the following functions or subroutines, the argument can be either real or complex but of double precision type.

1) .m. operator : multiplying two real/complex matrices

    REAL (kind=SELECTED_REAL_KIND(14)), DIMENSION(1000,1000)  :: a , b , c
    !!! suppose you have filled a and b with data and want c=a*b
    c= a .m. b

2) .d. operator : dividing two real/complex matrices

    COMPLEX (kind=SELECTED_REAL_KIND(14)), DIMENSION(1000,1000)  :: a , b , x
    !!! suppose you have filled a and b with data and want to solve ax=b
    x= b .d. a

3) LU of a matrix :

    lu(a)

4) Determinant of a real/complex matrix:

    det(a)

5) Inverse of a matrix:

    inv(a)

6) Are two matrices a and b equal?

    IF (a .meq. b) 

7) writing a matrix to a file in a nice, distinct, understandable way:

    CALL matwrite (a)

8) the almighty eig function!

    CLASS (eigen), POINTER :: ans
    COMPLEX(KIND=dp), DIMENSION(n,n) :: a, lvec, rvec
    COMPLEX(KIND=dp), DIMENSION(n) :: valc
    ! suppose you have filled a with data 
    ans => eig(a)
    lvec=ans%lgetleigvec() ! Get left eigen vectors 
    rvec=ans%getreigvec()  ! Get right eigen vectors
    valc=ans%geteigval()    ! Get the eigen values

If you don't know the eigen type (symmetric or non symmetric, complex or real), define the needed arrays/matrices and use the select case method:
    
    sol=>eig(ar)
    SELECT TYPE (sol)
        CLASS IS  (eigen_c)
            PRINT*, 'Complex eigenvalue, Nonsymmetric Complex eigenvectors.'
            valc=sol%geteigval()
            lvec=sol%getleigvec()
            rvec=sol%getreigvec()
        CLASS IS (eigen_c_h)
            PRINT*, 'Real eigenvalue, Symmetric Complex eigenvectors (Hermitian matrix).'
            valr=sol%geteigval()
            vecc=sol%geteigvec()
        CLASS IS (eigen_r_s)
            PRINT*, 'Real eigenvalue, Symmetric Real eigenvectors.'
            valr=sol%geteigval()
            vecr=sol%geteigvec()
        CLASS IS (eigen_r)
            PRINT*, 'Real eigenvalue, Nonsymmetric Real eigenvectors.'
            valr=sol%geteigval()
            lvecr=sol%getleigvec()
            rvecr=sol%getreigvec()
    END SELECT

C) How to compile?
MatModule relies on LAPACK.

    gfortran yourprogram.f90 matmodule.f90 -o yourprogram `pkg-config --libs lapack`
