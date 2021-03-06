# minila

"mini linear algebra"

A standalone Fortran implementation of common linear algebraic operations for dense matrices and vectors.

---

**implemented** (otherwise it has not yet been implemented)

(f) - function

(s) - subroutine

---

## Main Functions 

### Primitives

Simple linear algebra operations not included as Fortran intrinsics. Building blocks for other modules.

**`disp_mat(A)`**(s): displays matrix `A`

**`get_diag(A)`**(f): returns diagonal elements of matrix A as vector

**`set_diag(diag_vec)`**(f): returns square matrix of zeros with diagonal elements specified by `diag_vec`

**`identity(n)`**(f): returns the identity matrix of dimension `n`

**`rop(u1,u2)`**(f): rank one product of two vectors, produces matrix equal to `u1 u2^T`

**`diag_sub(D,b)`**(f): solve diagonal system `Dx=b` for `x=b/get_diag(D)`

**`bck_sub(T,b)`**(f): solve lower triangular system `Tx=b` for `x` using back substitution

**`fwd_sub(T,b)`**(f): solve upper triangular system `Tx=b` for `x` using forward substitution

**`strang(T,b)`**(f): solve triangular system `Tx=b` after detecting the triangularity type of T. Essentially wrapper around the `*sub` functions with O(n^2) checks and extra safeguards; if performance is critical use the `*_sub` functions directly. **What a serindiptously great name for the function.**


### Matrix Decompositions

**`lupp(A,L,U,p)`**(s): LU decomposition of square matrix with partial pivoting, results satisfy `PA=A(p,:)=LU` where `P` is the matrix corresponding to pivot vector `p` 

`ldl(A,L,D)`(s): LDL decomposition for symmetric matrix as `A=L D L^T`

`qr(A,Q,R)`(s): QR decomposition of matrix as `A=QR`

`eig(A,V,d)`(s): eigenvalue/spectral decomposition of matrix as `A = V diag(d) V^-1`

`svd(A,U,d,V)`(s): singular value decomposition of a matrix as `A = U diag(d) V^T`


### Linear Systems

**`linsol(A,B)`**(f): solves linear system `AX=B` via LU decomposition, the solution `X` will be of the same shape as the right hand side `B`, which can be a vector (single linear system) or a matrix (set of linear systems)

**`inverse(A)`**(f): computes the matrix inverse of `A`

**`sor(A,b,rel_tol,max_iter)`**(f): solves `Ax=b` via successive over-relaxation (an iterative method); terminates when relative error is less than `rel_tol` or when `max_iter` iterations have been performed


### Least Squares

`lssol(A,r)`(f): solve least squares problem `Ax \approx b`

---

## Helper Functions

Used internally by other functions, but not the meat and potatoes of a linear algebra library.
All included in `prim.f90` module currently.

**`is_square(A,n)`**(f): returns logical value on whether matrix A is square, if so sets `n` equal to dim(A)

**`is_triangular(A,type)`**(f): outputs a logical value, and assigns `type` to character corresponding to type of triangularity (`"q"`-not square, `"n"`-not triangular, `"u"`-upper triangular, `"l"`-lower triangular, `"d"`-diagonal)

**`get_triang(A,type)`**(f): returns lower (`type="l"`) or upper (`type="u"`) triangle (NOT including the diagonal) of matrix `A` as a matrix

---

*Notes*

I'm mainly using this library as a way to strengthen my programming
in a standard scientific computing language and to better understand
implementations of popular linear algebra methods.

This was originally supposed to be a C library, but implementing
even vector vector addition of arbitrarily sized arrays was
needlessly complicated (almost wholly because C does not
allow the determination of array size from pointers to arrays). I think I have figured our the proper way to get around this now, so maybe I'll write a C version in the future.

