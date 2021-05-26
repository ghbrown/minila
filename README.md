# minila

"mini linear algebra"

A light weight and standalone Fortran implementation of common linear algebraic operations for dense matrices and vectors.

---

**implemented** (otherwise it has not yet been implemented)

(f) - function

(s) - subroutine

---

## Main Functions 

Functions which will likely be directly useful to a user.

### Primitives

Simple linear algebra operations not included in Fortran built-ins.

**`dispmat(A)`**(s): displays matrix `A`

**`get_diag(A)`**(f): returns diagonal elements of matrix A as vector

**`set_diag(diag_vec)`**(f): returns square matrix of zeros with diagonal elements specified by `diag_vec`

**`rop(u1,u2)`**(f): rank one product of two vectors, produces matrix equal to `u1 u2^T`


**`diag_sub(D,b)`**(f): solve lower diagonal system `Dx=b` for `x=b/get_diag(D)`

**`bck_sub(T,b)`**(f): solve lower triangular system `Tx=b` for `x` using back substitution

**`fwd_sub(T,b)`**(f): solve upper triangular system `Tx=b` for `x` using forward substitution

**`strang(T,b)`**(f): solve triangular system `Tx=b` after detecting the triangularity type of T. Essentially wrapper around the `*sub` functions with O(n^2) checks and extra safeguards; if performance is critical use the `*sub` functions directly. **What a serindiptously great name for the function.**

### Matrix Decompositions

`lupp(A,L,U)`(s): LU decomposition of square matrix with partial pivoting

`qr(A,Q,R)`(s): QR decomposition of (nonsquare) matrix

### Linear Systems

`linsolve(A,b)`(f): solve (square) linear system `Ax=b` (combination of LU decomposition and triangular system solver)

`lssolve(A,r)`(f): solve least squares problem `Ax \approx b`

---

### Helper Functions

Functoins used internally by other functions. Not the meat and potatoes of a linear algebra library, but possibly useful.

**`is_square(A,n)`**(f): returns logical value on whether matrix A is square, if so sets `n` equal to dim(A)

**`is_triangular(A,type)`**(f): outputs a logical value, and assigns `type` to character corresponding to type of triangularity (`"q"`-not square, `"n"`-not triangular, `"u"`-upper triangular, `"l"`-lower triangular, `"d"`-diagonal)

---

*Notes*

I'm mainly using this library as a way to strengthen my programming
in a standard scientific computing language and to better understand
implementations of popular linear algebra methods.

This was originally supposed to be a C library, but implementing
even vector vector addition of arbitrarily sized arrays was
needlessly complicated. This is almost wholly because C does not
allow the determination of array size from pointers to arrays (I think I have figured our the proper way to get around this now, so maybe I'll write a C version in the future.)

