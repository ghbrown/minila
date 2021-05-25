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

**`dispmat(A)`**(s): displays matrix A

**`rop(u1,u2)`**(f): rank one product of two vectors, produces matrix equal to (u1)transpose(u2)

**`strang(T,b)`**(f): solve triangular system Tx=b given (upper *or* lower) triangular matrix T (what a serindiptously great name for the function)

### Matrix Decompositions

`lupp(A,L,U)`(s): LU decomposition of square matrix with partial pivoting

`qr(A,Q,R)`(s): QR decomposition of (nonsquare) matrix

### Linear Systems

`linsolve(A,b)`(f): solve (square) linear system `Ax=b` (combination of LU decomposition and triangular system solver)

`lssolve(A,r)`(f): solve least squares problem `Ax \approx b`

---

### Other Functions

Far more situational than the functions above, but used internally by other functions.

**`is_triangular(A,type)`**(f): outputs a logical value, and assigns `type` to character corresponding to type of triangularity (`n`-not triangular, `u`-upper triangular, `l`-lower triangular, `d`-diagonal)

---

*Notes*

I'm mainly using this library as a way to strengthen my programming
in a standard scientific computing language and to better understand
implementations of popular linear algebra methods.

This was originally supposed to be a C library, but implementing
even vector vector addition of arbitrarily sized arrays was
needlessly complicated. This is almost wholly because C does not
allow the determination of array size from pointers to arrays (I think I have figured our the proper way to get around this now, so maybe I'll write a C version in the future.)

