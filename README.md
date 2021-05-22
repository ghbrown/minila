# minila

"mini linear algebra"

A light weight and standalone Fortran implementation of common linear algebraic operations for dense matrices and vectors.

---

**implemented** (otherwise it has not yet been implemented)

(f) - function

(s) - subroutine

### Primitives

Simple linear algebra operations not included in Fortran built-ins.

**`dispmat(A)`**(s): displays matrix A

**`rop(u1,u2)`**(f): rank one product of two vectors, produces matrix equal to (u1)transpose(u2)

`strang(T,b)`(f): solve triangular system Tx=b given triangular matrix T (what a serindiptously great name for the function)

### Matrix Decompositions

`LU(A,L,U)`(s): LU decomposition of square matrix

`QR(A,Q,R)`(s): QR decomposition of (nonsquare) matrix

### Linear Systems

`linsolve(A,b)`(f): solve (square) linear system `Ax=b` (combination of LU decomposition and triangular system solver)

`lssolve(A,r)`(f): solve least squares problem `Ax \approx b`

---

*Notes*

This was originally supposed to be a C library, but implementing
even vector vector addition of arbitrarily sized arrays was
needlessly complicated. This is almost wholly because C does not
allow the determination of array size from pointers to arrays.

I'm mainly using this library as a way to strengthen my programming
in a standard scientific computing language and to better understand
implementations of popular linear algebra methods.
