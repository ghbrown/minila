# minila

"mini linear algebra"

A light weight and standalone Fortran implementation of common linear algebraic operations for dense matrices and vectors.

---

### Matrix Decompositions

`LU(A)`: LU decomposition of square matrix

`QR(A)`: QR decomposition of (nonsquare) matrix

### Linear Systems

`linsolve(A,b)`: solve (square) linear system `Ax=b`

`lssolve(A,r)`: solve least squares problem `Ax \approx b`


---

*Notes*

This was originally supposed to be a C library, but implementing
even vector vector addition of arbitrarily sized arrays was
needlessly complicated. This is almost wholly because C does not
allow the determination of array size from pointers to arrays.

I'm also using this library as a way to strengthen my programming
in a standard scientific computing language and to better understand
implementations of popular linear algebra methods.
