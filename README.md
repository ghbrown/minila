# minila

"mini linear algebra"

A light weight and standalone C implementation of common linear algebraic operations for dense matrices and vectors.

---

### Primitives

`vva`: vector vector addition

`mma`: matrix matrix addition

`mvp`: matrix vector product

`mmp`: matrix matrix product

`mme`: matrix matrix elementwise

---

*Notes*

This was originally supposed to be a C library, but implementing
even vector vector addition of arbitrarily sized arrays was
needlessly complicated. This is almost wholly because C does not
allow the determination of array size from pointers to arrays.

I'm also using this library as a way to strengthen my programming
in a standard scientific computing language and to better understand
implementations of popular linear algebra methods.
