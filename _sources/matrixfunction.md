# Matrix functions

We consider the solution of the linear matrix function equation
\begin{align}\label{eq:matrixfunc}
 \vec{x} = f(A) \vec{b}
\end{align}
for the vector $\vec{x}$ and $\vec{b}$ with a given square matrix $A\in \mathbb{R}^{n \times n}$. The matrix $A$ is assumed to be self-adjoint in the scalar product defined by the symmetric and positiv definit matrix $M$
\begin{align}
    \langle \vec v, \vec w\rangle_M := \vec v^\mathrm{T} M \cdot \vec w ,\quad
    ||\vec v||_M := \sqrt{ \langle \vec v, \vec v\rangle_M}
\end{align}

```{admonition} Include and Link
When you use the matrix-function extension you will need to have [boost](https://www.boost.org/) and lapack installed on your system. On linux, installing `libboost-dev` and `liblapacke-dev` should suffice. When configuring on a non-default system add lapack to the `LIBLAPACK` configuration variable in the Makefile.
```


```cpp
#include "dg/algorithm.h"
#include "dg/matrix/matrix.h" //requires boost and lapacke
```

to be continued ... (check out doxygen documentation for now)
