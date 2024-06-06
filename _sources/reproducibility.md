# Reproducibility and Accuracy

When computing the sum of floating point values as in the scalar product
\begin{align}
\vec x\cdot \vec y = \sum_k x_k y_k
\end{align}
round-off errors can (i) falsify the numerical result and (ii) make the result non-reproducible in a parallel setting. To see the first point imagine a naive implementation
like
```{code-block} cpp
double sum = 0;
for( unsigned u=0; u<size; u++)
    sum += x[u]*y[u];
```
 $ x\oplus y$ is numerically not exact, in fact $1\oplus \varepsilon = 1$ if $\varepsilon$ is the machine precision. This means that most accuracy is lost in $s \oplus x_u y_u$ at the end of the for loop, when
the variable `sum` gets increasingly larger in comparison to `x[u]*y[u]`.
The accuracy problem is somewhat alleviated in algorithms like tree-reduction but not entirely.

Another problem is that $(a \oplus b) \oplus c \neq a\oplus (b\oplus c)$, which leads to non-reproducible results in a parallel setting where the order of execution cannot be guaranteed.
```{admonition} Further reading
You can read more about this topic in our publication  M. Wiesenberger, L. Einkemmer, M. Held, A. Gutierrez-Milla, X. Sáez, R. Iakymchuk  [Reproducibility, accuracy and performance of the Feltor code and library on parallel computer architectures](https://doi.org/10.1016/j.cpc.2018.12.006) Computer Physics Communications 23
```

## The `dg::blas1::dot` function
We base our algorithm on the "exblas" approach, which computes sums with a long accumulator (a large fixed point number) that can accurately store the result and only perform one rounding operation at the end. This makes the `dg::blas1::dot` function both **as accurate as it can be** and **binary reproducible** across many different platforms.

## Is Feltor binary reproducible?
The short answer is no. The longer answer is "almost" or "at least if you don't change the compiler. In order to be binary reproducible **all computations** have to be implemented either with long accumulators or such that the order of execution is always guaranteed. Now the problem lies in the fact that we do not always have full control over what is implemented. Consider


``` cpp
#include "dg/exblas/exblas.h"
double x = 6.12610567450009658;
dg::exblas::udouble result;
result.d = sin(x);
//std::cout << result.d << " "<<result.i<<"\n";
std::cout << "Difference to correct result: "<<result.i  - -4628567870976535683<<"\n";
// Difference to correct result 1
```

Depending on your platform/compiler the result of the above example may or may not be zero. This is because the implementation of the sine function is not unique and rounding is not always done correctly for all inputs.

This simple example shows that from a developer perspective it is almost impossible (or at least highly impractical) to make an entire program binary reproducible across all platforms and compilers.
