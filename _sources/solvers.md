# Solvers in Feltor

In this chapter we will study how to solve
\begin{align}
f(x) = b
\end{align}
for both linear and nonlinear functions $f$ in Feltor.

## Overview

Feltor offers a growing range of solvers. Among others there are linear
- `dg::PCG` for preconditioned conjugate gradient
- `dg::LGMRES` for the restarted generalized minimal residual method
- `dg::BICGSTABL`
- `dg::ChebyshevIteration`

and non-linear solvers
- `dg::FixedPointSolver`
- `dg::AndersonAcceleration`

```{admonition} Matrix-free solvers
It is important to realize that Feltor's solvers are _matrix-free_ (and in consequence iterative). What this means is that the actual matrix does not need to be allocated. In fact, the matrix does not need to be known at all, in code or analytically. The only thing that is needed is the **application** of the matrix to a vector. In other words, the solver just needs a "black box" that takes a vector on one side and spits out a vector on the other side. Linearity is just a promise of the "black box" (which our solvers blindly trust).
```

These solvers have in common that they are used in some variation of
```cpp
Operator op; // some operator f to invert
Solver solver( copyable, construct_params...);
Vector x, b;
double rtol = 1e-8, atol = rtol;
solver.solve( op, x, b, solve_params...);
```
Importantly, many solvers take a preconditioner and the weights as solve parameters. The latter define the scalar product
\begin{align}
   \langle x, y \rangle := \sum_i x_i w_i y_i
\end{align}
in which in particular the error norm $|| v || = \sqrt{ \langle v, v\rangle }$ is computed.

Section [Operators](sec:operator) details the format and how to write the function $f$ in code, while [Preconditioners](sec:preconditioner) deal with the inverse operators.
```{seealso}
The [doxygen documentation](https://mwiesenberger.github.io/feltor/dg/html/group__invert.html) holds the complete list of solvers and details the parameters they use.
```

Furthermore, we have the **multigrid family** of solvers, which we designed as a customizable combination of the above solvers into
- `dg::nested_iterations`
- `dg::fmg_solve` (experimental).

Here, we refer to section [The multigrid solvers](sec:multigrid) for how to construct and use these solvers.
```{seealso}
The multigrid section in the [doxygen documentation](https://mwiesenberger.github.io/feltor/dg/html/group__multigrid.html)  holds the complete list of solvers and details the parameters they use.
```

```{admonition} Initial guess
Look at `dg::Extrapolation` to get an initial guess based on past solutions.
```

## A first example: inverting the Laplace operator

An often encountered problem is the discretization and inversion of elliptic equations, for example Poison's equation.
\begin{align}
-\Delta \phi = \rho
\end{align}
Using a local discontinuous Galerkin discretization this equation transforms to
\begin{align}
 M \vec \phi = \vec \rho
 \end{align}
where $M$ is self-adjoint in the dG weights $W$. See the [theory guide](https://www.overleaf.com/read/rpbjsqmmfzyj) on overleaf for more details on the discretization.
Since $M$ is self-adjoint we can use `dg::PCG` to invert the equation:
``` cpp
#include <iostream>
// include the basic dg header
#include "dg/algorithm.h"
```
We want to solve on the domain $[0,2\pi]\times [0,2\pi]$
``` cpp
const double lx = 2.*M_PI;
const double ly = 2.*M_PI;
```
With the somewhat trivial solution $\phi = \sin(x)\sin(y)$
``` cpp
double fct(double x, double y){
    return sin(y) * sin(x);
}
double laplace_fct( double x, double y) {
    return 2 * sin(y) * sin(x);
}
double initial( double x, double y) {
    return sin(0);
}

// The numbers in parentheses are our suggestions ...
unsigned n = 3, Nx = 10, Ny = 10;
const double eps = 1e-4;
// clang as script is unfortunately not very fast so don't make the numbers too big
//std::cout << "Type n (3), Nx (20) and Ny (20)! \n";
//std::cin >> n >> Nx >> Ny;
std::cout << "Computing on the Grid " <<n<<" x "<<Nx<<" x "<<Ny <<"\n";
dg::CartesianGrid2d grid( 0, lx, 0, ly, n, Nx, Ny, dg::DIR, dg::PER);
std::cout<<"Evaluate initial guess for iterative scheme\n";
dg::DVec x = dg::evaluate( initial, grid);
// create volume and inverse volume on previously defined grid
const dg::DVec vol2d = dg::create::volume( grid);

// Create negative, unnormalized, positive definite Laplacian
dg::Elliptic<dg::CartesianGrid2d, dg::DMatrix, dg::DVec> laplaceM( grid);
// allocate memory in conjugate gradient
unsigned max_iter = n*n*Nx*Ny;
dg::PCG<dg::DVec > pcg( x, max_iter);

// Evaluate right hand side and solution on the grid
dg::DVec b = dg::evaluate ( laplace_fct, grid);
const dg::DVec solution = dg::evaluate ( fct, grid);
dg::Timer t;
t.tic();
// we do not use a preconditioner in solution method
unsigned number = pcg.solve( laplaceM, x, b, 1., vol2d, eps);
t.toc();
std::cout << "Number of pcg iterations "<< number <<"\n";
std::cout << "For a precision of "<< eps<<"\n";
std::cout << "Took "<< t.diff()<<"s\n";

//Computing on the Grid 3 x 10 x 10
//Evaluate initial guess for iterative scheme
//Number of pcg iterations 15
//For a precision of 0.0001
//Took 1.82155s
```
This program executes on the device as is evident from the use of
 `dg::DVec`, a typedef for `thrust::device_vector<double>`, and `dg::DMatrix`.
 The new class that we encounter in this program is `dg::Elliptic`, a
 level 4 class. What this class does is to create and store the `dx` and `dy` matrices from the given grid. Furthermore, it allocates some
 internal workspace and defines a member function that uses `dg::blas1` and `dg::blas2` functions to discretize the Laplacian.
 The class depends on three template parameters,
 the geometry, the matrix and the vector class. The matrix and the
 vector class must fit together, that is, they must be useable in
 a `dg::blas2::symv` function.  
 ```{note}
 The `dg::Elliptic` class
 acts as a matrix itself. This becomes clear in the line `dg::blas2::symv( laplaceM, x, lap_x)` in the example below. The `M` in `laplaceM`
 reminds us that `dg::Elliptic` discretizes the negative Laplacian in order to
 generate a postive definite operator.
 ```

 The geometry `dg::CartesianGrid2d` defines a two-dimensional Euclidean metric
 tensor. Recall, that in general the Laplacian depends on metric coefficients.

Since the `dg::Elliptic` is a matrix, we can use it in a
conjugate gradient solver. The `dg::PCG` is a Level 2 class and
only depends on the vector class that we use. Recall here that the
conjugate gradient algorithm can be implemented with vector addition `dg::blas1::axpby`,
a scalar product (`dg::blas1::dot`) and a matrix-vector multiplication (`dg::blas2::symv`) alone.

``` cpp
//compute error
dg::DVec error( solution);
dg::blas1::axpby( 1.,x,-1.,error);

dg::DVec lap_x(x), residuum( b);
dg::blas2::symv(  laplaceM, x, lap_x);
dg::blas1::axpby( 1., lap_x, -1., residuum);

//global relative error in L2 norm is O(h^n)
double result;
result = sqrt(dg::blas2::dot( x, vol2d, x));
std::cout << "L2 Norm of x is               " << result << std::endl;
result = sqrt(dg::blas2::dot(solution, vol2d , solution));
std::cout << "L2 Norm of Solution is        " << result << std::endl;
result = sqrt(dg::blas2::dot(error, vol2d , error));
std::cout << "L2 Norm of Error is           " << result << std::endl;
result = sqrt(dg::blas2::dot( residuum, vol2d, residuum));
std::cout << "L2 Norm of Residuum is        " << result << std::endl;

//L2 Norm of x is               3.14152
//L2 Norm of Solution is        3.14159
//L2 Norm of Error is           0.00348843
//L2 Norm of Residuum is        0.000467984

```

(sec:operator)=
## Writing an operator to invert
Feltor has an assortment of already made Operators available that have special implementations, for example the
`dg::Elliptic` or `dg::Helmholtz` classes. Sometimes, however your own problem deviates slightly from what is already available or you want to start entirely from scratch.

For example in the first example we inverted the simple Laplacian. However, what if you want to solve the equation
\begin{align}
 \exp(\phi) - \Delta \phi = \rho
 \end{align}
which is non-linear. This means our function is now $f(x) = \exp(x)-\Delta x$ instead of $f(x) = \Delta x$. Feltor accepts several ways a function or matrix is implemented. One way you could implement the new $f$ based on the existing old $f$ (or in fact any $f$) is with **lambda functions**.
````{admonition} dg::apply
In Feltor the function that enables the implementation of a matrix-free algorithm is `dg::blas2::symv` or its alias `dg::apply`:
```cpp
    dg::apply( op, x,y); // compute y = f(x)
    dg::blas2::symv( op, x,y );//completely equivalent
```
`dg::apply` will accept many types, you will find the full list in the doxygen documentation.
Notably, it accepts functors and will try to call
```cpp
op(x,y); // compute y = f(x)
```
i.e. the apply call is equivalent to a functor call.
In order for this to work your own matrix class must be a functor/lambda with the signature `void operator()(const container&,container&)` where `container` is a placeholder for the vector class you use
````

### Lambda functions
In C++14 lambda functions are a very powerful tool.
They are basically compiler generated functors that
can be written in a very concise and short way.
```{seealso}
If you are unfamiliar with lambdas a good watch on youtube is for example [Lambdas from Scratch - Arthur O'Dwye](https://www.youtube.com/watch?v=3jCOwajNch0)
```

```cpp
auto op = [&lap = laplaceM, tmp=VectorType(copyble)]( const auto& x, auto& y) mutable
{
    dg::blas2::symv( lap, x, y);
    dg::blas1::transform( x, tmp, dg::EXP<double>());
    dg::blas1::axpby( 1., tmp, 1., y);
}
```
Several things are noteworthy here

  - Our lambda is a template or a "generic lambda" through the use of "auto" in the interface
  - Our lambda uses "init capture", a feature from C++14. We caught the elliptic object by reference (in order to avoid a deep copy but we should beware "dangling references") and allocated a "workspace" object tmp
  - Our lambda is declared "mutable" (because we write in the temporary tmp)


Note that the lambda now has the correct interface to be used in a solver. For example
```cpp
dg::AndersonAcceleration acc( copyable, 3);
acc.solve( op, phi, rho, weights, 1e-5, 1e-5, 1000, 1e-5, 3);
```

### Overload operator()
The equivalent of a lambda function is of course to write an operator yourself. The above lambda can be implemented like:
```cpp
template<class Geometry, class Matrix, class Container>
struct Operator
{
    Operator( dg::Elliptic<Geometry,Matrix,Container>& ell, const Container& copyable):
        m_tmp( copyable), m_lap(ell) {}
    template<class ContainerType0, class ContainerType1>
    void operator()( const ContainerType0& x, ContainerType1& y)
    {
        dg::blas2::symv( m_lap, x, y);
        dg::blas1::transform( x, m_tmp, dg::EXP<double>());
        dg::blas1::axpby( 1., m_tmp, 1., y);
    }
    private:
    Container m_tmp;
    dg::Elliptic<Geometry,Matrix,Container>& m_lap;
};
```
As you can see this is a lot more verbose than a lambda function. It is useful for complicated functors that hold
internal state. It is less useful for very small classes, where the parenthesis operator only consists of a line or two and when the operator is basically only
an adaptor of another type.

(sec:preconditioner)=
## Preconditioning and inverse operators

A preconditioner is a matrix that represents the approximate inverse of the matrix to invert.
From an
interface point of view, the preconditioner has the same conditions
as the matrix itself, i.e. must be usable in the `dg::blas2::symv` function
```cpp
Precond p(...);
dg::blas2::symv( p, x, y); // computes y = P(x)
```
This example works if `using Precond = Vector` is a
vector, which will be interpreted as a **diagonal preconditioner**.
For example the `dg::Elliptic` class comes with a default diagonal preconditioner:
```cpp
unsigned number = pcg.solve( laplaceM, x, b, laplaceM.precond(), vol2d, eps);
```
Here, the `precond` method returns a vector which contains the inverse of the volume per default.

But what if we want to use a more involved preconditioner? In fact, the interface of the preconditioner
must look exactly the same as that of the operator, which we discussed in the [previous section](sec:operator).
When constructing a preconditioner manually we thus have again two options: (i) write an entirely new class
and overload the function call operator or (ii) use lambda functions.

Consider the following piece of code
```cpp
dg::PCG<dg::DVec > pcg( x, max_iter);
auto inverse_op = [&] ( const auto& y, auto& x)
{
    pcg.solve( laplaceM, x, y, 1., vol2d, eps);
}
```
`inverse_op` now has the signature of a matrix in Feltor. That is we can write
```cpp
dg::blas2::symv( inverse_op, b, x);
```
which will compute $x = \Delta^{-1} b$. This means that `inverse_op` could be used as a preconditioner.
Of course, here the "preconditioner" is the exact inverse of the operator which means that any iterative scheme converges in one iteration.

```{admonition} Preconditioner
A preconditioner is nothing but an approximate inverse to a matrix.
```
We should remember that the preconditioner only needs to approximately invert the equation.
In the next example, we show
how to use a lambda to combine an elliptic operator with a Chebyshev solver that applies 10 iterations as a preconditioner that can be used in PCG (solves $-\Delta \Phi = \rho $)

```cpp
// The operator that we want to solve
dg::Elliptic<...> pol( grid);
pol.set_chi( chi);
//estimate EigenValue
dg::EVE<Vector> eve( copyable);
double eps_ev = 1e-2, ev;
// compute largest Eigenvalue in a few iterations
counter = eve.solve( pol, x, b, 1., pol.weights(), ev, eps_ev);
// Here we combine the Elliptic operator with the Chebyshev solver
dg::ChebyshevIterations<Vector> cheby( x);
auto precond = [nu=10, ev, &pol, &cheby](const auto& y, auto& x)
        {
            cheby.solve( pol, x, y, 1., ev/100., ev*1.1, nu+1, true);
        };
// Now we can call the solver
dg::PCG< Vector > pcg( x, 1000);
pcg.solve( pol, x, y, precond, pol.weights(), eps, 1, 1);
```

```{admonition} Is it faster?
Preconditioning aims at reducing the number of iterations. But the real question is:
Is the solve with preconditioner _faster_ than without preconditioner? We emphasize here that _faster_ is not the same as _less iterations_. In fact, the above example uses less PCG iterations than the unpreconditioned CG method, but is overall slower to execute (on desktop PCs).
Predicting why code is slower or faster is not always an easy task.
The hope in the above code for example is that the lack of dot products in the Chebyshev Iteration makes the code faster (because it lacks global communication) and for high enough process number the above code might actually perform better than unpreconditioned CG.
The only sure way to tell is to try and measure.
```

(sec:multigrid)=
## The multigrid solvers
The basic feature of a multigrid solver is that it (approximately) solves the equation on more than one grid.
The idea is that on a coarse grid a solver generally works faster than on a fine grids with more resolution.

The easiest multigrid solver is called _nested iterations_ or _full approximation scheme_ (for the non-linear case). Here the solution on a coarse grid is used as an improved initial guess for the fine grid solution.


### Construct Nested Grids
In Feltor the first thing to do for a multigrid solver is to construct multiple grids:
```cpp
unsigned stages = 3;
Geometry grid(...);
dg::NestedGrids<Geometry, Matrix, Container> nested( grid, stages);
```
The nested grids provide the functionality to project and interpolate between
various grids and also provides a workspace for a solve function.

### Construct Nested operators and solvers
We use the nested grids to construct an operator and a solver on each grid.
Like a preconditioner the solver needs to be represented as an inverse operator
```cpp
Container chi = ...;
std::vector<Container> multi_chi = nested.project( chi);
std::vector<dg::Elliptic<Geometry, Matrix, Container> >
    multi_pol( stages);
std::vector<std::function<void( const Container&, Container&)> >
    multi_inv_pol(stages);
for(unsigned u=0; u<stages; u++)
{
    multi_pol[u] = { nested.grid(u), dg::centered};
    multi_pol[u].set_chi( multi_chi[u]);
    multi_inv_pol[u] = [eps = multi_eps[u],
                        pcg = dg::PCG<Container>( multi_chi[u], 1000 ),
                        &pol = multi_pol[u]) ]
        ( const auto& y, auto& x) mutable
        {
            // Note that we do not need to provide an intial guess ourselves.
            //It is already contained in x and will be provided by the multigrid solver.
            pcg.solve( pol, x, y, pol.precond(), pol.weights(),
                eps, 1., 1 );
        };
}
```

Note how we init the solver in-place and capture the operator by reference in the lambda.
This
is very similar to how we constructed the preconditioner in section [preconditioner](sec:preconditioner).
```{admonition} Nested Type
Note how we store the lambdas in a vector of `std::function`, which erases the type of the lambda.
Therefore, it is possible to use different kind of solvers at each stage.
```
Often, the lambda body is used to provide more refined functionality.
```cpp
    multi_inv_pol[u] = [u=u, eps = multi_eps[u],
                        pcg = dg::PCG<Container>( multi_chi[u], 1000 ),
                        &pol = multi_pol[u]) ](
        const auto& y, auto& x)
    {
#ifdef MPI_VERSION
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif //MPI
        dg::Timer t;
        t.tic();
        if( u==0)
            number = pcg.solve( pol, x, y, pol.precond(), pol.weights(),
                eps, 1., 1 );
        else // check stopping criterion only every 10th iteration
            number = pcg.solve( pol, x, y, pol.precond(), pol.weights(),
                eps, 1., 10 );
        t.toc();
        DG_RANK0 std::cout << "# Nested iterations stage: " << u << ", iter: " << number << ", took "<<t.diff()<<"s\n";
    };
```
```{seealso}
Combining elliptic operators and PCG solves with nested iteration appears so often in application codes that we packaged it into the utilitiy class `dg::MultigridCG2d`. Check it out in the doxygen documentation.
```


### Solve
The final step is simply to call the appropriate function
```cpp
dg::nested_iterations( multi_pol, x, b, multi_inv_pol, nested);
```
or in the functional programming style make the multigrid itself an inverse operator
```cpp
auto inverse_op = [&] (const auto& y, auto& x)
{
    dg::nested_iterations( multi_pol, x, y, multi_inv_pol, nested);
};
```

```{admonition} Combine and build
The possibility to use lambdas in our solve methods opens up the possibility to combine various solvers arbitrarily in Feltor. For example we could use `dg::fmg_solve` as a preconditioner for a `dg::PCG` solver.
```
