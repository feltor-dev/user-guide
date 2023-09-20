# Grids and derivatives

In the first two chapters we saw how to deal with vectors and how to integrate
our first ordinary differential equation. But how about partial
differential equations? In order to define a partial differential
equation we need an actual physical space in two dimensions say and
define derivatives on it. However, there is such a multitude of possible
ways to discretize the physical space and the derivatives on it that we cannot
hope to cover them all in one library. We have chosen so-called
**discontinuous Galerkin** methods. You never heard of those before?
In simple words these methods use
a polynomial of arbitrary order in each grid cell instead of just one single
point as do finite difference methods. The order of the polynomial
defines the order of the method.
To help you get started we have written an
[introduction to dG methods](https://www.overleaf.com/read/rpbjsqmmfzyj).
## Function evaluation and integration
Let us look at a first example of what we can do with these methods. Let's
integrate a function:

```{code-block} cpp
#include <iostream>
//include the dg-library
#include "dg/algorithm.h"
//define a function to integrate
double function(double x, double y){
  return exp(x) * exp(y);
}

//create a discretization of [0,2]x[0,2] with
//3 polynomial coefficients and 20 cells in x and y
dg::CartesianGrid2d g2d( 0, 2, 0, 2, 3, 20, 20);
//discretize a function on this grid
const dg::DVec f = dg::evaluate( function, g2d);
//create the volume element
const dg::DVec vol2d = dg::create::volume( g2d);
//compute the integral on the device
double integral = dg::blas1::dot( vol2d, f);
// output approximates (exp(2)-exp(0))^2
std::cout << integral <<std::endl;
// outputs 40.82
```
In this program we encounter a new class, the `dg::CartesianGrid2d` and
two new functions, `dg::evaluate` and `dg::create::volume`.
The grid class represents the topology (what point is neighbor to what other
  point) and the geometry (the metric) that we use. In this
case it is a two-dimensional Cartesian geometry, that is, the metric tensor
  is the unit tensor.  Furthermore, it is a structured product space grid,
  which means that the topology is given implicitly by the grid coordinates.
  In the above example there are `60x60=3600` grid points in total.
This number comes from the 20 cells that we use in both x and y direction
and the 3 polynomial coefficients (that we use for both directions).
````{admonition} Choose polynomials
Per default the number of polynomials is the same in both x and y direction.
If you want the freedom to choose each individually you could use
```cpp
dg::CartesianGrid2d g2d( {0, 2, 3, 20}, {0, 2, 5, 18});
```
Now we have a grid with 5 coefficients and 18 cells in the y direction.
This highlights that the two-dimensional grid is made as a product
of 2 one-dimensional grids.
````
The `dg::evaluate` function is a template function that takes a **binary**
function or Functor as a first parameter.  In our example we aptly named
the function `function` and it depends on the x and y coordinate and
returns the result. The `dg::evaluate` function now takes our Cartesian grid,
constructs the grid coordinates from it and inserts the coordinate into the
function, one by one. The result is a vector containing the values of `function`
for all the grid coordinates. We effectively discretized our `function` on the
grid.
```{admonition} Evaluation direction, do you really need it?
If you think of the computational space as a box the first point in the result vector
corresponds to the lower left corner. The next one is the point to the right, that
is the x-direction is contiguous in memory. However, if you find that you need
this information, stop! You should not need it.
Use `dg::interpolate` or `dg::create::interpolation` if you want to access individual points
and the `dg::blas` functions to operate on them.
```
The last function is the `dg::create::volume` function. This function takes
our Cartesian grid and computes the volume form from the metric and also takes
the Gaussian weights (from our discontinuous Galerkin methods) into account.
The volume form is the `dxdy` in an integration of `Int f(x,y) dxdy`.

Finally, we use the familiar `dg::blas1::dot` function to compute the scalar
product between the function and the volume form and actually get
a very good approximation of the analytical result.

If we want to use MPI we have to change our code a little bit. When we use
MPI we have to distribute the vector `f` among all processes in the
communicator. We do this in the most direct way, namely simply allocating
an equally sized portion of the grid to every process in the communicator.
Then, each process can evaluate `function` only on the portion of the
grid that it was assigned to. The result looks like this:
```cpp
#include <iostream>
//activate MPI in FELTOR
#include "mpi.h"
#include "dg/algorithm.h"

double function(double x, double y){
    return exp(x)* exp(y);
}
int main(int argc, char* argv[])
{
    //init MPI and create a 2d Cartesian Communicator assuming 4 MPI threads
    MPI_Init( &argc, &argv);
    int periods[2] = {true, true}, np[2] = {2,2};
    MPI_Comm comm;
    MPI_Cart_create( MPI_COMM_WORLD, 2, np, periods, true, &comm);
    //create a 2d discretization of [0,2]x[0,2] with
    // 3 polynomial coefficients and 20 cells in x and y.
    // Each process gets 10 by 10 cells
    dg::CartesianMPIGrid2d g2d( 0, 2, 0, 2, 3, 20, 20, comm);
    //discretize a function on this grid
    const dg::MDVec f = dg::evaluate( function, g2d);
    //create the volume element
    const dg::MDVec vol2d = dg::create::volume( g2d);
    //compute the square L2 norm
    double norm = dg::blas1::dot( vol2d, f);
    //on every thread norm is now: (exp(2)-exp(0))^2
    //be a good MPI citizen and clean up
    MPI_Finalize();
    return 0;
}
```

## Derivatives
The next step after evaluating functions is to compute derivatives of course.
Consider this example code, which computes the Arakawa bracket
 $\{ f , g\} = \partial_x f \partial_y g - \partial_y f \partial_x g$

```cpp
 //define some test functions
double left( double x, double y) {
  return sin(x) * cos(y);
}
double right( double x, double y) {
  return sin(y) * cos(x);
}
// The analytical solution
double jacobian( double x, double y) {
  return cos(x) * cos(y) * cos(x) * cos(y) - sin(x) * sin(y) * sin(x) * sin(y);
}
double lx = 2.*M_PI, ly = lx;
unsigned n = 3, Nx = 40, Ny = 40;
//boundary conditions
dg::bc bcx = dg::PER, bcy = dg::PER;

//create a grid with Cartesian geometry
const dg::CartesianGrid2d grid( 0, lx, 0, ly, n, Nx, Ny, bcx, bcy);
//evaluate functions on the grid coordinates
const dg::DVec lhs = dg::evaluate( left, grid);
const dg::DVec rhs = dg::evaluate( right, grid);
const dg::DVec sol = dg::evaluate( jacobian, grid);


//allocate workspace
dg::DVec dxlhs(lhs), dxrhs(lhs), dylhs(lhs), dyrhs(lhs), jac(lhs);
//create the derivatives
dg::DMatrix dx = dg::create::dx( grid), dy = dg::create::dy(grid);
//apply the derivative to the functions
dg::blas2::symv( dx, lhs, dxlhs);
dg::blas2::symv( dy, lhs, dylhs);
dg::blas2::symv( dx, rhs, dxrhs);
dg::blas2::symv( dy, rhs, dyrhs);
//combine the results
dg::blas1::pointwiseDot( 1./3., dxlhs, dyrhs, -1./3., dylhs, dxrhs, 0., jac);
dg::blas1::pointwiseDot( 1./3.,   lhs, dyrhs, -1./3., dylhs,   rhs, 0., dylhs);
dg::blas1::pointwiseDot( 1./3., dxlhs,   rhs, -1./3.,   lhs, dxrhs, 0., dxrhs);
//add the remaining derivatives
dg::blas2::symv( 1., dx, dylhs, 1., jac);
dg::blas2::symv( 1., dy, dxrhs, 1., jac);

//create the volume form
dg::DVec vol = dg::create::volume(grid);
//test conservative property
double integral = dg::blas1::dot( vol, jac);
std::cout << "Integrated Jacobian is "<<integral<<"\n";
// should output 0 or something close to 1e-16

//now compute the distance to solution in L2 norm
dg::blas1::axpby( 1., sol, -1., jac);
double error = sqrt(dg::blas2::dot( jac, vol, jac));
std::cout << "Distance to solution "<<error<<"\n";
// Distance to solution 0.000754106
```

Here, we encounter the `dg::create::dx` and `dg::create::dy` functions that,
as the names suggest, create derivatives in x and y direction respectively.
The only parameter is the grid. Per default the boundary condition for the
derivative is periodic and we use a centered discretization.
The type of the matrix is a `dg::DMatrix` which is a typedef for an obscure
dg-specific data type that holds a block-strcutured sparse matrix. The important
thing is that a `dg::DMatrix` can be used together with a `dg::DVec`
in the `dg::blas2::symv` functions. These are level 1 functions that do
nothing but applying a given matrix to a vector and storing the result in
another vector. Just as the `dg::blas1` functions the `dg::blas2` functions
are templates that work for a variety of data types. In an MPI implementation
we would use a `dg::MDMatrix`.

## A note on boundary conditions
The matrices in the 'dg' library only know **homogeneous** boundary conditions.
For example when choosing `dg::DIR` as the boundary condition in `dg::create::dx`
the assumed boundary value is zero and when choosing `dg::NEU` the assumed
derivative on the boundary is also zero.
This has the advantage that our matrices are always linear. However, how do I
implement other boundary conditions then, you ask. For example when the
boundary value of my function is 1 and not 0?

The answer is that you'll have to manually subtract the value from the function
**before** you apply the derivative. For example if we assume that 'lhs' in the
previous example has non-homogeneous boundary conditions we would do something like

``` cpp
double function_with_dir_boundary( double x, double y){
    return 1.+sin(x)*sin(y);
}
dg::DVec fun = dg::evaluate( function_with_dir_boundary, grid), der(fun);
dg::blas1::plus( fun, -1.); // subtract boundary value before deriving
dg::blas2::symv( dx, fun, der);
```

This procedure can be done similarly for Neumann boundary conditions, where
you'll have to construct a function with your boundary conditions and subtract
it before applying the derivative.
