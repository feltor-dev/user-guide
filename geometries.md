# Curvilinear grids

## The Geometries extension
We provide an extension to the standard dg library in
`dg/geometries/geometries.h` with the goal to provide more versatile grids and
geometries and especially the utilities to model toroidal magnetic field geometries.
Logically, this extension enhances Levels 3 and 4 of the original dg library.

The following code demonstrates how to generate a structured orthogonal grid
between two contour lines of a two-dimensional flux-function that models the
magnetic flux in a magnetic confinement fusion device. We then use the
volume form of this grid to compute the area integral of a function on the
domain.

```{code-block} cpp
#include <iostream>
//include the dg-library
#include "dg/algorithm.h"
//include the geometries extension (introduces the dg::geo namespace)
#include "dg/geometries/geometries.h"

// define a function to integrate
double function( double x, double y){
    return x;
}

// parameters for the flux function
dg::geo::solovev::Parameters p;
p.A = 1.;
p.c = { 0.104225846989823495731985393840,
       -0.109784866431799048073838307831,
        0.172906442397641706343124761278,
       -0.0325732004904670117781012296259,
        0.0100841711884676126632200269819,
       -0.00334269397734777931041081168513,
       -0.000108019045920744348891466483526,
        0., 0., 0., 0., 0., 1.};
p.R_0 = 547.891714877869;
// indicate the two contour lines that bound the domain
double psi_0 = -20., psi_1 = -4.;
// create the flux function and its derivatives (up to 2nd order)
dg::geo::CylindricalFunctorsLvl2 psip = dg::geo::solovev::createPsip( p);
// create a grid generator for a simple orthogonal grid
dg::geo::SimpleOrthogonal generator( psip, psi_0, psi_1, p.R_0, 0., 0);
// create a grid with the help of the grid generator
unsigned n = 3, Nx = 8, Ny = 80;
dg::geo::CurvilinearGrid2d g2d( generator, n, Nx, Ny, dg::NEU);
// get an element of the metric tensor
dg::HVec g_xy = g2d.metric().value(0,1);
// let's check if the xy element is really 0
double zero = dg::blas1::dot( g_xy, g_xy);
// lo and behold
std::cout << "Norm of off-diagonal element " << zero << "\n";
// create the volume form
dg::HVec vol2d = dg::create::volume( g2d);
// pull back the function to the grid coordinates
dg::HVec f = dg::pullback( function, g2d);
// compute the area integral on the grid domain
double integral = dg::blas1::dot( vol2d, f);
// rather large ... ( 6.122e7)
std::cout << "The integral is " << integral << "\n";
```

The first part of the code just defines the parameters of our flux function and
is very specific to the type of function we use here. In an application code we
will typically write those parameters into an input file and then let the
program read in the parameters. With the help of the function and its derivatives
we can then construct a generator. This generator can then finally be used to
construct a grid. In fact, it constructs the complete coordinate transformation
including the Jacobian and metric tensors. Note that the whole process is
fully customizable, meaning you can create your own functions, create your
own grid generator or even your own grid by deriving from the relevant
abstract base class. The documentation will provide further detail on this.

The important thing to notice here is that the `dg::geo::CurvilinearGrid2d` is of
course fully compatible with the dg library, which means
we can evaluate functions, create derivatives and solve equations on the new grid
in the same way we did in the previous chapters.
In the example code above we just use some basic operations to demonstrate this point.

One thing to look out for is to use `dg::pullback` instead of `dg::evaluate` to
discretize functions on the grid. The distinction is the coordinate system
in which we define our function. In this case we define `function` in
Cartesian coordinates rather than in the transformed computational coordinates.
As the name suggests `dg::pullback` pulls the function from Cartesian back to
computational coordinates. In the example above `f` represents `x(zeta, eta)`, where `zeta` and `eta` are the computational space coordinates. See our
publication
 [Streamline integration as a method for two-dimensional elliptic grid generation](https://doi.org/10.1016/j.jcp.2017.03.056) for more details on this point and also
 on the grid generation method used.
