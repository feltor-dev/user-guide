# Applications
## Putting it all together

What remains to do is to
 put it all together and solve a "real life" problem.
 Let us look at the following code. What it does is to solve
 the heat diffusion equation with source term using a semi-implicit
 multistep scheme in time. The diffusive part is solved implicitly
 and the source term is contained in the explicit part. We
 manufactured a solution to test the implementation:

```{code-block} cpp
#include <iostream>

#include "dg/algorithm.h"

//method of manufactured solution
struct Solution{
  Solution(double t, double nu):t(t), nu(nu){}
  DG_DEVICE
     double operator()(double x, double y) const{
         return sin(t) * exp( -2*t*nu) * sin(x) * sin(y);
     }
  private:
  double t, nu;
};

struct Source{
  Source(double t, double nu):t(t), nu(nu){}
  DG_DEVICE
     double operator()(double x, double y) const{
         return sin(x) * sin(y) * cos(t) * exp(-2*t*nu) * (1-sin(t));
     }
  private:
  double t, nu;
};

//the explicit part contains the source Tp = S(x,y,t)
template<class container>
struct Explicit
{
  Explicit( const dg::Grid2d& g, double nu):
         m_nu( nu),
         m_x ( dg::evaluate(dg::cooX2d, g)),//x-coordinate
         m_y ( dg::evaluate(dg::cooY2d, g)) //y-coordinate
     {}
  void operator()( double t, const container& T, container& Tp) {
    dg::blas1::evaluate( Tp, dg::equals(), Source(t,m_nu), m_x, m_y);
  }
private:
  const double m_nu;
  const container m_x, m_y;

};

//the implicit part contains  Tp = nu Delta T(x,y,t) + cos(t) T(x,y,t)
template< class Matrix, class container>
struct Implicit
{
  Implicit( const dg::Grid2d& g, double nu):
         m_nu(nu),
         m_w2d( dg::create::weights(g)),
         m_v2d( dg::create::inv_weights(g)),
         m_LaplacianM( g, dg::normed)
         { }

  void operator()( double t, const container& T, container& Tp)
  {
    dg::blas2::gemv( m_LaplacianM, T, Tp);
    dg::blas1::axpby( cos(t), T, -m_nu, Tp);
  }
  //required by inversion in semi-implicit schemes
  const container& inv_weights(){return m_v2d;}
  const container& weights(){return m_w2d;}
  const container& precond(){return m_v2d;}
private:
  double m_nu;
  const container m_w2d, m_v2d, m_x, m_y;
  dg::Elliptic<dg::CartesianGrid2d, Matrix, container> m_LaplacianM;
};

const double lx = 2.*M_PI;
const double ly = 2.*M_PI;

int main()
{
  unsigned n = 3, Nx = 50 , Ny = 50;
  const double T = 0.1;
  const double NT= 40, eps = 1e-8;
  const double dt = (T/NT);
  const double nu = 0.01;
  //construct the grid and the explicit and implicit parts
  dg::Grid2d grid( 0, lx, 0, ly, n, Nx, Ny, dg::PER, dg::PER);
  Explicit<dg::DVec> ex( grid, nu);
  Implicit<dg::DMatrix, dg::DVec> im( grid, nu);

  //evaluate the initial condition
  const dg::DVec init( dg::evaluate(Solution(0.,nu), grid));
  dg::DVec y0(init);

  const dg::DVec sol = dg::evaluate( Solution(T,nu), grid);
  const dg::DVec w2d = dg::create::weights( grid);
  const double norm_sol = dg::blas2::dot( w2d, sol);
  double time = 0.;
  dg::DVec error( sol);

  //construct time stepper
  dg::Karniadakis< dg::DVec > karniadakis( y0, y0.size(), eps);
  time = 0., y0 = init;
  //initialize the timestepper
  karniadakis.init( std::tie(ex, im), time, y0, dt);
  //main time loop (NT = 20)
  for( unsigned i=0; i<NT; i++)
    karniadakis.step( std::tie(ex, im), time, y0); //inplace step

  dg::blas1::axpby( -1., sol, 1., y0);
  double res = sqrt(dg::blas2::dot( w2d, y0)/norm_sol);
  std::cout << "Relative error Karniadakis is "<< res<<std::endl;

  return 0;
}
```

We have now mastered the basic `dg` library. We hope it became clear
that it provides building blocks that you can use to solve your
very own problem. So how do I do Input and Output of fields in my code,
you might ask. The answer is that we do not do this. For input/output
you will have to use other libraries that specialize on these operations.
If we may recommend something,
then use JSON as a format for input files and NetCDF-4 for output files.
Please take a look at the `feltor/src/` folders to see examples of how
we ourseleves use our library in connection with JSON and NetCDF for a variety of physical models. 
