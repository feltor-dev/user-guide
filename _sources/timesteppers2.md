(sec:pdes)=
# Advanced Timesteppers

We now want to demonstrate how to use Feltor to solve partial differential equations.
We use the simple advection diffusion equation as a model equation
\begin{align}
    \frac{\partial \omega}{\partial t} &= -v\cdot \nabla\omega + D \Delta \omega \\
     -\Delta \phi &= \omega \\
     v_x &:= -\partial_y \phi \\
     v_y &:= \partial_x \phi
\end{align}

## Explicit stepping
As long as the diffusion coefficient is small enough to not influence the CFL condition we can compute everything explicitly. In Feltor we simply need a functor implementing the right hand side like we did in the
last chapter. We here repeat the building blocks that regard the timestepper:

```cpp
template<class Geometry, class Matrix, class Container>
struct Equations
{
    void operator()(double t, const Container& omega, Container& omegaDot)
    {
        // solve Poisson equation
        // implement advection term
        // implement diffusion term
    }
};
// Construct init condition
omega = myproject::initial_conditions(grid, js["init"] );

// Construct Equations
myproject::Equations<dg::x::CartesianGrid2d, dg::x::DMatrix,
    dg::x::DVec> rhs( grid, js);

// The timestepper
dg::Adaptive< dg::ERKStep< dg::x::DVec>> adapt(tableau, omega);

// The timeloop
dg::AdaptiveTimeloop<dg::x::DVec> timeloop( adapt, rhs,
                    dg::pid_control, dg::l2norm, rtol, atol);
for( unsigned u=1; u<=maxout; u++)
{

    timeloop.integrate( time, omega, u*deltaT, omega,
                          u < maxout ? dg::to::at_least : dg::to::exact);
    // ...
}
```
## Explicit advection - implicit diffusion
We want to split the PDE into two parts: $E(\phi, \omega) = -\vec v \cdot \nabla \omega$ and $I(\omega) = D\Delta\omega$ with $\phi = S(\omega) = \Delta^{-1}\omega$.
We intend to use a semi-implicit time integrator:
```cpp
template<class Geometry, class Matrix, class Container>
struct Explicit
{
    Explicit(...){}
    void operator()( double t, const Container& omega,
                    Container& k)
    {
        // Solve Phi = S(omega) with e.g. Multigrid
        // Compute k = E(phi,omega)
    }
    void implicit_part( double t, const Container& omega,
                       Container& k)
    {
        // Compute k = D Delta omega
    }
};
```
```{note}
There are several possibilities to partition the explicit, implicit and solver parts into functors (or lambdas in the main program). Generally, it is a good idea to keep all equation related functionality in one class. Since we cannot overload the `operator()` twice we reverted to a little trick, writing the `implicit_part` method that we later bind to the `operator()` of the Implicit class below.
```
```cpp
template<class Geometry, class Matrix, class Container>
struct Implicit
{
    Implicit( Explicit<...>& exp, ...): m_exp(exp) {...}
    void operator()( double t, const Container& omega, Container& k)
    {
        m_exp->implicit_part( t, omega, k);      
    }
```
```{note}
For the solve method we here chose a PCG solver since the Laplace is self-adjoint. Note, how we used  a small lambda wrapper to compute the implicit left hand side.

Typically, the solver would also write some information about its performance to `std::cout` so that a user is kept informed about the status of the integration.
```
```cpp
    void operator()( double alpha, double t, Container& omega, const Container& rhs)
    {
        auto wrapper = [=]( const auto& x, auto& y){
            // x - a I (x,t)
            operator()( t, x, y); // calls the above operator
            dg::blas1::axpby( 1., x, -alpha, y);
        };
        dg::blas1::copy( rhs, omega); // use rhs as initial guess
        unsigned number = m_pcg.solve( wrapper, omega, rhs, 1., m_weights, m_eps_time);
    }
    private:
    Explicit<Geometry,Matrix,Container>& m_exp;
    dg::PCG< Container> m_pcg;
    Container m_weights;
    value_type m_eps_time;
};

// Construct equations
Explicit<...> ex( ...);
Implicit<...> im(ex, ...);

// The timestepper
dg::Adaptive< dg::ARKStep< dg::x::DVec>> adapt(tableau, omega);

// The timeloop
dg::AdaptiveTimeloop<dg::x::DVec> timeloop( adapt, std::tie( ex, im, im),
                    dg::pid_control, dg::l2norm, rtol, atol);
```

## Implicit advection-diffusion solver
In order to solve the entire system implicitly we have to write both equations and Solvers.
To make it more clear let us reformulate the structure of equations that we have
\begin{align}
    \dot \omega &= I(\omega, \phi)\\
     0 &= R(\omega, \phi)
\end{align}
with
\begin{align}
I( \omega, \phi) &= -\vec v \cdot \nabla \omega + D\Delta\omega \\
R( \omega, \phi) &= \Delta \phi + \omega
\end{align}
To implement an implicit timestepper we need to solve the equation (the mass matrix is the identity)
\begin{align}
\begin{cases}
    \omega - \alpha I(\omega, \phi) &= \omega^*  \\
    \omega+\Delta\phi &= 0
\end{cases}
\end{align}
We can solve this equation for $\omega$ by first solving
\begin{align}
    -\Delta\phi - \alpha I(-\Delta\phi, \phi) = \omega^*
\end{align}
for $\phi$ and then using $\omega = -\Delta\phi$.
```{note}
The solution for $\phi$ is the same for both variants but the solution for $\omega$ is not (numerically). This can be seen by setting $\alpha=0$. Then in the first version $\omega=\omega^*$ but in the second version $\omega=-\Delta\phi = \omega^* + \mathcal O(\epsilon)$, depends on how well the equation is solved. For this reason the accuracy of the implicit solver should be well higher than the accuracy of the timestepper.
```

Since these equations are non-linear the solver needs to be a non-linear solver. Our idea is to use a multigrid FAS solver. For this solver we need to implement both the operator as well as its inverse on multiple grids.
```cpp
template<class Geometry, class Matrix, class Container>
struct Equations //corresponds to I
{
    void rhs( double t, const Container& omega, const Container& phi, Container& omegaDot)
    {
    // compute: omegaDot = I(omega, phi, t)
        // implement advection term
        // implement diffusion term
    }
    void compute_omega( const Container& phi, Container& omega)
    {
        // compute -Delta phi
        dg::blas2::symv( m_lapM, phi, omega);
    }
};
template<class Geometry, class Matrix, class Container>
struct Implicit
{
    Implicit(...)
    {
        // constructed nested grids
        ...
        for ( unsigned u=0; u<stages; u++)
        {
            // construct Equations
            m_eqs[u] = ...
            m_imp[u] = [&, m_omega = m_weights[u]]
                ( const auto& phi, auto& f) mutable
            {
                // omega - a I(omega, phi)
                m_eqs[u].compute_omega( phi, m_omega );
                m_eqs[u].rhs( m_time, m_omega, phi, f);
                dg::blas1::axpby( 1., m_omega, -m_alpha, f);
            }
            m_inv_imp[u] = [acc = dg::AndersonAcceleration( ...),
                           &I = m_imp[u],
                           &m_weights = m_weights[u]]
                ( const auto& omS, auto& phi) mutable
            {  
                // Solve Implicit( phi) = omS
                acc.solve( I, phi, omS, m_weights...);
            };
        }
    }

    void operator()(double t, double alpha, Container& omega, const Container& omS)
    {
        // Note how the operators and solvers hold these as references
        // so they know when they are updated
        m_time = t, m_alpha = alpha;

        // Solve the implicit equation
        nested_iterations( m_imp, m_phi, omS, m_inv_imp, m_nested);
        m_eqs[0].compute_omega( m_phi, omega);
    }
    private:
    // since we capture by reference we cannot easily copy this class
    Implicit( const Implicit&);
    Implicit( Implicit&&);
    double m_time, m_alpha;
    dg::NestedGrids< Geometry, Matrix, Container> m_nested;
    Container m_phi;
    std::vector<Container> m_weights;
    std::vector<Equations<Geometry, Matrix, Container>> m_eqs;
    std::vector< std::function<void( const Container&, Container&)>> m_imp, m_inv_imp;
};

// Construct init condition
omega = myproject::initial_conditions(grid, js["init"] );

// Construct Equations
myproject::Equations<dg::x::CartesianGrid2d, dg::x::DMatrix,
    dg::x::DVec> rhs( grid, js);

// The timestepper
dg::Adaptive< dg::DIRKStep< dg::x::DVec>> adapt(tableau, omega);

// The timeloop
dg::AdaptiveTimeloop<dg::x::DVec> timeloop( adapt, std::tie(imp,imp),
                    dg::pid_control, dg::l2norm, rtol, atol);
for( unsigned u=1; u<=maxout; u++)
{

    timeloop.integrate( time, omega, u*deltaT, omega,
                          u < maxout ? dg::to::at_least : dg::to::exact);
    // ...
}
```

## Implicit mass-matrix -  timestepper
 Consider the following general equation:
\begin{align}
    M(y,t)\frac{d y}{d t} &= F(y,t)
\end{align}
We intend to solve the resulting implicit equation with multigrid nested iteration. We thus have to implement and solve the following operators on multiply grids:
\begin{align}
     M(y,t) k &= F(y,t) \\
     M(y,t) ( y - y^* ) - \alpha F(y,t) &= 0
\end{align}
To simplify the implementation we restrict ourselves to solvers that
use the implicit solve only (i.e. Multistep, symplectic DIRK and SDIRK timesteppers).

```cpp
template<...>
struct Equations
{
    Equations(...){}
    void mass_matrix( double a, double t, const Container& y,
                     const Container& k,
                    double b, Container& result)
    {
        // Compute result = a M(y,t) k + b result
        //...
    }
    void rhs( double t, const Container& y,  Container& k)
    {
        // Compute k = F(y,t)
    }
};
template< ...>
class ImplicitSolver
{
    ImplicitSolver( ...){
        // Construct nested grids
        ...
        // Construct nested objects m_eqs
        ...
        // Construct nested Operators and Solvers:
        for( unsigend u=0; u<stages; u++)
        {
            m_imp[u] = [&ys = m_ys, &time = m_time, &alpha = m_alpha, k = m_weights[u]]
                ( const auto& y, auto& I) mutable
            {  
                //compute M(y,t) ( y- y^*) - a F(y,t)
                m_eqs[u].rhs(  time, y, I);
                dg::blas1::axpby( 1., y, -1., ys, k);
                m_eqs[u].mass_matrix(  1., time, y, k, -alpha,  I);
            };
            m_inv_imp[u] = [acc = dg::AndersonAcceleration( ...),
                                 &I = m_implicit[u],
                                 &m_weights = m_weights[u]]
                ( const auto& zero, auto& y) mutable
            {  
                // Solve Implicit( y) = 0
                acc.solve( I, y, zero, m_weights...);
            };
        }
    }
    void operator()( double t, const Container& y,  Container& k)
    {
        assert( false && "This should not be called");
    }
    void operator()( double alpha, double t, Container& y, const Container& ys)
    {
        // Note how the operators and solvers hold these as references
        // so they know when they are updated
        m_time = t, m_alpha = alpha;
        m_nested.project( ys, m_ys);
        // Solve the implicit equation
        nested_iterations( m_imp, y, 0., m_inv_imp, m_nested);
    }
    private:
    // since we capture by reference we cannot easily copy this class
    ImplicitSolver( const ImplicitSolver&);
    ImplicitSolver( ImplicitSolver&&);
    double m_time, m_alpha;
    dg::NestedGrids< ...> m_nested;
    std::vector<Container> m_weights, m_ys;
    std::vector<Equations<...>> m_eqs;
    std::vector< std::function<void( const Container&, Container&)>> m_imp, m_inv_imp;
};

// in main
ImplicitSolver<...> imp(...);
dg::ImplicitMultistep<Vector> multistep("BDF-3-3", y0);
multistep.init( std::tie( imp, imp), t0, y0, dt);
multistep.step( std::tie( imp, imp), t0, y0);
```
