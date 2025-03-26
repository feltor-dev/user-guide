(sec:timesteppers)=
# ODE integrators
This chapter deals with how to integrate differential equations of the form
\begin{align}
M(y,t) \frac{d}{dt}y = F(y,t)
\end{align}
where $M$ is the mass matrix and $F$ is the right hand side.

## Overview
Generally we need three things to solve an ODE numerically
 - **the stepper method** (there are several stepper classes each coming with a range of tableaus to chose from: Runge Kutta vs Multistep, explicit vs implicit and more) e.g. `dg::ExplicitMultistep<Vector>`
 - **the ode itself**, this has to be provided by the user as the right hand side functor or tuple of functors including the solver method if an implicit stepper is used and then tied to a single object typically using `std::tie`, e.g. `std::tie( ex, im, solve);`
 - **the timeloop** (adaptive timestep vs fixed timestep) -> can be implemented by the user or in Feltor represented by an instance of `dg::aTimeloop` class (useful if one wants to choose a method at runtime) e.g. `dg::AdaptiveTimeloop<Vector>`. The timeloop includes the **initial condition** and the integration boundaries (these are the paramters to the integrate function) `timeloop.integrate( t0, u0, t1, u1);`
 .

The `dg` library provides a wide selection of explicit, implicit and imex timesteppers:
- `dg::ERKStep` embedded Runge Kutta
- `dg::DIRKStep` diagonally implicit Runge Kutta
- `dg::ARKStep` additive Runge Kutta (imex)
- `dg::ExplicitMultistep`
- `dg::ImplicitMultistep`
- `dg::ImExMultistep`
-  see [doxygen documentation](https://mwiesenberger.github.io/feltor/dg/html/group__time.html) for a full list

The `dg::ERKStep`, `dg::DIRKStep` and `dg::ARKStep` (and in general any embedded method) can be used as a driver for an adaptive Timestepper
- `dg::Adaptive`
```{admonition} Vector type
All of the above classes are templates of the Vector type in use. Chosing an appropriate type for the integration variable(s) is usually the first decision to make when implementing a differential equation. Anything that works in the `dg::blas1` functions is allowed so check out the [Vectors](sec:vectors) chapter.
```


In order to implement any Runge-Kutta or multistep algorithm we need to be able to solve two general types of equations (see [theory guide on overleaf](https://www.overleaf.com/read/dfxncmnnpzfm)):
\begin{align}
    k =  M(y,t)^{-1} \cdot F(y,t) &\text{ given $t,y$ return $k$} \\
    M(y,t)\cdot ( y-y^*) - \alpha F(y,t) = 0 &\text{ given $\alpha, t, y^*$ return $y$}
\end{align}
For any explicit (part of the) equation we need to solve the first and for any implicit (part of the) equation  we need to solve both the first and the second equation.
````{note}
It is important to realize that the timestepper needs to know neither $F$ nor $M$, nor how they are implemented or how the implicit equation is solved.
The user just needs to provide oblique functor
objects with the signature
```cpp
// given t y write ydot
void operator()( value_type t, const ContainerType0& y, ContainerType1& ydot);
// given \alpha, t and y^* write a new y
void operator()( value_type alpha, value_type t, ContainerType& y, const ContainerType1& ystar);
```
````
More than one functor is needed for implicit or semi-implicit timesteppers. These are then expected as a `std::tuple` of functors.
Typically one would use `std::tie` to tie together two or more objects into a functor.

The timesteppers themselves are implemented using `dg::blas1` vector additions and all have a `step` method that advances the ode for one timestep
```cpp
dirk.step( std::tie( implicit_part, solver), t0, y0, t1, y1, dt);
```
This step method can be called in a loop to advance the ode for an arbitrary time. Or you can use one of the `dg::aTimeloop` classes
- `dg::AdaptiveTimeloop`
- `dg::SinglestepTimeloop`
- `dg::MultistepTimeloop`

These classes abstract the integration loop and can be called via a common interface
```cpp
timeloop.integrate( t0, y0, t1, y1);
```
This is useful especially if you want to choose various timesteppers at runtime.


We will now study a few case scenarios to clarify the above explanations.

## Integrate ODEs in Feltor - first timesteppers and simple time loops

As an example we solve the damped driven harmonic oscillator
\begin{align}
    \frac{d x}{d t} &= v \\
    \frac{d v}{d t} &= -2 \nu \omega_0 v - \omega_0^2 x + \sin (\omega_d t)
\end{align}
Since we have two scalar variables we will use a `std::array<double,2>` as the vector type to use.
Let us first choose some (somewhat random) parameters

```cpp
const double damping = 0.2, omega_0 = 1.0, omega_drive = 0.9;
```
We know that we can solve this ODE analytically. This comes in handy to verify our implementation.

```cpp
// We have an analytical solution
std::array<double,2> solution( double t)
{
    double tmp1 = (2.*omega_0*damping);
    double tmp2 = (omega_0*omega_0 - omega_drive*omega_drive)/omega_drive;
    double amp = 1./sqrt( tmp1*tmp1 + tmp2*tmp2);
    double phi = atan( 2.*omega_drive*omega_0*damping/(omega_drive*omega_drive-omega_0*omega_0));
    double x = amp*sin(omega_drive*t+phi)/omega_drive;
    double v = amp*cos(omega_drive*t+phi);
    return {x,v};
}
```
### Explicit Runge-Kutta - fixed step

First, show how to implement a simple timeloop with a fixed stepsize Runge Kutta integrator. We choose the classic 4-th order scheme, but consult the [documentation](https://mwiesenberger.github.io/feltor/dg/html/group__time.html) for an extensive list of available tableaus.

```cpp
// The right hand side needs to be a callable function in Feltor.
// In modern C++ this can for example be a lambda function:
auto rhs = [&]( double t, const std::array<double,2>& y,
            std::array<double,2>& yp)
{
    //damped driven harmonic oscillator
    // x -> y[0] , v -> y[1]
    yp[0] = y[1];
    yp[1] = -2.*damping*omega_0*y[1] - omega_0*omega_0*y[0]
            + sin(omega_drive*t);
};
// Let us choose an initial condition and the integration boundaries
double t0 = 0., t1 = 1.;
const std::array<double,2> u0 = solution(t0);

// Here, we choose the classic Runge-Kutta scheme to solve
dg::RungeKutta<std::array<double,2>> rk("Runge-Kutta-4-4", u0);
// Now we are ready to construct a time-loop by repeatedly stepping
// the Runge Kutta solve with a constant timestep
double t = t0;
std::array<double,2> u1( u0);
unsigned N = 20;
for( unsigned i=0; i<N; i++)
    rk.step( rhs, t, u1, t, u1, (t1-t0)/(double)N );

// Now let us compute the error
const std::array<double,2> sol = solution(t1);
dg::blas1::axpby( 1., sol , -1., u1);
std::cout << "Norm of error is " <<sqrt(dg::blas1::dot( u1, u1))<<"\n";
//Norm of error is 8.17315e-08
```

```{note}
Notice how the above program consists of the three mentioned ingredients. We wrote the ode as a single functor `rhs` expecting to use the explicit single step Runge Kutta method `dg::RungeKutta`. We then wrote a timeloop starting from the solution at time `t0`.
```
```{seealso}
Lambdas are so useful that they appear at several places in the user-guide and the `dg` library generally accepts them in many places. If you are not familiar with them, a good start is to watch youtube! For example [Lambdas from Scratch - Arthur O'Dwye](https://www.youtube.com/watch?v=3jCOwajNch0)
```

### Implicit Multistep  - fixed step
In the next example we want to solve the same ode with an implicit multistep method:

```cpp
// First we need to provide a solution method for the prototypical implicit equation
auto solve = [&]( double alpha, double t, std::array<double,2>& y,
            const std::array<double,2>& yp)
{
    // y - alpha RHS( t, y) = rho
    // can be solved analytically
    y[1] = ( yp[1] + alpha*sin(omega_drive*t) - alpha*omega_0*omega_0*yp[0])/
           (1.+2.*alpha*damping*omega_0+alpha*alpha*omega_0*omega_0);
    y[0] = yp[0] + alpha*y[1];
};
// Now we can construct a multistep method
dg::ImplicitMultistep<std::array<double,2>> multi("ImEx-BDF-3-3", u0);
// Let us choose the same initial conditions as before
t = t0; u1 = u0;
// Finally, we can construct a timeloop
multi.init( std::tie( rhs, solve), t, u0, (t1-t0)/(double)N);
for( unsigned i=0; i<N; i++)
    multi.step( std::tie(rhs, solve), t, u1);

dg::blas1::axpby( 1., sol , -1., u1);
std::cout << "Norm of error is " <<sqrt(dg::blas1::dot( u1, u1))<<"\n";
//dg::make_odeint( dirk, std::tie( rhs, solve), (t1-t0)/20.)->integrate( t0, u0, t1, u1);
//Norm of error is 3.60666e-05
```

```{note}
Notice the use of `std::tie` in the `init` and `step` methods. Since an implicit method needs both the implicit part `rhs` and an implicit solver `solve` to solve it.
```

### Embedded explicit Runge Kutta - adaptive step
As a last example we want to integrate the ode with an adaptive timestepper
```cpp
dg::Adaptive<dg::ERKStep<std::array<double,2>>> adapt("Tsitouras11-7-4-5", u0);
t = t0; u1 = u0;
double dt = 1e-6;
while( t < t1)
{
    if( t + dt > t1)
        dt = t1 - t;
    adapt.step( rhs, t, u1, t, u1, dt, dg::pid_control, dg::fast_l2norm, 1e-6, 1e-6);
}

dg::blas1::axpby( 1., sol , -1., u1);
std::cout << "Norm of error is " <<sqrt(dg::blas1::dot( u1, u1))
    <<" with "<<adapt.nsteps()<<" steps \n";
// Norm of error is 2.66678e-07 with 8 steps
```

## Abstract timeloops with dg::aTimeloop
We have seen that the main difference between the three previous methods was how to construct the timeloop. Runge-Kutta was just a basic for loop, Multistep needed to be initialized before usage, while the adaptive stepper needed to run in a while loop and needed a bunch of additional parameters. We will now introduce
an abstract interface that lets you choose at runtime which integration method to use:
```cpp
using Vec = std::array<double,2>;
const unsigned* nsteps;

std::string stepper = "Adaptive";
std::cin >> stepper;

// A pointer to an abstract integrator
auto odeint = std::unique_ptr<dg::aTimeloop<Vec>>();

if( "Adaptive" == stepper)
{
    dg::Adaptive<dg::ERKStep<Vec>> adapt("Tsitouras11-7-4-5", u0);
    odeint = std::make_unique<dg::AdaptiveTimeloop<Vec>>( adapt,
                rhs, dg::pid_control, dg::fast_l2norm, 1e-6, 1e-6);
    nsteps = &adapt.nsteps();
}
else if( "Singlestep" == stepper)
{
    dg::RungeKutta<Vec> rk("Runge-Kutta-4-4", u0);
    odeint = std::make_unique<dg::SinglestepTimeloop<Vec>>( rk,
                                         rhs, (t1-t0)/(double)N);
    nsteps = &N;
}
else if( "Multistep" == stepper)
{
    dg::ExplicitMultistep<Vec> multi( "ImEx-BDF-3-3", u0);
    odeint = std::make_unique<dg::MultistepTimeloop<Vec>>( multi,
                        rhs, t0, u0, (t1-t0)/(double)N);
    nsteps = &N;
}
else
    std::cerr << "Stepper "<<stepper<<" not known!\n";

odeint -> integrate(t0, u0, t1, u1);

dg::blas1::axpby( 1., sol , -1., u1);
std::cout << "Norm of error is " <<sqrt(dg::blas1::dot( u1, u1))
    <<" with "<<*nsteps<<" steps \n";
```
In order to build more complicated examples we will first need to learn how geometry and derivatives work and how we can build advanced solvers. We will revisit the timesteppers in [Section PDEs](sec:pdes).

## Integrating DAEs - solver only implicit integrators
Consider the van der Pol oscillator

\begin{align}
    \dot x &= v \\
    \dot v &= \sigma ( 1-x^2) v -x
 \end{align}
 where $\sigma$ represents the stiffness parameter. For
very large values these equations become **differential-algebraic equations of index 1** (which is easily seen if one divides by sigma and sets $1/\sigma = 0$, the remaining equation does not have derivatives, hence algebraic, index 1 means it is solvable).
It is then a singular perturbation problem.
This behaviour is difficult for time-integrators insofar it reduces the error of the method to its so-called _stage-order_.

Furthermore, for higher $\sigma$ it becomes paramount to use implicit integrators and sometimes it might be desirable to completely avoid any explicit stage altogether (recall that for example EDIRK schems have an explicit first stage). For this reason the doxygen documentation highlights implicit integrators of that sort:
- all of the implicit multistep methods avoid calling the rhs directly, however, currently only the `dg::BDF_X_X` schemes can avoid it in the initialization phase as well
- of the DIRK schemes only the SDIRK and symplectic DIRK schemes do not call the rhs.

This means for such schemes we only need to implement a single functor, namely one that solves
\begin{align}
    M(y,t)\cdot ( y-y^*) - \alpha F(y,t) = 0 &\text{ given $\alpha, t, y^*$ return $y$} 
\end{align}

Here we outline the structure of a possible implementation
```cpp

struct VanderPolSolver
{
     void operator()( double alpha, double t, std::array<double,2>& y, const std::array<double,2>& yp)
     {
         // solve y - alpha F(y,t) = yp
     }
}

auto doNotCall = [](double t, const auto& x, auto& y)
{
    // Just to be entirely certain
    assert(false && "This should never be called!");
}

VanderPolSolver solver( ...);
dg::DIRKStep<Vec> dirk( "SDIRK-4-2-3", y0);

odeint = std::make_unique<dg::SinglestepTimeloop<Vec>>( dirk,
                                                       std::tie( doNotCall, solver),
                                                       (t1-t0)/(double)N);
odeint -> integrate(t0, u0, t1, u1);
```
