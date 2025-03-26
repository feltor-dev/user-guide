(sec:project)=
# Create your own project

In this chapter we are going to learn how to write a code project implementing a set of equations from beginning to end including file I/O and diagnostics.
```{note}
In principle you can write your own codes however you want. Feltor does not force any specific structure onto you or force you to use more of it than you want to. This chapter is a reflection of best practices that we currently use for our own projects and that have proven themselves in scientific analysis. It is for example specifically designed to interoperate well with python. This design has evolved through the last 10 years and, for example when new C++ features become available, will likely continue to evolve in the future.
```

We use the simple advection diffusion equation as a model equation
\begin{align}
    \frac{\partial \omega}{\partial t} &= -v\cdot \nabla\omega + D \Delta \omega \\
     -\Delta \phi &= \omega \\
     v_x &:= -\partial_y \phi \\
     v_y &:= \partial_x \phi
\end{align}

## Basic file structure


- `init.h` Contains a list of initial conditions that can be chosen in the inputfile
- `diag.h` Contains a list of diagnostic functions that should be executed periodically in the main timeloop. This will determine the content of the output file
- `Makefile` How to compile the main program
- `parameters.h` [optional] Reads and stores parameters from the input file. The intention here is to provide quick (in terms of code length) access to often used parameters.
- `equations.h` Contains the main functor that is used in the timestepper. Implements the equations.
- `documentation.tex` Contains the documentation of the implemented equations, the numerical methods used, the initial conditions and the diagnostics. In particular it describes the parameters of the input file and what the output file contains.
- `input.json` The input file to a simulation
- `myprogram.cpp` Here is the main function. It reads in the input file, calls the relevant initial condition, constructs the Equations functor(s) and the Timestepper. Then, it implements the main time loop, the diagnostics and output strategy.

```{admonition} I/O file formats
We use and highly recommend [json](https://en.wikipedia.org/wiki/JSON) for the input file and [NetCDF-4](https://www.unidata.ucar.edu/software/netcdf/) for the output file.
The reason for this is that both of these formats interoperate well with many other programming languages, in particular python, which we use for data analysis.
```

## The init file
The goal of the init file is to provide various initial conditions that can be selected in the input file
```{admonition} Data layout
The return type of the initial conditions depends on what data layout you choose when implementing [The equations file](sec:equations). In this case it is `dg::x::DVec` because there is only one variable, but often used types are `std::vector<dg::x::DVec>`, `std::map<std::string, dg::x::DVec>` or `std::array<dg::x::DVec, 2>`
```

```cpp
// init.h
#pragma once

namespace myproject
{

dg::x::DVec initial_conditions(
    const dg::x::CartesianGrid2d& grid,
    dg::file::WrappedJsonValue init)
{
    dg::x::HVec omega;
    std::string initial = init.get("type", "lamb").asString();
    if( initial == "zero")
    {
        omega = dg::evaluate( dg::zero, grid);
    }
    else if( initial == "lamb")
    {
        double posX = init.get("posX", 0.5).asDouble();
        double posY = init.get("posY", 0.8).asDouble();
        double R =    init.get("sigma", 0.1).asDouble();
        double U =    init.get("velocity", 1).asDouble();
        dg::Lamb lamb( posX*grid.lx(), posY*grid.ly(), R, U);
        omega = dg::evaluate ( lamb, grid);
    }
    // ... implement more initial conditions here
    else
        throw dg::Error( dg::Message() << "Initial condition "
                        <<initial<<" not recognized!");
    return dg::construct<dg::x::DVec>(omega);
}

} //namespace myproject
```

The above example would require the following entry in the input file
```json
{
    "init" :
    {
        "type" : "lamb",
        "posX" : 0.5,
        "posY" : 0.8,
        "sigma" : 0.1,
        "velocity" : 1
    }
}
```
If the zero initial condition is chosen, the entry has to look like
```json
{
    "init" :
    {
        "type" : "zero"
    }
}
```
As you can see only the necessary parameters need to be present for each type of initial condition.
```{note}
You will have to document each initial condition and what parameters can be chosen in the documentation.tex file Here the "minted" package is highly recommended to document json parameters.
```

## The diag file

The desired list of variables in the output file typically grows (or sometimes shrinks) often during the lifetime of a project for example when new theoretical results are revealed or a new data-analysis needs to be performed.
The goal of our design here is therefore to provide a flexible and importantly an **extensible** approach to generating file output.

The purpose of the diag file is to provide one or more list of functions that generate output variables in the netcdf output file.
Each item in the list is a record consisting of a name, a description and a function that generates the data points of the variable. The functions' parameters are bundled in a struct `Variables`.
It is important to note that the list of records (the diagnostics\*\_list ) can be easily extended without side-effects of the existing records. The Variables can be extended, but the Variable construction line in the main file also needs to extend then (see [The main file](sec:main_file_construction) )

```cpp
// diag.h
#pragma once
#include "equations.h"

namespace myproject
{

struct Variables
{
    Equations<dg::x::CartesianGrid2d, dg::x::DMatrix, dg::x::DVec>& rhs;
    const dg::x::CartesianGrid2d& grid;
    const dg::x::DVec& omega;
};

struct Record
{
    std::string name; // variable name in the output file
    std::string long_name; // longer description as an attribute
    std::function<void(dg::x::DVec&, Variables&)> function;
    // function that generates the data points for the variable
};

// time - independent output (only called once)
std::vector<Record> diagnostics2d_static_list = {
    { "xc", "x-coordinate in Cartesian coordinate system",
        []( dg::x::DVec& result, Variables& v ) {
            result = dg::evaluate( dg::cooX2d, v.grid);
        }
    },
    { "yc", "y-coordinate in Cartesian coordinate system",
        []( dg::x::DVec& result, Variables& v ) {
            result = dg::evaluate( dg::cooY2d, v.grid);
        }
    },
    { "weights", "Gaussian integration weights",
        []( dg::x::DVec& result, Variables& v ) {
            result = dg::create::weights( v.grid);
        }
    },
    // ... extend here
};

// time - dependent output (called periodically)
std::vector<Record> diagnostics2d_list = {
    {"vorticity", "Vorticity in 2d",
        []( dg::x::DVec& result, Variables& v ) {
             dg::blas1::copy(v.omega, result);
        }
    },
    {"potential", "stream function",
        []( dg::x::DVec& result, Variables& v ) {
             dg::blas1::copy(v.rhs.potential(), result);
        }
    },
```
```{note}
You can write anything you want in the lambda, really.
At least anything you want to happen periodically in the main time loop.
```
```cpp
    {"enstrophy", "Squared vorticity",
        []( dg::x::DVec& result, Variables& v ) {
             // more complicated algorithms can be written here
             dg::blas1::pointwiseDot( v.omega, v.omega, result);
             dg::blas1::scal( result, 1./2.);
        }
    }
    // ... extend here
};

// ... write more lists, for example 1d diagnostics etc
// each list requires a loop in the main program


} //namespace myproject
```

(sec:equations)=
## The equations file

Our timesteppers require functors that implement the right hand side of the differential equation.
The purpose of the equations file is to provide such functor(s) for the specific equations at hand.

```cpp
// equations.h

#pragma once

#include "dg/algorithm.h"

namespace myproject
{

template<class Geometry, class Matrix, class Container>
struct Equations
{
    Equations( const Geometry& grid, dg::file::WrappedJsonValue& js):
        m_phi ( dg::evaluate(dg::zero, grid)),
        m_old_phi( 2, m_phi),
        m_adv( grid), // use grid's boundary conditions
        m_lapM( grid)
    {
        m_v = {m_phi, m_phi};
        m_pcg.construct(m_phi, grid.size());
        m_eps_pol = js["elliptic"].get("eps_pol", 1e-6).asDouble();
        m_nu = js["physical"].get( "nu", 1e-7).asDouble();
        m_centered[0] = dg::create::dx( grid, grid.bcx(), dg::centered);
        m_centered[1] = dg::create::dy( grid, grid.bcy(), dg::centered);
    }
    // accessors for diag.h
    const Container& potential() const {return m_phi;}
```
\begin{align}
    \frac{\partial \omega}{\partial t} &= -v\cdot \nabla\omega + D \Delta \omega \\
     -\Delta \phi &= \omega \\
     v_x &:= -\partial_y \phi \\
     v_y &:= \partial_x \phi
\end{align}
```{admonition} Data layout
Here, in the parameters for `operator()` you have to choose a suitable data-layout for the vector that is integrated in time. We chose a simple `Container` because there is only one dynamic variable in the equations. However, if there were three equations say then `std::array<Container,3>` could be a better type or even `std::vector<Container>` if you want to choose at runtime or the recently added `std::map<std::string, Container>` for a verbose access
```

```cpp
    // We implement advection diffusion equations:
    void operator()(double t, const Container& omega, Container& omegaDot)
    {

        // Solve potential equation
        m_old_phi.extrapolate( t, m_phi);
        // For demonstration we here use the simple unpreconditioned PCG
        // solver. In real code it is highly recommended to use nested_iterations instead
        // See "Solvers" chapter of this guide
        m_pcg.solve( m_lapM, m_phi, omega, 1., m_lapM.weights(), m_eps_pol);
        m_old_phi.update( t, m_phi);

        // add advection term
        dg::blas2::symv( -1., m_centered[1], m_phi, 0., m_v[0]);
        dg::blas2::symv( +1., m_centered[0], m_phi, 0., m_v[1]);
        m_adv.upwind( -1., m_v[0], m_v[1], omega, 0., omegaDot);

        // add diffusion
        dg::blas2::symv( -m_nu, m_lapM, omega, 1., omegaDot);
    }

    private:
    Container m_phi;
    dg::Extrapolation<Container> m_old_phi;
    dg::Advection<Geometry, Matrix, Container> m_adv;
    dg::Elliptic< Geometry, Matrix, Container> m_lapM;
    dg::PCG<Container> m_pcg;
    std::array<Matrix,2> m_centered;
    std::array<Container,2> m_v;
    double m_eps_pol;
    double m_nu;
};

}//namespace myproject
```
This necessitates two more parameters in the input file
```json
{
    "physical":
    {
        "nu" : 1e-7
    },
    "elliptic":
    {
        "eps_pol" : 1e-6
    }

}
```
We here showcase the simplest case where we implement everything explicitly (hence only one functor) and use a simple unpreconditioned PCG solver for the elliptic equation. This is just to show the basic structure. You can extend this however you want and make it arbitrarily complicated.

## The main file

The main function acts as a driver where everything comes together and has to be called in the right order.

### Front matter
We here decide to write a program that works for both shared and distributed memory. In order for this to work we need to initialize MPI, but only in case we actually need it. Typically, we would do the following:

```cpp
// myprogram.cpp
#include <iostream>
#include <iomanip>

#ifdef WITH_MPI
#include <mpi.h>
#endif //WITH_MPI

#include "dg/algorithm.h"
#include "dg/file/file.h"

#include "equations.h"
#include "init.h"
#include "diag.h"

int main( int argc, char* argv[])
{
#ifdef WITH_MPI
    dg::mpi_init( argc, argv);
    MPI_Comm comm;
    dg::mpi_init2d( dg::DIR, dg::PER, comm, std::cin, true);
    int rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
#endif //WITH_MPI
```
Here, we use a Macro `WITH_MPI` that is intended to be defined at compile-time.
```{note}
We therefore add `-DWITH_MPI` to the mpi compiler flags in [the Makefile](sec:Makefile)
```
Then in the main program we need to initialize MPI. The `dg::mpi_init2d` function generates a Cartesian Communicator that we will need to generate a grid and also will generate a prompt on `std::cin` asking the user to provide the partition of MPI processes.
```{admonition} MPI partitioning
We have not found a good way yet to automate the process of partitioning MPI processes amoung the dimensions. Therefore the user will be responsible for providing this information. For example if we have 8 processes overall we can provide `npx = 2`, `npy= 4` or `npx = 1`, `npy = 8` and so on.
```
Note that we here also need to distinguish between periodic and non-periodic boundary conditions in each direction (This influences the way communication needs to be set up).

### Read in the input file
```cpp
dg::file::WrappedJsonValue js( dg::file::error::is_throw);
if( argc != 3 )
{
    DG_RANK0 std::cerr << "ERROR: Wrong number of arguments!\nUsage: "
            << argv[0]<<" [input.json] [output.nc]\n \n"<<std::endl;
    dg::abort_program();
}
try{
    js = dg::file::file2Json( argv[1]);
} catch( std::exception& e) {
    DG_RANK0 std::cerr << "ERROR in input file "<<argv[1]<<std::endl;
    DG_RANK0 std::cerr << e.what()<<std::endl;
    dg::abort_program();
}
DG_RANK0 std::cout << js.toStyledString() << std::endl;
```
Here, we notice the `DG_RANK0` macro, which expands to `if(rank==0)` when compiled with MPI and else stays empty. Note that all processes read in the input file.
(sec:main_file_construction)=
### Construction phase
Now, we have to construct the essential objects for time integration: Grid, vector with initial condition, Equations, Timestepper
```cpp
// Construct grid
unsigned n, Nx, Ny;
double x0, x1, y0, y1;
dg::bc bcx, bcy;
try{
    dg::file::WrappedJsonValue grid = js["grid"];
    n = grid.get( "n", 3).asUInt();
    Nx = grid.get( "Nx", 48).asUInt();
    Ny = grid.get( "Ny", 48).asUInt();
    x0 = grid["x"].get( 0u, 0.).asDouble();
    x1 = grid["x"].get( 1u, 1.).asDouble();
    y0 = grid["y"].get( 0u, 0.).asDouble();
    y1 = grid["y"].get( 1u, 1.).asDouble();
    bcx = dg::str2bc ( grid["bc"].get( 0u, "DIR").asString());
    bcy = dg::str2bc ( grid["bc"].get( 1u, "PER").asString());
}catch ( std::exception& error){
    DG_RANK0 std::cerr << "Error in input file " << argv[1]<< std::endl;
    DG_RANK0 std::cerr << error.what() << std::endl;
    dg::abort_program();
}
dg::x::CartesianGrid2d grid( x0, x1, y0, y1 , n, Nx, Ny, bcx, bcy
    // The MPI version of CartesianGrid2d needs a communicator:
    #ifdef WITH_MPI
    , comm
    #endif //WITH_MPI
    );

// Construct Equations
myproject::Equations<dg::x::CartesianGrid2d, dg::x::DMatrix,
    dg::x::DVec> rhs( grid, js);

// Construct initial condition
dg::x::DVec omega;
try{

    omega = myproject::initial_conditions(grid, js["init"] );
}catch ( std::exception& error){
    DG_RANK0 std::cerr << "Error in input file " << argv[1]<< std::endl;
    DG_RANK0 std::cerr << error.what() << std::endl;
    dg::abort_program();
}

// Construct timestepper
std::string tableau;
double rtol, atol, time = 0.;
try{
    rtol = js["timestepper"].get("rtol", 1e-5).asDouble();
    atol = js["timestepper"].get("atol", 1e-5).asDouble();
    tableau = js[ "timestepper"].get( "tableau", "Bogacki-Shampine-4-2-3").asString();
}catch ( std::exception& error){
    DG_RANK0 std::cerr << "Error in input file " << argv[1]<< std::endl;
    DG_RANK0 std::cerr << error.what() << std::endl;
    dg::abort_program();
}
dg::Adaptive< dg::ERKStep< dg::x::DVec>> adapt(tableau, omega);
dg::AdaptiveTimeloop<dg::x::DVec> timeloop( adapt, rhs,
                    dg::pid_control, dg::l2norm, rtol, atol);

```

In this part we notice that we put more requirements on the input file.
Specifically
```json
{
    "grid":
    {
        "n" : 3,
        "Nx" : 20,
        "Ny" : 20,
        "x" : [0,1],
        "y" : [0,1],
        "bc" : ["DIR", "PER"]
    },
    "timestepper":
    {
        "tableau" : "Bogacki-Shampine-4-2-3",
        "rtol" : 1e-5,
        "atol" : 1e-6
    }
}
```
Finally, we need to initialize the diag Variables
```cpp
myproject::Variables var = {rhs, grid, omega};
// trigger first computation of potential
{
    dg::x::DVec temp = omega;
    rhs( 0., omega, temp);
}

```

(sec:main_file_output)=
### Initialize the file output
```cpp
// Create netcdf file
dg::file::NcFile file;
try{
    file.open( argv[2], dg::file::nc_clobber);
}catch( std::exception& e)
{
    DG_RANK0 std::cerr << "ERROR creating file "<<argv[1]<<std::endl;
    DG_RANK0 std::cerr << e.what() << std::endl;
    dg::abort_program();
}
```
```{admonition} Meta-data
In order to fulfill the CF conventions there is a set of standard fields that should go into any netcdf file that you produce:
```
```cpp
std::map<std::string, dg::file::nc_att_t> atts;
atts["title"] = "Output file of myproject/myprogram.cpp";
atts["Conventions"] = "CF-1.8";
atts["history"] = dg::file::timestamp( argc, argv);
atts["comment"] = "Find more info in myproject/documentation.tex";
atts["source"] = "FELTOR";
atts["references"] = "https://github.com/myname/myproject";
// Here we put the inputfile as a string without comments so that it can be read later by another parser
atts["inputfile"] = js.toStyledString();
file.put_atts(atts);
file.put_atts( dg::file::version_flags);
```
```{admonition} Compile time Macros
The relevant macros `GIT_HASH`, `GIT_BRANCH` and `COMPILE_TIME` used in `dg::file::version_flags` need to be defined at compile time. We provide the Makefile `feltor/config/version.mk`. Simply include it in your own Makefile and add `$(VERSION_FLAGS)` to your compilation command!
```

### File output
For two-dimensional and three-dimensional programs it is often necessary to save on storage space when writing files. For this reason not every timestep is written to file but only every 100th say. Of course, the compression in time can be combined with a compression in space as well. In order to do so you would typically interpolate the simulation results onto a lower resolution grid.
```{note}
In this example we use the serial netcdf approach where all mpi threads send their data to the master thread which will write it to file. Pay special attention to which functions are called for all threads and which ones are only called by the master thread.
```
Let us set this up via some input parameters:
```cpp

unsigned n_out     = js[ "output"]["n"].asUInt( 3);
unsigned Nx_out    = js[ "output"]["Nx"].asUInt( 48);
unsigned Ny_out    = js[ "output"]["Ny"].asUInt( 48);

dg::x::CartesianGrid2d grid_out( x0, x1, y0, y1,
            n_out, Nx_out, Ny_out, bcx, bcy
            #ifdef WITH_MPI
            , comm
            #endif //WITH_MPI
            );
dg::x::IHMatrix projection = dg::create::interpolation( grid_out, grid);
```
#### Defining dimensions and variables
```cpp
// Define unlimited time dimension and variables
file.def_dimvar_as<double>( "time", NC_UNLIMITED, {{"axis", "T"}});
// the dimensions are the ones of grid_out!
file.defput_dim( "x", {{"axis", "X"},
    {"long_name", "x-coordinate in Cartesian system"}},
    grid_out.abscissas(0));
file.defput_dim( "y", {{"axis", "Y"},
    {"long_name", "y-coordinate in Cartesian system"}},
    grid_out.abscissas(1));

for( auto& record : myproject::diagnostics2d_list)
{
    file.def_var_as<double>( record.name, {"time", "y","x"}, {{"long_name",
        record.long_name}});
    // and the 1d fields (our idea is to just volume integrate the 2d fields
    file.def_var_as<double>( record.name + "_1d", {"time"}, {{"long_name",
        record.long_name + " (Volume integrated)"}});
}
```
#### Output the static list

For the static (time independent) data we do not need to store the variable ids, we can directly write the output.
```cpp
dg::x::HVec resultH = dg::evaluate( dg::zero, grid);
dg::x::HVec transferH = dg::evaluate( dg::zero, grid_out);
dg::x::DVec resultD = transferH; // transfer to device
for( auto& record : myproject::diagnostics2d_static_list)
{
    record.function( resultD, var);
    dg::assign( resultD, resultH);
    dg::blas2::gemv( projection, resultH, transferH);
    file.defput_var( record.name, {"y", "x"}, {{"long_name",
        record.long_name}}, grid_out, transferH);
}
```
#### First file output
```cpp
dg::x::DVec volume = dg::create::volume( grid);
size_t start = 0;
for( auto& record : myproject::diagnostics2d_list)
{
    record.function( resultD, var);
    double result = dg::blas1::dot( volume, resultD);
    dg::assign( resultD, resultH);
    dg::blas2::gemv( projection, resultH, transferH);
    file.put_var( record.name, {start, grid_out}, transferH);
    file.put_var( record.name+"_1d", {start, grid_out}, result);
}
file.put_var( "time", {start}, time);
file.close();
```

### The timeloop
We are finally ready to construct the main timeloop.
Normally, a file and/or dignostic output is wanted only at a fixed interval, determined by two of the three parameters
\begin{align}
    \Delta T_{output} = \frac{T_{end}}{N_{output}}
\end{align}
where $T_{end}$ is the simulation end time, $N_{output}$ is the number
of outputs in the output file and $\Delta T_{output}$ is the time
interval between outputs. Care must be taken that $N_{output}$ is integer.
```{note}
Usually, the adaptive time-output
given by simply letting an adaptive time-stepper run for a fixed number
of steps is not useful as it makes comparing simulations (with different spatial resolution say) difficult and
takes away control from the user.
```
Optionally, there is also a second diagnostic frequency, where for example higher time-resolved (low spatial) diagnostics need to be computed.
This requires a second output frequency
\begin{align}
    \delta T_{diag} = \frac{\Delta T_{output}}{N_{diag}}
\end{align}
```{note}
Especially for this higher frequency it is usually not necessary to output at _exactly_  the intermediary timesteps, just _close to those timesteps_. This is relevant because it is better to let the timestepper chose its own timestep over forcing it to land on a given point in time. This is shown in the following loop.
```
```cpp
double Tend = js["output"].get("tend", 1.0).asDouble();
unsigned maxout = js["output"].get("maxout", 10).asUInt();
double deltaT = Tend/(double)maxout;
bool abort = false;
for( unsigned u=1; u<=maxout; u++)
{
    try{
        // the documentation of dg::aTimeloop holds more details about how this construct works ...
        timeloop.integrate( time, omega, u*deltaT, omega,
                          u < maxout ? dg::to::at_least : dg::to::exact);
    }catch ( std::exception& fail)
    {
        DG_RANK0 std::cerr << "ERROR in Timestepper\n";
        DG_RANK0 std::cerr << fail.what() << std::endl;
        DG_RANK0 std::cerr << "Writing last output and exit ..."<<std::endl;
        abort = true;
    }
    start = u;
    file.open( argv[2], dg::file::nc_write);
    // First write the time variable
    file.put_var( "time", {start}, time);
    for( auto& record : myproject::diagnostics2d_list)
    {
        record.function( resultD, var);
        double result = dg::blas1::dot( volume, resultD);
        dg::assign( resultD, resultH);
        dg::blas2::gemv( projection, resultH, transferH);
        file.put_var( record.name, {start, grid_out}, transferH);
        file.put_var( record.name+"_1d", {start, grid_out}, result);
    }
    file.close();
    if( abort) break;
}
```
```{admonition} Error handling mechanism
Here, a graceful exit strategy is
suggested, where an output is still written in case of failure in order to
analyse the reasons for failure.
```
```json
{
    "output":
    {
        "n" : 3,
        "Nx" : 20,
        "Ny" : 20,
        "tend" : 1.0,
        "maxout" : 10
    }
}
```

### Closing
```cpp
#ifdef WITH_MPI
    MPI_Finalize();
#endif // WITH_MPI
    return 0;
} //end of main
```

(sec:Makefile)=
## The Makefile
The easiest way to configure the Makefile variables is to simply include Feltor's configuration in your own Makefile
```make
# Makefile
device=omp
FELTOR_PATH=../feltor

#configure machine
include $(FELTOR_PATH)/config/default.mk
include $(FELTOR_PATH)/config/version.mk
include $(FELTOR_PATH)/config/*.mk
include $(FELTOR_PATH)/config/devices/devices.mk

INCLUDE+=-I$(FELTOR_PATH)/inc/

all: myprogram_hpc myprogram_mpi

# only necessary if you use the draw library
#myprogram: myprogram.cpp equations.h init.h diag.h
#    $(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(GLFLAGS) $(JSONLIB) -g

myprogram_hpc: myprogram.cpp equations.h init.h diag.h
    $(CC) $(OPT) $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) $(VERSION_FLAGS) -g

myprogram_mpi: myprogram.cpp equations.h init.h diag.h
    $(MPICC) $(OPT) $(MPICFLAGS) $< -o $@ $(INCLUDE) $(LIBS) $(JSONLIB) -DWITH_MPI $(VERSION_FLAGS) -g

.PHONY: clean

clean:
    rm -rf myprogram myprogram_hpc myprogram_mpi
```

```{admonition} git
If you change the `FELTOR_PATH` locally and you want git to ignore this line in the Makefile when you commit, you can follow this [stackoverflow question](https://stackoverflow.com/questions/6557467/can-git-ignore-a-specific-line)
```

Now you can compile your program like any feltor program. For example
```bash
make myprogram_hpc device=gpu # Compile for shared memory GPU platform
# or
make myprogram_mpi device=cpu # Compile pure mpi version
```
Then call the program with
```bash
./myprogram_hpc input.json output.nc
# or
echo 2 2 | mpirun -n 4 ./myprogram_mpi input.json output.nc
```
