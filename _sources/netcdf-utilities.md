# NetCDF utilities

As a binary output format we recommend [netCDF-4](https://www.unidata.ucar.edu/software/netcdf/). This format is widely used and bases on the HDF5 data format. A good start is to read about the [The NetCDF data model](https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_data_model.html) consisting of _Variables_, _Dimensions_ and _Attributes_.
```{admonition} Conventions
We usually provide Meta-data in the file based on the [CF-conventions](https://cfconventions.org/) and [netCDF-conventions](https://www.unidata.ucar.edu/software/netcdf/documentation/NUG/best_practices.html). See the [Create your own project](sec:project) page for more details on how to do this in practice.
```

```{seealso}
Read more about the [NetCDF user-guide](https://docs.unidata.ucar.edu/nug/current/best_practices.html) and [C documentation](https://www.unidata.ucar.edu/netcdf/docs) and [python documentation](http://unidata.github.io/netcdf4-python/)
```
Unfortunately, netCDF is a C-library with a from a C++ perspective frankly terrible interface.
In particular the task of "just write it to file" can be daunting with pure netCDF functions involving many low-level steps like creating and writing dimensions, managing variable ids, creating variables, and then only finally writing a variable. A bit further reaching is the question how to write to file in a MPI environment.

For this reason we provide the `dg::file::NcFile` class, which is a modern
C++-17 style implementation of the [NetCDF Data Model](https://docs.unidata.ucar.edu/netcdf-c/4.9.2/netcdf_data_model.html).
The goal of the class is to simplify reading and writing of data from/to NetCDF files to the degree that you can "just do it".

## Opening and closing a file
Opening and closing a file is as simple as
``` cpp
#include <iostream>
#include "dg/file/nc_utilities.h"

int main()
{
    // Create a file "outputfile.nc", overwrite ("clobber") if exists
    dg::file::NcFile file( "outputfile.nc", dg::file::nc_clobber);
    // Close the file
    file.close();
    return 0;
}
```
## Defining Attributes for a file
Attributes can be written with the `put_att` family of member functions
``` cpp
#include <iostream>
#include "dg/file/nc_utilities.h"

int main()
{
    // Create a file "outputfile.nc", overwrite ("clobber") if exists
    dg::file::NcFile file( "outputfile.nc", dg::file::nc_clobber);
    // A string attribute named "title" for the file
    file.put_att( {"title", "My title"});
    // A integer attribute named "truth"
    file.put_att( {"truth", 42});
    // Close the file
    file.close();
    return 0;
}
```

## Defining dimensions
We simplify the creation of dimensions given a `dg` grid. The product space nature of our grid maps directly to the netCDF data model. For example:
```cpp
    // Generate a grid
    const double x0 = 0., x1 = 2.*M_PI;
    dg::x::CartesianGrid2d grid( x0,x1,x0,x1,3,10,10);
    // and put dimensions to file
    file.defput_dim( "x", {{"axis", "X"},
        {"long_name", "x-coordinate in Cartesian system"}},
        grid.abscissas(0));
    file.defput_dim( "y", {{"axis", "Y"},
        {"long_name", "y-coordinate in Cartesian system"}},
        grid.abscissas(1));
```
creates two one-dimensional dimensions "x" and "y" and corresponding dimension variables with the same names. The data for "x" and "y" is written directly to file and is generated from the Gaussian nodes of the grid in x and y directly.
The two dimension variables are created together with two attritubes "axis" and "long_name" each.

## Defining and writing static variables
The netCDF-4 standard mandates that a variable should have dimensions and data and can have attributes.
```cpp
    // Generate some data and write to file
    dg::x::HVec data = dg::evaluate( function, grid);
    // Defne and write a variable in one go
    file.defput_var( "variable", {"y", "x"},
                {{"long_name", "A long explanation"}, {"unit", "m/s"}},
                grid, data);
```
Here, we defined a variable named "variable" with dimensions "x" and "y" and two attributes "long_name" and "unit". Note here that the "y" dimension comes first because our dg grids have the "x" dimension contiguous in memory.
Also note that grid appears in the defput member function. This is necessary to tell the library the extent of the
given data to write (Theoretically, you could just write a portion of the data, or in an MPI environment the grid tells the library which part of the data it should write).

## Defining and writing dynamic variables
The first thing to do to create a time series is to create an unlimited time dimension.
Afterwards we can define another variable, which depends on time and then write a time-series
```cpp
    // Generate an unlimited dimension and define another variable
    file.def_dimvar_as<double>( "time", NC_UNLIMITED, {{"axis", "T"}});
    file.def_var_as<double>( "dependent", {"time", "y", "x"},
        {{"long_name", "Really interesting"}});
    // Write timeseries
    for( unsigned u=0; u<=2; u++)
    {
        double time = u*0.01;
        file.put_var("time", {u}, time);
        // We can write directly from GPU
        dg::x::DVec data = dg::evaluate( function, grid);
        dg::blas1::scal( data, cos(time));
        file.put_var( "dependent", {u, grid}, data);
    }
```

## Reading variables
In order to read the data that we have just written we open the file in read mode
```cpp
    // Open file for reading
    file.open( "test.nc", dg::file::nc_nowrite);
    std::string title = file.get_att_as<std::string>( "title");
    // In MPI all ranks automatically get the right chunk of data
    file.get_var( "variable", grid, data);
    unsigned NT = file.get_dim_size( "time");
    // CHECK( NT == 3);
    double time;
    file.get_var( "time", {0}, time);
    file.close();
```

## What if my grid is not equidistant, curvilinear or even unstructured?
The `dg::file::NcFile` is just a wrapper around the underlying netcdf C-library so in principle
you should be able to do whatever you want.
For example let's view the example of a logarithmic grid (between $10^{-3}$ and $10^2$):
``` cpp
#include <iostream>
#include "dg/file/file.h"

int main()
{
    // create a netcdf file
    dg::file::NcFile file( "test.nc", dg::file::nc_clobber);

    // Computational grid
    dg::Grid2d g2d( -3,2 , -3,2, 1, 10, 20); //equidistant cell-centered grid
    // If you want to include the end points we need to massage the Grid constructor a bit:
    // double x0 = -3, x1 = 2;
    // double h = (x1-x0)/(N-1)
    // dg::Grid2d g2d( -3-h,2+h , -3-h,2+h, 1, 10, 20); //equidistant edge centered grid

    // Product space
    auto xpoints = dg::evaluate( [] (double x){ return pow( 10, x);}, g2d.gx());
    auto ypoints = dg::evaluate( [] (double y){ return pow( 10, y);}, g2d.gy());

    auto data = dg::kronecker( [](double x, double y){ return sin(x)*sin(y);},
                                     xpoints, ypoints);

    // Write the data and the actual grid points as variables
    file.defput_dim( "x", {}, xpoints);
    file.defput_dim( "y", {}, ypoints);
    file.defput_var( "Vector", {"y", "x"}, {}, g2d, data);

    // Close the file
    file.close();

    return 0;
}
```

## What to do for MPI distributed vectors

When writing data to file in an MPI setting data is distributed among processes. NetCDF offers a parallel writing backend, however in our tests this turned out to be slow and the compilation and linking of a program becomes more complicated.
A better approach is to use **serial netCDF**. Simply send all data to the master process, which then funnels the data into the output file.
Our `dg::file::NcFile` class automatically recognises if it is being used in an MPI environment and automatically dispatches all calls to the corresponding MPI call. This means that in MPI the class can be used exactly like in a shared memory environment.
