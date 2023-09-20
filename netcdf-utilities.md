# NetCDF utilities

As a binary output format we recommend [netCDF-4](https://www.unidata.ucar.edu/software/netcdf/). This format is widely used and bases on the HDF5 data format. Even though netCDF-4 offers an enhanced data model we usually do not create groups in files. In short, we use netcdf-4 for its HDF5 compatibility but only write the [The classic data model](https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_data_model.html#classic_model) consisting of _Variables_, _Dimensions_ and _Attributes_.
```{admonition} Conventions
We usually provide Meta-data in the file based on the [CF-conventions](https://cfconventions.org/) and [netCDF-conventions](https://www.unidata.ucar.edu/software/netcdf/documentation/NUG/best_practices.html). See the [Create your own project](sec:project) page for more details on how to do this in practice.
```

```{seealso}
Read more about the [NetCDF user-guide](https://docs.unidata.ucar.edu/nug/current/best_practices.html) and [C documentation](https://www.unidata.ucar.edu/netcdf/docs) and [python documentation](http://unidata.github.io/netcdf4-python/)
```
The goal of writing utility functions is to
- simplify error handling,
- simplify the creation of Dimensions for quantities that are defined on `dg` grids
- simplify the data management for writing `dg` data to file involved in a multi-threaded MPI program.

## Handling errors from netcdf calls
We provide a convenient class `dg.:file::NC_Error_Handle` which constructs from a netcdf return integer. If the integer is not zero it will throw.

``` cpp
#include "dg/file/nc_utilities.h"

dg::file::NC_Error_Handle err;
int ncid=-1;
try{
    err = nc_create( "outputfile.nc", NC_NETCDF4|NC_CLOBBER, &ncid);
}catch( std::exception& e)
{
    std::cerr << "ERROR creating file outputfile.nc"<<std::endl;
    std::cerr << e.what()<<std::endl;
}
```
## Defining dimensions
We simplify the creation of dimensions given a `dg` grid. The product space nature of our grid maps directly to the netCDF data model. For example
```cpp
int dim_ids[3], tvarID;
dg::CartesianGrid2d grid(...);
err = dg::file::define_dimensions( ncid, dim_ids, &tvarID, grid,
                {"time", "y", "x"});
```

creates three one-dimensional dimensions "time", "x" and "y" and corresponding dimension variables with the same names. The data for "x" and "y" is written directly to file and is generated from the Gaussian nodes of the grid in x and y directly. The time variable is an _unlimited_ variable expecting a time simulation with an unknown number of steps. Therefore, the function returns the id of the time variable for the user to write.
```{seealso}
See the `dg::file::define_dimensions` family of functions in the doxygen documentation. There is one for each grid.
```

## Defining and writing variables
The netCDF-4 standard mandates that a variable should have dimensions and data.
The netCDF C-interface already defines ready-to-use functions `nc_def_var`, `nc_put_var_double` (for writing variables in a single call) and `nc_put_vara_double` (for writing variables in chunks) for defining variables and writing data to file. The only issue when writing data to file is what to do in an MPI setting where the data is distributed among processes. NetCDF offers a parallel writing backend, however in our tests this turned out to be slow and the compilation and linking of a program becomes more complicated.

A better approach is to use **serial netCDF**. Simply send all data to the master process, which then funnels the data into the output file. The management of these data transfers and communication is hidden in the `dg::file::put_var_double` and `dg::file::put_vara_double` family of functions.
```cpp
int varID;
std::string name = "variable";
// only the master thread needs to define the variable
DG_RANK0 err = nc_def_var( ncid, name.data(), NC_DOUBLE, 3, dim_ids, &varID);
// generate data
dg::x::HVec transferH = ...;
int start = 0; // which timestep
// all threads need to call the writing function
dg::file::put_vara_double( ncid, varID, start, grid, transferH);
```
In a serial environment the `dg::file::put_var[a]_double` functions become simple wrappers around the corresponding netCDF function. They can therefore be used in a platform independent environment.
````{note}
If you actually do want to use the _parallel netcdf_ interface, you can use a hidden parameter in the function
```cpp
// now each process writes to the file in parallel
dg::file::put_vara_double( ncid, varID, start, grid, transferH, true);
```
````
```cpp
 err = nc_close(ncid);
 ```
