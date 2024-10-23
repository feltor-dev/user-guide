# A high level C++ NetCDF interface

It turns out that even with our netcdf-utilities the simple task of "just write the data to file" is still
not easy enough to become "simple". This is mainly because we still need to manually manage dimension and variable ids.

For this reason we provide the `dg::file::Reader` and the `dg::file::Writer` class with the express
purpose to make netCDF file reading and writing as convenient as possible.

## How to start

The high-level interface can be accessed by including
``` cpp
#include "dg/file/file.h"
```
This header automatically includes `"dg/file/json_utilities.h"`
and `"dg/file/netcdf_utilities.h"`.
```{note}
The high-level NetCDF interface depends on netcdf **and** json (mainly for convenient attribute reading and writing).
Read [](json-utilities) to decide which json library to use.
```

## A first example

Let us assume that we have a vector that lives on a 2d grid that we want to dump to a NetCDF file.
We can do this like so
``` cpp
#include <iostream>
#include "dg/file/file.h"


int main()
{
    // create a netcdf file
    int ncid;
    dg::file::NC_Error_Handle err;
    err = nc_create( "test.nc", NC_NETCDF4|NC_CLOBERR, &ncid);

    // create some data
    dg::Grid2d g2d( 0,1 , 0,1, 3, 20, 20);
    auto data = dg::evaluate( [](double x, double y){ return sin(x)*sin(y);}, g2d);

    // Create a writer and use it to write data
    dg::file::Writer<dg::Grid2d> writer( ncid, g2d, {"y", "x"});
    writer.def_and_put( "Vector", {}, data);

    // Close the file
    err = nc_close( ncid);

    return 0;
}
```

The design philosophy of the `dg::file::Writer` class is that it manages the dimensions and variable writing for variables that have equal type and shape. The type and dimensionality is given by the `dg::Grid2d` template parameter while the shape is inferred from the `g2d` parameter.
The two dimensions created are called "y" and "x" and the data "Vector" is created without attributes.

## A time dependent field

In this example we use show how to write a time dependent field to file and also make it agnostic of whether MPI is used or not
``` cpp
#include <iostream>
#include "dg/file/file.h"


int main()
{
    // create a netcdf file
    int ncid;
    dg::file::NC_Error_Handle err;
    DG_RANK0 err = nc_create( "test.nc", NC_NETCDF4|NC_CLOBERR, &ncid);

    // create some data
    dg::x::Grid2d g2d( 0,1 , 0,1, 3, 20, 20
#ifdef WITH_MPI
    , comm
#endif
    );
    dg::file::Writer<dg::x::Grid0d> write0d( ncid, {}, {"time"});
    write0d.def("time");
    write0d.def("Energy");
    // Note that a time dependent writer CANNOT write time-independent data
    dg::file::Writer<dg::x::Grid2d> write2d( ncid, g2d, {"time", "y", "x"});
    write2d.def( "Vector");

    double Tmax=2.*M_PI;
    double NT = 10;
    double h = Tmax/NT;
    for(unsigned i=0; i<=NT; i++)
    {
        // compute the data
        double time = i*h;
        auto data = dg::evaluate( [](double x, double y){ return sin(x)*sin(y);}, g2d);
        dg::blas1::scal( data, cos( time));
        double energy = dg::blas1::dot( data, data);

        // Write the data
        write0d.stack( "time", time);
        write0d.stack( "Energy", energy);
        write2d.stack( "Vector", data);
    }

    // Close the file
    DG_RANK0 err = nc_close( ncid);

    return 0;
}
```

## Writing and reading global attributes
Here we have the two functions `dg::file::json2nc_attrs` and `dg::file::nc_attrs2json`, which convert
json objects to netCDF attributes and back.

``` cpp
#include <iostream>
#include <ctime>
#include <iomanip>

#include "dg/file/file.h"


int main(int argc, char* argv[])
{
    dg::file::JsonType att;
    att["text"] = "Hello World!";
    att["number"] = 3e-4;
    att["int"] = -1;
    att["uint"] = 10;
    att["bool"] = true;
    att["array"] = dg::file::vec2json({-1.1, 42.3});

    int ncid;
    dg::file::NC_Error_Handle err;
    err = nc_create( "atts.nc", NC_NETCDF4|NC_CLOBBER, &ncid);

    // Write global attributes
    dg::file::json2nc_attrs( att, ncid, NC_GLOBAL);

    err = nc_close(ncid);

    std::cout << "Now test reading of attributes\n";
    err = nc_open( "atts.nc", 0, &ncid);

    // Read global attributes
    dg::file::JsonType read = dg::file::nc_attrs2json( ncid, NC_GLOBAL);

    std::cout << dg::file::WrappedJsonValue(read).toStyledString()<<"\n";
    err = nc_close(ncid);

    return 0;
}
```

## Caveats

### What if my grid is not equidistant, curvilinear or even unstructured?
Currently, our interface does not support the creation of non-equidistant dimensions (except for the non-equidistant dG nodes).

In the case that your grid is curvilinear or unstructured you need to pull the coordinates from physical space back to an equidistant computational space.
This is done by creating coordinate variables that can be used to map from the equidistant output grid to the actual grid of interest.

For example let's revisit the first example on a logarithmic grid (between $10^{-3}$ and $10^2$):
``` cpp
#include <iostream>
#include "dg/file/file.h"


int main()
{
    // create a netcdf file
    int ncid;
    dg::file::NC_Error_Handle err;
    err = nc_create( "test.nc", NC_NETCDF4|NC_CLOBERR, &ncid);

    // create some data
    unsigned N = 20;
    // Computational grid
    dg::Grid2d g2d( -3,2 , -3,2, 1, N, N); //equidistant cell-centered grid
    // If you want to include the end points we need to massage the Grid constructor a bit:
    // double x0 = -3, x1 = 2;
    // double h = (x1-x0)/(N-1)
    // dg::Grid2d g2d( -3-h,2+h , -3-h,2+h, 1, N, N); //equidistant edge centered grid
    // The curvilinear coordinate map X(x,y), Y(x,y)
    auto xpoints = dg::evaluate( [] (double x, double y){ return pow( 10, x);}, g2d);
    auto ypoints = dg::evaluate( [] (double x, double y){ return pow( 10, y);}, g2d);

    auto data = dg::blas1::evaluate( [](double x, double y){ return sin(x)*sin(y);},
                                     xpoints, ypoints);

    // Write the data and the actual grid points as variables
    dg::file::Writer<dg::Grid2d> writer( ncid, g2d, {"y", "x"});
    writer.def_and_put( "Vector", {}, data);
    writer.def_and_put( "YY", {}, ypoints);
    writer.def_and_put( "XX", {}, xpoints);

    // Close the file
    err = nc_close( ncid);

    return 0;
}
```

### What if I have data living on higher dimensional (>3) grids?
Currently, Feltor does not support more than three spatial dimensions, even though it is a somewhat small upgrade if it is ever needed. Write a support ticket if you really need it.

