(sec:welcome)=
# Welcome
This guide explains how to use Feltor.
It's structure follows the basic structure of Feltor itself as described
further below.
If you want, you can also follow a
[lecture by Matthias Wiesenberger](https://events.prace-ri.eu/event/989/sessions/3081/attachments/1199/2017/Wiesenberger_PRACE_high.mp4) held at the PRACE winter school on
GPU programming in Innsbruck (unfortunately takes some time to load).
If you haven't done so yet, please read the
[Quick Start Guide](https://github.com/feltor-dev/feltor) first, which
explains how to install the library and compile programs.

```{note}
This guide is generated using [jupyter-book](https://jupyterbook.org/intro.html) from [jupyter-notebooks](https://jupyter.org/) and Markdown files.
You can simply copy-paste the code examples into a textfile and compile the code yourself and we encourage you to do so and play around with the provided examples.
```

```{seealso}
You can look up any class or function beginning with `dg::` in the [doxygen documentation](https://mwiesenberger.github.io/feltor/dg/html/modules.html)
```

## What is FELTOR?

FELTOR (Full-F ELectromagnetic code in TORoidal geometry) is a modular
scientific software package. The following Figure shows its structure

![The structure of the FELTOR project](https://feltor-dev.github.io/images/FeltorStructure.png)


### Application codes

A collection of actual simulation projects and diagnostic
programs.  Execute one- two- and three-dimensional simulations or simply run a numerical algorithm for given parameters. Typically these programs read in input file(s), simulate, and either write results to disk or directly visualize them on screen. Some examples led to journal [publications](https://feltor-dev.github.io/publications) in the past.

### The dg library

The dg library is a **header-only** template C++ library separated into modules.
The core dg library header `dg/algorithm.h` includes the core modules in the above structure

**Advanced numerical algorithms**

Numerical schemes that are based on the existence of a geometry and/or a topology. These include notably the discretization of elliptic equations in arbitrary coordinates, multigrid algorithms and the flux coordinate independent approach in arbitrary coordinates (available through the _geometries_ extension `dg/geometries/geometries.h`).

**Topology and Geometry**

Here, we introduce data structures and functions that represent the concepts of Topology and Geometry and operations defined on them (for example the discontinuous Galerkin discretization of derivatives). The _geometries_ extension implements a large variety of grids and grid generation algorithms that can be used here.

**Basic numerical algorithms**

Algorithms like conjugate gradient (CG) or Runge-Kutta schemes that can be implemented with vector operations alone.

**Basic parallel operations**

In this "hardware abstraction" module we define the interface for various vector and matrix operations like additions, multiplications, scalar products and so on. These functions are then implemented  and optimized on a variety of hardware architectures and serve as building blocks for all higher
level algorithms.

**Geometries extension** `dg/geometries/geometries.h`

adds several grid generators, magnetic field structure and the flux-coordinate independent approach

**Matrix functions** `dg/matrix/matrix.h`

 adds matrix-function computations and exponential integrators (depends on boost and lapack libraries)

**File I/O operations** `dg/file/file.h`

simplifies common I/O operations in our programs (depends on NetCDF and jsoncpp libraries)

**Exblas** `dg/exblas/exblas.h`

is special because it is already included in the dg library  but can also be used as a standalone library. It provides binary reproducible and accurate scalar products on various architectures.
