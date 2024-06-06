(sec:vectors)=
# Vectors
Assume we want to compute $\vec y =  a\vec x+b\vec y$ where $\vec x$ and $\vec y$ are vectors and a and b are constants.
Furthermore, assume we wish to compute the dot product
of $\vec x$ and $\vec y$. And finally assume that we want to do this for both small
and large vectors, that is we want to eventually write parallel code using for example GPUs.
The following tutorial shows how to do exactly that
and more using the dg library.

##  A first example
First, we have to include the main feltor library header:

```cpp
// include the dg-library
#include "dg/algorithm.h"

double a = 0.5, b = 0.25;
std::array<double, 2> x ={2,2}, y = {4,4};
// compute a*x+b*y and store it in y
dg::blas1::axpby( a, x, b, y);
// compute Sum_i x_i y_i
double sum = dg::blas1::dot( x,y);
// output should be 8
std::cout << sum << std::endl;
// 8
```
In this code we encounter our first two dg functions, namely `dg::blas1::axpby`
and `dg::blas1::dot`. They perform very basic operations, namely adding vectors
 and computing scalar products respectively. (You can look up their formal documentation [here](https://mwiesenberger.github.io/feltor/dg/html/group__blas1.html)).
The remarkable thing about these two functions is that they are templates.
This means you can call them for many different vector classes. We can change the type of `x` and `y`;

```cpp
// works as well
std::vector<double> x(2,2), y(2,4);
dg::blas1::axpby( a, x, b, y);
double sum = dg::blas1::dot( x,y);
std::cout << sum << std::endl;
// 8
//or
thrust::host_vector<double> x(2,2), y(2,4);
dg::blas1::axpby( a, x, b, y);
double sum = dg::blas1::dot( x,y);
std::cout << sum << std::endl;
// 8
```
All three versions have the same result.

```{note}
All of these examples
execute on a single CPU thread, that is the compiler chooses the same
**serial** implementation of the vector addition and the scalar product.
```

So let us increase the vector size to say, a Million. Wouldn't it
be better to perform these operations in parallel? And a measurement
of the execution time would also be nice:

```cpp
//use the thrust library to allocate memory on the device
thrust::device_vector<double> x(1e6,2), y(1e6,4);
//create a Timer
dg::Timer t;
//start the clock
t.tic();
dg::blas1::axpby( a, x, b, y);
//stop the clock
t.toc();
std::cout << "Axpby took "<<t.diff()<<"s\n";
t.tic();
double sum = dg::blas1::dot( x,y);
t.toc();
std::cout << "Dot   took "<<t.diff()<<"s\n";
//output should be ... large
std::cout << sum << std::endl;

//Axpby took 0.913942s
//Dot   took 0.883924s
//4e+06

```

The first thing to notice is that we now use the thrust::device_vector<double> class. This is a vector class of the thrust library, which allocates memory on a GPU. The compiler recognizes that dg::blas1::axpby and dg::blas1::dot are now called with a GPU vector class and redirects the call to the corresponding CUDA implementation.

 When using a GPU a natural question is how the synchronization is done between the GPU kernel queue and the host execution. The good news is that you do not need to manually synchronize the GPUs. The library will do it for you whenever it is needed.

If you do not have a GPU you can also define the THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_OMP Macro, then the call redirects to a OpenMP parallelized version or THRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP (like it is done in this notebook because xeus-cling does not support OpenMP or cuda very well) then the calls redirect to the serial version again. This is a specialty of the thrust::device_vector class of the thrust library that is included in dg/algorithm.h.

 You could also use your own
 vector class in the `dg::blas1` functions. This is an advanced feature
 and requires you to provide a specialization of `dg::TensorTraits`
 for your class, where you specify the parallelization strategy that
 the libary should choose and how the data is laid out in memory.
 Please consult the [documentation](https://mwiesenberger.github.io/feltor/dg/html/index.html#dispatch) for further details on
how we dispatch the blas functions
and our template traits system.

The dg::Timer measures the time it took to execute the functions.

So, what if the vector size is even larger 1e8 say? Then an MPI implementation would be handy, wouldn't it:

```{note}
The following cell will not compile if you do not have MPI installed.. You will need to copy-paste the below code to a file and compile it using the instructions in the quick-start guide.  
```

```cpp
//activate MPI in FELTOR
#include "mpi.h" //needs to be included **before** dg/algorithm.h
#include "dg/algorithm.h"


int main(int argc, char* argv[])
{
    //init MPI
    MPI_Init( &argc, &argv);
    //let's take all processes
    MPI_Comm comm = MPI_COMM_WORLD;
    //get the number of MPI processes in the communicator
    int np,rank;
    MPI_Comm_size(comm, &np);
    //get the rank of the calling process
    MPI_Comm_rank(comm, &rank);
    //allocate and initialize local memory
    thrust::device_vector<double> x_local( 1e8/np, 2), y_local(1e8/np, 4);
    //combine the local vectors to a global MPI vector
    dg::MPI_Vector<thrust::device_vector<double>> x(x_local, comm);
    dg::MPI_Vector<thrust::device_vector<double>> y(y_local, comm);

    //now repeat the operations from before...
    double a = 0.5, b = 0.25;
    dg::Timer t;
    t.tic();
    dg::blas1::axpby( a, x, b, y);
    t.toc();
    if(rank==0)std::cout << "Axpby took "<<t.diff()<<"s\n";
    t.tic();
    double sum = dg::blas1::dot( x,y);
    t.toc();
    if(rank==0)std::cout << "Dot   took "<<t.diff()<<"s\n";
    if(rank==0)std::cout << sum << std::endl;
    //be a good MPI citizen and clean up
    MPI_Finalize();
    return 0;
}
```
```{note}
We have just written a hybrid MPI + X code, where X can be either serial, OpenMP or CUDA depending on how the code is compiled!
```
One remaining thing is that we quickly get tired
 of writing `thrust::device_vector<double>` and
especially `dg::MPI_Vector<thrust::device_vector<double>>`.
 So we invented convenient typedefs:
`
dg::DVec x_local( 1e8/np, 2), y_local(1e8/np, 4)
dg::x::DVec x(x_local, comm), y(y_local, comm);
`

 which is completely equivalent to the corresponding lines 18 and 20/ 21 above.

## Platform independent code

The remarkable thing in the above examples is that two lines of code never changed, even if we changed the class of vectors that we use it on

```cpp
// This is platform independent code
dg::blas1::axpby( a, x, b, y);
double sum = dg::blas1::dot( x,y);
```
This is a first example of **platform independent code**.
It gives us an abstract way to perform vector operations on many different types and with different parallelization. If we implement an algorithm only in terms of `dg::blas1` or similar functions, then
this algorithm is **automatically parallelized** for a variety of hardware and works with **many different vector types**.

```{note}
This is in fact exactly what the dg library does. It provides a set of algorithms all implemented as templates of (at least) the vector class using only the dg::blas1 (and later dg::blas2) functions. We will encounter a first example of such an algorithm in the next chapter of the tutorial: timesteppers.
```

## Arbitrarily parallel - The blas1 subroutine
One remaining question in this chapter is: what if we do not want to add vectors
but multiply them instead? Or take the exponential of each element?
There is a selection of predefined `dg::blas1` operations
you can choose from for example `dg::blas1::pointwiseDot` or `dg::blas1::scal`. Check out the
[documentation](https://mwiesenberger.github.io/feltor/dg/html/group__blas1.html).

If this still does not satisfy your needs, the answer probably is
the `dg::blas1::subroutine`.
```{admonition} Really blas1?
`blas` stands for basic linear algebra subroutins, the 1 stands for pure vector operations. In this sense `blas1::pointwiseDot` or `blas1::subroutine` are actually misnomers, since neither function is linear, nor are they part of the original blas library.
```


 As an example, let us assume that we have two vectors `v` and `w` and we want to compute the elementwise expression $v_i = v_i v_i + a\sin(w_i)$, where $a$ is a given scalar. The first thing is to write a lambda with our expression and then apply it to all elements in `v` and `w`
 ```{admonition} Lambdas are amazing
Lambda functions are very useful and will appear continuously throughout this guide, so it is a good idea to refresh your knowledge about them.
A good watch on youtube is for example [Lambdas from Scratch - Arthur O'Dwye](https://www.youtube.com/watch?v=3jCOwajNch0)
```

```cpp
double a = 8.;
auto lambda = [&]DG_DEVICE( double& vi, double wi ){
                              vi = vi*vi + a*sin(wi);
                          };
dg::HVec v( 1e6,2), w(1e6,4);
dg::blas1::subroutine( lambda, v,w);
std::cout << v[0] << std::endl;
// -2.05442
```

In the above example we have two host vectors `dg::HVec v` and `w`. In a shared
 memory code these will be declared as `dg::DVec`
In an MPI implementation we would simply write `dg::x::DVec` instead of `dg::DVec`.

`dg::blas1::subroutine` simply calls the given **custom Functor** (in the above case a lambda) for each element in `v` and `w`. There can be an **arbitrary** number of arguments to the subroutine function (in fact there have to be as many as arguments to the functor/lambda).

This in fact means that `dg::blas1::subroutine` is able to compute **any** trivially parallel
operation with **any** number of inputs and outputs.
```{note}
If you do not think this is pure magic, stop here and read again :)
(On a technical note, this neat behaviour is possible through `C++-11`-style
template parameter packs).
```

## A collection of vectors - recursive execution

In many practical problems of fluid dynamics we not only want to integrate one variable, `x` say in time, but several. For example we may want to integrate both density $n$ and velocity $v$. Then, the complete
set of variables should be tied together in a single entity.
In Feltor you do this simply by allocating "Vectors of vectors". The underlying engine automatically (and recursively) expands such constructions.
Let us see an example:

```cpp
double a = 0.5, b = 0.25;
dg::DVec two( 1e3, 2), four( 1e3, 4);
std::array<dg::DVec, 2> x = {two,two}, y = {four,four};
// compute a*x+b*y and store it in y
dg::blas1::axpby( a, x, b, y);
// compute Sum_i x_i y_i
double sum = dg::blas1::dot( x,y);
// output should be large ...
std::cout << sum << std::endl;
// 8000
```
In this example we have an array of two large vectors (that are not filled with too exciting values, but it is just to show that the code still works). It is also possible to have a vector or map of vectors

```cpp
double a = 0.5, b = 0.25;
dg::DVec two( 1e3, 2), four( 1e3, 4);
std::vector<dg::DVec> x {2,two}, y{2,four};
// compute a*x+b*y and store it in y
dg::blas1::axpby( a, x, b, y);
// compute Sum_i x_i y_i
double sum = dg::blas1::dot( x,y);
// output should be large ...
std::cout << sum << std::endl;
// 8000
```

```cpp
double a = 0.5, b = 0.25;
dg::DVec two( 1e3, 2), four( 1e3, 4);
std::map<std::string, dg::DVec> x {{ "density",two}, {"velocity", two}},
                                y {{ "density",four}, {"velocity", four}};
// compute a*x+b*y and store it in y
dg::blas1::axpby( a, x, b, y);
// compute Sum_i x_i y_i
double sum = dg::blas1::dot( x,y);
// output should be large ...
std::cout << sum << std::endl;
// 8000
```

Since the implementation works recursively, the "innver vector" can be a recusive vector again that is something like  `std::vector<std::array<dg::DVec,2>>` is also possible.
