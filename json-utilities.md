# JSON utilities

We use the [jsoncpp](https://github.com/open-source-parsers/jsoncpp) library
to parse json files. We chose it because it is widely used, very well maintained and even available as a linux package `libjsoncpp`.

## The `dg::file::file2Json` function
The first simplification regards how comments and errors are treated on file opening. We provide a shortcut:

```cpp
#include <iostream>
#include "dg/algorithm.h"
#include "dg/file/json_utilities.h"
dg::file::WrappedJsonValue js;
try{
    dg::file::file2Json( "does_not_exist.json", js.asJson(),
                    dg::file::comments::are_discarded,
                    dg::file::error::is_throw);
}catch( std::exception& e)
{
    std::cout << e.what();
}
```

An error occured while parsing does_not_exist.json
*** File does not exist! ***

In this example we open a file where C-style comments are allowed but discarded, while an error on opening or reading the file leads to a throw (which if not captured leads to immediate abortion of the program).

## The `dg::file::WrappedJsonValue` class
For our purposes the only downside of jsoncpp is that missing values do not trigger a throw and to manually check existence somewhat clutters the code. At the same time missing values often result from silly mistakes, which in a high performance computing environment result in real avoidable cost (for example because you ran a simulation with a default value that you did not really intend to).

For this reason we provide the `dg::file::WrappedJsonValue` class.
It basically wraps the
access to a `Json::Value` with guards that raise exceptions or display warnings in case an error occurs, for example when a key is misspelled,
missing or has the wrong type.
The goal is the composition of a good error message that helps a user
quickly debug the input file.

The Wrapper is necessary because Jsoncpp by default silently
generates a new key in case it is not present which in our scenario is an
invitation for stupid mistakes.
You can use the `WrappedJsonValue` like a `Json::Value` with read-only access:

```cpp
dg::file::WrappedJsonValue ws;
try{
    dg::file::file2Json( "test.json", ws.asJson(),
                    dg::file::comments::are_discarded,
                    dg::file::error::is_throw);
}catch( std::exception& e){ std::cerr << e.what();}
try{
    // Access a nested unsigned value
    unsigned n = ws["grid"]["n"].asUInt();
    // Access a nested unsigned value
    unsigned Nx = ws["grid"]["Nx"].asUInt();
    // Access a list item
    double x0 = ws["grid"]["x"][0].asDouble();
    // Acces using default initializer
    double x1 = ws["grid"]["x"].get( 1, 42.0).asDouble();
    // Access a string
    std::string bc = ws["bc"][0].asString();
    // Example of a throw
    std::string hello = ws[ "does not exist"].asString();
} catch ( std::exception& e){
    std::cout << "Error in file test.json\n";
    std::cout << e.what()<<std::endl;
}
```

Error in file test.json
*** Key error: "does not exist":  not found.

A feature of the class is that it keeps track of how a value is called.
For example
```cpp
void some_function( dg::file::WrappedJsonValue ws)
{
    int value = ws[ "some_non_existent_key"].asUInt();
    std::cout << value<<"\n";
}
try{
    some_function( ws["grid"]);
} catch ( std::exception& e){ std::cout << e.what()<<std::endl; }
// *** Key error: "grid": "some_non_existent_key":  not found.
```

The what string knows that "some_non_existent_key" is expected to be
contained in the "grid" key, which simplifies debugging to a great extent.
