# JSON utilities

We typically use and recommend [json](https://www.json.org/json-en.html) as a data format for input parameters to our codes. Json is easy for both humans and machines to read and write.

In order to parse json files in C++ we
recommend the [jsoncpp](https://github.com/open-source-parsers/jsoncpp)
or [nlohman/json](https://json.nlohmann.me)
library.
Both are widely used, very well maintained and even available as a linux package `libjsoncpp` and `nlohmann-json3-dev`.
The advantage of nlohmann-json is that it is header only, while jsoncpp needs to be precompiled and linked.

## How to start
Json parsers work quite well already out of the box. For some repetitive tasks specific to scientific simulations we have written utiltiy functions.
These can be accessed by influding
```cpp
#include "dg/file/json_utilities.h"
```
By default this header automatically includes jsoncpp's `"json/json.h"`
and thus incurs a dependency on jsoncpp. If one prefers to use nlohmann-json
one can define the macro
```cpp
#define DG_USE_JSONHPP
#include "dg/file/json_utilities.h"
```
Now, the file `<nlohmann/json.hpp>` is included instead of jsoncpp's header.
```{note}
Other json parsers are curretly not supported by our utility functions.
```

## File opening

Let us assume that in python we write an input file
```python
import json
params = {
    "grid" : {
        "n" : 3,
        "Nx" : 64,
        "x" : [0,1],
    },
    "bc" : ["DIR", "PER"],
}
with open("test.json") as f:
    json.dump( params, f, indent = 4)

```
Now, in C++ we open this file conveniently using

```cpp
#include <iostream>

#include "dg/file/json_utilities.h"

int main()
{
    dg::file::JsonType js; // either Json::Value or nlohmann::json
    try{
        js = dg::file::file2Json( "test.json");
    }catch( std::exception& e)
    {
        std::cout << e.what();
    }
    return 0;
}
```
In this example we open a file where C-style comments are allowed but discarded, while an error on opening or reading the file leads to a throw (which if we hadn't captured leads to immediate abortion of the program).

## The WrappedJsonValue class
For our purposes the only downside of jsoncpp or nlohmann-json is that missing values do not trigger a throw and to manually check existence somewhat clutters the code. At the same time missing values often result from silly mistakes, which in a high performance computing environment result in real avoidable cost (for example because you ran a simulation with a default value that you did not really intend to).

For this reason we provide the `dg::file::WrappedJsonValue` class.
It basically wraps the
access to a `Json::Value` with guards that raise exceptions or display warnings in case an error occurs, for example when a key is misspelled,
missing or has the wrong type.
The goal is the composition of a good error message that helps a user
quickly debug the input file.

The Wrapper is necessary because json parsers by default silently
generate a new key in case it is not present which in our scenario is an
invitation for stupid mistakes.
You can use the `WrappedJsonValue` like a `Json::Value` with read-only access:

```cpp
#include <iostream>

#include "dg/file/json_utilities.h"

int main()
{
    dg::file::WrappedJsonValue ws; // hide Json type ...
    try{
        ws = dg::file::file2Json( "test.json");
    }catch( std::exception& e)
    {
        std::cout << e.what();
    }
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
    // Error in file test.json
    // *** Key error: "does not exist":  not found.
    return 0;
}
```

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

The what string knows that `"some_non_existent_key"` is expected to be
contained in the `"grid"` key, which simplifies debugging to a great extent.
