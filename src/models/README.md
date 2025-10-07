# Connecting user-defined libraries to the `model` class

The outlined instructions present the steps necessary to compile the
user-defined model library and connect it to the `model` class __only when__
the compiler option is specified.

## Step 1: Creating a cmake compiler option

In the `src/models/CMakeLists.txt`, after the `model_library` is initialized, there are
conditional statements which depend on the value of `MODEL_TYPE` string
variable.
In order to create a new compile option, one should add a new condition such as

```
...
elseif("${MODEL_TYPE}" STREQUAL "<MODEL_NAME>")
...
```
Then, during compilation, if `<MODEL_NAME>` is specified in the `-DMODEL_TYPE`
flag ([as shown here](../../README.md)), only the user-specified library is
compiled.

Inside the conditional block, the order of the cmake commands will generally
be the same:
* Define an option variable that will be available to the compiler during the
  creation of the `model_library`:
```
...
target_compile_difinitions(model_library PUBLIC "<OPTION_NAME>")
...
```
Typically, `<OPTION_NAME>` should be (at least partially) the same as
`<MODEL_NAME>` and in all-caps.

* Add the subdirectory with the source files and `CMakeLists.txt` which are
used to compile the user-defined model library:

```
...
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/<library>)
...
```

* Include the directory that contains the header files for the user-defined
model library:

```
...
target_include_directories(model_library PUBLIC
                           ${INCLUDE_FRUSA_MODELS}/<library>)
...
```

* Finally, link the user-defined library to the main `model_library`

```
...
target_link_libraries(model_library PRIVATE <library_name>)
...
```
Note that `<library_name>` must be the same as defined in
`./<library>/CMakeLists.txt`!

## Step 2: Connecting the relevant libraries to `model` header

Now that there is an option for compiling the user-defined library, one can
call it in the header file in `frusa/include/models/model.h`.
In the header add the following lines on the top of the file (below previous
model ):

```
...
#elif defined <OPTION_NAME>
#include <library_headers>
...
using namespace <library_space>
...
```
Here, `<library_headers>` are all of the relevant headers for the user-defined
library.

__NOTE:__ This is not an optimal method since this space will grow fast with
more user-defined models.
The goal in the future is to generate a single header which is included in the
`model` header.

Typically, the functions in the user-defined library are contained in a
namespace `<library_space>`.
In c++ `using namespace ...` is not the best practice, so this may be changed
in the future.
