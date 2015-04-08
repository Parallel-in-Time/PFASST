# Style Guide                                                                    {#page_style_guide}

## Naming Convention

In general, abbreviations should only be used if they are absolutely
unambigious and would save a lot of horizontal space.

Always use underscores (`_`) to separate words.

* __Files__

  To accomodate systems which do not distinguish between upper and
  lower case in file and path names, always use lowercase for file
  names and directories.

  Header files should be named after the main/central class that they
  contain.  For example, a class `MyClass` should be in a header file
  called `my_header.hpp`.  Separate words by underscores.

* __Namespaces__

  Should always be lowercase and not too long.

* __Class Names__

  Should start with a capital letter and use _CamelCase_.

  Abstract classes should start with a capital `I` denoting an
  _interface_.  For example, `ISweeper`.

* __Function Names__

  Should always be lowercase.

* __Variable Names__

  Should always be lowercase.  Constants should all be uppercase with
  undersocres as word separation.

* __Template Parameters__

  When template parameters refer to a type name, append a capital `T`
  to the name of the template parameter.  For example
  ~~~{.cpp}
  template<typename typeT>
  int my_func(std::vector<typeT> &vec);
  ~~~

  Acceptable exceptions to this rule are `scalar` (which denotes a
  spatial type, e.g. `complex<double>`) and `time` (which denotes a
  temporal type, e.g. `double`).


## Techniques

### Includes

System headers and those from third party libraries should be included
using angular brackets.  Headers from our library should be included
using double quotes.

### Define Guards

Use `#%define` guards in header files to prevent multiple inclusion.
The variable defined should be derived from the header file name
itself.  For example, a header file `my_header.hpp` in the folder
`PFASST/include/subfolder` should have the following define guard:

~~~{.cpp}
#ifndef PFASST__SUBFOLDER__MY_HEADER_HPP
#define PFASST__SUBFOLDER__MY_HEADER_HPP
// file content
#endif // PFASST__SUBFOLDER__MY_HEADER_HPP
~~~

Always include the library name and at least one subfolder (if
applicable) to prevent multiple defines in files with the same name
but in different folders.

### Encapsulation / Clean API

We want to provide a clean API, thus we use encapsulation wherever
possible.  For classes, we prefer private member/data variables and
provide access to them via getter/setter methods.

### C++ Feature Set

We agreed on using the "most common" C++11 features.  These are:

* `auto` keyword
* `shared_ptr` and `unique_ptr`
* type traits
* static asserts
* Rvalue references and move semantics
* `decltype`
* delegating constructor
* `constexpr`
* right angle bracket (`vector<vector<T> >` vs. `<vector<vector<T>>`)

@todo Check the list of used/required C++11 features and possibly extend it by a list of minimum
  compiler versions supporting the features.


## Documentation

Document your code.

We use [Doxygen] for generating the documentation webpage. You can use pretty much all the features,
Doxygen provides, including [MathJAX] for formulas and [Markdown] for easy text formatting.

It is advised to use Doxygen's special commands wherever possible to aid readability of the 
generated documentation.
Especially, one should use `@tparam <T> <description>`, `@param[<in>,<out>] <param> <description>`,
`@returns <description>`, `@throws <exception> <description>`.

There is an additional custom defined block available to mark documentation of internals.
Therefore, sourround the respective block with `internals` and `endinternals` (as Doxygen commands).


## Formatting

In short: no tabs, 2 spaces, sane indent-continuation.

We provide a configuration file for the code beautifier [uncrustify]
which does a pretty good job to ensure consistent formatting of the
source files.


[Doxygen]: http://www.stack.nl/~dimitri/doxygen/index.html
[MathJAX]: http://www.stack.nl/~dimitri/doxygen/manual/formulas.html
[Markdown]: http://www.stack.nl/~dimitri/doxygen/manual/markdown.html
[uncrustify]: http://uncrustify.sourceforge.net/
