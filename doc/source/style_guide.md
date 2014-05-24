# Style Guide

To ease development we agree on a basic code style described as follows.

## Naming Convention

In general, abbreviations should only be used if they are absolutely unambigious and would save a
lot of horizontal space.
In case one hesitates whether to abbreviate a variable or not, rather not do it!

Always use underscores (`_`) as word separation.

* __Files__

  Due to systems which do not distinguish between upper and lower case in file and path names,
  always use lowercase for file names and directories.

  Header files should be named after the main/central class contained in.
  For example, a class `MyClass` should be in a header file called `my_header.hpp`.
  Separate words by underscores.

* __Namespaces__

  Should always be lowercase and not too long.

* __Class Names__

  Should start with a capital letter and use _CamelCase_.
  
  Abstract classes should start with a capital `I` denoting an _interface_. E.g. `ISweeper` in 
  `i_sweeper.hpp`.

* __Function Names__

  Should only consist of lowercase letters and underscores separating words.

* __Variable Names__

  Should always be lowercase using underscores as word separation.
  Member variables of classes should start with `m_`.
  Constants should all be uppercase with undersocres as word separation.

* __Template Parameter__

  In case, most template parameters denote a certain type, append a capital `T` to the name of the
  template parameter. E.g.
  ~~~{.cpp}
  template<typename typeT>
  int my_func(std::vector<typeT> &vec);
  ~~~


## Techniques

### Includes

System headers and those from third party libraries should be included using angular brackets.
Headers from our library should be included using double quotes.

For headers included via angular brackets, the preprocessor will search through all paths given via
the `-I` parameter.
While for includes with double quotes it will only look in the directory of the header including
the included file.
Only when there is no such header in there, the preprocessor will go on with the behaviour for 
angular brackets.

### Define Guards

Use `#define` guards in header files to prevent multiple inclusion.
The variable defined should be derived from the header file name itself.
For example, a header file `my_header.hpp` in the folder `PFASST/src/subfolder` should have the 
following define guard:

~~~{.cpp}
#ifndef PFASST__SUBFOLDER__MY_HEADER_HPP
#define PFASST__SUBFOLDER__MY_HEADER_HPP
// file content
#endif // PFASST__SUBFOLDER__MY_HEADER_HPP
~~~

Always include the library name and at least one subfolder (if applicable) to prevent possible 
multiple defines in files with the same name but in different folders.

### Encapsulation / Clean API

We want to provide a clean API, thus we use encapsulation wherever possible.
For classes, we prefer private member/data variables and provide access to them via getter/setter
methods.

### C++ Feature Set

We agreed on using the "most common" C++11 features.
These are:

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

We use [Doxygen] for generating the documentation webpage.
You can use pretty much all the features, Doxygen provides, including [MathJAX] for formulas and 
[Markdown] for easy text formatting.


## Formatting

In short: no tabs, 2 spaces, sane indent-continuation.

We provide a configuration file for the code beautifier [uncrustify] which does a pretty good job
to ensure consistent formatting of the source files.


[Doxygen]: http://www.stack.nl/~dimitri/doxygen/index.html
[MathJAX]: http://www.stack.nl/~dimitri/doxygen/manual/formulas.html
[Markdown]: http://www.stack.nl/~dimitri/doxygen/manual/markdown.html
[uncrustify]: http://uncrustify.sourceforge.net/
