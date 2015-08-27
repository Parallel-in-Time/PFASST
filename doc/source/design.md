# Design                                                                              {#page_design}

## Technical

### Default Constructors and Operator Overloading

To ensure the compiler is generating the default constructors (incl. copy- and move-ctors) as well
as the default copy- and move-assignment operators for our classes, these are explicitly stated as
`= default`.
