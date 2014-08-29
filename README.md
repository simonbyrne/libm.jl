# libm.jl

This is a very rough implementation of some of the standard math library
(libm) functions in pure julia. Much of it is based on
[openlibm](https://github.com/JuliaLang/openlibm), which in turn is based
around [fdlibm](http://www.netlib.org/fdlibm/).

You probably don't want to use this code directly, but you might find it
of interest to understand how math functions work under the hood.

[![Build Status](https://travis-ci.org/simonbyrne/libm.jl.svg?branch=master)](https://travis-ci.org/simonbyrne/libm.jl)
