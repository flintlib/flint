/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Tom Bachmann (C++ adaptation)

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Simple example demonstrating the use of the fmpq_polyxx module.
 */

#include "fmpq_polyxx.h"
#include <iostream>

using namespace flint;

int main(int argc, char* argv[])
{
    fmpq_polyxx f("2  1/2 3/5");
    fmpq_polyxx g("4  1/3 2 3/2 -1/2");

    std::cout << '(' << f.pretty("t") << ") * (" << g.pretty("t")
              << " = " << (f*g).pretty("t") << '\n';

    return 0;
}
