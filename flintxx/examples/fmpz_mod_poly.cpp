/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Tom Bachmann (C++ adaptation)

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Example program for the fmpz_mod_poly module.
*/

#include <cstdio>
#include <fmpz_mod_polyxx.h>

using namespace std;
using namespace flint;

int main(int argc, char* argv[])
{
    fmpzxx n(7);
    fmpz_modxx_ctx p(n);
    fmpz_mod_polyxx x(p);
    x.set_coeff(3, 5);
    x.set_coeff(0, 6);
    print(x);flint_printf("\n");
    print(x.sqr());flint_printf("\n");

    return 0;
}

