/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Tom Bachmann (C++ adaptation)

******************************************************************************/

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
    fmpz_mod_polyxx x(n);
    x.set_coeff(3, 5);
    x.set_coeff(0, 6);

    print(x);flint_printf("\n");
    print(x.sqr());flint_printf("\n");

    return 0;
}

