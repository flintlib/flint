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
    Simple example demonstrating the use of the fmpz_poly_q module.
 */

#include <iostream>
#include "fmpz_poly_qxx.h"

using namespace flint;
using namespace std;

int main(int argc, char* argv[])
{
    fmpz_poly_qxx f("2  1 3/1  2"), g("1  3/2  2 7");
    std::cout << f.pretty("t") << " * " << g.pretty("t")
              << " = " << (f*g).pretty("t") << '\n';
    return 0;
}

