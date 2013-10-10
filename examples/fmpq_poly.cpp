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

    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Tom Bachmann (C++ adaptation)

******************************************************************************/

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
