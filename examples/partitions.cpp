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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

/*
    FLINT program for demonstrating the Integer Partition function.
*/

#include <cstdlib>
#include <iostream>
#include "arithxx.h"
#include "fmpzxx.h"

using namespace flint;
using namespace std;

int
main(int argc, char * argv[])
{
    if (argc != 2)
    {
        std::cerr << "usage: partitions n\n";
        return 1;
    }

    ulong n;
    flint_sscanf(argv[1], "%wu", &n);

    std::cout << "p(" << n << ") =\n" << number_of_partitions(n) << '\n';

    return 0;
}
