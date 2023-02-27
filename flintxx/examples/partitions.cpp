/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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
