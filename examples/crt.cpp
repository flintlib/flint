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
    Copyright (C) 2013 Tom Bachmann (C++ adaptation)

******************************************************************************/

/*
    Demo FLINT program for incremental multimodular reduction and
    reconstruction using the Chinese Remainder Theorem.
*/

#include <iostream>
#include "fmpzxx.h"
#include "ulong_extras.h"

using namespace flint;

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        flint_printf("Syntax: crt <integer>\n");
        return EXIT_FAILURE;
    }
    
    fmpzxx x(argv[1]);
    slong bit_bound = bits(x) + 2;

    fmpzxx y(0);
    fmpzxx prod(1);

    mp_limb_t prime = 0;
    for (unsigned i = 0; bits(prod) < bit_bound; i++)
    {
        prime = n_nextprime(prime, 0);

        ulong res = (x % prime).to<ulong>();
        y = y.CRT(prod, res, prime, true);

        std::cout << "residue mod " << prime << " = " << res;
        std::cout << "; reconstruction = " << y << '\n';

        prod *= prime;
    }

    return 0;
}

