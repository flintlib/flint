/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Tom Bachmann (C++ adaptation)

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Demo FLINT program for balanced multimodular reduction and
    reconstruction using the Chinese Remainder Theorem.
*/

#include <iostream>
#include <vector>
#include "fmpzxx.h"
#include "ulong_extras.h"

using namespace flint;

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        std::cerr << "Syntax: crt <integer> <num_primes>\n";
        return 1;
    }

    slong num_primes = atoi(argv[2]);

    if (num_primes < 1)
    {
        std::cerr << "Requires num_primes >= 1\n";
        return 2;
    }

    fmpzxx x(argv[1]);

    std::vector<mp_limb_t> primes(num_primes), residues(num_primes);
    primes[0] = 2;
    for (unsigned i = 1; i < num_primes; i++)
        primes[i] = n_nextprime(primes[i-1], 0);

    fmpz_combxx comb(primes);
    multi_mod(residues, x, comb);

    for (unsigned i = 0; i < num_primes; i++)
        std::cout << "residue mod " << primes[i]
                  << " = " << residues[i] << '\n';

    std::cout << "reconstruction = " << multi_CRT(residues, comb, true)
              << '\n';

    return 0;
}
