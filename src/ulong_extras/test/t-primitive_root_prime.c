/*
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(n_primitive_root_prime, state)
{
    int i, j;

    {
        ulong n[] = { UWORD(15), UWORD(21), UWORD(33), UWORD(45) };

        for (i = 0; i < 4; i++)
        {
            ulong nonresidue = n_quadratic_nonresidue(n[i]);
            if (n_jacobi_unsigned(nonresidue, n[i]) != -1)
                TEST_FUNCTION_FAIL("%wu is not a Jacobi nonresidue mod %wu\n", nonresidue, n[i]);
        }
    }

    for (i = 0; i < 100; i++)
    {
        ulong p = n_randtest_prime(state, 1);
        ulong pinv = n_preinvert_limb(p);

        n_factor_t factors;
        n_factor_init(&factors);
        n_factor(&factors, p - 1, 1);

        ulong root;
        if (p != 2)
        {
            ulong nonresidue = n_quadratic_nonresidue(p);
            if (n_jacobi_unsigned(nonresidue, p) != -1)
                TEST_FUNCTION_FAIL("%wu is not a quadratic nonresidue mod %wu\n", nonresidue, p);
        }

        root = n_primitive_root_prime(p);

        for (j = 0; j < factors.num; j++)
        {
            if (n_powmod2_preinv(root, (p-1) / factors.p[j], p, pinv) == 1)
                TEST_FUNCTION_FAIL("%wu ** (%wu / %wu) == 1 mod %wu\n", root, p-1, factors.p[j], p);
        }
    }

    TEST_FUNCTION_END(state);
}
