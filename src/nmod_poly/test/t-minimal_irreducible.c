/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"

int
nmod_poly_irreducible_binomial_naive(nmod_poly_t res, ulong n)
{
    ulong p = res->mod.n;
    ulong a;

    if (n == 0)
        return 0;

    nmod_poly_zero(res);
    nmod_poly_set_coeff_ui(res, n, 1);

    for (a = 1; a < p; a++)
    {
        if (a > 5)  /* fine for the test code */
            return 0;

        nmod_poly_set_coeff_ui(res, 0, a);
        if (nmod_poly_is_irreducible(res))
            return 1;
    }

    return 0;
}

int
nmod_poly_irreducible_trinomial_naive(nmod_poly_t res, ulong n)
{
    ulong p = res->mod.n;
    ulong a, b, k;

    if (n <= 1)
        return 0;

    /* Iterate over all x^n + a x^k + b */
    for (a = 1; a < p; a++)
    {
        for (b = 1; b < p; b++)
        {
            for (k = 1; k < n; k++)
            {
                if (k > 2 * n)
                    continue;

                nmod_poly_zero(res);
                nmod_poly_set_coeff_ui(res, n, 1);
                nmod_poly_set_coeff_ui(res, 0, b);
                nmod_poly_set_coeff_ui(res, k, a);
                if (nmod_poly_is_irreducible(res))
                    return 1;
            }
        }
    }

    return 0;
}

int
nmod_poly_irreducible_tetranomial_naive(nmod_poly_t res, ulong n)
{
    ulong p = res->mod.n;
    ulong a, b, c, k, l;

    if (n <= 2)
        return 0;

    /* Iterate over all x^n + a x^k + b x^l + c */
    for (a = 1; a < p; a++)
    {
        for (b = 1; b < p; b++)
        {
            for (c = 1; c < p; c++)
            {
                for (k = 2; k < n; k++)
                {
                    for (l = 1; l < k; l++)
                    {
                        nmod_poly_zero(res);
                        nmod_poly_set_coeff_ui(res, n, 1);
                        nmod_poly_set_coeff_ui(res, k, a);
                        nmod_poly_set_coeff_ui(res, l, b);
                        nmod_poly_set_coeff_ui(res, 0, c);
                        if (nmod_poly_is_irreducible(res))
                            return 1;
                    }
                }
            }
        }
    }

    return 0;
}

int
nmod_poly_irreducible_pentanomial_naive(nmod_poly_t res, ulong n)
{
    ulong k, l, m;

    if (n <= 3)
        return 0;

    /* Iterate over all x^n + x^k + x^l + x^m + 1 */
    for (k = 3; k < n; k++)
    {
        for (l = 2; l < k; l++)
        {
            for (m = 1; m < l; m++)
            {
                nmod_poly_zero(res);
                nmod_poly_set_coeff_ui(res, n, 1);
                nmod_poly_set_coeff_ui(res, k, 1);
                nmod_poly_set_coeff_ui(res, l, 1);
                nmod_poly_set_coeff_ui(res, m, 1);
                nmod_poly_set_coeff_ui(res, 0, 1);
                if (nmod_poly_is_irreducible(res))
                    return 1;
            }
        }
    }

    return 0;
}

int nmod_poly_irreducible_binomial(nmod_poly_t res, ulong n);
int nmod_poly_irreducible_trinomial(nmod_poly_t res, ulong n);
int nmod_poly_irreducible_tetranomial(nmod_poly_t res, ulong n);
int nmod_poly_irreducible_pentanomial(nmod_poly_t res, ulong n);

int nontrinomials_tab[][3] = {
    {2, 1, 1},
    {2, 8, 5},
    {2, 13, 5},
    {2, 16, 5},
    {2, 19, 5},
    {2, 24, 5},
    {2, 26, 5},
    {2, 27, 5},
    {2, 32, 5},
    {2, 37, 5},
    {2, 38, 5},
    {2, 40, 5},
    {2, 43, 5},
    {2, 45, 5},
    {2, 48, 5},
    {2, 50, 5},
    {2, 51, 5},
    {2, 53, 5},
    {2, 56, 5},
    {2, 59, 5},
    {2, 61, 5},
    {2, 64, 5},
    {2, 67, 5},
    {2, 69, 5},
    {2, 70, 5},
    {2, 72, 5},
    {2, 75, 5},
    {2, 77, 5},
    {2, 78, 5},
    {2, 80, 5},
    {2, 82, 5},
    {2, 83, 5},
    {2, 85, 5},
    {2, 88, 5},
    {2, 91, 5},
    {2, 96, 5},
    {2, 99, 5},
    {2, 101, 5},
    {2, 104, 5},
    {2, 107, 5},
    {2, 109, 5},
    {2, 112, 5},
    {2, 114, 5},
    {2, 115, 5},
    {2, 116, 5},
    {2, 117, 5},
    {2, 120, 5},
    {2, 122, 5},
    {2, 125, 5},
    {2, 128, 5},
    {2, 131, 5},
    {2, 133, 5},
    {2, 136, 5},
    {2, 138, 5},
    {2, 139, 5},
    {2, 141, 5},
    {2, 143, 5},
    {2, 144, 5},
    {2, 149, 5},
    {2, 152, 5},
    {2, 157, 5},
    {2, 158, 5},
    {2, 160, 5},
    {2, 163, 5},
    {2, 164, 5},
    {2, 165, 5},
    {2, 168, 5},
    {2, 171, 5},
    {2, 173, 5},
    {2, 176, 5},
    {2, 179, 5},
    {2, 181, 5},
    {2, 184, 5},
    {2, 187, 5},
    {2, 188, 5},
    {2, 189, 5},
    {2, 190, 5},
    {2, 192, 5},
    {2, 195, 5},
    {2, 197, 5},
    {2, 200, 5},
    {2, 203, 5},
    {2, 205, 5},
    {2, 206, 5},
    {2, 208, 5},
    {2, 211, 5},
    {2, 213, 5},
    {2, 216, 5},
    {2, 219, 5},
    {2, 221, 5},
    {2, 222, 5},
    {2, 224, 5},
    {2, 226, 5},
    {2, 227, 5},
    {2, 229, 5},
    {2, 230, 5},
    {2, 232, 5},
    {2, 235, 5},
    {2, 237, 5},
    {2, 240, 5},
    {2, 243, 5},
    {2, 245, 5},
    {2, 246, 5},
    {2, 248, 5},
    {2, 251, 5},
    {2, 254, 5},
    {2, 256, 5},
    {2, 259, 5},
    {2, 261, 5},
    {2, 262, 5},
    {2, 264, 5},
    {2, 267, 5},
    {2, 269, 5},
    {2, 272, 5},
    {2, 275, 5},
    {2, 277, 5},
    {2, 280, 5},
    {2, 283, 5},
    {2, 285, 5},
    {2, 288, 5},
    {2, 290, 5},
    {2, 291, 5},
    {2, 293, 5},
    {2, 296, 5},
    {2, 298, 5},
    {2, 299, 5},
    {2, 301, 5},
    {2, 304, 5},
    {2, 306, 5},
    {2, 307, 5},
    {2, 309, 5},
    {2, 311, 5},
    {2, 312, 5},
    {2, 315, 5},
    {2, 317, 5},
    {2, 320, 5},
    {2, 323, 5},
    {2, 325, 5},
    {2, 326, 5},
    {2, 328, 5},
    {2, 331, 5},
    {2, 334, 5},
    {2, 335, 5},
    {2, 336, 5},
    {2, 338, 5},
    {2, 339, 5},
    {2, 341, 5},
    {2, 344, 5},
    {2, 347, 5},
    {2, 349, 5},
    {2, 352, 5},
    {2, 355, 5},
    {2, 356, 5},
    {2, 357, 5},
    {2, 360, 5},
    {2, 361, 5},
    {2, 363, 5},
    {2, 365, 5},
    {2, 368, 5},
    {2, 371, 5},
    {2, 373, 5},
    {2, 374, 5},
    {2, 376, 5},
    {2, 379, 5},
    {2, 381, 5},
    {2, 384, 5},
    {2, 387, 5},
    {2, 389, 5},
    {2, 392, 5},
    {2, 395, 5},
    {2, 397, 5},
    {2, 398, 5},
    {2, 400, 5},
    {2, 403, 5},
    {2, 405, 5},
    {2, 408, 5},
    {2, 410, 5},
    {2, 411, 5},
    {2, 413, 5},
    {2, 416, 5},
    {2, 419, 5},
    {2, 421, 5},
    {2, 424, 5},
    {2, 427, 5},
    {2, 429, 5},
    {2, 430, 5},
    {2, 432, 5},
    {2, 434, 5},
    {2, 435, 5},
    {2, 437, 5},
    {2, 440, 5},
    {2, 442, 5},
    {2, 443, 5},
    {2, 445, 5},
    {2, 448, 5},
    {2, 451, 5},
    {2, 452, 5},
    {2, 453, 5},
    {2, 454, 5},
    {2, 456, 5},
    {2, 459, 5},
    {2, 461, 5},
    {2, 464, 5},
    {2, 466, 5},
    {2, 467, 5},
    {2, 469, 5},
    {2, 472, 5},
    {2, 475, 5},
    {2, 477, 5},
    {2, 480, 5},
    {2, 482, 5},
    {2, 483, 5},
    {2, 485, 5},
    {2, 488, 5},
    {2, 491, 5},
    {2, 493, 5},
    {2, 496, 5},
    {2, 499, 5},
    {3, 1, 1},
    {3, 2, 2},
    {3, 49, 4},
    {3, 57, 4},
    {3, 65, 4},
    {3, 68, 4},
    {3, 75, 4},
    {3, 98, 4},
    {3, 105, 4},
    {3, 123, 4},
    {3, 129, 4},
    {3, 130, 4},
    {3, 132, 4},
    {3, 149, 4},
    {3, 161, 4},
    {3, 175, 4},
    {3, 189, 4},
    {3, 197, 4},
    {3, 207, 4},
    {3, 212, 4},
    {3, 213, 4},
    {3, 221, 4},
    {3, 223, 4},
    {3, 231, 4},
    {3, 233, 4},
    {3, 264, 4},
    {3, 267, 4},
    {3, 276, 4},
    {3, 281, 4},
    {3, 292, 4},
    {3, 297, 4},
    {3, 298, 4},
    {3, 303, 4},
    {3, 309, 4},
    {3, 311, 4},
    {3, 319, 4},
    {3, 332, 4},
    {3, 343, 4},
    {3, 391, 4},
    {3, 394, 4},
    {3, 397, 4},
    {3, 401, 4},
    {3, 404, 4},
    {3, 405, 4},
    {3, 411, 4},
    {3, 415, 4},
    {3, 426, 4},
    {3, 435, 4},
    {3, 436, 4},
    {3, 437, 4},
    {3, 439, 4},
    {3, 441, 4},
    {3, 442, 4},
    {3, 453, 4},
    {3, 459, 4},
    {3, 463, 4},
    {3, 489, 4},
    {3, 497, 4},
    {5, 1, 1},
    {5, 2, 2},
    {5, 4, 2},
    {5, 8, 2},
    {5, 16, 2},
    {5, 32, 2},
    {5, 35, 4},
    {5, 64, 2},
    {5, 70, 4},
    {5, 123, 4},
    {5, 125, 4},
    {5, 128, 2},
    {5, 140, 4},
    {5, 181, 4},
    {5, 191, 4},
    {5, 209, 4},
    {5, 213, 4},
    {5, 219, 4},
    {5, 237, 4},
    {5, 249, 4},
    {5, 250, 4},
    {5, 253, 4},
    {5, 256, 2},
    {5, 265, 4},
    {5, 273, 4},
    {5, 280, 4},
    {5, 285, 4},
    {5, 307, 4},
    {5, 319, 4},
    {5, 345, 4},
    {5, 351, 4},
    {5, 362, 4},
    {5, 375, 4},
    {5, 382, 4},
    {5, 407, 4},
    {5, 413, 4},
    {5, 415, 4},
    {5, 418, 4},
    {5, 421, 4},
    {5, 429, 4},
    {5, 433, 4},
    {5, 441, 4},
    {5, 445, 4},
    {5, 457, 4},
    {5, 461, 4},
    {5, 498, 4},
    {5, 500, 4},
    {7, 1, 1},
    {7, 2, 2},
    {7, 3, 2},
    {7, 6, 2},
    {7, 9, 2},
    {7, 18, 2},
    {7, 27, 2},
    {7, 54, 2},
    {7, 81, 2},
    {7, 124, 4},
    {7, 162, 2},
    {7, 163, 4},
    {7, 243, 2},
    {7, 268, 4},
    {7, 301, 4},
    {7, 359, 4},
    {7, 364, 4},
    {7, 372, 4},
    {7, 385, 4},
    {7, 389, 4},
    {7, 407, 4},
    {7, 431, 4},
    {7, 476, 4},
    {7, 486, 2},
    {7, 489, 4},
    {11, 1, 1},
    {11, 2, 2},
    {11, 5, 2},
    {11, 10, 2},
    {11, 25, 2},
    {11, 50, 2},
    {11, 125, 2},
    {11, 219, 4},
    {11, 250, 2},
    {11, 291, 4},
    {11, 467, 4},
    {13, 1, 1},
    {13, 2, 2},
    {13, 3, 2},
    {13, 4, 2},
    {13, 6, 2},
    {13, 8, 2},
    {13, 9, 2},
    {13, 12, 2},
    {13, 16, 2},
    {13, 18, 2},
    {13, 24, 2},
    {13, 27, 2},
    {13, 32, 2},
    {13, 36, 2},
    {13, 48, 2},
    {13, 54, 2},
    {13, 64, 2},
    {13, 72, 2},
    {13, 81, 2},
    {13, 96, 2},
    {13, 108, 2},
    {13, 128, 2},
    {13, 144, 2},
    {13, 162, 2},
    {13, 192, 2},
    {13, 216, 2},
    {13, 243, 2},
    {13, 256, 2},
    {13, 288, 2},
    {13, 324, 2},
    {13, 384, 2},
    {13, 432, 2},
    {13, 486, 2},
    {17, 1, 1},
    {17, 2, 2},
    {17, 4, 2},
    {17, 8, 2},
    {17, 16, 2},
    {17, 32, 2},
    {17, 64, 2},
    {17, 128, 2},
    {17, 231, 4},
    {17, 256, 2},
    {17, 375, 4},
    {19, 1, 1},
    {19, 2, 2},
    {19, 3, 2},
    {19, 6, 2},
    {19, 9, 2},
    {19, 18, 2},
    {19, 27, 2},
    {19, 54, 2},
    {19, 81, 2},
    {19, 162, 2},
    {19, 243, 2},
    {257, 1, 1},
    {257, 2, 2},
    {257, 4, 2},
    {257, 8, 2},
    {257, 16, 2},
    {257, 32, 2},
    {257, 64, 2},
    {367, 1, 1},
    {367, 2, 2},
    {367, 3, 2},
    {367, 6, 2},
    {367, 9, 2},
    {367, 18, 2},
    {367, 27, 2},
    {367, 54, 2},
    {367, 61, 2},
    {367, 81, 2},
    {1709, 1, 1},
    {1709, 2, 2},
    {1709, 4, 2},
    {1709, 7, 2},
    {1709, 8, 2},
    {1709, 14, 2},
    {1709, 16, 2},
    {1709, 28, 2},
    {1709, 32, 2},
    {1709, 49, 2},
    {1709, 56, 2},
    {1709, 61, 2},
    {1709, 64, 2},
    {1709, 98, 2},
    {13037, 1, 1},
    {13037, 2, 2},
    {13037, 4, 2},
    {13037, 8, 2},
    {13037, 16, 2},
    {13037, 32, 2},
    {13037, 64, 2},
    {15161, 1, 1},
    {15161, 2, 2},
    {15161, 4, 2},
    {15161, 5, 2},
    {15161, 8, 2},
    {15161, 10, 2},
    {15161, 16, 2},
    {15161, 20, 2},
    {15161, 25, 2},
    {15161, 32, 2},
    {15161, 40, 2},
    {15161, 50, 2},
    {15161, 64, 2},
    {15161, 80, 2},
    {15161, 100, 2},
    {69337, 1, 1},
    {69337, 2, 2},
    {69337, 3, 2},
    {69337, 4, 2},
    {69337, 6, 2},
    {69337, 8, 2},
    {69337, 9, 2},
    {69337, 12, 2},
    {69337, 16, 2},
    {69337, 18, 2},
    {69337, 24, 2},
    {69337, 27, 2},
    {69337, 32, 2},
    {69337, 36, 2},
    {69337, 48, 2},
    {69337, 54, 2},
    {69337, 64, 2},
    {69337, 72, 2},
    {69337, 81, 2},
    {69337, 96, 2},
    {438467, 1, 1},
    {438467, 2, 2},
    {438467, 7, 2},
    {438467, 14, 2},
    {438467, 49, 2},
    {438467, 98, 2},
    {0, 0, 0},
};


TEST_FUNCTION_START(nmod_poly_minimal_irreducible, state)
{
    ulong p, n, nmax;
    slong c1, c2, j, i, pi;
    nmod_poly_t f, g;
    int res1, res2;

    ulong ps[] = { 2, 3, 5, 7, 11, 17, 19, 29, 31, 37, 10007, UWORD(2147483659),
#if FLINT_BITS == 64
        UWORD(9223372036854775837),
#endif
        0 };

    for (i = 0; nontrinomials_tab[i][0] != 0; i++)
    {
        p = nontrinomials_tab[i][0];
        n = nontrinomials_tab[i][1];
        c2 = nontrinomials_tab[i][2];

        if (n > FLINT_MAX(50, flint_test_multiplier() * 10))
            continue;

        /* flint_printf("%wu %wu\n", p, n); */

        nmod_poly_init(f, p);
        nmod_poly_minimal_irreducible(f, n);
        c1 = 0;
        for (j = 0; j < f->length; j++)
            c1 += (f->coeffs[j] != 0);

        if (c1 != c2)
        {
            flint_printf("FAIL (wrong number of terms): p = %wu n = %wu  c1 = %wd c2 = %wd f = %{nmod_poly}\n",
                        p, n, c1, c2, f);
            flint_abort();
        }

        nmod_poly_clear(f);
    }

    for (pi = 0; (p = ps[pi]) != 0; pi++)
    {
        nmod_poly_init(f, p);
        nmod_poly_init(g, p);

        if (p == 2)
            nmax = 200;
        else if (p <= 5)
            nmax = 100;
        else
            nmax = 10;

        nmax = FLINT_MIN(nmax, FLINT_MAX(50, flint_test_multiplier() * 10));

        /* flint_printf("%wu\n", p); */

        for (n = 0; n <= nmax; n++)
        {
            /* flint_printf("n = %wd\n", n); */

            res1 = nmod_poly_irreducible_binomial(f, n);
            res2 = nmod_poly_irreducible_binomial_naive(g, n);

            if (res1 != res2 || (res1 && !nmod_poly_equal(f, g)))
            {
                flint_printf("FAIL (binomial): p = %wu n = %wu  res1 = %d res2 = %d\n", p, n, res1, res2);
                flint_abort();
            }

            res1 = nmod_poly_irreducible_trinomial(f, n);
            res2 = nmod_poly_irreducible_trinomial_naive(g, n);

            if (res1 != res2)
            {
                flint_printf("FAIL (trinomial): p = %wu n = %wu  res1 = %d res2 = %d\n", p, n, res1, res2);
                flint_abort();
            }

            res1 = nmod_poly_irreducible_tetranomial(f, n);
            res2 = nmod_poly_irreducible_tetranomial_naive(g, n);

            if (res1 != res2)
            {
                flint_printf("FAIL (tetranomial): p = %wu n = %wu  res1 = %d res2 = %d\n", p, n, res1, res2);
                flint_abort();
            }

            if (p == 2)
            {
                res1 = nmod_poly_irreducible_pentanomial(f, n);
                res2 = nmod_poly_irreducible_pentanomial_naive(g, n);

                if (res1 != res2)
                {
                    flint_printf("FAIL (pentanomial): p = %wu n = %wu  res1 = %d res2 = %d\n", p, n, res1, res2);
                    flint_abort();
                }
            }
        }

        nmod_poly_clear(f);
        nmod_poly_clear(g);
    }

    TEST_FUNCTION_END(state);
}
