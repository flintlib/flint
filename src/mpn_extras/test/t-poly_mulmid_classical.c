/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_poly.h"

TEST_FUNCTION_START(_flint_mpn_poly_mulmid_classical, state)
{
    slong iter;

    for (iter = 0; iter < 5000 * flint_test_multiplier(); iter++)
    {
        slong n1, n2, s, flen, glen, nlo, nhi, i, len, st1, st2;
        nn_ptr f, g, r1, r2;
        fmpz_poly_t F, G, H;
        fmpz_t t;
        int sgn, squaring, karatsuba;

        if (n_randint(state, 8) == 0)
        {
            /* large limb counts: exercises the Toom-Cook base case of
               Karatsuba (unsigned inputs only) */
            n1 = 28 + n_randint(state, 12);
            n2 = 28 + n_randint(state, 12);
            if (n1 < n2) FLINT_SWAP(slong, n1, n2);
            sgn = 0;
            squaring = (n1 == n2) && n_randint(state, 2);
            karatsuba = 1 + n_randint(state, 2);
        }
        else
        {
            n1 = 1 + n_randint(state, 5);
            n2 = 1 + n_randint(state, n1);
            /* 0: unsigned, 1: two's complement, 2: sign-magnitude (classical only) */
            sgn = n_randint(state, 3);
            squaring = (n1 == n2) && n_randint(state, 2);
            /* 0: classical, 1: full Karatsuba, 2: windowed Karatsuba
               (both also with unbalanced limb sizes) */
            karatsuba = (sgn == 2) ? 0 : n_randint(state, 3);
        }
        s = n1 + n2 - 1 + n_randint(state, 3);

        if (karatsuba == 2 && n_randint(state, 2))
        {
            /* unbalanced shapes, and the middle product shape 2n-1 x n */
            glen = 1 + n_randint(state, 30);
            flen = glen + n_randint(state, 3 * glen + 1);
            if (n_randint(state, 3) == 0)
                flen = 2 * glen - 1;
            squaring = 0;
        }
        else
        {
            flen = 1 + n_randint(state, karatsuba ? 40 : 12);
            glen = squaring ? flen : 1 + n_randint(state, karatsuba ? 40 : 12);
        }
        len = flen + glen - 1;

        if (karatsuba == 1)
        {
            nlo = 0;
            nhi = len;
        }
        else if (karatsuba == 2)
        {
            int shape = n_randint(state, 5);
            if (shape == 0)         /* low */
            {
                nlo = 0;
                nhi = 1 + n_randint(state, len);
            }
            else if (shape == 1)    /* high */
            {
                nlo = n_randint(state, len);
                nhi = len;
            }
            else if (shape == 2 && flen >= glen)    /* middle */
            {
                nlo = glen - 1;
                nhi = flen;
            }
            else
            {
                nlo = n_randint(state, len);
                nhi = nlo + 1 + n_randint(state, len - nlo);
            }
        }
        else
        {
            nlo = n_randint(state, len);
            nhi = nlo + 1 + n_randint(state, len - nlo);
        }

        /* sign-magnitude coefficients have an extra (sign) limb */
        st1 = n1 + (sgn == 2);
        st2 = n2 + (sgn == 2);

        f = flint_malloc(sizeof(ulong) * st1 * flen);
        g = squaring ? f : flint_malloc(sizeof(ulong) * st2 * glen);
        r1 = flint_malloc(sizeof(ulong) * s * len);
        r2 = flint_malloc(sizeof(ulong) * s * len);

        for (i = 0; i < st1 * flen; i++)
            f[i] = n_randtest(state);
        if (!squaring)
            for (i = 0; i < st2 * glen; i++)
                g[i] = n_randtest(state);
        for (i = 0; i < s * len; i++)
            r1[i] = n_randtest(state);

        if (sgn == 2)
        {
            for (i = 0; i < flen; i++)
                f[i * st1] &= 1;
            if (!squaring)
                for (i = 0; i < glen; i++)
                    g[i * st2] &= 1;
        }

        if (karatsuba)
        {
            /* the intermediate sums must not overflow: clear the top bits */
            int norm1 = 2 + n_randint(state, 4);
            int norm2 = squaring ? norm1 : 2 + n_randint(state, 4);
            for (i = n1 - 1; i < n1 * flen; i += n1)
                f[i] = sgn ? (ulong) (((slong) f[i]) >> norm1) : (f[i] >> norm1);
            if (!squaring)
                for (i = n2 - 1; i < n2 * glen; i += n2)
                    g[i] = sgn ? (ulong) (((slong) g[i]) >> norm2) : (g[i] >> norm2);
            /* enough output limbs for the growth; n1 + n2 + 1 is required
               when the coefficients get extended */
            s = n1 + n2 + 1;
            r1 = flint_realloc(r1, sizeof(ulong) * s * len);
            r2 = flint_realloc(r2, sizeof(ulong) * s * len);
            if (karatsuba == 1 && n1 == n2)
                _flint_mpn_poly_mul_karatsuba(r1, f, flen, g, glen, n1, s, 1 + n_randint(state, 6), norm1 - 1, sgn);
            else if (karatsuba == 1)
                _flint_mpn_poly_mulmid_karatsuba(r1, f, flen, n1, norm1 - 1, g, glen, n2, norm2 - 1, 0, len, s, 1 + n_randint(state, 6), sgn);
            else
                _flint_mpn_poly_mulmid_karatsuba(r1, f, flen, n1, norm1 - 1, g, glen, n2, norm2 - 1, nlo, nhi, s, 1 + n_randint(state, 6), sgn);
        }
        else
        {
            _flint_mpn_poly_mulmid_classical(r1, f, flen, n1, g, glen, n2, nlo, nhi, s, sgn);
        }

        /* reference */
        fmpz_poly_init2(F, flen);
        fmpz_poly_init2(G, glen);
        fmpz_poly_init(H);
        fmpz_init(t);

        for (i = 0; i < flen; i++)
        {
            if (sgn == 2) { fmpz_set_ui_array(t, f + i * st1 + 1, n1); if (f[i * st1]) fmpz_neg(t, t); }
            else if (sgn) fmpz_set_signed_ui_array(t, f + i * n1, n1);
            else fmpz_set_ui_array(t, f + i * n1, n1);
            fmpz_poly_set_coeff_fmpz(F, i, t);
        }
        for (i = 0; i < glen; i++)
        {
            if (sgn == 2) { fmpz_set_ui_array(t, g + i * st2 + 1, n2); if (g[i * st2]) fmpz_neg(t, t); }
            else if (sgn) fmpz_set_signed_ui_array(t, g + i * n2, n2);
            else fmpz_set_ui_array(t, g + i * n2, n2);
            fmpz_poly_set_coeff_fmpz(G, i, t);
        }

        fmpz_poly_mul_classical(H, F, G);

        for (i = nlo; i < nhi; i++)
        {
            fmpz_poly_get_coeff_fmpz(t, H, i);
            fmpz_fdiv_r_2exp(t, t, s * FLINT_BITS);
            fmpz_get_ui_array(r2 + (i - nlo) * s, s, t);
        }

        if (mpn_cmp(r1, r2, (nhi - nlo) * s) != 0)
            TEST_FUNCTION_FAIL(
                    "n1 = %wd, n2 = %wd, s = %wd, flen = %wd, glen = %wd, nlo = %wd, nhi = %wd, sgn = %d, squaring = %d, karatsuba = %d\n"
                    "F = %{fmpz_poly}\nG = %{fmpz_poly}\n"
                    "r1 = %{ulong*}\n"
                    "r2 = %{ulong*}\n",
                    n1, n2, s, flen, glen, nlo, nhi, sgn, squaring, karatsuba, F, G,
                    r1, (nhi - nlo) * s, r2, (nhi - nlo) * s);

        fmpz_poly_clear(F);
        fmpz_poly_clear(G);
        fmpz_poly_clear(H);
        fmpz_clear(t);
        flint_free(f);
        if (!squaring)
            flint_free(g);
        flint_free(r1);
        flint_free(r2);
    }

    TEST_FUNCTION_END(state);
}
