/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

/*
    Assumes that p is an odd prime, and that 0 <= a < p.
    Returns 1 if a is a quadratic residue and 0 otherwise.
    Does not support aliasing.
 */
static int _fmpz_sqrtmod(mpz_t rop, const mpz_t a, const mpz_t p) 
{
    slong i, r, m, iter;
    mpz_t p1, k, exp, b, g, bpow, gpow;

    if (mpz_jacobi(a, p) == -1)
        return 0;

    if (flint_mpz_congruent_ui_p(p, 3, 4))
    {
        mpz_init(exp);
        flint_mpz_add_ui(exp, p, 1);
        mpz_tdiv_q_2exp(exp, exp, 2);
        mpz_powm(rop, a, exp, p);
        mpz_clear(exp);

        return 1;
    }

    if (flint_mpz_congruent_ui_p(p, 5, 8))
    {
        mpz_init(exp);
        mpz_init(g);
        mpz_init(k);
        flint_mpz_add_ui(exp, p, 3);
        mpz_tdiv_q_2exp(exp, exp, 3);
        mpz_powm(rop, a, exp, p);
        mpz_set_ui(k, 2);
        mpz_powm(g, rop, k, p);
        if (mpz_cmp(g, a) == 0)
        {
           mpz_clear(exp);
           mpz_clear(g);
           mpz_clear(k);
           return 1;
        }
        flint_mpz_sub_ui(exp, p, 1);
        mpz_tdiv_q_2exp(exp, exp, 2);
        mpz_powm(g, k, exp, p);
        mpz_mul(rop, rop, g);
        mpz_mod(rop, rop, p);
        mpz_clear(exp);
        mpz_clear(g);
        mpz_clear(k);
        return 1;
    }

    mpz_init(p1);
    mpz_init(k);
    mpz_init(exp);
    mpz_init(b);
    mpz_init(g);
    mpz_init(bpow);
    mpz_init(gpow);

    r = 0;
    flint_mpz_sub_ui(p1, p, 1);
    do {
        mpz_tdiv_q_2exp(p1, p1, 1);
        r++;
    } while (mpz_even_p(p1));

    mpz_powm(b, a, p1, p);

    for (flint_mpz_set_ui(k, 3); ; flint_mpz_add_ui(k, k, 2)) /* 2 is a quadratic residue mod p = 8k + 1 */
    {
        if (mpz_jacobi(k, p) == -1) break;
    }

    mpz_powm(g, k, p1, p);

    flint_mpz_add_ui(exp, p1, 1);
    mpz_tdiv_q_2exp(exp, exp, 1);
    mpz_powm(rop, a, exp, p);

    iter = r - 1; /* maximum number of iterations if p is prime */

    while (flint_mpz_cmp_ui(b, 1))
    {
        mpz_set(bpow, b);
        m = 0;
        do
        {
            mpz_mul(bpow, bpow, bpow);
            mpz_mod(bpow, bpow, p);
            m++;
        } while (m < r && flint_mpz_cmp_ui(bpow, 1));

        mpz_set(gpow, g);
        for (i = 1; i < r - m; i++)
        {
            mpz_mul(gpow, gpow, gpow);
            mpz_mod(gpow, gpow, p);
        }

        mpz_mul(rop, rop, gpow);
        mpz_mod(rop, rop, p);

        mpz_mul(g, gpow, gpow);
        mpz_mod(g, g, p);

        mpz_mul(b, b, g);
        mpz_mod(b, b, p);

        r = m;
        if (iter-- == 0) /* too many iterations, p is not prime */
        {
            mpz_set_ui(rop, 0);
            break;
        }
    }

    mpz_clear(p1);
    mpz_clear(k);
    mpz_clear(exp);
    mpz_clear(b);
    mpz_clear(g);
    mpz_clear(bpow);
    mpz_clear(gpow);

    return mpz_sgn(rop) ? 1 : 0;
}

int fmpz_sqrtmod(fmpz_t b, const fmpz_t a, const fmpz_t p)
{
    if (b == a || b == p)
    {
        int ans;
        fmpz_t t;

        fmpz_init(t);
        ans = fmpz_sqrtmod(t, a, p);
        fmpz_swap(b, t);
        fmpz_clear(t);
        return ans;
    }

    fmpz_mod(b, a, p);

    if (fmpz_cmp_ui(b, 1) <= 0)
    {
        return 1;
    }

    if (!COEFF_IS_MPZ(*p))  /* p, and b are small */
    {
        mp_limb_t ans;

        ans = n_sqrtmod(*b, *p);
        if (ans)
            fmpz_set_ui(b, ans);
        return ans != 0;
    }
    else  /* p is large */
    {
        int ans;
        mpz_t t;
        __mpz_struct *bptr;
        
        bptr = _fmpz_promote_val(b);

        mpz_init(t);
        ans = _fmpz_sqrtmod(t, bptr, COEFF_TO_PTR(*p));
        mpz_swap(bptr, t);
        mpz_clear(t);

        _fmpz_demote_val(b);

        return ans;
    }
}

