/*
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint/fmpz.h"
#include "flint/profiler.h"

typedef struct
{
   slong bits;
}
info_t;

static ulong _fmpz_gcd_big_small_old(const fmpz_t g, ulong h)
{
    __mpz_struct * z = COEFF_TO_PTR(*g);

    return n_gcd(mpn_mod_1(z->_mp_d, FLINT_ABS(z->_mp_size), h), h);
}

static ulong _fmpz_gcd_small_old(const fmpz_t g, ulong h)
{
    if (!COEFF_IS_MPZ(*g))
        return n_gcd(FLINT_ABS(*g), h);
    else
        return _fmpz_gcd_big_small_old(g, h);
}

static void
fmpz_gcd3_small_old(fmpz_t res, const fmpz_t a, const fmpz_t b, ulong c)
{
    if (c <= 1)
    {
        if (c == 1)
            fmpz_one(res);
        else
            fmpz_gcd(res, a, b);
    }
    else
    {
        if (!COEFF_IS_MPZ(*a))
        {
            c = n_gcd(FLINT_ABS(*a), c);
            if (c != 1)
                c = _fmpz_gcd_small_old(b, c);
        }
        else
        {
            c = _fmpz_gcd_small_old(b, c);
            if (c != 1)
                c = _fmpz_gcd_big_small_old(a, c);
        }

        fmpz_set_ui(res, c);
    }
}

void
fmpz_gcd3_old(fmpz_t res, const fmpz_t a, const fmpz_t b, const fmpz_t c)
{
    if (!COEFF_IS_MPZ(*a))
    {
        fmpz_gcd3_small_old(res, b, c, FLINT_ABS(*a));
    }
    else if (!COEFF_IS_MPZ(*b))
    {
        fmpz_gcd3_small_old(res, a, c, FLINT_ABS(*b));
    }
    else if (!COEFF_IS_MPZ(*c))
    {
        fmpz_gcd3_small_old(res, a, b, FLINT_ABS(*c));
    }
    else
    {
        /* Three-way mpz_gcd. */
        __mpz_struct *rp, *ap, *bp, *cp, *tp;
        mp_size_t an, bn, cn, mn;

        /* If res is small, it cannot be aliased with a, b, c, so promoting is fine. */
        rp = _fmpz_promote(res);

        ap = COEFF_TO_PTR(*a);
        bp = COEFF_TO_PTR(*b);
        cp = COEFF_TO_PTR(*c);

        an = FLINT_ABS(ap->_mp_size);
        bn = FLINT_ABS(bp->_mp_size);
        cn = FLINT_ABS(cp->_mp_size);

        /* Select c to be the largest operand; we do the smaller gcd first. */
        mn = FLINT_MAX(FLINT_MAX(an, bn), cn);
        tp = cp;
        if (mn != cn)
        {
            if (mn == an)
            {
                cp = ap;
                ap = tp;
            }
            else
            {
                cp = bp;
                bp = tp;
            }

            cn = mn;
        }

        /* Handle aliasing */
        if (rp == cp)
        {
            mpz_t t;
            TMP_INIT;
            TMP_START;

            /* It would be more efficient to allocate temporary space for
               gcd(a, b), but we can't be sure that mpz_gcd never attempts
               to reallocate the output. */
            t->_mp_d = TMP_ALLOC(sizeof(mp_limb_t) * cn);
            t->_mp_size = t->_mp_alloc = cn;
            flint_mpn_copyi(t->_mp_d, cp->_mp_d, cn);

            mpz_gcd(rp, ap, bp);
            if (mpz_cmpabs_ui(rp, 1) != 0)
                mpz_gcd(rp, rp, t);

            TMP_END;
        }
        else
        {
            mpz_gcd(rp, ap, bp);
            if (mpz_cmpabs_ui(rp, 1) != 0)
                mpz_gcd(rp, rp, cp);
        }

        /* The result may be small */
        _fmpz_demote_val(res);
    }
}


void
sample_new(void * arg, ulong count)
{
    fmpz_t r, a, b, c;
    int ix;
    info_t * info = (info_t *) arg;
    slong bits = info->bits;

    fmpz_init(r);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);

    FLINT_TEST_INIT(state);

    prof_start();
    for (ix = 0; ix < 1000 * count; ix++)
    {
        fmpz_randtest(a, state, bits);
        fmpz_randtest(b, state, bits);
        fmpz_randtest(c, state, bits);
        fmpz_gcd3(r, a, b, c);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);
    fmpz_clear(b);

    flint_randclear(state);
}

void
sample_old(void * arg, ulong count)
{
    fmpz_t r, a, b, c;
    int ix;
    info_t * info = (info_t *) arg;
    slong bits = info->bits;

    fmpz_init(r);
    fmpz_init(a);
    fmpz_init(b);
    fmpz_init(c);

    FLINT_TEST_INIT(state);

    prof_start();
    for (ix = 0; ix < 1000 * count; ix++)
    {
        fmpz_randtest(a, state, bits);
        fmpz_randtest(b, state, bits);
        fmpz_randtest(c, state, bits);
        fmpz_gcd3_old(r, a, b, c);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);
    fmpz_clear(b);

    flint_randclear(state);
}

int
main()
{
    double minnew, maxnew, minold, maxold;
    int bits;
    info_t as;

    for (bits = 1; bits <= 200; bits += 3)
    {
        as.bits = bits;

        prof_repeat(&minnew, &maxnew, sample_new, &as);
        prof_repeat(&minold, &maxold, sample_old, &as);
        flint_printf("%3d bits:   (min) %.2fx speedup   (max) %.2fx speedup\n",
                bits, minold / minnew, maxold / maxnew);
    }
}
