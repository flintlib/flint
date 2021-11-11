/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint/fmpz.h"
#include "flint/ulong_extras.h"
#include "flint/profiler.h"

#define EXPBOUND64 34
#define EXPBOUND128 63

static const ulong mulfunc_bound_intersection[] = {
    0,                                          /* mul0            bound */
    0,4611686018427387903,                      /* mul1  and mul2  bound */
    2097151,55108,                              /* mul3  and mul4  bound */
    6208,1448,                                  /* mul5  and mul6  bound */
    511,234,                                    /* mul7  and mul8  bound */
    127,78,                                     /* mul9  and mul10 bound */
    52,38,                                      /* mul11 and mul12 bound */
    28,22,                                      /* mul13 and mul14 bound */
    18,15,                                      /* mul15 and mul16 bound */
    13,11,                                      /* mul17 and mul18 bound */
    9,8,                                        /* mul19 and mul20 bound */
    7,7,                                        /* mul21 and mul22 bound */
    6,6,                                        /* mul23 and mul24 bound */
    5,5,                                        /* mul25 and mul26 bound */
    5,4,                                        /* mul27 and mul28 bound */
    4,4,                                        /* mul29 and mul30 bound */
    4,3,                                        /* mul31 and mul32 bound */
    3,3                                         /* mul33 and mul34 bound */
};

static const ulong mulfunc_bound_newbound[] = {
    0xffffffffffffffff,                         /* mul0            bound */
    0xffffffffffffffff,4611686018427387903,     /* mul1  and mul2  bound */
    3037000499,3037000499,                      /* mul3  and mul4  bound */
    2097151,2097151,                            /* mul5  and mul6  bound */
    55108,55108,                                /* mul7  and mul8  bound */
    6208,6208,                                  /* mul9  and mul10 bound */
    1448,1448,                                  /* mul11 and mul12 bound */
    511,511,                                    /* mul13 and mul14 bound */
    234,234,                                    /* mul15 and mul16 bound */
    127,127,                                    /* mul17 and mul18 bound */
    78,78,                                      /* mul19 and mul20 bound */
    52,52,                                      /* mul21 and mul22 bound */
    38,38,                                      /* mul23 and mul24 bound */
    28,28,                                      /* mul25 and mul26 bound */
    22,22,                                      /* mul27 and mul28 bound */
    18,18,                                      /* mul29 and mul30 bound */
    15,15,                                      /* mul31 and mul32 bound */
    13,13,                                      /* mul33 and mul34 bound */
    11,11,                                      /* mul35 and mul36 bound */
    9,9,                                        /* mul37 and mul38 bound */
    8,8,                                        /* mul39 and mul40 bound */
    7,7,                                        /* mul41 and mul42 bound */
    7,7,                                        /* mul43 and mul44 bound */
    6,6,                                        /* mul45 and mul46 bound */
    6,6,                                        /* mul45 and mul46 bound */
    6,6,                                        /* mul47 and mul48 bound */
    5,5,                                        /* mul49 and mul50 bound */
    5,5,                                        /* mul49 and mul50 bound */
    5,5,                                        /* mul51 and mul52 bound */
    5,5,                                        /* mul53 and mul54 bound */
    4,4,                                        /* mul55 and mul56 bound */
    4,4,                                        /* mul57 and mul58 bound */
    4,4,                                        /* mul59 and mul60 bound */
    4,4,                                        /* mul61 and mul62 bound */
    3                                           /* mul63           bound */
};

void
fmpz_pow_ui_old(fmpz_t f, const fmpz_t g, ulong exp)
{
    fmpz c1;

    if (exp == WORD(0))
    {
        fmpz_one(f);
        return;
    }

    c1 = *g;

    if (!COEFF_IS_MPZ(c1))      /* g is small */
    {
        ulong u1 = FLINT_ABS(c1);
        ulong bits = FLINT_BIT_COUNT(u1);
        if (u1 <= UWORD(1))
        {
            fmpz_set_ui(f, u1);
        }
        else if (exp * bits <= FLINT_BITS - 2)
        {
            fmpz_set_ui(f, n_pow(u1, exp));
        }
        else
        {
            __mpz_struct *mpz_ptr = _fmpz_promote_val(f);

            flint_mpz_set_ui(mpz_ptr, u1);
            flint_mpz_pow_ui(mpz_ptr, mpz_ptr, exp);
            _fmpz_demote_val(f);    /* may actually fit into a small after all */
        }

        if ((c1 < WORD(0)) && (exp & 1)) /* sign is -ve if exp odd and g -ve */
            fmpz_neg(f, f);
    }
    else
    {
        __mpz_struct *mpz_ptr = _fmpz_promote_val(f);
        flint_mpz_pow_ui(mpz_ptr, COEFF_TO_PTR(c1), exp);
        /* no need to demote as it can't get smaller */
    }
}

typedef struct
{
   ulong exp;
}
info_t;

void
sample_new_trivial(void * arg, ulong count)
{
    fmpz_t r, a;
    int ix;
    info_t * info = (info_t *) arg;
    ulong exp = info->exp;

    fmpz_init(r);
    fmpz_init(a);

    FLINT_TEST_INIT(state);

    prof_start();
    for (ix = 0; ix < 10000 * count; ix++)
    {
        fmpz_set_ui(a, n_urandint(state, COEFF_MAX));
        fmpz_pow_ui(r, a, exp);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);

    flint_randclear(state);
}

void
sample_old_trivial(void * arg, ulong count)
{
    fmpz_t r, a;
    int ix;
    info_t * info = (info_t *) arg;
    ulong exp = info->exp;

    fmpz_init(r);
    fmpz_init(a);

    FLINT_TEST_INIT(state);

    prof_start();
    for (ix = 0; ix < 10000 * count; ix++)
    {
        fmpz_set_ui(a, n_urandint(state, COEFF_MAX));
        fmpz_pow_ui_old(r, a, exp);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);

    flint_randclear(state);
}

void
sample_new_64(void * arg, ulong count)
{
    fmpz_t r, a;
    int ix;
    info_t * info = (info_t *) arg;
    ulong exp = info->exp;

    fmpz_init(r);
    fmpz_init(a);

    FLINT_TEST_INIT(state);

    prof_start();
    for (ix = 0; ix < 10000 * count; ix++)
    {
        fmpz_set_ui(a, n_urandint(state, mulfunc_bound_intersection[exp]));
        fmpz_pow_ui(r, a, exp);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);

    flint_randclear(state);
}

void
sample_old_64(void * arg, ulong count)
{
    fmpz_t r, a;
    int ix;
    info_t * info = (info_t *) arg;
    ulong exp = info->exp;

    fmpz_init(r);
    fmpz_init(a);

    FLINT_TEST_INIT(state);

    prof_start();
    for (ix = 0; ix < 10000 * count; ix++)
    {
        fmpz_set_ui(a, n_urandint(state, mulfunc_bound_intersection[exp]));
        fmpz_pow_ui_old(r, a, exp);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);

    flint_randclear(state);
}

void
sample_new_128(void * arg, ulong count)
{
    fmpz_t r, a;
    int ix;
    info_t * info = (info_t *) arg;
    ulong exp = info->exp;

    fmpz_init(r);
    fmpz_init(a);

    FLINT_TEST_INIT(state);

    prof_start();
    for (ix = 0; ix < 10000 * count; ix++)
    {
        fmpz_set_ui(a, n_urandint(state, mulfunc_bound_newbound[exp]));
        fmpz_pow_ui(r, a, exp);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);

    flint_randclear(state);
}

void
sample_old_128(void * arg, ulong count)
{
    fmpz_t r, a;
    int ix;
    info_t * info = (info_t *) arg;
    ulong exp = info->exp;

    fmpz_init(r);
    fmpz_init(a);

    FLINT_TEST_INIT(state);

    prof_start();
    for (ix = 0; ix < 10000 * count; ix++)
    {
        fmpz_set_ui(a, n_urandint(state, mulfunc_bound_newbound[exp]));
        fmpz_pow_ui_old(r, a, exp);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);

    flint_randclear(state);
}

void
sample_new_large(void * arg, ulong count)
{
    fmpz_t r, a;
    int ix;
    info_t * info = (info_t *) arg;
    ulong exp = info->exp;

    fmpz_init(r);
    fmpz_init(a);

    FLINT_TEST_INIT(state);

    prof_start();
    for (ix = 0; ix < 10000 * count; ix++)
    {
        fmpz_randbits(a, state, 100);
        fmpz_pow_ui(r, a, exp);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);

    flint_randclear(state);
}

void
sample_old_large(void * arg, ulong count)
{
    fmpz_t r, a;
    int ix;
    info_t * info = (info_t *) arg;
    ulong exp = info->exp;

    fmpz_init(r);
    fmpz_init(a);

    FLINT_TEST_INIT(state);

    prof_start();
    for (ix = 0; ix < 10000 * count; ix++)
    {
        fmpz_randbits(a, state, 100);
        fmpz_pow_ui_old(r, a, exp);
    }
    prof_stop();

    fmpz_clear(r);
    fmpz_clear(a);

    flint_randclear(state);
}

int
main()
{
    double minnew, maxnew, minold, maxold;
    ulong exp;
    info_t as;

    flint_printf("For smaller base, trivial exponents:\n");
    for (exp = 0; exp <= 1; exp++)
    {
        as.exp = exp;

        prof_repeat(&minnew, &maxnew, sample_new_trivial, &as);
        prof_repeat(&minold, &maxold, sample_old_trivial, &as);
        flint_printf("exp = %lu:   (min) %.2fx speedup   (max) %.2fx speedup\n",
                exp, minold / minnew, maxold / maxnew);
    }

    flint_printf("\nComparing where the old small-method and the new small-method intersect:\n");
    for (exp = 2; exp <= EXPBOUND64; exp++)
    {
        as.exp = exp;

        prof_repeat(&minnew, &maxnew, sample_new_64, &as);
        prof_repeat(&minold, &maxold, sample_old_64, &as);
        flint_printf("exp = %lu:   (min) %.2fx speedup   (max) %.2fx speedup\n",
                exp, minold / minnew, maxold / maxnew);
    }

    flint_printf("\nWhen the new small-method is used:\n");
    for (exp = 2; exp <= EXPBOUND128; exp++)
    {
        as.exp = exp;

        prof_repeat(&minnew, &maxnew, sample_new_128, &as);
        prof_repeat(&minold, &maxold, sample_old_128, &as);
        flint_printf("exp = %lu:   (min) %.2fx speedup   (max) %.2fx speedup\n",
                exp, minold / minnew, maxold / maxnew);
    }

    flint_printf("\nFor larger base:\n");
    for (exp = 0; exp < 100; exp += 10)
    {
        as.exp = exp;

        prof_repeat(&minnew, &maxnew, sample_new_large, &as);
        prof_repeat(&minold, &maxold, sample_old_large, &as);
        flint_printf("exp = %lu:   (min) %.2fx speedup   (max) %.2fx speedup\n",
                exp, minold / minnew, maxold / maxnew);
    }
}
