/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_factor.h"

static slong trial_cutoff[15] = {4, 4, 4, 6, 11, 18, 31, 54, 97, 172, 309, 564, 1028, 1900, 3512};

/*
   Tuning values to give a roughly 1/3 chance of finding a factor of the given
   number of bits. Parameters are {bits, B1, curves}. B2 should be taken to be
   100*B1. Tuning values are a bit rough from 62 bits on.
*/
static slong ecm_tuning[][3] =
{
    {0, 0, 0}, {2, 1, 1}, {4, 3, 1}, {6, 5, 1}, {8, 7, 1},
    {10, 9, 2}, {12, 11, 2}, {14, 13, 2}, {16, 15, 2}, {18, 17, 2},
    {20, 19, 5}, {22, 25, 8}, {24, 32, 10}, {26, 47, 11}, {28, 67, 13},
    {30, 102, 13}, {32, 126, 13}, {34, 207, 15}, {36, 293, 16}, {38, 415, 17},
    {40, 610, 18}, {42, 920, 18}, {44, 1270, 20}, {46, 1800, 20}, {48, 2650, 20},
    {50, 3850, 21}, {52, 5300, 22}, {54, 8500, 22}, {56, 10000, 26}, {58, 12000, 33},
    {60, 14000, 42}, {62, 15000, 57}, {64, 16500, 72}, {66, 18000, 87}, {68, 22000, 102},
    {70, 26000, 117}, {72, 30000, 131}, {74, 40000, 146}, {76, 50000, 161}, {78, 60000, 175},
    {80, 70000, 190}, {82, 80000, 205}, {84, 100000, 220}, {86, 140000, 240}, {88, 190000, 255},
    {90, 240000, 291}, {92, 280000, 318}, {94, 320000, 345}, {96, 370000, 372}, {98, 420000, 400},
    {100, 470000, 430}
};

int _is_prime(const fmpz_t n, int proved)
{
    if (proved)
        return fmpz_is_prime(n);
    else
    	return fmpz_is_probabprime(n);
}

void remove_found_factors(fmpz_factor_t factor, fmpz_t n, fmpz_t f)
{
    slong i;
    fmpz_factor_t fac;

    fmpz_tdiv_q(n, n, f);

    fmpz_factor_init(fac);
    fmpz_factor_no_trial(fac, f);

    for (i = 0; i < fac->num; i++)
        fac->exp[i] += fmpz_remove(n, n, fac->p + i);

    _fmpz_factor_concat(factor, fac, 1);

    fmpz_factor_clear(fac);
}

int fmpz_factor_smooth(fmpz_factor_t factor, const fmpz_t n,
		                                        slong bits, int proved)
{
    ulong exp;
    mp_limb_t p;
    __mpz_struct * xsrc;
    mp_ptr xd;
    mp_size_t xsize;
    slong found;
    slong trial_stop;
    slong * idx;
    slong i, b, bits2, istride;
    const mp_limb_t * primes;
    int ret = 0;

    TMP_INIT;

    if (!COEFF_IS_MPZ(*n))
    {
        fmpz_factor_si(factor, *n);
        return 1;
    }

    _fmpz_factor_set_length(factor, 0);

    /* Get sign and size */
    xsrc = COEFF_TO_PTR(*n);
    if (xsrc->_mp_size < 0)
    {
        xsize = -(xsrc->_mp_size);
        factor->sign = -1;
    }
    else
    {
        xsize = xsrc->_mp_size;
        factor->sign = 1;
    }

    /* Just a single limb */
    if (xsize == 1)
    {
        _fmpz_factor_extend_factor_ui(factor, xsrc->_mp_d[0]);
        return 1;
    }

    /* Create a temporary copy to be mutated */
    TMP_START;
    xd = TMP_ALLOC(xsize * sizeof(mp_limb_t));
    flint_mpn_copyi(xd, xsrc->_mp_d, xsize);

    /* Factor out powers of two */
    xsize = flint_mpn_remove_2exp(xd, xsize, &exp);
    if (exp != 0)
        _fmpz_factor_append_ui(factor, UWORD(2), exp);

    if (bits <= 0)
    {
        flint_throw(FLINT_ERROR, "(fmpz_factor_smooth) Number of bits must be at least 1\n");
    }

    if (bits <= 15)
       trial_stop = trial_cutoff[bits - 1];
    else
       trial_stop = 3512;

    b = fmpz_sizeinbase(n, 2) - exp;
    idx = (slong *) flint_malloc((5 + b/4)*sizeof(slong));

    found = flint_mpn_factor_trial_tree(idx, xd, xsize, trial_stop);

    if (found)
    {
        primes = n_primes_arr_readonly(trial_stop);

        for (i = 0; i < found; i++)
        {
            p = primes[idx[i]];

            if (p == 2)
                continue;

            exp = 1;
            xsize = flint_mpn_divexact_1(xd, xsize, p);

            /* Check if p^2 divides n */
            if (flint_mpn_divisible_1_odd(xd, xsize, p))
            {
                xsize = flint_mpn_divexact_1(xd, xsize, p);
                exp = 2;
            }

            /* If we're up to cubes, then maybe there are higher powers */
            if (exp == 2 && flint_mpn_divisible_1_odd(xd, xsize, p))
            {
                xsize = flint_mpn_divexact_1(xd, xsize, p);
                xsize = flint_mpn_remove_power_ascending(xd, xsize, &p, 1, &exp);
                exp += 3;
            }

            _fmpz_factor_append_ui(factor, p, exp);
        }
    }

    if (xsize == 1)
    {
        /* Any single-limb factor left? */
        if (xd[0] != 1)
            _fmpz_factor_extend_factor_ui(factor, xd[0]);

        ret = 1;
    }
    else
    {
        fmpz_t n2, f;
        __mpz_struct * data;

        fmpz_init2(n2, xsize);

        data = _fmpz_promote(n2);
        flint_mpn_copyi(data->_mp_d, xd, xsize);
        data->_mp_size = xsize;

        if (proved != -1 && _is_prime(n2, proved))
        {
            _fmpz_factor_append(factor, n2, 1);
            ret = 1;
        }
        else
        {
            fmpz_t root;

            fmpz_init(root);

            exp = fmpz_is_perfect_power(root, n2);

            if (exp != 0)
            {
                fmpz_factor_t fac;

                fmpz_factor_init(fac);

                ret = fmpz_factor_smooth(fac, root, bits, proved);
                fmpz_set_ui(n2, 1);

                _fmpz_factor_concat(factor, fac, exp);

                fmpz_factor_clear(fac);
            }
            else if (bits >= 16) /* trial factored already up to 15 bits */
            {
                int found;
                flint_rand_t state;

                fmpz_init(f);
                flint_randinit(state);

                /* currently only tuning values up to factors of 100 bits */
                bits = FLINT_MIN(bits, 100);
                bits2 = (bits + 1)/2;

                /* tuning is in increments of 2 bits */
                istride = 3;
                /* start with 18-22 bits, advance by 6 bits at a time */
                for (i = 9 + (bits2 % 3); i <= bits2; i += istride)
                {
                    found = fmpz_factor_ecm(f, ecm_tuning[i][2],
                            ecm_tuning[i][1], ecm_tuning[i][1]*100, state, n2);

                    if (found != 0)
                    {
                        /* make sure all prime divisors in factor are removed from n2 */
                        remove_found_factors(factor, n2, f);

                        if (fmpz_is_one(n2))
                        {
                            ret = 1;

                            break;
                        }

                        /* if what remains is below the bound, just factor it */
                        if (fmpz_sizeinbase(n2, 2) < bits)
                        {
                            fmpz_factor_no_trial(factor, n2);

                            ret = 1;

                            break;
                        }

                        if (_is_prime(n2, proved))
                        {
                            _fmpz_factor_append(factor, n2, 1);

                            ret = 1;

                            break;
                        }

                        i -= istride; /* redo with the same parameters if factor found */
                    }
                }

                flint_randclear(state);
                fmpz_clear(f);
            }
        }

        if (ret != 1 && !fmpz_is_one(n2))
           _fmpz_factor_append(factor, n2, 1); /* place cofactor in factor struct */
        else
           ret = 1;

        fmpz_clear(n2);
    }

    flint_free(idx);

    TMP_END;
    return ret;
}

