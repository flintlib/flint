/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/* isolate the real roots of pol contained in [0,1] */
/* using VCA (Vincent-Collins-Akritas) method       */
/* the output are arrays of fmpz c and integers k so that the roots belong to */
/* [c*2^k, (c+1)*2^k[ */
/* if exact_roots is NULL only n_exact is updated */
/* similarly if c_array/k_array is NULL only n_intervals is updated (useful to */
/* count the roots) */
void _fmpz_poly_isolate_real_roots_0_1_vca(fmpq * exact_roots, slong * n_exact,
        fmpz * c_array, slong * k_array, slong * n_intervals,
        const fmpz * pol, slong len)
{
    fmpz_t c;
    slong k;
    fmpz * p;
    fmpz * p0;
    slong i;
    fmpz_t one;
    slong len0 = len;

    fmpz_init(one);
    fmpz_one(one);
    p0 = p = _fmpz_vec_init(len);
    _fmpz_vec_set(p, pol, len);

    fmpz_init(c);
    fmpz_zero(c);
    k = 0;

    while (1)
    {
        /* check for exact zero */
        while (fmpz_is_zero(p) && len)
        {
            if (exact_roots != NULL)
            {
                fmpz_set(fmpq_numref(exact_roots + *n_exact), c);
                fmpz_one(fmpq_denref(exact_roots + *n_exact));

                FLINT_ASSERT(k >= 0);
                fmpq_div_2exp(exact_roots + *n_exact, exact_roots + *n_exact, (ulong)k);
            }
            (*n_exact)++;
            p++;
            len--;
        }

        /* use Descartes bound */
        {
            const slong bound = _fmpz_poly_descartes_bound_0_1(p, len, 2);
            switch(bound)
            {
                case 2:
                case WORD_MAX:
                    /* unknown: go down */
                    k += 1;
                    fmpz_mul_2exp(c, c, 1);
                    _fmpz_poly_scale_2exp(p, len, -1);
                    continue;
                case 1:
                    /* got a root! */
                    if ((c_array != NULL) && (k_array != NULL))
                    {
                        fmpz_set(c_array + *n_intervals, c);
                        k_array[*n_intervals] = -k;
                    }
                    (*n_intervals)++;
                    break;
                case 0:
                    break;
                default:
                    flint_throw(FLINT_ERROR, "ERROR: expected 0,1,WORD_MAX as output from descartes_bound but got %wd\n", bound);
            }

            /* no root: go up */
            fmpz_add_ui(c, c, 1);
            i = (slong)fmpz_val2(c);
            if (k == i)
            {
                fmpz_clear(c);
                fmpz_clear(one);
                _fmpz_vec_clear(p0, len0);
                return;
            }

            /* go to the next node */
            _fmpz_poly_taylor_shift(p, one, len);
            if (i)
            {
                _fmpz_poly_scale_2exp(p, len, i);
                fmpz_fdiv_q_2exp(c, c, (ulong)i);

                FLINT_ASSERT(k >= i);
                k -= i;
            }
        }
    }
}

static void
_fmpz_poly_isolate_real_roots(fmpq * exact_roots, slong * n_exact, fmpz * c_array, slong * k_array, slong * n_interval, const fmpz_poly_t pol, int positive_only)
{
    slong i, k, n_neg, tmp, len, n_zeros, n_neg_exact;
    fmpz * p;

    n_neg = n_zeros = n_neg_exact = *n_exact = *n_interval = 0;
    len = pol->length;

    if (fmpz_poly_is_zero(pol))
        flint_throw(FLINT_ERROR, "ERROR (fmpz_poly_isolate_real_roots): zero polynomial\n");

    /* compute the number zero roots */
    /* (they will be inserted after the negative ones) */
    for (n_zeros = 0; (n_zeros < len) && fmpz_is_zero(pol->coeffs + n_zeros); n_zeros++);
    len -= n_zeros;
    p = _fmpz_vec_init(len);
    _fmpz_vec_set(p, pol->coeffs + n_zeros, len);

    if (!positive_only)
    {
        /* negative roots (use P(-x)) */
        for (i = 1; i < len; i += 2) fmpz_neg(p + i, p + i);
        k = _fmpz_poly_positive_root_upper_bound_2exp(p, len);

        if (k != WORD_MIN)
        {
            _fmpz_poly_scale_2exp(p, len, k);
            _fmpz_poly_isolate_real_roots_0_1_vca(exact_roots, n_exact, c_array, k_array, n_interval, p, len);
            n_neg = *n_interval;
            n_neg_exact = *n_exact;
            if ((c_array != NULL) && (k_array != NULL))
            {
                for (i = 0; i < *n_interval; i++)
                {
                    fmpz_add_ui(c_array + i, c_array + i, 1);
                    fmpz_neg(c_array + i, c_array + i);
                    k_array[i] += k;
                }
                for (i = 0; i < *n_interval / 2; i++)
                {
                    fmpz_swap(c_array + i, c_array + *n_interval - i - 1);
                    tmp = k_array[i];
                    k_array[i] = k_array[*n_interval - i - 1];
                    k_array[*n_interval - i - 1] = tmp;
                }
            }

            if (exact_roots != NULL)
            {
                for (i = 0; i < n_neg_exact; i++)
                {
                    fmpq_neg(exact_roots + i, exact_roots + i);
                    if (k > 0)
                        fmpq_mul_2exp(exact_roots + i, exact_roots + i, (ulong)k);
                    else if (k < 0)
                        fmpq_div_2exp(exact_roots + i, exact_roots + i, (ulong)-k);
                }
                for (i = 0; i < n_neg_exact/2; i++)
                {
                    fmpq_swap(exact_roots + i, exact_roots + *n_exact - i - 1);
                }
            }
        }
        else
            n_neg = 0;

        /* insert zero roots */
        if (exact_roots != NULL)
        {
            for (i = *n_exact; i < *n_exact+n_zeros; i++) fmpq_zero(exact_roots + i);
        }
        *n_exact += n_zeros;
    }

    /* positive roots */
    _fmpz_vec_set(p, pol->coeffs + n_zeros, len);
    k = _fmpz_poly_positive_root_upper_bound_2exp(p, len);
    if (k != WORD_MIN)
    {
        _fmpz_poly_scale_2exp(p, len, k);
        _fmpz_poly_isolate_real_roots_0_1_vca(exact_roots, n_exact, c_array, k_array, n_interval, p, len);

        if ((c_array != NULL) && (k_array != NULL))
        {
            for (i = n_neg; i < *n_interval; i++)
                k_array[i] += k;
        }

        if (exact_roots != NULL)
        {
            for (i = n_neg_exact + n_zeros; i < *n_exact; i++)
            {
                if (k > 0)
                    fmpq_mul_2exp(exact_roots + i, exact_roots + i, (ulong)k);
                else if (k < 0)
                    fmpq_div_2exp(exact_roots + i, exact_roots + i, (ulong)-k);
            }
        }
    }

    _fmpz_vec_clear(p, len);
}

void fmpz_poly_isolate_real_roots(fmpq * exact_roots, slong * n_exact, fmpz * c_array, slong * k_array, slong * n_interval, const fmpz_poly_t pol)
{
    _fmpz_poly_isolate_real_roots(exact_roots, n_exact, c_array, k_array, n_interval, pol, 0);
}

void fmpz_poly_isolate_positive_roots(fmpq * exact_roots, slong * n_exact, fmpz * c_array, slong * k_array, slong * n_interval, const fmpz_poly_t pol)
{
    _fmpz_poly_isolate_real_roots(exact_roots, n_exact, c_array, k_array, n_interval, pol, 1);
}

