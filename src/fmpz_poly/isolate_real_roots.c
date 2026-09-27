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
#include "fmpz_poly/impl.h"

/* isolate the real roots of pol contained in (0,1) */
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

    /* we only consider the open interval (0, 1), so ignore a root at 0 */
    while (len > 0 && fmpz_is_zero(p))
    {
        p++;
        len--;
    }

    if (len <= 1)
    {
        fmpz_clear(c);
        fmpz_clear(one);
        _fmpz_vec_clear(p0, len0);
        return;
    }

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

void
_fmpz_poly_isolate_real_roots_vca(fmpq * exact_roots, slong * n_exact, fmpz * c_array, slong * k_array, slong * n_interval, const fmpz_poly_t pol, int positive_only)
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
        k = _fmpz_poly_scale_positive_roots_0_1(p, len);

        if (k != WORD_MIN)
        {
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
    k = _fmpz_poly_scale_positive_roots_0_1(p, len);
    if (k != WORD_MIN)
    {
        /* the exact positive roots start here (note: roots at zero are
           not included when positive_only is set) */
        slong n_pos_exact_start = *n_exact;

        _fmpz_poly_isolate_real_roots_0_1_vca(exact_roots, n_exact, c_array, k_array, n_interval, p, len);

        if ((c_array != NULL) && (k_array != NULL))
        {
            for (i = n_neg; i < *n_interval; i++)
                k_array[i] += k;
        }

        if (exact_roots != NULL)
        {
            for (i = n_pos_exact_start; i < *n_exact; i++)
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

/* Width, in octaves, of the range spanned by the positive roots of
   (p, len) according to root bounds, or 0 if there are none. */
static slong
_positive_root_spread(const fmpz * p, slong len)
{
    fmpz * t;
    slong i, K, L;

    t = _fmpz_vec_init(len);
    _fmpz_vec_set(t, p, len);
    K = _fmpz_poly_scale_positive_roots_0_1(t, len);
    for (i = 0; i < len; i++)
        fmpz_set(t + i, p + len - 1 - i);
    L = _fmpz_poly_scale_positive_roots_0_1(t, len);
    _fmpz_vec_clear(t, len);

    if (K == WORD_MIN || L == WORD_MIN)
        return 0;

    return K + L;
}

/* Checks Newton's inequalities
   a_k^2 k (n - k) >= a_{k-1} a_{k+1} (k + 1) (n - k + 1),
   which hold if (p, len) has only real roots. */
static int
_newton_inequalities(const fmpz * p, slong len)
{
    slong k, n = len - 1;
    fmpz_t s, t;
    int ok = 1;

    fmpz_init(s);
    fmpz_init(t);

    for (k = 1; k < n && ok; k++)
    {
        fmpz_mul(t, p + k - 1, p + k + 1);
        if (fmpz_sgn(t) <= 0)
            continue;
        fmpz_mul_ui(t, t, k + 1);
        fmpz_mul_ui(t, t, n - k + 1);
        fmpz_mul(s, p + k, p + k);
        fmpz_mul_ui(s, s, k);
        fmpz_mul_ui(s, s, n - k);
        ok = (fmpz_cmp(s, t) >= 0);
    }

    fmpz_clear(s);
    fmpz_clear(t);
    return ok;
}

/* Try the faster methods for root isolation for polynomials where they
   apply. Returns 0 if the caller should use VCA instead. */
static int
_fmpz_poly_isolate_real_roots_try_fast(fmpq * exact_roots, slong * n_exact, fmpz * c_array, slong * k_array, slong * n_interval, const fmpz_poly_t pol, int positive_only)
{
    const fmpz * p = pol->coeffs;
    slong len = pol->length;
    slong j, n_pos, n_neg, bits, spread;
    int s, sp, sn;
    fmpz * q;

    while (len > 0 && fmpz_is_zero(p))
    {
        p++;
        len--;
    }

    if (len < FMPZ_POLY_ISOLATE_REAL_ROOTS_SIGNS_MIN_LEN)
        return 0;

    /* sign variations of p(x) and p(-x) */
    sp = sn = fmpz_sgn(p);
    n_pos = n_neg = 0;
    for (j = 1; j < len; j++)
    {
        s = fmpz_sgn(p + j);
        if (s != 0)
        {
            if (s != sp)
            {
                n_pos++;
                sp = s;
            }
            if (j % 2 == 1)
                s = -s;
            if (s != sn)
            {
                n_neg++;
                sn = s;
            }
        }
    }

    /* If Descartes' rule of signs already shows that there are at most one
       positive and one negative root, VCA terminates immediately. */
    if (n_pos <= 1 && (n_neg <= 1 || positive_only))
        return 0;

    bits = FLINT_ABS(_fmpz_vec_max_bits(p, len));

    /* Sturm sequence, if its coefficients remain small */
    if (len >= FMPZ_POLY_ISOLATE_REAL_ROOTS_STURM_MIN_LEN &&
        _fmpz_poly_isolate_real_roots_sturm(exact_roots, n_exact,
            c_array, k_array, n_interval, pol->coeffs, pol->length, positive_only,
            _fmpz_poly_num_real_roots_sturm_bound(len, bits)))
        return 1;

    /* Sign changes, for real-rooted polynomials whose roots are spread
       over many orders of magnitude. For real-rooted polynomials,
       n_pos + n_neg must equal the degree and Newton's inequalities must
       hold. */
    if (n_pos + n_neg == len - 1 && _newton_inequalities(p, len))
    {
        spread = (n_pos != 0) ? _positive_root_spread(p, len) : 0;

        if (!positive_only && n_neg != 0)
        {
            q = _fmpz_vec_init(len);
            for (j = 0; j < len; j++)
            {
                if (j % 2 == 1)
                    fmpz_neg(q + j, p + j);
                else
                    fmpz_set(q + j, p + j);
            }
            spread += _positive_root_spread(q, len);
            _fmpz_vec_clear(q, len);
        }

        /* we require that the root bounds span at least one octave per root
           on average (for example, the roots of the degree n Eulerian
           polynomial span about 2n octaves) */
        if (spread >= len - 1 &&
            _fmpz_poly_isolate_real_roots_signs(exact_roots, n_exact,
                c_array, k_array, n_interval, pol->coeffs, pol->length,
                positive_only, _fmpz_poly_isolate_real_roots_signs_budget(len)))
            return 1;
    }

    return 0;
}

void fmpz_poly_isolate_real_roots(fmpq * exact_roots, slong * n_exact, fmpz * c_array, slong * k_array, slong * n_interval, const fmpz_poly_t pol)
{
    if (!_fmpz_poly_isolate_real_roots_try_fast(exact_roots, n_exact, c_array, k_array, n_interval, pol, 0))
        _fmpz_poly_isolate_real_roots_vca(exact_roots, n_exact, c_array, k_array, n_interval, pol, 0);
}

void fmpz_poly_isolate_positive_roots(fmpq * exact_roots, slong * n_exact, fmpz * c_array, slong * k_array, slong * n_interval, const fmpz_poly_t pol)
{
    if (!_fmpz_poly_isolate_real_roots_try_fast(exact_roots, n_exact, c_array, k_array, n_interval, pol, 1))
        _fmpz_poly_isolate_real_roots_vca(exact_roots, n_exact, c_array, k_array, n_interval, pol, 1);
}

