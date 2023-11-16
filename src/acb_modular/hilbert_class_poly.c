/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "thread_support.h"
#include "ulong_extras.h"
#include "fmpz_poly.h"
#include "acb.h"
#include "arb_poly.h"
#include "acb_modular.h"

#ifdef __GNUC__
# define sqrt __builtin_sqrt
#else
# include <math.h>
#endif

static void
bsplit(arb_poly_t pol, const arb_t sqrtD,
            const slong * qbf, slong a, slong b, slong prec)
{
    if (b - a == 0)
    {
        arb_poly_one(pol);
    }
    else if (b - a == 1)
    {
        acb_t z;
        acb_init(z);

        /* j((-b+sqrt(-D))/(2a)) */
        arb_set_si(acb_realref(z), -FLINT_ABS(qbf[3 * a + 1]));
        arb_set(acb_imagref(z), sqrtD);
        acb_div_si(z, z, 2 * qbf[3 * a], prec);
        acb_modular_j(z, z, prec);

        if (qbf[3 * a + 1] < 0)
        {
            /* (x^2 - 2re(j) x + |j|^2) */
            arb_poly_fit_length(pol, 3);
            arb_mul(pol->coeffs, acb_realref(z), acb_realref(z), prec);
            arb_addmul(pol->coeffs, acb_imagref(z), acb_imagref(z), prec);
            arb_mul_2exp_si(pol->coeffs + 1, acb_realref(z), 1);
            arb_neg(pol->coeffs + 1, pol->coeffs + 1);
            arb_one(pol->coeffs + 2);
            _arb_poly_set_length(pol, 3);
        }
        else
        {
            /* (x-j) */
            arb_poly_fit_length(pol, 2);
            arb_neg(pol->coeffs, acb_realref(z));
            arb_one(pol->coeffs + 1);
            _arb_poly_set_length(pol, 2);
        }

        acb_clear(z);
    }
    else
    {
        arb_poly_t tmp;
        arb_poly_init(tmp);
        bsplit(pol, sqrtD, qbf, a, a + (b - a) / 2, prec);
        bsplit(tmp, sqrtD, qbf, a + (b - a) / 2, b, prec);
        arb_poly_mul(pol, pol, tmp, prec);
        arb_poly_clear(tmp);
    }
}

typedef struct
{
    const slong * qbf;
    arb_srcptr sqrtD;
    slong prec;
}
work_t;

static void
basecase(arb_poly_t res, slong a, slong b, work_t * work)
{
    bsplit(res, work->sqrtD, work->qbf, a, b, work->prec);
}

static void
merge(arb_poly_t res, const arb_poly_t a, const arb_poly_t b, work_t * work)
{
    arb_poly_mul(res, a, b, work->prec);
}

int
_acb_modular_hilbert_class_poly(fmpz_poly_t res, slong D,
        const slong * qbf, slong qbf_len, slong prec)
{
    arb_t sqrtD;
    arb_poly_t pol;
    int success;
    work_t work;

    arb_init(sqrtD);
    arb_poly_init(pol);

    arb_set_si(sqrtD, -D);
    arb_sqrt(sqrtD, sqrtD, prec);

    work.qbf = qbf;
    work.sqrtD = sqrtD;
    work.prec = prec;

    flint_parallel_binary_splitting(pol,
        (bsplit_basecase_func_t) basecase,
        (bsplit_merge_func_t) merge,
        sizeof(arb_poly_struct),
        (bsplit_init_func_t) arb_poly_init,
        (bsplit_clear_func_t) arb_poly_clear,
        &work, 0, qbf_len, 1, -1, 0);

/*    bsplit(pol, sqrtD, qbf, 0, qbf_len, prec); */

    success = arb_poly_get_unique_fmpz_poly(res, pol);

    arb_clear(sqrtD);
    arb_poly_clear(pol);

    return success;
}

void
acb_modular_hilbert_class_poly(fmpz_poly_t res, slong D)
{
    slong i, a, b, c, ac, qbf_alloc, qbf_len, prec;
    slong * qbf;
    double lgh;

    if (D >= 0 || ((D & 3) > 1))
    {
        fmpz_poly_zero(res);
        return;
    }

    qbf_alloc = qbf_len = 0;
    qbf = NULL;
    b = D & 1;

    /* Cohen algorithm 5.3.5 */
    do
    {
        ac = (b*b - D) / 4;
        a = FLINT_MAX(b, 1);

        do
        {
            if (ac % a == 0 && n_gcd(n_gcd(a, b), ac/a) == 1)
            {
                c = ac / a;

                if (qbf_len >= qbf_alloc)
                {
                    qbf_alloc = FLINT_MAX(4, FLINT_MAX(qbf_len + 1, qbf_alloc * 2));
                    qbf = flint_realloc(qbf, qbf_alloc * 3 * sizeof(slong));
                }

                if (a == b || a*a == ac || b == 0)
                {
                    qbf[3 * qbf_len + 0] = a;
                    qbf[3 * qbf_len + 1] = b;
                    qbf[3 * qbf_len + 2] = c;
                }
                else
                {
                    /* -b indicates that we have both b and -b */
                    qbf[3 * qbf_len + 0] = a;
                    qbf[3 * qbf_len + 1] = -b;
                    qbf[3 * qbf_len + 2] = c;
                }

                qbf_len++;
            }

            a++;
        }
        while (a*a <= ac);

        b += 2;
    }
    while (3*b*b <= -D);

    /* Estimate precision - see p.7 in http://hal.inria.fr/inria-00001040 */
    lgh = 0.0;
    for (i = 0; i < qbf_len; i++)
    {
        if (qbf[3 * i + 1] < 0)
            lgh += 2.0 / qbf[3 * i];
        else
            lgh += 1.0 / qbf[3 * i];
    }

    lgh = 3.141593 * sqrt(-D) * lgh;
#if 0
    lgh += 3.012 * h;
    prec = lgh * 1.442696;
    prec = prec + 10;
#else
    prec = lgh * 1.442696;     /* heuristic, but more accurate */
    prec = prec * 1.005 + 20;
#endif

    while (!_acb_modular_hilbert_class_poly(res, D, qbf, qbf_len, prec))
    {
        flint_printf("hilbert_class_poly failed at %wd bits of precision\n", prec);
        prec = prec * 1.2 + 10;
    }

    flint_free(qbf);
}

