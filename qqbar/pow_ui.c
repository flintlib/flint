/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/fmpz_poly_factor.h"
#include "arb_fmpz_poly.h"
#include "qqbar.h"

/* todo: use number field arithmetic, other optimisations ... */
void
qqbar_pow_ui(qqbar_t res, const qqbar_t x, ulong n)
{
    if (n == 0)
    {
        qqbar_one(res);
    }
    else if (n == 1)
    {
        qqbar_set(res, x);
    }
    else if (qqbar_is_rational(x))
    {
        fmpq_t t;
        fmpq_init(t);
        fmpz_neg(fmpq_numref(t), QQBAR_COEFFS(x));
        fmpz_set(fmpq_denref(t), QQBAR_COEFFS(x) + 1);
        fmpz_pow_ui(fmpq_numref(t), fmpq_numref(t), n);
        fmpz_pow_ui(fmpq_denref(t), fmpq_denref(t), n);
        qqbar_set_fmpq(res, t);
        fmpq_clear(t);
    }
    else
    {
        /* Fast detection of perfect powers */
        ulong f;
        f = arb_fmpz_poly_deflation(QQBAR_POLY(x));

        if (f % n == 0)
        {
            acb_t z, t, w;
            fmpz_poly_t H;
            slong prec;
            int pure_real, pure_imag;

            fmpz_poly_init(H);
            acb_init(z);
            acb_init(t);
            acb_init(w);

            arb_fmpz_poly_deflate(H, QQBAR_POLY(x), n);

            acb_set(z, QQBAR_ENCLOSURE(x));
            pure_real = (qqbar_sgn_im(x) == 0);
            pure_imag = (qqbar_sgn_re(x) == 0);

            for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
            {
                _qqbar_enclosure_raw(z, QQBAR_POLY(x), z, prec);
                if (pure_real)
                    arb_zero(acb_imagref(z));
                if (pure_imag)
                    arb_zero(acb_realref(z));

                acb_pow_ui(w, z, n, prec);

                if (_qqbar_validate_uniqueness(t, H, w, 2 * prec))
                {
                    fmpz_poly_set(QQBAR_POLY(res), H);
                    acb_set(QQBAR_ENCLOSURE(res), t);
                    break;
                }
            }

            fmpz_poly_clear(H);
            acb_clear(z);
            acb_clear(t);
            acb_clear(w);
            return;
        }

        if (n == 2)
        {
            qqbar_binary_op(res, x, x, 2);
            return;
        }

        /* todo: is left-to-right or right-to-left better? */
        if (n % 2 == 0)
        {
#if 0
            qqbar_pow_ui(res, x, n / 2);
            qqbar_binary_op(res, res, res, 2);
#else
            qqbar_binary_op(res, x, x, 2);
            qqbar_pow_ui(res, res, n / 2);
#endif
        }
        else
        {
            if (res == x)
            {
                qqbar_t t;
                qqbar_init(t);
                qqbar_pow_ui(t, x, n - 1);
                qqbar_binary_op(res, t, x, 2);
                qqbar_clear(t);
            }
            else
            {
                qqbar_pow_ui(res, x, n - 1);
                qqbar_mul(res, res, x);
            }
        }
    }
}

