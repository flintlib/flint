/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_factor.h"
#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"

#define FACTOR_SMOOTH_BOUND 32

/* Write |n| = A * B^2 */
static void
_fmpz_factor_square_root(fmpz_t A, fmpz_t B, const fmpz_t n, slong smooth_bound)
{
    fmpz_factor_t fac;
    fmpz_t t;
    slong i;

    fmpz_factor_init(fac);
    fmpz_factor_smooth(fac, n, smooth_bound, -1); /* -1  =>  no primality test */

    fmpz_one(A);
    fmpz_one(B);
    fmpz_init(t);

    for (i = 0; i < fac->num; i++)
    {
        if (fac->exp[i] == 1)
        {
            fmpz_mul(A, A, fac->p + i);
        }
        else if (fac->exp[i] == 2)
        {
            fmpz_mul(B, B, fac->p + i);
        }
        else
        {
            fmpz_pow_ui(t, fac->p + i, fac->exp[i] / 2);
            fmpz_mul(B, B, t);
            if (fac->exp[i] % 2)
                fmpz_mul(A, A, fac->p + i);
        }
    }

    fmpz_factor_clear(fac);
    fmpz_clear(t);
}

ca_field_ptr ca_field_cache_lookup_qqbar(ca_field_cache_t cache, const qqbar_t x, ca_ctx_t ctx);

ca_field_ptr ca_ctx_get_field_qqbar(ca_ctx_t ctx, const qqbar_t x)
{
    ca_ext_t ext;
    ca_ext_struct * ext_ptr[1];
    ca_field_ptr field;

    field = ca_field_cache_lookup_qqbar(CA_CTX_FIELD_CACHE(ctx), x, ctx);

    if (field == NULL)
    {
        /* todo: shallow copy */
        ca_ext_init_qqbar(ext, x, ctx);
        ext_ptr[0] = ca_ext_cache_insert(CA_CTX_EXT_CACHE(ctx), ext, ctx);
        field = ca_field_cache_insert_ext(CA_CTX_FIELD_CACHE(ctx), ext_ptr, 1, ctx);
        ca_ext_clear(ext, ctx);
    }

    return field;
}

#if 0

/* Find if sqrt(A) is cached; otherwise create new field */
/* Todo: fast search table */
slong ca_ctx_get_quadratic_field(ca_ctx_t ctx, const fmpz_t A)
{
    qqbar_t T;
    slong i;

    for (i = 0; i < ctx->fields_len; i++)
    {
        if (ctx->fields[i].type == CA_FIELD_TYPE_NF)
        {
            const qqbar_struct * v;

            v = CA_FIELD_NF_QQBAR(ctx->fields + i);

            if (qqbar_degree(v) == 2)
            {
                if (fmpz_is_one(QQBAR_COEFFS(v) + 2) &&
                    fmpz_is_zero(QQBAR_COEFFS(v) + 1) &&
                    (fmpz_cmpabs(QQBAR_COEFFS(v), A) == 0) &&
                    !fmpz_equal(QQBAR_COEFFS(v), A))
                {
                    if (fmpz_sgn(A) < 0)
                    {
                        if (qqbar_sgn_im(v) > 0)
                            return i;
                    }
                    else
                    {
                        if (qqbar_sgn_re(v) >= 0)
                            return i;
                    }
                }
            }
        }
    }

    i = ctx->fields_len;

    if (i >= ctx->fields_alloc)
    {
        ctx->fields = (ca_field_struct *) flint_realloc(ctx->fields, sizeof(ca_field_struct) * 2 * ctx->fields_alloc);
        ctx->fields_alloc = 2 * ctx->fields_alloc;
    }

    ctx->fields_len = i + 1;

    qqbar_init(T);
    qqbar_set_fmpz(T, A);
    qqbar_sqrt(T, T);

    ca_field_init_nf(ctx->fields + i, T);

    qqbar_clear(T);

    return i;
}

#else

ca_field_srcptr ca_ctx_get_quadratic_field(ca_ctx_t ctx, const fmpz_t A)
{
    ca_field_srcptr res;
    qqbar_t x;
    qqbar_init(x);

#if 0
    qqbar_set_fmpz(x, A);
    qqbar_sqrt(x, x);
#else
    fmpz_poly_fit_length(QQBAR_POLY(x), 3);
    _fmpz_poly_set_length(QQBAR_POLY(x), 3);
    fmpz_neg(QQBAR_COEFFS(x) + 0, A);
    fmpz_zero(QQBAR_COEFFS(x) + 1);
    fmpz_one(QQBAR_COEFFS(x) + 2);
    acb_set_fmpz(QQBAR_ENCLOSURE(x), A);
    acb_sqrt(QQBAR_ENCLOSURE(x), QQBAR_ENCLOSURE(x), QQBAR_DEFAULT_PREC);
#endif
    res = ca_ctx_get_field_qqbar(ctx, x);
    qqbar_clear(x);
    return res;
}

#endif

ca_field_srcptr ca_ctx_get_cyclotomic_field(ca_ctx_t ctx, ulong n)
{
    ca_field_srcptr res;
    qqbar_t x;
    qqbar_init(x);
    qqbar_root_of_unity(x, 1, n);
    res = ca_ctx_get_field_qqbar(ctx, x);
    qqbar_clear(x);
    return res;
}

static int
fmpz_discr_3(fmpz_t t, const fmpz_t D)
{
    if (fmpz_sgn(D) < 0)
        return 0;

    if (!fmpz_divisible_si(D, 3))
        return 0;

    fmpz_divexact_ui(t, D, 3);

    if (!fmpz_is_square(t))
        return 0;

    return 1;
}

void
ca_set_qqbar(ca_t res, const qqbar_t x, ca_ctx_t ctx)
{
    slong d;

    d = qqbar_degree(x);

    if (d == 1)
    {
        _ca_make_fmpq(res, ctx);
        qqbar_get_fmpq(CA_FMPQ(res), x);
    }
    else if (d == 2)
    {
        const fmpz *a, *b, *c;
        fmpz_t D, t;
        fmpz * res_num;
        fmpz * res_den;

        a = QQBAR_COEFFS(x) + 2;
        b = QQBAR_COEFFS(x) + 1;
        c = QQBAR_COEFFS(x) + 0;

        /* x = (-b +/- sqrt(b^2 - 4ac))/(2a) */

        fmpz_init(D);
        fmpz_init(t);
        fmpz_mul(D, a, c);
        fmpz_mul_2exp(D, D, 2);
        fmpz_submul(D, b, b);

        /* -D is a square <=> element of Q(i) */
        if (fmpz_is_square(D))
        {
            fmpz_sqrt(D, D);

            _ca_make_field_element(res, ctx->field_qq_i, ctx);

            res_num = QNF_ELEM_NUMREF(CA_NF_ELEM(res));
            res_den = QNF_ELEM_DENREF(CA_NF_ELEM(res));

            if (qqbar_sgn_im(x) > 0)
                fmpz_set(res_num + 1, D);
            else
                fmpz_neg(res_num + 1, D);

            fmpz_neg(res_num, b);
            fmpz_mul_2exp(res_den, a, 1);

            /* todo: is gcd avoidable? */
            /* todo: at least use fmpz_gcd3 */
            fmpz_gcd(D, res_num, res_num + 1);
            fmpz_gcd(D, D, res_den);

            if (!fmpz_is_one(D))
            {
                fmpz_divexact(res_num, res_num, D);
                fmpz_divexact(res_num + 1, res_num + 1, D);
                fmpz_divexact(res_den, res_den, D);
            }
        }
        else if (fmpz_discr_3(t, D))  /* Special case for Q(zeta_3) */
        {
            ca_field_srcptr K;

            fmpz_sqrt(D, t);

            K = ca_ctx_get_cyclotomic_field(ctx, 3);
            _ca_make_field_element(res, K, ctx);

            res_num = QNF_ELEM_NUMREF(CA_NF_ELEM(res));
            res_den = QNF_ELEM_DENREF(CA_NF_ELEM(res));

            if (qqbar_sgn_im(x) < 0)
            {
                fmpz_neg(D, D);
            }

            fmpz_sub(res_num, D, b);
            fmpz_mul_2exp(res_num + 1, D, 1);
            fmpz_mul_2exp(res_den, a, 1);

            /* todo: is gcd avoidable? */
            /* todo: at least use fmpz_gcd3 */
            fmpz_gcd(D, res_num, res_num + 1);
            fmpz_gcd(D, D, res_den);

            if (!fmpz_is_one(D))
            {
                fmpz_divexact(res_num, res_num, D);
                fmpz_divexact(res_num + 1, res_num + 1, D);
                fmpz_divexact(res_den, res_den, D);
            }
        }
        else
        {
            fmpz_t A, B;
            ca_field_srcptr field;

            fmpz_neg(D, D);

            fmpz_init(A);
            fmpz_init(B);

            /*
            sqrt(|D|) = A * B^2   =>   sqrt(D) = sqrt(A) * B     (D > 0)
                                                 sqrt(-A) * B    (D < 0)
            */
            _fmpz_factor_square_root(A, B, D, FACTOR_SMOOTH_BOUND);
            if (fmpz_sgn(D) < 0)
                fmpz_neg(A, A);

            field = ca_ctx_get_quadratic_field(ctx, A);
            _ca_make_field_element(res, field, ctx);

            res_num = QNF_ELEM_NUMREF(CA_NF_ELEM(res));
            res_den = QNF_ELEM_DENREF(CA_NF_ELEM(res));

            /* x = (-b +/- B*sqrt(A))/(2a) */

            fmpz_neg(res_num, b);
            fmpz_mul_2exp(res_den, a, 1);

            /* determine the correct sign */
            {
                if (fmpz_sgn(D) < 0)
                {
                    if (qqbar_sgn_im(x) > 0)
                        fmpz_set(res_num + 1, B);
                    else
                        fmpz_neg(res_num + 1, B);
                }
                else if (fmpz_is_zero(b))
                {
                    if (qqbar_sgn_re(x) > 0)
                        fmpz_set(res_num + 1, B);
                    else
                        fmpz_neg(res_num + 1, B);
                }
                else
                {
                    arb_t r1, r2;
                    slong prec;

                    arb_init(r1);
                    arb_init(r2);

                    for (prec = 64; ; prec *= 2)
                    {
                        arb_sqrt_fmpz(r1, A, prec);
                        arb_mul_fmpz(r1, r1, B, prec);
                        arb_add_fmpz(r2, r1, b, prec);
                        arb_neg(r2, r2);
                        arb_sub_fmpz(r1, r1, b, prec);
                        arb_div_fmpz(r1, r1, a, prec);
                        arb_div_fmpz(r2, r2, a, prec);
                        arb_mul_2exp_si(r1, r1, -1);
                        arb_mul_2exp_si(r2, r2, -1);

                        if (arb_overlaps(r1, acb_realref(QQBAR_ENCLOSURE(x))) &&
                            !arb_overlaps(r2, acb_realref(QQBAR_ENCLOSURE(x))))
                        {
                            fmpz_set(res_num + 1, B);
                            break;
                        }

                        if (arb_overlaps(r2, acb_realref(QQBAR_ENCLOSURE(x))) &&
                            !arb_overlaps(r1, acb_realref(QQBAR_ENCLOSURE(x))))
                        {
                            fmpz_neg(res_num + 1, B);
                            break;
                        }
                    }

                    arb_clear(r1);
                    arb_clear(r2);
                }
            }

            /* todo: use fmpz_gcd3 */
            fmpz_gcd(D, res_num, res_num + 1);
            fmpz_gcd(D, D, res_den);

            if (!fmpz_is_one(D))
            {
                fmpz_divexact(res_num, res_num, D);
                fmpz_divexact(res_num + 1, res_num + 1, D);
                fmpz_divexact(res_den, res_den, D);
            }

            fmpz_clear(A);
            fmpz_clear(B);
        }

        fmpz_clear(D);
        fmpz_clear(t);
    }
    else
    {
        ca_field_srcptr field;
        field = ca_ctx_get_field_qqbar(ctx, x);
        _ca_make_field_element(res, field, ctx);
        nf_elem_gen(CA_NF_ELEM(res), CA_FIELD_NF(field));
    }
}
