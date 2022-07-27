#include "ca.h"
#include "ca_mat.h"
#include "ca_poly.h"
#include "fexpr.h"
#include "gr.h"

#define GR_CA_CTX(ring_ctx) ((ca_ctx_struct *)((ring_ctx)->elem_ctx))

int
_gr_ca_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_RR_CA)
        gr_stream_write(out, "Real numbers (ca)");
    else if (ctx->which_ring == GR_CTX_CC_CA)
        gr_stream_write(out, "Complex numbers (ca)");
    else if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA)
        gr_stream_write(out, "Real algebraic numbers (ca)");
    else
        gr_stream_write(out, "Complex algebraic numbers (ca)");
    return GR_SUCCESS;
}

void
gr_ctx_ca_set_option(gr_ctx_t ctx, slong option, slong value)
{
    ca_ctx_set_option(GR_CA_CTX(ctx), option, value);
}

slong
gr_ctx_ca_get_option(gr_ctx_t ctx, slong option)
{
    return ca_ctx_get_option(GR_CA_CTX(ctx), option);
}

void
_gr_ca_init(ca_t x, gr_ctx_t ctx)
{
    ca_init(x, GR_CA_CTX(ctx));
}

void
_gr_ca_clear(ca_t x, gr_ctx_t ctx)
{
    ca_clear(x, GR_CA_CTX(ctx));
}

void
_gr_ca_swap(ca_t x, ca_t y, gr_ctx_t ctx)
{
    ca_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

/* todo: limits */
/* todo: faster real/algebraic constructions */
int
_gr_ca_randtest(ca_t res, flint_rand_t state, gr_ctx_t ctx)
{
    ca_randtest(res, state, 2, 10, GR_CA_CTX(ctx));

    if (ctx->which_ring == GR_CTX_RR_CA)
    {
        if (ca_check_is_real(res, GR_CA_CTX(ctx)) != T_TRUE)
        {
            ca_randtest_rational(res, state, 10, GR_CA_CTX(ctx));
        }
    }
    else if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA)
    {
        if (ca_check_is_real(res, GR_CA_CTX(ctx)) != T_TRUE ||
            ca_check_is_algebraic(res, GR_CA_CTX(ctx)) != T_TRUE)
        {
            ca_randtest_rational(res, state, 10, GR_CA_CTX(ctx));
        }
    }
    else if (ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        if (ca_check_is_algebraic(res, GR_CA_CTX(ctx)) != T_TRUE)
        {
            ca_randtest_rational(res, state, 10, GR_CA_CTX(ctx));
        }
    }

    return GR_SUCCESS;
}

/* todo */
int
_gr_ca_write(gr_stream_t out, const ca_t x, gr_ctx_t ctx)
{
    gr_stream_write_free(out, ca_get_str(x, GR_CA_CTX(ctx)));
    return GR_SUCCESS;
}

int
_gr_ca_zero(ca_t x, gr_ctx_t ctx)
{
    ca_zero(x, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_one(ca_t x, gr_ctx_t ctx)
{
    ca_one(x, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_set_si(ca_t res, slong v, gr_ctx_t ctx)
{
    ca_set_si(res, v, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_set_ui(ca_t res, ulong v, gr_ctx_t ctx)
{
    ca_set_ui(res, v, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_set_fmpz(ca_t res, const fmpz_t v, gr_ctx_t ctx)
{
    ca_set_fmpz(res, v, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_set_fmpq(ca_t res, const fmpq_t v, gr_ctx_t ctx)
{
    ca_set_fmpq(res, v, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_set_other(ca_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    slong target = ctx->which_ring;

    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            ca_set_fmpz(res, x, GR_CA_CTX(ctx));
            return GR_SUCCESS;

        case GR_CTX_FMPQ:
            ca_set_fmpq(res, x, GR_CA_CTX(ctx));
            return GR_SUCCESS;

        case GR_CTX_REAL_ALGEBRAIC_QQBAR:
            ca_set_qqbar(res, x, GR_CA_CTX(ctx));
            return GR_SUCCESS;

        case GR_CTX_COMPLEX_ALGEBRAIC_QQBAR:
            if (target == GR_CTX_CC_CA ||
                target == GR_CTX_COMPLEX_ALGEBRAIC_CA ||
                qqbar_is_real(x))
            {
                ca_set_qqbar(res, x, GR_CA_CTX(ctx));
                return GR_SUCCESS;
            }
            else
            {
                return GR_DOMAIN;
            }

        case GR_CTX_REAL_ALGEBRAIC_CA:
            ca_transfer(res, GR_CA_CTX(ctx), x, GR_CA_CTX(x_ctx));
            return GR_SUCCESS;

        case GR_CTX_COMPLEX_ALGEBRAIC_CA:
            {
                truth_t ok = T_UNKNOWN;

                if (target == GR_CTX_CC_CA || target == GR_CTX_COMPLEX_ALGEBRAIC_CA)
                {
                    ok = T_TRUE;
                }
                else if (target == GR_CTX_RR_CA)
                {
                    ok = ca_check_is_real(x, GR_CA_CTX(x_ctx));
                }
                else if (target == GR_CTX_REAL_ALGEBRAIC_CA)
                {
                    ok = truth_and(ca_check_is_algebraic(x, GR_CA_CTX(x_ctx)), ca_check_is_real(x, GR_CA_CTX(x_ctx)));
                }

                if (ok == T_TRUE)
                {
                    ca_transfer(res, GR_CA_CTX(ctx), x, GR_CA_CTX(x_ctx));
                    return GR_SUCCESS;
                }

                return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
            }

        case GR_CTX_RR_CA:
            {
                truth_t ok = T_UNKNOWN;

                if (target == GR_CTX_RR_CA || target == GR_CTX_CC_CA)
                {
                    ok = T_TRUE;
                }
                else if (target == GR_CTX_REAL_ALGEBRAIC_CA || target == GR_CTX_COMPLEX_ALGEBRAIC_CA)
                {
                    ok = ca_check_is_algebraic(x, GR_CA_CTX(x_ctx));
                }

                if (ok == T_TRUE)
                {
                    ca_transfer(res, GR_CA_CTX(ctx), x, GR_CA_CTX(x_ctx));
                    return GR_SUCCESS;
                }

                return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
            }

        case GR_CTX_CC_CA:
            {
                truth_t ok = T_UNKNOWN;

                if (ctx->which_ring == GR_CTX_CC_CA)
                {
                    ok = T_TRUE;
                }
                else if (ctx->which_ring == GR_CTX_RR_CA)
                {
                    ok = ca_check_is_real(x, GR_CA_CTX(x_ctx));
                }
                else if (ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
                {
                    ok = ca_check_is_algebraic(x, GR_CA_CTX(x_ctx));
                }
                else if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA)
                {
                    ok = truth_and(ca_check_is_algebraic(x, GR_CA_CTX(x_ctx)), ca_check_is_real(x, GR_CA_CTX(x_ctx)));
                }

                if (ok == T_TRUE)
                {
                    ca_transfer(res, GR_CA_CTX(ctx), x, GR_CA_CTX(x_ctx));
                    return GR_SUCCESS;
                }

                return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
            }
    }

    return GR_UNABLE;
}

int 
_gr_ca_get_arb_with_prec(arb_t res, gr_srcptr x, gr_ctx_t x_ctx, slong prec)
{
    int status = GR_UNABLE;
    truth_t ok;

    acb_t t;
    acb_init(t);

    /* todo: when to use accurate_parts? */
    ca_get_acb(t, x, prec, GR_CA_CTX(x_ctx));

    if (x_ctx->which_ring == GR_CTX_RR_CA || x_ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        (arb_is_zero(acb_imagref(t)) && arb_is_finite(acb_realref(t))))
    {
        status = GR_SUCCESS;
    }
    else
    {
        ok = ca_check_is_real(x, GR_CA_CTX(x_ctx));

        if (ok == T_TRUE)
            status = GR_SUCCESS;
        else
            status = (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
    }

    if (status == GR_SUCCESS)
        arb_set_round(res, acb_realref(t), prec);

    acb_clear(t);
    return status;
}

int 
_gr_ca_get_acb_with_prec(acb_t res, gr_srcptr x, gr_ctx_t x_ctx, slong prec)
{
    /* todo: when to use accurate_parts? */
    ca_get_acb(res, x, prec, GR_CA_CTX(x_ctx));
    acb_set_round(res, res, prec);
    return GR_SUCCESS;
}

int
_gr_ca_get_fmpz(fmpz_t res, const ca_t x, gr_ctx_t ctx)
{
    truth_t integer;

    /* todo: avoid duplicate computation */
    integer = ca_check_is_integer(x, GR_CA_CTX(ctx));

    if (integer == T_TRUE)
        return ca_get_fmpz(res, x, GR_CA_CTX(ctx)) ? GR_SUCCESS : GR_UNABLE;
    else if (integer == T_FALSE)
        return GR_DOMAIN;
    else
        return GR_UNABLE;
}

int
_gr_ca_get_si(slong * res, const ca_t x, gr_ctx_t ctx)
{
    fmpz_t n;
    int status;

    fmpz_init(n);
    status = _gr_ca_get_fmpz(n, x, ctx);

    if (status == GR_SUCCESS)
    {
        if (fmpz_fits_si(n))
            *res = fmpz_get_si(n);
        else
            status = GR_DOMAIN;
    }

    fmpz_clear(n);
    return status;
}

int
_gr_ca_get_ui(ulong * res, const ca_t x, gr_ctx_t ctx)
{
    fmpz_t n;
    int status;

    fmpz_init(n);
    status = _gr_ca_get_fmpz(n, x, ctx);

    if (status == GR_SUCCESS)
    {
        if (fmpz_sgn(n) >= 0 && fmpz_cmp_ui(n, UWORD_MAX) <= 0)
            *res = fmpz_get_ui(n);
        else
            status = GR_DOMAIN;
    }

    fmpz_clear(n);
    return status;
}

int 
_gr_ca_get_d(double * res, gr_srcptr x, gr_ctx_t ctx)
{
    arb_t t;
    int status = GR_UNABLE;

    arb_init(t);
    status = _gr_ca_get_arb_with_prec(t, x, ctx, 64);
    if (status == GR_SUCCESS)
        *res = arf_get_d(arb_midref(t), ARF_RND_NEAR);
    arb_clear(t);
    return status;
}

truth_t
_gr_ca_is_zero(const ca_t x, gr_ctx_t ctx)
{
    return ca_check_is_zero(x, GR_CA_CTX(ctx));
}

truth_t
_gr_ca_is_one(const ca_t x, gr_ctx_t ctx)
{
    return ca_check_is_one(x, GR_CA_CTX(ctx));
}

truth_t
_gr_ca_is_neg_one(const ca_t x, gr_ctx_t ctx)
{
    return ca_check_is_neg_one(x, GR_CA_CTX(ctx));
}

truth_t
_gr_ca_equal(const ca_t x, const ca_t y, gr_ctx_t ctx)
{
    return ca_check_equal(x, y, GR_CA_CTX(ctx));
}

int
_gr_ca_set(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_set(res, x, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_neg(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_neg(res, x, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_add(ca_t res, const ca_t x, const ca_t y, gr_ctx_t ctx)
{
    ca_add(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_add_si(ca_t res, const ca_t x, slong y, gr_ctx_t ctx)
{
    ca_add_si(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_add_ui(ca_t res, const ca_t x, ulong y, gr_ctx_t ctx)
{
    ca_add_ui(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_add_fmpz(ca_t res, const ca_t x, const fmpz_t y, gr_ctx_t ctx)
{
    ca_add_fmpz(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_add_fmpq(ca_t res, const ca_t x, const fmpq_t y, gr_ctx_t ctx)
{
    ca_add_fmpq(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_sub(ca_t res, const ca_t x, const ca_t y, gr_ctx_t ctx)
{
    ca_sub(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_sub_si(ca_t res, const ca_t x, slong y, gr_ctx_t ctx)
{
    ca_sub_si(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_sub_ui(ca_t res, const ca_t x, ulong y, gr_ctx_t ctx)
{
    ca_sub_ui(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_sub_fmpz(ca_t res, const ca_t x, const fmpz_t y, gr_ctx_t ctx)
{
    ca_sub_fmpz(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_sub_fmpq(ca_t res, const ca_t x, const fmpq_t y, gr_ctx_t ctx)
{
    ca_sub_fmpq(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mul(ca_t res, const ca_t x, const ca_t y, gr_ctx_t ctx)
{
    ca_mul(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mul_si(ca_t res, const ca_t x, slong y, gr_ctx_t ctx)
{
    ca_mul_si(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mul_ui(ca_t res, const ca_t x, ulong y, gr_ctx_t ctx)
{
    ca_mul_ui(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mul_fmpz(ca_t res, const ca_t x, const fmpz_t y, gr_ctx_t ctx)
{
    ca_mul_fmpz(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mul_fmpq(ca_t res, const ca_t x, const fmpq_t y, gr_ctx_t ctx)
{
    ca_mul_fmpq(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_inv(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_inv(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_ca_div(ca_t res, const ca_t x, const ca_t y, gr_ctx_t ctx)
{
    ca_div(res, x, y, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_ca_div_si(ca_t res, const ca_t x, slong y, gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        ca_div_si(res, x, y, GR_CA_CTX(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_ca_div_ui(ca_t res, const ca_t x, ulong y, gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        ca_div_ui(res, x, y, GR_CA_CTX(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_ca_div_fmpz(ca_t res, const ca_t x, const fmpz_t y, gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        ca_div_fmpz(res, x, y, GR_CA_CTX(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_ca_div_fmpq(ca_t res, const ca_t x, const fmpq_t y, gr_ctx_t ctx)
{
    if (fmpq_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        ca_div_fmpq(res, x, y, GR_CA_CTX(ctx));
        return GR_SUCCESS;
    }
}

truth_t
_gr_ca_is_invertible(const ca_t x, gr_ctx_t ctx)
{
    return truth_not(ca_check_is_zero(x, GR_CA_CTX(ctx)));
}

int
_gr_ca_pow_ui(ca_t res, const ca_t x, ulong exp, gr_ctx_t ctx)
{
    ca_pow_ui(res, x, exp, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_pow_si(ca_t res, const ca_t x, slong exp, gr_ctx_t ctx)
{
    ca_pow_si(res, x, exp, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_ca_pow_fmpz(ca_t res, const ca_t x, const fmpz_t exp, gr_ctx_t ctx)
{
    ca_pow_fmpz(res, x, exp, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_ca_pow_fmpq(ca_t res, const ca_t x, const fmpq_t exp, gr_ctx_t ctx)
{
    ca_pow_fmpq(res, x, exp, GR_CA_CTX(ctx));

    if (ctx->which_ring == GR_CTX_RR_CA || ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA)
    {
        truth_t real;

        real = ca_check_is_real(res, GR_CA_CTX(ctx));

        if (real == T_UNKNOWN)
            return GR_UNABLE;

        if (real == T_FALSE)
            return GR_DOMAIN;
    }

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_ca_pow(ca_t res, const ca_t x, const ca_t exp, gr_ctx_t ctx)
{
    ca_pow(res, x, exp, GR_CA_CTX(ctx));

    if (ctx->which_ring == GR_CTX_RR_CA || ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA)
    {
        truth_t real;

        real = ca_check_is_real(res, GR_CA_CTX(ctx));

        if (real == T_UNKNOWN)
            return GR_UNABLE;

        if (real == T_FALSE)
            return GR_DOMAIN;
    }

    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA || ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t algebraic;

        algebraic = ca_check_is_algebraic(res, GR_CA_CTX(ctx));

        if (algebraic == T_UNKNOWN)
            return GR_UNABLE;

        if (algebraic == T_FALSE)
            return GR_DOMAIN;
    }

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

truth_t
_gr_ca_is_square(const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_RR_CA || ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA)
    {
        return truth_not(ca_check_is_negative_real(x, GR_CA_CTX(ctx)));
    }
    else
    {
        return T_TRUE;
    }
}

int
_gr_ca_sqrt(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_sqrt(res, x, GR_CA_CTX(ctx));

    if (ctx->which_ring == GR_CTX_RR_CA || ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA)
    {
        truth_t real;

        real = ca_check_is_real(res, GR_CA_CTX(ctx));

        if (real == T_UNKNOWN)
            return GR_UNABLE;

        if (real == T_FALSE)
            return GR_DOMAIN;
    }

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_ca_rsqrt(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_sqrt(res, x, GR_CA_CTX(ctx));
    ca_inv(res, res, GR_CA_CTX(ctx));

    if (ctx->which_ring == GR_CTX_RR_CA || ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA)
    {
        truth_t real;

        real = ca_check_is_real(res, GR_CA_CTX(ctx));

        if (real == T_UNKNOWN)
            return GR_UNABLE;

        if (real == T_FALSE)
            return GR_DOMAIN;
    }

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    if (ca_is_special(res, GR_CA_CTX(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int
_gr_ca_floor(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_floor(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
_gr_ca_ceil(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_ceil(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

/* todo: trunc, nint in calcium */

int _arb_trunc(arb_t res, const arb_t x, slong prec);
int _arb_nint(arb_t res, const arb_t x, slong prec);

int
_gr_ca_trunc(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    acb_t t;
    int status;

    acb_init(t);

    ca_get_acb(t, x, 64, GR_CA_CTX(ctx));

    if (arf_cmpabs_2exp_si(arb_midref(acb_realref(t)), -1) < 0 && mag_cmp_2exp_si(arb_radref(acb_realref(t)), -1) < 0)
    {
        ca_zero(res, GR_CA_CTX(ctx));
        status = GR_SUCCESS;
    }
    else if (arb_is_positive(acb_realref(t)))
    {
        status = _gr_ca_floor(res, x, ctx);
    }
    else if (arb_is_negative(acb_realref(t)))
    {
        status = _gr_ca_ceil(res, x, ctx);
    }
    else
    {
        status = GR_UNABLE;
    }

    acb_clear(t);

    return status;
}

/* todo: fast numerical path */
int
_gr_ca_nint(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ca_check_is_integer(x, GR_CA_CTX(ctx)) == T_TRUE)
    {
        ca_set(res, x, GR_CA_CTX(ctx));
        return GR_SUCCESS;
    }
    else
    {
        ca_t t;
        truth_t integer;
        int status = GR_SUCCESS;

        ca_init(t, GR_CA_CTX(ctx));

        ca_set_d(t, 0.5, GR_CA_CTX(ctx));
        ca_add(t, x, t, GR_CA_CTX(ctx));
        ca_re(t, t, GR_CA_CTX(ctx));
        ca_floor(res, t, GR_CA_CTX(ctx));

        integer = ca_check_is_integer(t, GR_CA_CTX(ctx));

        if (integer == T_TRUE)
        {
            fmpz_t m;
            fmpz_init(m);
            if (ca_get_fmpz(m, t, GR_CA_CTX(ctx)))
            {
                if (fmpz_is_odd(m))
                    ca_sub_ui(res, res, 1, GR_CA_CTX(ctx));
            }
            else
            {
                status = GR_UNABLE;
            }
            fmpz_clear(m);
        }
        else if (integer == T_UNKNOWN)
        {
            status = GR_UNABLE;
        }

        ca_clear(t, GR_CA_CTX(ctx));
        return status;
    }
}

int
_gr_ca_i(ca_t res, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_RR_CA)
        return GR_DOMAIN;

    ca_i(res, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
_gr_ca_abs(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_abs(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
_gr_ca_conj(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_conj(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

/* todo: exploit when we know that the field is real */
int
_gr_ca_re(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_re(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

/* todo: exploit when we know that the field is real */
int
_gr_ca_im(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_im(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
_gr_ca_sgn(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_sgn(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
_gr_ca_csgn(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    ca_csgn(res, x, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}


#define CMP_UNDEFINED -2
#define CMP_UNKNOWN -3

int
_ca_cmp(const ca_t x, const ca_t y, ca_ctx_t ctx);

int
_gr_ca_cmp(int * res, const ca_t x, const ca_t y, gr_ctx_t ctx)
{
    int cmp = _ca_cmp(x, y, GR_CA_CTX(ctx));

    if (cmp == CMP_UNKNOWN)
        return GR_UNABLE;

    if (cmp == CMP_UNDEFINED)
        return GR_DOMAIN;

    *res = cmp;
    return GR_SUCCESS;
}

int
_gr_ca_pi(ca_t res, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA ||
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
        return GR_DOMAIN;

    ca_pi(res, GR_CA_CTX(ctx));

    if (ca_is_unknown(res, GR_CA_CTX(ctx)))
        return GR_UNABLE;

    return GR_SUCCESS;
}

int
_gr_ca_exp(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA || 
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_zero(x, GR_CA_CTX(ctx));

        if (ok == T_TRUE)
            return _gr_ca_one(res, ctx);

        return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
    }
    else
    {
        ca_exp(res, x, GR_CA_CTX(ctx));

        if (ca_is_unknown(res, GR_CA_CTX(ctx)))
            return GR_UNABLE;

        return GR_SUCCESS;
    }
}

int
_gr_ca_log(ca_t res, const ca_t x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA || 
        ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA)
    {
        truth_t ok = ca_check_is_one(x, GR_CA_CTX(ctx));

        if (ok == T_TRUE)
            return _gr_ca_zero(res, ctx);

        return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
    }
    else
    {
        ca_log(res, x, GR_CA_CTX(ctx));

        if (ctx->which_ring == GR_CTX_RR_CA)
        {
            truth_t ok = ca_check_is_real(res, GR_CA_CTX(ctx));

            if (ok == T_TRUE)
                return GR_SUCCESS;

            return (ok == T_FALSE) ? GR_DOMAIN : GR_UNABLE;
        }

        if (ca_check_is_infinity(res, GR_CA_CTX(ctx)) == T_TRUE)
            return GR_DOMAIN;

        if (ca_is_unknown(res, GR_CA_CTX(ctx)))
            return GR_UNABLE;

        return GR_SUCCESS;
    }
}

int
_gr_ca_poly_mullow(ca_ptr res,
    ca_srcptr poly1, slong len1,
    ca_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    _ca_poly_mullow(res, poly1, len1, poly2, len2, n, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mat_mul(ca_mat_t res, const ca_mat_t x, const ca_mat_t y, gr_ctx_t ctx)
{
    ca_mat_mul(res, x, y, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_mat_det(ca_t res, const ca_mat_t x, gr_ctx_t ctx)
{
    ca_mat_det(res, x, GR_CA_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ca_ctx_clear(gr_ctx_t ctx)
{
    ca_ctx_clear(GR_CA_CTX(ctx));
    flint_free(ctx->elem_ctx);
    return GR_SUCCESS;
}

truth_t
_gr_ca_ctx_is_algebraically_closed(gr_ctx_t ctx)
{
    return (ctx->which_ring == GR_CTX_CC_CA || ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_CA) ? T_TRUE : T_FALSE;
}

truth_t
_gr_ca_ctx_is_ordered_ring(gr_ctx_t ctx)
{
    return (ctx->which_ring == GR_CTX_RR_CA || ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_CA) ? T_TRUE : T_FALSE;
}

int
_gr_ca_mat_find_nonzero_pivot(slong * pivot_row, ca_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx)
{
    truth_t ok;

    ok = ca_mat_find_pivot(pivot_row, mat, start_row, end_row, column, GR_CA_CTX(ctx));

    if (ok == T_TRUE)
        return GR_SUCCESS;
    else if (ok == T_FALSE)
        return GR_DOMAIN;
    else
        return GR_UNABLE;
}

int _ca_methods_initialized = 0;

gr_static_method_table _ca_methods;

gr_method_tab_input _ca_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_ca_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_ca_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_true},

    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,
                                (gr_funcptr) _gr_ca_ctx_is_algebraically_closed},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) _gr_ca_ctx_is_ordered_ring},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_false},

    {GR_METHOD_INIT,            (gr_funcptr) _gr_ca_init},

    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_ca_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_ca_swap},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_ca_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_ca_write},

    {GR_METHOD_ZERO,            (gr_funcptr) _gr_ca_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_ca_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_ca_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_ca_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_ca_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_ca_equal},

    {GR_METHOD_SET,             (gr_funcptr) _gr_ca_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_ca_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_ca_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_ca_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_ca_set_fmpq},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_ca_set_other},

    {GR_METHOD_GET_SI,          (gr_funcptr) _gr_ca_get_si},
    {GR_METHOD_GET_UI,          (gr_funcptr) _gr_ca_get_ui},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_ca_get_fmpz},
    {GR_METHOD_GET_D,           (gr_funcptr) _gr_ca_get_d},

    {GR_METHOD_NEG,             (gr_funcptr) _gr_ca_neg},

    {GR_METHOD_ADD,             (gr_funcptr) _gr_ca_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_ca_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_ca_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_ca_add_fmpz},
    {GR_METHOD_ADD_FMPQ,        (gr_funcptr) _gr_ca_add_fmpq},

    {GR_METHOD_SUB,             (gr_funcptr) _gr_ca_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_ca_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_ca_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_ca_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,        (gr_funcptr) _gr_ca_sub_fmpq},

    {GR_METHOD_MUL,             (gr_funcptr) _gr_ca_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_ca_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_ca_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_ca_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,        (gr_funcptr) _gr_ca_mul_fmpq},

    {GR_METHOD_DIV,             (gr_funcptr) _gr_ca_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_ca_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_ca_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_ca_div_fmpz},
    {GR_METHOD_DIV_FMPQ,        (gr_funcptr) _gr_ca_div_fmpq},

    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_ca_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_ca_inv},

    {GR_METHOD_POW,             (gr_funcptr) _gr_ca_pow},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_ca_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_ca_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_ca_pow_fmpz},
    {GR_METHOD_POW_FMPQ,        (gr_funcptr) _gr_ca_pow_fmpq},

    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_ca_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_ca_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_ca_rsqrt},

    {GR_METHOD_FLOOR,           (gr_funcptr) _gr_ca_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) _gr_ca_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) _gr_ca_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) _gr_ca_nint},

    {GR_METHOD_I,               (gr_funcptr) _gr_ca_i},
    {GR_METHOD_ABS,             (gr_funcptr) _gr_ca_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) _gr_ca_conj},
    {GR_METHOD_RE,              (gr_funcptr) _gr_ca_re},
    {GR_METHOD_IM,              (gr_funcptr) _gr_ca_im},
    {GR_METHOD_SGN,             (gr_funcptr) _gr_ca_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) _gr_ca_csgn},

    {GR_METHOD_CMP,              (gr_funcptr) _gr_ca_cmp},

    {GR_METHOD_PI,              (gr_funcptr) _gr_ca_pi},
    {GR_METHOD_EXP,             (gr_funcptr) _gr_ca_exp},
    {GR_METHOD_LOG,             (gr_funcptr) _gr_ca_log},

    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_ca_poly_mullow},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_ca_mat_mul},
    {GR_METHOD_MAT_DET,         (gr_funcptr) _gr_ca_mat_det},
    {GR_METHOD_MAT_FIND_NONZERO_PIVOT,      (gr_funcptr) _gr_ca_mat_find_nonzero_pivot},

    {0,                         (gr_funcptr) NULL},
};

void
_gr_ctx_init_ca(gr_ctx_t ctx, int which_ring)
{
    ctx->flags = 0;
    ctx->which_ring = which_ring;
    ctx->sizeof_elem = sizeof(ca_struct);
    ctx->elem_ctx = flint_malloc(sizeof(ca_ctx_struct));
    ctx->size_limit = WORD_MAX;

    ca_ctx_init(GR_CA_CTX(ctx));

    ctx->methods = _ca_methods;

    if (!_ca_methods_initialized)
    {
        gr_method_tab_init(_ca_methods, _ca_methods_input);
        _ca_methods_initialized = 1;
    }
}

void
gr_ctx_init_real_ca(gr_ctx_t ctx)
{
    _gr_ctx_init_ca(ctx, GR_CTX_RR_CA);
}

void
gr_ctx_init_complex_ca(gr_ctx_t ctx)
{
    _gr_ctx_init_ca(ctx, GR_CTX_CC_CA);
}

void
gr_ctx_init_real_algebraic_ca(gr_ctx_t ctx)
{
    _gr_ctx_init_ca(ctx, GR_CTX_REAL_ALGEBRAIC_CA);
}

void
gr_ctx_init_complex_algebraic_ca(gr_ctx_t ctx)
{
    _gr_ctx_init_ca(ctx, GR_CTX_COMPLEX_ALGEBRAIC_CA);
}
