#include "acb.h"
#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_hypgeom.h"
#include "acf.h"
#include "qqbar.h"
#include "gr.h"

typedef struct
{
    slong prec;
}
gr_acb_ctx;

#define ACB_CTX_PREC(ring_ctx) (((gr_acb_ctx *)((ring_ctx)))->prec)

int
_gr_acb_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Complex numbers (acb, prec = ");
    gr_stream_write_si(out, ACB_CTX_PREC(ctx));
    gr_stream_write(out, ")");
    return GR_SUCCESS;
}

void
_gr_acb_init(acb_t x, const gr_ctx_t ctx)
{
    acb_init(x);
}

void
_gr_acb_clear(acb_t x, const gr_ctx_t ctx)
{
    acb_clear(x);
}

void
_gr_acb_swap(acb_t x, acb_t y, const gr_ctx_t ctx)
{
    acb_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

/* todo: limits */
int
_gr_acb_randtest(acb_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    acb_randtest(res, state, ACB_CTX_PREC(ctx), 10);
    return GR_SUCCESS;
}

/* todo */
int
_gr_acb_write(gr_stream_t out, const acb_t x, const gr_ctx_t ctx)
{
    slong digits = ACB_CTX_PREC(ctx) * 0.30102999566398 + 1;
    int flags = 0 /* ARB_STR_NO_RADIUS */;

    if (arb_is_zero(acb_imagref(x)))
    {
        gr_stream_write_free(out, arb_get_str(acb_realref(x), digits, flags));
    }
    else if (arb_is_zero(acb_realref(x)))
    {
        gr_stream_write_free(out, arb_get_str(acb_imagref(x), digits, flags));
        gr_stream_write(out, "*I");
    }
    else
    {
        gr_stream_write(out, "(");
        gr_stream_write_free(out, arb_get_str(acb_realref(x), digits, flags));

        if ((arb_is_exact(acb_imagref(x)) || (flags & ARB_STR_NO_RADIUS))
                && arf_sgn(arb_midref(acb_imagref(x))) < 0)
        {
            arb_t t;
            arb_init(t);
            arb_neg(t, acb_imagref(x));
            gr_stream_write(out, " - ");
            gr_stream_write_free(out, arb_get_str(t, digits, flags));
            arb_clear(t);
        }
        else
        {
            gr_stream_write(out, " + ");
            gr_stream_write_free(out, arb_get_str(acb_imagref(x), digits, flags));
        }

        gr_stream_write(out, "*I)");
    }

    return GR_SUCCESS;
}

int
_gr_acb_zero(acb_t x, const gr_ctx_t ctx)
{
    acb_zero(x);
    return GR_SUCCESS;
}

int
_gr_acb_one(acb_t x, const gr_ctx_t ctx)
{
    acb_one(x);
    return GR_SUCCESS;
}

int
_gr_acb_set_si(acb_t res, slong v, const gr_ctx_t ctx)
{
    acb_set_si(res, v);
    acb_set_round(res, res, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_set_ui(acb_t res, ulong v, const gr_ctx_t ctx)
{
    acb_set_ui(res, v);
    acb_set_round(res, res, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_set_fmpz(acb_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    acb_set_round_fmpz(res, v, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_set_fmpq(acb_t res, const fmpq_t v, const gr_ctx_t ctx)
{
    acb_set_fmpq(res, v, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_set_d(acb_t res, double x, const gr_ctx_t ctx)
{
    acb_set_d(res, x);
    arb_set_round(acb_realref(res), acb_realref(res), ACB_CTX_PREC(ctx));

    if (!arb_is_finite(acb_realref(res)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int 
_gr_ca_get_acb_with_prec(acb_t res, gr_srcptr x, gr_ctx_t x_ctx, slong prec);

int
_gr_acb_set_other(acb_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            return _gr_acb_set_fmpz(res, x, ctx);

        case GR_CTX_FMPQ:
            return _gr_acb_set_fmpq(res, x, ctx);

        case GR_CTX_REAL_ALGEBRAIC_QQBAR:
        case GR_CTX_COMPLEX_ALGEBRAIC_QQBAR:
            qqbar_get_acb(res, x, ACB_CTX_PREC(ctx));
            return GR_SUCCESS;

        case GR_CTX_RR_CA:
        case GR_CTX_REAL_ALGEBRAIC_CA:
        case GR_CTX_CC_CA:
        case GR_CTX_COMPLEX_ALGEBRAIC_CA:
            return _gr_ca_get_acb_with_prec(res, x, x_ctx, ACB_CTX_PREC(ctx));

        case GR_CTX_REAL_FLOAT_ARF:
            if (arf_is_finite(x))
            {
                arb_set_arf(acb_realref(res), x);
                arb_set_round(acb_realref(res), acb_realref(res), ACB_CTX_PREC(ctx));
                arb_zero(acb_imagref(res));
                return GR_SUCCESS;
            }
            else
            {
                return GR_DOMAIN;
            }

        case GR_CTX_COMPLEX_FLOAT_ACF:
            if (arf_is_finite(acf_realref((acf_srcptr) x)) && arf_is_finite(acf_imagref((acf_srcptr) x)))
            {
                arb_set_arf(acb_realref(res), acf_realref((acf_srcptr) x));
                arb_set_arf(acb_imagref(res), acf_imagref((acf_srcptr) x));
                acb_set_round(res, res, ACB_CTX_PREC(ctx));
                return GR_SUCCESS;
            }
            else
            {
                return GR_DOMAIN;
            }

        case GR_CTX_RR_ARB:
            arb_set_round(acb_realref(res), x, ACB_CTX_PREC(ctx));
            arb_zero(acb_imagref(res));
            return GR_SUCCESS;

        case GR_CTX_CC_ACB:
            acb_set_round(res, x, ACB_CTX_PREC(ctx));
            return GR_SUCCESS;
    }

    return GR_UNABLE;
}

/* xxx: assumes that ctx are not read */
int _gr_arf_get_fmpz(fmpz_t res, const arf_t x, const gr_ctx_t ctx);
int _gr_arf_get_si(slong * res, const arf_t x, const gr_ctx_t ctx);
int _gr_arf_get_ui(ulong * res, const arf_t x, const gr_ctx_t ctx);

int
_gr_acb_get_fmpz(fmpz_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (!acb_is_int(x))
    {
        if (acb_contains_int(x))
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }

    return _gr_arf_get_fmpz(res, arb_midref(acb_realref(x)), NULL);
}

int
_gr_acb_get_si(slong * res, const acb_t x, const gr_ctx_t ctx)
{
    if (!acb_is_int(x))
    {
        if (acb_contains_int(x))
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }

    return _gr_arf_get_si(res, arb_midref(acb_realref(x)), NULL);
}

int
_gr_acb_get_ui(ulong * res, const acb_t x, const gr_ctx_t ctx)
{
    if (!acb_is_int(x))
    {
        if (acb_contains_int(x))
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }

    return _gr_arf_get_ui(res, arb_midref(acb_realref(x)), NULL);
}

int
_gr_acb_get_d(double * res, const acb_t x, const gr_ctx_t ctx)
{
    *res = arf_get_d(arb_midref(acb_realref(x)), ARF_RND_NEAR);
    return GR_SUCCESS;
}

truth_t
_gr_acb_is_zero(const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_zero(x))
        return T_TRUE;

    if (acb_contains_zero(x))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
_gr_acb_is_one(const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_one(x))
        return T_TRUE;

    if (arb_contains_zero(acb_imagref(x)) && arb_contains_si(acb_realref(x), 1))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
_gr_acb_is_neg_one(const acb_t x, const gr_ctx_t ctx)
{
    if (acb_equal_si(x, -1))
        return T_TRUE;

    if (arb_contains_zero(acb_imagref(x)) && arb_contains_si(acb_realref(x), -1))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
_gr_acb_equal(const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    if (acb_is_exact(x) && acb_equal(x, y))
        return T_TRUE;

    if (acb_overlaps(x, y))
        return T_UNKNOWN;

    return T_FALSE;
}

int
_gr_acb_set(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_set(res, x);
    return GR_SUCCESS;
}

int
_gr_acb_neg(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_neg(res, x);
    return GR_SUCCESS;
}

int
_gr_acb_add(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_add(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_add_si(acb_t res, const acb_t x, slong y, const gr_ctx_t ctx)
{
    acb_add_si(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_add_ui(acb_t res, const acb_t x, ulong y, const gr_ctx_t ctx)
{
    acb_add_ui(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_add_fmpz(acb_t res, const acb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    acb_add_fmpz(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_sub(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_sub(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_sub_si(acb_t res, const acb_t x, slong y, const gr_ctx_t ctx)
{
    acb_sub_si(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_sub_ui(acb_t res, const acb_t x, ulong y, const gr_ctx_t ctx)
{
    acb_sub_ui(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_sub_fmpz(acb_t res, const acb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    acb_sub_fmpz(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_mul(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul_si(acb_t res, const acb_t x, slong y, const gr_ctx_t ctx)
{
    acb_mul_si(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul_ui(acb_t res, const acb_t x, ulong y, const gr_ctx_t ctx)
{
    acb_mul_ui(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul_fmpz(acb_t res, const acb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    acb_mul_fmpz(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_addmul(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_addmul(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_submul(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    acb_submul(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul_two(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_mul_2exp_si(res, x, 1);
    return GR_SUCCESS;
}

int
_gr_acb_sqr(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_sqr(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mul_2exp_si(acb_t res, const acb_t x, slong y, const gr_ctx_t ctx)
{
    acb_mul_2exp_si(res, x, y);
    return GR_SUCCESS;
}

int
_gr_acb_mul_2exp_fmpz(acb_t res, const acb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    acb_mul_2exp_fmpz(res, x, y);
    return GR_SUCCESS;
}

int
_gr_acb_set_fmpz_2exp_fmpz(acb_t res, const fmpz_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arb_set_fmpz_2exp(acb_realref(res), x, y);
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_get_fmpz_2exp_fmpz(fmpz_t res1, fmpz_t res2, const acb_t x, const gr_ctx_t ctx)
{
    if (!acb_is_exact(x) || !acb_is_real(x))
        return GR_UNABLE;

    if (!acb_is_finite(x))
        return GR_DOMAIN;

    arf_get_fmpz_2exp(res1, res2, arb_midref(acb_realref(x)));
    return GR_SUCCESS;
}

int
_gr_acb_inv(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_inv(res, x, ACB_CTX_PREC(ctx));
        if (acb_is_finite(res))
            return GR_SUCCESS;
        else
            return GR_UNABLE;
    }
}

int
_gr_acb_div(acb_t res, const acb_t x, const acb_t y, const gr_ctx_t ctx)
{
    if (acb_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_div(res, x, y, ACB_CTX_PREC(ctx));

        if (acb_is_finite(res))
            return GR_SUCCESS;
        else
            return GR_UNABLE;
    }
}

int
_gr_acb_div_si(acb_t res, const acb_t x, slong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_div_si(res, x, y, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_acb_div_ui(acb_t res, const acb_t x, ulong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_div_ui(res, x, y, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_acb_div_fmpz(acb_t res, const acb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_div_fmpz(res, x, y, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

truth_t
_gr_acb_is_invertible(const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_zero(x))
        return T_FALSE;

    if (acb_contains_zero(x))
        return T_UNKNOWN;

    return T_TRUE;
}

int
_gr_acb_pow_ui(acb_t res, const acb_t x, ulong exp, const gr_ctx_t ctx)
{
    acb_pow_ui(res, x, exp, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_pow_si(acb_t res, const acb_t x, slong exp, const gr_ctx_t ctx)
{
    if (exp < 0 && acb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else if (exp < 0 && acb_contains_zero(x))
    {
        return GR_UNABLE;
    }
    else
    {
        fmpz_t t;
        fmpz_init_set_si(t, exp);
        acb_pow_fmpz(res, x, t, ACB_CTX_PREC(ctx));
        fmpz_clear(t);
        return GR_SUCCESS;
    }
}

int
_gr_acb_pow_fmpz(acb_t res, const acb_t x, const fmpz_t exp, const gr_ctx_t ctx)
{
    if (fmpz_sgn(exp) < 0 && acb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else if (fmpz_sgn(exp) < 0 && acb_contains_zero(x))
    {
        return GR_UNABLE;
    }
    else
    {
        acb_pow_fmpz(res, x, exp, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_acb_pow(acb_t res, const acb_t x, const acb_t exp, const gr_ctx_t ctx)
{
    if (acb_is_int(exp))
    {
        if (arf_sgn(arb_midref(acb_realref(exp))) < 0)
        {
            if (acb_is_zero(x))
                return GR_DOMAIN;

            if (acb_contains_zero(x))
                return GR_UNABLE;
        }

        acb_pow(res, x, exp, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (!acb_contains_zero(x) || arb_is_positive(acb_realref(exp)))
    {
        acb_pow(res, x, exp, ACB_CTX_PREC(ctx));

        if (!acb_is_finite(res))
            return GR_UNABLE;

        return GR_SUCCESS;
    }
    else if (acb_is_zero(x) && arb_is_negative(acb_realref(exp)))
    {
        return GR_DOMAIN;
    }
    else
    {
        return GR_UNABLE;
    }
}

int
_gr_acb_pow_fmpq(acb_t res, const acb_t x, const fmpq_t exp, const gr_ctx_t ctx)
{
    int status;
    acb_t t;
    acb_init(t);
    acb_set_fmpq(t, exp, ACB_CTX_PREC(ctx) + 20);
    status = _gr_acb_pow(res, x, t, ctx);
    acb_clear(t);
    return status;
}

truth_t
_gr_acb_is_square(const acb_t x, const gr_ctx_t ctx)
{
    return T_TRUE;
}

int
_gr_acb_sqrt(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_sqrt(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_rsqrt(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (!acb_contains_zero(x))
    {
        acb_rsqrt(res, x, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (acb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        return GR_UNABLE;
    }
}

int _arb_trunc(arb_t res, const arb_t x, slong prec);
int _arb_nint(arb_t res, const arb_t x, slong prec);

int
_gr_acb_floor(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    arb_floor(acb_realref(res), acb_realref(x), ACB_CTX_PREC(ctx));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_ceil(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    arb_ceil(acb_realref(res), acb_realref(x), ACB_CTX_PREC(ctx));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_trunc(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    arb_zero(acb_imagref(res));
    return _arb_trunc(acb_realref(res), acb_realref(x), ACB_CTX_PREC(ctx));
}

int
_gr_acb_nint(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    arb_zero(acb_imagref(res));
    return _arb_nint(acb_realref(res), acb_realref(x), ACB_CTX_PREC(ctx));
}

int
_gr_acb_i(acb_t res, const gr_ctx_t ctx)
{
    acb_onei(res);
    return GR_SUCCESS;
}

int
_gr_acb_abs(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_abs(acb_realref(res), x, ACB_CTX_PREC(ctx));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_conj(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_conj(res, x);
    return GR_SUCCESS;
}

int
_gr_acb_re(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    arb_set(acb_realref(res), acb_realref(x));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_im(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    arb_set(acb_realref(res), acb_imagref(x));
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_sgn(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_sgn(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_csgn(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_csgn(acb_realref(res), x);
    arb_zero(acb_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acb_pi(acb_t res, const gr_ctx_t ctx)
{
    acb_const_pi(res, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_exp(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_exp(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_log(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (acb_contains_zero(x))
    {
        if (acb_is_zero(x))
            return GR_DOMAIN;
        return GR_UNABLE;
    }

    acb_log(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_sin(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_sin(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_cos(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_cos(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_tan(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_tan(res, x, ACB_CTX_PREC(ctx));
    return acb_is_finite(res) ? GR_SUCCESS : GR_DOMAIN;
}

int
_gr_acb_atan(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (!arb_is_zero(acb_imagref(x)) && arb_contains_zero(acb_realref(x)))
    {
        if (arb_is_zero(acb_realref(x)) && (arb_is_one(acb_imagref(x)) || arb_equal_si(acb_imagref(x), -1)))
            return GR_DOMAIN;

        if (arb_contains_si(acb_imagref(x), 1) || arb_contains_si(acb_imagref(x), -1))
            return GR_UNABLE;
    }

    acb_atan(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_erf(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_erf(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_erfc(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_erfc(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_erfi(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_hypgeom_erfi(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_gamma(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_int(x) && arb_is_nonpositive(acb_realref(x)))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_gamma(res, x, ACB_CTX_PREC(ctx));
        return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
    }
}

int
_gr_acb_rgamma(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    acb_rgamma(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_lgamma(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_int(x) && arb_is_nonpositive(acb_realref(x)))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_lgamma(res, x, ACB_CTX_PREC(ctx));
        return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
    }
}

int
_gr_acb_digamma(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (acb_is_int(x) && arb_is_nonpositive(acb_realref(x)))
    {
        return GR_DOMAIN;
    }
    else
    {
        acb_digamma(res, x, ACB_CTX_PREC(ctx));
        return acb_is_finite(res) ? GR_SUCCESS : GR_UNABLE;
    }
}

int
_gr_acb_zeta(acb_t res, const acb_t x, const gr_ctx_t ctx)
{
    if (arb_contains_si(acb_realref(x), 1) && arb_contains_zero(acb_imagref(x)))
    {
        if (acb_is_one(x))
            return GR_DOMAIN;
        else
            return GR_UNABLE;
    }
    else
    {
        acb_zeta(res, x, ACB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_acb_vec_dot(acb_t res, const acb_t initial, int subtract, acb_srcptr vec1, acb_srcptr vec2, slong len, gr_ctx_t ctx)
{
    acb_dot(res, initial, subtract, vec1, 1, vec2, 1, len, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_vec_dot_rev(acb_t res, const acb_t initial, int subtract, acb_srcptr vec1, acb_srcptr vec2, slong len, gr_ctx_t ctx)
{
    acb_dot(res, initial, subtract, vec1, 1, vec2 + len - 1, -1, len, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_poly_mullow(acb_ptr res,
    acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    _acb_poly_mullow(res, poly1, len1, poly2, len2, n, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mat_mul(acb_mat_t res, const acb_mat_t x, const acb_mat_t y, gr_ctx_t ctx)
{
    acb_mat_mul(res, x, y, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_acb_mat_det(acb_t res, const acb_mat_t x, gr_ctx_t ctx)
{
    acb_mat_det(res, x, ACB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int _acb_methods_initialized = 0;

gr_static_method_table _acb_methods;

gr_method_tab_input _acb_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_acb_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_acb_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_acb_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_acb_swap},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_acb_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_acb_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_acb_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_acb_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_acb_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_acb_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_acb_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_acb_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_acb_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_acb_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_acb_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_acb_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_acb_set_fmpq},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_acb_set_other},
    {GR_METHOD_SET_D,           (gr_funcptr) _gr_acb_set_d},
    {GR_METHOD_GET_SI,          (gr_funcptr) _gr_acb_get_si},
    {GR_METHOD_GET_UI,          (gr_funcptr) _gr_acb_get_ui},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_acb_get_fmpz},
    {GR_METHOD_GET_D,           (gr_funcptr) _gr_acb_get_d},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_acb_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_acb_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_acb_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_acb_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_acb_add_fmpz},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_acb_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_acb_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_acb_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_acb_sub_fmpz},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_acb_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_acb_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_acb_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_acb_mul_fmpz},
    {GR_METHOD_ADDMUL,          (gr_funcptr) _gr_acb_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) _gr_acb_submul},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_acb_mul_two},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_acb_sqr},
    {GR_METHOD_MUL_2EXP_SI,        (gr_funcptr) _gr_acb_mul_2exp_si},
    {GR_METHOD_MUL_2EXP_FMPZ,      (gr_funcptr) _gr_acb_mul_2exp_fmpz},
    {GR_METHOD_SET_FMPZ_2EXP_FMPZ, (gr_funcptr) _gr_acb_set_fmpz_2exp_fmpz},
    {GR_METHOD_GET_FMPZ_2EXP_FMPZ, (gr_funcptr) _gr_acb_get_fmpz_2exp_fmpz},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_acb_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_acb_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_acb_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_acb_div_fmpz},
    {GR_METHOD_INV,             (gr_funcptr) _gr_acb_inv},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_acb_is_invertible},
    {GR_METHOD_POW,             (gr_funcptr) _gr_acb_pow},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_acb_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_acb_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_acb_pow_fmpz},
    {GR_METHOD_POW_FMPQ,        (gr_funcptr) _gr_acb_pow_fmpq},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_acb_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_acb_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_acb_rsqrt},
    {GR_METHOD_FLOOR,           (gr_funcptr) _gr_acb_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) _gr_acb_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) _gr_acb_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) _gr_acb_nint},
    {GR_METHOD_I,               (gr_funcptr) _gr_acb_i},
    {GR_METHOD_ABS,             (gr_funcptr) _gr_acb_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) _gr_acb_conj},
    {GR_METHOD_RE,              (gr_funcptr) _gr_acb_re},
    {GR_METHOD_IM,              (gr_funcptr) _gr_acb_im},
    {GR_METHOD_SGN,             (gr_funcptr) _gr_acb_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) _gr_acb_csgn},
    {GR_METHOD_PI,              (gr_funcptr) _gr_acb_pi},
    {GR_METHOD_EXP,             (gr_funcptr) _gr_acb_exp},
    {GR_METHOD_LOG,             (gr_funcptr) _gr_acb_log},
    {GR_METHOD_SIN,             (gr_funcptr) _gr_acb_sin},
    {GR_METHOD_COS,             (gr_funcptr) _gr_acb_cos},
    {GR_METHOD_TAN,             (gr_funcptr) _gr_acb_tan},
    {GR_METHOD_ATAN,            (gr_funcptr) _gr_acb_atan},
    {GR_METHOD_ERF,             (gr_funcptr) _gr_acb_erf},
    {GR_METHOD_ERFI,            (gr_funcptr) _gr_acb_erfi},
    {GR_METHOD_ERFC,            (gr_funcptr) _gr_acb_erfc},
    {GR_METHOD_GAMMA,           (gr_funcptr) _gr_acb_gamma},
    {GR_METHOD_RGAMMA,          (gr_funcptr) _gr_acb_rgamma},
    {GR_METHOD_LGAMMA,          (gr_funcptr) _gr_acb_lgamma},
    {GR_METHOD_DIGAMMA,         (gr_funcptr) _gr_acb_digamma},
    {GR_METHOD_ZETA,            (gr_funcptr) _gr_acb_zeta},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) _gr_acb_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _gr_acb_vec_dot_rev},
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_acb_poly_mullow},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_acb_mat_mul},
    {GR_METHOD_MAT_DET,         (gr_funcptr) _gr_acb_mat_det},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_complex_acb(gr_ctx_t ctx, slong prec)
{
    ctx->which_ring = GR_CTX_CC_ACB;
    ctx->sizeof_elem = sizeof(acb_struct);
    ctx->size_limit = WORD_MAX;

    gr_ctx_arb_set_prec(ctx, prec);

    ACB_CTX_PREC(ctx) = prec;

    ctx->methods = _acb_methods;

    if (!_acb_methods_initialized)
    {
        gr_method_tab_init(_acb_methods, _acb_methods_input);
        _acb_methods_initialized = 1;
    }
}
