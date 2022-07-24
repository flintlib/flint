#include "arb.h"
#include "arb_poly.h"
#include "arb_mat.h"
#include "qqbar.h"
#include "gr.h"

typedef struct
{
    slong prec;
}
gr_arb_ctx;

#define ARB_CTX_PREC(ring_ctx) (((gr_arb_ctx *)((ring_ctx)->elem_ctx))->prec)

void gr_ctx_arb_set_prec(gr_ctx_t ctx, slong prec)
{
    prec = FLINT_MAX(prec, 2);
    prec = FLINT_MIN(prec, WORD_MAX / 8);

    ARB_CTX_PREC(ctx) = prec;
}

slong gr_ctx_arb_get_prec(gr_ctx_t ctx)
{
    return ARB_CTX_PREC(ctx);
}

int
_gr_arb_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Real numbers (arb, prec = ");
    gr_stream_write_si(out, ARB_CTX_PREC(ctx));
    gr_stream_write(out, ")");
    return GR_SUCCESS;
}

void
_gr_arb_init(arb_t x, const gr_ctx_t ctx)
{
    arb_init(x);
}

void
_gr_arb_clear(arb_t x, const gr_ctx_t ctx)
{
    arb_clear(x);
}

void
_gr_arb_swap(arb_t x, arb_t y, const gr_ctx_t ctx)
{
    arb_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

/* todo: limits */
int
_gr_arb_randtest(arb_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    arb_randtest(res, state, ARB_CTX_PREC(ctx), 10);
    return GR_SUCCESS;
}

/* todo */
int
_gr_arb_write(gr_stream_t out, const arb_t x, const gr_ctx_t ctx)
{
    gr_stream_write_free(out, arb_get_str(x, ARB_CTX_PREC(ctx) * 0.30102999566398 + 1, 0));
    return GR_SUCCESS;
}

int
_gr_arb_zero(arb_t x, const gr_ctx_t ctx)
{
    arb_zero(x);
    return GR_SUCCESS;
}

int
_gr_arb_one(arb_t x, const gr_ctx_t ctx)
{
    arb_one(x);
    return GR_SUCCESS;
}

int
_gr_arb_set_si(arb_t res, slong v, const gr_ctx_t ctx)
{
    arb_set_si(res, v);
    arb_set_round(res, res, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_set_ui(arb_t res, ulong v, const gr_ctx_t ctx)
{
    arb_set_ui(res, v);
    arb_set_round(res, res, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_set_fmpz(arb_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    arb_set_round_fmpz(res, v, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_set_fmpq(arb_t res, const fmpq_t v, const gr_ctx_t ctx)
{
    arb_set_fmpq(res, v, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_set_str(arb_t res, const char * x, const gr_ctx_t ctx)
{
    if (arb_set_str(res, x, ARB_CTX_PREC(ctx)))
        return GR_DOMAIN;

    return GR_SUCCESS;
}

int 
_gr_ca_get_arb_with_prec(arb_t res, gr_srcptr x, gr_ctx_t x_ctx, slong prec);

int
_gr_arb_set_other(arb_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            return _gr_arb_set_fmpz(res, x, ctx);

        case GR_CTX_FMPQ:
            return _gr_arb_set_fmpq(res, x, ctx);

        case GR_CTX_REAL_ALGEBRAIC_QQBAR:
            qqbar_get_arb(res, x, ARB_CTX_PREC(ctx));
            return GR_SUCCESS;

        case GR_CTX_COMPLEX_ALGEBRAIC_QQBAR:
            if (qqbar_is_real(x))
            {
                qqbar_get_arb(res, x, ARB_CTX_PREC(ctx));
                return GR_SUCCESS;
            }
            else
            {
                return GR_DOMAIN;
            }

        case GR_CTX_RR_CA:
        case GR_CTX_REAL_ALGEBRAIC_CA:
        case GR_CTX_CC_CA:
        case GR_CTX_COMPLEX_ALGEBRAIC_CA:
            return _gr_ca_get_arb_with_prec(res, x, x_ctx, ARB_CTX_PREC(ctx));

        case GR_CTX_REAL_FLOAT_ARF:
            if (arf_is_finite(x))
            {
                arb_set_arf(res, x);
                arb_set_round(res, res, ARB_CTX_PREC(ctx));
                return GR_SUCCESS;
            }
            else
            {
                return GR_DOMAIN;
            }

        case GR_CTX_RR_ARB:
            arb_set_round(res, x, ARB_CTX_PREC(ctx));
            return GR_SUCCESS;

        case GR_CTX_CC_ACB:
            if (arb_is_zero(acb_imagref((acb_srcptr) x)))
            {
                arb_set_round(res, x, ARB_CTX_PREC(ctx));
                return GR_SUCCESS;
            }
            else if (arb_contains_zero(acb_imagref((acb_srcptr) x)))
            {
                return GR_UNABLE;
            }
            else
            {
                return GR_DOMAIN;
            }
    }

    return GR_UNABLE;
}

/* xxx: assumes that ctx are not read */
int _gr_arf_get_fmpz(fmpz_t res, const arf_t x, const gr_ctx_t ctx);
int _gr_arf_get_si(slong * res, const arf_t x, const gr_ctx_t ctx);
int _gr_arf_get_ui(ulong * res, const arf_t x, const gr_ctx_t ctx);

int
_gr_arb_get_fmpz(fmpz_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (!arb_is_int(x))
    {
        if (arb_contains_int(x))
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }

    return _gr_arf_get_fmpz(res, arb_midref(x), NULL);
}

int
_gr_arb_get_si(slong * res, const arb_t x, const gr_ctx_t ctx)
{
    if (!arb_is_int(x))
    {
        if (arb_contains_int(x))
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }

    return _gr_arf_get_si(res, arb_midref(x), NULL);
}

int
_gr_arb_get_ui(ulong * res, const arb_t x, const gr_ctx_t ctx)
{
    if (!arb_is_int(x))
    {
        if (arb_contains_int(x))
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }

    return _gr_arf_get_ui(res, arb_midref(x), NULL);
}

int
_gr_arb_get_d(double * res, const arb_t x, const gr_ctx_t ctx)
{
    *res = arf_get_d(arb_midref(x), ARF_RND_NEAR);
    return GR_SUCCESS;
}

truth_t
_gr_arb_is_zero(const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_zero(x))
        return T_TRUE;

    if (arb_contains_zero(x))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
_gr_arb_is_one(const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_one(x))
        return T_TRUE;

    if (arb_contains_si(x, 1))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
_gr_arb_is_neg_one(const arb_t x, const gr_ctx_t ctx)
{
    if (arb_equal_si(x, -1))
        return T_TRUE;

    if (arb_contains_si(x, -1))
        return T_UNKNOWN;

    return T_FALSE;
}

truth_t
_gr_arb_equal(const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    if (arb_is_exact(x) && arb_equal(x, y))
        return T_TRUE;

    if (arb_overlaps(x, y))
        return T_UNKNOWN;

    return T_FALSE;
}

int
_gr_arb_set(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_set(res, x);
    return GR_SUCCESS;
}

int
_gr_arb_neg(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_neg(res, x);
    return GR_SUCCESS;
}

int
_gr_arb_add(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    arb_add(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_add_si(arb_t res, const arb_t x, slong y, const gr_ctx_t ctx)
{
    arb_add_si(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_add_ui(arb_t res, const arb_t x, ulong y, const gr_ctx_t ctx)
{
    arb_add_ui(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_add_fmpz(arb_t res, const arb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arb_add_fmpz(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_sub(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    arb_sub(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_sub_si(arb_t res, const arb_t x, slong y, const gr_ctx_t ctx)
{
    arb_sub_si(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_sub_ui(arb_t res, const arb_t x, ulong y, const gr_ctx_t ctx)
{
    arb_sub_ui(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_sub_fmpz(arb_t res, const arb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arb_sub_fmpz(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_mul(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    arb_mul(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_mul_si(arb_t res, const arb_t x, slong y, const gr_ctx_t ctx)
{
    arb_mul_si(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_mul_ui(arb_t res, const arb_t x, ulong y, const gr_ctx_t ctx)
{
    arb_mul_ui(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_mul_fmpz(arb_t res, const arb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arb_mul_fmpz(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_addmul(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    arb_addmul(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_submul(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    arb_submul(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_mul_two(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_mul_2exp_si(res, x, 1);
    return GR_SUCCESS;
}

int
_gr_arb_sqr(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_sqr(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_inv(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        arb_inv(res, x, ARB_CTX_PREC(ctx));
        if (arb_is_finite(res))
            return GR_SUCCESS;
        else
            return GR_UNABLE;
    }
}

int
_gr_arb_div(arb_t res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    if (arb_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        arb_div(res, x, y, ARB_CTX_PREC(ctx));

        if (arb_is_finite(res))
            return GR_SUCCESS;
        else
            return GR_UNABLE;
    }
}

int
_gr_arb_div_si(arb_t res, const arb_t x, slong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        arb_div_si(res, x, y, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_arb_div_ui(arb_t res, const arb_t x, ulong y, const gr_ctx_t ctx)
{
    if (y == 0)
    {
        return GR_DOMAIN;
    }
    else
    {
        arb_div_ui(res, x, y, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_arb_div_fmpz(arb_t res, const arb_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    if (fmpz_is_zero(y))
    {
        return GR_DOMAIN;
    }
    else
    {
        arb_div_fmpz(res, x, y, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

truth_t
_gr_arb_is_invertible(const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_zero(x))
        return T_FALSE;

    if (arb_contains_zero(x))
        return T_UNKNOWN;

    return T_TRUE;
}

int
_gr_arb_pow_ui(arb_t res, const arb_t x, ulong exp, const gr_ctx_t ctx)
{
    arb_pow_ui(res, x, exp, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_pow_si(arb_t res, const arb_t x, slong exp, const gr_ctx_t ctx)
{
    if (exp < 0 && arb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else if (exp < 0 && arb_contains_zero(x))
    {
        return GR_UNABLE;
    }
    else
    {
        fmpz_t t;
        fmpz_init_set_si(t, exp);
        arb_pow_fmpz(res, x, t, ARB_CTX_PREC(ctx));
        fmpz_clear(t);
        return GR_SUCCESS;
    }
}

int
_gr_arb_pow_fmpz(arb_t res, const arb_t x, const fmpz_t exp, const gr_ctx_t ctx)
{
    if (fmpz_sgn(exp) < 0 && arb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else if (fmpz_sgn(exp) < 0 && arb_contains_zero(x))
    {
        return GR_UNABLE;
    }
    else
    {
        arb_pow_fmpz(res, x, exp, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_arb_pow_fmpq(arb_t res, const arb_t x, const fmpq_t exp, const gr_ctx_t ctx)
{
    if (fmpq_sgn(exp) < 0 && arb_is_zero(x))
    {
        return GR_DOMAIN;
    }
    else if (fmpq_sgn(exp) < 0 && arb_contains_zero(x))
    {
        return GR_UNABLE;
    }
    else
    {
        if (fmpz_is_one(fmpq_denref(exp)) || arb_is_nonnegative(x))
        {
            arb_pow_fmpq(res, x, exp, ARB_CTX_PREC(ctx));
            return GR_SUCCESS;
        }
        else if (arb_is_negative(x))
        {
            return GR_DOMAIN;
        }
        else
        {
            return GR_UNABLE;
        }
    }
}

int
_gr_arb_pow(arb_t res, const arb_t x, const arb_t exp, const gr_ctx_t ctx)
{
    if (arb_is_int(exp))
    {
        if (arf_sgn(arb_midref(exp)) < 0)
        {
            if (arb_is_zero(x))
                return GR_DOMAIN;

            if (arb_contains_zero(x))
                return GR_UNABLE;
        }

        arb_pow(res, x, exp, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (arb_is_positive(x) || (arb_is_nonnegative(x) && arb_is_nonnegative(exp)))
    {
        arb_pow(res, x, exp, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (arb_is_zero(x) && arb_is_negative(exp))
    {
        return GR_DOMAIN;
    }
    else if (arb_is_negative(x) && !arb_contains_int(exp))
    {
        return GR_DOMAIN;
    }
    else
    {
        return GR_UNABLE;
    }
}

truth_t
_gr_arb_is_square(const arb_t x, const gr_ctx_t ctx)
{
    return T_TRUE;
}

int
_gr_arb_sqrt(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_nonnegative(x))
    {
        arb_sqrt(res, x, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (arb_is_negative(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        return GR_UNABLE;
    }
}

int
_gr_arb_rsqrt(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_positive(x))
    {
        arb_rsqrt(res, x, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }
    else if (arb_is_nonpositive(x))
    {
        return GR_DOMAIN;
    }
    else
    {
        return GR_UNABLE;
    }
}

int
_gr_arb_floor(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_floor(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_ceil(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_ceil(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_arb_trunc(arb_t res, const arb_t x, slong prec)
{
    if (arb_contains_zero(x))
    {
        arb_t a;
        arb_init(a);

        mag_one(arb_radref(a));

        if (arb_contains(a, x))
        {
            arb_zero(res);
        }
        else
        {
            arb_t b;
            arb_init(b);
            arb_floor(a, x, prec);
            arb_ceil(b, x, prec);
            arb_union(res, a, b, prec);
            arb_clear(b);
        }

        arb_clear(a);
    }
    else if (arf_sgn(arb_midref(x)) > 0)
    {
        arb_floor(res, x, prec);
    }
    else
    {
        arb_ceil(res, x, prec);
    }

    return GR_SUCCESS;
}

int
_arb_nint(arb_t res, const arb_t x, slong prec)
{
    if (arb_is_int(x))
    {
        arb_set(res, x);
    }
    else
    {
        arb_t t, u;
        arb_init(t);
        arb_init(u);

        arb_set_d(t, 0.5);
        arb_add(t, x, t, prec);

        arb_mul_2exp_si(u, x, 1);
        arb_sub_ui(u, u, 1, prec);
        arb_mul_2exp_si(u, u, -2);

        arb_floor(res, t, prec);

        /* nint(x) = floor(x+0.5) - isint((2*x-1)/4) */

        if (arb_is_int(u))
        {
            arb_sub_ui(res, res, 1, prec);
        }
        else if (arb_contains_int(u))
        {
            arf_one(arb_midref(u));
            mag_one(arb_radref(u));
            arb_mul_2exp_si(u, u, -1);
            arb_sub_ui(res, res, 1, prec);
        }
    }

    return GR_SUCCESS;
}

int
_gr_arb_trunc(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    return _arb_trunc(res, x, ARB_CTX_PREC(ctx));
}

int
_gr_arb_nint(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    return _arb_nint(res, x, ARB_CTX_PREC(ctx));
}

int
_gr_arb_abs(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_abs(res, x);
    return GR_SUCCESS;
}

int
_gr_arb_conj(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_set(res, x);
    return GR_SUCCESS;
}

int
_gr_arb_im(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_zero(res);
    return GR_SUCCESS;
}

int
_gr_arb_sgn(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_sgn(res, x);
    return GR_SUCCESS;
}

int
_gr_arb_cmp(int * res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    if ((arb_is_exact(x) && arb_is_exact(y)) || !arb_overlaps(x, y))
    {
        *res = arf_cmp(arb_midref(x), arb_midref(y));
        return GR_SUCCESS;
    }
    else
    {
        *res = 0;
        return GR_UNABLE;
    }
}

int
_gr_arb_cmpabs(int * res, const arb_t x, const arb_t y, const gr_ctx_t ctx)
{
    arb_t t, u;

    *t = *x;
    *u = *y;

    if (arf_sgn(arb_midref(t)) < 0)
        ARF_NEG(arb_midref(t));

    if (arf_sgn(arb_midref(u)) < 0)
        ARF_NEG(arb_midref(u));

    return _gr_arb_cmp(res, t, u, ctx);
}

int
_gr_arb_pi(arb_t res, const gr_ctx_t ctx)
{
    arb_const_pi(res, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_exp(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_exp(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_log(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    if (arb_is_positive(x))
    {
        arb_log(res, x, ARB_CTX_PREC(ctx));
        return GR_SUCCESS;
    }

    if (arb_is_nonpositive(x))
        return GR_DOMAIN;

    return GR_UNABLE;
}

int
_gr_arb_sin(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_sin(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_cos(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_cos(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_atan(arb_t res, const arb_t x, const gr_ctx_t ctx)
{
    arb_atan(res, x, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_vec_dot(arb_t res, const arb_t initial, int subtract, arb_srcptr vec1, arb_srcptr vec2, slong len, gr_ctx_t ctx)
{
    arb_dot(res, initial, subtract, vec1, 1, vec2, 1, len, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_vec_dot_rev(arb_t res, const arb_t initial, int subtract, arb_srcptr vec1, arb_srcptr vec2, slong len, gr_ctx_t ctx)
{
    arb_dot(res, initial, subtract, vec1, 1, vec2 + len - 1, -1, len, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_poly_mullow(arb_ptr res,
    arb_srcptr poly1, slong len1,
    arb_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    _arb_poly_mullow(res, poly1, len1, poly2, len2, n, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_mat_mul(arb_mat_t res, const arb_mat_t x, const arb_mat_t y, gr_ctx_t ctx)
{
    arb_mat_mul(res, x, y, ARB_CTX_PREC(ctx));
    return GR_SUCCESS;
}

int
_gr_arb_ctx_clear(gr_ctx_t ctx)
{
    flint_free(ctx->elem_ctx);
    return GR_SUCCESS;
}

int _arb_methods_initialized = 0;

gr_static_method_table _arb_methods;

gr_method_tab_input _arb_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_arb_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_arb_ctx_write},
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
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_arb_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_arb_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_arb_swap},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_arb_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_arb_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_arb_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_arb_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_arb_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_arb_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_arb_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_arb_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_arb_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_arb_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_arb_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_arb_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_arb_set_fmpq},
    {GR_METHOD_SET_STR,         (gr_funcptr) _gr_arb_set_str},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_arb_set_other},
    {GR_METHOD_GET_SI,          (gr_funcptr) _gr_arb_get_si},
    {GR_METHOD_GET_UI,          (gr_funcptr) _gr_arb_get_ui},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_arb_get_fmpz},
    {GR_METHOD_GET_D,           (gr_funcptr) _gr_arb_get_d},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_arb_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_arb_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_arb_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_arb_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_arb_add_fmpz},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_arb_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_arb_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_arb_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_arb_sub_fmpz},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_arb_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_arb_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_arb_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_arb_mul_fmpz},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_arb_mul_two},
    {GR_METHOD_ADDMUL,          (gr_funcptr) _gr_arb_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) _gr_arb_submul},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_arb_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_arb_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_arb_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_arb_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_arb_div_fmpz},
    {GR_METHOD_INV,             (gr_funcptr) _gr_arb_inv},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_arb_is_invertible},
    {GR_METHOD_POW,             (gr_funcptr) _gr_arb_pow},
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_arb_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_arb_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_arb_pow_fmpz},
    {GR_METHOD_POW_FMPQ,        (gr_funcptr) _gr_arb_pow_fmpq},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_arb_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_arb_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_arb_rsqrt},
    {GR_METHOD_FLOOR,           (gr_funcptr) _gr_arb_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) _gr_arb_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) _gr_arb_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) _gr_arb_nint},
    {GR_METHOD_ABS,             (gr_funcptr) _gr_arb_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) _gr_arb_conj},
    {GR_METHOD_RE,              (gr_funcptr) _gr_arb_set},
    {GR_METHOD_IM,              (gr_funcptr) _gr_arb_im},
    {GR_METHOD_SGN,             (gr_funcptr) _gr_arb_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) _gr_arb_sgn},
    {GR_METHOD_CMP,             (gr_funcptr) _gr_arb_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) _gr_arb_cmpabs},
    {GR_METHOD_I,               (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_PI,              (gr_funcptr) _gr_arb_pi},
    {GR_METHOD_EXP,             (gr_funcptr) _gr_arb_exp},
    {GR_METHOD_LOG,             (gr_funcptr) _gr_arb_log},
    {GR_METHOD_SIN,             (gr_funcptr) _gr_arb_sin},
    {GR_METHOD_COS,             (gr_funcptr) _gr_arb_cos},
    {GR_METHOD_ATAN,            (gr_funcptr) _gr_arb_atan},
    {GR_METHOD_VEC_DOT,         (gr_funcptr) _gr_arb_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _gr_arb_vec_dot_rev},
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_arb_poly_mullow},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_arb_mat_mul},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_real_arb(gr_ctx_t ctx, slong prec)
{
    ctx->flags = 0;
    ctx->which_ring = GR_CTX_RR_ARB;
    ctx->sizeof_elem = sizeof(arb_struct);
    ctx->elem_ctx = flint_malloc(sizeof(gr_arb_ctx));
    ctx->size_limit = WORD_MAX;

    gr_ctx_arb_set_prec(ctx, prec);

    ctx->methods = _arb_methods;

    if (!_arb_methods_initialized)
    {
        gr_method_tab_init(_arb_methods, _arb_methods_input);
        _arb_methods_initialized = 1;
    }
}
