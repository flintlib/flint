#include "arf.h"
#include "acb.h"
#include "gr.h"

typedef struct
{
    slong prec;
    arf_rnd_t rnd;
}
gr_arf_ctx;

#define ARF_CTX_PREC(ring_ctx) (((gr_arf_ctx *)((ring_ctx)->elem_ctx))->prec)
#define ARF_CTX_RND(ring_ctx) (((gr_arf_ctx *)((ring_ctx)->elem_ctx))->rnd)

void gr_ctx_arf_set_prec(gr_ctx_t ctx, slong prec)
{
    prec = FLINT_MAX(prec, 2);
    prec = FLINT_MIN(prec, WORD_MAX / 8);

    ARF_CTX_PREC(ctx) = prec;
}

slong gr_ctx_arf_get_prec(gr_ctx_t ctx)
{
    return ARF_CTX_PREC(ctx);
}

int
_gr_arf_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Floating-point numbers (arf, prec = ");
    gr_stream_write_si(out, ARF_CTX_PREC(ctx));
    gr_stream_write(out, ", radix = 2, +inf, -inf, nan)");
    return GR_SUCCESS;
}

void
_gr_arf_init(arf_t x, const gr_ctx_t ctx)
{
    arf_init(x);
}

void
_gr_arf_clear(arf_t x, const gr_ctx_t ctx)
{
    arf_clear(x);
}

void
_gr_arf_swap(arf_t x, arf_t y, const gr_ctx_t ctx)
{
    arf_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

/* todo: limits */
int
_gr_arf_randtest(arf_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    arf_randtest(res, state, ARF_CTX_PREC(ctx), 10);
    return GR_SUCCESS;
}

/* todo */
int
_gr_arf_write(gr_stream_t out, const arf_t x, const gr_ctx_t ctx)
{
    gr_stream_write_free(out, arf_get_str(x, ARF_CTX_PREC(ctx) * 0.30102999566398 + 1));
    return GR_SUCCESS;
}

int
_gr_arf_zero(arf_t x, const gr_ctx_t ctx)
{
    arf_zero(x);
    return GR_SUCCESS;
}

int
_gr_arf_one(arf_t x, const gr_ctx_t ctx)
{
    arf_one(x);
    return GR_SUCCESS;
}

int
_gr_arf_set_si(arf_t res, slong v, const gr_ctx_t ctx)
{
    arf_set_round_si(res, v, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_set_ui(arf_t res, ulong v, const gr_ctx_t ctx)
{
    arf_set_round_ui(res, v, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_set_fmpz(arf_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    arf_set_round_fmpz(res, v, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_set_fmpq(arf_t res, const fmpq_t v, const gr_ctx_t ctx)
{
    arf_set_fmpq(res, v, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_set_d(arf_t res, double x, const gr_ctx_t ctx)
{
    arf_set_d(res, x);
    return GR_SUCCESS;
}

/* todo: set_round? */
int
_gr_arf_set(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    arf_set(res, x);
    return GR_SUCCESS;
}

int
_gr_arf_set_str(arf_t res, const char * x, const gr_ctx_t ctx)
{
    arb_t t;
    arb_init(t);
    arb_set_str(t, x, ARF_CTX_PREC(ctx) + 20);
    arf_set_round(res, arb_midref(t), ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    arb_clear(t);
    return GR_SUCCESS;
}


int
_gr_arf_set_other(arf_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            return _gr_arf_set_fmpz(res, x, ctx);

        case GR_CTX_FMPQ:
            return _gr_arf_set_fmpq(res, x, ctx);

        case GR_CTX_REAL_FLOAT_ARF:
            return _gr_arf_set(res, x, ctx);

        case GR_CTX_RR_ARB:
            return _gr_arf_set(res, arb_midref((arb_srcptr) x), ctx);

        default:
            {
                gr_ctx_t cctx;
                acb_t z;
                int status;

                gr_ctx_init_complex_acb(cctx, ARF_CTX_PREC(ctx) + 20);
                acb_init(z);

                status = gr_set_other(z, x, x_ctx, cctx);

                if (status == GR_SUCCESS)
                {
                    if (acb_is_real(z))
                        status = _gr_arf_set(res, arb_midref(acb_realref(z)), ctx);
                    else
                        status = GR_DOMAIN;
                }

                acb_clear(z);
                gr_ctx_clear(cctx);

                return status;
            }
    }
}

int
_gr_arf_get_fmpz(fmpz_t res, const arf_t x, const gr_ctx_t ctx)
{
    if (!arf_is_int(x))
        return GR_DOMAIN;

    /* todo: detect mpz overflow */
    if (arf_cmpabs_2exp_si(x, WORD_MAX) >= 0)
        return GR_UNABLE;

    arf_get_fmpz(res, x, ARF_RND_DOWN);
    return GR_SUCCESS;
}

int
_gr_arf_get_si(slong * res, const arf_t x, const gr_ctx_t ctx)
{    fmpz_t t;

    if (!arf_is_int(x))
        return GR_DOMAIN;

    if (arf_cmp_si(x, WORD_MIN) < 0 || arf_cmp_si(x, WORD_MAX) > 0)
        return GR_DOMAIN;
    fmpz_init(t);
    arf_get_fmpz(t, x, ARF_RND_DOWN);
    *res = fmpz_get_si(t);
    fmpz_clear(t);

    return GR_SUCCESS;
}

int
_gr_arf_get_ui(ulong * res, const arf_t x, const gr_ctx_t ctx)
{
    fmpz_t t;

    if (!arf_is_int(x))
        return GR_DOMAIN;

    if (arf_sgn(x) < 0 || arf_cmp_ui(x, UWORD_MAX) > 0)
        return GR_DOMAIN;

    fmpz_init(t);
    arf_get_fmpz(t, x, ARF_RND_DOWN);
    *res = fmpz_get_ui(t);
    fmpz_clear(t);

    return GR_SUCCESS;
}

int
_gr_arf_get_d(double * res, const arf_t x, const gr_ctx_t ctx)
{
    *res = arf_get_d(x, ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_arf_is_zero(const arf_t x, const gr_ctx_t ctx)
{
    return arf_is_zero(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_arf_is_one(const arf_t x, const gr_ctx_t ctx)
{
    return arf_is_one(x) ? T_TRUE : T_FALSE;
}

truth_t
_gr_arf_is_neg_one(const arf_t x, const gr_ctx_t ctx)
{
     return arf_equal_si(x, -1) ? T_TRUE : T_FALSE;
}

truth_t
_gr_arf_equal(const arf_t x, const arf_t y, const gr_ctx_t ctx)
{
    if (arf_is_nan(x) || arf_is_nan(y))
        return T_FALSE;

    return arf_equal(x, y) ? T_TRUE : T_FALSE;
}

/* todo: neg_round? */
int
_gr_arf_neg(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    arf_neg(res, x);
    return GR_SUCCESS;
}

int
_gr_arf_add(arf_t res, const arf_t x, const arf_t y, const gr_ctx_t ctx)
{
    arf_add(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_add_si(arf_t res, const arf_t x, slong y, const gr_ctx_t ctx)
{
    arf_add_si(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_add_ui(arf_t res, const arf_t x, ulong y, const gr_ctx_t ctx)
{
    arf_add_ui(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_add_fmpz(arf_t res, const arf_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arf_add_fmpz(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_sub(arf_t res, const arf_t x, const arf_t y, const gr_ctx_t ctx)
{
    arf_sub(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_sub_si(arf_t res, const arf_t x, slong y, const gr_ctx_t ctx)
{
    arf_sub_si(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_sub_ui(arf_t res, const arf_t x, ulong y, const gr_ctx_t ctx)
{
    arf_sub_ui(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_sub_fmpz(arf_t res, const arf_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arf_sub_fmpz(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_mul(arf_t res, const arf_t x, const arf_t y, const gr_ctx_t ctx)
{
    arf_mul(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_mul_si(arf_t res, const arf_t x, slong y, const gr_ctx_t ctx)
{
    arf_mul_si(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_mul_ui(arf_t res, const arf_t x, ulong y, const gr_ctx_t ctx)
{
    arf_mul_ui(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_mul_fmpz(arf_t res, const arf_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arf_mul_fmpz(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_addmul(arf_t res, const arf_t x, const arf_t y, const gr_ctx_t ctx)
{
    arf_addmul(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_submul(arf_t res, const arf_t x, const arf_t y, const gr_ctx_t ctx)
{
    arf_submul(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_mul_two(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    arf_mul_2exp_si(res, x, 1);
    return GR_SUCCESS;
}

int
_gr_arf_sqr(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    arf_mul(res, x, x, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_inv(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    /* todo */
    arf_ui_div(res, 1, x, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_div(arf_t res, const arf_t x, const arf_t y, const gr_ctx_t ctx)
{
    arf_div(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_div_si(arf_t res, const arf_t x, slong y, const gr_ctx_t ctx)
{
    arf_div_si(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_div_ui(arf_t res, const arf_t x, ulong y, const gr_ctx_t ctx)
{
    arf_div_ui(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_div_fmpz(arf_t res, const arf_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arf_div_fmpz(res, x, y, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_sqrt(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    arf_sqrt(res, x, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_abs(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    arf_abs(res, x);
    return GR_SUCCESS;
}

int
_gr_arf_conj(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    arf_set(res, x);
    return GR_SUCCESS;
}

int
_gr_arf_im(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    arf_zero(res);
    return GR_SUCCESS;
}

/* todo: sign of nan? */
int
_gr_arf_sgn(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    arf_set_si(res, arf_sgn(x));
    return GR_SUCCESS;
}

int
_gr_arf_rsqrt(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    arf_rsqrt(res, x, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_floor(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    arf_floor(res, x);
    return GR_SUCCESS;
}

int
_gr_arf_ceil(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    arf_ceil(res, x);
    return GR_SUCCESS;
}

int
_gr_arf_trunc(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    if (arf_is_int(x) || arf_is_special(x))
    {
        arf_set(res, x);
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        arf_get_fmpz(t, x, ARF_RND_DOWN);
        arf_set_fmpz(res, t);
        fmpz_clear(t);
    }

    return GR_SUCCESS;
}

int
_gr_arf_nint(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    if (arf_is_int(x) || arf_is_special(x))
    {
        arf_set(res, x);
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        arf_get_fmpz(t, x, ARF_RND_NEAR);
        arf_set_fmpz(res, t);
        fmpz_clear(t);
    }

    return GR_SUCCESS;
}

/* todo: handling nan */
int
_gr_arf_cmp(int * res, const arf_t x, const arf_t y, const gr_ctx_t ctx)
{
    *res = arf_cmp(x, y);
    return GR_SUCCESS;
}

int
_gr_arf_cmpabs(int * res, const arf_t x, const arf_t y, const gr_ctx_t ctx)
{
    *res = arf_cmpabs(x, y);
    return GR_SUCCESS;
}

#define ARF_FUNC_VIA_ARB(res, arb_func, x) \
    arb_t r, t; \
    slong prec, wp, extra; \
    int status = GR_SUCCESS; \
    prec = ARF_CTX_PREC(ctx); \
    arb_init(r); \
    *arb_midref(t) = *x; \
    mag_init(arb_radref(t)); \
    for (extra = 10 + prec * 0.01; ; extra += FLINT_MAX(extra, 32)) \
    { \
        wp = prec + extra; \
        if (wp > 10 * prec + 1000) \
        { \
            status = GR_UNABLE; \
            arf_nan(res); \
            break; \
        } \
        arb_func(r, t, wp); \
        if (arb_rel_accuracy_bits(r) >= prec) \
        { \
            arf_set_round(res, arb_midref(r), prec, ARF_CTX_RND(ctx)); \
            break; \
        } \
    } \
    arb_clear(r); \
    return status; \

#define ARF_FUNC2_VIA_ARB(res, arb_func, x, y) \
    arb_t r, t, u; \
    slong prec, wp, extra; \
    int status = GR_SUCCESS; \
    prec = ARF_CTX_PREC(ctx); \
    arb_init(r); \
    *arb_midref(t) = *x; \
    mag_init(arb_radref(t)); \
    *arb_midref(u) = *y; \
    mag_init(arb_radref(u)); \
    for (extra = 10 + prec * 0.01; ; extra += FLINT_MAX(extra, 32)) \
    { \
        wp = prec + extra; \
        if (wp > 10 * prec + 1000) \
        { \
            status = GR_UNABLE; \
            arf_nan(res); \
            break; \
        } \
        arb_func(r, t, u, wp); \
        if (arb_rel_accuracy_bits(r) >= prec) \
        { \
            arf_set_round(res, arb_midref(r), prec, ARF_CTX_RND(ctx)); \
            break; \
        } \
    } \
    arb_clear(r); \
    return status; \

/* todo: lots of special cases */
int
_gr_arf_pow(arf_t res, const arf_t x, const arf_t y, const gr_ctx_t ctx)
{
    if (!arf_is_finite(x) || !arf_is_finite(y)
        || (arf_is_zero(x) && arf_sgn(y) < 0)
        || (arf_sgn(x) < 0 && !arf_is_int(y)))
    {
        arf_nan(res);
        return GR_SUCCESS;
    }
    else
    {
        ARF_FUNC2_VIA_ARB(res, arb_pow, x, y)
    }
}

int
_gr_arf_pi(arf_t res, const gr_ctx_t ctx)
{
    arb_t t;
    arb_init(t);
    arb_const_pi(t, ARF_CTX_PREC(ctx) + 30);
    arf_set_round(res, arb_midref(t), ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    arb_clear(t);
    return GR_SUCCESS;
}

int
_gr_arf_exp(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_one(res);
        else if (arf_is_pos_inf(x))
            arf_pos_inf(res);
        else if (arf_is_neg_inf(x))
            arf_zero(res);
        else
            arf_nan(res);
        return GR_SUCCESS;
    }
    else
    {
        ARF_FUNC_VIA_ARB(res, arb_exp, x)
    }
}

int
_gr_arf_log(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_neg_inf(res);
        else if (arf_is_pos_inf(x))
            arf_pos_inf(res);
        else
            arf_nan(res);
        return GR_SUCCESS;
    }
    else
    {
        ARF_FUNC_VIA_ARB(res, arb_log, x)
    }
}

int
_gr_arf_sin(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_zero(res);
        else
            arf_nan(res);
        return GR_SUCCESS;
    }
    else
    {
        ARF_FUNC_VIA_ARB(res, arb_sin, x)
    }
}

int
_gr_arf_cos(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_one(res);
        else
            arf_nan(res);
        return GR_SUCCESS;
    }
    else
    {
        ARF_FUNC_VIA_ARB(res, arb_cos, x)
    }
}

int
_gr_arf_tan(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_zero(res);
        else
            arf_nan(res);
        return GR_SUCCESS;
    }
    else
    {
        ARF_FUNC_VIA_ARB(res, arb_tan, x)
    }
}

int
_gr_arf_atan(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            arf_zero(res);
        }
        else if (arf_is_pos_inf(x))
        {
            _gr_arf_pi(res, ctx);
            arf_mul_2exp_si(res, res, -1);
        }
        else if (arf_is_neg_inf(x))
        {
            _gr_arf_pi(res, ctx);
            arf_mul_2exp_si(res, res, -1);
            arf_neg(res, res);
        }
        else
        {
            arf_nan(res);
        }
        return GR_SUCCESS;
    }
    else
    {
        ARF_FUNC_VIA_ARB(res, arb_atan, x)
    }
}

int
_gr_arf_sinh(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_zero(res);
        else if (arf_is_inf(x))
            arf_set(res, x);
        else
            arf_nan(res);
        return GR_SUCCESS;
    }
    else
    {
        ARF_FUNC_VIA_ARB(res, arb_sinh, x)
    }
}

int
_gr_arf_cosh(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_one(res);
        else if (arf_is_inf(x))
            arf_pos_inf(res);
        else
            arf_nan(res);
        return GR_SUCCESS;
    }
    else
    {
        ARF_FUNC_VIA_ARB(res, arb_cosh, x)
    }
}

int
_gr_arf_tanh(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_zero(res);
        else if (arf_is_inf(x))
            arf_set_si(res, arf_sgn(x));
        else
            arf_nan(res);
        return GR_SUCCESS;
    }
    else
    {
        ARF_FUNC_VIA_ARB(res, arb_tanh, x)
    }
}

/* todo: configurable function to return pole */

int
_gr_arf_gamma(arf_t res, const arf_t x, const gr_ctx_t ctx)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
            arf_pos_inf(res);
        else if (arf_is_inf(x))
            arf_pos_inf(res);
        else
            arf_nan(res);
        return GR_SUCCESS;
    }
    else if (arf_sgn(x) < 0 && arf_is_int(x))
    {
        arf_pos_inf(res);
        return GR_SUCCESS;
    }
    else
    {
        ARF_FUNC_VIA_ARB(res, arb_gamma, x)
    }
}


int
_gr_arf_vec_dot(arf_t res, const arf_t initial, int subtract, arf_srcptr vec1, arf_srcptr vec2, slong len, gr_ctx_t ctx)
{
    arf_approx_dot(res, initial, subtract, vec1, 1, vec2, 1, len, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_arf_vec_dot_rev(arf_t res, const arf_t initial, int subtract, arf_srcptr vec1, arf_srcptr vec2, slong len, gr_ctx_t ctx)
{
    arf_approx_dot(res, initial, subtract, vec1, 1, vec2 + len - 1, -1, len, ARF_CTX_PREC(ctx), ARF_CTX_RND(ctx));
    return GR_SUCCESS;
}

/*
int
_gr_arf_poly_mullow(arf_ptr res,
    arf_srcptr poly1, slong len1,
    arf_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    _arf_poly_mullow(res, poly1, len1, poly2, len2, n, ARF_CTX_PREC(ctx));
    return GR_SUCCESS;
}
*/

#include "gr_mat.h"
#include "arb_mat.h"

int
_gr_arf_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong prec;
    slong cutoff;

    prec = ARF_CTX_PREC(ctx);

    /* todo: detect small-integer matrices */
    if (prec <= 2 * FLINT_BITS)
        cutoff = 120;
    else if (prec <= 16 * FLINT_BITS)
        cutoff = 60;
    else
        cutoff = 40;

    if (A->r <= cutoff || A->c <= cutoff || B->c <= cutoff)
    {
        return gr_mat_mul_classical(C, A, B, ctx);
    }
    else
    {
        /* todo: direct algorithm, avoiding copies */
        arb_mat_t RC, RB, RA;
        slong i, j;
        arf_t zero;

        arb_mat_init(RA, A->r, A->c);
        arb_mat_init(RB, B->r, B->c);
        arb_mat_init(RC, C->r, C->c);

        arf_init(zero);

        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                *arb_midref(arb_mat_entry(RA, i, j)) = ((arf_srcptr) A->rows[i])[j];

        for (i = 0; i < B->r; i++)
            for (j = 0; j < B->c; j++)
                *arb_midref(arb_mat_entry(RB, i, j)) = ((arf_srcptr) B->rows[i])[j];

        arb_mat_approx_mul(RC, RA, RB, prec);

        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                *arb_midref(arb_mat_entry(RA, i, j)) = *zero;

        for (i = 0; i < B->r; i++)
            for (j = 0; j < B->c; j++)
                *arb_midref(arb_mat_entry(RB, i, j)) = *zero;

        for (i = 0; i < C->r; i++)
            for (j = 0; j < C->c; j++)
                arf_swap(((arf_ptr) C->rows[i]) + j, arb_midref(arb_mat_entry(RC, i, j)));

        arb_mat_clear(RA);
        arb_mat_clear(RB);
        arb_mat_clear(RC);

        return GR_SUCCESS;
    }
}

int
_gr_arf_ctx_clear(gr_ctx_t ctx)
{
    flint_free(ctx->elem_ctx);
    return GR_SUCCESS;
}

int _arf_methods_initialized = 0;

gr_static_method_table _arf_methods;

gr_method_tab_input _arf_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_arf_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_arf_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_arf_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_arf_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_arf_swap},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_arf_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_arf_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_arf_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_arf_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_arf_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_arf_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_arf_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_arf_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_arf_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_arf_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_arf_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_arf_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_arf_set_fmpq},
    {GR_METHOD_SET_D,           (gr_funcptr) _gr_arf_set_d},
    {GR_METHOD_SET_STR,         (gr_funcptr) _gr_arf_set_str},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_arf_set_other},

    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_arf_get_fmpz},
    {GR_METHOD_GET_UI,          (gr_funcptr) _gr_arf_get_ui},
    {GR_METHOD_GET_SI,          (gr_funcptr) _gr_arf_get_si},
    {GR_METHOD_GET_D,           (gr_funcptr) _gr_arf_get_d},

    {GR_METHOD_NEG,             (gr_funcptr) _gr_arf_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_arf_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_arf_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_arf_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_arf_add_fmpz},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_arf_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_arf_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_arf_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_arf_sub_fmpz},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_arf_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_arf_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_arf_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_arf_mul_fmpz},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_arf_mul_two},
    {GR_METHOD_ADDMUL,          (gr_funcptr) _gr_arf_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) _gr_arf_submul},
    {GR_METHOD_SQR,             (gr_funcptr) _gr_arf_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_arf_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_arf_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_arf_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_arf_div_fmpz},
    {GR_METHOD_INV,             (gr_funcptr) _gr_arf_inv},

    {GR_METHOD_POW,             (gr_funcptr) _gr_arf_pow},
/*
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_arf_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_arf_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_arf_pow_fmpz},
    {GR_METHOD_POW_FMPQ,        (gr_funcptr) _gr_arf_pow_fmpq},
*/
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_arf_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_arf_rsqrt},

    {GR_METHOD_FLOOR,           (gr_funcptr) _gr_arf_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) _gr_arf_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) _gr_arf_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) _gr_arf_nint},

    {GR_METHOD_ABS,             (gr_funcptr) _gr_arf_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) _gr_arf_conj},
    {GR_METHOD_RE,              (gr_funcptr) _gr_arf_set},
    {GR_METHOD_IM,              (gr_funcptr) _gr_arf_im},
    {GR_METHOD_SGN,             (gr_funcptr) _gr_arf_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) _gr_arf_sgn},
    {GR_METHOD_CMP,             (gr_funcptr) _gr_arf_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) _gr_arf_cmpabs},

    {GR_METHOD_I,               (gr_funcptr) gr_not_in_domain},
    {GR_METHOD_PI,              (gr_funcptr) _gr_arf_pi},
    {GR_METHOD_EXP,             (gr_funcptr) _gr_arf_exp},
    {GR_METHOD_LOG,             (gr_funcptr) _gr_arf_log},
    {GR_METHOD_SIN,             (gr_funcptr) _gr_arf_sin},
    {GR_METHOD_COS,             (gr_funcptr) _gr_arf_cos},
    {GR_METHOD_TAN,             (gr_funcptr) _gr_arf_tan},

    {GR_METHOD_SINH,            (gr_funcptr) _gr_arf_sinh},
    {GR_METHOD_COSH,            (gr_funcptr) _gr_arf_cosh},
    {GR_METHOD_TANH,            (gr_funcptr) _gr_arf_tanh},

    {GR_METHOD_ATAN,            (gr_funcptr) _gr_arf_atan},

    {GR_METHOD_GAMMA,            (gr_funcptr) _gr_arf_gamma},


    {GR_METHOD_VEC_DOT,         (gr_funcptr) _gr_arf_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _gr_arf_vec_dot_rev},
/*
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_arf_poly_mullow},
*/
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_arf_mat_mul},
    {GR_METHOD_MAT_DET,         (gr_funcptr) gr_mat_det_generic_field},
    {GR_METHOD_MAT_FIND_NONZERO_PIVOT,     (gr_funcptr) gr_mat_find_nonzero_pivot_large_abs},

    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_real_float_arf(gr_ctx_t ctx, slong prec)
{
    ctx->flags = 0;
    ctx->which_ring = GR_CTX_REAL_FLOAT_ARF;
    ctx->sizeof_elem = sizeof(arf_struct);
    ctx->elem_ctx = flint_malloc(sizeof(gr_arf_ctx));
    ctx->size_limit = WORD_MAX;

    gr_ctx_arf_set_prec(ctx, prec);
    ARF_CTX_RND(ctx) = ARF_RND_NEAR;

    ctx->methods = _arf_methods;

    if (!_arf_methods_initialized)
    {
        gr_method_tab_init(_arf_methods, _arf_methods_input);
        _arf_methods_initialized = 1;
    }
}
