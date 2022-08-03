#include "acf.h"
#include "acb.h"
#include "gr.h"

typedef struct
{
    slong prec;
    arf_rnd_t rnd;
}
gr_acf_ctx;

#define ACF_CTX_PREC(ring_ctx) (((gr_acf_ctx *)((ring_ctx)))->prec)
#define ACF_CTX_RND(ring_ctx) (((gr_acf_ctx *)((ring_ctx)))->rnd)

void gr_ctx_acf_set_prec(gr_ctx_t ctx, slong prec)
{
    prec = FLINT_MAX(prec, 2);
    prec = FLINT_MIN(prec, WORD_MAX / 8);

    ACF_CTX_PREC(ctx) = prec;
}

slong gr_ctx_acf_get_prec(gr_ctx_t ctx)
{
    return ACF_CTX_PREC(ctx);
}

int
_gr_acf_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Complex floating-point numbers (acf, prec = ");
    gr_stream_write_si(out, ACF_CTX_PREC(ctx));
    gr_stream_write(out, ")");
    return GR_SUCCESS;
}

void
_gr_acf_init(acf_t x, const gr_ctx_t ctx)
{
    acf_init(x);
}

void
_gr_acf_clear(acf_t x, const gr_ctx_t ctx)
{
    acf_clear(x);
}

void
_gr_acf_swap(acf_t x, acf_t y, const gr_ctx_t ctx)
{
    acf_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

/* todo: limits */
int
_gr_acf_randtest(acf_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    arf_randtest(acf_realref(res), state, ACF_CTX_PREC(ctx), 10);
    arf_randtest(acf_imagref(res), state, ACF_CTX_PREC(ctx), 10);
    return GR_SUCCESS;
}

/* todo */
int
_gr_acf_write(gr_stream_t out, const acf_t x, const gr_ctx_t ctx)
{
    slong digits = ACF_CTX_PREC(ctx) * 0.30102999566398 + 1;

    if (arf_is_zero(acf_imagref(x)))
    {
        gr_stream_write_free(out, arf_get_str(acf_realref(x), digits));
    }
    else if (arf_is_zero(acf_realref(x)))
    {
        gr_stream_write_free(out, arf_get_str(acf_imagref(x), digits));
        gr_stream_write(out, "*I");
    }
    else
    {
        gr_stream_write(out, "(");
        gr_stream_write_free(out, arf_get_str(acf_realref(x), digits));

        if (arf_sgn(acf_imagref(x)) < 0)
        {
            arf_t t;
            arf_init_neg_shallow(t, acf_imagref(x));
            gr_stream_write(out, " - ");
            gr_stream_write_free(out, arf_get_str(t, digits));
        }
        else
        {
            gr_stream_write(out, " + ");
            gr_stream_write_free(out, arf_get_str(acf_imagref(x), digits));
        }

        gr_stream_write(out, "*I)");
    }

    return GR_SUCCESS;
}

int
_gr_acf_zero(acf_t x, const gr_ctx_t ctx)
{
    arf_zero(acf_realref(x));
    arf_zero(acf_imagref(x));
    return GR_SUCCESS;
}

int
_gr_acf_one(acf_t x, const gr_ctx_t ctx)
{
    arf_one(acf_realref(x));
    arf_zero(acf_imagref(x));
    return GR_SUCCESS;
}

int
_gr_acf_set_si(acf_t res, slong v, const gr_ctx_t ctx)
{
    arf_set_round_si(acf_realref(res), v, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_zero(acf_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acf_set_ui(acf_t res, ulong v, const gr_ctx_t ctx)
{
    arf_set_round_ui(acf_realref(res), v, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_zero(acf_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acf_set_fmpz(acf_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    arf_set_round_fmpz(acf_realref(res), v, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_zero(acf_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acf_set_fmpq(acf_t res, const fmpq_t v, const gr_ctx_t ctx)
{
    arf_set_fmpq(acf_realref(res), v, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_zero(acf_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acf_set_d(acf_t res, double x, const gr_ctx_t ctx)
{
    arf_set_d(acf_realref(res), x);
    arf_zero(acf_imagref(res));
    return GR_SUCCESS;
}

/* todo: set_round? */
int
_gr_acf_set(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
    acf_set(res, x);
    return GR_SUCCESS;
}

/* todo
int
_gr_acf_set_str(acf_t res, const char * x, const gr_ctx_t ctx)
{
}
*/

int
_gr_acf_set_other(acf_t res, gr_srcptr x, gr_ctx_t x_ctx, const gr_ctx_t ctx)
{
    switch (x_ctx->which_ring)
    {
        case GR_CTX_FMPZ:
            return _gr_acf_set_fmpz(res, x, ctx);

        case GR_CTX_FMPQ:
            return _gr_acf_set_fmpq(res, x, ctx);

        case GR_CTX_REAL_FLOAT_ARF:
            arf_set_round(acf_realref(res), x, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
            arf_zero(acf_imagref(res));
            return GR_SUCCESS;

        case GR_CTX_COMPLEX_FLOAT_ACF:
            return _gr_acf_set(res, x, ctx);

        case GR_CTX_RR_ARB:
            arf_set_round(acf_realref(res), arb_midref((arb_srcptr) x), ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
            arf_zero(acf_imagref(res));
            return GR_SUCCESS;

        case GR_CTX_CC_ACB:
            arf_set_round(acf_realref(res), arb_midref(acb_realref((acb_srcptr) x)), ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
            arf_set_round(acf_imagref(res), arb_midref(acb_imagref((acb_srcptr) x)), ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
            return GR_SUCCESS;

        default:
            {
                gr_ctx_t cctx;
                acb_t z;
                int status;

                gr_ctx_init_complex_acb(cctx, ACF_CTX_PREC(ctx) + 20);
                acb_init(z);

                status = gr_set_other(z, x, x_ctx, cctx);

                if (status == GR_SUCCESS)
                {
                    arf_set_round(acf_realref(res), arb_midref(acb_realref(z)), ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
                    arf_set_round(acf_imagref(res), arb_midref(acb_imagref(z)), ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
                }

                acb_clear(z);
                gr_ctx_clear(cctx);

                return status;
            }
    }
}

int
_gr_acf_get_fmpz(fmpz_t res, const acf_t x, const gr_ctx_t ctx)
{
    if (!arf_is_zero(acf_imagref(x)))
        return GR_DOMAIN;

    if (!arf_is_int(acf_realref(x)))
        return GR_DOMAIN;

    /* todo: detect mpz overflow */
    if (arf_cmpabs_2exp_si(acf_realref(x), WORD_MAX) >= 0)
        return GR_UNABLE;

    arf_get_fmpz(res, acf_realref(x), ARF_RND_DOWN);
    return GR_SUCCESS;
}

int
_gr_acf_get_si(slong * res, const acf_t x, const gr_ctx_t ctx)
{
    fmpz_t t;

    if (!arf_is_zero(acf_imagref(x)))
        return GR_DOMAIN;

    if (!arf_is_int(acf_realref(x)))
        return GR_DOMAIN;

    if (arf_cmp_si(acf_realref(x), WORD_MIN) < 0 || arf_cmp_si(acf_realref(x), WORD_MAX) > 0)
        return GR_DOMAIN;
    fmpz_init(t);
    arf_get_fmpz(t, acf_realref(x), ARF_RND_DOWN);
    *res = fmpz_get_si(t);
    fmpz_clear(t);

    return GR_SUCCESS;
}

int
_gr_acf_get_ui(ulong * res, const acf_t x, const gr_ctx_t ctx)
{
    fmpz_t t;

    if (!arf_is_zero(acf_imagref(x)))
        return GR_DOMAIN;

    if (!arf_is_int(acf_realref(x)))
        return GR_DOMAIN;

    if (arf_sgn(acf_realref(x)) < 0 || arf_cmp_ui(acf_realref(x), UWORD_MAX) > 0)
        return GR_DOMAIN;

    fmpz_init(t);
    arf_get_fmpz(t, acf_realref(x), ARF_RND_DOWN);
    *res = fmpz_get_ui(t);
    fmpz_clear(t);

    return GR_SUCCESS;
}

int
_gr_acf_get_d(double * res, const acf_t x, const gr_ctx_t ctx)
{
    if (!arf_is_zero(acf_imagref(x)))
        return GR_DOMAIN;

    *res = arf_get_d(acf_realref(x), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_acf_is_zero(const acf_t x, const gr_ctx_t ctx)
{
    return (arf_is_zero(acf_realref(x)) && arf_is_zero(acf_imagref(x))) ? T_TRUE : T_FALSE;
}

truth_t
_gr_acf_is_one(const acf_t x, const gr_ctx_t ctx)
{
    return (arf_is_one(acf_realref(x)) && arf_is_zero(acf_imagref(x))) ? T_TRUE : T_FALSE;
}

truth_t
_gr_acf_is_neg_one(const acf_t x, const gr_ctx_t ctx)
{
    return (arf_equal_si(acf_realref(x), -1) && arf_is_zero(acf_imagref(x))) ? T_TRUE : T_FALSE;
}

truth_t
_gr_acf_equal(const acf_t x, const acf_t y, const gr_ctx_t ctx)
{
    if (arf_is_nan(acf_realref(x)) || arf_is_nan(acf_imagref(x)) ||
        arf_is_nan(acf_realref(y)) || arf_is_nan(acf_imagref(y)))
        return T_FALSE;

    return acf_equal(x, y) ? T_TRUE : T_FALSE;
}

/* todo: neg_round? */
int
_gr_acf_neg(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
    arf_neg(acf_realref(res), acf_realref(x));
    arf_neg(acf_imagref(res), acf_imagref(x));
    return GR_SUCCESS;
}

int
_gr_acf_add(acf_t res, const acf_t x, const acf_t y, const gr_ctx_t ctx)
{
    acf_add(res, x, y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_add_si(acf_t res, const acf_t x, slong y, const gr_ctx_t ctx)
{
    arf_add_si(acf_realref(res), acf_realref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_set_round(acf_realref(res), acf_realref(x), ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_add_ui(acf_t res, const acf_t x, ulong y, const gr_ctx_t ctx)
{
    arf_add_ui(acf_realref(res), acf_realref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_set_round(acf_realref(res), acf_realref(x), ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_add_fmpz(acf_t res, const acf_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arf_add_fmpz(acf_realref(res), acf_realref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_set_round(acf_realref(res), acf_realref(x), ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_sub(acf_t res, const acf_t x, const acf_t y, const gr_ctx_t ctx)
{
    acf_sub(res, x, y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_sub_si(acf_t res, const acf_t x, slong y, const gr_ctx_t ctx)
{
    arf_sub_si(acf_realref(res), acf_realref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_set_round(acf_realref(res), acf_realref(x), ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_sub_ui(acf_t res, const acf_t x, ulong y, const gr_ctx_t ctx)
{
    arf_sub_ui(acf_realref(res), acf_realref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_set_round(acf_realref(res), acf_realref(x), ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_sub_fmpz(acf_t res, const acf_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arf_sub_fmpz(acf_realref(res), acf_realref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_set_round(acf_realref(res), acf_realref(x), ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_mul(acf_t res, const acf_t x, const acf_t y, const gr_ctx_t ctx)
{
    acf_mul(res, x, y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_mul_si(acf_t res, const acf_t x, slong y, const gr_ctx_t ctx)
{
    arf_mul_si(acf_realref(res), acf_realref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_mul_si(acf_imagref(res), acf_imagref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_mul_ui(acf_t res, const acf_t x, ulong y, const gr_ctx_t ctx)
{
    arf_mul_ui(acf_realref(res), acf_realref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_mul_ui(acf_imagref(res), acf_imagref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_mul_fmpz(acf_t res, const acf_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arf_mul_fmpz(acf_realref(res), acf_realref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_mul_fmpz(acf_imagref(res), acf_imagref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_mul_two(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
    arf_mul_2exp_si(acf_realref(res), acf_realref(x), 1);
    arf_mul_2exp_si(acf_imagref(res), acf_imagref(x), 1);
    return GR_SUCCESS;
}

int
_gr_acf_sqr(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
    acf_mul(res, x, x, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_inv(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
    acf_approx_inv(res, x, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_div(acf_t res, const acf_t x, const acf_t y, const gr_ctx_t ctx)
{
    acf_approx_div(res, x, y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_div_si(acf_t res, const acf_t x, slong y, const gr_ctx_t ctx)
{
    arf_div_si(acf_realref(res), acf_realref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_div_si(acf_imagref(res), acf_imagref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_div_ui(acf_t res, const acf_t x, ulong y, const gr_ctx_t ctx)
{
    arf_div_ui(acf_realref(res), acf_realref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_div_ui(acf_imagref(res), acf_imagref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_div_fmpz(acf_t res, const acf_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    arf_div_fmpz(acf_realref(res), acf_realref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_div_fmpz(acf_imagref(res), acf_imagref(x), y, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_sqrt(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
    acf_approx_sqrt(res, x, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_abs(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
    if (arf_is_zero(acf_imagref(x)))
    {
        arf_abs(acf_realref(res), acf_realref(x));
    }
    else if (arf_is_zero(acf_realref(x)))
    {
        arf_abs(acf_realref(res), acf_imagref(x));
    }
    else
    {
        arf_sosq(acf_realref(res), acf_realref(x), acf_imagref(x), ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
        arf_sqrt(acf_realref(res), acf_realref(res), ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    }

    arf_zero(acf_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acf_conj(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
    arf_set(acf_realref(res), acf_realref(x));
    arf_neg(acf_imagref(res), acf_imagref(x));
    return GR_SUCCESS;
}

int
_gr_acf_re(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
    arf_set(acf_realref(res), acf_realref(x));
    arf_zero(acf_imagref(res));
    return GR_SUCCESS;
}


int
_gr_acf_im(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
    arf_set(acf_realref(res), acf_imagref(x));
    arf_zero(acf_imagref(res));
    return GR_SUCCESS;
}

/*
int
_gr_acf_sgn(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
    return GR_SUCCESS;
}

int
_gr_acf_csgn(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
}

int
_gr_acf_rsqrt(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
}

int
_gr_acf_floor(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
}

int
_gr_acf_ceil(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
}

int
_gr_acf_trunc(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
}

int
_gr_acf_nint(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
}
*/

/* todo: handling nan */
int
_gr_acf_cmp(int * res, const acf_t x, const acf_t y, const gr_ctx_t ctx)
{
    if (!arf_is_zero(acf_imagref(x)) || !arf_is_zero(acf_imagref(y)))
        return GR_DOMAIN;

    *res = arf_cmp(acf_realref(x), acf_realref(y));
    return GR_SUCCESS;
}

int
_gr_acf_cmpabs(int * res, const acf_t x, const acf_t y, const gr_ctx_t ctx)
{
    if (!arf_is_zero(acf_imagref(x)) || !arf_is_zero(acf_imagref(y)))
        return GR_UNABLE;

    *res = arf_cmp(acf_realref(x), acf_realref(y));
    return GR_SUCCESS;
}

#define ACF_FUNC_VIA_ACB(res, acb_func, x) \
    acb_t r, t; \
    slong prec, wp, extra; \
    int status = GR_SUCCESS; \
    prec = ACF_CTX_PREC(ctx); \
    acb_init(r); \
    *arb_midref(acb_realref(t)) = *acf_realref(x); \
    *arb_midref(acb_imagref(t)) = *acf_imagref(x); \
    mag_init(arb_radref(acb_realref(t))); \
    mag_init(arb_radref(acb_imagref(t))); \
    for (extra = 10 + prec * 0.01; ; extra += FLINT_MAX(extra, 32)) \
    { \
        wp = prec + extra; \
        if (wp > 10 * prec + 1000) \
        { \
            status = GR_UNABLE; \
            arf_nan(acf_realref(res)); \
            arf_nan(acf_imagref(res)); \
            break; \
        } \
        acb_func(r, t, wp); \
        if (acb_rel_accuracy_bits(r) >= prec) \
        { \
            arf_set_round(acf_realref(res), arb_midref(acb_realref(r)), prec, ACF_CTX_RND(ctx)); \
            arf_set_round(acf_imagref(res), arb_midref(acb_imagref(r)), prec, ACF_CTX_RND(ctx)); \
            break; \
        } \
    } \
    acb_clear(r); \
    return status; \

#define ACF_FUNC2_VIA_ACB(res, acb_func, x, y) \
    acb_t r, t, u; \
    slong prec, wp, extra; \
    int status = GR_SUCCESS; \
    prec = ACF_CTX_PREC(ctx); \
    acb_init(r); \
    *arb_midref(acb_realref(t)) = *acf_realref(x); \
    *arb_midref(acb_imagref(t)) = *acf_imagref(x); \
    mag_init(arb_radref(acb_realref(t))); \
    mag_init(arb_radref(acb_imagref(t))); \
    *arb_midref(acb_realref(u)) = *acf_realref(y); \
    *arb_midref(acb_imagref(u)) = *acf_imagref(y); \
    mag_init(arb_radref(acb_realref(u))); \
    mag_init(arb_radref(acb_imagref(u))); \
    for (extra = 10 + prec * 0.01; ; extra += FLINT_MAX(extra, 32)) \
    { \
        wp = prec + extra; \
        if (wp > 10 * prec + 1000) \
        { \
            status = GR_UNABLE; \
            arf_nan(acf_realref(res)); \
            arf_nan(acf_imagref(res)); \
            break; \
        } \
        acb_func(r, t, u, wp); \
        if (acb_rel_accuracy_bits(r) >= prec) \
        { \
            arf_set_round(acf_realref(res), arb_midref(acb_realref(r)), prec, ACF_CTX_RND(ctx)); \
            arf_set_round(acf_imagref(res), arb_midref(acb_imagref(r)), prec, ACF_CTX_RND(ctx)); \
            break; \
        } \
    } \
    acb_clear(r); \
    return status; \

/* todo: lots of special cases */
int
_gr_acf_pow(acf_t res, const acf_t x, const acf_t y, const gr_ctx_t ctx)
{
    ACF_FUNC2_VIA_ACB(res, acb_pow, x, y)
}

int
_gr_acf_i(acf_t res, const gr_ctx_t ctx)
{
    arf_zero(acf_realref(res));
    arf_one(acf_imagref(res));
    return GR_SUCCESS;
}

int
_gr_acf_pi(acf_t res, const gr_ctx_t ctx)
{
    arb_t t;
    arb_init(t);
    arb_const_pi(t, ACF_CTX_PREC(ctx) + 30);
    arf_set_round(acf_realref(res), arb_midref(t), ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    arf_zero(acf_imagref(res));
    arb_clear(t);
    return GR_SUCCESS;
}

int
_gr_acf_exp(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
    ACF_FUNC_VIA_ACB(res, acb_exp, x)
}

int
_gr_acf_log(acf_t res, const acf_t x, const gr_ctx_t ctx)
{
    ACF_FUNC_VIA_ACB(res, acb_log, x)
}

int
_gr_acf_vec_dot(acf_t res, const acf_t initial, int subtract, acf_srcptr vec1, acf_srcptr vec2, slong len, gr_ctx_t ctx)
{
    acf_approx_dot(res, initial, subtract, vec1, 1, vec2, 1, len, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

int
_gr_acf_vec_dot_rev(acf_t res, const acf_t initial, int subtract, acf_srcptr vec1, acf_srcptr vec2, slong len, gr_ctx_t ctx)
{
    acf_approx_dot(res, initial, subtract, vec1, 1, vec2 + len - 1, -1, len, ACF_CTX_PREC(ctx), ACF_CTX_RND(ctx));
    return GR_SUCCESS;
}

/*
int
_gr_acf_poly_mullow(acf_ptr res,
    acf_srcptr poly1, slong len1,
    acf_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    _acf_poly_mullow(res, poly1, len1, poly2, len2, n, ACF_CTX_PREC(ctx));
    return GR_SUCCESS;
}
*/

#include "gr_mat.h"
#include "acb_mat.h"

int
_gr_acf_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong prec;
    slong cutoff;

    prec = ACF_CTX_PREC(ctx);

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
        acb_mat_t RC, RB, RA;
        slong i, j;
        acf_t zero;

        acb_mat_init(RA, A->r, A->c);
        acb_mat_init(RB, B->r, B->c);
        acb_mat_init(RC, C->r, C->c);

        acf_init(zero);

        for (i = 0; i < A->r; i++)
        {
            for (j = 0; j < A->c; j++)
            {
                *arb_midref(acb_realref(acb_mat_entry(RA, i, j))) = *acf_realref(((acf_srcptr) A->rows[i]) + j);
                *arb_midref(acb_imagref(acb_mat_entry(RA, i, j))) = *acf_imagref(((acf_srcptr) A->rows[i]) + j);
            }
        }

        for (i = 0; i < B->r; i++)
        {
            for (j = 0; j < B->c; j++)
            {
                *arb_midref(acb_realref(acb_mat_entry(RB, i, j))) = *acf_realref(((acf_srcptr) B->rows[i]) + j);
                *arb_midref(acb_imagref(acb_mat_entry(RB, i, j))) = *acf_imagref(((acf_srcptr) B->rows[i]) + j);
            }
        }

        acb_mat_approx_mul(RC, RA, RB, prec);

        for (i = 0; i < A->r; i++)
        {
            for (j = 0; j < A->c; j++)
            {
                *arb_midref(acb_realref(acb_mat_entry(RA, i, j))) = *acf_realref(zero);
                *arb_midref(acb_imagref(acb_mat_entry(RA, i, j))) = *acf_imagref(zero);
            }
        }

        for (i = 0; i < B->r; i++)
        {
            for (j = 0; j < B->c; j++)
            {
                *arb_midref(acb_realref(acb_mat_entry(RB, i, j))) = *acf_realref(zero);
                *arb_midref(acb_imagref(acb_mat_entry(RB, i, j))) = *acf_imagref(zero);
            }
        }

        for (i = 0; i < C->r; i++)
        {
            for (j = 0; j < C->c; j++)
            {
                arf_swap(acf_realref(((acf_ptr) C->rows[i]) + j), arb_midref(acb_realref(acb_mat_entry(RC, i, j))));
                arf_swap(acf_imagref(((acf_ptr) C->rows[i]) + j), arb_midref(acb_imagref(acb_mat_entry(RC, i, j))));
            }
        }

        acb_mat_clear(RA);
        acb_mat_clear(RB);
        acb_mat_clear(RC);

        return GR_SUCCESS;
    }
}


int _acf_methods_initialized = 0;

gr_static_method_table _acf_methods;

gr_method_tab_input _acf_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_acf_ctx_write},
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
    {GR_METHOD_INIT,            (gr_funcptr) _gr_acf_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_acf_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_acf_swap},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_acf_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_acf_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_acf_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_acf_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_acf_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_acf_is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _gr_acf_is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_acf_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_acf_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_acf_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_acf_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_acf_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _gr_acf_set_fmpq},
    {GR_METHOD_SET_D,           (gr_funcptr) _gr_acf_set_d},
/*    {GR_METHOD_SET_STR,         (gr_funcptr) _gr_acf_set_str}, */
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _gr_acf_set_other},

    {GR_METHOD_GET_FMPZ,        (gr_funcptr) _gr_acf_get_fmpz},
    {GR_METHOD_GET_UI,          (gr_funcptr) _gr_acf_get_ui},
    {GR_METHOD_GET_SI,          (gr_funcptr) _gr_acf_get_si},
    {GR_METHOD_GET_D,           (gr_funcptr) _gr_acf_get_d},

    {GR_METHOD_NEG,             (gr_funcptr) _gr_acf_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_acf_add},
    {GR_METHOD_ADD_UI,          (gr_funcptr) _gr_acf_add_ui},
    {GR_METHOD_ADD_SI,          (gr_funcptr) _gr_acf_add_si},
    {GR_METHOD_ADD_FMPZ,        (gr_funcptr) _gr_acf_add_fmpz},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_acf_sub},
    {GR_METHOD_SUB_UI,          (gr_funcptr) _gr_acf_sub_ui},
    {GR_METHOD_SUB_SI,          (gr_funcptr) _gr_acf_sub_si},
    {GR_METHOD_SUB_FMPZ,        (gr_funcptr) _gr_acf_sub_fmpz},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_acf_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_acf_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_acf_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_acf_mul_fmpz},
    {GR_METHOD_MUL_TWO,         (gr_funcptr) _gr_acf_mul_two},
/*
    {GR_METHOD_ADDMUL,          (gr_funcptr) _gr_acf_addmul},
    {GR_METHOD_SUBMUL,          (gr_funcptr) _gr_acf_submul},
*/
    {GR_METHOD_SQR,             (gr_funcptr) _gr_acf_sqr},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_acf_div},
    {GR_METHOD_DIV_UI,          (gr_funcptr) _gr_acf_div_ui},
    {GR_METHOD_DIV_SI,          (gr_funcptr) _gr_acf_div_si},
    {GR_METHOD_DIV_FMPZ,        (gr_funcptr) _gr_acf_div_fmpz},
    {GR_METHOD_INV,             (gr_funcptr) _gr_acf_inv},

    {GR_METHOD_POW,             (gr_funcptr) _gr_acf_pow},
/*
    {GR_METHOD_POW_UI,          (gr_funcptr) _gr_acf_pow_ui},
    {GR_METHOD_POW_SI,          (gr_funcptr) _gr_acf_pow_si},
    {GR_METHOD_POW_FMPZ,        (gr_funcptr) _gr_acf_pow_fmpz},
    {GR_METHOD_POW_FMPQ,        (gr_funcptr) _gr_acf_pow_fmpq},
*/
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_acf_sqrt},
/*
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_acf_rsqrt},

    {GR_METHOD_FLOOR,           (gr_funcptr) _gr_acf_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) _gr_acf_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) _gr_acf_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) _gr_acf_nint},
*/

    {GR_METHOD_ABS,             (gr_funcptr) _gr_acf_abs},
    {GR_METHOD_CONJ,            (gr_funcptr) _gr_acf_conj},
    {GR_METHOD_RE,              (gr_funcptr) _gr_acf_re},
    {GR_METHOD_IM,              (gr_funcptr) _gr_acf_im},
/*
    {GR_METHOD_SGN,             (gr_funcptr) _gr_acf_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) _gr_acf_sgn},
*/
    {GR_METHOD_CMP,             (gr_funcptr) _gr_acf_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) _gr_acf_cmpabs},

    {GR_METHOD_I,               (gr_funcptr) _gr_acf_i},
    {GR_METHOD_PI,              (gr_funcptr) _gr_acf_pi},
    {GR_METHOD_EXP,             (gr_funcptr) _gr_acf_exp},
    {GR_METHOD_LOG,             (gr_funcptr) _gr_acf_log},

    {GR_METHOD_VEC_DOT,         (gr_funcptr) _gr_acf_vec_dot},
    {GR_METHOD_VEC_DOT_REV,     (gr_funcptr) _gr_acf_vec_dot_rev},
/*
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_acf_poly_mullow},
*/
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_acf_mat_mul},
    {GR_METHOD_MAT_DET,         (gr_funcptr) gr_mat_det_generic_field},
    /* todo */
    {GR_METHOD_MAT_FIND_NONZERO_PIVOT,     (gr_funcptr) gr_mat_find_nonzero_pivot_large_abs},

    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_complex_float_acf(gr_ctx_t ctx, slong prec)
{
    ctx->which_ring = GR_CTX_COMPLEX_FLOAT_ACF;
    ctx->sizeof_elem = sizeof(acf_struct);
    ctx->size_limit = WORD_MAX;

    gr_ctx_acf_set_prec(ctx, prec);
    ACF_CTX_RND(ctx) = ARF_RND_NEAR;

    ctx->methods = _acf_methods;

    if (!_acf_methods_initialized)
    {
        gr_method_tab_init(_acf_methods, _acf_methods_input);
        _acf_methods_initialized = 1;
    }
}
