#include "gr_poly.h"
#include "gr_ore_poly.h"

/* gr_poly/compose.c */
int _gr_poly_inflate(gr_ptr poly, slong len, slong n, gr_ctx_t ctx);

int
sigma_delta_unable(gr_ptr sigma, gr_ptr delta,
                   gr_srcptr a,
                   gr_ore_poly_ctx_struct * ctx)
{
    return GR_UNABLE;
}

int
sigma_delta_commutative(gr_ptr sigma, gr_ptr delta,
                        gr_srcptr a,
                        gr_ore_poly_ctx_struct * ctx)
{
    gr_ctx_ptr cctx = GR_ORE_POLY_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    if (sigma != NULL && sigma != a)
        status |= gr_set(sigma, a, cctx);

    if (delta != NULL)
        status |= gr_zero(delta, cctx);

    return status;
}

int
sigma_delta_derivative(gr_ptr sigma, gr_ptr delta,
                       gr_srcptr a,
                       gr_ore_poly_ctx_struct * ctx)
{
    gr_ctx_struct * cctx = GR_ORE_POLY_ELEM_CTX(ctx);
    /* todo: gr_derivative taking a generator index? */
    if (cctx->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;

    int status = GR_SUCCESS;

    if (sigma != NULL && sigma != a)
        status |= gr_set(sigma, a, cctx);

    if (delta != NULL)
        status |= gr_poly_derivative(delta, a, POLYNOMIAL_ELEM_CTX(cctx));

    return status;
}

int
sigma_delta_euler_derivative(gr_ptr sigma, gr_ptr delta,
                             gr_srcptr a,
                             gr_ore_poly_ctx_struct * ctx)
{
    int status = GR_SUCCESS;
    gr_ctx_struct * cctx = GR_ORE_POLY_ELEM_CTX(ctx);
    if (cctx->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;

    /* - for polynomial base rings, could call _gr_poly_derivative and avoid the
     *   post-shift
     * - gr_euler_derivative for arbitrary rings??? */

    status |= sigma_delta_derivative(sigma, delta, a, ctx);

    if (status == GR_SUCCESS && delta != NULL)
        status |= gr_poly_shift_left(delta, delta, 1, POLYNOMIAL_ELEM_CTX(cctx));

    return status;
}

static int
sigma_delta_shift_si(gr_ptr sigma, gr_ptr delta,
                     gr_srcptr a, slong shift, slong difference,
                     gr_ore_poly_ctx_struct * ctx)
{
    gr_ctx_struct * cctx = GR_ORE_POLY_ELEM_CTX(ctx);
    if (cctx->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;
    gr_ctx_struct * sctx = POLYNOMIAL_ELEM_CTX(cctx);

    int status = GR_SUCCESS;

    gr_ptr _shift;
    gr_poly_t shifted;
    GR_TMP_INIT(_shift, sctx);
    gr_poly_init(shifted, sctx);

    status |= gr_set_si(_shift, shift, sctx);
    status |= gr_poly_taylor_shift(shifted, a, _shift, sctx);

    if (delta != NULL)
    {
        if (difference == 0)
            status |= gr_zero(delta, cctx);
        else if (difference == 1)
            status |= gr_poly_sub(delta, shifted, a, sctx);
        else if (difference == -1)
            status |= gr_poly_sub(delta, a, shifted, sctx);
    }

    if (sigma != NULL)
        gr_swap(sigma, shifted, cctx);

    gr_poly_clear(shifted, sctx);
    GR_TMP_CLEAR(_shift, sctx);

    return status;
}

int
sigma_delta_forward_shift(gr_ptr sigma, gr_ptr delta,
                          gr_srcptr a,
                          gr_ore_poly_ctx_struct * ctx)
{
    return sigma_delta_shift_si(sigma, delta, a, +1, 0, ctx);
}

int
sigma_delta_backward_shift(gr_ptr sigma, gr_ptr delta,
                           gr_srcptr a,
                           gr_ore_poly_ctx_struct * ctx)
{
    return sigma_delta_shift_si(sigma, delta, a, -1, 0, ctx);
}

int
sigma_delta_forward_difference(gr_ptr sigma, gr_ptr delta,
                               gr_srcptr a,
                               gr_ore_poly_ctx_struct * ctx)
{
    return sigma_delta_shift_si(sigma, delta, a, +1, +1, ctx);
}

int
sigma_delta_backward_difference(gr_ptr sigma, gr_ptr delta,
                                gr_srcptr a,
                                gr_ore_poly_ctx_struct * ctx)
{
    return sigma_delta_shift_si(sigma, delta, a, -1, -1, ctx);
}

int
sigma_delta_compose(gr_ptr sigma, gr_ptr delta,
                    gr_srcptr a,
                    gr_ore_poly_ctx_struct * ctx)
{
    gr_ctx_struct * cctx = GR_ORE_POLY_ELEM_CTX(ctx);
    if (cctx->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;
    gr_ctx_struct * sctx = POLYNOMIAL_ELEM_CTX(cctx);

    int status = GR_SUCCESS;

    if (sigma != NULL)
        status |= gr_poly_compose(sigma, a, GR_ORE_POLY_ORE_DATA(ctx)->sigma_x, sctx);

    if (delta != NULL)
        status |= gr_zero(delta, cctx);

    return status;
}

int
sigma_delta_mahler(gr_ptr sigma, gr_ptr delta,
                   gr_srcptr a,
                   gr_ore_poly_ctx_struct * ctx)
{
    gr_ctx_struct * cctx = GR_ORE_POLY_ELEM_CTX(ctx);
    if (cctx->which_ring != GR_CTX_GR_POLY)
        return GR_UNABLE;
    gr_ctx_struct * sctx = POLYNOMIAL_ELEM_CTX(cctx);

    int status = GR_SUCCESS;

    if (sigma != NULL)
    {
        gr_poly_struct * _sigma = sigma;
        slong b = GR_ORE_POLY_ORE_DATA(ctx)->mahler_base;
        status |= gr_set(sigma, a, cctx);
        if (_sigma->length > 1)
        {
            slong newlen = (_sigma->length - 1) * b + 1;
            gr_poly_fit_length(sigma, newlen, sctx);
            status |= _gr_poly_inflate(_sigma->coeffs, _sigma->length, b, sctx);
            gr_poly_fit_length(sigma, newlen, sctx);
            _gr_poly_set_length(sigma, newlen, sctx);
            _gr_poly_normalise(sigma, sctx);
        }
    }

    if (delta != NULL)
        status |= gr_zero(delta, cctx);

    return status;
}

int
sigma_delta_frobenius(gr_ptr sigma, gr_ptr delta,
                      gr_srcptr a,
                      gr_ore_poly_ctx_struct * ctx)
{
    gr_ctx_ptr cctx = GR_ORE_POLY_ELEM_CTX(ctx);
    int status = GR_SUCCESS;

    if (sigma != NULL)
        status |= gr_fq_frobenius(sigma, a, 1, cctx);

    if (delta != NULL)
        status |= gr_zero(delta, cctx);

    return status;
}

const gr_ore_poly_sigma_delta_t _gr_ore_poly_default_sigma_delta[] =
{
    [ORE_ALGEBRA_CUSTOM]              = sigma_delta_unable,
    [ORE_ALGEBRA_COMMUTATIVE]         = sigma_delta_commutative,
    [ORE_ALGEBRA_DERIVATIVE]          = sigma_delta_derivative,
    [ORE_ALGEBRA_EULER_DERIVATIVE]    = sigma_delta_euler_derivative,
    [ORE_ALGEBRA_FORWARD_SHIFT]       = sigma_delta_forward_shift,
    [ORE_ALGEBRA_FORWARD_DIFFERENCE]  = sigma_delta_forward_difference,
    [ORE_ALGEBRA_BACKWARD_SHIFT]      = sigma_delta_backward_shift,
    [ORE_ALGEBRA_BACKWARD_DIFFERENCE] = sigma_delta_backward_difference,
    [ORE_ALGEBRA_Q_SHIFT]             = sigma_delta_compose,
    [ORE_ALGEBRA_MAHLER]              = sigma_delta_mahler,
    [ORE_ALGEBRA_FROBENIUS]           = sigma_delta_frobenius,
    [ORE_POLY_NUM_ALGEBRAS]           = NULL,
};
