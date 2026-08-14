#include "gr_ore_poly.h"
#include "gr_poly.h"

int
_gr_ore_poly_indicial_polynomial_euler_derivative(gr_ptr ind, gr_srcptr op,
                                                  slong len, gr_ctx_t ctx)
{
    gr_ctx_struct * Cst, * Pol;

    int status = GR_SUCCESS;

    GR_MUST_SUCCEED (gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Cst, &Pol, ctx));

    for (slong i = 0; i < len; i++)
        status |= gr_poly_get_coeff_scalar(
                GR_ENTRY(ind, i, Cst->sizeof_elem),
                GR_ENTRY(op, i, Pol->sizeof_elem),
                0, Cst);

    return status;
}

int
_gr_ore_poly_indicial_polynomial(gr_ptr ind,
                                 gr_srcptr op, slong len, gr_ctx_t ctx)
{
    switch (GR_ORE_POLY_CTX(ctx)->which_algebra)
    {
        case ORE_ALGEBRA_EULER_DERIVATIVE:
            return _gr_ore_poly_indicial_polynomial_euler_derivative(ind, op, len, ctx);
        default:
            return GR_UNABLE;
    }
}

int
gr_ore_poly_indicial_polynomial(gr_poly_t ind,
                                const gr_ore_poly_t op, gr_ctx_t ctx)
{
    gr_ctx_struct * Cst;

    if (gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Cst, NULL, ctx) != GR_SUCCESS)
        return GR_UNABLE;

    gr_poly_fit_length(ind, op->length, Cst);

    int status = _gr_ore_poly_indicial_polynomial(ind->coeffs, op->coeffs,
                                                  op->length, ctx);

    _gr_poly_set_length(ind, op->length, Cst);
    _gr_poly_normalise(ind, Cst);

    if (gr_poly_is_zero(ind, Cst) != T_FALSE)
        return GR_UNABLE;

    return status;
}

