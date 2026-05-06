#include "gr_ore_poly.h"
#include "gr_poly.h"

int
gr_ore_poly_ctx_over_gr_poly_base_ptrs(gr_ctx_struct ** Scalars,
                                    gr_ctx_struct ** Pol, const gr_ctx_t Ore)
{
    gr_ctx_struct * _Pol = GR_ORE_POLY_ELEM_CTX(Ore);

    ulong which_base = _Pol->which_ring;
    if (which_base != GR_CTX_GR_POLY)
        return GR_UNABLE;

    if (Pol != NULL)
        * Pol = _Pol;
    if (Scalars != NULL)
        * Scalars = POLYNOMIAL_ELEM_CTX(_Pol);

    return GR_SUCCESS;
}

int
_gr_ore_poly_indicial_polynomial_euler_derivative(gr_ptr ind, gr_srcptr op,
                                                  slong len, gr_ctx_t Ore)
{
    gr_ctx_struct * Scalars, * Pol;

    int status = GR_SUCCESS;

    GR_MUST_SUCCEED (gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Scalars, &Pol, Ore));

    for (slong i = 0; i < len; i++)
        status = gr_poly_get_coeff_scalar(
                GR_ENTRY(ind, i, Scalars->sizeof_elem),
                GR_ENTRY(op, i, Pol->sizeof_elem),
                0, Scalars);

    return status;
}

int
_gr_ore_poly_indicial_polynomial(gr_ptr ind, gr_srcptr op, slong len, gr_ctx_t Ore)
{
    switch (GR_ORE_POLY_CTX(Ore)->which_algebra)
    {
        case ORE_ALGEBRA_EULER_DERIVATIVE:
            return _gr_ore_poly_indicial_polynomial_euler_derivative(ind, op, len, Ore);
        default:
            return GR_UNABLE;
    }
}

int
gr_ore_poly_indicial_polynomial(gr_poly_t ind, const gr_ore_poly_t op, gr_ctx_t Ore)
{
    gr_ctx_struct * Scalars;

    if (gr_ore_poly_ctx_over_gr_poly_base_ptrs(&Scalars, NULL, Ore) != GR_SUCCESS)
        return GR_UNABLE;

    // XXX move to level depending on operator type?
    gr_poly_fit_length(ind, op->length, Scalars);

    int status = _gr_ore_poly_indicial_polynomial(ind->coeffs, op->coeffs, op->length, Ore);

    _gr_poly_set_length(ind, op->length, Scalars);
    _gr_poly_normalise(ind, Scalars);

    return status;
}
