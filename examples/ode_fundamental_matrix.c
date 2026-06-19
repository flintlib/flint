#include "acb_types.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_ode.h"
#include "gr.h"
#include "gr_ore_poly.h"

int main(void)
{
    gr_ctx_t ZZ, Pol, Dop;
    gr_ptr dop;

    int status = GR_SUCCESS;

    slong prec = 50;

    /* Differential operator */

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ctx_init_gr_ore_poly(Dop, Pol, 0, ORE_ALGEBRA_EULER_DERIVATIVE);
    status |= gr_ctx_set_gen_name(Pol, "z");
    status |= gr_ctx_set_gen_name(Dop, "Tz");

    GR_TMP_INIT(dop, Dop);

    status |= gr_ore_poly_set_str(dop,  "(z^2 + 1)*Tz^2 + (z^2 - 1)*Tz", Dop);

    flint_printf("dop = %{gr}\n", dop, Dop);

    /* Evaluation point */

    acb_t pt;
    acb_init(pt);
    acb_const_pi(pt, prec + 4);
    arb_one(acb_imagref(pt));
    acb_inv(pt, pt, prec + 4);

    /* Output matrix */

    acb_mat_t mat;
    slong dop_order = gr_ore_poly_length(dop, Dop) - 1;
    acb_mat_init(mat, dop_order, dop_order);

    /* Numerical solution */

    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt, 0, prec);

    acb_mat_printd(mat, prec * 0.3);

    acb_clear(pt);
    acb_mat_clear(mat);
    GR_TMP_CLEAR(dop, Dop);

    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);

    flint_cleanup_master();

    return (status == GR_SUCCESS);
}
