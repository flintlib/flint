#include <math.h>

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

    // Differential operator

    // XXX fmpz_poly? fmpq_poly??
    // XXX std derivative
    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_gr_poly(Pol, ZZ);
    gr_ctx_init_gr_ore_poly(Dop, Pol, 0, ORE_ALGEBRA_EULER_DERIVATIVE);
    status |= gr_ctx_set_gen_name(Pol, "z");
    status |= gr_ctx_set_gen_name(Dop, "Tz");
    GR_TMP_INIT(dop, Dop);
    status |= gr_ore_poly_set_str(dop,  "(z^2 + 1)*Tz^2 + (z^2 - 1)*Tz", Dop);
    // status |= gr_ore_poly_set_str(dop,  "4*Tz^2 - 4*Tz - z^2 + 8*z - 11", Dop);
    // status |= gr_ore_poly_set_str(dop, "Tz^6 - 6*Tz^5 + 12*Tz^4 - 10*Tz^3 + 3*Tz^2 + z^2", Dop);
    // status |= gr_ore_poly_set_str(dop, "4*Tz^2 - 4*Tz - z^2 + 8*z - 11", Dop);
    GR_MUST_SUCCEED(status);

    // flint_printf("dop = %{gr}\n", dop, Dop);

    // Evaluation point

    acb_t pt;
    acb_init(pt);
    // acb_set_d(pt, .5);
    acb_const_pi(pt, prec + 4);
    arb_one(acb_imagref(pt));
    acb_inv(pt, pt, prec + 4);

    // Output matrix

    acb_mat_t mat;
    slong dop_order = gr_ore_poly_length(dop, Dop) - 1;
    acb_mat_init(mat, dop_order, dop_order);

    // Numerical solution

    status |= acb_ode_fundamental_matrix(mat, dop, Dop, NULL, NULL, pt, 0, prec);
    GR_MUST_SUCCEED(status);

    flint_printf("\n------\n");
    acb_mat_printd(mat, prec+4);
    flint_printf("\n%f\n", prec/log2(10.));

    acb_clear(pt);
    acb_mat_clear(mat);
    GR_TMP_CLEAR(dop, Dop);

    gr_ctx_clear(Dop);
    gr_ctx_clear(Pol);

    flint_cleanup_master();
}
