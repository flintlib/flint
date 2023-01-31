
#include "acb_theta.h"

void
arb_mat_randtest_sym_pos(arb_mat_t mat, flint_rand_t state, slong prec,
                         slong mag_bits)
{
    slong g = arb_mat_nrows(mat);
    arb_mat_t tp;

    arb_mat_init(tp, g, g);
    arb_mat_randtest_cho(mat, state, prec, mag_bits);
    arb_mat_transpose(tp, mat);
    arb_mat_mul(mat, tp, mat, prec);

    arb_mat_clear(tp);
}
