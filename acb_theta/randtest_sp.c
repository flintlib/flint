
#include "acb_theta.h"

static void
randtest_trig_sp(fmpz_mat_t mat, flint_rand_t state, slong bits)
{
    slong g = fmpz_mat_nrows(mat) / 2;
    fmpz_mat_t b, bt;

    fmpz_mat_init(b, g, g);
    fmpz_mat_init(bt, g, g);
    bits = FLINT_MAX(bits, 1);

    fmpz_mat_randbits(b, state, bits);
    fmpz_mat_transpose(bt, b);
    fmpz_mat_add(b, b, bt);
    fmpz_mat_scalar_tdiv_q_2exp(b, b, 1);
    fmpz_mat_trig_sp(mat, b);

    fmpz_mat_clear(b);
    fmpz_mat_clear(bt);
}

static void
randtest_diag_sp(fmpz_mat_t mat, flint_rand_t state, slong bits)
{
    slong g = fmpz_mat_nrows(mat) / 2;
    fmpz_mat_t u;

    fmpz_mat_init(u, g, g);
    bits = FLINT_MAX(bits, 1);

    fmpz_mat_one(u);
    fmpz_mat_randops(u, state, 2 * bits * g);
    fmpz_mat_diag_sp(mat, u);

    fmpz_mat_clear(u);
}

void
fmpz_mat_randtest_sp(fmpz_mat_t mat, flint_rand_t state, slong bits)
{
    slong g = fmpz_mat_nrows(mat) / 2;
    fmpz_mat_t aux;

    fmpz_mat_init(aux, 2 * g, 2 * g);

    randtest_trig_sp(mat, state, bits);
    randtest_diag_sp(aux, state, bits);
    fmpz_mat_mul(mat, mat, aux);
    fmpz_mat_J(aux);
    fmpz_mat_mul(mat, mat, aux);
    randtest_trig_sp(aux, state, bits);
    fmpz_mat_mul(mat, mat, aux);
    fmpz_mat_J(aux);
    fmpz_mat_mul(mat, mat, aux);
    randtest_diag_sp(aux, state, bits);
    fmpz_mat_mul(mat, mat, aux);
    fmpz_mat_J(aux);
    fmpz_mat_mul(mat, mat, aux);
    randtest_trig_sp(aux, state, bits);
    fmpz_mat_mul(mat, mat, aux);

    fmpz_mat_clear(aux);
}
