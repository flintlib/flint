#include <stdlib.h>
#include <stdio.h>

#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpq_vec.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include "flint/gr.h"
#include "flint/gr_poly.h"
#include "flint/gr_vec.h"
#include <flint/profiler.h>


int main(int argc, char* argv[])
{
    gr_ctx_t ctx;
    gr_ctx_init_fmpq(ctx);

    flint_rand_t state;
    flint_rand_init(state);

    char *strP;
    gr_poly_t P;
    gr_vec_t x;
    gr_vec_t y;
    slong n, j;


    n = 200;

    gr_vec_init(x, n, ctx);
    gr_vec_init(y, n, ctx);
    gr_poly_init(P, ctx);
    for (j = 0; j < n; j++) {
            gr_set_si(gr_vec_entry_ptr(x, j, ctx), n_randint(state, 10000000000), ctx);
            gr_set_si(gr_vec_entry_ptr(y, j, ctx), n_randint(state, 10000000000), ctx);
    }
    //fmpq_poly_interpolate_newton_fmpz_vec(P, x, y, n);


    timeit_t t0;

    timeit_start(t0);
    int plop = gr_poly_interpolate_fast(P, x, y, ctx);
    timeit_stop(t0);
    flint_printf("%d, %d\n", plop, GR_SUCCESS);
    flint_printf("cpu = %ld ms wall = %ld ms\n", t0->cpu, t0->wall);
    /*strP = fmpq_poly_get_str_pretty(P, "t");
    flint_printf("Poly: %s\n", strP);*/
    //gr_poly_print(P, ctx);
    //flint_printf("\n");

    flint_free(strP);
    gr_poly_clear(P, ctx);
    gr_vec_clear(x, ctx);
    gr_vec_clear(y, ctx);

    return 0;
}