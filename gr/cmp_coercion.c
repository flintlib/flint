#include "gr.h"

/* todo: deal correctly with nested structures */

int gr_ctx_cmp_coercion(gr_ctx_t ctx1, gr_ctx_t ctx2)
{
    if (ctx1->which_ring < ctx2->which_ring)
        return -1;
    if (ctx1->which_ring > ctx2->which_ring)
        return 1;

    if (ctx1->which_ring == GR_CTX_GR_POLY)
    {
        return gr_ctx_cmp_coercion(POLYNOMIAL_ELEM_CTX(ctx1), POLYNOMIAL_ELEM_CTX(ctx2));
    }

    if (ctx1->which_ring == GR_CTX_GR_MAT)
    {
        return gr_ctx_cmp_coercion(MATRIX_CTX(ctx1)->base_ring, MATRIX_CTX(ctx2)->base_ring);
    }

    return 1;
}
