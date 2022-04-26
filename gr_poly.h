#ifndef GR_POLY_H
#define GR_POLY_H

#ifdef GR_POLY_INLINES_C
#define GR_POLY_INLINE FLINT_DLL
#else
#define GR_POLY_INLINE static __inline__
#endif

#include "flint/fmpz_poly.h"
#include "gr.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
    gr_ptr coeffs;
    slong length;
    slong alloc;
}
gr_poly_struct;

typedef gr_poly_struct gr_poly_t[1];

void gr_poly_init(gr_poly_t poly, gr_ctx_t ctx);
void gr_poly_init2(gr_poly_t poly, slong len, gr_ctx_t ctx);
void gr_poly_clear(gr_poly_t poly, gr_ctx_t ctx);

GR_POLY_INLINE gr_ptr
gr_poly_entry_ptr(gr_poly_t poly, slong i, gr_ctx_t ctx)
{
    return GR_ENTRY(poly->coeffs, i, ctx->sizeof_elem);
}

GR_POLY_INLINE slong gr_poly_length(const gr_poly_t poly, gr_ctx_t ctx)
{
    return poly->length;
}

GR_POLY_INLINE void
gr_poly_swap(gr_poly_t poly1, gr_poly_t poly2, gr_ctx_t ctx)
{
    gr_poly_struct t = *poly1;
    *poly1 = *poly2;
    *poly2 = t;
}

void gr_poly_fit_length(gr_poly_t poly, slong len, gr_ctx_t ctx);
void _gr_poly_set_length(gr_poly_t poly, slong len, gr_ctx_t ctx);
void _gr_poly_normalise(gr_poly_t poly, gr_ctx_t ctx);

int gr_poly_set(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx);

int _gr_poly_reverse(gr_ptr res, gr_srcptr poly, slong len, slong n, gr_ctx_t ctx);
int gr_poly_reverse(gr_poly_t res, const gr_poly_t poly, slong n, gr_ctx_t ctx);

GR_POLY_INLINE int
gr_poly_zero(gr_poly_t poly, gr_ctx_t ctx)
{
    _gr_poly_set_length(poly, 0, ctx);
    return GR_SUCCESS;
}

int gr_poly_one(gr_poly_t poly, gr_ctx_t ctx);
int gr_poly_neg_one(gr_poly_t poly, gr_ctx_t ctx);

int gr_poly_write(gr_stream_t out, const gr_poly_t poly, gr_ctx_t ctx);
int gr_poly_print(const gr_poly_t poly, gr_ctx_t ctx);

int gr_poly_randtest(gr_poly_t poly, flint_rand_t state, slong len, gr_ctx_t ctx);

truth_t _gr_poly_equal(gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
truth_t gr_poly_equal(const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx);

int gr_poly_set_scalar(gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx);
int gr_poly_set_si(gr_poly_t poly, slong x, gr_ctx_t ctx);
int gr_poly_set_ui(gr_poly_t poly, slong x, gr_ctx_t ctx);
int gr_poly_set_fmpz(gr_poly_t poly, const fmpz_t x, gr_ctx_t ctx);
int gr_poly_set_fmpq(gr_poly_t poly, const fmpq_t x, gr_ctx_t ctx);

/* todo: document/test */
int gr_poly_get_fmpz_poly(gr_poly_t res, const fmpz_poly_t src, gr_ctx_t ctx);
int gr_poly_set_fmpq_poly(gr_poly_t res, const fmpq_poly_t src, gr_ctx_t ctx);

int gr_poly_set_coeff_scalar(gr_poly_t poly, slong n, gr_srcptr x, gr_ctx_t ctx);
int gr_poly_set_coeff_si(gr_poly_t poly, slong n, slong x, gr_ctx_t ctx);
int gr_poly_set_coeff_ui(gr_poly_t poly, slong n, ulong x, gr_ctx_t ctx);
int gr_poly_set_coeff_fmpz(gr_poly_t poly, slong n, const fmpz_t x, gr_ctx_t ctx);
int gr_poly_set_coeff_fmpq(gr_poly_t poly, slong n, const fmpq_t x, gr_ctx_t ctx);

int gr_poly_neg(gr_poly_t res, const gr_poly_t src, gr_ctx_t ctx);

int _gr_poly_add(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
int gr_poly_add(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx);

int _gr_poly_sub(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
int gr_poly_sub(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx);

int _gr_poly_mul(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
int gr_poly_mul(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx);

int _gr_poly_mullow_generic(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx);
int gr_poly_mullow(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, slong n, gr_ctx_t ctx);

int gr_poly_mul_scalar(gr_poly_t res, const gr_poly_t poly, gr_srcptr c, gr_ctx_t ctx);

int _gr_poly_inv_series(gr_ptr Qinv, gr_srcptr Q, slong Qlen, slong len, gr_ctx_t ctx);
int gr_poly_inv_series(gr_poly_t Qinv, const gr_poly_t Q, slong len, gr_ctx_t ctx);

int _gr_poly_evaluate_rectangular(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx);
int gr_poly_evaluate_rectangular(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx);

int _gr_poly_evaluate_horner(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx);
int gr_poly_evaluate_horner(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx);

int _gr_poly_evaluate(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr x, gr_ctx_t ctx);
int gr_poly_evaluate(gr_ptr res, const gr_poly_t poly, gr_srcptr x, gr_ctx_t ctx);

int _gr_poly_taylor_shift_horner(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx);
int gr_poly_taylor_shift_horner(gr_poly_t res, const gr_poly_t f, gr_srcptr c, gr_ctx_t ctx);

int _gr_poly_taylor_shift_divconquer(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx);
int gr_poly_taylor_shift_divconquer(gr_poly_t res, const gr_poly_t f, gr_srcptr c, gr_ctx_t ctx);

int _gr_poly_taylor_shift(gr_ptr res, gr_srcptr poly, slong len, gr_srcptr c, gr_ctx_t ctx);
int gr_poly_taylor_shift(gr_poly_t res, const gr_poly_t f, gr_srcptr c, gr_ctx_t ctx);

int _gr_poly_compose_horner(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
int gr_poly_compose_horner(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx);

int _gr_poly_compose_divconquer(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
int gr_poly_compose_divconquer(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx);

int _gr_poly_compose(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, gr_ctx_t ctx);
int gr_poly_compose(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx);

int _gr_poly_derivative(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx);
int gr_poly_derivative(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx);

int _gr_poly_integral(gr_ptr res, gr_srcptr poly, slong len, gr_ctx_t ctx);
int gr_poly_integral(gr_poly_t res, const gr_poly_t poly, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
