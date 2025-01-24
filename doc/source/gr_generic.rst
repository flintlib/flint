.. _gr-generic:

**gr_generic.h** -- basic algorithms and fallback implementations for generic elements
======================================================================================

.. function:: void gr_generic_init(void)
              void gr_generic_clear(void)
              void gr_generic_swap(void)
              void gr_generic_randtest(void)
              void gr_generic_write(void)
              void gr_generic_zero(void)
              void gr_generic_one(void)
              void gr_generic_equal(void)
              void gr_generic_set(void)
              void gr_generic_set_si(void)
              void gr_generic_set_ui(void)
              void gr_generic_set_fmpz(void)
              void gr_generic_neg(void)
              void gr_generic_add(void)
              void gr_generic_sub(void)
              void gr_generic_mul(void)

.. function:: int gr_generic_ctx_clear(gr_ctx_t ctx)

.. function:: void gr_generic_set_shallow(gr_ptr res, gr_srcptr x, const gr_ctx_t ctx)

.. function:: int gr_generic_write_n(gr_stream_t out, gr_srcptr x, slong n, gr_ctx_t ctx)

.. function:: int gr_generic_randtest_not_zero(gr_ptr x, flint_rand_t state, gr_ctx_t ctx)

.. function:: int gr_generic_randtest_small(gr_ptr x, flint_rand_t state, gr_ctx_t ctx)

.. function:: truth_t gr_generic_is_zero(gr_srcptr x, gr_ctx_t ctx)
              truth_t gr_generic_is_one(gr_srcptr x, gr_ctx_t ctx)
              truth_t gr_generic_is_neg_one(gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_generic_neg_one(gr_ptr res, gr_ctx_t ctx)

.. function:: int gr_generic_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t xctx, gr_ctx_t ctx)
              int gr_generic_set_fmpq(gr_ptr res, const fmpq_t y, gr_ctx_t ctx)

Generic string parsing
-----------------------------------------------------------------------------------------

.. macro :: GR_PARSE_BALANCE_ADDITIONS

.. macro :: GR_PARSE_RING_EXPONENTS

.. function:: int gr_generic_set_str_expr(gr_ptr res, const char * s, int flags, gr_ctx_t ctx)
              int gr_generic_set_str(gr_ptr res, const char * s, gr_ctx_t ctx)
              int gr_generic_set_str_balance_additions(gr_ptr res, const char * s, gr_ctx_t ctx)
              int gr_generic_set_str_ring_exponents(gr_ptr res, const char * s, gr_ctx_t ctx)

    Parses expression string. Generators returned by :func:`gr_gens_recursive` are handled
    automatically. We have the following flags:

    * ``GR_PARSE_RING_EXPONENTS`` - by default, only (nonnegative) integer literals are allowed
      for exponents. If this flag is set, exponents are parsed as arbitrary subexpressions
      within the same ring.
    * ``GR_PARSE_BALANCE_ADDITIONS`` - attempt to improve performance for huge sums
      by reordering additions (useful for polynomials)

Generic arithmetic
-----------------------------------------------------------------------------------------

.. function:: int gr_generic_add_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_generic_add_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_generic_add_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_generic_add_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
              int gr_generic_add_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
              int gr_generic_other_add(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

.. function:: int gr_generic_sub_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_generic_sub_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_generic_sub_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_generic_sub_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
              int gr_generic_sub_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
              int gr_generic_other_sub(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

.. function:: int gr_generic_mul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_generic_mul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_generic_mul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_generic_mul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
              int gr_generic_mul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
              int gr_generic_other_mul(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

.. function:: int gr_generic_addmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_generic_addmul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_generic_addmul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_generic_addmul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_generic_addmul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
              int gr_generic_addmul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)

.. function:: int gr_generic_submul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_generic_submul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_generic_submul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_generic_submul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_generic_submul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
              int gr_generic_submul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)

.. function:: int gr_generic_mul_two(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_generic_sqr(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_generic_mul_2exp_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_generic_mul_2exp_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)

.. function:: int gr_generic_set_fmpz_2exp_fmpz(gr_ptr res, const fmpz_t x, const fmpz_t y, gr_ctx_t ctx)

.. function:: int gr_generic_get_fmpz_2exp_fmpz(fmpz_t res1, fmpz_t res2, gr_ptr x, gr_ctx_t ctx)

.. function:: int gr_generic_inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

.. function:: truth_t gr_generic_is_invertible(gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_generic_div_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
              int gr_generic_div_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
              int gr_generic_div_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
              int gr_generic_div_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
              int gr_generic_div_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
              int gr_generic_other_div(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

.. function:: int gr_generic_divexact(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

.. function:: int gr_generic_pow_fmpz_sliding(gr_ptr f, gr_srcptr g, const fmpz_t pow, gr_ctx_t ctx)
              int gr_generic_pow_ui_sliding(gr_ptr f, gr_srcptr g, ulong pow, gr_ctx_t ctx)
              int gr_generic_pow_fmpz_binexp(gr_ptr res, gr_srcptr x, const fmpz_t exp, gr_ctx_t ctx)
              int gr_generic_pow_ui_binexp(gr_ptr res, gr_srcptr x, ulong e, gr_ctx_t ctx)

.. function:: int gr_generic_pow_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t e, gr_ctx_t ctx)
              int gr_generic_pow_si(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx)
              int gr_generic_pow_ui(gr_ptr res, gr_srcptr x, ulong e, gr_ctx_t ctx)
              int gr_generic_pow_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
              int gr_generic_pow_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
              int gr_generic_other_pow(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

.. function:: int _gr_fmpz_poly_evaluate_horner(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx)
              int gr_fmpz_poly_evaluate_horner(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx)
              int _gr_fmpz_poly_evaluate_rectangular(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx)
              int gr_fmpz_poly_evaluate_rectangular(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx)
              int _gr_fmpz_poly_evaluate(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx)
              int gr_fmpz_poly_evaluate(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx)

    Sets *res* to the value of the integer polynomial *f* evaluated
    at the argument *x*.

.. function:: int gr_fmpz_mpoly_evaluate_iter(gr_ptr res, const fmpz_mpoly_t f, gr_srcptr x, const fmpz_mpoly_ctx_t mctx, gr_ctx_t ctx)
              int gr_fmpz_mpoly_evaluate_horner(gr_ptr res, const fmpz_mpoly_t f, gr_srcptr x, const fmpz_mpoly_ctx_t mctx, gr_ctx_t ctx)
              int gr_fmpz_mpoly_evaluate(gr_ptr res, const fmpz_mpoly_t f, gr_srcptr x, const fmpz_mpoly_ctx_t mctx, gr_ctx_t ctx)

    Sets *res* to value of the multivariate polynomial *f* (with
    corresponding context object *mctx*) evaluated at the vector
    of arguments in *x*.

.. function:: truth_t gr_generic_is_square(gr_srcptr x, gr_ctx_t ctx)
              int gr_generic_sqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_generic_rsqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

    Currently these methods check for the special values 0 and 1.

.. function:: int gr_generic_numerator(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
              int gr_generic_denominator(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

.. function:: int gr_generic_cmp(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_generic_cmpabs(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
              int gr_generic_cmp_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
              int gr_generic_cmpabs_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)

Generic special functions
-----------------------------------------------------------------------------------------

To do: move to ``gr_special``

.. function:: int gr_generic_bernoulli_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
              int gr_generic_bernoulli_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
              int gr_generic_bernoulli_vec(gr_ptr res, slong len, gr_ctx_t ctx)
              int gr_generic_eulernum_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
              int gr_generic_eulernum_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
              int gr_generic_eulernum_vec(gr_ptr res, slong len, gr_ctx_t ctx)
              int gr_generic_stirling_s1u_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx)
              int gr_generic_stirling_s1_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx)
              int gr_generic_stirling_s2_uiui(gr_ptr res, ulong x, ulong y, gr_ctx_t ctx)
              int gr_generic_stirling_s1u_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx)
              int gr_generic_stirling_s1_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx)
              int gr_generic_stirling_s2_ui_vec(gr_ptr res, ulong x, slong len, gr_ctx_t ctx)


Generic vector methods
-----------------------------------------------------------------------------------------

To do: move to ``gr_vec``

.. function:: void gr_generic_vec_init(gr_ptr vec, slong len, gr_ctx_t ctx)

.. function:: void gr_generic_vec_clear(gr_ptr vec, slong len, gr_ctx_t ctx)

.. function:: void gr_generic_vec_swap(gr_ptr vec1, gr_ptr vec2, slong len, gr_ctx_t ctx)

.. function:: int gr_generic_vec_zero(gr_ptr vec, slong len, gr_ctx_t ctx)

.. function:: int gr_generic_vec_set(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)

.. function:: int gr_generic_vec_neg(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)

.. function:: int gr_generic_vec_normalise(slong * res, gr_srcptr vec, slong len, gr_ctx_t ctx)

.. function:: slong gr_generic_vec_normalise_weak(gr_srcptr vec, slong len, gr_ctx_t ctx)

.. function:: int gr_generic_vec_mul_scalar_2exp_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)

.. function:: int gr_generic_vec_scalar_addmul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)

.. function:: int gr_generic_vec_scalar_submul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)

.. function:: int gr_generic_vec_scalar_addmul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)

.. function:: int gr_generic_vec_scalar_submul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)

.. function:: truth_t gr_generic_vec_equal(gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)

.. function:: int gr_generic_vec_is_zero(gr_srcptr vec, slong len, gr_ctx_t ctx)

.. function:: int gr_generic_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)

.. function:: int gr_generic_vec_dot_rev(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)

.. function:: int gr_generic_vec_dot_ui(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const ulong * vec2, slong len, gr_ctx_t ctx)

.. function:: int gr_generic_vec_dot_si(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const slong * vec2, slong len, gr_ctx_t ctx)

.. function:: int gr_generic_vec_dot_fmpz(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const fmpz * vec2, slong len, gr_ctx_t ctx)

.. function:: int gr_generic_vec_set_powers(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)

.. function:: int gr_generic_vec_reciprocals(gr_ptr res, slong len, gr_ctx_t ctx)

.. function:: int gr_generic_vec_add(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
              int gr_generic_vec_sub(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
              int gr_generic_vec_mul(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
              int gr_generic_vec_div(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
              int gr_generic_vec_divexact(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
              int gr_generic_vec_pow(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
              int gr_generic_vec_add_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
              int gr_generic_vec_sub_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
              int gr_generic_vec_mul_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
              int gr_generic_vec_div_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
              int gr_generic_vec_divexact_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
              int gr_generic_vec_pow_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
              int gr_generic_vec_add_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
              int gr_generic_vec_sub_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
              int gr_generic_vec_mul_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
              int gr_generic_vec_div_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
              int gr_generic_vec_divexact_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
              int gr_generic_vec_pow_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
              int gr_generic_vec_add_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
              int gr_generic_vec_sub_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
              int gr_generic_vec_mul_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
              int gr_generic_vec_div_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
              int gr_generic_vec_divexact_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
              int gr_generic_vec_pow_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
              int gr_generic_vec_add_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
              int gr_generic_vec_sub_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
              int gr_generic_vec_mul_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
              int gr_generic_vec_div_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
              int gr_generic_vec_divexact_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
              int gr_generic_vec_pow_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
              int gr_generic_vec_add_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
              int gr_generic_vec_sub_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
              int gr_generic_vec_mul_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
              int gr_generic_vec_div_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
              int gr_generic_vec_divexact_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
              int gr_generic_vec_pow_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
              int gr_generic_scalar_add_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int gr_generic_scalar_sub_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int gr_generic_scalar_mul_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int gr_generic_scalar_div_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int gr_generic_scalar_divexact_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int gr_generic_scalar_pow_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int gr_generic_vec_add_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
              int gr_generic_vec_sub_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
              int gr_generic_vec_mul_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
              int gr_generic_vec_div_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
              int gr_generic_vec_divexact_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
              int gr_generic_vec_pow_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
              int gr_generic_other_add_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
              int gr_generic_other_sub_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
              int gr_generic_other_mul_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
              int gr_generic_other_div_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
              int gr_generic_other_divexact_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
              int gr_generic_other_pow_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
              int gr_generic_vec_add_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
              int gr_generic_vec_sub_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
              int gr_generic_vec_mul_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
              int gr_generic_vec_div_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
              int gr_generic_vec_divexact_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
              int gr_generic_vec_pow_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
              int gr_generic_scalar_other_add_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int gr_generic_scalar_other_sub_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int gr_generic_scalar_other_mul_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int gr_generic_scalar_other_div_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int gr_generic_scalar_other_divexact_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
              int gr_generic_scalar_other_pow_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)


.. raw:: latex

    \newpage
