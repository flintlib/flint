.. _gr-domains:

**gr.h (continued)** -- builtin domains and types
===============================================================================

Coercions
-------------------------------------------------------------------------------

.. function:: int gr_ctx_cmp_coercion(gr_ctx_t ctx1, gr_ctx_t ctx2)

    Returns 1 if coercing elements into *ctx1* is more meaningful,
    and returns -1 otherwise.



Domain properties
-------------------------------------------------------------------------------

.. function:: truth_t gr_ctx_is_finite(gr_ctx_t ctx)
              truth_t gr_ctx_is_multiplicative_group(gr_ctx_t ctx)
              truth_t gr_ctx_is_ring(gr_ctx_t ctx)
              truth_t gr_ctx_is_commutative_ring(gr_ctx_t ctx)
              truth_t gr_ctx_is_integral_domain(gr_ctx_t ctx)
              truth_t gr_ctx_is_unique_factorization_domain(gr_ctx_t ctx)
              truth_t gr_ctx_is_field(gr_ctx_t ctx)
              truth_t gr_ctx_is_algebraically_closed(gr_ctx_t ctx)
              truth_t gr_ctx_is_finite_characteristic(gr_ctx_t ctx)
              truth_t gr_ctx_is_ordered_ring(gr_ctx_t ctx)
              truth_t gr_ctx_is_zero_ring(gr_ctx_t ctx)

    Returns whether the structure satisfies the respective
    mathematical property.
    The result can be ``T_UNKNOWN``.

.. function:: truth_t gr_ctx_is_exact(gr_ctx_t ctx)

    Returns whether the representation of elements is always exact.

.. function:: truth_t gr_ctx_is_canonical(gr_ctx_t ctx)

    Returns whether the representation of elements is always canonical.

.. function:: truth_t gr_ctx_has_real_prec(gr_ctx_t ctx)

    Returns whether *ctx* or a base field thereof represents real or complex
    numbers using finite-precision approximations.
    This returns ``T_TRUE`` both for floating-point approximate
    fields and for rigorous fields based on ball or interval arithmetic.

.. function:: int gr_ctx_set_real_prec(gr_ctx_t ctx, slong prec)
              int gr_ctx_get_real_prec(slong * prec, gr_ctx_t ctx)

    Sets or retrieves the floating-point precision in bits.


Groups
-------------------------------------------------------------------------------

.. function:: void gr_ctx_init_perm(gr_ctx_t ctx, ulong n)

    Initializes *ctx* to the symmetric group `S_n` representing
    permutations of `[0, 1, \ldots, n - 1]`.
    Elements are currently represented as pointers (the representation
    may change in the future).

.. function:: void gr_ctx_init_psl2z(gr_ctx_t ctx)

    Initializes *ctx* to the modular group `\text{PSL}(2, \mathbb{Z})`
    with elements of type :type:`psl2z_t`.

.. function:: int gr_ctx_init_dirichlet_group(gr_ctx_t ctx, ulong q)

    Initializes *ctx* to the Dirichlet group `G_q`
    with elements of type :type:`dirichlet_char_t`.
    Fails and returns ``GR_DOMAIN`` if *q* is zero.
    Fails and returns ``GR_UNABLE`` if *q* has a prime factor
    larger than `10^{16}`, which is currently unsupported
    by the implementation.

Base rings and fields
-------------------------------------------------------------------------------

.. function:: void gr_ctx_init_random(gr_ctx_t ctx, flint_rand_t state)

    Initializes *ctx* to a random ring. This will currently
    only generate base rings and composite rings over certain
    simple base rings.

.. function:: void gr_ctx_init_fmpz(gr_ctx_t ctx)

    Initializes *ctx* to the ring of integers
    `\mathbb{Z}` with elements of type :type:`fmpz`.

.. function:: void gr_ctx_init_fmpq(gr_ctx_t ctx)

    Initializes *ctx* to the field of rational numbers
    `\mathbb{Q}` with elements of type :type:`fmpq`.

.. function:: void gr_ctx_init_fmpzi(gr_ctx_t ctx)

    Initializes *ctx* to the ring of Gaussian integers
    `\mathbb{Z}[i]` with elements of type :type:`fmpzi_t`.

.. function:: void gr_ctx_init_nmod8(gr_ctx_t ctx, unsigned char n)
              void gr_ctx_init_nmod32(gr_ctx_t ctx, unsigned int n)

    Initializes *ctx* to the ring `\mathbb{Z}/n\mathbb{Z}`
    of integers modulo *n* where
    elements have type :type:`uint8` or :type:`uint32`. The modulus must be
    nonzero.

.. function:: void gr_ctx_init_nmod(gr_ctx_t ctx, ulong n)

    Initializes *ctx* to the ring `\mathbb{Z}/n\mathbb{Z}`
    of integers modulo *n* where
    elements have type :type:`ulong`. We require `n \ne 0`.

.. function:: void gr_ctx_init_fmpz_mod(gr_ctx_t ctx, const fmpz_t n)

    Initializes *ctx* to the ring `\mathbb{Z}/n\mathbb{Z}`
    of integers modulo *n* where
    elements have type :type:`fmpz`. The modulus must be positive.

.. function:: void gr_ctx_fmpz_mod_set_primality(gr_ctx_t ctx, truth_t is_prime)

    For a ring initialized with :func:`gr_ctx_init_fmpz_mod`,
    indicate whether the modulus is prime. This can speed up
    some computations.

.. function:: void gr_ctx_init_fq(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var)
              void gr_ctx_init_fq_nmod(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var)
              void gr_ctx_init_fq_zech(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var)

    Initializes *ctx* to the finite field `\mathbb{F}_q`
    where `q = p^d`. It is assumed (not checked) that *p* is prime.
    The variable name *var* can be ``NULL`` to use a default.

    The corresponding element types are ``fq_t``, ``fq_nmod_t``, ``fq_zech_t``.
    The ``fq_nmod`` context requires `p < 2^{64}` while ``fq_zech``
    requires `q < 2^{64}` (and in practice a much smaller value
    than this).

.. function:: void gr_ctx_init_nf(gr_ctx_t ctx, const fmpq_poly_t poly)
              void gr_ctx_init_nf_fmpz_poly(gr_ctx_t ctx, const fmpz_poly_t poly)

    Initializes *ctx* to the number field with defining polynomial
    ``poly`` which *must* be irreducible (this is not checked).
    The elements have type :type:`nf_elem_t`.

.. function:: void gr_ctx_init_real_qqbar(gr_ctx_t ctx)
              void gr_ctx_init_complex_qqbar(gr_ctx_t ctx)

    Initializes *ctx* to the field of real or complex algebraic
    numbers with elements of type :type:`qqbar_t`.

.. function:: void gr_ctx_init_real_arb(gr_ctx_t ctx, slong prec)
              void gr_ctx_init_complex_acb(gr_ctx_t ctx, slong prec)

    Initializes *ctx* to the field of real or complex
    numbers represented by elements of type :type:`arb_t`
    and  :type:`acb_t`.

.. function:: void gr_ctx_arb_set_prec(gr_ctx_t ctx, slong prec)
              slong gr_ctx_arb_get_prec(gr_ctx_t ctx)

    Sets or retrieves the bit precision of *ctx* which must be
    an Arb context (this is currently not checked).

.. function:: void gr_ctx_init_real_ca(gr_ctx_t ctx)
              void gr_ctx_init_complex_ca(gr_ctx_t ctx)
              void gr_ctx_init_real_algebraic_ca(gr_ctx_t ctx)
              void gr_ctx_init_complex_algebraic_ca(gr_ctx_t ctx)

    Initializes *ctx* to the field of real, complex, real algebraic
    or complex algebraic numbers represented by elements of type
    :type:`ca_t`.

.. function:: void gr_ctx_ca_set_option(gr_ctx_t ctx, slong option, slong value)
              slong gr_ctx_ca_get_option(gr_ctx_t ctx, slong option)

    Sets or retrieves options of a Calcium context object.

Extended number sets
-------------------------------------------------------------------------------

.. function:: void gr_ctx_init_complex_extended_ca(gr_ctx_t ctx)

    Like :func:`gr_ctx_init_complex_ca` but allows special values
    (infinities, undefined).

Floating-point arithmetic
-------------------------------------------------------------------------------

Although domains of floating-point numbers approximate
real and complex fields, they are not rings or fields.
Floating-point arithmetic can be used in many places where a ring
or field is normally assumed, but predicates like "is field"
return false.

* Equality compares equality of floating-point numbers,
  with the special rule that NaN is not equal to itself.
* In general, the following implementations do not currently
  guarantee correct rounding except for atomic arithmetic operations
  (add, sub, mul, div, sqrt) on real floating-point numbers.

.. function:: void gr_ctx_init_real_float_arf(gr_ctx_t ctx, slong prec)

    Initializes *ctx* to the floating-point arithmetic with elements
    of type :type:`arf_t` and a default precision of *prec* bits.

.. function:: void gr_ctx_init_complex_float_acf(gr_ctx_t ctx, slong prec)

    Initializes *ctx* to the complex floating-point arithmetic with elements
    of type :type:`acf_t` and a default precision of *prec* bits.

Vectors
-------------------------------------------------------------------------------

.. function:: void gr_ctx_init_vector_gr_vec(gr_ctx_t ctx, gr_ctx_t base_type)

    Initializes *ctx* to the domain of all vectors (of any length)
    over the given *base_type*.
    Elements have type :type:`gr_vec_struct`.

.. function:: void gr_ctx_init_vector_space_gr_vec(gr_ctx_t ctx, gr_ctx_t base_type, slong n)

    Initializes *ctx* to the space of all vectors of length *n*
    over the given *base_type*.
    Elements have type :type:`gr_vec_struct`.

Matrices
-------------------------------------------------------------------------------

.. function:: void gr_ctx_init_matrix_domain(gr_ctx_t ctx, gr_ctx_t base_ring)

    Initializes *ctx* to the domain of all matrices (of any shape)
    over the given *base_ring*.
    Elements have type :type:`gr_mat_struct`.

.. function:: void gr_ctx_init_matrix_space(gr_ctx_t ctx, gr_ctx_t base_ring, slong n, slong m)

    Initializes *ctx* to the space of matrices over *base_ring*
    with *n* rows and *m* columns.
    Elements have type :type:`gr_mat_struct`.

.. function:: void gr_ctx_init_matrix_ring(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)

    Initializes *ctx* to the ring of matrices over *base_ring*
    with *n* rows columns.
    Elements have type :type:`gr_mat_struct`.

Polynomial rings
-------------------------------------------------------------------------------

.. function:: void gr_ctx_init_fmpz_poly(gr_ctx_t ctx)

    Initializes *ctx* to a ring of integer polynomials of
    type :type:`fmpz_poly_struct`.

.. function:: void gr_ctx_init_fmpq_poly(gr_ctx_t ctx)

    Initializes *ctx* to a ring of rational polynomials of
    type :type:`fmpq_poly_struct`.

.. function:: void gr_ctx_init_gr_poly(gr_ctx_t ctx, gr_ctx_t base_ring)

    Initializes *ctx* to a ring of densely represented univariate polynomials
    over the given *base_ring*.
    Elements have type :type:`gr_poly_struct`.

.. function:: void gr_ctx_init_fmpz_mpoly(gr_ctx_t ctx, slong nvars, const ordering_t ord)

    Initializes *ctx* to a ring of sparsely represented multivariate
    polynomials in *nvars* variables over the integers,
    with monomial ordering *ord*.
    Elements have type :type:`fmpz_mpoly_struct`.

.. function:: void gr_ctx_init_gr_mpoly(gr_ctx_t ctx, gr_ctx_t base_ring, slong nvars, const ordering_t ord)

    Initializes *ctx* to a ring of sparsely represented multivariate
    polynomials in *nvars* variables over the given *base_ring*,
    with monomial ordering *ord*.
    Elements have type :type:`gr_mpoly_struct`.

Fraction fields
-------------------------------------------------------------------------------

.. function:: void gr_ctx_init_fmpz_mpoly_q(gr_ctx_t ctx, slong nvars, const ordering_t ord)

    Initializes *ctx* to a ring of sparsely represented multivariate
    fractions in *nvars* variables over the integers (equivalently, rationals),
    with monomial ordering *ord*.
    Elements have type :type:`fmpz_mpoly_q_struct`.

Symbolic expressions
-------------------------------------------------------------------------------

.. function:: void gr_ctx_init_fexpr(gr_ctx_t ctx)

    Initializes *ctx* to handle symbolic expressions.
    Elements have type :type:`fexpr_struct`.

.. raw:: latex

    \newpage
