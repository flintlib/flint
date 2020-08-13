.. _aprcl:

**aprcl.h** -- APRCL primality testing
========================================================================================

This module implements the rigorous APRCL primality test, suitable for integers
up to a few thousand digits.

Authors:

* Vladimir Glazachev (Google Summer of Code, 2015)

Primality test functions
--------------------------------------------------------------------------------

.. function:: int aprcl_is_prime(const fmpz_t n)

    Tests `n` for primality using the APRCL test.
    This is the same as :func:`aprcl_is_prime_jacobi`.

.. function:: int aprcl_is_prime_jacobi(const fmpz_t n)

    If `n` prime returns 1; otherwise returns 0. The algorithm is well described
    in "Implementation of a New Primality Test" by H. Cohen and A.K. Lenstra and
    "A Course in Computational Algebraic Number Theory" by H. Cohen.

    It is theoretically possible that this function fails to prove that
    `n` is prime. In this event, :func:`flint_abort` is called.
    To handle this condition, the :func:`_aprcl_is_prime_jacobi` function
    can be used.

.. function:: int aprcl_is_prime_gauss(const fmpz_t n)

    If `n` is prime returns 1; otherwise returns 0.
    Uses the Cyclotomic primality testing algorithm described in
    "Four primality testing algorithms" by Rene Schoof.
    The minimum required numbers `s` and `R` are computed automatically.

    By default `R \ge 180`. In some cases this function fails to prove
    that `n` is prime. This means that we select a too small `R` value.
    In this event, :func:`flint_abort` is called.
    To handle this condition, the :func:`_aprcl_is_prime_jacobi` function
    can be used.

.. function:: primality_test_status _aprcl_is_prime_jacobi(const fmpz_t n, const aprcl_config config)

    Jacobi sum test for `n`. Possible return values:
    ``PRIME``, ``COMPOSITE`` and ``UNKNOWN`` (if we cannot
    prove primality).

.. function:: primality_test_status _aprcl_is_prime_gauss(const fmpz_t n, const aprcl_config config)

    Tests `n` for primality with fixed ``config``. Possible return values:
    ``PRIME``, ``COMPOSITE`` and ``PROBABPRIME``
    (if we cannot prove primality).

.. function:: void aprcl_is_prime_gauss_min_R(const fmpz_t n, ulong R)

    Same as :func:`aprcl_is_prime_gauss` with fixed minimum value of `R`.

.. function:: int aprcl_is_prime_final_division(const fmpz_t n, const fmpz_t s, ulong r)

    Returns 0 if for some `a = n^k \bmod s`, where `k \in [1, r - 1]`, 
    we have that `a | n`; otherwise returns 1.

Configuration functions
--------------------------------------------------------------------------------

.. type:: _aprcl_config

.. type:: aprcl_config

    Holds precomputed parameters.

.. function:: void aprcl_config_gauss_init(aprcl_config conf, const fmpz_t n)

    Computes the `s` and `R` values used in the cyclotomic primality test,
    `s^2 > n` and `s=\prod\limits_{\substack{q-1|R \\ q \text{ prime}}}q`.
    Also stores factors of `R` and `s`.

.. function:: void aprcl_config_gauss_init_min_R(aprcl_config conf, const fmpz_t n, ulong R)

    Computes the `s` with fixed minimum `R` such that `a^R \equiv 1 \mod{s}`
    for all integer `a` coprime to `s`. 

.. function:: void aprcl_config_gauss_clear(aprcl_config conf)

    Clears the given ``aprcl_config`` element. It must be reinitialised in
    order to be used again.

.. function:: ulong aprcl_R_value(const fmpz_t n)

    Returns a precomputed `R` value for APRCL, such that the
    corresponding `s` value is greater than `\sqrt{n}`. The maximum
    stored value `6983776800` allows to test numbers up to `6000` digits.

.. function:: void aprcl_config_jacobi_init(aprcl_config conf, const fmpz_t n)

    Computes the `s` and `R` values used in the cyclotomic primality test,
    `s^2 > n` and `a^R \equiv 1 \mod{s}` for all `a` coprime to `s`.
    Also stores factors of `R` and `s`.

.. function:: void aprcl_config_jacobi_clear(aprcl_config conf)

    Clears the given ``aprcl_config`` element. It must be reinitialised in
    order to be used again.

Cyclotomic arithmetic
--------------------------------------------------------------------------------

This code implements arithmetic in cyclotomic rings.

Types
................................................................................

.. type:: _unity_zp

.. type:: unity_zp

    Represents an element of `\mathbb{Z}[\zeta_{p^{exp}}]/(n)` as an
    :type:`fmpz_mod_poly_t` reduced modulo a cyclotomic polynomial.

.. type:: _unity_zpq

.. type:: unity_zpq

    Represents an element of `\mathbb{Z}[\zeta_q, \zeta_p]/(n)`
    as an array of :type:`fmpz_mod_poly_t`.

Memory management
................................................................................

.. function:: void unity_zp_init(unity_zp f, ulong p, ulong exp, const fmpz_t n)

    Initializes `f` as an element of `\mathbb{Z}[\zeta_{p^{exp}}]/(n)`.

.. function:: void unity_zp_clear(unity_zp f)

    Clears the given element. It must be reinitialised in
    order to be used again.

.. function:: void unity_zp_copy(unity_zp f, const unity_zp g)

    Sets `f` to `g`. `f` and `g` must be initialized with same `p` and `n`.

.. function:: void unity_zp_swap(unity_zp f, unity_zp q)

    Swaps `f` and `g`. `f` and `g` must be initialized with same `p` and `n`.

.. function:: void unity_zp_set_zero(unity_zp f)

    Sets `f` to zero.

Comparison
................................................................................

.. function:: slong unity_zp_is_unity(const unity_zp f)

    If `f = \zeta^h` returns h; otherwise returns -1.

.. function:: int unity_zp_equal(const unity_zp f, const unity_zp g)

    Returns nonzero if `f = g` reduced by the `p^{exp}`-th cyclotomic
    polynomial.

Output
................................................................................

.. function:: void unity_zp_print(const unity_zp f)

    Prints the contents of the `f`.

Coefficient management
................................................................................

.. function:: void unity_zp_coeff_set_fmpz(unity_zp f, ulong ind, const fmpz_t x)
              void unity_zp_coeff_set_ui(unity_zp f, ulong ind, ulong x)

    Sets the coefficient of `\zeta^{ind}` to `x`.
    `ind` must be less than `p^{exp}`.

.. function:: void unity_zp_coeff_add_fmpz(unity_zp f, ulong ind, const fmpz_t x)
              void unity_zp_coeff_add_ui(unity_zp f, ulong ind, ulong x)

    Adds `x` to the coefficient of `\zeta^{ind}`.
    `x` must be less than `n`.
    `ind` must be less than `p^{exp}`.

.. function:: void unity_zp_coeff_inc(unity_zp f, ulong ind)

    Increments the coefficient of `\zeta^{ind}`.
    `ind` must be less than `p^{exp}`.

.. function:: void unity_zp_coeff_dec(unity_zp f, ulong ind)

    Decrements the coefficient of `\zeta^{ind}`.
    `ind` must be less than `p^{exp}`.

Scalar multiplication
................................................................................

.. function:: void unity_zp_mul_scalar_fmpz(unity_zp f, const unity_zp g, const fmpz_t s)

    Sets `f` to `s \cdot g`. `f` and `g` must be initialized with
    same `p`, `exp` and `n`.

.. function:: void unity_zp_mul_scalar_ui(unity_zp f, const unity_zp g, ulong s)

    Sets `f` to `s \cdot g`. `f` and `g` must be initialized with
    same `p`, `exp` and `n`.

Addition and multiplication
................................................................................

.. function:: void unity_zp_add(unity_zp f, const unity_zp g, const unity_zp h)

    Sets `f` to `g + h`.
    `f`, `g` and `h` must be initialized with same `p`, `exp` and `n`.

.. function:: void unity_zp_mul(unity_zp f, const unity_zp g, const unity_zp h)

    Sets `f` to `g \cdot h`.
    `f`, `g` and `h` must be initialized with same `p`, `exp` and `n`.

.. function:: void unity_zp_sqr(unity_zp f, const unity_zp g)

    Sets `f` to `g \cdot g`.
    `f`, `g` and `h` must be initialized with same `p`, `exp` and `n`.

.. function:: void unity_zp_mul_inplace(unity_zp f, const unity_zp g, const untiy_zp h, fmpz_t * t)

    Sets `f` to `g \cdot h`. If `p^{exp} = 3, 4, 5, 7, 8, 9, 11, 16` special
    multiplication functions are used. The preallocated array `t` of ``fmpz_t`` is
    used for all computations in this case.
    `f`, `g` and `h` must be initialized with same `p`, `exp` and `n`.

.. function:: void unity_zp_sqr_inplace(unity_zp f, const unity_zp g, fmpz_t * t)

    Sets `f` to `g \cdot g`. If `p^{exp} = 3, 4, 5, 7, 8, 9, 11, 16` special
    multiplication functions are used. The preallocated array `t` of ``fmpz_t`` is
    used for all computations in this case.
    `f` and `g` must be initialized with same `p`, `exp` and `n`.

Powering functions
................................................................................

.. function:: void unity_zp_pow_fmpz(unity_zp f, unity_zp g, const fmpz_t pow)

    Sets `f` to `g^{pow}`. `f` and `g` must be initialized with
    same `p`, `exp` and `n`.

.. function:: void unity_zp_pow_ui(unity_zp f, unity_zp g, ulong pow)

    Sets `f` to `g^{pow}`. `f` and `g` must be initialized with
    same `p`, `exp` and `n`.

.. function:: ulong _unity_zp_pow_select_k(const fmpz_t n)

    Returns the smallest integer `k` satisfying
    `\log (n) < (k(k + 1)2^{2k}) / (2^{k + 1} - k - 2) + 1`

.. function:: void unity_zp_pow_2k_fmpz(unity_zp f, unity_zp g, const fmpz_t pow)

    Sets `f` to `g^{pow}` using the `2^k`-ary exponentiation method.
    `f` and `g` must be initialized with same `p`, `exp` and `n`.

.. function:: void unity_zp_pow_2k_ui(unity_zp f, const unity_zp g, ulong pow)

    Sets `f` to `g^{pow}` using the `2^k`-ary exponentiation method.
    `f` and `g` must be initialized with same `p`, `exp` and `n`.

.. function:: void unity_zp_pow_sliding_fmpz(unity_zp f, unity_zp g, const fmpz_t pow)

    Sets `f` to `g^{pow}` using the sliding window exponentiation method.
    `f` and `g` must be initialized with same `p`, `exp` and `n`.


Cyclotomic reduction
................................................................................

.. function:: void _unity_zp_reduce_cyclotomic_divmod(unity_zp f)
              void _unity_zp_reduce_cyclotomic(unity_zp f)

    Sets `f = f \bmod \Phi_{p^{exp}}`. `\Phi_{p^{exp}}` is the `p^{exp}`-th
    cyclotomic polynomial. `g` must be reduced by `x^{p^{exp}}-1` poly.
    `f` and `g` must be initialized with same `p`, `exp` and `n`.

.. function:: void unity_zp_reduce_cyclotomic(unity_zp f, const unity_zp g)

    Sets `f = g \bmod \Phi_{p^{exp}}`. `\Phi_{p^{exp}}` is the `p^{exp}`-th
    cyclotomic polynomial.

Automorphism and inverse
................................................................................

.. function:: void unity_zp_aut(unity_zp f, const unity_zp g, ulong x)

    Sets `f = \sigma_x(g)`, the automorphism `\sigma_x(\zeta)=\zeta^x`.
    `f` and `g` must be initialized with the same `p`, `exp` and `n`.

.. function:: void unity_zp_aut_inv(unity_zp f, const unity_zp g, ulong x)

    Sets `f = \sigma_x^{-1}(g)`, so `\sigma_x(f) = g`.
    `g` must be reduced by `\Phi_{p^{exp}}`.
    `f` and `g` must be initialized with the same `p`, `exp` and `n`.

Jacobi sum
................................................................................

Here `\chi_{p, q}` is the character defined by
`\chi_{p, q}(g^x) = \zeta_{p^k}^x`, where `g` is
a primitive root modulo `q`.

.. function:: void unity_zp_jacobi_sum_pq(unity_zp f, ulong q, ulong p)

    Sets `f` to the Jacobi sum `J(p, q) = j(\chi_{p, q}, \chi_{p, q})`.

.. function:: void unity_zp_jacobi_sum_2q_one(unity_zp f, ulong q)

    Sets `f` to the Jacobi sum
    `J_2(q) = j(\chi_{2, q}^{2^{k - 3}}, \chi_{2, q}^{3 \cdot 2^{k - 3}}))^2`.

.. function:: void unity_zp_jacobi_sum_2q_two(unity_zp f, ulong q)

    Sets `f` to the Jacobi sum
    `J_3(1) = j(\chi_{2, q}, \chi_{2, q}, \chi_{2, q}) =
    J(2, q) \cdot j(\chi_{2, q}^2, \chi_{2, q})`.

Extended rings
................................................................................

.. function:: void unity_zpq_init(unity_zpq f, ulong q, ulong p, const fmpz_t n)

    Initializes `f` as an element of `\mathbb{Z}[\zeta_q, \zeta_p]/(n)`.

.. function:: void unity_zpq_clear(unity_zpq f)

    Clears the given element. It must be reinitialized in
    order to be used again.

.. function:: void unity_zpq_copy(unity_zpq f, const unity_zpq g)

    Sets `f` to `g`. `f` and `g` must be initialized with
    same `p`, `q` and `n`.

.. function:: void unity_zpq_swap(unity_zpq f, unity_zpq q)

    Swaps `f` and `g`. `f` and `g` must be initialized with
    same `p`, `q` and `n`.

.. function:: int unity_zpq_equal(const unity_zpq f, const unity_zpq g)

    Returns nonzero if `f = g`.

.. function:: ulong unity_zpq_p_unity(const unity_zpq f)

    If `f = \zeta_p^x` returns `x in [0, p - 1]`; otherwise returns `p`.

.. function:: int unity_zpq_is_p_unity(const unity_zpq f)

    Returns nonzero if `f = \zeta_p^x`.

.. function:: int unity_zpq_is_p_unity_generator(const unity_zpq f)

    Returns nonzero if `f` is a generator of the cyclic group `<\zeta_p>`.

.. function:: void unity_zpq_coeff_set_fmpz(unity_zpq f, ulong i, ulong j, const fmpz_t x)

    Sets the coefficient of `\zeta_q^i \zeta_p^j` to `x`.
    `i` must be less than `q` and `j` must be less than `p`.

.. function:: void unity_zpq_coeff_set_ui(unity_zpq f, ulong i, ulong j, ulong x)

    Sets the coefficient of `\zeta_q^i \zeta_p^j` to `x`.
    `i` must be less than `q` and `j` must be less then `p`.

.. function:: void unity_zpq_coeff_add(unity_zpq f, ulong i, ulong j, const fmpz_t x)

    Adds `x` to the coefficient of `\zeta_p^i \zeta_q^j`. `x` must be less than `n`.

.. function:: void unity_zpq_add(unity_zpq f, const unity_zpq g, const unity_zpq h)

    Sets `f` to `g + h`.
    `f`, `g` and `h` must be initialized with same
    `q`, `p` and `n`.

.. function:: void unity_zpq_mul(unity_zpq f, const unity_zpq g, const unity_zpq h)

    Sets the `f` to `g \cdot h`.
    `f`, `g` and `h` must be initialized with same
    `q`, `p` and `n`.

.. function:: void _unity_zpq_mul_unity_p(unity_zpq f)

    Sets `f = f \cdot \zeta_p`.

.. function:: void unity_zpq_mul_unity_p_pow(unity_zpq f, const unity_zpq g, ulong k)

    Sets `f` to `g \cdot \zeta_p^k`.

.. function:: void unity_zpq_pow(unity_zpq f, unity_zpq g, const fmpz_t p)

    Sets `f` to `g^p`. `f` and `g` must be initialized with same `p`, `q` and `n`.

.. function:: void unity_zpq_pow_ui(unity_zpq f, unity_zpq g, ulong p)

    Sets `f` to `g^p`. `f` and `g` must be initialized with same `p`, `q` and `n`.

.. function:: void unity_zpq_gauss_sum(unity_zpq f, ulong q, ulong p)

    Sets `f = \tau(\chi_{p, q})`.

.. function:: void unity_zpq_gauss_sum_sigma_pow(unity_zpq f, ulong q, ulong p)

    Sets `f = \tau^{\sigma_n}(\chi_{p, q})`.
