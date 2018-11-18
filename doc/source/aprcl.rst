.. _aprcl:

**aprcl.h** -- APRCL primality testing
========================================================================================


Primality test functions
--------------------------------------------------------------------------------


.. function:: int is_prime_aprcl(const fmpz_t n)

    Test `n` for primality using APRCL test. Calls ``is_prime_jacobi()``
    function. Returns non zero value if `n` is prime.

.. function:: primality_test_status _is_prime_jacobi(const fmpz_t n, const aprcl_config config)

    Jacobi sum test for ``fmpz`` value `n`. Possible return values:
    ``PRIME``, ``COMPOSITE`` and ``UNKNOWN`` (if we can't 
    prove primality). For implementation details see ``is_prime_jacobi.c``
    source code.

.. function:: int is_prime_jacobi(const fmpz_t n)

    If `n` prime returns 1; otherwise returns 0. The algorithm is well described
    in "Implementation of a New Primality Test" by H. Cohen and A.K. Lenstra and 
    "A Course in Computational Algebraic Number Theory" by H. Cohen.
    It is theoretically possible that returns 0 for prime number. 
    It means that we can't check `L_p` condition (see ``is_prime_jacobi.c`` 
    source code for implementation details), this case can be find
    using ``_is_prime_jacobi()`` function.

.. function:: primality_test_status _is_prime_gauss(const fmpz_t n, const aprcl_config config)

    Test `n` for primality with fixed ``config``. Possible return values:
    ``PRIME``, ``COMPOSITE`` and ``PROBABPRIME`` 
    (if we can't prove primality).

.. function:: int is_prime_gauss(const fmpz_t n)

    If `n` exactly prime returns 1; otherwise returns 0.
    In function used Cyclotomic primality testing algorithm discribed in 
    "Four primality testing algorithms" by Rene Schoof. 
    The minimum required numbers `s` and `R` computed automatically. 
    By default `R >= 180`. In some cases returns 0 for prime numbers. 
    It means that we select too small R value. To find this case 
    ``_is_prime_gauss()`` function can be used.

.. function:: is_prime_gauss_min_R(const fmpz_t n, ulong R)

    Same as ``is_prime_gauss`` function with fixed minimum value of `R`.

.. function:: int is_prime_final_division(const fmpz_t n, const fmpz_t s, ulong r)

    Returns 0 if for some `a = n^k \bmod s`, where `k \in [1, r - 1]`, 
    we have that `a | n`; otherwise returns 1.



Configuration functions
--------------------------------------------------------------------------------


.. function:: void config_gauss_init(aprcl_config conf, const fmpz_t n)

    Compute the `s` and `R` values used in cyclotomic primality test.
    `s^2 > n` and `s=\prod\limits_{\substack{q-1|R \\ q \text{ prime}}}q`.
    Also store factors of `R` and `s`.

.. function:: void config_gauss_init_min_R(aprcl_config conf, const fmpz_t n, ulong R)

    Compute the `s` with fixed minimum `R` such that `a^R \equiv 1 \mod{s}`
    for all integer `a` coprime to `s`. 

.. function:: void config_gauss_clear(aprcl_config conf)

    Clears the given ``aprcl_config`` element. It must be reinitialised in
    order to be used again.

.. function:: ulong _R_value(const fmpz_t n)

    Returns precomputed `R` value for APR-CL configration. `R` such that the 
    corresponding `s` value greater then `\sqrt{n}`. Maximum stored value:
    `6983776800` allows to test numbers up to `6000` digit.

.. function:: void config_jacobi_init(aprcl_config conf, const fmpz_t n)

    Compute the `s` and `R` values used in cyclotomic primality test.
    `s^2 > n` and `a^R \equiv 1 \mod{s}` for all `a` coprime to `s`.
    Also store factors of `R` and `s`.

.. function:: void config_jacobi_clear(aprcl_config conf)

    Clears the given ``aprcl_config`` element. It must be reinitialised in
    order to be used again.



$\mathbb{Z}[\zeta_p]/(n)$. Memory management
--------------------------------------------------------------------------------

``unity_zp`` represents as ``fmpz_mod_poly`` reduced modulo
cyclotomic polynomial.

.. function:: void unity_zp_init(unity_zp f, ulong p, const fmpz_t n)

    Initialized ``unity_zp`` element of `\mathbb{Z}[\zeta_p]/(n)`.

.. function:: void unity_zp_clear(unity_zp f)

    Clears the given ``unity_zp`` element. It must be reinitialised in
    order to be used again.

.. function:: void unity_zp_copy(unity_zp f, const unity_zp g)

    Sets `f` to `g`. `f` and `g` must be initialized with same `p` and `n`.

.. function:: void unity_zp_swap(unity_zp f, unity_zp q)

    Swaps `f` and `g`. `f` and `g` must be initialized with same `p` and `n`.

.. function:: void unity_zp_set_zero(unity_zp f)

    Sets `f` to zero value.



$\mathbb{Z}[\zeta_p]/(n)$. Comparision
--------------------------------------------------------------------------------


.. function:: slong unity_zp_is_unity(const unity_zp f)

    If `f = \zeta^h` returns h; otherwise returns -1.

.. function:: int unity_zp_equal(const unity_zp f, const unity_zp g)

    Returns non zero value if `f == g` reduced by `p^{exp}`-th cyclotomic
    polynomial.



$\mathbb{Z}[\zeta_p]/(n)$. Output
--------------------------------------------------------------------------------


.. function:: void unity_zp_print(const unity_zp f)

    Prints the contents of the `f`.



$\mathbb{Z}[\zeta_p]/(n)$. Coefficient management
--------------------------------------------------------------------------------


.. function:: void unity_zp_coeff_set_fmpz(unity_zp f, ulong ind, const fmpz_t x)

    Sets the coefficient of the `\zeta^{ind}` equal to x.
    `ind` must be less then `p^{exp}`.

.. function:: void unity_zp_coeff_set_ui(unity_zp f, ulong ind, ulong x)

    Sets the coefficient of the `\zeta^{ind}` equal to x.
    `ind` must be less then `p^{exp}`.

.. function:: void unity_zp_coeff_add_fmpz(unity_zp f, ulong ind, const fmpz_t x)

    Sets the `a` in `a*\zeta^{ind}` to `a + x`. `x` must be less then `n`.
    `ind` must be less then `p^{exp}`.

.. function:: void unity_zp_coeff_add_ui(unity_zp f, ulong ind, ulong x)

    Sets the `a` in `a*\zeta^{ind}` to `a + x`. `x` must be less then `n`.
    `ind` must be less then `p^{exp}`.

.. function:: void unity_zp_coeff_inc(unity_zp f, ulong ind)

    Increase the coefficient at `\zeta^{ind}`. 
    `ind` must be less then `p^{exp}`.

.. function:: void unity_zp_coeff_dec(unity_zp f, ulong ind)

    Decrease the coefficient at `\zeta^{ind}`. 
    `ind` must be less then `p^{exp}`.



$\mathbb{Z}[\zeta_p]/(n)$. Scalar multiplication
--------------------------------------------------------------------------------


.. function:: void unity_zp_mul_scalar_fmpz(unity_zp f, const unity_zp g, const fmpz_t s)

    Sets the `f` to `s * g`. `f` and `g` must be initialized with
    same `p`, `exp` and `n`.

.. function:: void unity_zp_mul_scalar_ui(unity_zp f, const unity_zp g, ulong s)

    Sets the `f` to `s * g`. `f` and `g` must be initialized with
    same `p`, `exp` and `n`.



$\mathbb{Z}[\zeta_p]/(n)$. Addition and multiplication
--------------------------------------------------------------------------------


.. function:: void unity_zp_add(unity_zp f, const unity_zp g, const unity_zp h)

    Sets the `f` to `g + h`.
    `f`, `g` and `h` must be initialized with same `p`, `exp` and `n`.

.. function:: void unity_zp_mul(unity_zp f, const unity_zp g, const unity_zp h)

    Sets `f` to `g * h`.
    `f`, `g` and `h` must be initialized with same `p`, `exp` and `n`.

.. function:: void unity_zp_sqr(unity_zp f, const unity_zp g)

    Sets `f` to `g * g`.
    `f`, `g` and `h` must be initialized with same `p`, `exp` and `n`.

.. function:: void untiy_zp_mul_inplace(unity_zp f, const unity_zp g, const untiy_zp h, fmpz_t * t)

    Sets `f` to `g * h`. If `p^{exp} = 3, 4, 5, 7, 8, 9, 11, 16` special
    multiplication functions are used. Allocated array `t` of ``fmpz_t`` are
    used for all computations in this case.
    `f`, `g` and `h` must be initialized with same `p`, `exp` and `n`.

.. function:: void unity_zp_sqr_inplace(unity_zp f, const unity_zp g, fmpz_t * t)

    Sets `f` to `g * g`. If `p^{exp} = 3, 4, 5, 7, 8, 9, 11, 16` special
    multiplication functions are used. Allocated array `t` of ``fmpz_t`` are
    used for all computations in this case.
    `f` and `g` must be initialized with same `p`, `exp` and `n`.



$\mathbb{Z}[\zeta_p]/(n)$. Powering functions
--------------------------------------------------------------------------------


.. function:: void unity_zp_pow_fmpz(unity_zp f, unity_zp g, const fmpz_t pow)

    Sets the `f` to `g^{pow}`. `f` and `g` must be initialized with
    same `p`, `exp` and `n`.

.. function:: void unity_zp_pow_ui(unity_zp f, unity_zp g, ulong pow)

    Sets the `f` to `g^{pow}`. `f` and `g` must be initialized with 
    same `p`, `exp` and `n`.

.. function:: ulong _unity_zp_pow_select_k(const fmpz_t n)

    Returns smallest integer `k` satisfies: 
    `\log (n) < (k(k + 1)2^{2k}) / (2^{k + 1} - k - 2) + 1`

.. function:: void unity_zp_pow_2k_fmpz(unity_zp f, unity_zp g, const fmpz_t pow)

    Sets the `f` to `g^{pow}` using `2^k`-ary exponentiation method.
    `f` and `g` must be initialized with same `p`, `exp` and `n`.

.. function:: void unity_zp_pow_2k_ui(unity_zp f, const unity_zp g, ulong pow)

    Sets the `f` to `g^{pow}` using `2^k`-ary exponentiation method.
    `f` and `g` must be initialized with same `p`, `exp` and `n`.

.. function:: void unity_zp_pow_sliding_fmpz(unity_zp f, unity_zp g, const fmpz_t pow)

    Sets the `f` to `g^{pow}` using sliding window exponentiation method.
    `f` and `g` must be initialized with same `p`, `exp` and `n`.



$\mathbb{Z}[\zeta_p]/(n)$. Cyclotomic reduction
--------------------------------------------------------------------------------


.. function:: void _unity_zp_reduce_cyclotomic_divmod(unity_zp f)

    Sets `f = \sigma_x(g)`, there automorphism `\sigma_x(\zeta)=\zeta^x`.
    `f` and `g` must be initialized with same `p`, `exp` and `n`.

.. function:: void _unity_zp_reduce_cyclotomic(unity_zp f)

    Sets `f = f \bmod \Phi_{p^{exp}}`. `\Phi_{p^{exp}}` is the `p^{exp}`-th
    cyclotomic polynomial. `g` must be reduced by `x^{p^{exp}}-1` poly.
    `f` and `g` must be initialized with same `p`, `exp` and `n`.

.. function:: void unity_zp_reduce_cyclotomic(unity_zp f, const unity_zp g)

    Sets `f = g \bmod \Phi_{p^{exp}}`. `\Phi_{p^{exp}}` is the `p^{exp}`-th
    cyclotomic polynomial.



$\mathbb{Z}[\zeta_p]/(n)$. Automorphism and inverse
--------------------------------------------------------------------------------


.. function:: void unity_zp_aut(unity_zp f, const unity_zp g, ulong x)

    Sets `f = \sigma_x(g)`, there automorphism `\sigma_x(\zeta)=\zeta^x`.
    `f` and `g` must be initialized with same `p`, `exp` and `n`.

.. function:: void unity_zp_aut_inv(unity_zp f, const unity_zp g, ulong x)

    Sets `f = \sigma_x^{-1}(g)`, so `\sigma_x(f) = g`. 
    `g` must be reduced by `\Phi_{p^{exp}}`.
    `f` and `g` must be initialized with same `p`, `exp` and `n`.



$\mathbb{Z}[\zeta_p]/(n)$. Jacobi sum
--------------------------------------------------------------------------------

In this part `\chi_{p, q}` is the character defined by 
`\chi_{p, q}(g^x) = \zeta_{p^k}^x`, there `g` is
a primitive root modulo `q`.

.. function:: void unity_zp_jacobi_sum_pq(unity_zp f, ulong q, ulong p)

    Sets `f` to Jacobi sum `J(p, q) = j(\chi_{p, q}, \chi_{p, q})`.

.. function:: void unity_zp_jacobi_sum_2q_one(unity_zp f, ulong q)

    Sets `f` to Jacobi sum 
    `J_2(q) = j(\chi_{2, q}^{2^{k - 3}}, \chi_{2, q}^{3 * 2^{k - 3}}))^2`

.. function:: void unity_zp_jacobi_sum_2q_two(unity_zp f, ulong q)

    Sets `f` to Jacobi sum
    `J_3(1) = j(\chi_{2, q}, \chi_{2, q}, \chi_{2, q}) = 
    J(2, q) * j(\chi_{2, q}^2, \chi_{2, q})`



Operations in $\mathbb{Z}[\zeta_q, \zeta_p]/(n)$.
--------------------------------------------------------------------------------

``unity_zpq`` represents as array of ``fmpz_mod_poly``.

.. function:: void unity_zpq_init(unity_zpq f, ulong q, ulong p, const fmpz_t n)

    Initialized ``unity_zpq`` element of `\mathbb{Z}[\zeta_q, \zeta_p]/(n)`.
    ``unity_zpq`` is an array of ``fmpz_mod_poly_t`` elements.

.. function:: void unity_zpq_clear(unity_zpq f)

    Clears the given ``unity_zpq`` element. It must be reinitialised in
    order to be used again.

.. function:: void unity_zpq_copy(unity_zpq f, const unity_zpq g)

    Sets `f` to `g`. `f` and `g` must be initialized with
    same `p`, `q` and `n`.

.. function:: void unity_zpq_swap(unity_zpq f, unity_zpq q)

    Swaps `f` and `g`. `f` and `g` must be initialized with
    same `p`, `q` and `n`.

.. function:: int unity_zpq_equal(const unity_zpq f, const unity_zpq g)

    Returns non zero value if `f == g`.

.. function:: ulong unity_zpq_p_unity(const unity_zpq f)

    If `f = \zeta_p^x` returns `x in [0, p - 1]`; otherwise returns `f->p`. 

.. function:: int unity_zpq_is_p_unity(const unity_zpq f)

    Returns non zero value if `f = \zeta_p^x`.

.. function:: int unity_zpq_is_p_unity_generator(const unity_zpq f)

    Returns non zero value if `f` is a generator of `<\zeta_p>` cyclic group.

.. function:: void unity_zpq_coeff_set_fmpz(unity_zpq f, ulong i, ulong j, const fmpz_t x)

    Sets the coefficient of the `\zeta_q^i*\zeta_p^j` equal to x.
    ``i`` must be less then {q} and ``j`` must be less then {p}.

.. function:: void unity_zpq_coeff_set_ui(unity_zpq f, ulong i, ulong j, ulong x)

    Sets the coefficient of the `\zeta_q^i*\zeta_p^j` equal to x.
    ``i`` must be less then {q} and ``j`` must be less then {p}.

.. function:: void unity_zpq_coeff_add(unity_zpq f, ulong i, ulong j, const fmpz_t x)

    Sets `a` in `a*\zeta_p^i*\zeta_q^j` to `a + x`. `x` must be less then `n`.

.. function:: void unity_zpq_add(unity_zpq f, const unity_zpq g, const unity_zpq h)

    Sets the ``f`` to the ``g``+``h``.
    ``f``, ``g`` and ``h`` must be initialized with same
    ``q``, ``p`` and ``n``.

.. function:: void unity_zpq_mul(unity_zpq f, const unity_zpq g, const unity_zpq h)

    Sets the ``f`` to the ``g``*``h``.
    ``f``, ``g`` and ``h`` must be initialized with same
    ``q``, ``p`` and ``n``.

.. function:: void _unity_zpq_mul_unity_p(unity_zpq f)

    Sets `f = f * \zeta_p`.

.. function:: void unity_zpq_mul_unity_p_pow(unity_zpq f, const unity_zpq g, ulong k)

    Sets `f` to `g * \zeta_p^k`.

.. function:: void unity_zpq_pow(unity_zpq f, unity_zpq g, const fmpz_t p)

    Sets the `f` to `g^p`. `f` and `g` must be initialized with same `p`, `q` and `n`.

.. function:: void unity_zpq_pow_ui(unity_zpq f, unity_zpq g, ulong p)

    Sets the `f` to `g^p`. `f` and `g` must be initialized with same `p`, `q` and `n`.

.. function:: void unity_zpq_gauss_sum(unity_zpq f, ulong q, ulong p)

    Sets `f = \tau(\chi_{p, q})`.

.. function:: void unity_zpq_gauss_sum_sigma_pow(unity_zpq f, ulong q, ulong p)

    Sets `f = \tau^{\sigma_n}(\chi_{p, q})`.
