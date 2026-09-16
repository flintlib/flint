.. _ecpp:

**ecpp.h** -- elliptic curve primality proving
===============================================================================

Elliptic curve primality proving (ECPP) in the Atkin-Morain form. Unlike
APRCL, the proof produces a certificate that can be verified independently
and much faster than it was produced, and the running time scales better
with the size of the input.

A certificate is a chain of steps `(n_i, a_i, b_i, m_i, q_i, P_i)`. Step `i`
shows that `n_i` is prime provided `q_i` is prime: the point `P_i` on the
curve `y^2 = x^3 + a_i x + b_i` modulo `n_i` satisfies `m_i P_i = O` and
`(m_i / q_i) P_i \neq O` with `q_i \mid m_i` and `q_i > (n_i^{1/4} + 1)^2`.
The chain continues with `n_{i+1} = q_i` until `q < 2^{64}`, which is
checked directly.

The prover uses the organisation of the discriminant search of fastECPP
(Franke, Kleinjung, Morain, Wirth). A pool of fundamental discriminants
`D < 0` that are smooth over a fixed set of small primes (the set, the
class number bound and the trial division bound depend on the size of the
input) is built once, sorted by the cost of realising a step with `D`:
number of units, class number `h(D)`, odd part of `h(D)`. For each `n`
the pool is scanned in that order. A discriminant is only used when all
its genus characters `(p^*/n)` are 1, since otherwise `n` cannot be
represented by the principal form; then `\sqrt{D} \bmod n` is the product
of the square roots of the `p^* \mid D`, which are computed once per `n`,
and Cornacchia writes `4n = t^2 + |D| v^2`. The possible curve orders
`n + 1 \pm t` (and the extra twists for `D = -3, -4`) are collected in
batches; one remainder tree removes all their prime factors below the
trial division bound, and the cofactors `q` with `(n^{1/4} + 1)^2 < q
\le n/2` are tested for probable primality. The candidate with the best
estimated cost per bit gained (the cost being that of finding a root of
the factor of `H_D` over the genus field, see ``ecpp_class_poly_genus``,
of degree `h / 2^{g-1}` where `g` is the number of genus characters,
by radicals up to degree 4) is realised: a root modulo `n` of the
Hilbert class polynomial `H_D` gives the curve, the right twist and a point
are found by testing, and the prover recurses on `q`. Dead ends continue
with the next candidate of the batch, then with the next batch.

``fmpz_is_prime`` uses ECPP instead of APRCL for its final proof from a
size that depends on the number of available threads: 250 bits with one
thread, 550 with two, 1000 with three, 2000 with four and 8000 with up
to eight (APRCL parallelises better than the largely serial ECPP chain;
these are crossovers measured on random primes). With threads, the
square roots, the Cornacchia steps, the probable prime tests and the
reduction of the primorial are distributed, the two twists of a curve are
tried at once, and the realisation of a step (class polynomial, root,
curve and point) runs on a pool thread while the search for the next
step proceeds.

Prior art. The organisation of the search follows the fastECPP of
Franke, Kleinjung, Morain and Wirth [FKMW2004]_ and Morain [Morain2007]_,
in the form of the FastECPP implementations by Jared Asuncion in PARI/GP
(``primecert``, [PARI]_) and by Andreas Enge in the CM library
[EngeCM]_, [Enge2024]_: the discriminants smooth over a set of small
primes, the square roots of the primes computed once per `n`, the batch
factoring of the curve orders with a remainder tree and the choice of
the trial division bound and of the minimum gain per step; the rounds
of the search extending the prime set until a few prime cofactors are
expected, the trial division bound fixed for the whole chain, the
cache-blocked class number sieve and the class field tower with the
Hecke representation of its coefficients are as in CM, which also
provided the reference for the normalisations of Weber's function as a
class invariant ([Schertz2002]_, [YuiZagier1997]_). The class field tower
is the decomposition of Enge and Morain [EngeMorain2003]_; the factor of
`H_D` over the genus field and the descent by radicals and by Kummer
theory through the tower are original to this implementation. The
timings of both PARI/GP and CM on the same inputs were used throughout
to tune the parameters.

``ecpp_set_verbose(1)`` reports the steps of a proof as they are found;
compiling ``prove.c`` with ``-DECPP_PROFILE`` prints a breakdown of the
time at the end, and ``profile/p-prove.c`` times ECPP against APRCL.

Example. ``build/examples/factor_integer -certify n`` proves each prime
factor above 64 bits with ECPP and prints the certificate; ``-pari``
prints it in the ``primecert`` syntax of PARI/GP instead, ``-verbose``
reports the steps::

    $ build/examples/factor_integer -certify "2^128+1"
    340282366920938463463374607431768211457 =
    59649589127497217 * 5704689200685129054721

    ECPP certificate for 5704689200685129054721 (verified):
    ecpp certificate, 1 steps
    [0] n = 5704689200685129054721
        D = -40
        a = 580466530013776483294
        b = 4831761323372247440219
        m = 5704689200761597106060
        q = 636165875659
        x = 2451296810544034633782
        y = 5469387309100357564736

    $ build/examples/factor_integer -pari "2^128+1"
    ...
    [[5704689200685129054721, -76468051338, 8967298340, 580466530013776483294, [2451296810544034633782, 5469387309100357564736]]]

The step reads: the curve `y^2 = x^3 + a x + b` over `\mathbb{Z}/n`
has `m = n + 1 - t = s q` points with `q = 636165875659` prime and
`(n^{1/4} + 1)^2 < q`, and `P = (x, y)` satisfies `m P = O`,
`(m/q) P \ne O`; hence `n` is prime. In the PARI/GP form the entries
are `[n, t, s, a, [x, y]]`, with `b` implied by the point, and each
step's `q = (n + 1 - t)/s` is the `n` of the next one (the last `q` is
a word-size prime). ``primecertisvalid`` in PARI/GP accepts the output
of ``-pari``.

Point arithmetic modulo `n` is done in Jacobian coordinates; every quantity
whose invertibility a formula assumes is accumulated and checked with a
single gcd at the end, so that a composite `n` is detected rather than
producing a wrong proof.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: ecpp_step_struct

    One step of a certificate: `n`, `D`, curve `(a, b)`, `m`, `q` and the
    point `(x, y)`.

.. type:: ecpp_cert_struct

.. type:: ecpp_cert_t

    A certificate: an array of steps.

.. type:: ecpp_point_struct

.. type:: ecpp_point_t

    A point `(X : Y : Z)` in Jacobian coordinates; `Z = 0` is the point at
    infinity.

Certificates
-------------------------------------------------------------------------------

.. function:: void ecpp_cert_init(ecpp_cert_t cert)
              void ecpp_cert_clear(ecpp_cert_t cert)

.. function:: ecpp_step_struct * ecpp_cert_push(ecpp_cert_t cert)
              void ecpp_cert_pop(ecpp_cert_t cert)

    Appends an initialised step, respectively removes the last step.

.. macro:: ECPP_CERT_FORMAT_FLINT
           ECPP_CERT_FORMAT_PARI

    Text formats of a certificate: the labelled form shown above, and the
    ``primecert`` syntax of PARI/GP, ``[[n, t, s, a, [x, y]], ...]`` with
    `n + 1 - t = s q`.

.. function:: void ecpp_cert_print(const ecpp_cert_t cert, int format)
              void ecpp_cert_fprint(FILE * file, const ecpp_cert_t cert, int format)
              char * ecpp_cert_get_str(const ecpp_cert_t cert, int format)

    Prints, respectively returns as a string (to be freed with
    ``flint_free``), the certificate in the given format.

.. function:: int ecpp_cert_set_str(ecpp_cert_t cert, const char * str)

    Reads a certificate in either format (recognised by the first
    character). Returns 1 on success, 0 (with an empty certificate) if
    the string is malformed. In the PARI/GP format `b` is recovered from
    the point and the discriminant is not available (set to 0); the
    certificate verifies all the same.

.. function:: void ecpp_set_verbose(int verbose)

    With ``verbose`` nonzero, ``ecpp_prove`` prints each step as it is
    found (the sizes, the discriminant, the class number and the degree
    of the polynomial whose root is needed).

Proving and verifying
-------------------------------------------------------------------------------

.. function:: int ecpp_prove(ecpp_cert_t cert, const fmpz_t n)

    Attempts to prove that `n` is prime. Returns 1 and fills ``cert`` with
    a certificate if `n` is prime, 0 if `n` is composite, and -1 if no proof
    could be found. For `n < 2^{64}` the certificate is empty. The input
    is expected to be a probable prime; composites are usually detected
    quickly by the initial BPSW test.

.. function:: int ecpp_verify(const ecpp_cert_t cert, const fmpz_t n)

    Returns 1 if ``cert`` is a valid certificate proving that `n` is prime,
    0 otherwise. An empty certificate is valid exactly for `n < 2^{64}`
    prime.

.. function:: int ecpp_verify_step(const ecpp_step_struct * s)

    Checks one step, i.e. that `n` is prime if `q` is prime.

.. function:: int ecpp_is_prime(const fmpz_t n)

    Runs ``ecpp_prove`` and discards the certificate.

Curve arithmetic
-------------------------------------------------------------------------------

.. function:: void ecpp_point_init(ecpp_point_t P)
              void ecpp_point_clear(ecpp_point_t P)
              void ecpp_point_set_affine(ecpp_point_t P, const fmpz_t x, const fmpz_t y)
              int ecpp_point_is_zero(const ecpp_point_t P)

.. function:: void ecpp_gr_ctx_init(gr_ctx_t gctx, const fmpz_mod_ctx_t ctx)

    Initialises ``gctx`` as the ring of integers modulo the modulus of
    ``ctx``: ``mpn_mod`` for moduli of 2 to ``MPN_MOD_MAX_LIMBS`` limbs,
    ``fmpz_mod`` otherwise.

.. function:: int ecpp_point_mul_gr(gr_ptr R, gr_srcptr P, const fmpz_t k, gr_srcptr a, gr_ptr acc, gr_ctx_t ctx)

    The scalar multiplication over a generic ring (a point is three
    consecutive elements X, Y, Z), the implementation behind
    ``ecpp_point_mul``, which converts to the ring given by
    ``ecpp_gr_ctx_init``; about 1.7 times faster than the ``fmpz_mod``
    arithmetic at 200 bits and 1.6 at 400.

.. function:: int ecpp_point_mul(ecpp_point_t R, const ecpp_point_t P, const fmpz_t k, const fmpz_t a, fmpz_t acc, const fmpz_mod_ctx_t ctx)

    Sets `R = k P` on `y^2 = x^3 + a x + b` modulo the modulus of ``ctx``,
    for an affine point `P` (`Z = 1`) and `k \ge 0`. Every quantity whose
    invertibility the formulas assume is multiplied into ``acc``; the
    result is only meaningful if `\gcd(\mathrm{acc}, n) = 1` at the end,
    otherwise `n` is composite. Returns 0 if ``acc`` became zero.

Discriminants and Cornacchia
-------------------------------------------------------------------------------

.. function:: slong ecpp_disc_table(ecpp_disc_struct ** table, const ulong * primes, slong nprimes, slong Dmax, slong hmax, slong omax, double costmax)

    Fills ``*table`` (allocated with ``flint_malloc``) with the fundamental
    discriminants `-D_{\max} \le D < 0` of class number at most ``hmax``
    whose odd prime factors are among the ``nprimes`` odd primes in
    ``primes``, sorted by number of units, class number, odd part of the
    class number and `|D|`, and returns their number. Each entry records
    the class number, the number `g` of genus characters, the odd part
    `o = h / 2^{g-1}`, the factor `q_0 \in \{1, -4, 8, -8\}` of `D` at 2
    and the indices of the odd primes dividing `D`.

.. function:: int ecpp_root_radicals(fmpz_t x, const fmpz_mod_poly_t f, flint_rand_t state, const fmpz_mod_ctx_t ctx)

    Sets ``x`` to a root of the monic polynomial ``f`` of degree 2, 3 or 4,
    assumed to split completely modulo the prime modulus `n` of ``ctx``,
    by radicals (Cardano, Ferrari): a few square and cube roots modulo
    `n`, some ten times cheaper than a powering of a polynomial modulo
    ``f``. Degrees 3 and 4 require `n \equiv 1 \pmod 3`. Returns 1 on
    success and 0 if the degree or `n` is not handled or the computation
    failed.

.. function:: int ecpp_poly_root(fmpz_t x, const fmpz_mod_poly_t f, flint_rand_t state, const fmpz_mod_ctx_t ctx)

    Sets ``x`` to a root of the monic polynomial ``f``, assumed to split
    completely modulo the prime modulus of ``ctx``: by radicals when
    possible, else by random splitting with `(x + r)^{(n-1)/2}`. Returns
    1 on success, 0 otherwise.

.. function:: int ecpp_class_poly_tower(fmpz_t j, slong D, flint_rand_t state, const fmpz_mod_ctx_t ctx)

    Sets ``j`` to a root modulo the prime modulus `n` of ``ctx`` of the
    Hilbert class polynomial `H_D`, for a fundamental discriminant `D`
    such that `n` splits completely in the Hilbert class field (`4n = t^2
    + |D| v^2`), through the class field tower (Enge, Morain): a
    composition series of the class group with prime quotients `p_i` gives
    a tower of fields with relative polynomials whose coefficients are
    stored in Hecke form (numerators over the derivative of the polynomial
    of the level below, so that all coefficients are integers), computed
    numerically from the `j`-values ordered by cosets, and modulo `n` a
    root is found level by level in degrees `p_i` instead of `h`. Returns
    1 on success, 0 on failure. For even `D` the values of Weber's
    function `f` (a class invariant of height `72/e` times smaller than
    `j`'s, with the normalisations of Enge's CM library) are used instead
    of `j`, and `j` is recovered from the root modulo `n` by a closed
    formula. The Lagrange numerators are built by a subproduct tree, so
    that about 1.5 times the height of the class polynomial suffices:
    `h = 256` takes 0.13 s with the Weber invariant and 3.4 s with `j`.
    For odd `D` the caller may pass `4D`: the order of conductor 2 has the
    same class field (for `D \equiv 1 \pmod 8`, or three times the class
    number for `D \equiv 5 \pmod 8`) and admits the Weber invariant; its
    curves have the right order when `v` is even. Levels of prime degree
    `p \ge 5` also carry Kummer data (the `p`-th power of the Lagrange
    resolvent and the ratios of the resolvents, as elements of the level
    below tensored with `\mathbb{Q}(\zeta_p)`), so that when `p \mid n - 1`
    the level is descended with one `p`-th root instead of a root of a
    degree-`p` polynomial.

.. function:: int ecpp_class_poly_genus(fmpz_mod_poly_t F, slong D, const slong * pstar, slong g, const fmpz * sqrts, const fmpz_mod_ctx_t ctx)

    Sets ``F`` to the factor of the Hilbert class polynomial `H_D` modulo
    `n` (the modulus of ``ctx``) corresponding to the principal genus. Here
    `D = p_1^* \cdots p_g^*` with the `p_i^*` given in ``pstar`` (the
    signed odd primes `p^* = (-1)^{(p-1)/2} p` dividing `D` and, for even
    `D`, the factor `q_0 \in \{-4, 8, -8\}`), and ``sqrts[i]`` is a square
    root of ``pstar[i]`` modulo `n`, which must exist. `H_D` factors over
    the genus field `K(\sqrt{p_1^*}, \ldots, \sqrt{p_g^*})` into `2^{g-1}`
    conjugate factors of degree `h(D) / 2^{g-1}`, whose coefficients are
    identified numerically from the values of `j` at the reduced forms
    grouped by genus, and ``F`` is the image of one of them modulo `n`. Its
    roots are `j`-invariants of curves with CM by the maximal order of
    `\mathbb{Q}(\sqrt{D})`, and finding one costs a fraction `2^{-(g-1)}`
    of a root of `H_D`. Returns 1 on success, 0 if `g < 2` or the
    coefficients could not be identified. Requires `D` fundamental.

.. function:: int ecpp_cornacchia(fmpz_t t, fmpz_t v, const fmpz_t n, slong D, const fmpz_t sqrtD)

    Given `D < 0`, `n` odd with `\gcd(n, D) = 1` and a square root of `D`
    modulo `n`, finds `t, v \ge 0` with `t^2 + |D| v^2 = 4n` and returns 1,
    or returns 0 if there is no solution.
