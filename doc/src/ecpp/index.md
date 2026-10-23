---
title: ecpp.h
header: ecpp.h
description: Elliptic curve primality proving with verifiable certificates.
category: Integers
order: 10
---

This module implements elliptic curve primality proving (ECPP) in the
Atkin–Morain form [@AtkMor1993], with the organisation of the
discriminant search of fastECPP [@FKMW2004; @Morain2007]. Unlike APRCL,
ECPP produces a *certificate*: a chain of elliptic curves and points from
which anybody can check the primality of $n$ much faster than it was
proved, with nothing more than arithmetic modulo the numbers in the
chain. Its running time also grows more slowly with the size of $n$
(heuristically $\tilde O(\log^4 n)$ against $(\log n)^{O(\log \log \log n)}$
for APRCL). `fmpz_is_prime` uses ECPP for its final proof from about 250
bits (more with several threads, see below).

[](ecpp/prove) shows how to prove, verify, print and read
certificates. The mathematics (the theorem behind a step, the construction
of curves by complex multiplication and the choice of discriminants) is
described in [](algorithms/number-theory/primality-proving).

::: summary algorithms/number-theory/primality-proving :::

## Proving and verifying {#proving}

::: function
int ecpp_prove(ecpp_cert_t cert, const fmpz_t n)

Attempts to prove that $n$ is prime. Returns 1 and fills `cert` with a
certificate if $n$ is prime, 0 if $n$ is composite (or $n \le 1$), and
$-1$ if no proof could be found. For $n < 2^{64}$ the certificate is
empty and primality is decided by `fmpz_is_prime`. The input is first
subjected to a BPSW test, so composites are detected quickly.
:::

::: function
int ecpp_verify(const ecpp_cert_t cert, const fmpz_t n)

Returns 1 if `cert` is a valid certificate proving that $n$ is prime, and
0 otherwise. An empty certificate is valid exactly for primes
$n < 2^{64}$.
:::

::: function
int ecpp_verify_step(const ecpp_step_struct * s)

Checks one step: returns 1 if it proves that its $n$ is prime provided
that its $q$ is prime, and 0 otherwise.
:::

::: function
int ecpp_is_prime(const fmpz_t n)

Runs `ecpp_prove` and discards the certificate; same return values.
:::

::: function
void ecpp_set_verbose(int verbose)

With `verbose` nonzero, `ecpp_prove` prints each step as it is found
(the sizes, the discriminant, the class number and the degree of the
polynomial whose root is needed). The setting is global.
:::

## Certificates {#certificates}

A certificate is a chain of steps $(n_i, a_i, b_i, m_i, q_i, P_i)$ with
$n_0 = n$. Step $i$ shows that $n_i$ is prime provided that $q_i$ is
prime: the point $P_i$ on the curve $y^2 = x^3 + a_i x + b_i$ modulo
$n_i$ satisfies $m_i P_i = O$ and $(m_i/q_i) P_i \ne O$, with
$q_i \mid m_i$ and $q_i > (n_i^{1/4} + 1)^2$. The chain continues with
$n_{i+1} = q_i$ until $q < 2^{64}$, whose primality is checked directly.

::: type
ecpp_step_struct

One step of a certificate, with the fields `n`, `D` (the discriminant
used, for information; 0 if unknown), `a`, `b`, `m`, `q`, `x` and `y`
(the point $P$), all of type `fmpz_t` except `D`.
:::

::: type
ecpp_cert_struct
ecpp_cert_t

A certificate: an array `steps` of `num` steps (with `alloc` allocated).
:::

::: function
void ecpp_cert_init(ecpp_cert_t cert)
void ecpp_cert_clear(ecpp_cert_t cert)

Initialise an empty certificate, and free it.
:::

::: function
ecpp_step_struct * ecpp_cert_push(ecpp_cert_t cert)
void ecpp_cert_pop(ecpp_cert_t cert)

Append a new step (with initialised fields) and return a pointer to it,
respectively remove the last step.
:::

::: macro
ECPP_CERT_FORMAT_FLINT
ECPP_CERT_FORMAT_PARI

The text formats of a certificate: a labelled form listing $n$, $D$, $a$,
$b$, $m$, $q$, $x$ and $y$ for each step, and the `primecert` syntax of
PARI/GP [@PARI], `[[n, t, s, a, [x, y]], ...]`, where $n + 1 - t = sq$
and $b$ is implied by the point. `primecertisvalid` in PARI/GP accepts
the second form.
:::

::: function
void ecpp_cert_print(const ecpp_cert_t cert, int format)
void ecpp_cert_fprint(FILE * file, const ecpp_cert_t cert, int format)
char * ecpp_cert_get_str(const ecpp_cert_t cert, int format)

Print the certificate in the given format to `stdout` or to `file`, or
return it as a string to be freed with `flint_free`.
:::

::: function
int ecpp_cert_set_str(ecpp_cert_t cert, const char * str)

Reads a certificate in either format (recognised by its first
character). Returns 1 on success, and 0, leaving an empty certificate, if
the string is malformed. In the PARI/GP format, $b$ is recovered from the
point and the discriminant is not available (it is set to 0); the
certificate verifies all the same.
:::

## Curve arithmetic {#curves}

Points on $y^2 = x^3 + ax + b$ modulo $n$ are kept in Jacobian
coordinates. Since $n$ is not known to be prime during the proof, every
quantity whose invertibility the formulas assume is multiplied into an
accumulator, which is checked with a single GCD at the end: a composite
$n$ is then detected instead of producing a wrong proof.

::: type
ecpp_point_struct
ecpp_point_t

A point $(X : Y : Z)$ in Jacobian coordinates (fields `X`, `Y`, `Z`);
$Z = 0$ is the point at infinity.
:::

::: function
void ecpp_point_init(ecpp_point_t P)
void ecpp_point_clear(ecpp_point_t P)
void ecpp_point_set_affine(ecpp_point_t P, const fmpz_t x, const fmpz_t y)
int ecpp_point_is_zero(const ecpp_point_t P)

Initialise and free a point, set it to the affine point $(x, y)$ (with
$Z = 1$), and test whether it is the point at infinity.
:::

::: function
int ecpp_point_mul(ecpp_point_t R, const ecpp_point_t P, const fmpz_t k, const fmpz_t a, fmpz_t acc, const fmpz_mod_ctx_t ctx)

Sets $R = kP$ on the curve with coefficient $a$ modulo the modulus $n$ of
`ctx`, for an affine point $P$ and $k \ge 0$. Every quantity whose
invertibility the formulas assume is multiplied into `acc`; the result is
meaningful only if $\gcd(\mathit{acc}, n) = 1$ at the end, and otherwise
$n$ is composite. Returns 0 if `acc` became zero and 1 otherwise.
:::

::: function
void ecpp_gr_ctx_init(gr_ctx_t gctx, const fmpz_mod_ctx_t ctx)

Initialises `gctx` as the ring of integers modulo the modulus of `ctx`:
`mpn_mod` for moduli of 2 to `MPN_MOD_MAX_LIMBS` limbs, and `fmpz_mod`
otherwise.
:::

::: function
int ecpp_point_mul_gr(gr_ptr R, gr_srcptr P, const fmpz_t k, gr_srcptr a, gr_ptr acc, gr_ctx_t ctx)

The scalar multiplication over a generic ring (a point being three
consecutive elements $X$, $Y$, $Z$), on which `ecpp_point_mul` is
built after converting to the ring given by `ecpp_gr_ctx_init`.
:::

## Discriminants {#discriminants}

The prover uses a pool of fundamental discriminants $D < 0$ whose odd
prime factors lie in a fixed set of small primes, so that $\sqrt D \bmod n$
is a product of square roots computed once per $n$. They are ordered by
the estimated cost of realising a step with them.

::: type
ecpp_disc_struct

A discriminant of the pool: `D`, the class number `h`, the number `g` of
genus characters, the odd part `o` $= h/2^{g-1}$ of the class number,
the factor `q0` $\in \{1, -4, 8, -8\}$ of $D$ at 2 (so that
$D = q_0 \prod p^*$), the number `nfac` of odd primes dividing $D$ with
their indices `fac` in the prime list, and the estimated cost `cost`.
:::

::: macro
ECPP_DISC_MAXFAC
ECPP_DISC_COST_MAX

The maximum number (8) of odd primes dividing a discriminant of the pool,
and the maximum estimated cost (20) of a discriminant admitted to it.
:::

::: function
slong ecpp_disc_table(ecpp_disc_struct ** table, const ulong * primes, slong nprimes, slong Dmax, slong hmax, slong omax, double costmax)

Sets `*table` (allocated with `flint_malloc`) to the fundamental
discriminants $-D_{\max} \le D < 0$ of class number at most `hmax`, odd
part of the class number at most `omax` and cost at most `costmax`, whose
odd prime factors are among the `nprimes` odd primes in `primes`, and
returns their number. They are sorted by number of units, odd part of the
class number, cost and $|D|$. The class numbers are computed with a
cache-blocked sieve over reduced forms.
:::

::: function
double ecpp_disc_cost(const ecpp_disc_struct * d)
double ecpp_disc_cost_v(const ecpp_disc_struct * d, int veven)
double ecpp_disc_cost_n(const ecpp_disc_struct * d, int veven, const ulong * nmodp, slong bits)

The estimated cost of realising a step with the discriminant `d`, in
units of one scalar multiplication on the curve at a few thousand bits:
the twists, the numerical computation of the class polynomials of the
class field tower (cheap with Weber's function, which is usable for even
$D$ and, through the order of conductor 2, for odd $D$ when the $v$ of
$4n = t^2 + |D| v^2$ is even, as `veven` says), and the descent through
the prime factors $p$ of $h$ (a square root for $p = 2$, radicals or a
small powering for $p = 3$, a $p$-th root or a powering of a degree $p$
polynomial otherwise), or the cost of a root of the genus field factor if
that is smaller. `ecpp_disc_cost_n` uses $n$ modulo 5, 7, 11 and 13
(`nmodp`, or the average over $n$ if `NULL`) and scales the cost of the
class polynomial numerics, which does not depend on $n$, to the bit size
of $n$. `ecpp_disc_cost` is the optimistic estimate (with `veven` = 1)
used to order and admit the pool. Class numbers with a prime factor
above 13 get a cost above `ECPP_DISC_COST_MAX`.
:::

::: function
int ecpp_disc_use_tower(const ecpp_disc_struct * d, int veven, slong * Dt)

Returns whether the class field tower is the cheaper way to a curve for
`d` (class number at least 2, and Weber's function usable or $h \le 48$),
and sets `Dt` to the discriminant to build the tower on ($D$, or $4D$
for odd $D$ with `veven` set).
:::

::: function
int ecpp_cornacchia(fmpz_t t, fmpz_t v, const fmpz_t n, slong D, const fmpz_t sqrtD)

Given $D < 0$, an odd $n$ with $\gcd(n, D) = 1$ and a square root of $D$
modulo $n$, finds $t, v \ge 0$ with $t^2 + |D| v^2 = 4n$ and returns 1,
or returns 0 if there is no solution (Cornacchia's algorithm in the
form of Cohen [@Coh1996, Algorithm 1.5.3]).
:::

## Class polynomials and roots {#class-polynomials}

::: function
int ecpp_class_poly_tower(fmpz_t j, slong D, flint_rand_t state, const fmpz_mod_ctx_t ctx)
int _ecpp_class_poly_tower(fmpz_t j, slong D, int flags, flint_rand_t state, const fmpz_mod_ctx_t ctx)

Set `j` to a root modulo the prime modulus $n$ of `ctx` of the Hilbert
class polynomial $H_D$, for a fundamental discriminant $D$ such that $n$
splits completely in the Hilbert class field (that is,
$4n = t^2 + |D| v^2$), through the class field tower of Enge and Morain
[@EngeMorain2003]. Return 1 on success and 0 on failure. The flags of
the underscore version are `ECPP_TOWER_KUMMER` (compute the Kummer data,
which `ecpp_class_poly_tower` does from `ECPP_KUMMER_BITS` bits) and
`ECPP_TOWER_J` (use $j$ rather than Weber's function).

::: details
A composition series of the class group with prime quotients $p_i$ gives
a tower of fields with relative polynomials, whose coefficients are
stored in Hecke form (numerators over the derivative of the polynomial of
the level below, so that all coefficients are integers) and computed
numerically from the $j$-values ordered by cosets. Modulo $n$, a root is
found level by level in degrees $p_i$ instead of $h$. For even $D$, the
values of Weber's function $\mathfrak f$ (a class invariant whose class
polynomial has height $72/e$ times smaller than that of $j$, with the
normalisations of Enge's CM library [@EngeCM; @Schertz2002;
@YuiZagier1997]) are used instead of $j$, and $j$ is recovered from the
root modulo $n$ by a closed formula. The Lagrange numerators are built by
a subproduct tree, so that about 1.5 times the height of the class
polynomial suffices as precision. For odd $D$ the caller may pass $4D$: the
order of conductor 2 has the same class field (for $D \equiv 1 \pmod 8$;
three times the class number for $D \equiv 5 \pmod 8$) and admits
Weber's function, and its curves have the right order when $v$ is even.
Levels of prime degree $p \ge 5$ also carry Kummer data (the $p$-th power
of the Lagrange resolvent and the ratios of the resolvents, as elements
of the level below tensored with $\mathbb{Q}(\zeta_p)$), so that when
$p \mid n \pm 1$ the level is descended with one $p$-th root instead of a
root of a polynomial of degree $p$.
:::
:::

::: macro
ECPP_TOWER_KUMMER
ECPP_TOWER_J

Flags of `_ecpp_class_poly_tower` (1 and 2).
:::

::: function
int ecpp_class_poly_genus(fmpz_mod_poly_t F, slong D, const slong * pstar, slong g, const fmpz * sqrts, const fmpz_mod_ctx_t ctx)

Sets `F` to the factor of the Hilbert class polynomial $H_D$ modulo $n$
(the modulus of `ctx`) that corresponds to the principal genus. Here
$D = p_1^* \cdots p_g^*$ with the $p_i^*$ given in `pstar` (the signed odd
primes $p^* = (-1)^{(p-1)/2} p$ dividing $D$ and, for even $D$, the
factor $q_0 \in \{-4, 8, -8\}$), and `sqrts[i]` is a square root of
`pstar[i]` modulo $n$, which must exist. Over the genus field
$K(\sqrt{p_1^*}, \ldots, \sqrt{p_g^*})$, $H_D$ factors into $2^{g-1}$
conjugate factors of degree $h(D)/2^{g-1}$, whose coefficients are
identified numerically from the values of $j$ at the reduced forms
grouped by genus; `F` is the image of one of them modulo $n$, so that
finding a root costs a fraction $2^{-(g-1)}$ of a root of $H_D$. Returns 1
on success and 0 if $g < 2$ or the coefficients could not be identified.
Requires $D$ fundamental.
:::

::: function
int ecpp_root_radicals(fmpz_t x, const fmpz_mod_poly_t f, flint_rand_t state, const fmpz_mod_ctx_t ctx)

Sets `x` to a root of the monic polynomial `f` of degree 2, 3 or 4,
assumed to split completely modulo the prime modulus $n$ of `ctx`, by
radicals (Cardano, Ferrari): a few square and cube roots modulo $n$, some
ten times cheaper than a powering of a polynomial modulo `f`. Degrees 3
and 4 require $n \equiv 1 \pmod 3$. Returns 1 on success and 0 if the
degree or $n$ is not handled or the computation failed.
:::

::: function
int ecpp_poly_root(fmpz_t x, const fmpz_mod_poly_t f, flint_rand_t state, const fmpz_mod_ctx_t ctx)

Sets `x` to a root of the monic polynomial `f`, assumed to split
completely modulo the prime modulus of `ctx`: by radicals when possible,
and otherwise by random splitting with $\gcd(f, (x + r)^{(n-1)/2} - 1)$.
Returns 1 on success and 0 otherwise.
:::

## Tuning parameters {#tuning}

::: macro
ECPP_MIN_PRIME
ECPP_BIGPOOL_BITS
ECPP_POOL_HMAX
ECPP_POOL_OMAX
ECPP_DMAX_BITS(bits)
ECPP_PRIMORIAL_CHUNKS
ECPP_KUMMER_BITS

Parameters of the prover: the number of prime cofactors (4) expected per
round of the discriminant search; the size (2500 bits) from which the
pool admits all discriminants with a cheap class field tower rather than
only those with a small genus factor; the bounds on the class number
(400) and its odd part (16) in the pool; the bound $2^{20}$ to $2^{23}$
on $|D|$, depending on the bit size; the number of chunks (8) in which the
primorial for the batch factoring is reduced when threads are used; and
the size (1500 bits) from which the class field tower carries Kummer
data.
:::
