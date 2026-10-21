---
title: nmod.h
header: nmod.h
description: Arithmetic in $\Z/n\Z$ for a word-size modulus $n$, with Montgomery arithmetic, square roots and discrete logarithms.
extra-headers: gr
category: Integers mod n
order: 1
---

This module implements arithmetic in $\Z/n\Z$ for a word-size modulus
$n \ge 1$. Elements are plain `ulong` values in $[0, n)$, and the modulus is
passed as an `nmod_t` holding $n$ together with a precomputed inverse.
Most functions are inline wrappers around the `_preinv` functions of
`ulong_extras`. They are the scalar building blocks of `nmod_vec`,
`nmod_poly`, `nmod_mat` and the other modules over $\Z/n\Z$. The module
also provides square roots modulo primes, discrete logarithms, Montgomery
arithmetic, and the generic ring (`gr`) contexts for $\Z/n\Z$.

The guides below explain how to use the module; start with
[](nmod/arithmetic). The conventions common to all FLINT modules, including the
naming of functions and the use of `nmod_t` as a context, are described
in [](guide/conventions/types-and-naming).

## Conventions {#conventions}

Throughout, $B = 2^{\bits}$.

::: type
nmod_t

A modulus $n$ with precomputed data, set up by `nmod_init`. The fields are
`n` (the modulus), `norm` (the number of leading zero bits of $n$) and
`ninv` (the value `n_preinvert_limb(n)` used for
[Möller–Granlund division](algorithms/arithmetic/moller-granlund)).
An `nmod_t` is a small struct and is passed by value.
:::

- *Reduced* means $0 \le a < n$. Unless stated otherwise, operands must be
  reduced and results are reduced. Functions named `nmod_set_*` and
  `nmod_ui_*` accept any input.
- None of these assumptions is checked, except by `FLINT_ASSERT` in builds
  configured with `--enable-assert`.
- A leading underscore marks a function with an extra restriction on the
  modulus or the input (for example `_nmod_add`), or an internal helper.

::: summary algorithms/arithmetic/choosing :::

## Initialisation and conversion {#init}

::: function
void nmod_init(nmod_t * mod, ulong n)

Sets `mod` to the modulus $n$, computing `norm` and `ninv`. Requires
$n \ne 0$.
:::

::: macro
NMOD_BITS(mod)

The number of bits of `mod.n`, that is, `FLINT_BITS - mod.norm`.
:::

::: function
ulong nmod_set_ui(ulong x, nmod_t mod)
ulong nmod_set_si(slong x, nmod_t mod)

Return $x \bmod n$ in $[0, n)$, for any $x$.
:::

## Arithmetic {#arithmetic}

Addition, multiplication, inversion and division of residues. The guide
[](nmod/arithmetic) explains which of them need reduced operands and what
to do when the divisor is not invertible.

::: function
ulong nmod_add(ulong a, ulong b, nmod_t mod)
ulong nmod_sub(ulong a, ulong b, nmod_t mod)
ulong nmod_neg(ulong a, nmod_t mod)

Return $a + b$, $a - b$ and $-a$ modulo $n$, for reduced operands and any
$n$.
:::

::: function
ulong _nmod_add(ulong a, ulong b, nmod_t mod)
ulong _nmod_sub(ulong a, ulong b, nmod_t mod)

Branch-free versions of `nmod_add` and `nmod_sub`. Require
$n \le 2^{\bits - 1}$; for larger $n$ the result is wrong.
:::

::: function
ulong nmod_ui_add_ui(ulong a, ulong b, nmod_t mod)

Returns $(a + b) \bmod n$ for any $a$ and $b$.
:::

::: function
ulong nmod_mul(ulong a, ulong b, nmod_t mod)

Returns $ab \bmod n$ for reduced operands and any $n$.
:::

::: function
ulong _nmod_mul_fullword(ulong a, ulong b, nmod_t mod)

Returns $ab \bmod n$ for reduced operands. Requires `mod.norm` $= 0$, that
is, $n \ge 2^{\bits - 1}$. It saves the normalisation shift of `nmod_mul`.
:::

::: function
ulong nmod_ui_mul_ui(ulong a, ulong b, nmod_t mod)

Returns $ab \bmod n$ for any $a$ and $b$.
:::

::: function
ulong nmod_addmul(ulong s, ulong a, ulong b, nmod_t mod)
ulong nmod_fmma(ulong a, ulong b, ulong c, ulong d, nmod_t mod)

Return $s + ab$ and $ab + cd$ modulo $n$, for reduced operands.
`nmod_addmul` forms $s + ab$ as a two-word integer and reduces it once.
:::

::: macro
NMOD_ADDMUL(r, a, b, mod)

Sets `r` to $(r + ab) \bmod n$, for reduced `r`, `a` and `b`. This is
`r = nmod_addmul(r, a, b, mod)`.
:::

::: function
ulong nmod_inv(ulong a, nmod_t mod)

Returns $a^{-1} \bmod n$ for reduced $a$. Throws an exception if
$\gcd(a, n) \ne 1$.
:::

::: function
ulong nmod_div(ulong a, ulong b, nmod_t mod)

Returns $a b^{-1} \bmod n$ for reduced $a$ and $b$. Throws an exception if
$\gcd(b, n) \ne 1$.
:::

::: function
int nmod_divides(ulong * a, ulong b, ulong c, nmod_t mod)

If $ac \equiv b \pmod n$ has a solution, sets `*a` to one solution and
returns 1; otherwise returns 0. Unlike `nmod_div`, this works when $c$ is
not invertible. Requires reduced $b$ and $c$.

::: details
The congruence $ac \equiv b \pmod n$ is solvable if and only if
$g = \gcd(c, n)$ divides $b$, and then it has exactly $g$ solutions
modulo $n$. `nmod_divides` computes $u$ with $u c \equiv g \pmod n$ by
the extended Euclidean algorithm and returns $a = (b/g)\, u \bmod n$,
which is not necessarily the least solution. If $c = 0$, then `*a` is set
to 0, and the return value is 1 if $b = 0$ and 0 otherwise. When there is
no solution and $c \ne 0$, `*a` is not modified.
:::
:::

## Powering {#powering}

::: function
ulong nmod_pow_ui(ulong a, ulong exp, nmod_t mod)
ulong nmod_pow_fmpz(ulong a, const fmpz_t exp, nmod_t mod)

Return $a^{\mathit{exp}} \bmod n$ for reduced $a$, with $a^0 = 1$ (which is
0 when $n = 1$). For a negative exponent, `nmod_pow_fmpz` returns
$(a^{-1})^{-\mathit{exp}}$ and throws `FLINT_IMPINV` if $a$ is not
invertible modulo $n$.

::: details
`nmod_pow_ui` uses binary exponentiation with `n_mulmod_preinv` on the
normalised modulus. For odd $n$ and $\mathit{exp} \ge 2^{11}$ it uses
[Montgomery multiplication](algorithms/arithmetic/montgomery) instead.
`nmod_pow_fmpz` calls `n_powmod2_fmpz_preinv`, after inverting $a$ with
`n_gcdinv` if the exponent is negative, as `n_powmod2_preinv` does.
:::
:::

::: function
ulong nmod_ui_pow_ui(ulong a, ulong exp, nmod_t mod)

Returns $a^{\mathit{exp}} \bmod n$ for any $a$.
:::

::: function
ulong nmod_2_pow_ui(ulong exp, nmod_t mod)

Returns $2^{\mathit{exp}} \bmod n$.

::: details
Exponents below $2^{20}$, and all exponents for even $n$, use
left-to-right binary exponentiation, in which multiplication by the
base 2 is a modular doubling. Larger exponents with odd $n$ use the
same method in Montgomery arithmetic.
:::
:::

::: function
ulong _nmod_pow_ui_binexp(ulong a, ulong exp, nmod_t mod)
ulong _nmod_pow_ui_redc(ulong a, ulong exp, nmod_t mod)
ulong _nmod_2_pow_ui_binexp(ulong exp, nmod_t mod)

The algorithms behind `nmod_pow_ui` and `nmod_2_pow_ui`, with the same
results. The first two require $\mathit{exp} \ge 1$ (they do not terminate
for $\mathit{exp} = 0$), and `_nmod_pow_ui_redc` requires odd $n$.
`_nmod_2_pow_ui_binexp` accepts any $n$ and any exponent.
:::

## Square roots {#sqrt}

::: function
int nmod_sqrt(ulong * res, ulong a, nmod_t mod)

For prime $n$ and reduced $a$: if $a$ is a square modulo $n$, sets `*res`
to the square root in $[0, n/2]$ and returns 1; otherwise sets `*res` to 0
and returns 0. Primality is not checked. For composite $n$, a return value
of 1 still gives a genuine square root, but 0 does not prove that none
exists.
:::

::: function
int _nmod_sqrt_small(ulong * res, ulong a, ulong n)

For odd $n$ and $a < n$, finds the least square root of $a$ modulo $n$ by
exhaustive search in $O(n)$ steps. Returns 1 and sets `*res` to the root,
or returns 0 and sets `*res` to 0 if there is none. The answer is correct
for composite $n$ too.
:::

::: macro
NMOD_SQRT_SMALL_CUTOFF

The bound 128 below which `nmod_sqrt` uses `_nmod_sqrt_small` for odd $n$.
:::

::: summary algorithms/number-theory/sqrt-mod-p :::

## Discrete logarithms {#discrete-log}

These functions compute discrete logarithms modulo a prime $p$ to the base
of the least primitive root, by the Pohlig–Hellman method with
baby-step giant-step. They are efficient when $p - 1$ has only small prime
factors. The guide [](nmod/discrete_log) shows how to use them, including
logarithms to other bases.

::: type
nmod_discrete_log_pohlig_hellman_struct
nmod_discrete_log_pohlig_hellman_t

Precomputed data for discrete logarithms modulo one prime: the
factorisation of $p - 1$, the primitive root and lookup tables. The field
`mod` is an `nmod_t` for $p$.
:::

::: function
void nmod_discrete_log_pohlig_hellman_init(nmod_discrete_log_pohlig_hellman_t L)
void nmod_discrete_log_pohlig_hellman_clear(nmod_discrete_log_pohlig_hellman_t L)

Initialise `L` and free its memory. An initialised `L` must be set up with
`nmod_discrete_log_pohlig_hellman_precompute_prime` before it is used.
:::

::: function
double nmod_discrete_log_pohlig_hellman_precompute_prime(nmod_discrete_log_pohlig_hellman_t L, ulong p)

Sets up `L` for the prime $p$, which is not checked, freeing any earlier
data. Returns an estimate of the number of multiplications modulo $p$
needed by one call to `nmod_discrete_log_pohlig_hellman_run`.
:::

::: function
ulong nmod_discrete_log_pohlig_hellman_primitive_root(const nmod_discrete_log_pohlig_hellman_t L)

Returns the base $g$ of the logarithms, the least primitive root modulo $p$.
:::

::: function
ulong nmod_discrete_log_pohlig_hellman_run(const nmod_discrete_log_pohlig_hellman_t L, ulong y)

Returns the $x \in [0, p - 1)$ with $g^x \equiv y \pmod p$. Requires
$0 < y < p$.

::: details
`L` is not modified, so several threads can call
`nmod_discrete_log_pohlig_hellman_run` with the same `L` at the same time.
For $y = 0$, which has no logarithm, the function throws an exception when
$p$ is odd.
:::
:::

::: summary algorithms/number-theory/discrete-log :::

## Montgomery arithmetic {#redc}

In Montgomery arithmetic modulo an odd $n$, a residue $x$ is represented by
$x R \bmod n$, where $R = B$ for the standard functions and
$R = 2^{\bits/2}$ for the `half` functions. Multiplication then needs no
division. Values are converted into this *Montgomery form* once and back at
the end. The prefix of a function name tells which representation it uses:

| Prefix | Requires | $R$ | Residues |
|---|---|---|---|
| `nmod_redc_` | odd $n$ | $2^{\bits}$ | $[0, n)$ |
| `nmod_redc_fast_` | odd $n < 2^{\bits - 2}$ | $2^{\bits}$ | $[0, 2n)$ |
| `nmod_redc_half_` | odd $n < 2^{\bits/2 - 1}$ | $2^{\bits/2}$ | $[0, n)$ |
| `nmod_redc_half_fast_` | odd $n < 2^{\bits/2 - 2}$ | $2^{\bits/2}$ | $[0, 2n)$ |

The `fast` functions take and return *lazy* residues in $[0, 2n)$, which
saves a conditional subtraction per multiplication. The guide
[](nmod/redc) explains when Montgomery form pays off and the rules for
mixing standard, `fast` and `half` functions.

::: summary algorithms/arithmetic/montgomery :::

::: type
nmod_redc_ctx_struct
nmod_redc_ctx_t

A context for Montgomery arithmetic. The field `mod` is an ordinary
`nmod_t` for $n$, used for conversions and for addition, and `nred` is
$-n^{-1} \bmod R$. The context needs no clearing.
:::

::: function
void nmod_redc_ctx_init_nmod(nmod_redc_ctx_t ctx, nmod_t mod)
void nmod_redc_ctx_init_ui(nmod_redc_ctx_t ctx, ulong n)

Initialise `ctx` for the standard and `fast` functions modulo $n$, which
must be odd.
:::

::: function
int nmod_redc_can_use_fast(const nmod_redc_ctx_t ctx)

Returns whether $n < 2^{\bits - 2}$, the condition for the `fast`
functions.
:::

::: function
ulong nmod_redc_set_nmod(ulong x, const nmod_redc_ctx_t ctx)
ulong nmod_redc_set_ui(ulong x, const nmod_redc_ctx_t ctx)

Return the Montgomery form $x R \bmod n$ of $x$. `nmod_redc_set_nmod`
requires $x < n$; `nmod_redc_set_ui` accepts any $x$.
:::

::: function
ulong nmod_redc_get_nmod(ulong x, const nmod_redc_ctx_t ctx)

Returns $x R^{-1} \bmod n$ in $[0, n)$, converting from Montgomery form.
Any word $x$ is accepted, in particular a lazy residue.
:::

::: function
ulong nmod_redc_neg(ulong x, const nmod_redc_ctx_t ctx)
ulong nmod_redc_add(ulong x, ulong y, const nmod_redc_ctx_t ctx)
ulong nmod_redc_sub(ulong x, ulong y, const nmod_redc_ctx_t ctx)
ulong nmod_redc_mul(ulong x, ulong y, const nmod_redc_ctx_t ctx)

Negation, addition, subtraction and multiplication in Montgomery form, for
inputs and output in $[0, n)$. Negation, addition and subtraction are the
same as for ordinary residues.
:::

::: function
ulong nmod_redc_fast_neg(ulong x, const nmod_redc_ctx_t ctx)
ulong nmod_redc_fast_add(ulong x, ulong y, const nmod_redc_ctx_t ctx)
ulong nmod_redc_fast_sub(ulong x, ulong y, const nmod_redc_ctx_t ctx)
ulong nmod_redc_fast_mul(ulong x, ulong y, const nmod_redc_ctx_t ctx)

The same operations on lazy residues: inputs and output are in $[0, 2n)$.
Require $n < 2^{\bits - 2}$.
:::

::: function
ulong nmod_redc_fast_mul_two(ulong x, const nmod_redc_ctx_t ctx)

Returns a lazy residue in $[0, 2n)$ congruent to $2x$, for $x \in [0, 2n)$.
:::

::: function
ulong nmod_redc_fast_normalise(ulong x, const nmod_redc_ctx_t ctx)

Reduces a lazy residue $x \in [0, 2n)$ to $[0, n)$.
:::

::: function
ulong _nmod_redc_pow_ui(ulong a, ulong exp, const nmod_redc_ctx_t ctx)
ulong _nmod_redc_fast_pow_ui(ulong a, ulong exp, const nmod_redc_ctx_t ctx)

Return $a^{\mathit{exp}}$ with $a$ and the result in Montgomery form, by
binary exponentiation. Require $\mathit{exp} \ge 1$. The `fast` version
takes and returns lazy residues.
:::

::: function
ulong _nmod_redc_2_pow_ui(ulong exp, const nmod_redc_ctx_t ctx)
ulong _nmod_redc_fast_2_pow_ui(ulong exp, const nmod_redc_ctx_t ctx)

Return the Montgomery form of $2^{\mathit{exp}}$, for any $\mathit{exp}$.
The `fast` version returns a lazy residue.
:::

### Half-word Montgomery arithmetic {#redc-half}

These functions use $R = 2^{\bits/2}$, so that every product fits in one
word. They take the same context type, but it must be initialised with
the functions below, and the resulting Montgomery form differs from that
of the full-word functions, so the two cannot be mixed.

::: function
void nmod_redc_half_ctx_init_nmod(nmod_redc_ctx_t ctx, nmod_t mod)
void nmod_redc_half_ctx_init_ui(nmod_redc_ctx_t ctx, ulong n)

Initialise `ctx` for the `half` functions modulo $n$, which must be odd
with $n < 2^{\bits/2 - 1}$.
:::

::: function
int nmod_redc_half_can_use_fast(const nmod_redc_ctx_t ctx)

Returns whether $n < 2^{\bits/2 - 2}$, the condition for the `half_fast`
functions.
:::

::: function
ulong nmod_redc_half_set_nmod(ulong x, const nmod_redc_ctx_t ctx)
ulong nmod_redc_half_set_ui(ulong x, const nmod_redc_ctx_t ctx)
ulong nmod_redc_half_get_nmod(ulong x, const nmod_redc_ctx_t ctx)

Convert to and from Montgomery form with $R = 2^{\bits/2}$.
`nmod_redc_half_set_nmod` requires $x < n$ and `nmod_redc_half_set_ui`
accepts any $x$. `nmod_redc_half_get_nmod` returns a value in $[0, n)$ for
any $x < nR$, in particular for lazy residues.
:::

::: function
ulong nmod_redc_half_add(ulong x, ulong y, const nmod_redc_ctx_t ctx)
ulong nmod_redc_half_sub(ulong x, ulong y, const nmod_redc_ctx_t ctx)
ulong nmod_redc_half_mul(ulong x, ulong y, const nmod_redc_ctx_t ctx)

Addition, subtraction and multiplication in Montgomery form, for inputs
and output in $[0, n)$.
:::

::: function
ulong nmod_redc_half_fast_add(ulong x, ulong y, const nmod_redc_ctx_t ctx)
ulong nmod_redc_half_fast_sub(ulong x, ulong y, const nmod_redc_ctx_t ctx)
ulong nmod_redc_half_fast_mul(ulong x, ulong y, const nmod_redc_ctx_t ctx)

The same on lazy residues in $[0, 2n)$. Require $n < 2^{\bits/2 - 2}$.
:::

### Low-level functions {#redc-low}

The context functions are wrappers of the following ones, which take the
odd modulus $n$ and `nred` $= -n^{-1} \bmod R$ directly. Here $R = B$,
except for the `half` functions, where $R = 2^{\bits/2}$.

::: function
ulong n_redc(ulong x, ulong n, ulong nred)
ulong n_redc_fast(ulong x, ulong n, ulong nred)

Return $x R^{-1} \bmod n$ for any word $x$. `n_redc_fast` returns a value
in $[0, n]$, which may need one more subtraction of $n$.
:::

::: function
ulong n_ll_redc(ull_t x, ulong n, ulong nred)
ulong n_ll_redc_fast(ull_t x, ulong n, ulong nred)

Return a value congruent to $x R^{-1}$ modulo $n$, for a two-word $x$.
`n_ll_redc` requires $x < nB$ and returns a value in $[0, n)$.
`n_ll_redc_fast` requires $x + (B - 1) n < B^2$ and returns a value below
$x/B + n$.
:::

::: function
ulong n_mulmod_redc(ulong a, ulong b, ulong n, ulong nred)
ulong n_mulmod_redc_fast(ulong a, ulong b, ulong n, ulong nred)

Return a value congruent to $ab R^{-1}$ modulo $n$. `n_mulmod_redc`
requires $ab < nB$ (for example $a, b < n$) and returns a value in
$[0, n)$. `n_mulmod_redc_fast` maps $a, b \in [0, 2n)$ to $[0, 2n)$ when
$n < 2^{\bits - 2}$.
:::

::: function
ulong n_to_redc_preinv(ulong a, ulong n, ulong ninv, ulong norm)

Returns $a B \bmod n$ for $a < n$, given `ninv` and `norm` as stored in an
`nmod_t` for $n$. This is one division step with the precomputed inverse.
:::

::: function
ulong n_redc_half(ulong x, ulong n, ulong nred)
ulong n_redc_half_fast(ulong x, ulong n, ulong nred)
ulong n_mulmod_redc_half(ulong a, ulong b, ulong n, ulong nred)
ulong n_mulmod_redc_half_fast(ulong a, ulong b, ulong n, ulong nred)
ulong n_to_redc_half_preinv(ulong a, ulong n, ulong ninv, ulong norm)

The same with $R = 2^{\bits/2}$. The reductions require
$x + (R - 1) n < B$, and the ranges of inputs and outputs are those of the
corresponding `nmod_redc_half` functions. `n_to_redc_half_preinv` ignores
`norm`.
:::

::: type
ull_t

An unsigned two-word integer: `__uint128_t` when compiling for a 64-bit
machine with a GCC-compatible compiler, and otherwise a struct of two
words `lo` and `hi`. It exists for the Montgomery helpers above and may
change.
:::

::: function
ull_t ull(ulong hi, ulong lo)
ulong ull_hi(ull_t x)
ulong ull_lo(ull_t x)

Construct the two-word integer $\mathit{hi} B + \mathit{lo}$, and return
the high or the low word of $x$.
:::

::: function
ull_t ull_add(ull_t x, ull_t y)
ull_t ull_add_u(ull_t x, ulong y)
ull_t ull_u_mul_u(ulong x, ulong y)

Return $(x + y) \bmod B^2$, and the full product $xy$ of two words.
:::

## Generic rings {#gr}

These functions create `gr` contexts for $\Z/n\Z$. Those without `redc` in
the name are declared in `gr.h`.

::: function
int gr_ctx_init_nmod(gr_ctx_t ctx, ulong n)
int gr_ctx_init_nmod8(gr_ctx_t ctx, ulong n)
int gr_ctx_init_nmod32(gr_ctx_t ctx, ulong n)
int gr_ctx_init_nmod_redc(gr_ctx_t ctx, ulong n)
int gr_ctx_init_nmod_redc_fast(gr_ctx_t ctx, ulong n)

Initialise `ctx` for the ring $\Z/n\Z$ in the given representation and
return `GR_SUCCESS`. If $n = 0$, return `GR_DOMAIN`; if $n$ violates the
restrictions of the representation, return `GR_UNABLE`. In both cases
`ctx` is not initialised.

::: details
| Function | Element type | Representation | Allowed $n$ |
|---|---|---|---|
| `gr_ctx_init_nmod` | `ulong` | reduced residue | $n \ge 1$ |
| `gr_ctx_init_nmod8` | `uint8_t` | reduced residue | $1 \le n \le 255$ |
| `gr_ctx_init_nmod32` | `uint32_t` | reduced residue | $1 \le n < 2^{32}$ |
| `gr_ctx_init_nmod_redc` | `ulong` | Montgomery form in $[0, n)$ | odd $n \ge 3$ |
| `gr_ctx_init_nmod_redc_fast` | `ulong` | Montgomery form in $[0, 2n)$ | odd $3 \le n < 2^{\bits - 2}$ |

Some generic algorithms ask whether the ring is a field. For the `nmod`,
`nmod_redc` and `nmod_redc_fast` contexts, the answer is computed with
`n_is_prime` on every query, unless it has been set with
`gr_ctx_set_is_field(ctx, T_TRUE)` or `gr_ctx_set_is_field(ctx, T_FALSE)`.
Setting it is recommended when the primality of $n$ is known. It is not
done automatically, as it would slow down creating a context that is only
used for basic arithmetic. The `nmod8` and `nmod32` contexts behave in the
same way.

Many operations of the `nmod8`, `nmod32`, `nmod_redc` and `nmod_redc_fast`
contexts are not as optimised as those of `nmod`. Use `nmod8` and `nmod32`
mainly to save memory.
:::
:::

::: function
void _gr_ctx_init_nmod(gr_ctx_t ctx, void * nmod_t_ref)

Initialises `ctx` like `gr_ctx_init_nmod`, copying the `nmod_t` that
`nmod_t_ref` points to instead of recomputing it.
:::

::: macro
NMOD_CTX(ctx)
NMOD_CTX_REF(ctx)

The `nmod_t` stored in a context created by `gr_ctx_init_nmod`, and a
pointer to it.
:::

::: macro
GR_NMOD_REDC_CTX(ctx)

A pointer to the `nmod_redc_ctx_struct` stored in a context created by
`gr_ctx_init_nmod_redc` or `gr_ctx_init_nmod_redc_fast`.
:::

## Low-level macros {#reduction}

These macros are the building blocks of the arithmetic functions above
and of the vector, polynomial and matrix modules. Most reduce one, two or
three words modulo $n$ with one or two steps of
[Möller–Granlund division](algorithms/arithmetic/moller-granlund).
Unless stated otherwise, the first argument
receives the result and `mod` is an `nmod_t`. A pair `a_hi, a_lo` stands
for the two-word integer $a_{hi} B + a_{lo}$.

::: macro
NMOD_RED(r, a, mod)

Sets `r` to $a \bmod n$ for any word $a$.
:::

::: macro
NMOD_RED2(r, a_hi, a_lo, mod)

Sets `r` to $(a_{hi} B + a_{lo}) \bmod n$. Requires $a_{hi} < n$.
:::

::: macro
NMOD2_RED2(r, a_hi, a_lo, mod)

Sets `r` to $(a_{hi} B + a_{lo}) \bmod n$ with no restriction on
$a_{hi}$. It reduces $a_{hi}$ first, so it costs two division steps.
:::

::: macro
NMOD_RED3(r, a_hi, a_me, a_lo, mod)

Sets `r` to $(a_{hi} B^2 + a_{me} B + a_{lo}) \bmod n$. Requires
$a_{hi} < n$.
:::

::: macro
NMOD_RED2_FULLWORD(r, a_hi, a_lo, mod)
NMOD_RED2_NONFULLWORD(r, a_hi, a_lo, mod)

Versions of `NMOD_RED2` without the test on `mod.norm`.
`NMOD_RED2_FULLWORD` requires `mod.norm` $= 0$ (the top bit of $n$ is set)
and `NMOD_RED2_NONFULLWORD` requires `mod.norm` $> 0$. Both require
$a_{hi} < n$.
:::

::: macro
NMOD_MUL_PRENORM(res, a, b, mod)
NMOD_MUL_FULLWORD(res, a, b, mod)

Multiply without the operand shift done by `nmod_mul`.
`NMOD_MUL_PRENORM` requires $a < n$ and $b = c\, 2^{\mathit{norm}}$ with
$c < n$, that is, $b$ is a reduced value already shifted left by
`mod.norm` bits, and sets `res` to the reduced residue $ac \bmod n$;
`nmod_mul(a, c, mod)` is this macro applied to `a` and `c << mod.norm`.
`NMOD_MUL_FULLWORD` requires `mod.norm` $= 0$ and reduced operands, and
sets `res` to $ab \bmod n$.
:::

::: macro
NMOD_CAN_USE_SHOUP(mod)

Nonzero if $n < 2^{\bits - 1}$, the condition for multiplication by a
fixed operand with `n_mulmod_shoup`.
:::
