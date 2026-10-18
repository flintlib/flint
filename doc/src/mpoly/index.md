---
title: mpoly.h
header: mpoly.h
description: Infrastructure shared by all sparse multivariate polynomial modules: monomial orderings, packed exponent vectors, monomial arithmetic, heaps, and helpers for GCD and factorisation.
category: Multivariate polynomials
order: 1
---

This module contains no polynomial type of its own. It implements what the
sparse multivariate modules (`fmpz_mpoly`, `fmpq_mpoly`, `nmod_mpoly`,
`fmpz_mod_mpoly`, `fq_nmod_mpoly`, `fq_zech_mpoly` and `gr_mpoly`) have in
common: the context `mpoly_ctx_t` that fixes the number of variables and
the monomial ordering, the packing of exponent vectors into machine words,
arithmetic and comparison on packed monomials, the heaps used for sparse
multiplication and division, and helpers for the GCD and factoring code.
Most users never call these functions directly; they are needed to write
new algorithms on the internal representation, and they explain how that
representation works.

The representation is described in detail, with the proofs that make it
work, on the algorithm page [](algorithms/multivariate/representation).
The short version: a polynomial is an array of coefficients together with
an array of monomials, sorted from the largest monomial to the smallest.
A monomial is stored as $N$ words holding the exponents (and, for degree
orderings, the total degree) as bit fields of a common width *bits*,
arranged so that comparing monomials amounts to comparing words, possibly
after an exclusive or with a mask, and multiplying monomials amounts to
adding words. Each field has a spare high bit so that overflow of an
addition can be detected.

::: summary algorithms/multivariate/representation :::

## Conventions {#conventions}

- *User exponent vectors* are arrays of length `nvars` of type `ulong *`,
  `slong *`, `fmpz *` (written *ffmpz* in function names: a flat array of
  `fmpz`) or `fmpz **` (written *pfmpz*: an array of pointers). Entry $i$
  is the exponent of the variable of index $i$, and the variable of index
  0 is the most significant in the ordering. With variables $x > y > z$,
  the polynomial $7x^2y + 8yz + 9$ has degrees $(2, 1, 1)$, its term of
  index 0 has exponent vector $(2, 1, 0)$ and its term of index 2 has
  exponent vector $(0, 0, 0)$.
- *Packed monomials* (often called `exps`, `poly_exps` or `Aexps`) are
  arrays of $N$ words per monomial, where `N = mpoly_words_per_exp(bits,
  mctx)`. Functions on single packed monomials take $N$ (and sometimes
  `bits` or a mask) instead of the context. Functions on arrays of
  monomials take the array, its length and `bits`.
- *Unpacked fields* are the `nfields` integers stored in a packed monomial,
  in storage order (least significant field first). They differ from user
  exponents in order, and for degree orderings by the extra total degree
  field.
- A packing width `bits` is valid if $8 \le \mathit{bits} \le$
  `FLINT_BITS` (*single precision*, functions with suffix `_sp`) or if it
  is a multiple of `FLINT_BITS` (*multiprecision*, suffix `_mp`). Every
  stored field must be less than $2^{\mathit{bits} - 1}$: the top bit of
  each field is a guard bit that is zero in every valid monomial.
- *Masks.* `cmpmask` is the array of $N$ words returned by
  `mpoly_get_cmpmask`; monomials compare like the words `exp[i] ^
  cmpmask[i]` from $i = N - 1$ down to $0$. The overflow mask returned by
  `mpoly_overflow_mask_sp` has the guard bit of each field set.

## Orderings and contexts {#context}

::: type
ordering_t

An enumeration of the supported monomial orderings: `ORD_LEX`,
`ORD_DEGLEX` and `ORD_DEGREVLEX`.

::: details
With variables $x_0 > x_1 > \cdots > x_{n-1}$ and exponent vectors $a$,
$b$:

- `ORD_LEX`: $x^a > x^b$ if the first nonzero entry of $a - b$ is positive.
- `ORD_DEGLEX`: $x^a > x^b$ if $|a| > |b|$, or $|a| = |b|$ and $x^a > x^b$
  in lex order, where $|a| = a_0 + \cdots + a_{n-1}$.
- `ORD_DEGREVLEX`: $x^a > x^b$ if $|a| > |b|$, or $|a| = |b|$ and the last
  nonzero entry of $a - b$ is negative.

All three are monomial orderings: total orders compatible with
multiplication in which $1$ is the smallest monomial.
:::
:::

::: macro
ORD_LEX
ORD_DEGLEX
ORD_DEGREVLEX

The three values of `ordering_t`: lexicographic, degree lexicographic and
degree reverse lexicographic ordering.
:::

::: macro
MPOLY_NUM_ORDERINGS

The number of supported orderings, 3.
:::

::: type
mpoly_ctx_struct
mpoly_ctx_t

The monomial part of a polynomial ring context. The fields are `nvars`
(the number of variables), `nfields` (the number of fields of a packed
monomial: `nvars`, plus one for degree orderings), `ord`, the flags `deg`
(the ordering compares total degrees first) and `rev` (the variables are
stored in reverse order and compared with reversed sign), and two lookup
tables indexed by `bits - 1` for $1 \le \mathit{bits} \le$ `FLINT_BITS`:
`lut_words_per_exp` and `lut_fix_bits`. The contexts of the polynomial
modules, such as `fmpz_mpoly_ctx_t`, contain an `mpoly_ctx_t` called
`minfo`.
:::

::: function
void mpoly_ctx_init(mpoly_ctx_t ctx, slong nvars, const ordering_t ord)

Initialises `ctx` for `nvars` variables and the ordering `ord`. With
`nvars` equal to zero the context describes constants; it then behaves
like a degree ordering with one field that is always zero. An invalid
ordering throws `FLINT_ERROR`.
:::

::: function
void mpoly_ctx_init_rand(mpoly_ctx_t mctx, flint_rand_t state, slong max_nvars)

Initialises `mctx` with a random ordering and a random number of
variables between 0 and `max_nvars`, for testing.
:::

::: function
void mpoly_ctx_clear(mpoly_ctx_t mctx)

Releases the context. The context owns no memory, so this does nothing,
but it should be called for symmetry with `mpoly_ctx_init`.
:::

::: function
ordering_t mpoly_ordering_randtest(flint_rand_t state)

Returns a random ordering.
:::

::: function
int mpoly_ordering_isdeg(const mpoly_ctx_t mctx)
int mpoly_ordering_isrev(const mpoly_ctx_t mctx)

Return whether the ordering of `mctx` is a degree ordering (`ORD_DEGLEX`
or `ORD_DEGREVLEX`), respectively a reversed ordering (`ORD_DEGREVLEX`).
:::

::: function
void mpoly_ordering_print(ordering_t ord)

Prints `lex`, `deglex` or `degrevlex` to standard output.
:::

## Packing parameters and masks {#packing}

::: macro
MPOLY_MIN_BITS

The smallest packing width used by the polynomial modules, 8 bits. Freshly
initialised polynomials use this width, and it grows as needed.
:::

::: function
slong mpoly_words_per_exp_sp(flint_bitcnt_t bits, const mpoly_ctx_t mctx)
slong mpoly_words_per_exp_mp(flint_bitcnt_t bits, const mpoly_ctx_t mctx)
slong mpoly_words_per_exp(flint_bitcnt_t bits, const mpoly_ctx_t mctx)

Return the number $N$ of words of a monomial packed with fields of `bits`
bits. In the single precision case $1 \le \mathit{bits} \le$ `FLINT_BITS`
this is $\lceil \mathit{nfields} / \lfloor \text{FLINT\_BITS}/\mathit{bits}
\rfloor \rceil$ (read from a table); in the multiprecision case `bits`
must be a multiple of `FLINT_BITS` and $N = \mathit{nfields} \cdot
\mathit{bits} / \text{FLINT\_BITS}$. The unsuffixed version dispatches on
`bits`.
:::

::: function
flint_bitcnt_t mpoly_fix_bits(flint_bitcnt_t bits, const mpoly_ctx_t mctx)

Rounds a required number of bits up to a valid packing width. For $\mathit{bits}
\le$ `FLINT_BITS` the result is the largest width that is at least
`MPOLY_MIN_BITS` and at least `bits` and needs no more words than `bits`
itself, so no space is wasted; for larger `bits` it is the next multiple
of `FLINT_BITS`. For example, with three variables and lex ordering, any
`bits` up to 21 is rounded to 21 on a 64-bit machine.
:::

::: function
void mpoly_get_cmpmask(ulong * cmpmask, slong N, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

Sets `(cmpmask, N)` to the comparison mask for monomials packed with
`bits` bits: zero for `ORD_LEX` and `ORD_DEGLEX`, and for `ORD_DEGREVLEX`
all ones in every field except the total degree field. Two packed
monomials compare in the monomial ordering as the multiword integers
obtained by an exclusive or with this mask. Every function that compares
monomials takes such a mask.
:::

::: function
ulong mpoly_overflow_mask_sp(flint_bitcnt_t bits)

Returns the word with the high (guard) bit of each `bits`-bit field set,
for $\mathit{bits} \le$ `FLINT_BITS`. A packed word is valid if and only
if its bitwise and with this mask is zero.
:::

::: macro
MPOLY_GET_CMPMASK_TMP_ALLOC(cmpmask, N, bits, mctx)
MPOLY_GET_CMPMASK_FLINT_MALLOC(cmpmask, N, bits, mctx)

Set `N` to `mpoly_words_per_exp(bits, mctx)`, allocate `N` words for
`cmpmask` and fill them with `mpoly_get_cmpmask`. The first macro
allocates with `TMP_ALLOC` and may only be used between `TMP_START` and
`TMP_END`; the second allocates with `flint_malloc`, and the caller must
free the mask with `flint_free` (use it when the mask must outlive the
current `TMP` scope, for instance when it is passed to worker threads).
:::

## Converting monomials {#convert}

These functions translate between user exponent vectors, unpacked fields
and packed monomials.

::: function
flint_bitcnt_t mpoly_exp_bits_required_ui(const ulong * user_exp, const mpoly_ctx_t mctx)
flint_bitcnt_t mpoly_exp_bits_required_ffmpz(const fmpz * user_exp, const mpoly_ctx_t mctx)
flint_bitcnt_t mpoly_exp_bits_required_pfmpz(fmpz * const * user_exp, const mpoly_ctx_t mctx)

Return the number of bits needed to pack the user exponent vector
`user_exp`, including the guard bit: one more than the bit length of the
largest exponent, or of the total degree for degree orderings. The result
is not rounded with `mpoly_fix_bits`. For `_ui`, a total degree that
overflows a word gives `2*FLINT_BITS`.
:::

::: function
flint_bitcnt_t mpoly_gen_pow_exp_bits_required(slong v, ulong e, const mpoly_ctx_t mctx)

Returns the number of bits needed to pack the monomial $x_v^e$, that is
$1 + \lceil \log_2(e + 1) \rceil$.
:::

::: function
void mpoly_set_monomial_ui(ulong * exp1, const ulong * exp2, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_set_monomial_ffmpz(ulong * exp1, const fmpz * exp2, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_set_monomial_pfmpz(ulong * exp1, fmpz * const * exp2, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

Pack the user exponent vector `exp2` into the monomial `exp1` with
`bits` bits per field, computing the total degree field for degree
orderings. It is not checked that the exponents fit; use
`mpoly_exp_bits_required_ui` and friends first.
:::

::: function
void mpoly_get_monomial_ui(ulong * user_exps, const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_get_monomial_ui_sp(ulong * user_exps, const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_get_monomial_ui_mp(ulong * user_exps, const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_get_monomial_si(slong * user_exps, const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_get_monomial_si_mp(slong * user_exps, const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

Unpack the packed monomial `poly_exps` into the user exponent vector
`user_exps`. The `_sp` versions require $\mathit{bits} \le$ `FLINT_BITS`
and the `_mp` versions a multiple of `FLINT_BITS` greater than
`FLINT_BITS`; the unsuffixed versions dispatch. The `_mp` versions throw
`FLINT_ERROR` if an exponent does not fit a `ulong` (respectively an
`slong`).
:::

::: function
void mpoly_get_monomial_ffmpz(fmpz * exps, const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_get_monomial_pfmpz(fmpz ** exps, const ulong * poly_exps, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

Unpack the packed monomial `poly_exps` into a user exponent vector of
`fmpz`, for any packing width.
:::

::: function
ulong mpoly_get_monomial_var_exp_ui(const ulong * poly_exps, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
ulong mpoly_get_monomial_var_exp_ui_sp(const ulong * poly_exps, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
ulong mpoly_get_monomial_var_exp_ui_mp(const ulong * poly_exps, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
slong mpoly_get_monomial_var_exp_si(const ulong * poly_exps, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
slong mpoly_get_monomial_var_exp_si_mp(const ulong * poly_exps, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

Return the exponent of the variable of index `var` in the packed monomial
`poly_exps`, without unpacking the other fields. The multiprecision
versions throw `FLINT_ERROR` if the exponent does not fit the return
type.
:::

::: function
void mpoly_get_monomial_ui_unpacked_ui(ulong * user_exps, const ulong * poly_exps, const mpoly_ctx_t mctx)
void mpoly_get_monomial_ui_unpacked_ffmpz(ulong * user_exps, const fmpz * poly_exps, const mpoly_ctx_t mctx)
void mpoly_get_monomial_ffmpz_unpacked_ffmpz(fmpz * user_exps, const fmpz * poly_exps, const mpoly_ctx_t mctx)
void mpoly_get_monomial_pfmpz_unpacked_ffmpz(fmpz ** user_exps, const fmpz * poly_exps, const mpoly_ctx_t mctx)

Convert an already unpacked array of `nfields` fields into a user
exponent vector, by reordering the fields and dropping the total degree
field. They are used together with `mpoly_max_fields_fmpz` to compute
degrees. The `_ui` versions assume that the fields fit a `ulong`.
:::

::: function
void mpoly_pack_vec_ui(ulong * exp1, const ulong * exp2, flint_bitcnt_t bits, slong nfields, slong len)
void mpoly_pack_vec_fmpz(ulong * exp1, const fmpz * exp2, flint_bitcnt_t bits, slong nfields, slong len)

Pack `len` consecutive vectors of `nfields` fields from `exp2` into
`exp1`, with `bits` bits per field. In the single precision case the
fields of one vector are written from the least significant bits of the
first word up, starting a new word whenever the next field does not fit,
and the unused high bits are zero. It is not checked that the fields fit.
:::

::: function
void mpoly_unpack_vec_ui(ulong * exp1, const ulong * exp2, flint_bitcnt_t bits, slong nfields, slong len)
void mpoly_unpack_vec_fmpz(fmpz * exp1, const ulong * exp2, flint_bitcnt_t bits, slong nfields, slong len)

Unpack `len` packed vectors of `nfields` fields of `bits` bits from
`exp2` into `exp1`. The `_ui` version requires the fields to fit a
`ulong`.
:::

::: function
int mpoly_repack_monomials(ulong * exps1, flint_bitcnt_t bits1, const ulong * exps2, flint_bitcnt_t bits2, slong len, const mpoly_ctx_t mctx)

Converts `len` monomials packed with `bits2` bits per field in `exps2`
into monomials packed with `bits1` bits in `exps1`, which must have room
for `len*mpoly_words_per_exp(bits1, mctx)` words. Returns 1 for success.
Widening always succeeds; narrowing returns 0 as soon as a field does not
fit in `bits1 - 1` bits, leaving `exps1` partially written.
:::

::: function
ulong * mpoly_monomials_repack_wider(flint_bitcnt_t * bits_ptr, slong * N_ptr, const ulong * old_exps, flint_bitcnt_t old_bits, slong len, slong alloc, const mpoly_ctx_t mctx)
ulong * mpoly_monomials_repack_wider_cmpmask(flint_bitcnt_t * bits_ptr, slong * N_ptr, ulong ** cmpmask_ptr, const ulong * old_exps, flint_bitcnt_t old_bits, slong len, slong alloc, const mpoly_ctx_t mctx)

Helpers for algorithms that detect exponent overflow while producing
their output. They replace `*bits_ptr` by `mpoly_fix_bits(*bits_ptr + 1,
mctx)`, update `*N_ptr`, and return a new array, allocated with
`flint_malloc` with room for `alloc` monomials, containing the first
`len` monomials of `old_exps` (packed with `old_bits` bits) repacked to
the new width. `old_exps` is not freed. The second version also
reallocates `*cmpmask_ptr` with `flint_realloc` and refills it, so it must
have been allocated with `flint_malloc`.
:::

::: function
int mpoly_term_exp_fits_ui(ulong * exps, flint_bitcnt_t bits, slong n, const mpoly_ctx_t mctx)
int mpoly_term_exp_fits_si(ulong * exps, flint_bitcnt_t bits, slong n, const mpoly_ctx_t mctx)

Return whether every exponent of the monomial of index $n$ in the array
`exps` fits a `ulong`, respectively an `slong`.
:::

::: function
void mpoly_monomial_randbits_fmpz(fmpz * exp, flint_rand_t state, flint_bitcnt_t exp_bits, const mpoly_ctx_t mctx)

Sets `exp` to a random user exponent vector that can be packed into
`exp_bits` bits (including the guard bit), for testing.
:::

## Arithmetic on packed monomials {#arithmetic}

All functions in this section are inline. They operate word by word on
monomials of $N$ words and rely on the guard bits: as long as the result
of each field fits below its guard bit, word arithmetic is field
arithmetic. Except where stated, the output may alias the inputs.

::: function
void mpoly_monomial_zero(ulong * exp_ptr, slong N)
void mpoly_monomial_set(ulong * exp2, const ulong * exp3, slong N)
void mpoly_monomial_swap(ulong * exp2, ulong * exp3, slong N)
void mpoly_monomial_set_extra(ulong * exp2, const ulong * exp3, slong N, slong offset, ulong extra)

Set a monomial to $1$ (all fields zero), copy a monomial, swap two
monomials. `mpoly_monomial_set_extra` copies `exp3` and adds `extra` to
the word of index `offset`.
:::

::: function
void mpoly_copy_monomials(ulong * exp1, const ulong * exp2, slong len, slong N)

Copies `len` monomials of $N$ words. With GCC this is a macro for
`memcpy`, so the arrays must not overlap.
:::

::: function
void mpoly_monomial_add(ulong * exp_ptr, const ulong * exp2, const ulong * exp3, slong N)
void mpoly_monomial_add_mp(ulong * exp_ptr, const ulong * exp2, const ulong * exp3, slong N)
void mpoly_monomial_sub(ulong * exp_ptr, const ulong * exp2, const ulong * exp3, slong N)
void mpoly_monomial_sub_mp(ulong * exp_ptr, const ulong * exp2, const ulong * exp3, slong N)

Set `exp_ptr` to the product (sum of fields) or quotient (difference of
fields) of the monomials `exp2` and `exp3`. The plain versions work
word by word and are for single precision fields; the `_mp` versions
propagate carries across words and are for multiprecision fields. The
caller is responsible for overflow (add) and for divisibility (sub).
:::

::: function
void mpoly_monomial_madd(ulong * exp1, const ulong * exp2, ulong scalar, const ulong * exp3, slong N)
void mpoly_monomial_madd_mp(ulong * exp1, const ulong * exp2, ulong scalar, const ulong * exp3, slong N)
void mpoly_monomial_madd_inplace_mp(ulong * exp12, ulong scalar, const ulong * exp3, slong N)
void mpoly_monomial_msub(ulong * exp1, const ulong * exp2, ulong scalar, const ulong * exp3, slong N)
void mpoly_monomial_msub_mp(ulong * exp1, const ulong * exp2, ulong scalar, const ulong * exp3, slong N)

Set `exp1` to `exp2` plus (minus) `scalar` times `exp3`, that is, to
$\mathit{exp2} \cdot \mathit{exp3}^{\pm \mathit{scalar}}$ as monomials;
`madd_inplace_mp` adds to `exp12`. The `_mp` versions use `mpn`
arithmetic over all $N$ words.
:::

::: function
void mpoly_monomial_madd_ui_array(ulong * exp1, const ulong * exp2, const ulong * scalar, slong scalar_limbs, const ulong * exp3, slong N)
void mpoly_monomial_msub_ui_array(ulong * exp1, const ulong * exp2, const ulong * scalar, slong scalar_limbs, const ulong * exp3, slong N)
void mpoly_monomial_madd_fmpz(ulong * exp1, const ulong * exp2, const fmpz_t scalar, const ulong * exp3, slong N)

Multiprecision versions of `madd` and `msub` with a scalar of
`scalar_limbs` words (at most $N$) or a nonnegative `fmpz`, computed modulo
$2^{N \cdot \text{FLINT\_BITS}}$. `exp1` must not alias `exp3`.
:::

::: function
void mpoly_monomial_mul_ui(ulong * exp2, const ulong * exp3, slong N, ulong c)
void mpoly_monomial_mul_ui_mp(ulong * exp2, const ulong * exp3, slong N, ulong c)
void mpoly_monomial_mul_fmpz(ulong * exp2, const ulong * exp3, slong N, const fmpz_t c)

Set `exp2` to the $c$-th power of the monomial `exp3`, that is, multiply
every field by $c$. The plain version multiplies each word and is only
correct if no field overflows; the `_mp` and `fmpz` versions multiply the
$N$-word integer, which is correct for multiprecision fields. In
`mpoly_monomial_mul_fmpz`, $c \ge 0$ must have at most $N$ words.
:::

::: function
ulong mpoly_monomial_max1(ulong exp2, ulong exp3, flint_bitcnt_t bits, ulong mask)
ulong mpoly_monomial_min1(ulong exp2, ulong exp3, flint_bitcnt_t bits, ulong mask)
void mpoly_monomial_max(ulong * exp1, const ulong * exp2, const ulong * exp3, flint_bitcnt_t bits, slong N, ulong mask)
void mpoly_monomial_min(ulong * exp1, const ulong * exp2, const ulong * exp3, flint_bitcnt_t bits, slong N, ulong mask)
void mpoly_monomial_max_mp(ulong * exp1, const ulong * exp2, const ulong * exp3, flint_bitcnt_t bits, slong N)
void mpoly_monomial_min_mp(ulong * exp1, const ulong * exp2, const ulong * exp3, flint_bitcnt_t bits, slong N)

Set `exp1` to the field-wise maximum (least common multiple) or minimum
(greatest common divisor) of `exp2` and `exp3`. The single precision
versions take the overflow mask and compute all fields of a word at once
without branches: the guard bit of each field of `mask + exp2 - exp3`
records whether that field of `exp2` is at least that of `exp3`.
:::

::: function
int mpoly_monomial_overflows(ulong * exp2, slong N, ulong mask)
int mpoly_monomial_overflows1(ulong exp, ulong mask)
int mpoly_monomial_overflows_mp(ulong * exp_ptr, slong N, flint_bitcnt_t bits)

Return whether a field of the monomial has its guard bit set, that is,
whether an addition overflowed or a subtraction went negative. The single
precision versions take the overflow mask; the `_mp` version tests the
top bit of each multiword field.
:::

::: function
int mpoly_monomial_divides(ulong * exp_ptr, const ulong * exp2, const ulong * exp3, slong N, ulong mask)
int mpoly_monomial_divides1(ulong * exp_ptr, const ulong exp2, const ulong exp3, ulong mask)
int mpoly_monomial_divides_mp(ulong * exp_ptr, const ulong * exp2, const ulong * exp3, slong N, flint_bitcnt_t bits)

Return whether the monomial `exp3` divides `exp2`, and if so set
`exp_ptr` to the quotient. If not, `exp_ptr` holds garbage.
:::

::: function
int mpoly_monomial_divides_test(const ulong * exp2, const ulong * exp3, slong N, ulong mask)
int mpoly_monomial_divides_mp_test(const ulong * exp2, const ulong * exp3, slong N, flint_bitcnt_t bits)

Return whether `exp3` divides `exp2` without computing the quotient.
:::

::: function
int mpoly_monomial_halves(ulong * exp_ptr, const ulong * exp2, slong N, ulong mask)
int mpoly_monomial_halves1(ulong * exp_ptr, const ulong exp2, ulong mask)
int mpoly_monomial_halves_mp(ulong * exp_ptr, const ulong * exp2, slong N, flint_bitcnt_t bits)

Return whether the monomial `exp2` is a square, that is, every field is
even, and if so set `exp_ptr` to its square root. The single precision
versions shift whole words right by one bit, so that the low bit of each
field lands in the guard bit of the field below, where the overflow test
detects it.
:::

::: function
int mpoly_monomial_valid_any_bits(const ulong * exp, slong N, flint_bitcnt_t bits)
void mpoly_monomial_add_any_bits(ulong * exp1, const ulong * exp2, const ulong * exp3, slong N, flint_bitcnt_t bits)
void mpoly_monomial_sub_any_bits(ulong * exp1, const ulong * exp2, const ulong * exp3, slong N, flint_bitcnt_t bits)
void mpoly_monomial_madd_any_bits(ulong * exp1, const ulong * exp2, ulong scalar, const ulong * exp3, slong N, flint_bitcnt_t bits)
void mpoly_monomial_msub_any_bits(ulong * exp1, const ulong * exp2, ulong scalar, const ulong * exp3, slong N, flint_bitcnt_t bits)
void mpoly_monomial_max_any_bits(ulong * exp1, const ulong * exp2, const ulong * exp3, flint_bitcnt_t bits, slong N, ulong mask)
void mpoly_monomial_min_any_bits(ulong * exp1, const ulong * exp2, const ulong * exp3, flint_bitcnt_t bits, slong N, ulong mask)
int mpoly_monomial_overflows_any_bits(ulong * exp, slong N, ulong mask, flint_bitcnt_t bits)
int mpoly_monomial_divides_any_bits(ulong * exp_ptr, const ulong * exp2, const ulong * exp3, slong N, ulong mask, flint_bitcnt_t bits)
int mpoly_monomial_halves_any_bits(ulong * exp_ptr, const ulong * exp2, slong N, ulong mask, flint_bitcnt_t bits)
void mpoly_monomial_mul_ui_any_bits(ulong * exp2, const ulong * exp3, slong N, ulong c, flint_bitcnt_t bits)

Versions of the functions above that dispatch to the single precision or
the `_mp` version according to `bits`, for code that handles both packing
widths. The `mask` argument is only used in the single precision case.
`mpoly_monomial_valid_any_bits` returns whether no guard bit is set; in
builds with assertions the other functions assert that their inputs are
valid.
:::

## Comparing monomials {#compare}

::: function
int mpoly_monomial_is_zero(const ulong * exp, slong N)
int mpoly_monomial_equal(const ulong * exp2, const ulong * exp3, slong N)
int mpoly_monomial_equal_extra(const ulong * exp2, const ulong * exp3, slong N, slong offset, ulong extra)

Return whether a monomial is $1$, whether two monomials are equal, and
whether `exp2` equals `exp3` with `extra` added to its word of index
`offset`.
:::

::: function
int mpoly_monomial_cmp(const ulong * exp2, const ulong * exp3, slong N, const ulong * cmpmask)
int mpoly_monomial_lt(const ulong * exp3, const ulong * exp2, slong N, const ulong * cmpmask)
int mpoly_monomial_gt(const ulong * exp3, const ulong * exp2, slong N, const ulong * cmpmask)

Compare two packed monomials in the monomial ordering given by the mask
from `mpoly_get_cmpmask`. `mpoly_monomial_cmp` returns $1$, $0$ or $-1$
as `exp2` is greater than, equal to or less than `exp3`;
`mpoly_monomial_lt` and `mpoly_monomial_gt` return whether their *first*
argument is less, respectively greater, than the second.
:::

::: function
int mpoly_monomial_cmp1(ulong a, ulong b, ulong cmpmask)
int mpoly_monomial_gt1(ulong a, ulong b, ulong cmpmask)
int mpoly_monomial_ge1(ulong a, ulong b, ulong cmpmask)

Versions for one-word monomials: compare, $a > b$, $a \ge b$.
:::

::: function
int mpoly_monomial_cmp_nomask(const ulong * exp2, const ulong * exp3, slong N)
int mpoly_monomial_lt_nomask(const ulong * exp2, const ulong * exp3, slong N)
int mpoly_monomial_gt_nomask(const ulong * exp2, const ulong * exp3, slong N)

Compare the monomials as $N$-word unsigned integers, which is the
monomial ordering for `ORD_LEX` and `ORD_DEGLEX` (and is used for lex
orderings in internal code).
:::

::: function
int mpoly_monomial_cmp_nomask_extra(const ulong * exp2, const ulong * exp3, slong N, slong offset, ulong extra)
int mpoly_monomial_lt_nomask_extra(const ulong * exp2, const ulong * exp3, slong N, slong offset, ulong extra)
int mpoly_monomial_gt_nomask_extra(const ulong * exp2, const ulong * exp3, slong N, slong offset, ulong extra)

As the `_nomask` functions, comparing `exp2` with `exp3` plus `extra`
added to the word of index `offset`. Used for recursive representations
in which one variable is kept out of the packed monomial.
:::

::: function
int mpoly_monomial_cmp_general(ulong * Aexp, flint_bitcnt_t Abits, ulong * Bexp, flint_bitcnt_t Bbits, const mpoly_ctx_t mctx)

Compares two packed monomials with possibly different packing widths in
the ordering of `mctx`, returning $1$, $0$ or $-1$.
:::

::: function
int mpoly_monomials_cmp(const ulong * Aexps, flint_bitcnt_t Abits, const ulong * Bexps, flint_bitcnt_t Bbits, slong length, const mpoly_ctx_t mctx)

Compares two arrays of `length` monomials lexicographically term by term,
returning the comparison of the first pair that differs, or 0. The
packing widths may differ.
:::

## Generators {#generators}

::: function
slong mpoly_gen_index(slong v, const mpoly_ctx_t mctx)

Returns the index of the field that holds the exponent of the variable of
index `v`: `v` for reversed orderings and `nvars - 1 - v` otherwise.
:::

::: function
void mpoly_gen_fields_ui(ulong * exp, slong var, const mpoly_ctx_t mctx)
void mpoly_gen_fields_fmpz(fmpz * exp, slong var, const mpoly_ctx_t mctx)

Set the `nfields` unpacked fields of the generator of index `var`.
:::

::: function
flint_bitcnt_t mpoly_gen_bits_required(slong var, const mpoly_ctx_t mctx)

Returns the number of bits needed to pack a generator, which is
`MPOLY_MIN_BITS`.
:::

::: function
void mpoly_gen_monomial_sp(ulong * oneexp, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_gen_offset_shift_sp(slong * offset, slong * shift, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_gen_monomial_offset_shift_sp(ulong * mexp, slong * offset, slong * shift, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

For a single precision packing: set `oneexp` (or `mexp`) to the packed
generator of index `var`, and set `offset` and `shift` such that the
exponent of `var` in a packed monomial `e` is `(e[offset] >> shift)` with
the low `bits` bits kept.
:::

::: function
slong mpoly_gen_offset_mp(slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
slong mpoly_gen_monomial_offset_mp(ulong * mexp, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

For a multiprecision packing: return the word offset of the field of
`var`; the second version also sets `mexp` to the packed generator.
:::

::: function
int mpoly_is_gen(ulong * exps, slong var, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

Returns whether the packed monomial `exps` is the generator of index
`var`, or, if `var` is negative, whether it is any generator.
:::

::: function
ulong _mpoly_bidegree(const ulong * Aexps, flint_bitcnt_t Abits, const mpoly_ctx_t mctx)

For a lex ordering and single precision packing, returns the exponents of
the variables of index 0 and 1 in the monomial `Aexps`, packed into one
word with `pack_exp2`.
:::

::: function
void mpoly_compose_mat_gen(fmpz_mat_t M, const slong * c, const mpoly_ctx_t mctxB, const mpoly_ctx_t mctxAC)
void mpoly_compose_mat_fill_column(fmpz_mat_t M, const ulong * Cexp, flint_bitcnt_t Cbits, slong Bvar, const mpoly_ctx_t mctxB, const mpoly_ctx_t mctxAC)

Fill the matrix of a monomial substitution, used by the
`*_compose_*_gen` functions and the matrix method of composition. `M`
has `mctxAC->nfields + 1` rows and `mctxB->nfields` columns; the column
of a variable of $B$ holds the unpacked fields of its image in the ring
of $A$ and $C$, and a 1 in the last row means that the variable is mapped
to zero. The first function maps the variable of index $i$ to the
variable of index `c[i]` (or to zero if `c[i]` is negative or out of
range). The second fills the column of `Bvar` with the monomial
`(Cexp, Cbits)`, or marks it as zero if `Cexp` is `NULL`.
:::

## Arrays of monomials {#arrays}

These functions work on the monomial array of a polynomial: `len`
monomials of `N = mpoly_words_per_exp(bits, mctx)` words each, sorted in
decreasing order where stated.

::: function
int mpoly_monomial_exists(slong * index, const ulong * poly_exps, const ulong * exp, slong len, slong N, const ulong * cmpmask)

Searches the sorted array `(poly_exps, len)` for the monomial `exp` by
binary search. Returns 1 and sets `index` to its position if it is
present; otherwise returns 0 and sets `index` to the position where it
would have to be inserted to keep the array sorted, in $[0,
\mathit{len}]$.
:::

::: function
slong mpoly_monomial_index1_nomask(ulong * Aexps, slong Alen, ulong e)

Returns the index of the one-word monomial `e` in the array `(Aexps,
Alen)` sorted in decreasing order as unsigned integers, or $-1$ if it does
not occur.
:::

::: function
slong mpoly_monomial_index_ui(const ulong * Aexp, flint_bitcnt_t Abits, slong Alength, const ulong * exp, const mpoly_ctx_t mctx)
slong mpoly_monomial_index_pfmpz(const ulong * Aexp, flint_bitcnt_t Abits, slong Alength, fmpz * const * exp, const mpoly_ctx_t mctx)
slong mpoly_monomial_index_monomial(const ulong * Aexp, flint_bitcnt_t Abits, slong Alength, const ulong * Mexp, flint_bitcnt_t Mbits, const mpoly_ctx_t mctx)

Return the index of the monomial with user exponent vector `exp`
(respectively the packed monomial `Mexp` packed with `Mbits` bits) in
the sorted array `(Aexp, Alength)` packed with `Abits` bits, or $-1$ if it
does not occur.
:::

::: function
void mpoly_max_fields_ui_sp(ulong * max_fields, const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_max_fields_fmpz(fmpz * max_fields, const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_min_fields_ui_sp(ulong * min_fields, const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_min_fields_fmpz(fmpz * min_fields, const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

Set `max_fields` (`min_fields`) to the `nfields` unpacked field-wise
maxima (minima) of the monomials. The `_ui_sp` versions require single
precision packing. The minimum functions require $\mathit{len} > 0$; the
maximum of an empty array is zero.
:::

::: function
int mpoly_degrees_fit_si(const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

Returns whether the degree of the monomial array in each variable fits an
`slong`.
:::

::: function
void mpoly_degrees_si(slong * user_degs, const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_degrees_si_threaded(slong * user_degs, const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx, const thread_pool_handle * handles, slong num_handles)
void mpoly_degrees_ffmpz(fmpz * user_degs, const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_degrees_pfmpz(fmpz ** user_degs, const ulong * poly_exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

Set `user_degs` to the degree in each variable, all $-1$ if `len` is
zero. The `_si` versions assume that the degrees fit an `slong`. The
threaded version splits the array among the given thread pool handles.
:::

::: function
slong mpoly_degree_si(const ulong * poly_exps, slong len, flint_bitcnt_t bits, slong var, const mpoly_ctx_t mctx)
void mpoly_degree_fmpz(fmpz_t deg, const ulong * poly_exps, slong len, flint_bitcnt_t bits, slong var, const mpoly_ctx_t mctx)

Return or set the degree in the variable of index `var` ($-1$ if `len`
is zero). For a lex ordering and `var` equal to 0 only the first
monomial is read.
:::

::: function
int mpoly_total_degree_fits_si(const ulong * exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
slong mpoly_total_degree_si(const ulong * exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_total_degree_fmpz(fmpz_t totdeg, const ulong * exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
void mpoly_total_degree_fmpz_ref(fmpz_t totdeg, const ulong * exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

Compute the total degree, the maximum over the monomials of the sum of
the exponents, which is $-1$ for an empty array. For degree orderings
`mpoly_total_degree_fmpz` reads it from the degree field of the first
monomial; `mpoly_total_degree_fmpz_ref` always scans all monomials and
serves as a reference implementation for testing.
`mpoly_total_degree_si` assumes that the result fits.
:::

::: function
void mpoly_used_vars_or(int * used, const ulong * exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

For each variable of index $i$ that occurs with a positive exponent in
one of the monomials, sets `used[i]` to a nonzero value. Other entries
are left unchanged, so the results for several arrays can be combined.
:::

::: function
int mpoly_is_poly(const ulong * Aexps, slong Alen, flint_bitcnt_t Abits, slong var, const mpoly_ctx_t mctx)

Returns whether no variable other than the one of index `var` occurs in
the monomials.
:::

::: function
int mpoly_monomials_valid_test(ulong * exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
int mpoly_monomials_overflow_test(ulong * exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)
int mpoly_monomials_inorder_test(ulong * exps, slong len, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

Consistency checks used by the `assert_canonical` functions. The first
returns whether every total degree field (for degree orderings) equals the
sum of the exponents, the second whether some field has its guard bit
set, and the third whether the monomials are strictly decreasing.
:::

::: function
void mpoly_reverse(ulong * Aexp, const ulong * Bexp, slong len, slong N)

Sets `Aexp` to the monomials of `Bexp` in reverse order. In place
operation is allowed.
:::

::: function
void mpoly_monomials_deflation(fmpz * shift, fmpz * stride, const ulong * Aexps, flint_bitcnt_t Abits, slong Alength, const mpoly_ctx_t mctx)

For each variable $v$, with $S_v$ the set of exponents of $v$ in the
monomials, sets `shift[v]` to $\min S_v$ and `stride[v]` to $\gcd\{e -
\min S_v : e \in S_v\}$, which is zero if $S_v$ has one element. Every
monomial then has exponent $\mathit{shift}_v + k \cdot \mathit{stride}_v$
in $v$. For an empty array all shifts and strides are zero.
:::

::: function
void mpoly_monomials_deflate(ulong * Aexps, flint_bitcnt_t Abits, const ulong * Bexps, flint_bitcnt_t Bbits, slong Blength, const fmpz * shift, const fmpz * stride, const mpoly_ctx_t mctx)
void mpoly_monomials_inflate(ulong * Aexps, flint_bitcnt_t Abits, const ulong * Bexps, flint_bitcnt_t Bbits, slong Blength, const fmpz * shift, const fmpz * stride, const mpoly_ctx_t mctx)

Apply $e \mapsto (e - \mathit{shift}_v)/\mathit{stride}_v$ (deflate) or
$e \mapsto \mathit{shift}_v + \mathit{stride}_v e$ (inflate) to the
exponent of each variable $v$ of the monomials `Bexps`, writing the
result to `Aexps` with `Abits` bits, which must be large enough. When
deflating, a zero stride maps the exponent to zero, and every $e -
\mathit{shift}_v$ must be divisible by the stride. The order of the
monomials may change.
:::

::: function
void _mpoly_gen_shift_right(ulong * Aexp, flint_bitcnt_t Abits, slong Alength, slong var, ulong amount, const mpoly_ctx_t mctx)
void _mpoly_gen_shift_right_fmpz(ulong * Aexp, flint_bitcnt_t Abits, slong Alength, slong var, const fmpz_t amount, const mpoly_ctx_t mctx)
void _mpoly_gen_shift_left(ulong * Aexp, flint_bitcnt_t Abits, slong Alength, slong var, ulong amount, const mpoly_ctx_t mctx)

Divide (shift right) or multiply (shift left) every monomial in place by
$x_{\mathit{var}}^{\mathit{amount}}$, without changing the packing width.
The division must be exact and the product must fit; `_mpoly_gen_shift_right`
and `_mpoly_gen_shift_left` require single precision packing.
:::

::: function
void mpoly_monomials_shift_right_ui(ulong * Aexps, flint_bitcnt_t Abits, slong Alength, const ulong * user_exps, const mpoly_ctx_t mctx)
void mpoly_monomials_shift_right_ffmpz(ulong * Aexps, flint_bitcnt_t Abits, slong Alength, const fmpz * user_exps, const mpoly_ctx_t mctx)

Divide every monomial in place by the monomial with user exponent vector
`user_exps`, which must divide all of them.
:::

::: function
void mpoly_remove_var_powers(fmpz * var_powers, ulong * Aexps, flint_bitcnt_t Abits, slong Alen, const mpoly_ctx_t mctx)

Sets `var_powers` to the user exponent vector of the greatest common
divisor of the monomials ($\mathit{Alen} > 0$) and divides every monomial
by it in place.
:::

::: function
void mpoly1_fill_marks(ulong ** Dcoeffs, slong * Dlen, slong * Dalloc, const ulong * Aexps, slong Alen, flint_bitcnt_t Abits, const mpoly_ctx_t mctx)
void mpoly2_fill_marks(ulong ** Dcoeffs, slong * Dlen, slong * Dalloc, const ulong * Aexps, slong Alen, flint_bitcnt_t Abits, const mpoly_ctx_t mctx)

For a lex ordering and single precision packing, partition the sorted
monomials into runs with the same exponent of the variable of index 0
(respectively the same exponents of the variables of index 0 and 1). The
array `(*Dcoeffs, *Dlen)` receives the start index of each run, and is
reallocated as needed (`*Dalloc` is its allocation) so that
`(*Dcoeffs)[*Dlen]` equals `Alen`.
:::

::: function
void mpoly_to_mpolyl_perm_deflate(ulong * Aexps, flint_bitcnt_t Abits, const mpoly_ctx_t Actx, ulong * Bexps, flint_bitcnt_t Bbits, const mpoly_ctx_t Bctx, slong length, const slong * perm, const ulong * shift, const ulong * stride)
void mpoly_from_mpolyl_perm_inflate(ulong * Bexps, flint_bitcnt_t Bbits, const mpoly_ctx_t Bctx, ulong * Aexps, flint_bitcnt_t Abits, const mpoly_ctx_t Actx, slong length, const slong * perm, const ulong * shift, const ulong * stride)

Convert monomials between a context `Bctx` with $n$ variables and a lex
context `Actx` with $m \le n$ variables, both with single precision
packing. Deflation sets the exponent of variable $k$ of $A$ to
$(e_l - \mathit{shift}_l)/\mathit{stride}_l$ where $l = \mathit{perm}[k]$
and $e_l$ is the exponent of variable $l$ of $B$; inflation is the
inverse, setting the exponents of variables of $B$ not in the image of
`perm` to their shift. The GCD and factoring code use these to move to a
lex ordering in the variables that actually occur. The order of the
monomials is not maintained; the caller sorts afterwards.
:::

## Packing into a single word {#tight}

The array (dense) multiplication and division methods map each monomial
to an integer index in a dense array. With bounds $b_0, b_1, \ldots$ on
the exponents of the fields, a monomial with fields $e_0, e_1, \ldots$ is
encoded in the mixed radix representation $e_0 + b_0 (e_1 + b_1 (e_2 +
\cdots))$, called *tight* or *factorial* packing.

::: function
void mpoly_pack_monomials_tight(ulong * exp1, const ulong * exp2, slong len, const slong * mults, slong num, slong bits)
void mpoly_unpack_monomials_tight(ulong * e1, ulong * e2, slong len, slong * mults, slong num, slong bits)

Convert `len` one-word monomials between the usual packing with `bits`
bits per field and the tight packing with the bases `mults` (the first
entry corresponding to the least significant field). Only the least
significant `num` fields are converted; each must be less than its base.
:::

::: function
int mpoly_monomial_divides_tight(slong e1, slong e2, slong * prods, slong num)

Returns whether the tightly packed monomial `e2` divides `e1`. The array
`(prods, num + 1)` holds the partial products $1, b_0, b_0 b_1, \ldots$
of the bases.
:::

::: function
void mpoly_max_degrees_tight(slong * max_exp, ulong * exps, slong len, slong * prods, slong num)

Sets `max_exp` to the maximum of each of the `num` fields over the
tightly packed monomials `(exps, len)`, with `prods` as for
`mpoly_monomial_divides_tight`.
:::

::: function
void mpoly_main_variable_split_LEX(slong * ind, ulong * pexp, const ulong * Aexp, slong l1, slong Alen, const ulong * mults, slong num, slong Abits)
void mpoly_main_variable_split_DEG(slong * ind, ulong * pexp, const ulong * Aexp, slong l1, slong Alen, ulong deg, slong num, slong Abits)

Split one-word monomials `(Aexp, Alen)` by the value of their most
significant field (the main variable for lex, the total degree for degree
orderings), which takes at most `l1` values. The array `ind` of length
$\mathit{l1} + 1$ receives the start of each chunk, and `pexp` the
remaining `num` fields of each monomial in tight packing with the bases
`mults` (for lex) or all equal to `deg` (for degree orderings, where the
lowest field is dropped because the total degree determines it). Used by
the array multiplication.
:::

::: function
void mpoly_main_variable_terms1(slong * i1, slong * n1, const ulong * exp1, slong l1, slong len1, slong k, slong num, slong bits)

Given one-word monomials `(exp1, len1)` with fields of `bits` bits, sets
the `l1` start indices `i1` and lengths `n1` (which may be zero) of the
chunks of monomials with exponent $\mathit{l1} - 1, \mathit{l1} - 2,
\ldots, 0$ in the field of index $k - 1$, which must be the most
significant nonzero field. The fields are not masked.
:::

::: function
ulong pack_exp2(ulong e0, ulong e1)
ulong pack_exp3(ulong e0, ulong e1, ulong e2)
ulong extract_exp(ulong e, int idx, int nvars)

Pack two (three) small exponents into one word with fields of
$\text{FLINT\_BITS}/2$ ($\text{FLINT\_BITS}/3$) bits, `e0` in the most
significant field, so that comparing the words compares the exponent
tuples lexicographically. `extract_exp` returns the field of index `idx`
(counted from the least significant) of a word packed with `nvars` equal
fields. Used for the bivariate and trivariate base cases of the GCD and
factoring code.
:::

## Heaps {#heaps}

Johnson's heap method for multiplication and the Monagan–Pearce division
algorithms keep the candidate terms of the result in a binary max-heap
ordered by monomial. Nodes with equal monomials are *chained* instead of
being inserted separately, so that the heap holds each monomial at most
once. See [](algorithms/multivariate/multiplication).

::: type
mpoly_heap_t
mpoly_nheap_t

A chain node: the indices `i` and `j` of the two input terms whose
product the node represents, and a pointer `next` to the next node with
the same monomial. `mpoly_nheap_t` has an extra field `p` for the index of
the divisor in the ideal division functions.
:::

::: type
mpoly_heap1_s
mpoly_heap_s

A heap entry: a monomial `exp` (one word, respectively a pointer to $N$
words) and a pointer `next` to the first node of the chain of products
with this monomial.
:::

::: macro
HEAP_LEFT(i)
HEAP_RIGHT(i)
HEAP_PARENT(i)
HEAP_ASSIGN(h, c1, c2)

Index arithmetic for a binary heap stored from index 1 ($2i$, $2i + 1$,
$\lfloor i/2 \rfloor$), and assignment of the monomial `c1` and chain
`c2` to the heap entry `h`.
:::

::: function
int _mpoly_heap_insert(mpoly_heap_s * heap, ulong * exp, void * x, slong * next_loc, slong * heap_len, slong N, const ulong * cmpmask)
void _mpoly_heap_insert1(mpoly_heap1_s * heap, ulong exp, void * x, slong * next_loc, slong * heap_len, ulong maskhi)

Insert the chain node `x` with monomial `exp` into the heap. The entry of
index 0 is unused, so `*heap_len` is one more than the number of entries.
If an entry with the same monomial is found at the top, at index
`*next_loc`, or on the path from the insertion point to the root, `x` is
chained to it (and `*next_loc` remembers its position); otherwise a new
entry is created and `*heap_len` is incremented. `_mpoly_heap_insert`
returns 1 if a new entry was created and 0 if the node was chained.
The one-word version takes `maskhi = cmpmask[0]`. Initialise `*next_loc`
to a value larger than any heap length.
:::

::: function
void * _mpoly_heap_pop(mpoly_heap_s * heap, slong * heap_len, slong N, const ulong * cmpmask)
void * _mpoly_heap_pop1(mpoly_heap1_s * heap, slong * heap_len, ulong maskhi)

Remove the entry with the largest monomial and return its chain (cast to
`void *`), decrementing `*heap_len`. The heap must not be empty.
:::

::: function
void mpoly_search_monomials(slong ** e_ind, ulong * e, slong * e_score, slong * t1, slong * t2, slong * t3, slong lower, slong upper, const ulong * a, slong a_len, const ulong * b, slong b_len, slong N, const ulong * cmpmask)

Splits the products of two sorted monomial arrays for parallel
multiplication. Finds a monomial `e` such that the number of products
$a_i b_j \le e$ is between `lower` and `upper` (or as close as possible),
stores this number in `e_score`, and sets `*e_ind` to one of the work
arrays `t1`, `t2`, `t3` (each of length `a_len`): entry $i$ of it is the
first index $j$ with $a_i b_j \le e$. Assumes that $1$ is the smallest
monomial.
:::

::: function
slong mpoly_divide_threads(slong n, double la, double lb)

Returns $m$ with $0 \le m < n$ such that $(m + 1)/(n - m)$ is close to
$\mathit{la}/\mathit{lb}$: the number of threads (besides the calling
one) to give to the first of two independent tasks of sizes `la` and
`lb` when $n$ threads are available.
:::

## Red-black trees {#rbtree}

Maps from `ulong` or `fmpz` keys to fixed-size data blocks, used to
collect terms by exponent.

::: type
mpoly_rbnode_ui_struct
mpoly_rbtree_ui_struct
mpoly_rbtree_ui_t
mpoly_rbnode_fmpz_struct
mpoly_rbtree_fmpz_struct
mpoly_rbtree_fmpz_t

A red-black tree with keys of type `ulong` or `fmpz_t`. Nodes are stored
in an array and refer to each other by index; the data blocks, each of
`data_size` bytes, are stored in a separate array in the same order.
:::

::: function
void mpoly_rbtree_ui_init(mpoly_rbtree_ui_t T, slong data_size)
void mpoly_rbtree_fmpz_init(mpoly_rbtree_fmpz_t T, slong data_size)
void mpoly_rbtree_ui_clear(mpoly_rbtree_ui_t T)
void mpoly_rbtree_fmpz_clear(mpoly_rbtree_fmpz_t T)

Initialise an empty tree with data blocks of `data_size` bytes, or free
it. Clearing does not clear the data blocks themselves.
:::

::: function
void * mpoly_rbtree_ui_lookup(mpoly_rbtree_ui_t T, int * its_new, ulong key)
void * mpoly_rbtree_fmpz_lookup(mpoly_rbtree_fmpz_t T, int * its_new, const fmpz_t key)

Return a pointer to the data block for `key`, inserting a new node if the
key is not present; `*its_new` is set to whether this happened (the new
block is uninitialised). The pointer is valid until the next insertion.
:::

::: function
slong mpoly_rbtree_ui_head(const mpoly_rbtree_ui_t T)
slong mpoly_rbtree_fmpz_head(const mpoly_rbtree_fmpz_t T)

Return the index of the root node, from which the tree can be traversed
through the `left` and `right` fields of the nodes.
:::

## Helpers for GCD computations {#gcd}

The GCD functions of all coefficient rings first analyse the exponents of
both inputs and record the results in an `mpoly_gcd_info_t`, which is
used to dispose of easy cases and to choose an algorithm and a variable
order. See [](algorithms/multivariate/gcd).

::: macro
MPOLY_GCD_USE_HENSEL
MPOLY_GCD_USE_BROWN
MPOLY_GCD_USE_ZIPPEL
MPOLY_GCD_USE_ZIPPEL2
MPOLY_GCD_USE_PRS
MPOLY_GCD_USE_ALL

Bit flags (1, 2, 4, 8, 16 and their union 31) selecting the GCD
algorithms that may be used: Hensel lifting, Brown's dense algorithm,
Zippel's sparse algorithm, the sparse algorithm based on Berlekamp–Massey
interpolation, and pseudo-remainder sequences.
:::

::: type
mpoly_gcd_info_struct
mpoly_gcd_info_t

Data about a GCD problem $\gcd(A, B)$ in $n$ variables: for each variable
the maximal and minimal exponent in $A$ and $B$ and the number of terms
attaining them, the common stride of the exponents, the degrees after
deflation, bounds and term count estimates for the GCD, the number
`mvars` of variables that actually occur, the densities of $A$ and $B$,
for each algorithm an estimated cost and a variable permutation, and the
bit set `can_use` of algorithms that may be applied.
:::

::: function
void mpoly_gcd_info_init(mpoly_gcd_info_t Iv, slong nvars)
void mpoly_gcd_info_clear(mpoly_gcd_info_t Iv)

Initialise for $\mathit{nvars} > 0$ variables, and free.
:::

::: function
void mpoly_gcd_info_limits(ulong * Amax_exp, ulong * Amin_exp, slong * Amax_exp_count, slong * Amin_exp_count, const ulong * Aexps, flint_bitcnt_t Abits, slong Alength, const mpoly_ctx_t mctx)

For each variable, sets the maximal and minimal exponent in the monomials
`(Aexps, Alength)` ($\mathit{Alength} > 0$, single precision packing) and
the number of monomials attaining each.
:::

::: function
void mpoly_gcd_info_stride(ulong * strides, const ulong * Aexps, flint_bitcnt_t Abits, slong Alength, const ulong * Amax_exp, const ulong * Amin_exp, const ulong * Bexps, flint_bitcnt_t Bbits, slong Blength, const ulong * Bmax_exp, const ulong * Bmin_exp, const mpoly_ctx_t mctx)

For each variable $v$, sets `strides[v]` to the GCD of all $e -
\min_A$ for exponents $e$ of $v$ in $A$ and all $e - \min_B$ for
exponents in $B$, given the limits computed by `mpoly_gcd_info_limits`.
Deflating by these strides does not change the GCD.
:::

::: function
void mpoly_gcd_info_set_perm(mpoly_gcd_info_t Iv, slong Alength, slong Blength, const mpoly_ctx_t mctx)

Sets `Iv->mvars` to the number of variables that occur in $A$ (after the
monomial content is removed, they occur in both), initialises the
variable permutations of all algorithms to these variables in their
natural order, sets the densities (the number of terms divided by the
product of the deflated degrees plus one) and clears `can_use`.
:::

::: function
slong mpoly_gcd_info_get_brown_upper_limit(const mpoly_gcd_info_t Iv, slong var, slong bound)

Returns the number of interpolation steps in the variable of index `var`
of the Brown permutation after which Brown's algorithm should stop
attempting early divisibility tests, or 0 if the degree bounds of `Iv` are
not reliable.
:::

::: function
void mpoly_gcd_info_measure_hensel(mpoly_gcd_info_t Iv, slong Alength, slong Blength, const mpoly_ctx_t mctx)
void mpoly_gcd_info_measure_brown(mpoly_gcd_info_t Iv, slong Alength, slong Blength, const mpoly_ctx_t mctx)
void mpoly_gcd_info_measure_bma(mpoly_gcd_info_t Iv, slong Alength, slong Blength, const mpoly_ctx_t mctx)
void mpoly_gcd_info_measure_zippel(mpoly_gcd_info_t Iv, slong Alength, slong Blength, const mpoly_ctx_t mctx)
void mpoly_gcd_info_measure_zippel2(mpoly_gcd_info_t Iv, slong Alength, slong Blength, const mpoly_ctx_t mctx)

Decide whether an algorithm is applicable and, if so, set its bit in
`Iv->can_use`, choose its variable order and estimate its cost. Brown and
Hensel need at least two variables whose dense representations fit, and
their cost is modelled from the products of the degrees of $G$, $A/G$
and $B/G$; Zippel chooses as main variable one with few leading and
trailing terms; the Berlekamp–Massey based method (`measure_bma`, used
for `MPOLY_GCD_USE_ZIPPEL2`) needs at least three variables and chooses
two main variables of small degree that fit in half a word.
`measure_zippel2` is the analogue for the Zippel2 variant used by the
`nmod` and `fq_nmod` modules.
:::

::: function
int mpoly_monomial_cofactors(fmpz * Abarexps, fmpz * Bbarexps, const ulong * Aexps, flint_bitcnt_t Abits, const ulong * Bexps, flint_bitcnt_t Bbits, slong length, const mpoly_ctx_t mctx)

Tests whether two monomial arrays of the same length differ by a monomial
factor term by term up to a common monomial, that is, whether $A = a
\cdot H$ and $B = b \cdot H$ term by term for monomials $a$, $b$ with
$\gcd(a, b) = 1$. If so, sets `Abarexps` and `Bbarexps` to the user
exponent vectors of $a$ and $b$ and returns 1; otherwise returns 0. Used
to detect GCDs whose cofactors are monomials.
:::

## Helpers for factorisation {#factor}

::: macro
MPOLY_FACTOR_USE_ZAS
MPOLY_FACTOR_USE_WANG
MPOLY_FACTOR_USE_ZIP
MPOLY_FACTOR_USE_ALL

Bit flags (1, 2, 4 and their union 7) selecting the factoring algorithms:
Zassenhaus-style recombination, Wang's algorithm with dense Hensel
lifting, and Zippel's sparse Hensel lifting. See
[](algorithms/multivariate/factorisation).
:::

::: function
int mpoly_is_proved_not_square(const ulong * Aexps, slong Alen, flint_bitcnt_t Abits, slong N, ulong * t)

Tries to prove that a polynomial with monomials `(Aexps, Alen)` is not a
square, using `t` ($N$ words) as scratch space. Returns 1 if the
field-wise maximum of the monomials is not a square monomial, that is,
some variable (or, for degree orderings, the total degree) has odd
maximal degree, since the degrees of a square are even. Returns 0 if this
test is inconclusive, and for an empty array.
:::

::: type
mpoly_compression_struct
mpoly_compression_t

Describes a change of variables that reduces the number of variables of
a polynomial: an invertible integer matrix `umat`, a shift `deltas`, the
degrees `degs` in the new variables, the exponent matrix `exps` in the
new variables, the numbers `nvars` and `mvars` of old and new variables,
and flags `is_trivial` (the matrix is the identity and the ordering lex),
`is_perm` (the matrix is a permutation) and `is_irred` (the Newton
polytope proves the polynomial irreducible).
:::

::: function
void mpoly_compression_init(mpoly_compression_t M)
void mpoly_compression_clear(mpoly_compression_t M)

Initialise and free.
:::

::: function
void mpoly_compression_set(mpoly_compression_t M, const ulong * Aexps, flint_bitcnt_t Abits, slong Alen, const mpoly_ctx_t mctx)

Computes a compression of the polynomial with monomials `(Aexps, Alen)`
(single precision packing): a unimodular change of exponents $e = V e' +
D$ after which the exponents $e'$ are nonnegative and the polynomial
depends on as few variables with as small degrees as the heuristic finds,
and tests irreducibility from the Newton polytope.
:::

::: function
slong _mpoly_compress_exps(slong * V, slong * D, slong * deg, slong * S, slong n, slong l)

The core of `mpoly_compression_set`: given the $l \times n$ exponent
matrix `S` (rows are monomials, $l > 1$), finds a unimodular $n \times n$
matrix `V`, a vector `D` and new exponents (overwriting `S`) with $V
S'_i + D = S_i$ for every row, $\min_i S'_{ij} = 0$ and `deg[j]` one more
than the maximum of column $j$, by greedy unimodular column operations
that reduce the product of the degrees. Returns the number of columns of
$S'$ that are not identically zero, which come first.
:::

::: function
int mpoly_test_irreducible(ulong * Aexps, flint_bitcnt_t Abits, slong Alen, const mpoly_ctx_t ctx)
int _mpoly_test_irreducible(slong * Aexps, slong stride, slong Alen, slong nvars, flint_rand_t state, slong tries_left)

Try to prove that a polynomial with the given monomials is absolutely
irreducible from its Newton polytope alone, independently of the
coefficients, following Gao and Lauder [@GaoLau2001]. The exponents are
projected onto random planes, and the projection of the polytope is tested
for integral indecomposability, at most `tries_left` times. Return 1 if
irreducibility was proved and 0 if not; a return of 0 says nothing. The
underscore version takes the exponent matrix with row stride `stride`;
the public version returns 0 for multiprecision packing or fewer than two
terms.
:::

::: function
void mpoly_bivar_cld_bounds(slong * l, slong n)

For $f \in F[y][x]$ of degree $n - 1$ in $x$, transforms the input
$l_i = 1 + \deg_y(\operatorname{coeff}(f, x^i))$ for $0 \le i < n$ into
bounds $l_i > \deg_y \operatorname{coeff}(f g_x/g, x^i)$ valid for every
divisor $g$ of $f$, using the upper convex hull of the points $(i, l_i)$.
These logarithmic derivative bounds are used by the bivariate
factoring code over finite fields.
:::

::: function
void _slong_array_fit_length(slong ** array, slong * alloc, slong len)

Reallocates the `slong` array `*array` with allocation `*alloc` to hold at
least `len` entries, growing geometrically.
:::

## Generic univariates and parsing {#generic}

A few algorithms are written once for all coefficient rings, which are
passed as a table of function pointers.

::: type
mpoly_void_ring_t

A ring described by function pointers: `elem_size` (the size of an
element in bytes), the context pointer `ctx`, and functions `init`,
`clear`, `is_zero`, `zero`, `one`, `set_fmpz`, `set`, `swap`, `neg`,
`add`, `sub`, `mul_fmpz`, `mul`, `divexact`, `divides`, `pow_fmpz` and
`length`. Each polynomial module provides a function such as
`mpoly_void_ring_init_fmpz_mpoly_ctx` that fills it for its polynomial
ring.
:::

::: function
void * mpoly_void_ring_elem_init(mpoly_void_ring_t R)
void mpoly_void_ring_elem_clear(void * a, mpoly_void_ring_t R)

Allocate and initialise an element of $R$, and clear and free it.
:::

::: type
mpoly_univar_struct
mpoly_univar_t

A sparse univariate polynomial over a ring given by an
`mpoly_void_ring_t`: `length` terms with coefficients stored contiguously
in `coeffs` (each `elem_size` bytes) and `fmpz` exponents in `exps`, in
decreasing order of exponent.
:::

::: function
void mpoly_univar_init(mpoly_univar_t A, mpoly_void_ring_t R)
void mpoly_univar_init2(mpoly_univar_t A, slong len, mpoly_void_ring_t R)
void mpoly_univar_clear(mpoly_univar_t A, mpoly_void_ring_t R)
void mpoly_univar_swap(mpoly_univar_t A, mpoly_univar_t B)
void mpoly_univar_fit_length(mpoly_univar_t A, slong len, mpoly_void_ring_t R)

Memory management: initialise (with room for `len` terms), clear, swap,
and ensure room for `len` terms.
:::

::: function
int mpoly_univar_pseudo_gcd_ducos(mpoly_univar_t G, mpoly_univar_t B, mpoly_univar_t A, mpoly_void_ring_t R)

Sets $G$ to the last nonzero subresultant of $B$ and $A$, computed with
Ducos' algorithm [@Duc2000], which divides out the predictable factors
of the subresultant sequence early. Requires $\deg B \ge \deg A \ge 0$
and nonzero inputs; $A$ and $B$ are destroyed. Returns 0 if a ring
operation failed (for example an exponent overflow in a power), 1
otherwise.
:::

::: function
int mpoly_univar_resultant(void * r, mpoly_univar_t fx, mpoly_univar_t gx, mpoly_void_ring_t R)
int mpoly_univar_discriminant(void * d, mpoly_univar_t fx, mpoly_void_ring_t R)

Set $r$ to the resultant of $f$ and $g$, respectively $d$ to the
discriminant of $f$, using the subresultant algorithm; the inputs are
destroyed. The resultant of a zero polynomial is zero, the discriminant
of a constant $a$ is $1/a^2$ (so it can fail), and the discriminant of a
linear polynomial is 1. Return 1 for success and 0 if a ring operation
failed.
:::

::: type
string_with_length_struct
mpoly_parse_struct
mpoly_parse_t

An operator precedence parser for polynomial expressions over a ring given
by an `mpoly_void_ring_t`, with a stack of operators and expressions and a
table of named terminals (variables and other named constants).
:::

::: function
void mpoly_parse_init(mpoly_parse_t E)
void mpoly_parse_clear(mpoly_parse_t E)

Initialise and free a parser. The ring `E->R` must be set before the
first use.
:::

::: function
void mpoly_parse_add_terminal(mpoly_parse_t E, const char * s, const void * v)

Adds the name `s` with value `v` (an element of `E->R`, copied).
Terminals are kept sorted by decreasing length, so that the longest
matching name is used.
:::

::: function
int mpoly_parse_parse(mpoly_parse_t E, void * res, const char * s, slong len)

Parses the string `(s, len)` and sets `res` to its value, returning 0 for
success and $-1$ for failure; on failure `res` holds an unspecified value.
The grammar is described at `fmpz_mpoly_set_str_pretty`: nonnegative
integer constants, terminals, `+`, `-` (binary and unary), `*`, `/` (exact
division), `^` followed immediately by a nonnegative integer, parentheses
and spaces.
:::
