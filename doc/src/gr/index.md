---
title: gr.h
header: gr.h
description: Generic rings: context objects, elements, status codes and the generic interface to arithmetic in any ring.
category: Generic rings
order: 1
---

The `gr` module is FLINT's interface for generic programming. A
mathematical structure $R$, usually a ring, is represented at run time by
a *context object* of type `gr_ctx_t`, which records the size of an
element and a table of methods; elements are untyped blocks of memory
passed by pointer. A function written against this interface, such as
`exp_partial_sum` in [](gr/getting-started), runs unchanged over $\Z$,
$\Z/n\Z$, finite fields, number fields, real and complex balls, exact
algebraic numbers, polynomials, matrices, power series and other rings
built from these. The generic algorithms of `gr_vec`, `gr_poly`,
`gr_mat`, `gr_series` and `gr_mpoly` are written this way, and many
type-specific FLINT functions are thin wrappers around them. Almost every
function returns a status code that reports whether the result is
mathematically undefined (`GR_DOMAIN`) or could not be computed
(`GR_UNABLE`), and predicates return a three-valued `truth_t`, so that
inexact and partially decidable rings are handled rigorously.

The constructors of the built-in rings are listed in
[`gr_domains`](gr_domains), the default implementations of methods in
[`gr_generic`](gr_generic), special functions in
[`gr_special`](gr_special) and vector operations in [`gr_vec`](gr_vec). Start with the guide [](gr/getting-started); the
guide [](gr/implementing) explains how to define a new ring. The design
of the framework (method tables, inexact elements, error handling and
testing) is described in [](algorithms/generic/overview) and in
[@Joh2025gr].

::: summary algorithms/generic/overview :::

## Types {#types}

::: type
gr_ctx_struct
gr_ctx_t

A context object representing a mathematical structure $R$. A
`gr_ctx_t` is an array of length one of `gr_ctx_struct`, so that it is
passed by reference. The structure has the fields

* `data`: `GR_CTX_STRUCT_DATA_BYTES` (six words) of storage for the
  parameters of the ring, for example a modulus, a precision, or a
  pointer to a base ring or to heap-allocated data;
* `which_ring`: a code of type `gr_which_structure` identifying the
  implementation;
* `sizeof_elem`: the size of an element in bytes;
* `methods`: a pointer to the method table;
* `size_limit`: a bound on the size of elements (for example a number of
  terms), beyond which some rings refuse operations with `GR_UNABLE`;
  `WORD_MAX` for most rings.

Contexts are not passed as `const`, so that a ring can keep mutable data
such as caches or a table of discovered zero divisors.
:::

::: type
gr_ctx_ptr

A pointer to a context object (`void *`).
:::

::: type
gr_ptr
gr_srcptr

A pointer to an element or to an array of contiguous elements, and a
pointer to a read-only element or array. These are `void *` and
`const void *`, so a pointer to any element type can be passed.
:::

There is no type for a single generic element, since its size is only
known at run time. Single elements are allocated on the stack with the
`GR_TMP_INIT` macros, on the heap with `gr_heap_init`, or declared with
their concrete type when it is known: an element of the context created
by `gr_ctx_init_fmpz` is an `fmpz`, an element of `gr_ctx_init_nmod` is
a `ulong`, an element of `gr_ctx_init_complex_acb` is an `acb_struct`.

::: macro
GR_ENTRY(vec, i, size)

The address of entry $i$ of the array `vec` whose elements have `size`
bytes, that is, `(char *) vec + i * size` as a `void *`. Typically
`size` is `ctx->sizeof_elem`.
:::

::: function
slong gr_ctx_sizeof_elem(gr_ctx_t ctx)

Returns `ctx->sizeof_elem`, the size in bytes of an element of `ctx`.
:::

::: function
slong gr_ctx_sizeof_ctx(void)

Returns `sizeof(gr_ctx_struct)`, for interfaces that allocate context
objects without access to the C type.
:::

::: function
void * gr_ctx_data_ptr(gr_ctx_t ctx)
void * gr_ctx_data_as_ptr(gr_ctx_t ctx)

Return a pointer to the `data` field of `ctx`, respectively the pointer
stored at the start of the `data` field (the value of
`GR_CTX_DATA_AS_PTR(ctx)`), for rings whose parameters are kept in a
separately allocated structure.
:::

::: macro
GR_CTX_DATA_AS_PTR(ctx)

The first word of the `data` field of `ctx`, read as a `void *` lvalue.
:::

## Status codes {#status}

Functions that can fail return an `int` combining the following flags.
When the result is not `GR_SUCCESS`, the values of the output variables
are unspecified (but valid, so they can be cleared or overwritten).
Statuses are combined with `|`, so a sequence of operations can be
checked once at the end:

```c
status = GR_SUCCESS;
status |= gr_add(res, a, b, ctx);
status |= gr_pow_ui(res, res, 2, ctx);
```

The two failure flags have distinct meanings, explained with examples in
[](algorithms/generic/errors). Functions that return a status code are
declared with the `WARN_UNUSED_RESULT` attribute, so that compilers warn
when a status is ignored.

::: macro
GR_SUCCESS

The value 0: the operation succeeded and the output is a correct element
of the target structure (for inexact representations, a correct
enclosure).
:::

::: macro
GR_DOMAIN

The value 1: the result is not defined in the target structure, for
example division by zero, $1/2$ in $\Z$, $\sqrt{-1}$ in $\mathbb{R}$,
the inverse of a zero divisor, or a conversion of an element that does
not belong to the target (such as `gr_get_ui` of a negative integer).
:::

::: macro
GR_UNABLE

The value 2: the result may exist but the implementation could not
compute it: the operation is not implemented for this ring, the result
would be too large, the inputs are inexact and an exact decision is
needed, or the computation would take too long. When this flag is set,
`GR_DOMAIN` may also apply without being known; `GR_DOMAIN | GR_UNABLE`
means "unable".
:::

::: macro
GR_TEST_FAIL

The value 4, used only by test code to report a wrong result.
:::

::: macro
GR_MUST_SUCCEED(expr)

Evaluates `expr` and throws an exception with `flint_throw` (abort by
default) if the result is not `GR_SUCCESS`. Use it when failure is
impossible or fatal.
:::

::: macro
GR_IGNORE(expr)

Evaluates `expr` and discards the status, avoiding a warning about an
unused result.
:::

## Truth values {#truth}

Predicates return a `truth_t` instead of a C boolean, so that they can
answer "unknown" when the question cannot be decided, for example
whether two overlapping balls are equal.

::: type
truth_t

An enumeration with the values `T_TRUE` (0), `T_FALSE` (1) and
`T_UNKNOWN` (2). Since `T_TRUE` is zero, a `truth_t` must never be used
as a C condition: write `t == T_TRUE` or `t != T_FALSE` according to
whether the unknown case should be excluded or included.
:::

::: macro
T_TRUE
T_FALSE
T_UNKNOWN

The three truth values.
:::

::: function
truth_t truth_and(truth_t x, truth_t y)
truth_t truth_or(truth_t x, truth_t y)
truth_t truth_not(truth_t x)

Kleene's three-valued conjunction, disjunction and negation: for
example `truth_and` returns `T_FALSE` if either argument is `T_FALSE`
and `T_UNKNOWN` if neither is `T_FALSE` and one is `T_UNKNOWN`.
:::

::: function
void truth_println(truth_t x)

Prints `T_TRUE`, `T_FALSE` or `T_UNKNOWN` followed by a newline.
:::

::: function
truth_t gr_in_domain(int status)

Converts a status into the truth of "the result is defined": `T_TRUE`
for `GR_SUCCESS`, `T_FALSE` when `GR_DOMAIN` is set and `GR_UNABLE` is
not, and `T_UNKNOWN` otherwise.
:::

::: function
int gr_check(truth_t t)

Converts a truth value into a status: `T_TRUE` to `GR_SUCCESS`,
`T_FALSE` to `GR_DOMAIN` and `T_UNKNOWN` to `GR_UNABLE`.
:::

## Contexts {#contexts}

::: function
void gr_ctx_clear(gr_ctx_t ctx)

Clears the context object, freeing any memory it owns. All elements of
the ring must be cleared before the context; a ring built on top of a
base ring (for example a polynomial ring) must be cleared before its
base ring.
:::

::: function
int gr_ctx_write(gr_stream_t out, gr_ctx_t ctx)
int gr_ctx_print(gr_ctx_t ctx)
int gr_ctx_println(gr_ctx_t ctx)
int gr_ctx_get_str(char ** s, gr_ctx_t ctx)

Write a description of the ring, such as `Integers mod 7 (_gr_nmod)`,
to a stream, print it to standard output (with a newline for
`gr_ctx_println`), or set `*s` to a newly allocated string that must be
freed with `flint_free`. The format specifier `%{gr_ctx}` of
`flint_printf` does the same ([](flint/printf)).
:::

::: function
gr_ptr gr_ctx_base(gr_ctx_t ctx)

Returns a pointer to the context of the base ring if `ctx` is
constructed over one (polynomials, matrices, series, fractions,
quotients); otherwise returns `NULL`.
:::

::: function
int gr_ctx_set_gen_name(gr_ctx_t ctx, const char * s)
int gr_ctx_set_gen_names(gr_ctx_t ctx, const char ** s)

Set the name of the generator, for rings with one generator (polynomial
rings, finite fields, series), or the names of the generators of a
multivariate ring. The names are used for printing and parsing. Return
`GR_UNABLE` if the ring has no named generators.
:::

::: function
int gr_ctx_ngens(slong * ngens, gr_ctx_t ctx)
int gr_ctx_gen_name(char ** name, slong i, gr_ctx_t ctx)

Set `*ngens` to the number of generators (0 by default), respectively
set `*name` to a newly allocated copy of the name of generator $i$,
which must be freed with `flint_free`.
:::

::: function
void gr_ctx_uninitialized(gr_ctx_t ctx)

Sets `ctx` to a placeholder context that can be printed and cleared but
supports no other operation; initialising an element throws an
exception. Wrapper code uses it for a context whose constructor failed
(for example `gr_ctx_init_mpn_mod` with an unsupported modulus), so that
it can still be cleared.
:::

## Properties of a ring {#properties}

These predicates read flags of the context or compute a property once;
they are cheap and generic algorithms use them to choose an algorithm,
for example fraction-free elimination over an integral domain that is
not a field. The answer `T_UNKNOWN`, the default for most properties,
means that the implementation does not know.

::: function
truth_t gr_ctx_is_ring(gr_ctx_t ctx)
truth_t gr_ctx_is_commutative_ring(gr_ctx_t ctx)
truth_t gr_ctx_is_integral_domain(gr_ctx_t ctx)
truth_t gr_ctx_is_field(gr_ctx_t ctx)
truth_t gr_ctx_is_unique_factorization_domain(gr_ctx_t ctx)
truth_t gr_ctx_is_finite(gr_ctx_t ctx)
truth_t gr_ctx_is_finite_characteristic(gr_ctx_t ctx)
truth_t gr_ctx_is_algebraically_closed(gr_ctx_t ctx)
truth_t gr_ctx_is_ordered_ring(gr_ctx_t ctx)
truth_t gr_ctx_is_zero_ring(gr_ctx_t ctx)
truth_t gr_ctx_is_rational_vector_space(gr_ctx_t ctx)
truth_t gr_ctx_is_real_vector_space(gr_ctx_t ctx)
truth_t gr_ctx_is_complex_vector_space(gr_ctx_t ctx)
truth_t gr_ctx_is_multiplicative_group(gr_ctx_t ctx)

Return whether the structure has the mathematical property. Rings with
inexact elements answer for the ring they represent, so real balls form
a field. Floating-point contexts are not rings and return `T_FALSE` for
the ring properties.
:::

::: function
truth_t gr_ctx_is_approx_commutative_ring(gr_ctx_t ctx)

Returns whether the structure $\tilde R$ implemented by `ctx` is an
approximate version of a commutative ring $R$: the axioms need only hold
up to a small perturbation, $|\tilde a - \tilde b| < \varepsilon$ in
place of $a = b$. Polynomials with floating-point coefficients form an
approximate commutative ring; $2 \times 2$ matrices with floating-point
entries do not. Every commutative ring is an approximate commutative
ring.
:::

::: function
truth_t gr_ctx_is_exact(gr_ctx_t ctx)
truth_t gr_ctx_is_canonical(gr_ctx_t ctx)

Return whether every element is represented exactly (no balls, no
truncated series, no rounding), respectively whether every element has
a unique representation, so that equality can be decided by comparing
representations.
:::

::: function
truth_t gr_ctx_is_threadsafe(gr_ctx_t ctx)

Returns whether different threads may operate simultaneously on
different elements of the ring with the same context; `T_TRUE` by
default. Contexts that cache data without locking (for example
`gr_ctx_init_real_ca`) return `T_FALSE`, and parallel generic algorithms
then run serially.
:::

::: function
truth_t gr_ctx_has_real_prec(gr_ctx_t ctx)
int gr_ctx_set_real_prec(gr_ctx_t ctx, slong prec)
int gr_ctx_get_real_prec(slong * prec, gr_ctx_t ctx)

Return whether `ctx`, or a base ring of it, approximates real or complex
numbers with a working precision (ball arithmetic as well as
floating-point), and set or get that precision in bits. The setters
return `GR_UNABLE` for rings without a precision.
:::

::: function
int gr_ctx_set_is_field(gr_ctx_t ctx, truth_t is_field)

Declares whether the ring is a field. For $\Z/n\Z$ this declares whether
$n$ is prime: `gr_ctx_init_nmod` tests primality when asked, while the
contexts with large moduli (`gr_ctx_init_fmpz_mod`,
`gr_ctx_init_mpn_mod`) report `T_UNKNOWN` until the user declares it,
which enables algorithms for fields. For `gr_ctx_init_gr_poly_quotient`
it declares that the modulus is irreducible. The declaration is not
checked. Returns `GR_UNABLE` for rings that do not support it.
:::

::: function
truth_t gr_ctx_is_pretend_field(gr_ctx_t ctx)
int gr_ctx_set_is_pretend_field(gr_ctx_t ctx, truth_t is_pretend_field)
int gr_ctx_recover_zero_divisor(gr_ptr res, gr_ctx_t ctx)

Let a ring that is not known to be a field compute as if it were one, in
the style of dynamic evaluation [@DDDD1985]. While the pretence is on,
an inversion or division that meets a nonzero non-invertible element
fails with `GR_UNABLE` (not `GR_DOMAIN`, which keeps its meaning of
division by zero) and records that zero divisor in the context;
`gr_ctx_recover_zero_divisor` sets `res` to the first recorded zero
divisor, or returns `GR_UNABLE` if there is none.

::: details
For $\Z/n\Z$ the recorded zero divisor $a$ gives the nontrivial factor
$\gcd(a, n)$ of $n$, so a typical use is to run a field algorithm modulo
an $n$ that is probably prime and to split $n$ when it fails. For
`gr_ctx_init_gr_poly_quotient` the recorded element is a nontrivial
monic factor of the modulus; see [`gr_domains`](gr_domains#poly-quotient).

`gr_ctx_is_pretend_field` returns `T_TRUE` if the ring is a field or
pretends to be one; `gr_ctx_is_field` is not affected by the pretence.
Algorithms that only need inversions to fail cleanly (the Euclidean and
half-GCD algorithms for polynomials, LU decomposition, the default
`gr_gcd` and `gr_canonical_associate`) check `gr_ctx_is_pretend_field`;
algorithms whose correctness depends on the ring being a field (root
finding and factorisation, for example) check `gr_ctx_is_field`.
Setting the pretence to `T_FALSE` forgets the recorded zero divisor.

Supported by `gr_ctx_init_nmod`, `gr_ctx_init_fmpz_mod`,
`gr_ctx_init_mpn_mod` and `gr_ctx_init_gr_poly_quotient`; the recording
is thread-safe (see `_gr_ctx_zero_divisor_lock`). The other rings return
`GR_UNABLE` from the setter and from the recovery function.
:::
:::

## Memory management {#memory}

An element must be initialised before it is used, which sets it to zero,
and cleared afterwards. The guide [](gr/getting-started#elements)
compares the ways of allocating elements.

::: function
void gr_init(gr_ptr res, gr_ctx_t ctx)
void gr_clear(gr_ptr res, gr_ctx_t ctx)

Initialise `res`, which must point to `ctx->sizeof_elem` bytes of memory,
to the zero element of the ring; respectively clear it, freeing any
memory the element owns but not the memory it occupies.
:::

::: function
void gr_swap(gr_ptr x, gr_ptr y, gr_ctx_t ctx)

Swaps `x` and `y` efficiently, without copying owned data.
:::

::: function
void gr_set_shallow(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Copies the bytes of `x` to `res` (by default with `memcpy`), so that
both share any memory owned by `x`. `res` must not be initialised or
cleared independently afterwards; this is for building temporary views.
:::

::: function
gr_ptr gr_heap_init(gr_ctx_t ctx)
void gr_heap_clear(gr_ptr x, gr_ctx_t ctx)

Allocate one element on the heap with `flint_malloc` and initialise it
to zero; clear and free an element allocated this way.
:::

::: function
gr_ptr gr_heap_init_vec(slong len, gr_ctx_t ctx)
void gr_heap_clear_vec(gr_ptr x, slong len, gr_ctx_t ctx)

Allocate an array of `len` contiguous elements on the heap, initialised
to zero; clear the `len` elements and free the array.
:::

::: macro
GR_TMP_INIT(x1, ctx)
GR_TMP_INIT2(x1, x2, ctx)
GR_TMP_INIT3(x1, x2, x3, ctx)
GR_TMP_INIT4(x1, x2, x3, x4, ctx)
GR_TMP_INIT5(x1, x2, x3, x4, x5, ctx)
GR_TMP_CLEAR(x1, ctx)
GR_TMP_CLEAR2(x1, x2, ctx)
GR_TMP_CLEAR3(x1, x2, x3, ctx)
GR_TMP_CLEAR4(x1, x2, x3, x4, ctx)
GR_TMP_CLEAR5(x1, x2, x3, x4, x5, ctx)

Allocate one to five temporary elements contiguously on the stack with
`alloca`, initialised to zero, and assign their addresses to the
`gr_ptr` variables `x1`, ...; the matching `GR_TMP_CLEAR` macro clears
them. An `INIT` and a `CLEAR` macro must match exactly.

::: details
The memory is reclaimed only when the enclosing function returns, so
these macros must not be used inside a loop: allocate before the loop,
move the loop body to a separate function, or use the heap. Elements are
expected to be small (a few kilobytes at most), which holds for all
built-in rings.
:::
:::

::: macro
GR_TMP_INIT_VEC(vec, len, ctx)
GR_TMP_CLEAR_VEC(vec, len, ctx)

Allocate a temporary array of `len` elements initialised to zero and
assign it to `vec`; clear and free it. The array is placed on the stack
if it takes at most `GR_TMP_VEC_ALLOC_MAX_STACK` (1024) bytes and on the
heap otherwise. The same rules about loops apply.
:::

::: macro
GR_TMP_VEC_ALLOC_MAX_STACK
GR_TMP_ALLOC(size)
GR_TMP_FREE(ptr, size)
GR_TMP_ALLOC_SMALL(size)

The allocation helpers behind the temporary macros: `GR_TMP_ALLOC`
allocates on the stack up to `GR_TMP_VEC_ALLOC_MAX_STACK` bytes and with
`flint_malloc` above, `GR_TMP_FREE` frees what `GR_TMP_ALLOC` allocated
on the heap, and `GR_TMP_ALLOC_SMALL` always uses `alloca`.
:::

## Random elements {#random}

These functions are intended for testing. The distributions are chosen by
each ring and are deliberately non-uniform, favouring special values
such as $0$, $\pm 1$, values close to the boundaries of a
representation, and elements with structure (sparse polynomials,
integers with long runs of ones), so that corner cases are found
quickly.

::: function
int gr_randtest(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)

Sets `res` to a random element of the ring.
:::

::: function
int gr_randtest_not_zero(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)

Sets `res` to a random element that is not zero. Returns `GR_DOMAIN` in
the zero ring.
:::

::: function
int gr_randtest_invertible(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)

Sets `res` to a random invertible element. The default tries
`gr_randtest` and `gr_randtest_small` a few times each and falls back
to $\pm 1$.
:::

::: function
int gr_randtest_small(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)

Sets `res` to a small random element, for tests where large elements
would make operations expensive. The default is $a g + b$ for small
integers $a, b$ and the generator $g$ (if any), sometimes divided by a
small integer.
:::

## Input and output {#io}

Output goes through a stream, which is either a C `FILE` or a growing
string; `gr_print` and `gr_get_str` are the usual entry points, and
`flint_printf` accepts `%{gr}` (element and context) and `%{gr*}`
(array, length and context), see [](flint/printf). The output of
`gr_write` can be read back with `gr_set_str` in the rings that
`gr_test_ring` checks (balls are written in the `[m +/- r]` notation).

::: type
gr_stream_struct
gr_stream_t

An output stream: either a `FILE *` (field `fp`) or a heap-allocated
string with fields `s`, `len` and `alloc`.
:::

::: function
void gr_stream_init_file(gr_stream_t out, FILE * fp)
void gr_stream_init_str(gr_stream_t out)

Initialise a stream writing to `fp`, respectively to a new string that
the caller takes over as `out->s` and frees with `flint_free`.
`gr_stream_init_file` is declared only if `stdio.h` was included before
`gr.h`.
:::

::: function
int gr_stream_write(gr_stream_t out, const char * s)
int gr_stream_write_si(gr_stream_t out, slong x)
int gr_stream_write_ui(gr_stream_t out, ulong x)
int gr_stream_write_free(gr_stream_t out, char * s)
int gr_stream_write_fmpz(gr_stream_t out, const fmpz_t x)

Write a string, an integer in decimal, or a string allocated with
`flint_malloc` that is freed afterwards.
:::

::: function
int gr_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx)
int gr_print(gr_srcptr x, gr_ctx_t ctx)
int gr_println(gr_srcptr x, gr_ctx_t ctx)
int gr_get_str(char ** s, gr_srcptr x, gr_ctx_t ctx)

Write `x` to a stream, print it to standard output (with a newline for
`gr_println`), or set `*s` to a newly allocated string that must be
freed with `flint_free`.
:::

::: function
int gr_write_n(gr_stream_t out, gr_srcptr x, slong n, gr_ctx_t ctx)
int gr_get_str_n(char ** s, gr_srcptr x, slong n, gr_ctx_t ctx)

As `gr_write` and `gr_get_str`, rounding real and complex numbers (also
inside polynomials, matrices, ...) to $n$ significant digits. Exact rings
ignore $n$.
:::

::: function
int gr_set_str(gr_ptr res, const char * x, gr_ctx_t ctx)

Sets `res` to the value of the expression `x`, parsed by
`gr_generic_set_str` unless the ring overrides it. Expressions may use
integer and decimal literals (`2.5e-3`, converted with
`gr_set_fmpz_10exp_fmpz`), `+ - * / ^`, parentheses, the generators of
the ring and its base rings by name, the constants `pi`, `i`, `inf` and
`nan`, functions such as `sqrt`, `exp`, `log` and `gamma`, ball
notation `m +/- r`, and `O(x^n)` in series rings; see
`gr_generic_set_str_expr`. Returns `GR_UNABLE` if the string cannot be
parsed and propagates the status of the arithmetic (for example
`GR_DOMAIN` for `"1/2"` in $\Z$).
:::

::: function
int gr_get_fexpr(fexpr_t res, gr_srcptr x, gr_ctx_t ctx)
int gr_get_fexpr_serialize(fexpr_t res, gr_srcptr x, gr_ctx_t ctx)
int gr_set_fexpr(gr_ptr res, fexpr_vec_t inputs, gr_vec_t outputs, const fexpr_t x, gr_ctx_t ctx)

Convert between elements and symbolic expressions (`fexpr`). The
*serialize* version may produce an expression of the internal
representation that is not meant to be read by humans, but that
`gr_set_fexpr` converts back exactly. `gr_set_fexpr` evaluates the
expression in the ring; `inputs` and `outputs` are a dictionary of
symbols and their values (initially empty, also used as scratch space).
These functions are declared only if `fexpr.h` was included before
`gr.h`.
:::

::: macro
GR_FEXPR_SERIALIZE
GR_FEXPR_COMPACT

Flag values for conversions to symbolic expressions.
:::

## Assignment and conversions {#conversions}

::: function
int gr_set(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Sets `res` to a copy of `x`.
:::

::: function
int gr_set_ui(gr_ptr res, ulong x, gr_ctx_t ctx)
int gr_set_si(gr_ptr res, slong x, gr_ctx_t ctx)
int gr_set_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx)
int gr_set_fmpq(gr_ptr res, const fmpq_t x, gr_ctx_t ctx)
int gr_set_d(gr_ptr res, double x, gr_ctx_t ctx)

Set `res` to the image of the integer, rational or double `x` in the
ring.
:::

::: function
int gr_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)

Sets `res` to the element `x` of another ring `x_ctx` converted to
`ctx`.
:::

::: function
int gr_get_si(slong * res, gr_srcptr x, gr_ctx_t ctx)
int gr_get_ui(ulong * res, gr_srcptr x, gr_ctx_t ctx)
int gr_get_fmpz(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
int gr_get_fmpq(fmpq_t res, gr_srcptr x, gr_ctx_t ctx)
int gr_get_d(double * res, gr_srcptr x, gr_ctx_t ctx)

Convert `x` to a C or FLINT type: an integer, a rational number, or a
double (rounded).
:::

::: function
int gr_set_fmpz_2exp_fmpz(gr_ptr res, const fmpz_t a, const fmpz_t b, gr_ctx_t ctx)
int gr_get_fmpz_2exp_fmpz(fmpz_t res1, fmpz_t res2, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the dyadic number $a \cdot 2^b$, respectively write a dyadic
`x` as $m \cdot 2^e$ with $m$ odd (or $m = e = 0$).
:::

::: function
int gr_set_fmpz_10exp_fmpz(gr_ptr res, const fmpz_t a, const fmpz_t b, gr_ctx_t ctx)

Sets `res` to the decimal number $a \cdot 10^b$.
:::

::: function
int gr_get_d_2exp_si(double * res, slong * exp, gr_srcptr x, gr_ctx_t ctx)

Writes a real `x` approximately as $m 2^e$ with $0.5 \le |m| < 1$ (or
$m = e = 0$), which works when `x` is outside the range of a `double`.
:::

## Special values and generators {#special-values}

::: function
int gr_zero(gr_ptr res, gr_ctx_t ctx)
int gr_one(gr_ptr res, gr_ctx_t ctx)
int gr_neg_one(gr_ptr res, gr_ctx_t ctx)

Set `res` to $0$, $1$ or $-1$. `gr_one` fails with `GR_DOMAIN` in
structures without a unit element (for example the domain of matrices
of all shapes created by `gr_ctx_init_matrix_domain`).
:::

::: function
int gr_gen(gr_ptr res, gr_ctx_t ctx)

Sets `res` to the generator of the ring: $x$ in $R[x]$ or $R[[x]]$, the
generator of the defining polynomial in a finite field or number field,
the image of $x$ in a quotient ring $R[x]/(m)$. Returns `GR_UNABLE` for
rings without a distinguished generator.
:::

::: function
int gr_gens(gr_vec_t res, gr_ctx_t ctx)
int gr_gens_recursive(gr_vec_t res, gr_ctx_t ctx)

Set `res` to the vector of generators of the ring (empty by default),
for example $x_1, \ldots, x_n$ in a multivariate polynomial ring. The
*recursive* version also includes the generators of the base rings,
recursively, embedded in `ctx`: in $(\Z[x])[y]$ it gives $x, y$.
`gr_set_str` recognises the names of all these generators.
:::

## Equality {#predicates}

::: function
truth_t gr_is_zero(gr_srcptr x, gr_ctx_t ctx)
truth_t gr_is_one(gr_srcptr x, gr_ctx_t ctx)
truth_t gr_is_neg_one(gr_srcptr x, gr_ctx_t ctx)

Return whether `x` is equal to $0$, $1$ or $-1$.
:::

::: function
truth_t gr_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
truth_t gr_not_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

Return whether `x` and `y` are equal, respectively different, as
mathematical elements of the ring; for inexact representations the
answer is `T_UNKNOWN` unless it can be decided from the enclosures.
:::

## Arithmetic {#arithmetic}

Every ring implements `neg`, `add`, `sub` and `mul`. The variants with
an operand of type `ulong`, `slong`, `fmpz_t` or `fmpq_t` convert that
operand to the ring unless the ring has a faster method; a `fmpq` operand
fails with `GR_DOMAIN` if its denominator is not invertible. The
`_other` variants take one operand from another ring `x_ctx` or `y_ctx`
and convert it with `gr_set_other` (the order of the operands matters
for `gr_other_mul` in noncommutative rings). All functions allow `res`
to alias any input in the same ring.

::: function
int gr_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Sets `res` to $-x$.
:::

::: function
int gr_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_add_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_add_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_add_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int gr_add_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
int gr_add_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
int gr_other_add(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

Set `res` to $x + y$.
:::

::: function
int gr_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_sub_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_sub_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_sub_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int gr_sub_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
int gr_sub_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
int gr_other_sub(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

Set `res` to $x - y$.
:::

::: function
int gr_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_mul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_mul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_mul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int gr_mul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
int gr_mul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
int gr_other_mul(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

Set `res` to $x y$.
:::

::: function
int gr_addmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_addmul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_addmul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_addmul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int gr_addmul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
int gr_addmul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
int gr_submul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_submul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_submul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_submul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int gr_submul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
int gr_submul_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)

Set `res` to $\mathit{res} + x y$, respectively $\mathit{res} - x y$.
The default forms the product in a temporary; rings may fuse the
operation, avoiding the temporary or an intermediate rounding.
For sums of many products use the dot products of `gr_vec`
(`_gr_vec_dot`), which rings optimise further.
:::

::: function
int gr_mul_two(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_sqr(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to $2x$, respectively $x^2$; by default $x + x$ and $x \cdot x$.
:::

::: function
int gr_mul_2exp_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_mul_2exp_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)

Set `res` to $x \cdot 2^y$. For negative $y$ this divides by $2^{-y}$
(failing with `GR_DOMAIN` if the division is not possible), so it also
works in rings where $2^y$ is not representable as an element but the
quotient is.
:::

## Division {#division}

The function page [](gr/division) explains the differences between these
functions, with examples.

::: function
truth_t gr_is_invertible(gr_srcptr x, gr_ctx_t ctx)
int gr_inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Return whether `x` is a unit, respectively set `res` to $x^{-1}$,
failing with `GR_DOMAIN` if `x` is not invertible.
:::

::: function
int gr_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_div_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_div_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_div_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int gr_div_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
int gr_div_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
int gr_other_div(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

Set `res` to the quotient $q = x / y$, the unique $q$ with $q y = x$ in
the ring. Fail with `GR_DOMAIN` if there is no such $q$ or it is not
unique (for example $0/0$, $1/2$ in $\Z$, $3/2$ in $\Z/6\Z$).
:::

::: function
int gr_div_nonunique(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
truth_t gr_divides(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

`gr_div_nonunique` sets `res` to some $q$ with $q y = x$, failing with
`GR_DOMAIN` only if no such $q$ exists; `gr_divides` returns whether
$x \mid y$, that is, whether $y = q x$ for some $q$.
:::

::: function
int gr_divexact(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_divexact_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_divexact_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_divexact_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int gr_divexact_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
int gr_divexact_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
int gr_other_divexact(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

Set `res` to $x / y$ assuming that the division is exact. The ring may
skip the verification: if the division is not exact, the result may be
meaningless even when `GR_SUCCESS` is returned.
:::

::: function
int gr_euclidean_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_euclidean_rem(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_euclidean_divrem(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

Euclidean division $x = q y + r$ with remainder, setting `res` (or
`res1`) to $q$ and `res` (or `res2`) to $r$; division by zero fails with
`GR_DOMAIN`. The choice of quotient is defined by the ring: the floor
quotient in $\Z$ (as `fmpz_fdiv_qr`), a nearest quotient in $\Z[i]$,
polynomial division in $F[x]$. Rings without a Euclidean structure
return `GR_UNABLE` (no default implementation).
:::

## Powers {#powering}

::: function
int gr_pow(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_pow_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx)
int gr_pow_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx)
int gr_pow_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx)
int gr_pow_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx)
int gr_pow_other(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)
int gr_other_pow(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_srcptr y, gr_ctx_t ctx)

Set `res` to $x^y$. Fail with `GR_DOMAIN` if the power has no value in
the ring (for example $0^{-1}$, $2^{-1}$ in $\Z$, $(-1)^{1/2}$ in
$\mathbb{R}$).

::: details
For integer exponents the default implementations use binary
exponentiation, inverting `x` first when the exponent is negative
(`gr_generic_pow_ui`, `gr_generic_pow_fmpz`). For $y = p/2$ the default
of `gr_pow_fmpq` uses `gr_sqrt` or `gr_rsqrt` and then an integer power;
other rational exponents are converted to the ring and passed to
`gr_pow`. In subrings of $\mathbb{C}$ a non-integer power is the
principal branch $x^y = \exp(y \log x)$ for $x \ne 0$. In other rings
the meaning of a non-integer power is defined by the ring; most return
`GR_UNABLE`.
:::
:::

::: function
int gr_derivative_gen(gr_ptr res, gr_srcptr x, slong var, gr_ctx_t ctx)

Sets `res` to the derivative of `x` with respect to the generator of
index `var` of the ring (for polynomial and series rings).
:::

## Square roots {#roots}

In subrings of $\mathbb{C}$ the principal square root is computed; in
other rings (for example finite fields) the ring chooses a root. The
default implementations handle $0$ and $1$ and return `GR_UNABLE`
otherwise.

::: function
truth_t gr_is_square(gr_srcptr x, gr_ctx_t ctx)

Returns whether `x` is a square in the ring.
:::

::: function
int gr_sqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_rsqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to a square root $\sqrt{x}$, respectively $1/\sqrt{x}$. Fail
with `GR_DOMAIN` if `x` is not a square (and for `gr_rsqrt` if $x = 0$)
and with `GR_UNABLE` if this cannot be decided, for example for a real
ball containing zero.
:::

## Divisibility and factorisation {#gcd}

::: function
int gr_gcd(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_lcm(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

Set `res` to a greatest common divisor, respectively least common
multiple, of `x` and `y`. These are unique only up to units; the ring
chooses a representative, typically the canonical associate (nonnegative
in $\Z$, monic in $F[x]$, 1 or 0 in a field).
:::

::: function
int gr_canonical_associate(gr_ptr res1, gr_ptr res2, gr_srcptr x, gr_ctx_t ctx)

Sets `res1` to a canonical representative $u x$ of the class of `x`
under multiplication by units and `res2` to the unit $u$. For $x = 0$
sets `res1` to 0 and `res2` to 1. Fraction fields use it to normalise
denominators.
:::

::: function
int gr_factor(gr_ptr c, gr_vec_t factors, fmpz_vec_t exponents, gr_srcptr x, int flags, gr_ctx_t ctx)

Computes a factorisation $x = c f_1^{e_1} \cdots f_n^{e_n}$ into
distinct irreducible (or prime) factors $f_k$ with positive exponents.
`factors` and `exponents` are resized to length $n$.

::: details
The unit or content $c$ is an element of the same ring as the input: the
sign $-1$, $0$ or $1$ in $\Z$, the sign times the content (as a constant
polynomial) in $\Z[x]$, the leading coefficient in $F[x]$. The factors
are not sorted in any particular order. `flags` selects options of the
ring's algorithm; 0 gives the default. Implemented for $\Z$, $\mathbb{Q}$,
multivariate polynomials and rational functions over $\Z$ and $\mathbb{Q}$,
and univariate polynomials (`gr_ctx_init_gr_poly`, `gr_ctx_init_fmpz_poly`)
over $\Z$, $\Z/p\Z$, finite fields, algebraic numbers (`qqbar`, `ca`)
and complex balls and floating-point numbers, the last ones by
root-finding. Other rings return `GR_UNABLE`.
:::
:::

::: function
int gr_numerator(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_denominator(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to a numerator $p$, respectively a denominator $q$, with
$x = p / q$, both in the same ring as `x`. Fraction fields return the
reduced numerator and the canonical denominator; the default is $p = x$
and $q = 1$.
:::

## Integer, real and complex parts {#parts}

These functions are defined for subrings of $\mathbb{C}$ (and some for
ordered rings); other rings return `GR_DOMAIN` or `GR_UNABLE`.

::: function
int gr_floor(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_ceil(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_trunc(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_nint(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the integer obtained by rounding `x` towards $-\infty$,
$+\infty$, zero, or to the nearest integer with ties to even. For
complex balls (`acb`) the real part is rounded and the result is real.
:::

::: function
int gr_abs(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Sets `res` to the absolute value of `x`, in an ordered ring or a subring
of $\mathbb{C}$.
:::

::: function
int gr_i(gr_ptr res, gr_ctx_t ctx)

Sets `res` to the imaginary unit $i$; `GR_DOMAIN` in real rings.
:::

::: function
int gr_conj(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_re(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_im(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_sgn(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_csgn(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
int gr_arg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the complex conjugate, the real part, the imaginary part,
the sign $x / |x|$ (0 for $x = 0$), the real-valued sign used for branch
cuts (`csgn`: the sign of $\operatorname{Re}(x)$, or of
$\operatorname{Im}(x)$ when the real part is zero), and the argument
in $(-\pi, \pi]$.
:::

## Infinities and special values {#infinities}

::: function
int gr_pos_inf(gr_ptr res, gr_ctx_t ctx)
int gr_neg_inf(gr_ptr res, gr_ctx_t ctx)
int gr_uinf(gr_ptr res, gr_ctx_t ctx)
int gr_undefined(gr_ptr res, gr_ctx_t ctx)
int gr_unknown(gr_ptr res, gr_ctx_t ctx)

Set `res` to $+\infty$, $-\infty$, the unsigned infinity
$\tilde \infty$, the value *Undefined* or the value *Unknown*. Only
extended structures such as `gr_ctx_init_complex_extended_ca` and the
floating-point contexts have such elements; other rings return
`GR_DOMAIN` or `GR_UNABLE`.
:::

## Ordering {#ordering}

::: function
int gr_cmp(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_cmp_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)

Set `*res` to $-1$, $0$ or $1$ according to whether $x < y$, $x = y$ or
$x > y$. Fail with `GR_DOMAIN` if the ring is not ordered and with
`GR_UNABLE` if the comparison cannot be decided (overlapping balls).
:::

::: function
int gr_cmpabs(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_cmpabs_other(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t y_ctx, gr_ctx_t ctx)

As `gr_cmp`, comparing the absolute values $|x|$ and $|y|$; this is
defined also for complex numbers.
:::

::: function
truth_t gr_le(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
truth_t gr_lt(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
truth_t gr_ge(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
truth_t gr_gt(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
truth_t gr_abs_le(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
truth_t gr_abs_lt(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
truth_t gr_abs_ge(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
truth_t gr_abs_gt(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

Return the truth of $x \le y$, $x < y$, $x \ge y$, $x > y$ and of the
same relations between $|x|$ and $|y|$. They call `gr_cmp` or
`gr_cmpabs` and return `T_UNKNOWN` if it fails.
:::

::: function
int gr_min(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int gr_max(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

Set `res` to the minimum or maximum of `x` and `y`. Ball arithmetic
returns a ball enclosing the minimum or maximum of all points in the
input balls.
:::

## Enclosures and error terms {#enclosures}

::: function
int gr_set_interval_mid_rad(gr_ptr res, gr_srcptr m, gr_srcptr r, gr_ctx_t ctx)

In ball representations of real numbers, sets `res` to the ball
$m \pm r$.

::: details
In vector spaces over the real numbers represented with balls, the radii
are applied to the components separately: in the complex numbers,
$a + bi \pm (0.1 + 0.2i)$ means $(a \pm 0.1) + (b \pm 0.2) i$, and in a
polynomial ring each coefficient gets its own radius. `gr_set_str`
parses `m +/- r` with this function. Exact rings return `GR_UNABLE`.
:::
:::

::: function
int gr_set_interval_mid_inf(gr_ptr res, gr_srcptr m, gr_ctx_t ctx)

Sets `res` to $m \pm \infty$, the whole real line for a real ball. In
vector spaces over the reals the infinite radius is a real constant, so
only the real part of a complex number, or the constant coefficient of
a polynomial, gets an infinite radius, as with
`gr_set_interval_mid_rad`.

::: details
`gr_set_str` parses `m +/- inf` (and `+/- inf` with $m = 0$) by passing
$m$ to this function without evaluating `inf`, so that such balls can
be read back in rings where infinity is not an element. The default
uses `gr_pos_inf` and `gr_set_interval_mid_rad` and returns `GR_UNABLE`
when the ring has no infinity.
:::
:::

::: function
int gr_get_interval_mid_rad(gr_ptr m, gr_ptr r, gr_srcptr x, gr_ctx_t ctx)

Sets `m` and `r` to the midpoint and radius of `x` as exact elements of
the ring, so that `gr_set_interval_mid_rad` with `m` and `r`
reproduces `x`; radii of vector spaces are split as in the setter. The
default, used by exact rings, sets $m = x$ and $r = 0$.
:::

::: function
int gr_big_o_base_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t n, gr_ctx_t ctx)

Sets `res` to the error term $O(x^n)$, where `x` is the generator of a
power series ring or the base $p$ of a $p$-adic ring. Returns
`GR_UNABLE` in rings that cannot represent such a term (for example an
exact polynomial ring).

::: details
If `x` is not the generator of `ctx`, the ring searches its base rings:
in $(R[[x]])[[y]]$, $O(x^2)$ means $O(x^2) y^0$, so that nested error
terms such as $(x + O(x^5)) + (3 + x + O(x^4)) y^2 + O(y^5)$ can be
built. `gr_set_str` parses `O(x^n)` by calling this function with `x`
and `n`, without evaluating the power $x^n$.
:::
:::

## Finite fields {#finite-fields}

These methods are implemented by the finite field contexts
(`gr_ctx_init_fq`, `gr_ctx_init_fq_nmod`, `gr_ctx_init_fq_zech`). The
contexts for $\Z/n\Z$ (`nmod`, `fmpz_mod`, `mpn_mod`) implement
`gr_ctx_fq_prime`, `gr_ctx_fq_degree` and `gr_ctx_fq_order` when $n$ is
known to be prime (see `gr_ctx_set_is_field`), with degree 1 and order
$n$. Other rings return `GR_UNABLE`.

::: function
int gr_ctx_fq_prime(fmpz_t p, gr_ctx_t ctx)
int gr_ctx_fq_degree(slong * deg, gr_ctx_t ctx)
int gr_ctx_fq_order(fmpz_t q, gr_ctx_t ctx)

Set the characteristic $p$, the degree $d$ over the prime field and the
order $q = p^d$ of the finite field $\mathbb{F}_q$.
:::

::: function
int gr_fq_frobenius(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx)

Sets `res` to $x^{p^e}$, the image of `x` under the $e$-th power of the
Frobenius automorphism.
:::

::: function
int gr_fq_multiplicative_order(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
truth_t gr_fq_is_primitive(gr_srcptr x, gr_ctx_t ctx)

Set `res` to the multiplicative order of `x` (`GR_DOMAIN` for
$x = 0$), respectively return whether `x` generates the multiplicative
group $\mathbb{F}_q^\times$.
:::

::: function
int gr_fq_norm(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
int gr_fq_trace(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)

Set `res` to the norm or the trace of `x` over the prime field, as an
integer in $[0, p)$.
:::

::: function
int gr_fq_pth_root(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)

Sets `res` to the unique $p$-th root of `x`.
:::

## Transformed representations {#transformed}

A *transformed representation* holds ring elements in a form in which
multiplication is a cheap pointwise operation, at the price of
conversions in and out: for example polynomials as their values under a
number-theoretic transform. Algorithms performing many multiplications
between few conversions (matrix products, half-GCD, middle products) run
substantially faster in such a representation. The functions below
require FLINT built with `fft_small`; otherwise the constructors return
`GR_UNABLE` and callers fall back to their ordinary code paths. The
design of the transform space is described in
[](fft_small/transform-space).

::: type
gr_transformed_poly_workload_struct
gr_transformed_poly_workload_t

An advisory estimate of the intended workload: the number of operands
that will be converted in (`num_inputs`), of pointwise multiplications
(`num_muls`), of results converted out (`num_outputs`), of elements
alive at the same time (`num_live`, 0 to derive
`num_inputs + num_outputs + 2`), a memory bound in bytes (`mem_limit`, 0
for the implementation default), and a flag `force` that skips the
profitability model and the storage budget, declining only on
implementation bounds (for tests, where small unprofitable sizes catch
bugs most easily). Implementations use it to decide whether switching
representation pays off and to size their tables.
:::

::: function
int gr_ctx_init_gr_poly_transformed_repr(gr_ctx_t ctx, gr_ctx_t base, slong len_bound, slong terms_bound, const gr_transformed_poly_workload_struct * workload)

Constructs a ring of transformed polynomials over `base` admitting
lengths up to `len_bound` and up to `terms_bound` accumulated elementary
products per element, for the estimated workload (which may be `NULL`,
though implementations then usually judge the switch unprofitable).
Returns `GR_UNABLE` if the base ring provides no transformed
representation (the method `GR_METHOD_CTX_INIT_TRANSFORMED_POLY_REPR`,
implemented for `nmod` and `mpn_mod` base rings) or judges the switch
unprofitable.

::: details
Conversions are performed by the `SET_GR_POLY`, `GET_GR_POLY` and
`GET_GR_POLY_WINDOW` methods of the constructed context; the window
variant writes the coefficients $[z_l, z_h)$ of the represented
polynomial, with zeros beyond its length. `_gr_get_gr_poly_destructive`
converts out through the `GET_GR_POLY_DESTRUCTIVE` method where the
representation provides it: the inverse transforms then run on the
element's own storage, skipping a transform-sized copy, and the element
may only be cleared or fully overwritten afterwards; representations
without the method fall back to the copying conversion.
:::
:::

::: function
int gr_ctx_init_gr_poly_transformed_cyclic_repr(gr_ctx_t ctx, gr_ctx_t base, slong * len, slong terms_bound, const gr_transformed_poly_workload_struct * workload)

Constructs a ring of transformed polynomials modulo $x^L - 1$ over
`base`, in which products wrap around. The implementation may round the
requested cyclic length `*len` up (for example to a power of two) and
writes the actual $L$ back. Conversions in require length at most $L$;
conversions out give the coefficients of the reduced representative.
This is the `CTX_INIT_TRANSFORMED_POLY_CYCLIC_REPR` method of the base
ring (`GR_UNABLE` by default); it is used for middle products and for
products whose high part is known, such as the remainder step of
division by a precomputed modulus.
:::

::: function
int _gr_nmod_ctx_init_transformed_poly_repr(gr_ctx_t ctx, gr_ctx_t base, slong len_bound, slong terms_bound, const gr_transformed_poly_workload_struct * workload)
int _gr_nmod_ctx_init_transformed_poly_cyclic_repr(gr_ctx_t ctx, gr_ctx_t base, slong * len, slong terms_bound, const gr_transformed_poly_workload_struct * workload)
int _gr_nmod_tpoly_get_gr_poly_window_destructive(nn_ptr cc, gr_ptr x, slong zl, slong zh, gr_ctx_t ctx)

The implementations of the two constructors for `nmod` base rings,
installed in their method table, and the destructive windowed conversion
of the resulting ring, which writes the coefficients $[z_l, z_h)$ of `x`
to `cc` and consumes `x`.
:::

::: function
int gr_ctx_init_transformed_mpn(gr_ctx_t ctx, slong bits_bound, slong terms_bound, int is_signed, slong num_live, int alloc_strategy)

Constructs the ring of transformed big integers: bilinear expressions
over $\Z$ whose value is below $2^{\mathit{bits\_bound}}$ in absolute
value, with at most `terms_bound` accumulated elementary products and
multiplicative depth two. Elements carry a sign bit when `is_signed` is
nonzero; mixed-sign accumulations switch the pointwise additions and
subtractions, and the sign of a mixed result is resolved at conversion
out.

::: details
`bits_bound` bounds the magnitude of any single represented product; the
representation provisions the accumulation capacity from `terms_bound`
and the sign headroom from `is_signed` itself, so callers must not add
either to the bound. `num_live` declares how many elements are expected
alive simultaneously; the constructor declines (`GR_UNABLE`) when their
combined transform storage would exceed
`flint_fft_small_max_transformed_ring_size`.

`alloc_strategy` selects how element storage is obtained.
`GR_TRANSFORMED_MPN_ALLOC_MALLOC` allocates each element separately.
`GR_TRANSFORMED_MPN_ALLOC_FIT_BUFFER` reserves a stable region of the
thread's `fft_small` scratch buffer holding `num_live` element slabs,
from which `gr_init` draws and to which `gr_clear` returns;
initialisations beyond `num_live` simultaneously live elements fall back
to separate allocation. The reservation makes repeated context creation
cheap, since the region persists in the thread's scratch buffer across
contexts. It is designed for a context that is the only user of
`fft_small` on its thread while it lives (construct, compute, destroy),
which is how every library caller uses it. This is an assumption rather
than a requirement: while a reservation is live the thread's scratch
buffer can neither grow nor move, and an interleaved `fft_small`
operation, such as an integer multiplication large enough to reach it,
is served from a secondary buffer instead (see `mpn_ctx_fit_buffer`). A
context held across unrelated work therefore stays correct. One
reservation exists per thread, so a context constructed while another
holds it silently uses separate allocation. Elements must be initialised
with `gr_init` or `GR_TMP_INIT_VEC` in every strategy.

Conversions out whose destination window is too short to hold the
reconstruction are staged inside the ring, and that staging is
provisioned with the reservation (well under one element's storage,
reserved past the slabs). A caller that will always have an element
dead by the time it converts may add the flag
`GR_TRANSFORMED_MPN_SCRATCH_FROM_SLAB` to `alloc_strategy`, and the ring
then takes a free slab for the staging instead of reserving anything.
The flag declares a footprint, not a correctness contract: if no slab is
free the ring falls back to an allocation cached for the lifetime of the
context. One conversion may be in flight at a time.

A context may be shared by several threads, each working on its own
elements: element storage and staging are managed by an operand pool
(shared with the complex products of `fft_small`) whose bookkeeping is
locked, and every conversion takes independent scratch. A context
holding a reservation pins the creating thread's scratch buffer and must
be cleared on that thread, after all its elements.
:::
:::

::: macro
GR_TRANSFORMED_MPN_ALLOC_MALLOC
GR_TRANSFORMED_MPN_ALLOC_FIT_BUFFER
GR_TRANSFORMED_MPN_ALLOC_STRATEGY_MASK
GR_TRANSFORMED_MPN_SCRATCH_FROM_SLAB

Values of the `alloc_strategy` argument: the two strategies (0 and 1),
the mask selecting the strategy bit, and the flag (2) that takes
conversion staging from a free slab.
:::

::: function
int gr_ctx_init_transformed_mpn_window(gr_ctx_t ctx, slong bits_bound, slong terms_bound, int is_signed, slong num_live, int alloc_strategy, slong lo, slong hi)

As `gr_ctx_init_transformed_mpn`, for products of which only the limbs
$[\mathit{lo}, \mathit{hi})$ are wanted. The plan is the window plan of
`fft_small_plan_init_mpn_window`: when the window allows it, the
transforms are cyclic of a length shorter than the products, and a
product whose slots then alias can be converted only through
`gr_transformed_mpn_get_window` with a window inside
$[\mathit{lo}, \mathit{hi})$; the other conversions of such a product
return `GR_UNABLE`. Elements that do not alias (operands, and products
short enough) convert as usual.
:::

::: function
int gr_transformed_mpn_use_fit_buffer(gr_ctx_t ctx, slong num_live)

Switches a context created with `GR_TRANSFORMED_MPN_ALLOC_MALLOC` to the
reserved-slab strategy for `num_live` simultaneously live elements, as
if it had been created that way. This serves callers whose live count
depends on quantities derived from the constructed context, such as a
blocking decision based on the element size. Returns `GR_UNABLE`,
leaving separate allocation in place, if the reservation cannot be made.
:::

::: function
int gr_transformed_mpn_set(gr_ptr res, nn_srcptr a, slong an, int sign, gr_ctx_t ctx)
int gr_transformed_mpn_get(nn_ptr z, slong zn, slong * zn_out, int * sign, gr_srcptr x, gr_ctx_t ctx)
int gr_transformed_mpn_get_destructive(nn_ptr z, slong zn, slong * zn_out, int * sign, gr_ptr x, gr_ctx_t ctx)
slong gr_transformed_mpn_get_limbs(gr_ctx_t ctx, gr_srcptr x)
slong gr_transformed_mpn_get_limbs_bound(gr_ctx_t ctx)

Conversions by limb arrays with an explicit sign. `gr_transformed_mpn_set`
converts the integer $\pm(a_0, \ldots, a_{an-1})$ in, failing with
`GR_DOMAIN` if it does not fit the context. The `get` functions write
the absolute value to `z`, its length to `*zn_out` and the sign to
`*sign`; the *destructive* form may consume the element, skipping a
copy of the transform. `gr_transformed_mpn_get_limbs` returns the number
of limbs the reconstruction of `x` needs, which exceeds the length of
the value by the length of the top CRT coefficient, and
`gr_transformed_mpn_get_limbs_bound` an upper bound over every element of
the context.

::: details
A window of at least `gr_transformed_mpn_get_limbs` limbs receives the
reconstruction directly, without staging or copying. A shorter window is
also accepted: the ring stages the conversion in its own scratch and
copies the result out, succeeding whenever the value fits in `zn` limbs
and returning `GR_DOMAIN` otherwise. For the destructive form that
refusal comes with the element already consumed, so it cannot be
retried. Limbs of `z` above `*zn_out` are zeroed. Callers converting
into an `fmpz` or into a window of `get_limbs` limbs need not consult
`get_limbs_bound`; a caller that wants a buffer of its own should size it
from that bound rather than from its own model of the representation.
:::
:::

::: function
int gr_transformed_mpn_get_trunc(nn_ptr z, slong zn, slong * zn_out, int * sign, slong lo, gr_srcptr x, gr_ctx_t ctx)
int gr_transformed_mpn_get_trunc_destructive(nn_ptr z, slong zn, slong * zn_out, int * sign, slong lo, gr_ptr x, gr_ctx_t ctx)
slong gr_transformed_mpn_get_limbs_trunc(gr_ctx_t ctx, gr_srcptr x, slong lo)

Truncated conversion out: the limbs of the value starting at limb `lo`,
with an error against the exact value within $(-1.5, +0.5)$ units of
the lowest returned limb, equivalently at most 1 from the floor of the
truncated value. The export includes every CRT slot whose coefficient
span can reach the first returned limb, with carries propagated, so the
error is the truncation itself (one-sided, below one unit) plus the
total of the wholly dropped slots (under half a unit either way).
`gr_transformed_mpn_get_limbs_trunc` returns the number of limbs needed.
:::

::: function
int gr_transformed_mpn_get_low(nn_ptr z, slong zn, gr_srcptr x, gr_ctx_t ctx)
int gr_transformed_mpn_get_low_destructive(nn_ptr z, slong zn, gr_ptr x, gr_ctx_t ctx)

Residue conversion out: `z` receives the value modulo $2^{\bits \cdot
\mathit{zn}}$ as a two's-complement residue of exactly `zn` limbs, with
no sign, whatever the size of the value. Only the CRT slots starting
below limb `zn` are reconstructed, so this costs a fraction of a full
conversion when `zn` is a fraction of the length of the value.
:::

::: function
int gr_transformed_mpn_get_window(nn_ptr z, slong zn, slong lo, gr_srcptr x, gr_ctx_t ctx)
int gr_transformed_mpn_get_window_destructive(nn_ptr z, slong zn, slong lo, gr_ptr x, gr_ctx_t ctx)

Window conversion out, for signed contexts: `z` receives the limbs
$[\mathit{lo}, \mathit{lo} + \mathit{zn})$ of the value as a
two's-complement window, within one unit of $\lfloor x / B^{lo}
\rfloor$ (two when the element carries a small side value, which is
folded in exactly only for $\mathit{lo} = 0$) and exact for
$\mathit{lo} = 0$; the value need not fit. The sign is applied in
transform space, so the window is that of the signed value. This is the
conversion for windowed contexts with cyclic plans. Returns `GR_UNABLE`
for unsigned contexts.
:::

::: function
int gr_transformed_mpn_get_fmpz_destructive(fmpz_t f, slong lo_limbs, gr_ptr x, gr_ctx_t ctx)

Converts `x` into `f`, consuming it. The width the reconstruction needs
is taken from the element and `f` is grown to match, so the conversion
writes straight into the limbs of `f`, with neither staging nor a copy.
A nonzero `lo_limbs` selects the truncated conversion, with the error
described for `gr_transformed_mpn_get_trunc`. This is the conversion to
use whenever the destination is an `fmpz`.
:::

::: function
ulong gr_transformed_mpn_sizeof_data(gr_ctx_t ctx)

Returns the size in bytes of the transform storage of one element. A
context with reserved slabs occupies `num_live` times this, plus the
conversion staging region when it is not taken from a slab.
:::

## Method tables {#method-tables}

The guide [](gr/implementing) explains how to define a ring with these
types and functions.

::: type
gr_funcptr

The type `int (*)(void)` of a method table entry; methods are stored
with this type and cast back to their actual signature by the wrapper
functions.
:::

::: type
gr_method

The enumeration indexing method tables: the method of `gr_add` is at
index `GR_METHOD_ADD`, that of `gr_ctx_is_field` at
`GR_METHOD_CTX_IS_FIELD`, and so on. `GR_METHOD_TAB_SIZE` is the number
of methods. Some methods (for example `GR_METHOD_IS_INTEGER`,
`GR_METHOD_ROOT_OF_UNITY_UI`) have no wrapper function yet.
:::

::: type
gr_static_method_table

An array of `GR_METHOD_TAB_SIZE` entries of type `gr_funcptr`.
:::

::: type
gr_method_tab_input

A pair (`index`, `function`) of a `gr_method` and a `gr_funcptr`; arrays
of such pairs, terminated by an entry with a `NULL` function, describe
the methods a ring implements.
:::

::: function
void gr_method_tab_init(gr_funcptr * methods, gr_method_tab_input * tab)

Initialises the method table `methods`: every entry is first set to
`gr_not_implemented`, then the defaults of [`gr_generic`](gr_generic)
are installed, and finally the entries listed in `tab` overwrite them.
The table is built in a temporary and copied at the end, so that
concurrent initialisation of the same table by several threads is
harmless.
:::

::: function
void gr_method_tab_extend(gr_funcptr * methods, gr_method_tab_input * tab)

Overwrites the entries listed in `tab` in an initialised table, for
rings that share most methods with another table.
:::

::: function
int gr_not_implemented(void)
int gr_not_in_domain(void)

Do nothing and return `GR_UNABLE`, respectively `GR_DOMAIN`. The first
fills the entries of methods without an implementation; the second can
be installed for an operation that never makes sense in the ring, such
as the constant $\pi$ in $\mathbb{Q}$.
:::

::: function
truth_t gr_generic_ctx_predicate(gr_ctx_t ctx)
truth_t gr_generic_ctx_predicate_true(gr_ctx_t ctx)
truth_t gr_generic_ctx_predicate_false(gr_ctx_t ctx)

Context predicates that return `T_UNKNOWN`, `T_TRUE` and `T_FALSE`, for
filling the `GR_METHOD_CTX_IS_*` entries of a method table.
:::

::: type
gr_which_structure

The enumeration of implementation codes stored in
`ctx->which_ring`: `GR_CTX_FMPZ`, `GR_CTX_NMOD`, `GR_CTX_RR_ARB`,
`GR_CTX_GR_POLY` and so on, with `GR_CTX_UNKNOWN_DOMAIN` for rings
defined outside FLINT. Generic code uses it for conversions between
specific rings and for shortcuts; `gr_ctx_cmp_coercion` compares these
codes.
:::

## Testing {#testing}

These functions test a ring implementation with random elements; they
are used by FLINT's test suite and are useful for checking a new ring.
A test fails only on a definite contradiction: an identity such as
$(xy)z = x(yz)$ is reported as wrong only if `gr_equal` returns
`T_FALSE` and all operations succeeded, so inexact rings and partial
implementations pass as long as they never return a wrong result. The
algorithm page [](algorithms/generic/testing) lists the properties
checked.

::: macro
GR_TEST_VERBOSE
GR_TEST_ALWAYS_ABLE

Flags for the `test_flags` arguments: print each test with the counts
of successful, domain and unable outcomes and the timing; and treat
`GR_UNABLE` as a failure, for rings that are expected to compute every
operation.
:::

::: function
void gr_test_ring(gr_ctx_t R, slong iters, int test_flags)

Tests the ring `R`, running each test `iters` times (vector tests
`iters / 10 + 1` times): conversions and printing, the ring axioms
(associativity, commutativity if `gr_ctx_is_commutative_ring` is
`T_TRUE`, distributivity), aliasing of the arithmetic functions,
consistency of the type variants (`_ui`, `_si`, `_fmpz`, `_fmpq`,
`_other`), division, inversion and exact division, powering, square
roots, ordering, GCD, LCM and factorisation when the ring is a unique
factorisation domain, vector operations, dot products and classical
matrix multiplication. Throws an exception on the first failure.
:::

::: function
void gr_test_multiplicative_group(gr_ctx_t R, slong iters, int test_flags)

Tests `R` as a multiplicative group: `gr_one`, associativity and
aliasing of `gr_mul`, commutativity if the group is commutative,
division, inversion and powering.
:::

::: function
void gr_test_floating_point(gr_ctx_t R, slong iters, int test_flags)

Tests a floating-point context, for which the ring axioms hold only
approximately: initialisation, printing and parsing, commutativity of
addition, aliasing of the arithmetic functions, vector operations and
conversions.
:::

::: function
int gr_test_cmp_fun(gr_ctx_t R, gr_method_binary_op_get_int op, gr_ctx_t R_ref, flint_rand_t state, int test_flags)

Compares the comparison function `op` (such as `gr_cmp`) on random
elements of `R` with the same function on the elements converted to the
reference ring `R_ref`.
:::

::: function
int gr_test_approx_unary_op(gr_ctx_t R, gr_method_unary_op op, gr_ctx_t R_ref, gr_srcptr rel_tol, flint_rand_t state, int test_flags)
int gr_test_approx_binary_op(gr_ctx_t R, gr_method_binary_op op, gr_ctx_t R_ref, gr_srcptr rel_tol, flint_rand_t state, int test_flags)
int gr_test_approx_binary_op_type_variants(gr_ctx_t R, const char * opname, gr_method_binary_op gr_op, gr_method_binary_op_ui gr_op_ui, gr_method_binary_op_si gr_op_si, gr_method_binary_op_fmpz gr_op_fmpz, gr_method_binary_op_fmpq gr_op_fmpq, int fused, int small_test_values, gr_srcptr rel_tol, flint_rand_t state, int test_flags)
int gr_test_approx_dot(gr_ctx_t R, gr_ctx_t R_ref, slong maxlen, gr_srcptr rel_tol, flint_rand_t state, int test_flags)

Tests for approximate (floating-point) rings: compute an operation, its
type variants, or a dot product of length at most `maxlen` in `R` and
in a more precise reference ring `R_ref` (for example balls), and check
that the relative error is at most `rel_tol`, an element of `R_ref`.
For the type variants, `fused` selects the `addmul`-type signature and
`small_test_values` restricts the scalar operands to small values.
These functions abort on failure.
:::

::: function
int gr_test_factor(gr_ctx_t R, flint_rand_t state, int test_flags)

Tests `gr_factor` on a small random element, sometimes a product of
two: the lengths of `factors` and `exponents` agree and
$c \prod f_k^{e_k}$ equals the input. `GR_DOMAIN` is accepted only for
the input zero.
:::

::: function
int gr_test_pow_ui_exponent_addition(gr_ctx_t R, flint_rand_t state, int test_flags)
int gr_test_pow_ui_base_scalar_multiplication(gr_ctx_t R, flint_rand_t state, int test_flags)
int gr_test_pow_ui_base_multiplication(gr_ctx_t R, flint_rand_t state, int test_flags)
int gr_test_pow_ui_aliasing(gr_ctx_t R, flint_rand_t state, int test_flags)
int gr_test_pow_fmpz_exponent_addition(gr_ctx_t R, flint_rand_t state, int test_flags)

Single tests of `gr_pow_ui` and `gr_pow_fmpz`: $x^{a+b} = x^a x^b$,
$(c x)^a = c^a x^a$ for an integer $c$, $(xy)^a = x^a y^a$ (for
commutative rings), and aliasing of input and output. They return
`GR_TEST_FAIL` on failure and otherwise the status of the operations.
:::

## Internal functions {#internal}

::: function
slong _gr_ctx_get_real_prec(gr_ctx_t ctx)

Returns the working precision of `ctx`, or 0 if it has none.
:::

::: function
void _gr_length(gr_srcptr x, gr_ctx_t ctx)

Calls the internal `_GR_METHOD_LENGTH` method of the ring (a size
measure of `x`, 0 by default) and discards the result.
:::

::: function
void _gr_ctx_zero_divisor_lock(void)
void _gr_ctx_zero_divisor_unlock(void)

Lock and unlock the global mutex that serialises the recording of zero
divisors by rings pretending to be fields (no-ops in builds without
threads).
:::
