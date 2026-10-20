---
title: ca_ext.h
header: ca_ext.h
description: Extension numbers, the generators of the fields in which exact real and complex numbers are represented: algebraic numbers, constants and function values, with cached enclosures.
category: Exact real and complex numbers
order: 6
---

An *extension number* is a fixed real or complex number $a$ that serves
as a generator of the fields $\mathbb{Q}(a_1, \ldots, a_n)$ in which
`ca_t` numbers are represented (see [](ca) and [](ca_field)). A
`ca_ext_t` holds one of the following:

- an algebraic number, given in canonical form by a `qqbar_t` (see
  [](qqbar)); for example $i$ is the root of $x^2 + 1$ with positive
  imaginary part;
- a function value $f(x_1, \ldots, x_n)$, where $f$ is one of the
  builtin functions identified by a `calcium_func_code` (see
  [](calcium)) and $x_1, \ldots, x_n$ are `ca_t` numbers, for example
  $e^{x}$ or $x^y$;
- a builtin constant such as $\pi$, which is the case $n = 0$ of a
  function value.

An extension number caches data for repeated computations with the
same number: a numerical enclosure that is refined on demand and, for
algebraic numbers, the number field data used for arithmetic
in $\mathbb{Q}(a)$. Because the arguments of a function value are
themselves elements of fields with their own extension numbers,
extension numbers form a directed acyclic graph, and their operations
need the context object `ca_ctx_t`.

Users do not normally create extension numbers. The intended way to
work with an extension number $a$ is through the `ca_t` number $1 \cdot
a$, the generator of the field $\mathbb{Q}(a)$; given a `ca_t`,
`ca_is_gen_as_ext` returns its extension number if it is a generator.
The context keeps every extension number in a cache without
duplicates, so that structurally equal extension numbers are the same
object and can be compared by pointer.

## Types and macros {#types}

::: type
ca_ext_struct
ca_ext_t

An extension number. The structure holds the head (a
`calcium_func_code`), a hash value, the depth, and either a `qqbar_t`
together with a pointer to an `nf_struct` for arithmetic in
$\mathbb{Q}(a)$ (when the head is `CA_QQBar`), or the array of
arguments, their number, a cached enclosure with its precision and an
optional cached `qqbar_t` value (for a function value). A `ca_ext_t` is
an array of one structure, so it is passed by reference.
:::

::: type
ca_ext_ptr
ca_ext_srcptr

Aliases for `ca_ext_struct *` and `const ca_ext_struct *`.
:::

::: macro
CA_EXT_HEAD(x)
CA_EXT_HASH(x)
CA_EXT_DEPTH(x)

Access the head, the hash value and the depth of `x`. The head is
`CA_QQBar` for an algebraic number and the function code (`CA_Exp`,
`CA_Pi`, ...) for a function value. The depth is 0 for an algebraic
number or a constant, and one more than the largest depth of the
extension numbers involved in the arguments for a function value.
:::

::: macro
CA_EXT_IS_QQBAR(x)
CA_EXT_QQBAR(x)
CA_EXT_QQBAR_NF(x)

`CA_EXT_IS_QQBAR` tests whether `x` is an algebraic number in canonical
form. For such an `x`, the other macros give its `qqbar_t` value and a
pointer to the number field object `nf_struct` of $\mathbb{Q}(x)$ (see
[](nf)).
:::

::: macro
CA_EXT_FUNC_ARGS(x)
CA_EXT_FUNC_NARGS(x)
CA_EXT_FUNC_ENCLOSURE(x)
CA_EXT_FUNC_PREC(x)

For a function value or constant `x`: the array of arguments (a
`ca_ptr`), the number of arguments, the cached enclosure (an `acb_t`)
and the precision at which it was computed (0 if none has been).
:::

## Memory management {#memory}

::: function
void ca_ext_init_qqbar(ca_ext_t res, const qqbar_t x, ca_ctx_t ctx)

Initialises `res` as the algebraic number `x`, and initialises the
number field data of $\mathbb{Q}(x)$ from the minimal polynomial of
`x`.
:::

::: function
void ca_ext_init_const(ca_ext_t res, calcium_func_code func, ca_ctx_t ctx)

Initialises `res` as the constant identified by `func`, for example
`CA_Pi` for $\pi$ or `CA_Euler` for $\gamma$.
:::

::: function
void ca_ext_init_fx(ca_ext_t res, calcium_func_code func, const ca_t x, ca_ctx_t ctx)
void ca_ext_init_fxy(ca_ext_t res, calcium_func_code func, const ca_t x, const ca_t y, ca_ctx_t ctx)
void ca_ext_init_fxn(ca_ext_t res, calcium_func_code func, ca_srcptr x, slong nargs, ca_ctx_t ctx)

Initialises `res` as the function value $f(x)$, $f(x, y)$ or $f(x_1,
\ldots, x_n)$ with $n$ = `nargs`, where $f$ is identified by `func`
(for example `CA_Exp` for $e^x$ or `CA_Pow` for $x^y$). The arguments
are copied. Nothing is evaluated or simplified, and the arguments
should be numbers for which $f$ is defined.
:::

::: function
void ca_ext_init_set(ca_ext_t res, const ca_ext_t x, ca_ctx_t ctx)

Initialises `res` as a copy of `x` (without the cached data).
:::

::: function
void ca_ext_clear(ca_ext_t res, ca_ctx_t ctx)

Clears `res`, freeing its arguments or number field data and caches.
:::

## Structure {#structure}

::: function
slong ca_ext_nargs(const ca_ext_t x, ca_ctx_t ctx)

Returns the number of arguments of `x`: 0 for an algebraic number and
for a constant such as $\pi$.
:::

::: function
void ca_ext_get_arg(ca_t res, const ca_ext_t x, slong i, ca_ctx_t ctx)

Sets `res` to the argument with index `i` (starting from 0) of `x`.
Throws an exception if `x` is an algebraic number or `i` is out of
range.
:::

::: function
ulong ca_ext_hash(const ca_ext_t x, ca_ctx_t ctx)

Returns the hash of the structural representation of `x`, computed at
initialisation from the minimal polynomial (for an algebraic number, so
that conjugates have equal hashes) or from the head and the hashes of
the arguments.
:::

::: function
int ca_ext_equal_repr(const ca_ext_t x, const ca_ext_t y, ca_ctx_t ctx)

Returns whether `x` and `y` are structurally equal: the same algebraic
number, or the same head with arguments of equal representation (see
`ca_equal_repr`).
:::

::: function
int ca_ext_cmp_repr(const ca_ext_t x, const ca_ext_t y, ca_ctx_t ctx)

Compares `x` and `y` in a total structural order, returning $-1$, 0 or
1. The result says nothing about the values of the numbers.

::: details
Algebraic numbers come first, ordered by degree, then by the
coefficients of the minimal polynomial (by absolute value from the
leading coefficient down, then by sign), then by real and imaginary
parts. Function values follow, ordered by depth, then by head, then by
number of arguments, then by the arguments with `ca_cmp_repr`. Ordering
by depth makes $x$ precede $f(x)$; the fields sort their generators by
decreasing order, so that a lexicographic monomial order eliminates the
more complex generators first (see [](algorithms/exact/field-towers)).
:::
:::

## Input and output {#io}

::: function
void ca_ext_print(const ca_ext_t x, ca_ctx_t ctx)

Prints a description of `x` to standard output, without a newline:
`Algebraic` followed by a decimal approximation of an algebraic number
(or `I`), or the function name followed by its arguments in parentheses,
printed with `ca_print`.
:::

## Numerical evaluation {#numerical}

::: function
void ca_ext_get_acb_raw(acb_t res, ca_ext_t x, slong prec, ca_ctx_t ctx)

Sets `res` to an enclosure of the value of `x` computed with working
precision `prec`, without adaptive refinement. An algebraic number is
refined to the requested precision with `qqbar_cache_enclosure`; for a
function value, the arguments are evaluated with `ca_get_acb_raw` and
the corresponding [](acb) function is applied. The enclosure is cached
in `x` and reused for any later request at a precision not larger than
the cached one.
:::

## Cache {#cache}

::: type
ca_ext_cache_struct
ca_ext_cache_t

A set of structurally distinct extension numbers: an array of pointers
to heap-allocated `ca_ext_struct` objects together with a hash table
for lookup. Each context object contains one such cache.
:::

::: function
void ca_ext_cache_init(ca_ext_cache_t cache, ca_ctx_t ctx)
void ca_ext_cache_clear(ca_ext_cache_t cache, ca_ctx_t ctx)

Initialises `cache` as empty, or clears it, freeing all the extension
numbers it holds.
:::

::: function
ca_ext_ptr ca_ext_cache_insert(ca_ext_cache_t cache, const ca_ext_t x, ca_ctx_t ctx)

Returns a pointer to the extension number in `cache` that is
structurally equal to `x`, inserting a copy of `x` first if there is
none. The returned object is owned by the cache and lives as long as
the context.
:::
