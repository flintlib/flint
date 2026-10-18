---
title: qfb.h
header: qfb.h
description: Positive definite binary quadratic forms — reduction, enumeration of the reduced forms of a discriminant, NUCOMP composition, powering, and the order and exponent of the form class group.
category: Number theory
order: 6
---

This module implements binary quadratic forms $f(x, y) = ax^2 + bxy +
cy^2$ with integer coefficients, written $(a, b, c)$, of discriminant $D
= b^2 - 4ac$. It supports positive definite forms ($D < 0$, $a > 0$):
reduction to the unique reduced form of an equivalence class, the list of
all reduced primitive forms of a discriminant (whose number is the class
number $h(D)$), composition by Shanks' NUCOMP algorithm, powering, prime
forms, and the computation of the order of a form and of the exponent of
the class group. Indefinite forms are not supported.

The guide [](qfb/class-group) shows how to compute in the class group,
and [](algorithms/number-theory/binary-quadratic-forms) describes the
algorithms. The original module was written by William Hart; the
documentation was later revised by Håvard Damm-Johnsen.

::: summary algorithms/number-theory/binary-quadratic-forms :::

## Types and memory management {#types}

::: type
qfb
qfb_t

A binary quadratic form: a struct with three `fmpz_t` fields `a`, `b`,
`c`. A `qfb_t` is an array of length one of `qfb`.
:::

::: function
void qfb_init(qfb_t q)
void qfb_clear(qfb_t q)

Initialise $q$ (to $(0, 0, 0)$) and free its memory.
:::

::: function
void qfb_array_clear(qfb ** forms, slong num)

Clears the first `num` forms of the array `*forms` and frees the array,
as allocated by `qfb_reduced_forms`.
:::

## Basic operations {#basic}

::: function
void qfb_set(qfb_t f, qfb_t g)

Sets $f$ to a copy of $g$.
:::

::: function
int qfb_equal(qfb_t f, qfb_t g)

Returns 1 if $f$ and $g$ have the same coefficients, otherwise 0. For
reduced forms this tests equivalence.
:::

::: function
void qfb_print(qfb_t q)

Prints $q$ as `(a, b, c)`.
:::

::: function
void qfb_discriminant(fmpz_t D, qfb_t f)

Sets $D$ to the discriminant $b^2 - 4ac$ of $f = (a, b, c)$.
:::

::: function
int qfb_is_primitive(qfb_t f)

Returns 1 if $\gcd(a, b, c) = 1$, otherwise 0.
:::

## Reduction and enumeration {#reduction}

A positive definite form is *reduced* if $|b| \le a \le c$, and $b \ge 0$
if $|b| = a$ or $a = c$. Every positive definite form is equivalent under
$\mathrm{SL}_2(\Z)$ to exactly one reduced form.

::: function
void qfb_reduce(qfb_t r, qfb_t f, fmpz_t D)

Sets $r$ to the reduced form equivalent to the positive definite form
$f$, whose discriminant must be $D$ (used to recompute $c$). Aliasing is
allowed.
:::

::: function
int qfb_is_reduced(qfb_t r)

Returns 1 if the positive definite form $r$ is reduced, otherwise 0.
:::

::: function
slong qfb_reduced_forms(qfb ** forms, slong d)
slong qfb_reduced_forms_large(qfb ** forms, slong d)

Compute all reduced primitive positive definite forms of discriminant $d
< 0$, allocate an array for them, store it in `*forms` and return their
number, the class number $h(d)$. The forms are in no particular order.
The array must be freed with `qfb_array_clear(forms, num)`. Returns 0
(with `*forms` set to `NULL`) if $d \not\equiv 0, 1 \pmod 4$, and throws
for $d \ge 0$.

::: details
`qfb_reduced_forms` loops over $b$ and sieves the factorisations of $(b^2
- d)/4$ with the square roots of $d$ modulo small primes; the forms are
the factorisations into $ac$ with $|b| \le a \le c$. It needs the primes
up to about $\sqrt{|d|}$ and calls `qfb_reduced_forms_large` when that
would require too large a table. `qfb_reduced_forms_large` loops over $a
\le \sqrt{|d|/3}$ instead, with the square roots of $d$ modulo $4a$, and
needs only the primes up to $|d|^{1/4}$. Both take $O(|d|^{1/2 + o(1)})$
time and memory.
:::
:::

## The class group {#group}

The primitive forms of discriminant $D$ modulo equivalence form a finite
abelian group under Gauss composition, the *form class group*, whose
identity is the class of the principal form and in which the inverse of
the class of $(a, b, c)$ is the class of $(a, -b, c)$.

::: function
void qfb_principal_form(qfb_t f, fmpz_t D)

Sets $f$ to the principal form of discriminant $D$, $(1, 0, -D/4)$ for $D
\equiv 0 \pmod 4$ and $(1, 1, (1 - D)/4)$ for $D \equiv 1 \pmod 4$.
:::

::: function
int qfb_is_principal_form(qfb_t f, fmpz_t D)

Returns 1 if $f$ is the reduced principal form of discriminant $D$ (the
identity of the class group, for reduced $f$), otherwise 0.
:::

::: function
void qfb_inverse(qfb_t r, qfb_t f)

Sets $r$ to the reduced inverse of the reduced form $f$: $(a, -b, c)$, or
$f$ itself if $a = b$ or $a = c$ (such forms have order at most 2).
:::

::: function
void qfb_prime_form(qfb_t r, fmpz_t D, fmpz_t p)

Sets $r$ to a form $(p, b, c)$ of discriminant $D$ with $0 \le b \le p$,
for a prime $p$ such that $D$ is a square modulo $4p$ (that is, $p$
splits or ramifies). The result is not specified otherwise. If $p^2 \mid
D$ the form may be imprimitive.
:::

::: function
void qfb_nucomp(qfb_t r, const qfb_t f, const qfb_t g, fmpz_t D, fmpz_t L)

Sets $r$ to a form in the class of the composition of $f$ and $g$, both
of discriminant $D$, by Shanks' NUCOMP algorithm in the version of
Jacobson and van der Poorten [@JvdP2002], with a partial
[Lehmer](algorithms/integers/gcd#lehmer) extended GCD (`fmpz_xgcd_partial`).
The parameter must be $L = \lfloor |D|^{1/4} \rfloor$. The result is *nearly reduced* (its
coefficients are of the size of $\sqrt{|D|}$) but must be passed to
`qfb_reduce` to obtain the reduced form. $f$ must be primitive.
:::

::: function
void qfb_nudupl(qfb_t r, const qfb_t f, fmpz_t D, fmpz_t L)

As `qfb_nucomp` for the composition of $f$ with itself (NUDUPL), which
needs fewer operations. $f$ must be primitive.
:::

::: function
void qfb_pow_ui(qfb_t r, qfb_t f, fmpz_t D, ulong exp)
void qfb_pow(qfb_t r, qfb_t f, fmpz_t D, fmpz_t exp)
void qfb_pow_with_root(qfb_t r, qfb_t f, fmpz_t D, fmpz_t exp, fmpz_t L)

Set $r$ to the reduced form in the class of $f^{\mathit{exp}}$, by binary
powering with NUDUPL and NUCOMP and a reduction after each step. $D$ must
be the discriminant of $f$, which must be primitive. For exponent 0 the
result is the principal form and for exponent 1 a copy of $f$ (not
reduced if $f$ is not). `qfb_pow_with_root` takes $L = \lfloor
|D|^{1/4} \rfloor$ as an argument instead of computing it. A negative
`fmpz_t` exponent is allowed: the result is then the power $-\mathit{exp}$
of the inverse form (`qfb_inverse`).
:::

## Order and exponent {#exponent}

::: function
int qfb_exponent_element(fmpz_t exponent, qfb_t f, fmpz_t n, ulong B1, ulong B2_sqrt)

Tries to find the order of the class of the reduced primitive form $f$ in
the class group of discriminant $n$. Returns 1 and sets `exponent` to the
order on success, and returns 0 if the order is not of the form $s \cdot
q$ with $s$ a product of prime powers below about `B1` and $q$ a prime
below about `B2_sqrt` squared.

::: details
The function raises $f$ to products of the primes up to `B1`, in batches
of increasing size, until the principal form appears, then determines
the exact power of each prime (stage 1). If this fails it looks for a
single larger prime factor with a baby-step giant-step search (stage 2),
which stores about `B2_sqrt` forms in a hash table. See
[](algorithms/number-theory/binary-quadratic-forms#order).
:::
:::

::: function
int qfb_exponent(fmpz_t exponent, fmpz_t n, ulong B1, ulong B2_sqrt, slong c)

Tries to compute the exponent of the class group of discriminant $n <
0$, the least common multiple of the orders of its elements, as the least
common multiple of the orders of prime forms for small primes. Returns 1
and sets `exponent` on success, and 0 if `qfb_exponent_element` fails.
The parameter $c$ controls the number of prime forms tried; the result
is the exponent with high probability (heuristically, at least $1 -
2^{-c}$), but it is not proved. The bounds `B1` and `B2_sqrt` are as for
`qfb_exponent_element`.

::: details
The algorithm is Algorithm 8.1 of [@Sut2007]: the prime forms are raised
to the exponent found so far, so that each further order computation
only finds the missing part. Non-fundamental discriminants are handled by
dividing imprimitive prime forms by their content.
:::
:::

::: function
int qfb_exponent_grh(fmpz_t exponent, fmpz_t n, ulong B1, ulong B2_sqrt)

As `qfb_exponent`, but uses the prime forms of all primes $p < 6 \log^2
|n|$ that split or ramify. Assuming the generalised Riemann hypothesis,
these generate the class group [@Bach1990], [@BD1992], so the result is
the exponent of the class group if the function returns 1.
:::

## Hash tables {#hash}

These functions implement a simple hash table of forms with open
addressing, used for the baby-step giant-step stage of
`qfb_exponent_element`. Forms are hashed by $a$ modulo $2^{\mathit{depth}}$.

::: type
qfb_hash_t

An entry of the table: the form `q`, a second form `q2` and an integer
`iter`.
:::

::: function
qfb_hash_t * qfb_hash_init(slong depth)
void qfb_hash_clear(qfb_hash_t * qhash, slong depth)

Allocate a table of $2^{\mathit{depth}}$ empty entries (with `q->a` zero)
and free it.
:::

::: function
void qfb_hash_insert(qfb_hash_t * qhash, qfb_t q, qfb_t q2, slong iter, slong depth)

Stores $q$, `q2` (unless it is `NULL`) and `iter` in the first free entry
at or after the position given by the hash of $q$. An exception is
raised if the table is full.
:::

::: function
slong qfb_hash_find(qfb_hash_t * qhash, qfb_t q, slong depth)

Returns the index of an entry containing $q$ or its inverse (the same $a$
and $|b|$), or $-1$ if there is none (also when the table is full and
does not contain it).
:::
