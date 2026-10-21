---
title: Binary quadratic forms and class groups
summary: Reduction of positive definite binary quadratic forms, enumeration of the reduced forms of a discriminant (the class number), composition by Shanks' NUCOMP, and the computation of orders and of the exponent of the form class group by a two-stage method in the style of Pollard and Sutherland.
complexity: Reduction and NUCOMP cost $O(M(\log |D|) \log\log |D|)$ or less per operation; listing the reduced forms costs $O(|D|^{1/2+o(1)})$; the order of an element whose order is $B_1$-smooth apart from one prime below $B_2$ costs $O(\pi(B_1) \log B_1 + \sqrt{B_2})$ compositions.
order: 90
---

A binary quadratic form $f = (a, b, c)$ is the polynomial $ax^2 + bxy +
cy^2$ with integer coefficients; its discriminant is $D = b^2 - 4ac$,
which is congruent to 0 or 1 modulo 4. The group
$\mathrm{SL}_2(\Z)$ acts on forms by linear substitution of the
variables, preserving the discriminant, and two forms in the same orbit
are *equivalent*. Only positive definite forms ($D < 0$, $a > 0$) are
supported by `qfb` and treated here. The standard references are
[@Buell1989], [@Coh1996, Chapter 5] and [@Cox2013].

## Reduction {#reduction}

::: definition {#def-reduced}
A positive definite form $(a, b, c)$ is *reduced* if $|b| \le a \le c$,
and $b \ge 0$ whenever $|b| = a$ or $a = c$.
:::

::: theorem {#thm-reduced}
Every positive definite form is equivalent to exactly one reduced form.
A reduced form of discriminant $D$ has $a \le \sqrt{|D|/3}$.
:::

The bound follows from $|D| = 4ac - b^2 \ge 4a^2 - a^2 = 3a^2$. The
existence and uniqueness are classical [@Coh1996, §5.3].
The reduction algorithm alternates two transformations:

::: algorithm Reduction of a positive definite form {#alg-reduce}
**Input:** $(a, b, c)$ with $a > 0$ and discriminant $D < 0$.

**Output:** the equivalent reduced form.

1. If $|b| > a$: replace $b$ by its residue modulo $2a$ in $(-a, a]$ and
   set $c \gets (b^2 - D)/(4a)$ (a substitution $x \mapsto x + ky$).
2. If $c < a$: exchange $a$ and $c$ and negate $b$ (the substitution
   $(x, y) \mapsto (-y, x)$), and go to step 1.
3. If $|b| = a$ or $a = c$, replace $b$ by $|b|$. Return $(a, b, c)$.
:::

Each exchange in step 2 decreases $a$, and after the first step the
algorithm behaves like a Euclidean algorithm on the coefficients, with
$O(\log(a/\sqrt{|D|}) + 1)$ iterations. `qfb_reduce` implements it with
`fmpz` arithmetic; the result of a composition is already nearly reduced,
so in the class group computations below a reduction costs only a few
steps.

## Reduced forms and the class number {#enumeration}

The primitive forms of discriminant $D$ (those with $\gcd(a, b, c) = 1$)
modulo equivalence form a finite abelian group under composition, the
*form class group* $\mathrm{Cl}(D)$, of order the class number $h(D)$.
By [](#thm-reduced), $h(D)$ is the number of primitive reduced forms, and
these can be listed by searching $|b| \le a \le \sqrt{|D|/3}$ with $b
\equiv D \pmod 2$ and $a \mid (b^2 - D)/4$, $c = (b^2 - D)/(4a) \ge a$.
FLINT does this in two ways, both in $O(|D|^{1/2 + o(1)})$ operations:

- by $b$ (`qfb_reduced_forms`): sieve the factorisations of $(b^2 -
  D)/4$ for all $b$ at once, using the square roots of $D$ modulo each
  small prime $p$ (the $b$ for which $p$ divides $(b^2 - D)/4$ form
  arithmetic progressions), and read off the divisors $a$;
- by $a$ (`qfb_reduced_forms_large`): factor $4a$ for all $a$ by
  sieving, and for each $a$ enumerate the square roots $b$ of $D$ modulo
  $4a$ with `n_sqrtmodn`. This needs only the primes up to $|D|^{1/4}$
  and is used for large $|D|$.

Asymptotically faster methods for $h(D)$, of complexity about $|D|^{1/4}$
(Shanks' baby-step giant-step on the analytic class number formula) or
subexponential (index calculus), are not implemented.

## Composition {#composition}

Gauss composition of two primitive forms $f_1 = (a_1, b_1, c_1)$ and $f_2
= (a_2, b_2, c_2)$ of discriminant $D$ is as follows (Dirichlet's
composition, see [@Cox2013, §3]). Let $g = \gcd(a_1, a_2,
(b_1 + b_2)/2)$ and find $B$ with
$$
B \equiv b_1 \pmod{2a_1/g}, \qquad B \equiv b_2 \pmod{2a_2/g}, \qquad
B^2 \equiv D \pmod{4 a_1 a_2/g^2}.
$$
Then $f_1 f_2$ is the class of $(a_1 a_2/g^2, B, (B^2 - D) g^2/(4 a_1
a_2))$. The coefficients of this form are of the size of $|D|$, and a
reduction brings them back to $\sqrt{|D|}$.

Shanks' NUCOMP [@Sha1989] avoids the large intermediate form. It
computes the same quantities, but stops the extended Euclidean algorithm
that determines $B$ halfway, when the remainders reach $L = \lfloor
|D|^{1/4} \rfloor$, and assembles from the partial cofactors a form that
is already nearly reduced. All intermediate numbers then have about half
the size, which roughly halves the cost of a composition followed by a
reduction. FLINT follows the presentation of Jacobson and van der Poorten
[@JvdP2002] (`qfb_nucomp`), with the specialisation NUDUPL for squaring
(`qfb_nudupl`), and uses `fmpz_xgcd_partial` for the partial Euclidean
algorithm, a variant of
[Lehmer's algorithm](algorithms/integers/gcd#lehmer). Powers are computed by binary powering with NUDUPL, NUCOMP
and a reduction after each step (`qfb_pow_ui`, `qfb_pow`).

## Prime forms {#prime-forms}

For a prime $p$ such that $D$ is a square modulo $4p$, there is a form
$(p, b, c)$ of discriminant $D$ with $0 \le b \le p$, obtained from a
square root of $D$ modulo $p$ (`qfb_prime_form`). These *prime forms*
correspond to the prime ideals of norm $p$, and the prime forms with $p$
below a bound of order $\log^2 |D|$ generate the class group under the
generalised Riemann hypothesis: Bach [@Bach1990] proved that $p \le 6
\log^2 |D|$ suffices; see also [@BD1992].

## Orders and the exponent {#order}

The order of an element $f$ of a finite group can be found by an
analogue of Pollard's $p - 1$ method when it is smooth [@Sut2007].

**Stage 1.** Let $E = \prod_{p \le B_1} p^{\lfloor \log_p B_1 \rfloor}$.
If $f^E = 1$, then the order of $f$ divides $E$, and it is found by
removing the prime factors of $E$ one at a time: for each prime $p$,
raise $f^{E/p^{k}}$ to successive powers of $p$ until the identity
appears. FLINT accumulates $E$ in products of four primes, doubles the
number of primes in each round, and saves intermediate powers so that the
exact order can be found without starting again.

**Stage 2.** If $f^E \ne 1$, the order of $h = f^E$ may still be a single
prime $q \le B_2$. With $m \approx \sqrt{B_2}$, store the baby steps $h^j$
for odd $j < m$ in a hash table and compute the giant steps $h^{2im}$.
Since the table is searched for a form *or its inverse* (equal $a$ and
$|b|$), a match $h^{2im} = h^{\pm j}$ gives $h^{2im \mp j} = 1$, so the
baby steps need only cover half of the residues, and odd $j$ suffice
because $q$ is odd. The candidate $2im \mp j$ is accepted if it is prime.

**Exponent.** The exponent $\lambda$ of $\mathrm{Cl}(D)$ is the least
common multiple of the orders of its elements, a multiple of the order of
every element and a divisor of $h(D)$. Sutherland's Algorithm 8.1
[@Sut2007] computes it from random elements, here prime forms: each new
form is first raised to the exponent $\lambda'$ found so far, so that only
the missing part of its order has to be found, and $\lambda' \gets
\lambda' \cdot \operatorname{ord}(f^{\lambda'})$. With prime forms for all
primes up to Bach's bound, the result is the exponent assuming GRH; with
a fixed number of random primes, it is the exponent with high probability.
When $D$ is not fundamental, a prime form may be imprimitive; it is then
divided by its content, which changes the discriminant to $D/p^2$.

## Use in FLINT {#use}

| Function | Method |
|---|---|
| `qfb_reduce`, `qfb_is_reduced` | [](#alg-reduce) |
| `qfb_reduced_forms`, `qfb_reduced_forms_large` | enumeration by $b$ or by $a$ with sieving |
| `qfb_nucomp`, `qfb_nudupl` | NUCOMP and NUDUPL [@JvdP2002] |
| `qfb_pow_ui`, `qfb_pow`, `qfb_pow_with_root` | binary powering with reduction |
| `qfb_prime_form` | square root of $D$ modulo $p$ |
| `qfb_exponent_element` | stage 1 and stage 2 with `qfb_hash_t` |
| `qfb_exponent`, `qfb_exponent_grh` | Sutherland's Algorithm 8.1 on prime forms |

::: functions qfb_reduce qfb_reduced_forms qfb_nucomp qfb_exponent_element qfb_exponent_grh :::
