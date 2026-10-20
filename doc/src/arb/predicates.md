---
title: Comparisons and predicates
functions: arb_is_positive arb_contains_zero arb_lt arb_equal arb_overlaps arb_contains arb_is_zero arb_is_exact
---

A ball represents an unknown real number in an interval, so a question
such as "is $x < y$?" can have three answers: certainly true, certainly
false, or unknown. FLINT's predicates return an `int` that is 1 when the
answer is certainly true and 0 otherwise.

## Three-valued logic with two values {#logic}

Because 0 means "false or unknown", the negation of a predicate is not
obtained by negating its return value. To branch on the sign of $x$,
test both properties:

```c
if (arb_is_positive(x))
    ...;    /* x > 0 */
else if (arb_is_nonpositive(x))
    ...;    /* x <= 0 */
else
    ...;    /* not known: increase the precision, or give up */
```

The pairs of *pointwise* predicates that are each other's negations are
`arb_is_positive` / `arb_is_nonpositive`, `arb_is_negative` /
`arb_is_nonnegative`, and `arb_is_zero` / `arb_is_nonzero`, and for two
balls `arb_lt` / `arb_ge`, `arb_le` / `arb_gt` and `arb_eq` / `arb_ne`.
The `arb_contains_` predicates are their duals: `arb_contains_zero(x)` is
exactly `!arb_is_nonzero(x)`, `arb_contains_positive(x)` is
`!arb_is_nonpositive(x)`, and so on. They ask whether *some* point of the
ball has the property, which is what an algorithm needs to know before
it may assume the opposite.

## Pointwise predicates and set predicates {#kinds}

Three different kinds of question can be asked about balls.

- *Pointwise predicates* such as `arb_is_positive`, `arb_lt` and `arb_eq`
  are the ball versions of predicates on real numbers: $x < y$ is true if
  $t < u$ for all $t \in x$ and $u \in y$. They are the ones to use when
  the balls stand for unknown numbers. `arb_eq(x, y)` can only be true if
  both are the same exact number.
- *Set predicates* such as `arb_overlaps` and `arb_contains` treat the
  balls as sets. `arb_overlaps(x, y)` returns 1 if the two numbers *might*
  be equal (it is the negation of `arb_ne`), and `arb_contains(x, y)`
  returns 1 if every number in `y` is in `x`, which is how a test program
  checks that an enclosure contains a known value.
- *Representation tests* such as `arb_equal`, `arb_is_exact` and
  `arb_is_finite` look at the midpoint and the radius. `arb_equal(x, x)`
  is 1 for every ball, while `arb_eq(x, x)` is 0 for an inexact ball,
  because two numbers chosen in the same ball need not be equal.

```c example name=arb_predicates
#include <flint/flint.h>
#include <flint/arb.h>

static const char * sign(const arb_t x)
{
    if (arb_is_positive(x))
        return "positive";
    else if (arb_is_negative(x))
        return "negative";
    else if (arb_is_zero(x))
        return "zero";
    else
        return "unknown";
}

int main(void)
{
    arb_t a, b, c;
    arb_init(a); arb_init(b); arb_init(c);

    /* a = [4 +/- 1] = [3, 5], b = [6 +/- 1] = [5, 7], c = [6 +/- 2] = [4, 8] */
    arb_set_ui(a, 4); mag_set_ui(arb_radref(a), 1);
    arb_set_ui(b, 6); mag_set_ui(arb_radref(b), 1);
    arb_set_ui(c, 6); mag_set_ui(arb_radref(c), 2);

    flint_printf("a < b: %d   a <= b: %d   b < a: %d\n",
        arb_lt(a, b), arb_le(a, b), arb_lt(b, a));
    flint_printf("a < c: %d   c <= a: %d   overlaps: %d\n",
        arb_lt(a, c), arb_le(c, a), arb_overlaps(a, c));
    flint_printf("c contains b: %d   b contains c: %d\n",
        arb_contains(c, b), arb_contains(b, c));

    arb_sub(c, b, a, 53);
    flint_printf("b - a: %s\n", sign(c));
    arb_sub(c, a, a, 53);
    flint_printf("a - a: %s, eq(a, a): %d, equal(a, a): %d\n",
        sign(c), arb_eq(a, a), arb_equal(a, a));
    arb_zero(c);
    flint_printf("0: %s, eq(0, 0): %d\n", sign(c), arb_eq(c, c));

    arb_indeterminate(c);
    flint_printf("nan: positive %d, nonpositive %d, contains_zero %d, overlaps %d\n",
        arb_is_positive(c), arb_is_nonpositive(c), arb_contains_zero(c),
        arb_overlaps(c, a));

    arb_clear(a); arb_clear(b); arb_clear(c);
    flint_cleanup();
    return 0;
}
```
```output
a < b: 0   a <= b: 1   b < a: 0
a < c: 0   c <= a: 0   overlaps: 1
c contains b: 1   b contains c: 0
b - a: unknown
a - a: unknown, eq(a, a): 0, equal(a, a): 1
0: zero, eq(0, 0): 1
nan: positive 0, nonpositive 0, contains_zero 1, overlaps 1
```

The balls $a = [3, 5]$ and $b = [5, 7]$ touch, so $a \le b$ is certain
but $a < b$ is not. $b - a = [0, 4]$ is certainly nonnegative but not
certainly positive. An indeterminate ball satisfies no pointwise
predicate and every set predicate that asks whether something is
possible.

## NaN and infinities {#special}

An indeterminate ball (NaN midpoint) may represent any number, so all
pointwise predicates and `arb_contains(y, x)` for an ordinary `y` return
0, while the "possibly" predicates (`arb_contains_zero`,
`arb_contains_positive`, ..., `arb_overlaps`, and `arb_contains(x, y)`
with the NaN ball as `x`) return 1. Infinities are points of the extended
real line: $[+\infty] > [10^{100}]$ is certain, $[-\infty \pm 3]$ and
$[-\infty]$ compare equal, and a ball with an infinite radius is the whole
line, about which nothing is certain.

## Deciding by increasing precision {#precision}

When a predicate returns 0 for both a property and its negation, the
usual remedy is to recompute at a higher precision. This terminates if
the answer is determined by a nonzero quantity, but loops forever when
the quantity is exactly zero, since ball arithmetic cannot prove that a
number is zero unless the computation is exact (see
[](guide/types/ball-arithmetic#guess-and-verify)). Programs should limit
the precision, and treat the remaining case explicitly.
