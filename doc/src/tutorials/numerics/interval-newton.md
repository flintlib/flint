---
title: Proving and computing a root with the interval Newton method
summary: A complete rigorous computation from scratch: the interval Newton operator, a proof of existence and uniqueness of a root, exclusion of roots by bisection, and quadratic refinement to 1000 digits with precision doubling, for the omega constant x e^x = 1.
order: 6
---

This tutorial solves an equation rigorously from scratch, using nothing
but ball arithmetic and one theorem. The equation is
$$
x e^x = 1,
$$
whose unique real solution is the omega constant
$\Omega = W(1) = 0.5671432904\ldots$, the value of the Lambert $W$
function at 1. We will (1) prove that $[0, 2]$ contains exactly one
solution, locating it in a subinterval; (2) refine it to 1000 digits with
quadratic convergence; (3) check the result against the library function
`arb_lambertw`. The same program works for any function whose value and
derivative can be evaluated in ball arithmetic.

## The interval Newton operator {#operator}

Let $f$ be continuously differentiable on an interval $X = [a, b]$, let
$m \in X$ (usually the midpoint) and let $F'$ be an interval containing
$f'(x)$ for all $x \in X$. The *interval Newton operator* is
$$
N(X) = m - \frac{f(m)}{F'},
$$
the set of all $m - f(m)/d$ with $d \in F'$, assuming $0 \notin F'$.

::: theorem Interval Newton {#thm-newton}
Suppose that $0 \notin F'$.

1. Every root of $f$ in $X$ lies in $N(X)$.
2. If $N(X) \cap X = \emptyset$, then $f$ has no root in $X$.
3. If $N(X) \subseteq X$, then $f$ has exactly one root in $X$.
:::

::: proof
(1) If $\xi \in X$ and $f(\xi) = 0$, the mean value theorem gives
$0 = f(m) + f'(\eta)(\xi - m)$ with $\eta$ between $m$ and $\xi$, so
$\xi = m - f(m)/f'(\eta) \in N(X)$ since $f'(\eta) \in F'$.
(2) follows from (1).
(3) Since $f'$ is continuous and does not vanish on $X$, $f$ is strictly
monotonic on $X$, so it has at most one root there. For existence,
assume $f' > 0$ on $X$ (the other case is symmetric). The mean value
theorem gives $f(a) = f(m) - f'(\eta)(m - a)$ for some $\eta \in X$.
Since $m - f(m)/f'(\eta) \in N(X) \subseteq X$, we have
$m - f(m)/f'(\eta) \ge a$, that is $f(m) \le f'(\eta)(m - a)$, so
$f(a) \le 0$. In the same way $f(b) \ge 0$, and the intermediate value
theorem gives a root.
:::

The theorem goes back to Moore [@Moo1966]; see [@Neu1990], [@Tuc2011] or
[@Rum2010] for its many variants, including the multivariate Krawczyk
operator. In ball arithmetic, we cannot compute $N(X)$ exactly, but we
can compute a ball $\tilde N$ that contains it: evaluate $f(m)$ at the
exact point $m$ (as a ball, since $f(m)$ is rounded), $F' = f'(X)$ on the
ball $X$, and the quotient and difference. Then $\tilde N \subseteq X$
implies $N(X) \subseteq X$, and $\tilde N \cap X = \emptyset$ implies
$N(X) \cap X = \emptyset$, so both conclusions of the theorem can be
drawn from the computed enclosure. Moreover, by (1), replacing $X$ by
$\tilde N \cap X$ never loses the root, and repeating this converges
quadratically: the width of $N(X)$ is about $|f(m)| \cdot
\operatorname{width}(1/F') = O(\operatorname{width}(X)^2)$, since $f(m) =
O(\operatorname{width}(X))$ near the root.

## Wide balls {#wide-balls}

The operator needs a reasonably narrow starting interval. Ball arithmetic
represents $[a, b]$ as a midpoint and a radius, and operations on wide
balls overestimate. For $X = [0, 1]$, the program below computes
$f'(X) = (X + 1) e^X$ as the product of $[1, 2]$ and an enclosure of
$[1, e]$; a product of balls $[m_1 \pm r_1][m_2 \pm r_2]$ has radius
$|m_1| r_2 + |m_2| r_1 + r_1 r_2$, which gives about $[0.14, 5.44]$
instead of $[1, 2e]$. Then $N([0, 1])$ is about $[0.53, 1.75]$, which is
not contained in $[0, 1]$, and the test fails although the root is
there. The cure is to bisect: on smaller intervals the overestimation is
smaller, and the test succeeds. The function `search` below bisects
$[0, 2]$ until each part is either proved to contain no root or proved
to contain exactly one.

## The program {#program}

```c example name=numerics_interval_newton
#include <flint/arb.h>

/* f(x) = x exp(x) - 1 and f'(x) = (x + 1) exp(x) */
void f_eval(arb_t y, const arb_t x, slong prec)
{
    arb_exp(y, x, prec);
    arb_mul(y, y, x, prec);
    arb_sub_ui(y, y, 1, prec);
}

void fprime_eval(arb_t y, const arb_t x, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_exp(t, x, prec);
    arb_add_ui(y, x, 1, prec);
    arb_mul(y, y, t, prec);
    arb_clear(t);
}

/* N = m - f(m) / f'(X); returns 0 if f'(X) contains zero */
int newton_operator(arb_t N, const arb_t X, slong prec)
{
    arb_t m, fm, fpX;
    int ok;
    arb_init(m); arb_init(fm); arb_init(fpX);

    arb_get_mid_arb(m, X);              /* exact midpoint */
    f_eval(fm, m, prec);
    fprime_eval(fpX, X, prec);
    ok = !arb_contains_zero(fpX);
    arb_div(N, fm, fpX, prec);
    arb_sub(N, m, N, prec);

    arb_clear(m); arb_clear(fm); arb_clear(fpX);
    return ok;
}

/* search [a, b] for roots; returns 1 and sets X if one root is proved */
int search(arb_t X, double a, double b)
{
    arb_t N;
    int found = 0;
    arb_init(N);

    arb_set_d(X, (a + b) / 2);
    mag_set_d(arb_radref(X), (b - a) / 2);    /* X = [a, b] (or slightly wider) */

    if (newton_operator(N, X, 53) && arb_contains(X, N))
    {
        flint_printf("[%.4f, %.4f]: exactly one root\n", a, b);
        found = 1;
    }
    else if (!arb_overlaps(N, X))
    {
        flint_printf("[%.4f, %.4f]: no root\n", a, b);
    }
    else
    {
        flint_printf("[%.4f, %.4f]: undecided, bisect\n", a, b);
        found = search(X, a, (a + b) / 2);
        if (!found)
            found = search(X, (a + b) / 2, b);
        else
            search(N, (a + b) / 2, b);     /* report the rest, X is kept */
    }

    arb_clear(N);
    return found;
}

int main(void)
{
    arb_t X, N, W;
    slong prec, target = 3330, step;

    arb_init(X); arb_init(N); arb_init(W);

    if (!search(X, 0.0, 2.0))
        return 1;

    /* iterate X = N(X) intersected with X, doubling the precision */
    for (step = 1; ; step++)
    {
        /* twice the current accuracy, plus guard bits */
        prec = FLINT_MAX(64, 2 * arb_rel_accuracy_bits(X) + 30);
        prec = FLINT_MIN(prec, target + 30);
        newton_operator(N, X, prec);
        arb_intersection(X, X, N, prec);
        flint_printf("step %wd: prec = %4wd, accuracy = %4wd bits\n",
            step, prec, arb_rel_accuracy_bits(X));
        if (arb_rel_accuracy_bits(X) >= target)
            break;
    }

    arb_printn(X, 1000, ARB_STR_CONDENSE * 10);
    flint_printf("\n");

    arb_one(W);
    arb_lambertw(W, W, 0, target + 30);
    flint_printf("overlaps W(1): %d\n", arb_overlaps(X, W));

    arb_clear(X); arb_clear(N); arb_clear(W);
    flint_cleanup();
    return 0;
}
```

```output
[0.0000, 2.0000]: undecided, bisect
[0.0000, 1.0000]: undecided, bisect
[0.0000, 0.5000]: no root
[0.5000, 1.0000]: undecided, bisect
[0.5000, 0.7500]: exactly one root
[0.7500, 1.0000]: no root
[1.0000, 2.0000]: undecided, bisect
[1.0000, 1.5000]: no root
[1.5000, 2.0000]: no root
step 1: prec =   64, accuracy =    5 bits
step 2: prec =   64, accuracy =   15 bits
step 3: prec =   64, accuracy =   39 bits
step 4: prec =  108, accuracy =   87 bits
step 5: prec =  204, accuracy =  184 bits
step 6: prec =  398, accuracy =  380 bits
step 7: prec =  790, accuracy =  772 bits
step 8: prec = 1574, accuracy = 1558 bits
step 9: prec = 3146, accuracy = 3129 bits
step 10: prec = 3360, accuracy = 3357 bits
[0.5671432904{...980 digits...}6227976096 +/- 4.92e-1001]
overlaps W(1): 1
```

The output is a proof: $f$ has no root in $[0, 0.5]$, $[0.75, 1]$ and
$[1, 2]$, and exactly one in $[0.5, 0.75]$, and that root is enclosed in
the final ball. (On $[1, 2]$ even the exact operator does not decide:
$N([1, 2]) \approx [0.45, 1.24]$ overlaps $[1, 2]$, so one more
bisection is needed.)

## Notes on the precision {#precision}

The accuracy roughly doubles in each step, as expected from quadratic
convergence: $5, 15, 39, 87, 184, \ldots$ bits. The working precision
only needs to be somewhat larger than the accuracy that the step can
produce, so the program sets it to twice the current accuracy plus 30
guard bits. With this schedule, the total cost is dominated by the last
step or two: the cost of a step at precision $p$ is at least linear in
$p$, and the precisions form a geometric sequence. A naive program that
runs every step at the full target precision would cost about
$\log_2(\text{target})$ times as much.

Three details keep the iteration correct and convergent:

- The midpoint $m$ is taken exactly (`arb_get_mid_arb`), and $f(m)$ is
  evaluated on this exact point, so its radius is only the rounding
  error. Evaluating $f$ on the whole ball $X$ instead would give a ball
  of width proportional to the width of $X$, and no convergence.
- The intersection with the previous $X$ guarantees that the enclosure
  never grows. In exact arithmetic $N(X) \subseteq X$ once the method
  converges; with rounding, the intersection also removes overestimation.
- The stopping criterion uses the accuracy of the enclosure, not the
  size of the last correction. The result is therefore a proved
  enclosure regardless of when the iteration is stopped.

The library uses the same ideas internally: `arb_fmpz_poly_refine_root_arb`
refines a real root of an integer polynomial with interval Newton steps
(with ordinary Newton steps and bisection as fallbacks);
`arb_calc_refine_root_newton` uses the closely related error bound of
[](tutorials/numerics/real-roots#refinement); and `acb_mat_eig_simple`
certifies eigenvalues with a matrix version of the method due to Rump
[@Rum2010].
