---
title: Exponential and logarithm
functions: padic_exp padic_log
---

The $p$-adic exponential and logarithm are defined by their power series

$$
\exp(x) = \sum_{i \ge 0} \frac{x^i}{i!}, \qquad
\log(x) = \sum_{i \ge 1} (-1)^{i - 1} \frac{(x - 1)^i}{i},
$$

which converge for $\operatorname{ord}_p(x) > 1/(p - 1)$ and for
$\operatorname{ord}_p(x - 1) > 0$, respectively. For $x \in \mathbf{Q}_p$
this means $\operatorname{ord}_p(x) \ge 1$ ($\ge 2$ if $p = 2$) for the
exponential and $x \equiv 1 \pmod p$ for the logarithm. On these
domains they are inverse isomorphisms between $p\mathbf{Z}_p$ (or
$4\mathbf{Z}_2$) and $1 + p\mathbf{Z}_p$ (or $1 + 4\mathbf{Z}_2$). The
functions return 0 and leave the output unchanged when the series
diverges, and otherwise return 1 and set the output to the value reduced
modulo $p^N$. Since the result of the exponential is a unit and that of
the logarithm has valuation at least that of $x - 1$, the precision $N$
of the output is also, up to the valuation, its number of correct
digits.

::: summary algorithms/p-adic/exp-log :::

## Algorithms {#algorithms}

Both series are truncated after the first term whose valuation, and that
of every later term, is at least $N$ (`_padic_exp_bound` and
`_padic_log_bound`). The tail bounds use
$\operatorname{ord}_p(i!) \le (i - 1)/(p - 1)$ and
$\operatorname{ord}_p(i) \le \log_p i$; for an argument of valuation $v$
this gives about $N(p - 1)/((p - 1)v - 1)$ terms for the exponential
and about $N/v$ for the logarithm.

- *Rectangular splitting* evaluates the truncated series with about
  $2\sqrt n$ multiplications of full-size integers modulo $p^{N + k}$,
  the $k$ extra digits absorbing the powers of $p$ in the denominators
  $i!$ or $i$. It is used by `padic_exp` for $N < 1024$ and by
  `padic_log` for $N < 512/b$, where $b$ is the bit length of $p$.
- The *balanced* algorithms split the argument into pieces
  $x = x_0 + x_1 + \cdots$ with $x_j$ having valuation at least $2^j$
  and $2^j$ digits, and use $\exp(x) = \prod \exp(x_j)$, respectively
  $\log(x) = \sum \log(1 - y_j)$ for a matching factorisation
  $x = \prod (1 - y_j)$. Each short series is summed exactly by binary
  splitting. This is the $p$-adic analogue of Brent's bit-burst method
  and costs $O(M(N \log p) \log^2 N)$, quasi-linear in $N$.
- *Satoh's method* (`padic_log_satoh`) uses
  $\log(x) = p^{-k} \log(x^{p^k})$: since $x^{p^k} - 1$ has valuation
  at least $k + 1$, the series for the right-hand side has about
  $N/(k + 1)$ terms, and with $k = \lfloor \sqrt N \rfloor$ the power
  (about $k \log_2 p$ squarings) and the shorter series together cost
  $O(\sqrt N \log p)$ multiplications [@SST2003]. It is
  never selected automatically, and for $N < 16$ it falls back to
  rectangular splitting.

## Example {#example}

With $p = 5$ and precision 20, the logarithm inverts the exponential at
$15 = 3 \cdot 5$, by each of the three algorithms. Neither series
converges at the unit 3.

```c example name=padic_exp_log
#include <flint/flint.h>
#include <flint/padic.h>

int main(void)
{
    fmpz_t p;
    padic_ctx_t ctx;
    padic_t x, y, z;
    int ok;

    fmpz_init_set_ui(p, 5);
    padic_ctx_init(ctx, p, 0, 20, PADIC_TERSE);
    padic_init2(x, 20);
    padic_init2(y, 20);
    padic_init2(z, 20);

    padic_set_ui(x, 15, ctx);            /* x = 3 * 5, valuation 1 */
    ok = padic_exp(y, x, ctx);
    flint_printf("exp(15) converges: %d\ny = ", ok);
    padic_print(y, ctx); flint_printf("\n");

    ok = padic_log(z, y, ctx);
    flint_printf("log(y) converges: %d, log(y) == 15: %d\n",
                 ok, padic_equal(z, x));

    padic_log_balanced(z, y, ctx);
    flint_printf("balanced agrees: %d\n", padic_equal(z, x));
    padic_log_satoh(z, y, ctx);
    flint_printf("Satoh agrees: %d\n", padic_equal(z, x));

    padic_set_ui(x, 3, ctx);             /* a unit, not 1 mod 5 */
    flint_printf("exp(3) converges: %d, log(3) converges: %d\n",
                 padic_exp(y, x, ctx), padic_log(z, x, ctx));

    padic_clear(x); padic_clear(y); padic_clear(z);
    padic_ctx_clear(ctx);
    fmpz_clear(p);
    return 0;
}
```
```output
exp(15) converges: 1
y = 87481552191941
log(y) converges: 1, log(y) == 15: 1
balanced agrees: 1
Satoh agrees: 1
exp(3) converges: 0, log(3) converges: 0
```

The value printed is the representative of $\exp(15)$ modulo $5^{20}$.
Applying the logarithm to it recovers 15 exactly, because
$\log(\exp(x)) = x$ on the domain of convergence and both functions are
computed to the full precision of the output.
