---
title: Bell numbers
functions: arith_bell_number arith_bell_number_vec arith_bell_number_nmod arith_bell_number_nmod_vec arith_bell_number_size
---

The Bell number $B_n$ is the number of partitions of a set of $n$
elements: $B_0, B_1, \ldots = 1, 1, 2, 5, 15, 52, 203, \ldots$ In terms
of Stirling numbers of the second kind, $B_n = \sum_{k=0}^n S(n,k)$, and
the exponential generating function is $\sum_{n \ge 0} B_n x^n/n! =
\exp(e^x - 1)$ [@Comtet1974]. The numbers grow slightly slower than $n!$:
$\log B_n = n(\log n - \log\log n - 1 + o(1))$, 6403 bits for
$n = 1000$.

The module computes an isolated $B_n$ (`arith_bell_number`), the vector
$B_0, \ldots, B_{n-1}$ (`arith_bell_number_vec`), and the same modulo a
word-size integer (`arith_bell_number_nmod`,
`arith_bell_number_nmod_vec`), which is much cheaper when only the
residue is needed.

::: summary algorithms/number-theory/stirling-bell :::

## Isolated values {#single}

`arith_bell_number` reads $B_n$ from a table for $n \le 25$ (64-bit),
runs the Bell triangle in fixed-size arrays of up to three words for $n
\le 58$, uses Dobiński's formula
$$
B_n = \frac{1}{e} \sum_{k=0}^{\infty} \frac{k^n}{k!}
$$
for $n \le 3400$, in an exact form without floating-point arithmetic
(`arith_bell_number_dobinski`), and the multimodular algorithm above
(`arith_bell_number_multi_mod`). The latter computes $B_n$ modulo primes
$p > n$ of $\bits - 3$ bits with the identity
$$
B_n = \sum_{k=0}^n \frac{(n-k)^n}{(n-k)!} \sum_{j=0}^k
\frac{(-1)^j}{j!},
$$
which needs $O(n)$ multiplications modulo $p$ per prime, and combines the
residues by Chinese remaindering; the primes are processed in parallel by
the threads of `flint_set_num_threads`. About $n \log_2 n / \bits$ primes
are needed, so the residues cost $O(n^2 \log n)$ bit operations.

## Vectors {#vectors}

The Bell triangle (Aitken's array) computes $B_0, \ldots, B_{n-1}$ with
$n^2/2$ additions: each row starts with the last entry of the previous
row, and each further entry is the sum of its left neighbour and the
entry above it; the rows start with the Bell numbers. Over $\Z$ this
costs $O(n^3 \log n)$ bit operations and is used for $n < 5000$; for
longer vectors `arith_bell_number_vec_multi_mod` computes the vector
modulo word-size primes and reconstructs each entry with only as many
primes as it needs.

Modulo $m$, `arith_bell_number_nmod_vec` chooses between three methods:
the triangle ($O(n^2)$ operations, any $m$, used for $n < 300$), the
power series $\exp(e^x - 1)$ ($O(M(n))$ operations, requires $(n-1)!$ to
be invertible modulo $m$, so in particular $m \ge n$ for prime $m$), and
the ordinary generating function
$$
\sum_{n \ge 0} B_n x^n = \sum_{k \ge 0} \frac{x^k}{(1 - x)(1 - 2x)
\cdots (1 - kx)}
$$
summed by binary splitting ($O(M(n) \log n)$ operations, any $m$).

## Example {#example}

The program prints the first Bell numbers, the size of $B_{1000}$, and
checks `arith_bell_number_nmod` against the residue of the exact value
for the prime $m = 10^9 + 7$ (where the series method applies) and for
the modulus $m = 2^{10} = 1024$ (where it does not).

```c example name=arith_bell_basic
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/nmod.h>
#include <flint/arith.h>

int main(void)
{
    fmpz * b;
    fmpz_t B;
    nmod_t mod;
    slong i;

    b = _fmpz_vec_init(10);
    arith_bell_number_vec(b, 10);
    flint_printf("B_0..B_9:");
    for (i = 0; i < 10; i++)
    {
        flint_printf(" ");
        fmpz_print(b + i);
    }
    flint_printf("\n");

    fmpz_init(B);
    arith_bell_number(B, 1000);
    flint_printf("B_1000 has %wu bits\n", fmpz_bits(B));

    nmod_init(&mod, 1000000007);
    flint_printf("mod 10^9+7: %wu %wu\n",
        arith_bell_number_nmod(1000, mod), fmpz_fdiv_ui(B, mod.n));
    nmod_init(&mod, 1024);
    flint_printf("mod 1024:   %wu %wu\n",
        arith_bell_number_nmod(1000, mod), fmpz_fdiv_ui(B, mod.n));

    fmpz_clear(B);
    _fmpz_vec_clear(b, 10);
    return 0;
}
```
```output
B_0..B_9: 1 1 2 5 15 52 203 877 4140 21147
B_1000 has 6403 bits
mod 10^9+7: 465231251 465231251
mod 1024:   443 443
```
