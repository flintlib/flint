---
title: Bernoulli numbers
summary: Computes exact Bernoulli numbers from the zeta function (numerically, with the denominator from von Staudt–Clausen), from Harvey's multimodular algorithm, by a reverse zeta iterator with a table of powers for ranges, and by Ramanujan's recursion or power series modulo primes for short vectors.
complexity: $O(n^2 \log^{O(1)} n)$ bit operations for an isolated $B_n$ by either the zeta or the multimodular method; $O(n^3)$ with a small constant for all of $B_0, \ldots, B_n$ with the reverse iterator.
order: 65
---

The Bernoulli numbers $B_n$ are the rational numbers defined by
$x/(e^x - 1) = \sum_{n \ge 0} B_n x^n/n!$. They vanish for odd $n \ge 3$,
alternate in sign for even $n$, and grow quickly:

$$
B_{2m} = (-1)^{m+1} \frac{2\, (2m)!}{(2\pi)^{2m}} \zeta(2m), \qquad
\log_2 |B_n| = n \log_2 \frac{n}{2\pi e} + O(\log n).
$$

The numerator of $B_n$ thus has about $n \log_2 n$ bits; that of $B_{10^6}$
has about 15.8 million. For the classical theory see [@IR1990, Chapter
15] and [@Coh2007b, Chapter 9].

## The denominator {#denominator}

::: theorem von Staudt–Clausen {#thm-vsc}
For even $n \ge 2$,
$$
B_n + \sum_{(p - 1) \mid n} \frac{1}{p} \in \Z,
$$
the sum running over the primes $p$ with $(p - 1) \mid n$. In particular
the denominator of $B_n$ in lowest terms is $\prod_{(p-1) \mid n} p$,
squarefree and divisible by 6.
:::

The denominator is therefore found from the divisors of $n$
(`arith_bernoulli_number_denom` tests the primes up to $n + 1$, with a
table for $n \le 178$). It is small: at most $\prod_{d \mid n} (d + 1)$,
so it has at most $\sigma_0(n) \log_2(n+1)$ bits, and usually far
fewer. All methods below therefore compute the *numerator* $N_n = d_n B_n$,
an integer.

## The zeta function method {#zeta}

Since $d_n$ is known, it suffices to compute $B_n$ numerically with an
absolute error below $1/(2 d_n)$, that is, to about
$\log_2 |B_n| + \log_2 d_n + 1$ significant bits. With $p$ bits of
precision, the formula above needs $n!$ and $\pi^n$ (cheap: binary
splitting and a power) and $\zeta(n)$, which is $1 + O(2^{-n})$ and for
which very few terms are needed:

$$
\zeta(n) = \sum_{j=1}^{J} j^{-n} + \epsilon, \qquad 0 < \epsilon <
\int_J^\infty t^{-n}\, dt < J^{1-n},
$$

so $J \approx 2^{p/(n-1)}$ terms suffice, which is $O(n/\log n)$ for $p
\approx n \log_2 n$. Equivalently one evaluates the Euler product
$\zeta(n)^{-1} = \prod_{q \le J} (1 - q^{-n})$ over primes, with
$O(n/\log^2 n)$ factors. FLINT's `arb_bernoulli_ui_zeta` uses the Euler
product when $n$ is small compared with the precision and an asymptotic
expansion of $\zeta(n)$ otherwise, in ball arithmetic. `_bernoulli_fmpq_ui_zeta`
multiplies the ball by $d_n$ and returns the unique integer it contains
(increasing the precision in the unlikely event that the ball contains
two). Each power $q^{-n}$ costs $O(M(n \log n) \log n)$, so the total is
$O(n M(n \log n))$, about $n^2 \log^{2+o(1)} n$ bit operations.

### The reverse iterator {#reverse}

To compute all of $B_0, B_2, \ldots, B_n$, evaluating each from scratch
would cost $n^3$. Bloemen's observation [@Blo2009] is that the power
table can be shared: if the powers $j^{-n}$ are stored in fixed point
with $p$ fractional bits, then $j^{-(n-2)} = j^{-n} \cdot j^2$ is one
multiplication by a word. Going *down* in $n$ thus updates the table
with cheap exact multiplications (going up would need divisions), and
the precision and the number of terms needed decrease along the way.

::: algorithm Reverse iteration {#alg-rev}
**Input:** an even $n_{\max}$.

**Output:** $B_{n_{\max}}, B_{n_{\max} - 2}, \ldots, B_0$.

1. Choose the precision $p$ from the size of $B_{n_{\max}}$ and of its
   denominator, and $J \approx 2^{(p+1)/(n_{\max} - 1)}$.
2. For odd $j = 3, 5, \ldots, J$ compute $P_j = \lfloor 2^p j^{-n_{\max}}
   \rfloor$ with a bound on the rounding error. Compute $F = 2\, n_{\max}!
   / (2\pi)^{n_{\max}}$.
3. For $n = n_{\max}, n_{\max} - 2, \ldots$:
   1. $S = 2^{-p} \sum_j P_j$ encloses $\sum_{j \ge 3 \text{ odd}} j^{-n}$
      up to the accumulated rounding errors and the tail of the series.
   2. $\delta = \zeta(n) - 1 = S + (1 + S)/(2^n - 1)$; $B_n = \pm F (1 +
      \delta)$; round $d_n B_n$ to the integer, if the error bound allows
      it.
   3. Update $P_j \gets j^2 P_j$, $F \gets F (2\pi)^2 / (n(n-1))$.
   4. Now and then, lower $p$ and $J$ to the values needed for the
      current $n$, truncating the $P_j$.
:::

Each step costs $O(J)$ multiplications of a $p$-bit number by a word,
$O(n/\log n \cdot n \log n) = O(n^2)$ bit operations, plus two
multiplications of $p$-bit numbers, so all the even Bernoulli numbers up
to $n$ cost $O(n^3)$ in the worst case but with a very small constant.

FLINT's implementation (`bernoulli_rev_init`, `bernoulli_rev_next`)
organizes this as follows.

* **Batches.** The numbers are produced in batches of $k = 16$
  consecutive indices ($k = 8$ below $n = 4000$) and yielded from a
  buffer. Step 3.1 for the whole batch is one pass over the table: for each
  $j$, the $k$ multiplications by $j^2$ and the $k$ additions into the sums
  are done while the number $P_j$ is in cache. The frame $p$ and the number
  of terms $J$ are lowered at every batch.
* **The table.** With AVX-512 IFMA, the $P_j$ are stored as digits of
  50 bits (redundant, below $2^{51}$) in 64-bit lanes, multiplied by $j^2$
  with `vpmadd52luq`/`vpmadd52huq` and summed without carries, the digit
  sums being converted to limbs every 4096 powers; this needs $j <
  2^{24}$. Otherwise the table is in limbs, updated by `mpn_mul_1` and
  summed by `mpn_add_n` directly into the accumulators, eight powers at a
  time by blocks of 256 limbs (the carries of the products between blocks
  being kept), so that a block of a power is stepped $k$ times while in
  cache and the blocks of the $k$ sums serve all eight powers; streaming
  each power and the $k$ sums $k$ times over is bound by the memory
  bandwidth once the powers exceed the L2 cache (half the speed for powers
  of $2 \cdot 10^5$ limbs on one thread). In step 2 the
  powers of the primes are computed by powering in ball arithmetic and
  those of the composites $j = ab$ ($a$ the smallest prime factor) as
  $\lfloor P_a P_b / 2^p \rfloor$ from the top limbs of the entries of $a$
  and $b$, so that nothing beyond the table is held; the primes, and then
  the composites in $[L, 3L)$ for $L = 9, 27, \ldots$ (whose cofactors are
  below $L$), are done in parallel. The table takes about
  $\sum_{j \le J} (p - n \log_2 j) \approx n^2/24$ bits; it shrinks as $n$
  decreases (the frame is lowered by about $2 \log_2(n/2\pi)$ bits per
  step while $P_j$ grows by $2 \log_2 j$ bits, $j \le J \approx n/(2\pi e)$),
  and a buffer is reallocated when it is more than about 1/8 empty.
* **The sieve.** For $n \ge 10^5$, the table only holds the $j$ prime to
  $3$ and $5$ and, apart, the $\{3,5\}$-smooth $j$ (3, 5, 9, 15, 25, ...,
  a few dozen): by unique factorization
  $1 + \sum_{j \ge 3 \text{ odd}} j^{-n} = (1 + G)(1 + S)$ with $G$ and
  $S$ the sums over the smooth and the prime-to-15 $j \le J$, up to terms
  $j^{-n}$ with distinct odd $j > J$, which the tail bound covers. This
  costs one product $GS$ per number, of the size of $15^{-n}$ in the
  frame, and removes $1 - (2/3)(4/5) = 7/15$ of the $j$ and about as much
  of the table (more of the small $j$, which hold the most limbs, are
  removed than of the large ones), saving more in the table pass than the
  product costs at these sizes. The smooth powers stay in limbs (with
  digits, their frame is lowered by shifts) and are summed by one task of
  the pass.
* **The products.** $F$ is a floating-point number of $W \approx p/64 + 2$
  limbs. The two products of a step, $F (2\pi)^2$ and $F \delta$, both
  only need the high half: $\delta$ is aligned so that both read the
  window of limbs $[W - 1, 2W)$ of a product of $W$-limb operands, and with
  `fft_small` the forward transform of $F$ is shared by the two products
  and that of $(2\pi)^2$ is cached, so that a step costs four transforms.
  $W$ is lowered in steps of about 3%, so that the transform plan serves
  many batches. The division by $2^n - 1$ in step 3.2 is done in linear
  time, by summing the $n$-bit chunks of the dividend from the top.
* **Rounding.** All the error terms are tracked, and the fraction of $d_n
  F (1 + \delta)$ must stay away from $1/2$ by more than the error bound;
  otherwise (which is not expected to happen) $B_n$ is computed by
  `_bernoulli_fmpq_ui`.
* **Threads.** The numbers of a batch become independent if the
  prefactors are obtained as $F_{n_0 - 2t} = F_{n_0} (2\pi)^{2t} /
  (n_0 (n_0 - 1) \cdots (n_0 - 2t + 1))$ from the transform of $F_{n_0}$
  and cached transforms of $(2\pi)^{2t}$, at the cost of a division by a
  short multi-limb integer per number. With threads, the numbers of a
  batch are tasks taken by the threads (the expensive ones first), and the
  threads then take chunks of powers for the table pass of the next batch,
  summing into their own accumulators, which the tasks of the next batch
  add up. The threads are requested when a batch starts, and the serial
  scheme is used if none are free.

Threading inside each number (for example in the transforms) does not
pay at these sizes, while separate iterators on disjoint ranges of $n$,
which parallelise perfectly, each need their own table and initialisation
and, to yield the numbers in order, the buffering of a whole range.
`bernoulli_fmpq_vec_no_cache` (which fills the cache used by the `arb`
special functions) stores the whole range anyway: with up to three
threads it runs one iterator, and with more it runs independent iterators
on 32 blocks of the range (largest blocks first), which cost about 1.45
times as much in total but parallelise perfectly.

## Harvey's multimodular algorithm {#multimodular}

Harvey [@Har2010] computes $B_n \bmod p$ for a prime $p$ in time
$O(p \log^{1+o(1)} p)$, independently of $n$:

1. If $(p - 1) \mid n$, $p$ divides the denominator. Otherwise, by
   **Kummer's congruence** $B_n/n \equiv B_m/m \pmod p$ for $m = n \bmod
   (p-1)$, so it suffices to compute $B_m/m \bmod p$ with $2 \le m \le p
   - 3$.
2. **Voronoi's congruence** expresses $(c^m - 1) B_m / m$, for any $c$
   coprime to $p$, as a sum over $x \in (\Z/p\Z)^\times$ of $x^{m-1}$
   times a simple function of $\lfloor cx/p \rfloor$. Enumerating $x$ as
   the powers of a primitive root $g$, and using the symmetry $x
   \leftrightarrow p - x$, gives a sum of $(p-1)/2$ terms, each one
   multiplication modulo $p$.
3. With $c = 2$ (when $2^m \not\equiv 1 \pmod p$), the floor terms are
   the bits of the binary expansion of $x/p$. Harvey reorganises the sum
   so that the inner loop reads the expansion of $1/p$ word by word and
   accumulates into tables indexed by 8-bit chunks; the tables are
   combined with precomputed weights at the end. This makes the constant
   very small. If $2^m \equiv 1$, the slower generic sum with $c = g$ is
   used.

With all primes $5 \le p \le X$, $(p-1) \nmid n$, up to $X \approx n
\log n$, the product exceeds $|N_n|$, and the Chinese remainder theorem
gives the numerator: $O(n^2 \log^{O(1)} n)$ operations in total, but
fully parallel and with smaller constants than the zeta method for large
$n$. FLINT combines the two methods: the residues modulo a product $M$ of
primes give the low part $r = N_n \bmod M$, and the zeta method only
needs to find the integer $t = (N_n - r)/M$, a number with a fraction $1
- \alpha$ of the bits of $N_n$. The default $\alpha$ grows from 0 at $n =
18000$ to $0.28$ for large $n$.

Asymptotically faster methods exist: Harvey's later algorithm computes
$B_n$ in $O(n^{4/3 + o(1)})$ bit operations [@HarveyBern2014], and
Brent and Harvey compute all of $B_0, \ldots, B_n$ in quasi-optimal time
with power series [@BrentHarvey2013]; they are not implemented.

## Short vectors {#vectors}

For a vector $B_0, \ldots, B_{n-1}$ with moderate $n$, `arith` uses:

- **Ramanujan's recursion** (`_arith_bernoulli_number_vec_recursive`,
  $n < 700$). Ramanujan's identity [@Ram1911]
  $$
  \binom{m+3}{m} B_m = \begin{cases} \frac{m+3}{3} \\ -\frac{m+3}{6}
  \end{cases} - \sum_{j=1}^{\lfloor m/6 \rfloor} \binom{m+3}{m - 6j}
  B_{m-6j}, \quad m \equiv 0, 2 \text{ or } 4 \pmod 6,
  $$
  (with $(m+3)/3$ for $m \equiv 0, 2$ and $-(m+3)/6$ for $m \equiv 4$)
  needs only a sixth of the previous values. FLINT runs it on integers
  over the common denominator $(n + 1)\#$, so that no GCD is needed until
  the final division by the true denominators.
- **Power series modulo primes**
  (`_arith_bernoulli_number_vec_multi_mod`, $n \ge 700$). The even
  Bernoulli numbers are the coefficients of
  $$
  \frac{x^2}{\cosh x - 1} = \sum_{k \ge 0} \frac{(2 - 4k)\, B_{2k}}{(2k)!}
  x^{2k},
  $$
  obtained by one power series inversion of $(\cosh x - 1)/x^2$, which has
  only the even terms [@BuhlerCrandallSompolski1992]. Multiplying by the
  known denominators and factorials modulo each prime gives the
  numerators modulo the prime, and Chinese remaindering the numerators.

## Euler numbers {#euler}

The Euler numbers $E_n$, $1/\cosh t = \sum E_n t^n/n!$, are integers with
$|E_{2m}| = 2^{2m+2} (2m)!\, \beta(2m+1)/\pi^{2m+1}$, where $\beta(s) =
\sum_{k \ge 0} (-1)^k (2k+1)^{-s}$ is the Dirichlet beta function. FLINT
computes an isolated $E_n$ from this formula for $n < 2000$ and with a
multimodular method combined with the beta function above
(`arb_fmpz_euler_number_ui`), and a vector by inverting the power series
of $\cosh$ modulo primes (`arith_euler_number_vec`).

## Use in FLINT {#use}

| Function | Method |
|---|---|
| `bernoulli_fmpq_ui` | the cache if possible; zeta for $n < 18000$; zeta with multimodular low part above |
| `_bernoulli_fmpq_ui_zeta` | the zeta function via `arb_bernoulli_ui_zeta` |
| `_bernoulli_fmpq_ui_multi_mod` | `bernoulli_mod_p_harvey` on many primes in parallel, plus the zeta function for the high part |
| `bernoulli_mod_p_harvey` | Harvey's algorithm with Kummer's congruence |
| `bernoulli_rev_init`, `bernoulli_rev_next` | the reverse iterator [](#alg-rev), in batches, threaded |
| `bernoulli_fmpq_vec_no_cache`, `bernoulli_cache_compute` | one reverse iterator, or iterators on blocks in parallel (more than three threads) |
| `arith_bernoulli_number_vec` | Ramanujan's recursion ($n < 700$) or power series modulo primes |
| `arith_euler_number`, `arith_euler_number_vec` | the beta function and multimodular methods; power series modulo primes |

The generic functions call these: `gr_bernoulli_ui` uses
`bernoulli_fmpq_ui` over $\Q$ (and in rings without a real precision,
followed by a conversion), `gr_bernoulli_vec` uses
`bernoulli_fmpq_vec_no_cache` over $\Q$, and
`gr_eulernum_vec` uses `arith_euler_number_vec` over $\Z$; in real and
complex rings they compute numerically with `arb`.

::: functions bernoulli_fmpq_ui bernoulli_mod_p_harvey bernoulli_fmpq_vec_no_cache arith_bernoulli_number_vec :::
