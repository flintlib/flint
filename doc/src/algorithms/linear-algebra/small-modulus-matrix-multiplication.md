---
title: Matrix multiplication modulo small integers
summary: How `nmod_mat_mul_u8` multiplies matrices modulo $N \le 255$ with entries in bytes, replacing multiplications by in-register table lookups for $N \le 16$, and with bit-packed Strassen–Winograd and the method of four Russians for $N = 2$ and $N = 3$.
complexity: For $N \le 16$, one table lookup and one byte addition per product, with a reduction every $\lfloor (256 - N)/(N - 1) \rfloor$ steps; for $N = 2$, about $mkn/(64 t)$ word operations with $t = 8$, plus Strassen–Winograd above dimension 8192.
order: 25
---

Moduli below 256 occur in coding theory, combinatorics and modular
algorithms over small finite fields. Storing such residues in 64-bit
words wastes most of the memory bandwidth, and multiplying them with
64-bit or floating-point instructions wastes most of the arithmetic
width. `nmod_mat_mul_u8` converts the entries to bytes (or bits) and
uses methods that exploit their size. `_nmod_mat_mul_u8` is the same
multiplication on byte matrices, used directly by the `gr` contexts for
`nmod8`. Throughout, $N$ is the modulus and $A$, $B$ have dimensions
$m \times k$ and $k \times n$.

## Byte kernel for $N \le 16$ {#lookup}

If $N \le 16$, every residue is below 16, so it can index a 16-byte
table held in a vector register. Let $\lambda_c$ be the table with
$\lambda_c[j] = cj \bmod N$ for $0 \le j < 16$. For an entry $c$ of $A$
and a vector $b$ of 16 or 32 bytes of a row of $B$, a byte shuffle
instruction (`vpshufb` on x86-64 with AVX2, `tbl` on AArch64) computes
all products $(c \, b_j) \bmod N$ at once. The multiplication is
replaced by a lookup, and the products are already reduced.

The micro-kernel keeps a tile of $5 \times 64$ entries of $C$ in byte
accumulators. Each step adds the looked-up products of one entry of $A$
per row to one row of 64 bytes of $B$.

::: proposition {#prop-cadence}
If each accumulator starts below $N$ and receives at most
$s = \lfloor (255 - (N - 1))/(N - 1) \rfloor$ products before it is
reduced, then it never exceeds 255.
:::

::: proof
Each product from the table is at most $N - 1$, so after $s$ steps the
accumulator is at most $(N - 1) + s (N - 1) \le 255$.
:::

The *cadence* $s$ is 254 for $N = 2$, 41 for $N = 7$, 20 for $N = 13$
and 16 for $N = 16$. An accumulator byte $x = 16 x_1 + x_0$ is reduced
with two more lookups, $y = (16 x_1 \bmod N) + (x_0 \bmod N) < 2N$,
followed by $\min(y, y - N)$ computed with byte wraparound, which is
$y \bmod N$. A product of a byte would contribute up to $(N - 1)^2$ to
the accumulator, so the lookup also makes the cadence $N - 1$ times
longer than with multiplications.

Without a byte shuffle instruction, the portable code multiplies in
16-bit lanes and reduces every
$\lfloor (65535 - (N - 1))/(N - 1)^2 \rfloor$ steps, which is still 291
steps at $N = 16$.

Above $N = 16$ the residues no longer fit a single 16-byte table, and
the lookups would have to be combined.

## Larger moduli {#sgemm}

For $17 \le N \le 255$, `nmod_mat_mul_u8` uses single-precision
floating point. By
[](algorithms/linear-algebra/matrix-multiplication#prop-exact), a block
of $k_{\max} = \lfloor (2^{24} - 1)/(N - 1)^2 \rfloor$ columns of $A$ can be multiplied exactly by `flint_sgemm` with residues in
$[0, N)$; $k_{\max} \ge 260$. The blocks of the inner dimension are converted from bytes to
floats, multiplied, and the exact block results are added in 32-bit
integers, which are reduced modulo $N$ only before they could exceed
$2^{31}$ and at the end. When several threads are available, $m, n \ge 256$
and the whole product fits in one block, a single threaded `sgemm` on
the whole matrices is used instead.

## Strassen–Winograd on bytes {#strassen}

Above the leaves, `nmod_mat_mul_u8` applies
[Strassen–Winograd](algorithms/linear-algebra/strassen-winograd) to the
byte matrices, with byte additions modulo $N$ and temporaries in bytes,
an eighth of the size of `ulong` temporaries. The recursion splits each
dimension $d$ into $\lceil d/2 \rceil$ and $\lfloor d/2 \rfloor$, extends
short operands by zeros in the additions, and stops when a dimension is
below 256 for $N \le 16$ and below 2048 for larger $N$, whose `sgemm`
leaves need large blocks to amortise the conversions. With several
threads, $C$ is divided into a grid of tiles of at least $256 \times 256$
entries, each computed by an independent recursion.

## Bit-packed kernels for $N = 2$ and $N = 3$ {#packed}

For $N = 2$ an entry is one bit, and 64 entries of a row are packed into
a word, so that adding rows is an exclusive or of words. FLINT multiplies
by the *method of four Russians* [@ADKF1970], as in the M4RI library
[@ABH2010]:

::: algorithm Method of four Russians over GF(2) {#alg-m4r}
**Input:** $A$ ($m \times k$) and $B$ ($k \times n$) over GF(2), rows
packed into words; a parameter $t$.

**Output:** $C = AB$.

1. For each group of $t$ consecutive rows $B_l, \ldots, B_{l + t - 1}$
   of $B$:
   1. build the table $T$ of all $2^t$ sums
      $\sum_{i \in S} B_{l + i}$, $S \subseteq \{0, \ldots, t - 1\}$,
      visiting the subsets in Gray code order, so that each entry is the
      previous one plus a single row;
   2. for each row $r$ of $A$, read the $t$ bits
      $A_{r, l}, \ldots, A_{r, l + t - 1}$ as an index $s$ and add
      $T[s]$ to row $r$ of $C$.
:::

Each table costs $2^t$ row additions and saves $t - 1$ of the $t$ row
additions per row of $A$, so with $t = 8$ the algorithm needs about
$mkn / (64 t) + 2^t k n / (64 t)$ word operations. FLINT builds eight
tables at a time and absorbs them in one pass over $C$, so that each row
of $C$ is read and written once per 64 columns of $A$.

For $N = 3$ an entry is represented by two bits in two separate words
(*bit planes*), with $0 = (0, 0)$, $1 = (1, 0)$ and $2 = (0, 1)$, so a
pair of words holds 64 entries [@BB2009]. Addition modulo 3 of two such
pairs takes six bitwise operations, and negation exchanges the two
planes, so subtraction costs the same as addition. The method of four
Russians works as for GF(2), with tables of the $3^t$ combinations
$\sum_{i} c_i B_{l + i}$, $c_i \in \{0, 1, 2\}$, for $t = 5$ (243
entries), two tables at a time, and an index formed from the base-3
digits of the entries of $A$.

Above these leaves, both cases use Strassen–Winograd on the packed
matrices, whose additions cost one (GF(2)) or six (GF(3)) operations per
64 entries. Because the leaves are so cheap, the recursion only starts
when all dimensions are at least 8192 for $N = 2$ and 16384 for $N = 3$,
with at most 6 levels. The entries are packed directly from the `ulong`
matrices, without an intermediate byte image. With several threads, the
packed operands are shared, and $C$ is divided into a grid of tiles,
splitting rows before columns since each extra strip of rows only costs
rebuilding the tables. When $A$ has fewer than 64 ($N = 2$) or 192
($N = 3$) rows, the tables do not pay off and the byte kernel is used
instead.
