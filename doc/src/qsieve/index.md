---
title: qsieve.h
header: qsieve.h
description: The self-initialising multiple polynomial quadratic sieve for factoring integers of 40 to 90 digits.
category: Integers
order: 8
---

This module implements the self-initialising multiple polynomial quadratic
sieve (SIQS) with the large prime variation, the method `fmpz_factor`
uses for the composite cofactors left after trial division and the
elliptic curve method. It is the method of choice for integers without
small factors of about 40 to 100 digits. The sieve runs in parallel on the
threads made available with `flint_set_num_threads`.

The only function meant for general use is `qsieve_factor`; normally one
calls `fmpz_factor`, which also handles primes, perfect powers and small
factors, and finishes the factorisation. The remaining sections
document the stages of the sieve. The algorithm, its parameters and its
complexity are described in [](algorithms/number-theory/quadratic-sieve).

::: summary algorithms/number-theory/quadratic-sieve :::

## Factoring {#factoring}

::: function
void qsieve_factor(fmpz_factor_t factors, const fmpz_t n)

Finds a nontrivial factorisation of $n$, which must be composite and not a
perfect power; the sign of $n$ is multiplied into `factors->sign`. The bases written to `factors` (which
should be freshly initialised) multiply to $|n|$ with the stored
exponents, but they need not be prime: the sieve splits $n$ along the
square roots it finds and stops as soon as all parts are probable primes
or it runs out of square roots. Use `fmpz_factor` for a complete
factorisation into proven primes.

::: details
If $n$ has a small factor that turns up while the factor base is computed
(a prime up to the largest factor base prime), that factor and the
complete factorisation of the cofactor (by `fmpz_factor_no_trial`) are
returned instead. The tuning parameters are read from `qsieve_tune` for
the bit size of $n$.
:::
:::

::: function
void qsieve_factor_with_tune(fmpz_factor_t factors, const fmpz_t n, ulong ks_primes, slong fb_primes, slong small_primes, slong sieve_size, ulong sieve_bits)

`qsieve_factor` with the given tuning parameters instead of those from
the table (see `qsieve_tune` for their meaning). Used for tuning.
:::

## Parameters {#parameters}

::: variable
const ulong qsieve_tune

The table of tuning parameters `qsieve_tune[][6]`, a static array in the
header with `QS_TUNE_SIZE` rows, each of the form `{bits, ks_primes,
fb_primes, small_primes, sieve_size, sieve_bits}`: the bit size of $n$
from which the row applies, the number of primes considered in the
Knuth–Schroeppel function, the number of factor base primes (including
the multiplier and 2), the number of small primes not sieved with, the
size of the sieve interval in bytes, and the sieve threshold (in units of
$\log_2$, after biasing). The rows go from 10 to 300 bits; for example,
for 200-bit numbers (61 digits) the factor base has 4500 primes and the
sieve 131072 entries, and at 300 bits (91 digits) 140000 primes and
2097152 entries.
:::

::: macro
QS_TUNE_SIZE

The number of rows of `qsieve_tune`.
:::

::: macro
BLOCK_SIZE
QS_L1_BITS
QS_L1_SIZE
QS_SMALL_BOUND

The sieve interval is processed in blocks of `BLOCK_SIZE` ($2^{19}$)
bytes, which fit in the L2 cache. Within a block, the primes below
`QS_SMALL_BOUND` (4096), which hit it many times, are sieved over
sub-blocks of `QS_L1_SIZE` ($2^{15}$) bytes, which fit in the L1 cache,
and the primes up to `BLOCK_SIZE` directly over the block. The primes
above `BLOCK_SIZE` hit a block at most a few times and are handled by a
bucket sieve once per polynomial.
:::

::: macro
QS_LP_MULT
QS_TRIAL_CHUNK
QS_HASH_BITS_MIN
QS_HASH_BITS_MAX
BITS_ADJUST
QS_DEBUG

Further parameters: large primes are accepted below `QS_LP_MULT` (60)
times the largest factor base prime; candidates are trial divided by
`QS_TRIAL_CHUNK` (64) primes at a time; the hash table of large primes
has between $2^{12}$ and $2^{20}$ buckets, depending on the size of the
factor base; `BITS_ADJUST` (25) compensates for the approximate logarithms
of the sieve; and `QS_DEBUG` selects the amount of debugging output (0,
none, by default).
:::

## Types {#types}

::: type
qs_s
qs_t

The state of a factorisation: $n$, the multiplier $k$ and $kn$, the factor
base and its data, the sieve parameters, the data of the current
polynomials, the relations found so far with the hash table of their
large primes, and the matrix for the linear algebra. Initialised by
`qsieve_init` and freed by `qsieve_clear`, `qsieve_poly_clear` and
`qsieve_linalg_clear`.
:::

::: type
qs_poly_s
qs_poly_t

The per-thread polynomial data: the coefficients $B$ and $C$ of the
current polynomial, the two sieve start positions for each factor base
prime, the bucket sieve buffers, and space for the factorisation of a
candidate.
:::

::: type
prime_t

A factor base prime `p`, with two precomputed inverses (`pinv`, `pinv2`)
and its size in bits (`size`), used as its approximate logarithm.
:::

::: type
fac_t

A factor of a relation: the index `ind` of a factor base prime and its
exponent `exp`.
:::

::: type
relation_t
qs_rel_s

A relation as returned by `qsieve_get_relation` and as held in the
relation store: its large primes `lp` and `lp2` (both 1 for a full
relation, `lp2` = 1 for a partial relation with one large prime), the
exponents of the small primes, the other factor base primes dividing it
with their exponents, and the value $Y$ whose square is congruent to the
product of its primes modulo $kn$.
:::

::: type
hash_t

An entry of the hash table of large primes: the prime, the index of the
next entry in its chain, its number of occurrences, and its parent in the
union-find forest of the graph of large primes.
:::

::: type
la_col_t

A column of the sparse matrix over $\mathbb{F}_2$: the list `data` of the
`weight` rows in which it is 1, and the index `orig` of the relation it
comes from.
:::

## Initialisation {#init}

::: function
void qsieve_init(qs_t qs_inf, const fmpz_t n)
void qsieve_init_with_tune(qs_t qs_inf, const fmpz_t n, ulong ks_primes, slong fb_primes, slong small_primes, slong sieve_size, ulong sieve_bits)

Initialise `qs_inf` for factoring $n$, with the tuning parameters of the
table row for the bit size of $n$, respectively the given ones. The sieve
values are biased so that an entry is over the threshold exactly when its
top bit is set, which lets the scan of the sieve test a word at a time.
:::

::: function
void qsieve_clear(qs_t qs_inf)

Frees $n$, $kn$, the factor base and the relation store of `qs_inf`.
:::

## Multiplier and factor base {#factor-base}

::: function
ulong qsieve_knuth_schroeppel(qs_t qs_inf)

Chooses the Knuth–Schroeppel multiplier $k$ for $n$, among the 29
squarefree numbers from 1 to 47, maximising the expected contribution of
the small primes to the sieve: $kn \bmod 8$ and the primes $p$ (up to
`ks_primes` of them) for which $kn$ is a square modulo $p$, minus
$\frac12 \log k$ [@Sil1987]. Stores it in `qs_inf->k` and returns 0, or
returns a small prime factor of $n$ if one is found (2 if $n$ is even).
:::

::: function
prime_t * compute_factor_base(ulong * small_factor, qs_t qs_inf, slong num_primes)

Extends the factor base to `num_primes` primes: the multiplier, 2, and
the primes $p$ for which $kn$ is a nonzero square modulo $p$, with their
inverses and the square roots of $kn$ modulo $p$. If a prime dividing $n$
is met, it is stored in `small_factor` and the extension stops. Returns
the factor base.
:::

::: function
ulong qsieve_primes_init(qs_t qs_inf)
ulong qsieve_primes_increment(qs_t qs_inf, ulong delta)

Compute the initial factor base (of `fb_primes` primes, plus `ks_primes`
extra primes held in reserve), respectively extend it by `delta` primes
when sieving has to continue with a larger factor base. Return a small
prime factor of $n$ if one is found, and 0 otherwise.
:::

## Polynomials {#polynomials}

SIQS uses polynomials $Q(x) = (Ax + B)^2 - kn = A (Ax^2 + 2Bx + C)$,
where $A$ is a product of $s$ factor base primes close to its optimal size
$\sqrt{2kn}/M$ for a sieve interval $[-M, M]$, and $B$ runs over the
$2^{s-1}$ solutions of $B^2 \equiv kn \pmod A$ up to sign, which are
enumerated in Gray code order so that each new $B$ and its sieve start
positions cost only additions.

::: function
ulong qsieve_poly_init(qs_t qs_inf)
void qsieve_poly_clear(qs_t qs_inf)

Allocate the polynomial data (for $s$ factors of $A$, the factor base and
one `qs_poly_s` per thread), set up the parameters of the bucket sieve,
and free them.
:::

::: function
int qsieve_init_A(qs_t qs_inf)

Chooses the number $s$ of prime factors of $A$ and the range of the factor
base from which they are taken (primes of about $(\log_2 A)/s$ bits,
following msieve), allocates the polynomial data and computes the first
$A$. Returns 0 if no suitable range exists, in which case the factor base
must be extended.
:::

::: function
int qsieve_next_A(qs_t qs_inf)

Finds the next value of $A$. If $A$ has $s > 3$ factors, the first $s - 1$
are chosen at random among the smallest primes of the allowed range whose
index is not divisible by 4 (the number of candidates being adjusted so
that the products have about the right size), and the last factor is found
by binary search among those with index divisible by 4, so that $A$ is
within a factor of two of its optimal value; values of $A$ used before are
rejected. With at most 3 factors, the next subset in lexicographic order
is used. Returns 0 if no new value is found.
:::

::: function
void qsieve_reinit_A(qs_t qs_inf)
void qsieve_A_used_clear(qs_t qs_inf)

Restart the enumeration of the values of $A$ after a failed linear algebra
step, and free the list of the values of $A$ used so far.
:::

::: function
void qsieve_init_poly_first(qs_t qs_inf)

For the current $A$, computes the data for generating the values of $B$
(for each prime $q_i \mid A$, the term $B_i = (A/q_i) \gamma_i$ with
$\gamma_i \equiv \sqrt{kn} (A/q_i)^{-1} \pmod{q_i}$, and
$2 B_i A^{-1} \bmod p$ for every factor base prime $p$), the first $B$ and
the sieve start positions $(\pm \sqrt{kn} - B) A^{-1} \bmod p$ for every
factor base prime.
:::

::: function
void qsieve_init_poly_next(qs_t qs_inf, slong i)

Moves to the $(i+1)$-th value of $B$ in Gray code order, updating $B$ and
the sieve start positions by adding or subtracting $2 B_\nu$ and its
precomputed multiples.
:::

::: function
void qsieve_compute_C(fmpz_t C, qs_t qs_inf, qs_poly_t poly)

Sets $C = (B^2 - kn)/A$ for the current polynomial (an exact division).
:::

::: function
void qsieve_poly_copy(qs_poly_t poly, qs_t qs_inf)

Copies the current $B$ and sieve start positions to the per-thread data
`poly`.
:::

## Sieving {#sieving}

::: function
void qsieve_do_sieving(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)
void qsieve_do_sieving2(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)
void qsieve_do_sieving3(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)

Sieve the interval for the polynomial `poly`: initialise the sieve and add
the approximate logarithm of $p$ at the positions congruent to either
start position modulo $p$, for every factor base prime $p$ beyond the
small primes. The first version sieves the whole array prime by prime;
the second works block by block, with the primes below `QS_SMALL_BOUND`
sieved over L1-sized sub-blocks; the third (used unless the sieve is
small) adds a bucket sieve for the primes above `BLOCK_SIZE`: each of them
is visited once per polynomial and its hits are appended to a bucket per
sub-block, and the buckets are applied to the sub-blocks in turn and kept
for `qsieve_evaluate_candidate`.
:::

::: function
slong qsieve_evaluate_sieve(qs_t qs_inf, unsigned char * sieve, qs_poly_t poly)

Scans the sieve for entries over the threshold, testing a word at a time,
and calls `qsieve_evaluate_candidate` for each. Returns the number of
relations found.
:::

::: function
slong qsieve_evaluate_candidate(qs_t qs_inf, ulong i, unsigned char * sieve, qs_poly_t poly)

Checks whether $Q(x)$ at sieve position $i$ factors over the factor base,
apart from at most one large prime below the large prime bound. The small
primes are divided out by trial division, the medium primes by checking
the sieve start positions, and the large factor base primes are read off
the buckets. Full and partial relations are stored with
`qsieve_add_relation`. Returns the number of relations found (0 or 1).
:::

::: function
slong qsieve_collect_relations(qs_t qs_inf, unsigned char * sieve)

Sieves with all the polynomials for the current $A$ (distributing the
values of $B$ over the threads) and returns the number of relations
found.
:::

## Relations and large primes {#relations}

The relations are kept in memory. A partial relation with large prime
$\ell$ is an edge between $\ell$ and the vertex 1 of a graph whose
vertices are the large primes; every independent cycle of this graph
gives a full relation, by multiplying the partial relations around it
(their large primes then occur squared).

::: function
void qsieve_relations_init(qs_t qs_inf)
void qsieve_relations_clear(qs_t qs_inf)
void qsieve_relations_reset(qs_t qs_inf)

Initialise the store of relations (called by `qsieve_init_with_tune`),
free it (called by `qsieve_clear`), and discard all relations while
keeping the allocations (used when the linear algebra failed to split $n$
and sieving starts again with a larger factor base).
:::

::: function
void qsieve_add_relation(qs_t qs_inf, ulong prime, ulong prime2, const fmpz_t Y, const qs_poly_t poly)

Appends a relation to the store, with large primes `prime` and `prime2`
(`prime2` = 1 for at most one large prime, both 1 for a full relation),
the square root `Y`, and the small prime exponents and remaining factors
taken from `poly`. The store supports two large primes, but the sieve
currently only produces partial relations with one. Called from the
sieving threads, so the caller must hold `qs_inf->mutex`.
:::

::: function
void qsieve_add_partial_edge(qs_t qs_inf, ulong prime, ulong prime2)

Adds the edge `prime` – `prime2` (with `prime2` = 1 for a single large
prime) to the graph of large primes, incrementing `qs_inf->num_cycles` if
it closes a new independent cycle, as detected by a union-find
structure. The caller must hold `qs_inf->mutex`.
:::

::: function
hash_t * qsieve_get_table_entry(qs_t qs_inf, ulong prime)
void qsieve_add_to_hashtable(qs_t qs_inf, ulong prime)

Return the entry of `prime` in the hash table of large primes, creating
it if necessary, and add an occurrence of `prime` to the table.
:::

::: function
relation_t qsieve_get_relation(qs_t qs_inf, slong i)

Returns a copy of the $i$-th relation of the store, which the caller must
free.
:::

::: function
relation_t qsieve_merge_relation(qs_t qs_inf, relation_t a, relation_t b)

Combines two partial relations with the same large prime $\ell$ into a
full relation (the product, with $\ell^2$ removed from the factorisation
and $Y$ adjusted accordingly).
:::

::: function
int qsieve_compare_relation(const void * a, const void * b)
int qsieve_remove_duplicates(relation_t * rel_list, slong num_relations)

Comparison of relations (by large prime, then number of factors, then the
factors), for sorting, and removal of duplicates from a list of relations
by sorting it; returns the new length.
:::

::: function
void qsieve_insert_relation(qs_t qs_inf, relation_t * rel_list, slong num_relations)

Inserts the given full relations as columns of the matrix, recording
their factors and $Y$ values for the square root step.
:::

::: function
int qsieve_process_relation(qs_t qs_inf)

Processes the relations once enough are expected: removes duplicates,
combines the partial relations along the independent cycles of the graph
of large primes (found from a spanning forest) into full relations, and
if there are at least `extra_rels` more than factor base primes inserts
them in the matrix and returns 1. Otherwise lowers the estimate
`qs_inf->num_cycles` by the shortfall and returns 0, so that sieving
continues. Returns $-1$ if a large prime turned out to divide $n$ (stored
in `qs_inf->small_factor`).
:::

## Linear algebra and square root {#linalg}

::: function
void qsieve_linalg_init(qs_t qs_inf)
void qsieve_linalg_realloc(qs_t qs_inf)
void qsieve_linalg_clear(qs_t qs_inf)

Allocate, resize after the factor base was extended, and free the matrix
and the data for the square root step.
:::

::: function
void reduce_matrix(qs_t qs_inf, slong * nrows, slong * ncols, la_col_t * cols)

Light filtering of the sparse matrix: repeatedly deletes the columns that
contain a row with a single nonzero entry, then trims the matrix to a few
more columns than rows, updating the dimensions.
:::

::: function
uint64_t * block_lanczos(flint_rand_t state, slong nrows, slong dense_rows, slong ncols, la_col_t * B)

Montgomery's block Lanczos algorithm over $\mathbb{F}_2$ with blocks of
64 vectors [@Mon1995]: returns an array of `ncols` words whose bit $l$ in
entry $i$ gives the $i$-th coordinate of the $l$-th vector of a basis of
(part of) the null space of the matrix with columns `B`, or `NULL` if the
iteration failed (it is then retried with another random start). Adapted
from msieve by Jason Papadopoulos.
:::

::: function
uint64_t get_null_entry(uint64_t * nullrows, slong i, slong l)

Returns the entry $i$ of the $l$-th null space vector (nonzero if it is
1).
:::

::: function
void qsieve_square_root(fmpz_t X, fmpz_t Y, qs_t qs_inf, uint64_t * nullrows, slong ncols, slong l, fmpz_t N)

For the $l$-th null space vector, sets $Y$ to the product of the $Y$
values of the relations in it and $X$ to the square root of the product
of their factor base primes (from the halved exponent counts), both modulo
$N = kn$. Then $X^2 \equiv Y^2 \pmod{kn}$, and $\gcd(X - Y, n)$ is a
factor of $n$, nontrivial with probability at least about $1/2$ if $n$ is
not a prime power.
:::

::: function
void insert_col_entry(la_col_t * col, slong entry)
void swap_cols(la_col_t * col2, la_col_t * col1)
void clear_col(la_col_t * col)
void free_col(la_col_t * col)

Inline helpers for matrix columns: append the row `entry`, exchange two
columns, set the weight to zero, and free the column's data.
:::
