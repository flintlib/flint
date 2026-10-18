---
title: Hermite and Smith normal forms
summary: The Hermite and Smith normal forms of integer matrices and the algorithms FLINT implements for them — classical and gcd-based elimination, Kannan–Bachem, Domich–Kannan–Trotter and other computations modulo a determinant, the Pernet–Stein algorithm, Iliopoulos' algorithm, iterated Hermite forms, and Lübeck's $p$-adic method for elementary divisors — with the choices made by fmpz_mat_hnf and fmpz_mat_snf.
complexity: Polynomial in the dimension and the bit size for all algorithms except the classical one; computations modulo $D$ keep entries below $D$ and cost $O(n^3)$ operations on integers of $\log_2 D$ bits for an $n \times n$ matrix.
order: 65
---

The Hermite normal form (HNF) is to integer matrices what the reduced row
echelon form is to matrices over a field: a canonical form under the row
operations that are invertible over $\Z$. The Smith normal form (SNF)
allows column operations as well and describes a matrix up to
equivalence; it gives the structure of finitely generated abelian groups.
Both are computed by elimination with gcds, and the main difficulty is the
growth of the entries, which the algorithms below control in different
ways. There are no generic `gr_mat` versions of these algorithms; they are
specific to $\Z$ (`fmpz_mat`). Over $\Z/N\Z$ the Howell form plays the
role of the HNF ([](algorithms/linear-algebra/gaussian-elimination#howell)).

## Definitions {#definitions}

A matrix $U \in \Z^{m \times m}$ is *unimodular* if $\det U = \pm 1$,
equivalently if $U^{-1}$ is an integer matrix. Multiplying $A$ on the left
by $U$ performs a sequence of elementary row operations over $\Z$:
exchanging rows, negating a row, adding an integer multiple of a row to
another. The rows of $A$ and of $UA$ generate the same lattice
$L(A) = \{x^T A : x \in \Z^m\} \subseteq \Z^n$.

::: definition {#def-hnf}
An $m \times n$ integer matrix $H$ of rank $r$ is in *Hermite normal form*
if its first $r$ rows are nonzero and the others zero, the first nonzero
entry of row $i < r$ is positive and lies in column $j_i$ with
$j_0 < j_1 < \cdots < j_{r-1}$, and $0 \le H_{k j_i} < H_{i j_i}$ for all
$k < i$.
:::

::: theorem {#thm-hnf}
For every integer matrix $A$ there is exactly one $H$ in Hermite normal
form with $H = UA$ for a unimodular $U$. Two matrices of the same shape
have the same Hermite normal form if and only if their rows generate the
same lattice.
:::

::: proof
Existence follows from the algorithms below. For uniqueness, let $H$ and
$H'$ be in Hermite normal form with the same row lattice. The pivot
columns are determined by the lattice ($j_i$ is the least $j$ such that
the lattice vectors vanishing in columns before $j$ have rank $r - i$),
so they agree, and the pivots agree since $H_{i j_i}$ is the positive
generator of the ideal of values in column $j_i$ of lattice vectors
vanishing before $j_i$. If $H \ne H'$, let $i$ be the last row in which
they differ; then $H_i - H'_i$ is a lattice vector vanishing up to column
$j_i$, so it is a combination of rows $i + 1, \ldots$ of $H$; its first
nonzero entry is in some pivot column $j_k$ with $k > i$ and is then a
multiple of $H_{k j_k}$, which contradicts the reduction
$0 \le H_{i j_k}, H'_{i j_k} < H_{k j_k}$.
:::

For a nonsingular $n \times n$ matrix, $H$ is upper triangular and
$\prod_i H_{ii} = |\det A| = [\Z^n : L(A)]$.

::: theorem Smith normal form {#thm-smith}
For every $m \times n$ integer matrix $A$ of rank $r$ there are
unimodular $U$, $V$ with $UAV = \operatorname{diag}(d_1, \ldots, d_r,
0, \ldots, 0)$, $d_i > 0$ and $d_1 \mid d_2 \mid \cdots \mid d_r$. The
$d_i$ are unique: $d_1 d_2 \cdots d_k$ is the gcd $g_k$ of all
$k \times k$ minors of $A$.
:::

::: proof
Existence follows from the algorithms below. The gcd of the $k \times k$
minors is invariant under multiplication by unimodular matrices, since by
the Cauchy–Binet formula the minors of $UA$ are integer combinations of
those of $A$ and conversely. For the diagonal form the gcd is
$d_1 \cdots d_k$ because of the divisibility chain.
:::

The $d_i$ are the *elementary divisors* (or *invariant factors*) of $A$,
and $\Z^n/L(A) \cong \bigoplus_i \Z/d_i\Z \oplus \Z^{n - r}$. The largest
one, $d_r$, is the exponent of the torsion part.

## Hermite normal form algorithms {#hnf}

### Classical elimination {#classical}

`fmpz_mat_hnf_classical` processes the columns from left to right. In the
current column it repeatedly picks the nonzero entry of least absolute
value at or below the current row, moves it into place and reduces the
entries below it with integer division, until they are zero; then it
makes the pivot positive and reduces the entries above it
[@Coh1996, Algorithm 2.4.4]. `fmpz_mat_hnf_xgcd` instead combines two
rows with the cofactors $u, v$ of $g = \gcd(a, b) = ua + vb$ by the
unimodular transformation
$$
\begin{pmatrix} u & v \\ -b/g & a/g \end{pmatrix},
$$
which clears an entry in one step [@Coh1996, Algorithm 2.4.5]. Both are
fast for small matrices, but the entries of the rows not yet processed
can grow exponentially with the dimension [@HM1991].

### Kannan–Bachem {#kannan-bachem}

`fmpz_mat_hnf_minors` [@KanBac1979] brings the leading principal
$k \times k$ submatrix into Hermite normal form for $k = 1, 2, \ldots$:
row $k$ is eliminated against rows $0, \ldots, k - 1$ with the gcd
transformation above, and then the entries above the new diagonal are
reduced modulo the pivots. Because all rows processed so far are kept
reduced, Kannan and Bachem prove that every intermediate entry has a
number of bits polynomial in $n$ and in the size of $A$, which gives the
first polynomial-time Hermite normal form algorithm. FLINT's version
requires $A$ to have full column rank; it brings a zero diagonal entry
into place by exchanging with a later row. `fmpz_mat_hnf_minors_transform`
applies the same operations to an identity matrix to obtain $U$.

### Computing modulo a determinant {#modular}

The entries of an HNF are bounded by its pivots, so it is natural to
compute modulo a multiple of them.

::: lemma {#lem-dzn}
Let $L \subseteq \Z^n$ be a lattice of rank $n$ with elementary divisors
$d_1 \mid \cdots \mid d_n$. Then $d_n \Z^n \subseteq L$, and hence
$D \Z^n \subseteq L$ for every multiple $D$ of $d_n$, in particular for
$D = \det L = d_1 \cdots d_n$.
:::

::: proof
The group $\Z^n / L \cong \bigoplus_i \Z/d_i\Z$ is annihilated by $d_n$.
:::

Hence $L$ is generated by the rows of $A$ together with the rows of
$D I_n$, and all row operations can be done modulo $D$, which bounds the
entries by $D$. `fmpz_mat_hnf_modular` is the algorithm of Domich, Kannan
and Trotter [@DomKanTro1987; @Coh1996, Algorithm 2.4.8]: it eliminates
column $k$ with gcd transformations, reducing modulo the current $R$,
puts $\gcd(H_{kk}, R)$ on the diagonal, and divides $R$ by the pivot for
the next column, since the lattice generated by the remaining columns
has determinant $R/H_{kk}$. `fmpz_mat_hnf_modular_eldiv` uses the
weaker requirement $d_n \mid D$: it computes the strong echelon form
modulo $D$ (`nmod_mat_strong_echelon_form` or
`fmpz_mat_strong_echelon_form_mod`), which by the result of Fieker and
Hofmann [@FieHof2014] is the Hermite normal form of the lattice
$L + D\Z^n = L$ reduced modulo $D$, and replaces zero diagonal entries,
which stand for pivots equal to $D$, by $D$. Hafner and McCurley
[@HM1991] and Storjohann and Labahn [@SL1996] give asymptotically faster
modular algorithms based on matrix multiplication.

### Pernet–Stein {#pernet-stein}

For a nonsingular square matrix, $D = |\det A|$ can be huge, but for
random matrices the lattice generated by $n - 1$ of the rows together
with a combination of the last two has a small determinant. Pernet and
Stein [@PernetStein2010], after Micciancio and Warinschi [@MW2001],
exploit this. For an $m \times n$ matrix of rank $r$,
`fmpz_mat_hnf_pernet_stein` does the following.

1. Find the pivot rows and columns modulo a random prime; let $B$ be the
   first $r - 2$ pivot rows and $c$, $d$ the next two, all restricted to
   the first $r - 1$ pivot columns.
2. *Double determinant*: compute $d_1 = \det \binom{B}{c}$ and
   $d_2 = \det \binom{B}{d}$, both from one solution of a linear system
   over $\Q$ (`fmpq_mat_solve_fmpz_mat`) and a few modular determinants.
   Let $g = \gcd(d_1, d_2) = s d_1 + t d_2$.
3. Replace the pair $c, d$ by $sc + td$; the resulting
   $(r - 1) \times (r - 1)$ matrix $C$ has determinant $\pm g$, which is
   small with high probability for random input. Compute its HNF modulo
   $g$ (`fmpz_mat_hnf_modular`), or recursively if $g$ is large.
4. *Add columns*: the HNF $H_C$ of $C$ is $UC$ with $U = H_C C^{-1}$
   unimodular, so the HNF of the same $r - 1$ rows on all columns is $U$
   times these rows. The other columns are obtained by solving a linear
   system with $C$.
5. *Add rows*: insert the remaining rows of $A$, including the original
   rows $c$ and $d$, one at a time into the HNF by gcd elimination and
   reduction.
6. If the pivot structure obtained differs from the one predicted modulo
   the prime, the prime was bad: restart.

The expected cost for random matrices is dominated by linear system
solving, $O^\sim(n^{\omega + 1} b)$ bit operations. FLINT falls back to
`fmpz_mat_hnf_minors` when $A$ has full column rank and fewer than about
35 to 52 columns, and when the rank is too small for the decomposition.
Faster algorithms that certify unimodularity and handle structured input
better are given in [@PauSto2013].

### Choice and transformation {#choice}

`fmpz_mat_hnf` calls `fmpz_mat_hnf_classical` for matrices with few rows
(fewer than 52 rows for entries of 1 or 2 bits, decreasing to 2 rows for
entries above 512 bits) and `fmpz_mat_hnf_pernet_stein` otherwise.
`fmpz_mat_hnf_transform` computes $U$ with `fmpz_mat_hnf_minors_transform`
for matrices of full column rank, and otherwise as the right block of the
HNF of $(A \mid I_m)$: that matrix has full row rank, so its HNF is
$U(A \mid I_m) = (UA \mid U)$ with $U$ unimodular.

## Smith normal form algorithms {#smith}

### Elimination {#smith-elimination}

`fmpz_mat_snf_kannan_bachem` [@KanBac1979] alternately clears the first
column with row gcd transformations and the first row with column gcd
transformations until both are zero apart from the corner, and repeats if
the corner does not divide the remaining entries (adding a row to the
first row and continuing). Then it proceeds with the submatrix.
`fmpz_mat_snf_iliopoulos` [@Iliopoulos1989] does the same eliminations
modulo a multiple $D$ of the largest elementary divisor of a nonsingular
square matrix, which by [](#lem-dzn) does not change the result and bounds
the entries by $D$; the diagonal is finally normalised by
`fmpz_mat_snf_diagonal`, which replaces pairs of diagonal entries by their
gcd and lcm.

### Iterated Hermite forms {#iterated}

`fmpz_mat_snf` uses a different approach: it alternates the row HNF of
the matrix and the row HNF of its transpose (that is, the column HNF)
until the matrix is diagonal. Each round is a unimodular transformation,
and the leading entries can only become divisors of what they were;
Kannan and Bachem show that the number of rounds is bounded by a
polynomial in the dimension and the size of the entries [@KanBac1979].
The diagonal is then turned into a divisibility chain by gcd and lcm
steps. Since `fmpz_mat_hnf` is fast (Pernet–Stein for
large matrices), this inherits its performance; two or three rounds
suffice for typical matrices. `fmpz_mat_snf_transform` does the same with
`fmpz_mat_hnf_transform` and accumulates $U$ and $V$; the final gcd steps
are done with extended gcds on pairs of rows and columns.

### Elementary divisors by p-adic ranks {#luebeck}

When only the elementary divisors are needed and their prime factors are
known, Lübeck's method [@Lüb2002] computes the valuations from ranks
modulo primes. Let $R$ be an $r \times n$ basis of the row lattice $L$ (the
nonzero rows of the HNF), with elementary divisors $e_1, \ldots, e_r$.

::: proposition {#prop-luebeck}
Let $p$ be a prime, let $y_1, \ldots, y_k$ be a basis of the left kernel
of $R$ modulo $p$ (so $k = r - \operatorname{rank}_p R$), and let $R'$ be
obtained from $R$ by replacing, for each $j$, a row $i_j$ in which $y_j$
has entry 1 (and the other $y_l$ have entry 0) by $y_j R / p$. Then $R'$ is
a basis of $L' = L + \frac1p (L \cap p\Z^n)$, and the elementary divisors
of $L'$ are $e_i/p$ for the $k$ indices with $p \mid e_i$ and $e_i$ for
the others.
:::

::: proof
Choose bases in which $L = \bigoplus_i e_i \Z f_i$ (Smith form). Then
$L \cap p\Z^n = \bigoplus_i \operatorname{lcm}(e_i, p) \Z f_i$ and
$L' = \bigoplus_i e_i' \Z f_i$ with $e_i' = e_i/p$ if $p \mid e_i$ and
$e_i' = e_i$ otherwise; in particular $\operatorname{rank}_p R = \#\{i :
p \nmid e_i\}$, so $k = \#\{i : p \mid e_i\}$. The vectors $y_j R$ span
$L \cap p\Z^n$ modulo $pL$, so the rows of $R'$ generate $L'$, and the
replaced rows are recovered as $p \cdot (y_j R/p)$ minus the other terms
of $y_j R$, so nothing is lost; $R'$ has $r$ rows, so it is a basis.
:::

Iterating until $R$ has full rank modulo $p$, the nullities $k_1 \ge k_2
\ge \cdots$ count the elementary divisors with $v_p(e_i) \ge 1, 2,
\ldots$, which determines all valuations. `fmpz_mat_elementary_divisors`
computes the HNF, factors its pivots (whose product is divisible by every
prime dividing an $e_i$), applies this for each prime, and falls back to
`fmpz_mat_snf` when a pivot is too large to factor or has a prime factor
beyond a word. Faster Smith form algorithms based on linear system
solving are given in [@Gie1995; @EGV2000].

## Use in FLINT {#use}

| Task | Functions |
|---|---|
| HNF, default | `fmpz_mat_hnf` |
| HNF, specific algorithm | `fmpz_mat_hnf_classical`, `fmpz_mat_hnf_xgcd`, `fmpz_mat_hnf_minors`, `fmpz_mat_hnf_modular`, `fmpz_mat_hnf_modular_eldiv`, `fmpz_mat_hnf_pernet_stein` |
| HNF with transformation | `fmpz_mat_hnf_transform`, `fmpz_mat_hnf_minors_transform` |
| SNF | `fmpz_mat_snf`, `fmpz_mat_snf_kannan_bachem`, `fmpz_mat_snf_iliopoulos`, `fmpz_mat_snf_diagonal`, `fmpz_mat_snf_transform` |
| Elementary divisors | `fmpz_mat_elementary_divisors` |
| Checks | `fmpz_mat_is_in_hnf`, `fmpz_mat_is_in_snf` |
| Modulo $N$ | `fmpz_mat_howell_form_mod`, `fmpz_mat_strong_echelon_form_mod`, `nmod_mat_howell_form`, `nmod_mat_strong_echelon_form` |
