---
title: The example programs
summary: What each program in FLINT's examples directory does, how to build and run it, and which ones are good starting points for one's own programs.
order: 3
---

The directory `examples` of the FLINT source tree contains about fifty
complete programs. Some are short demonstrations of a module, others are
command-line tools and benchmarks that were used for record computations.
They are good starting points for one's own programs.

## Building and running {#building}

In a configured FLINT source tree,

```
make examples
```

builds the library and all example programs, placing the executables in
`build/examples`; for example, `build/examples/partitions 1000000000`
computes the number of partitions of $10^9$. `make checkexamples` runs
some of them and checks their output. A single example can also be
compiled against an installed FLINT, since the examples include the
headers as `<flint/...>`:

```
gcc -O2 examples/crt.c -o crt -lflint -lmpfr -lgmp -lm
```

(with `-I`, `-L` and `-Wl,-rpath` options for a nonstandard prefix, as
in [](tutorials/getting-started/first-program)). Most programs print a
usage message when called without arguments, and many accept
`-threads n` to use $n$ threads. Several print their running time and
memory usage with the functions of [](profiler).

## Integers and number theory {#integers}

| Program | What it does | Example |
|---|---|---|
| `crt` | Reconstructs an integer from its residues modulo $2, 3, 5, \ldots$ with `fmpz_CRT_ui`, printing each step | `crt 10382788` |
| `multi_crt` | The same with all residues at once, using `fmpz_multi_mod_ui` and `fmpz_multi_CRT_ui` with a precomputed `fmpz_comb_t` | `multi_crt 100493287498239 13` |
| `factor_integer` | Factors an integer, given as an expression such as `2^128+1`, with `fmpz_factor`; `-certify` proves the primality of the large factors with ECPP and prints the certificates | `factor_integer -timing "2^128+1"` |
| `AKS` | Tests the primality of an integer with the AKS algorithm (`fmpz_is_prime_AKS`) and compares with `fmpz_is_prime` | `AKS 1000003` |
| `primegen` | Prints the primes up to $N$ using the prime iterator of [](ulong_extras), or with `-c` only counts them | `primegen -c 1000000` |
| `check_n_is_prime` | Checks `n_is_prime` against a file of all base-2 Fermat pseudoprimes below $2^{64}$ (not included) | `check_n_is_prime psps.txt` |
| `partitions` | Computes the number of partitions $p(n)$ with `partitions_fmpz_fmpz` | `partitions 1000000000 -threads 4` |
| `bernoulli` | Computes the Bernoulli number $B_n$ exactly | `bernoulli 1000000 -threads 8` |
| `delta_qexp` | Computes the coefficient $\tau(n)$ of $q^n$ in the $q$-expansion of the modular discriminant $\Delta$, with `arith_ramanujan_tau` | `delta_qexp 1000000` |
| `class_poly` | Computes the Hilbert class polynomial of a discriminant | `class_poly -1000004 -threads 8` |
| `mfcoefs` | Computes many Fourier coefficients of a few built-in modular forms (named `11`, `23`, `41`, `131` and `13_4` after their level) from products of Eisenstein series; needs [](fft_small) | `mfcoefs 23 100000` |

For example:

```
> build/examples/crt 10382788
residue mod 2 = 0; reconstruction = 0
residue mod 3 = 1; reconstruction = -2
residue mod 5 = 3; reconstruction = -2
residue mod 7 = 3; reconstruction = -32
residue mod 11 = 9; reconstruction = -662
residue mod 13 = 0; reconstruction = -7592
residue mod 17 = 4; reconstruction = 172588
residue mod 19 = 10; reconstruction = 683098
residue mod 23 = 13; reconstruction = 10382788
> build/examples/factor_integer "2^128+1"
340282366920938463463374607431768211457 =
59649589127497217 * 5704689200685129054721
```

## Polynomials {#polynomials}

| Program | What it does | Example |
|---|---|---|
| `fmpq_poly` | Multiplies two polynomials with rational coefficients; a minimal demonstration of [](fmpq_poly) | `fmpq_poly` |
| `fmpz_poly_q` | Multiplies two rational functions with [](fmpz_poly_q) | `fmpz_poly_q` |
| `fmpz_mod_poly` | Squares a polynomial modulo 7 with [](fmpz_mod_poly) | `fmpz_mod_poly` |
| `fmpz_poly_factor_zassenhaus` | Factors a polynomial of degree 62 over $\Z$ with the Zassenhaus algorithm; a small change makes it read a larger polynomial from the file `fmpz_poly_hensel_P1` | `fmpz_poly_factor_zassenhaus` |
| `factor_polynomial` | Factors a polynomial in $x$, $y$, $z$ over $\Z$, given as an expression, with [](fmpz_mpoly_factor) | `factor_polynomial "(1+x+2*y+3*z)^3+1"` |
| `minimal_irreducibles` | Prints the minimal irreducible polynomials (`nmod_poly_minimal_irreducible`) of degrees $n_{\min}$ to $n_{\max}$ over $\mathbb{F}_p$, in parallel with `flint_parallel_do` | `minimal_irreducibles 2 1 8` |
| `fq_poly` | Times the multiplication of polynomials over three finite fields with the classical, reordering and Kronecker substitution algorithms | `fq_poly` |
| `radix` | Times the radix conversion of a polynomial modulo $17^{26}$ with `fmpz_mod_poly_radix` | `radix` |
| `groebner` | Computes Gröbner bases over $\Z$ or $\Z/p\Z$ with Buchberger's algorithm for built-in or user-given systems, with timings and verification | `groebner -vars x,y -input "x^2+y^2-1,x-y"` |
| `swinnerton_dyer_poly` | Computes the Swinnerton-Dyer polynomial $S_n$ exactly, using exact real numbers | `swinnerton_dyer_poly 10` |
| `poly_roots` | Isolates all complex roots of an integer polynomial (`arb_fmpz_poly_complex_roots`) | `poly_roots -print 15 t 100` |

For example, `minimal_irreducibles 2 1 8` prints the polynomials
`x`, `x^2 + x + 1`, `x^3 + x + 1`, `x^4 + x + 1`, `x^5 + x^2 + 1`,
`x^6 + x + 1`, `x^7 + x + 1` and `x^8 + x^4 + x^3 + x + 1`, and

```
> build/examples/factor_polynomial "(1+x+2*y+3*z)^3+1"
x^3+6*x^2*y+9*x^2*z+3*x^2+12*x*y^2+36*x*y*z+12*x*y+27*x*z^2+18*x*z+3*x+8*y^3+36*y^2*z+12*y^2+54*y*z^2+36*y*z+6*y+27*z^3+27*z^2+9*z+2 =
1*(x+2*y+3*z+2)^1*(x^2+4*x*y+6*x*z+x+4*y^2+12*y*z+2*y+9*z^2+3*z+1)^1
```

## Matrices {#matrices}

| Program | What it does | Example |
|---|---|---|
| `stirling_matrix` | Builds the $n \times n$ matrices of Stirling numbers of the first and second kind and multiplies them, which gives the identity matrix | `stirling_matrix 10` |
| `hilbert_matrix` | Computes the determinant (or with `-eig` the eigenvalues) of the $n \times n$ Hilbert matrix in ball arithmetic, doubling the precision until the result is accurate | `hilbert_matrix 200` |
| `hilbert_matrix_ca` | Computes the eigenvalues of the Hilbert matrix as exact algebraic numbers and verifies the trace and determinant | `hilbert_matrix_ca 5` |

## $p$-adic numbers {#p-adic}

| Program | What it does | Example |
|---|---|---|
| `padic` | Demonstrates conversion and printing of $p$-adic numbers in the different print modes, and the logarithm, with [](padic) | `padic` |
| `qadic` | Demonstrates powers, sums, Teichmüller lifts, inverses, Frobenius images and square roots in unramified extensions with [](qadic) | `qadic` |

## Real and complex numbers {#numerics}

These programs use ball arithmetic ([](arb), [](acb) and the modules
built on them) and are described in [](tutorials/numerics).

| Program | What it does | Example |
|---|---|---|
| `pi` | Computes $\pi$ (or $e$, $\gamma$ and other constants with `-constant`) to $n$ digits | `pi 1000000 -threads 4` |
| `pi_agm` | Computes $\pi$ with the Brent–Salamin AGM iteration | `pi_agm 100000` |
| `logistic` | Iterates the logistic map $x_{k+1} = r x_k (1 - x_k)$, restarting with doubled precision until the result is accurate | `logistic 10000` |
| `real_roots` | Isolates and refines the real roots of a function on an interval with [](arb_calc) | `real_roots 0 0.0 50.0 -verbose` |
| `zeta_zeros` | Computes consecutive zeros of the Riemann zeta function on the critical line | `zeta_zeros -n 1 -count 10` |
| `keiper_li` | Computes the Keiper–Li coefficients $\lambda_0, \ldots, \lambda_n$ | `keiper_li 1000 -threads 2` |
| `lvalue` | Evaluates Dirichlet $L$-functions and their derivatives | `lvalue -character 4 3 -re 2` |
| `lcentral` | Computes the central values $L(1/2, \chi)$ for all characters $\chi$ modulo $q$ in a range | `lcentral 1 8` |
| `complex_plot` | Plots a function of a complex variable by domain colouring to a PPM image | `complex_plot zeta -range -10 10 0 40` |
| `integrals` | Computes a list of difficult integrals with rigorous error bounds using `acb_calc_integrate` | `integrals -i 4` |
| `integrals_double_exp` | Computes integrals with the double exponential (tanh-sinh) quadrature | `integrals_double_exp -p 333` |
| `taylor_integrals` | Computes integrals (an elliptic integral and others) with the Taylor series method of [](acb_calc) | `taylor_integrals 50` |
| `fpwrap` | Demonstrates the floating-point wrappers of [](arb_fpwrap) | `fpwrap` |
| `functions_benchmark` | Times standard elementary and special functions at several precisions | `functions_benchmark` |

## Exact real and complex numbers {#exact}

These programs use [](ca) and [](qqbar) and are described in
[](tutorials/exact): `elementary` and `machin` in
[](tutorials/exact/identities), `binet` (with `swinnerton_dyer_poly`) in
[](tutorials/exact/number-fields), `huge_expr` in
[](tutorials/exact/expressions), and `dft` (with `hilbert_matrix_ca`) in
[](tutorials/exact/linear-algebra).

| Program | What it does | Example |
|---|---|---|
| `elementary` | Evaluates and simplifies elementary expressions exactly, and decides whether some are zero | `elementary` |
| `binet` | Computes Fibonacci numbers with Binet's formula in exact arithmetic | `binet 250` |
| `machin` | Verifies variations of Machin's formula for $\pi$ | `machin` |
| `huge_expr` | Proves the equality of two algebraic numbers given by large expressions with nested square roots | `huge_expr` |
| `dft` | Verifies the identity $\mathrm{DFT}^{-1}(\mathrm{DFT}(x)) = x$ exactly for vectors of exact real numbers such as $x_k = \sqrt{k + 2}$ or $\log(k + 2)$ | `dft 8` |
