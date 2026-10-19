---
title: Parallel loops and binary splitting
summary: How flint_parallel_do schedules independent calls over threads, why strided scheduling balances monotone costs to within one call, and how parallel binary splitting divides threads over a recursion tree.
order: 2
---

Two generic helpers distribute work over the threads of the pool
([](algorithms/parallel/thread-pool)): `flint_parallel_do` runs
independent calls $f(0), \ldots, f(n-1)$, and
`flint_parallel_binary_splitting` runs the two halves of a binary
splitting recursion in parallel.

## Static scheduling of a loop {#static}

`flint_parallel_do` evaluates $f(0), \ldots, f(n-1)$ with $T$ threads,
where $T$ is the number of workers obtained plus the calling thread. The
indices are divided once, before the threads start; there is no work
stealing. Two schedules are available.

**Blocks** (`FLINT_PARALLEL_UNIFORM`): with $s = \lceil n/T \rceil$, thread
$t$ gets the indices $ts \le i < \min((t+1)s, n)$. When all calls cost the
same, the loads differ by at most $s$ calls, and each thread
works on a contiguous range, which is good for memory locality.

**Strides** (`FLINT_PARALLEL_STRIDED`): thread $t$ gets the indices
$i \equiv t \pmod T$. This is the better choice when the cost $c_i$ of
call $i$ varies monotonically, as it does when $i$ indexes coefficients or
subproblems of increasing size.

::: theorem Balance of strided scheduling {#thm-strided}
Let $c_0 \le c_1 \le \cdots \le c_{n-1}$ and let
$S_t = \sum_{i \equiv t \bmod T} c_i$ be the load of thread $t$ under
strided scheduling. Then $\lvert S_t - S_u \rvert \le c_{n-1}$ for all
$0 \le t, u < T$.
:::

::: proof
Let $t < u$ and let $m_t$ and $m_u$ be the numbers of indices of the two
threads. The index sets are $\{t + kT : 0 \le k < m_t\}$ and
$\{u + kT : 0 \le k < m_u\}$, and $m_t - 1 \le m_u \le m_t$ because
$t < u < t + T$.

Since $t + kT < u + kT$, monotonicity gives $c_{t+kT} \le c_{u+kT}$ for
$k < m_u$. Hence $S_t \le S_u + (m_t - m_u) c_{n-1} \le S_u + c_{n-1}$.

Since $u + kT < t + (k+1)T$, also $c_{u+kT} \le c_{t+(k+1)T}$ whenever
$k + 1 < m_t$. If $m_u = m_t - 1$, summing over $k < m_u$ gives
$S_u \le S_t - c_t \le S_t$. If $m_u = m_t$, the same sum covers all but the
last index of thread $u$, so $S_u \le S_t - c_t + c_{u + (m_u - 1)T} \le
S_t + c_{n-1}$. The case of nonincreasing costs is symmetric.
:::

With blocks, monotone costs can make the imbalance as large as
$\lceil n/T \rceil (c_{n-1} - c_0)$, roughly the difference between the
first and last block.

The flag `FLINT_PARALLEL_DYNAMIC` is reserved for dynamic scheduling, in
which threads take new indices from a shared counter as they finish. It is
not implemented and currently gives blocks.

## Parallel binary splitting {#bsplit}

Binary splitting evaluates $P(a, b) = f(a) \circ \cdots \circ f(b-1)$ as
$P(a, m) \circ P(m, b)$ with $m = a + \lfloor (b - a)/2 \rfloor$, recursing
until $b - a$ is at most a cutoff [@HP1998; @Ber2008]. The two recursive
calls are independent, so the recursion tree has natural parallelism near
its root.

::: algorithm flint_parallel_binary_splitting {#alg-bsplit}
**Input:** a range $[a, b)$, a cutoff $c \ge 1$ and a thread budget $t$
(`thread_limit`, or `flint_get_num_threads()` if not positive).

1. If $b - a \le c$, compute $P(a, b)$ with the base case and return.
2. Request one worker, subject to the caller's own budget.
3. If one was obtained: wake it on $P(m, b)$ with budget
   $\lfloor t/2 \rfloor - 1$, lower the caller's budget to
   $t - \lfloor t/2 \rfloor - 1$ and compute $P(a, m)$, then restore the
   budget and wait for the worker.
4. Otherwise compute $P(a, m)$ and $P(m, b)$ one after the other.
5. Return the merge of the two results.
:::

Each level halves the budget, so with $T = 2^k$ threads the
top $k$ levels of the tree run in parallel and every thread then works on
one subtree of size about $(b - a)/T$. The two recursive calls receive
their shares as `thread_limit`, so no more than $t$ threads run base
cases at the same time. For products whose operands grow along the
recursion, most of the time is spent in the last merges near the root, so
the speedup of binary splitting by itself is limited; it is most effective
when the merge operation, typically a large multiplication, is itself
parallel.

::: functions flint_parallel_do flint_parallel_binary_splitting :::
