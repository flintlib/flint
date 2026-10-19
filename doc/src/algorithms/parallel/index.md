---
title: Parallel computation
short: Parallelism
summary: How FLINT distributes work over threads: the thread pool, parallel loops and parallel binary splitting.
order: 16
---

FLINT uses a single process-wide pool of threads, created on demand. The
sections below describe how the pool hands out threads, how the generic
parallel helpers divide work among them, and how read-only data is shared
between threads without locks.
