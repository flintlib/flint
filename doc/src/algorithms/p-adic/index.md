---
title: p-adic numbers
short: p-adic
summary: Arithmetic with p-adic numbers and unramified extensions: precision models, Newton and Hensel lifting, exponentials and logarithms, Frobenius and Teichmüller lifts, norms and traces.
order: 10
---

These pages describe how the modules `padic`, `padic_poly`, `padic_mat`,
`qadic` and `padic_radix` compute in the field $\mathbf{Q}_p$ of
$p$-adic numbers and its unramified extensions $\mathbf{Q}_q$. The
overview compares the modules and lists the method used for each
operation. The page on precision states how much precision each
operation preserves and how the modules represent elements; the other
pages give the algorithms with proofs: Hensel lifting for inverses,
square roots and Teichmüller lifts, the exponential and logarithm, and
the arithmetic of unramified extensions including the Frobenius
automorphism, the norm and the trace. Standard introductions to
$p$-adic numbers are [@Gou1997] and [@Rob2000]; [@Car2017] surveys the
computational models.
