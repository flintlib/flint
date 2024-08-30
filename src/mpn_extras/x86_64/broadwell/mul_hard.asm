dnl
dnl Copyright (C) 2023 Albin Ahlbäck
dnl
dnl This file is part of FLINT.
dnl
dnl FLINT is free software: you can redistribute it and/or modify it under
dnl the terms of the GNU Lesser General Public License (LGPL) as published
dnl by the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.  See <https://www.gnu.org/licenses/>.
dnl

include(`config.m4')

dnl TODO: A lot to fix here...
dnl * Instead of flint_mpn_mul_M_N for hardcoded M and N, do flint_mpn_mul_M_n,
dnl   where n is a variable instead. This will reduce the amount of code, and
dnl   probably be around the same speed, although one register has to go to n
dnl   (%rcx).
dnl * Fix latencies.
dnl * Minimize 32-bit instructions by registers.
dnl * Investigate partially storing in stack(?) instead of storing in registers.

	TEXT

.macro	m3_str res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5
	mulx	(0+\ap_offset)*8(\ap), \r0, \r3
	mulx	(1+\ap_offset)*8(\ap), \r1, \r4
	mulx	(2+\ap_offset)*8(\ap), \r2, \r5
	add	\r3, \r1
	adc	\r4, \r2
	mov	\r0, (0+\res_offset)*8(\res)
	mov	\r1, (1+\res_offset)*8(\res)
	mov	\r2, (2+\res_offset)*8(\res)
	adc	$0, \r5
.endm

.macro	m4_str res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5
	mulx	(0+\ap_offset)*8(\ap), \r0, \r3
	mulx	(1+\ap_offset)*8(\ap), \r1, \r5
	mulx	(2+\ap_offset)*8(\ap), \r2, \r4
	add	\r3, \r1
	adc	\r5, \r2
	mov	\r0, (0+\res_offset)*8(\res)
	mov	\r1, (1+\res_offset)*8(\res)
	mov	\r2, (2+\res_offset)*8(\res)
	mulx	(3+\ap_offset)*8(\ap), \r3, \r5
	adc	\r3, \r4
	adc	$0, \r5
	mov	\r4, (3+\res_offset)*8(\res)
.endm

.macro	m5_str res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5
	mulx	(0+\ap_offset)*8(\ap), \r0, \r3
	mulx	(1+\ap_offset)*8(\ap), \r5, \r4
	mulx	(2+\ap_offset)*8(\ap), \r2, \r1
	add	\r3, \r5
	adc	\r4, \r2
	mov	\r0, (0+\res_offset)*8(\res)
	mov	\r5, (1+\res_offset)*8(\res)
	mov	\r2, (2+\res_offset)*8(\res)
	mulx	(3+\ap_offset)*8(\ap), \r3, \r4
	mulx	(4+\ap_offset)*8(\ap), \r0, \r5
	adc	\r3, \r1
	adc	\r4, \r0
	mov	\r1, (3+\res_offset)*8(\res)
	mov	\r0, (4+\res_offset)*8(\res)
	adc	$0, \r5
.endm

.macro	m6_str res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5
	mulx	(0+\ap_offset)*8(\ap), \r0, \r3
	mulx	(1+\ap_offset)*8(\ap), \r1, \r5
	mulx	(2+\ap_offset)*8(\ap), \r2, \r4
	add	\r3, \r1
	adc	\r5, \r2
	mov	\r0, (0+\res_offset)*8(\res)
	mov	\r1, (1+\res_offset)*8(\res)
	mov	\r2, (2+\res_offset)*8(\res)
	mulx	(3+\ap_offset)*8(\ap), \r3, \r5
	mulx	(4+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r4
	adc	\r5, \r0
	mov	\r4, (3+\res_offset)*8(\res)
	mov	\r0, (4+\res_offset)*8(\res)
	mulx	(5+\ap_offset)*8(\ap), \r3, \r5
	adc	\r3, \r1
	adc	$0, \r5
	mov	\r1, (5+\res_offset)*8(\res)
.endm

.macro	m7_str res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5
	mulx	(0+\ap_offset)*8(\ap), \r0, \r3
	mulx	(1+\ap_offset)*8(\ap), \r1, \r4
	mulx	(2+\ap_offset)*8(\ap), \r2, \r5
	add	\r3, \r1
	adc	\r4, \r2
	mov	\r0, (0+\res_offset)*8(\res)
	mov	\r1, (1+\res_offset)*8(\res)
	mov	\r2, (2+\res_offset)*8(\res)
	mulx	(3+\ap_offset)*8(\ap), \r3, \r4
	mulx	(4+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r5
	adc	\r4, \r0
	mov	\r5, (3+\res_offset)*8(\res)
	mov	\r0, (4+\res_offset)*8(\res)
	mulx	(5+\ap_offset)*8(\ap), \r3, \r4
	mulx	(6+\ap_offset)*8(\ap), \r2, \r5
	adc	\r3, \r1
	adc	\r4, \r2
	mov	\r1, (5+\res_offset)*8(\res)
	mov	\r2, (6+\res_offset)*8(\res)
	adc	$0, \r5
.endm

.macro	m8_str res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5
	mulx	(0+\ap_offset)*8(\ap), \r0, \r3
	mulx	(1+\ap_offset)*8(\ap), \r1, \r5
	mulx	(2+\ap_offset)*8(\ap), \r2, \r4
	add	\r3, \r1
	adc	\r5, \r2
	mov	\r0, (0+\res_offset)*8(\res)
	mov	\r1, (1+\res_offset)*8(\res)
	mov	\r2, (2+\res_offset)*8(\res)
	mulx	(3+\ap_offset)*8(\ap), \r3, \r5
	mulx	(4+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r4
	adc	\r5, \r0
	mov	\r4, (3+\res_offset)*8(\res)
	mov	\r0, (4+\res_offset)*8(\res)
	mulx	(5+\ap_offset)*8(\ap), \r3, \r5
	mulx	(6+\ap_offset)*8(\ap), \r2, \r4
	adc	\r3, \r1
	adc	\r5, \r2
	mov	\r1, (5+\res_offset)*8(\res)
	mov	\r2, (6+\res_offset)*8(\res)
	mulx	(7+\ap_offset)*8(\ap), \r3, \r5
	adc	\r3, \r4
	adc	$0, \r5
	mov	\r4, (7+\res_offset)*8(\res)
.endm

.macro	m9_str res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5
	mulx	(0+\ap_offset)*8(\ap), \r0, \r3
	mulx	(1+\ap_offset)*8(\ap), \r5, \r4
	mulx	(2+\ap_offset)*8(\ap), \r2, \r1
	add	\r3, \r5
	adc	\r4, \r2
	mov	\r0, (0+\res_offset)*8(\res)
	mov	\r5, (1+\res_offset)*8(\res)
	mov	\r2, (2+\res_offset)*8(\res)
	mulx	(3+\ap_offset)*8(\ap), \r3, \r4
	mulx	(4+\ap_offset)*8(\ap), \r0, \r5
	adc	\r3, \r1
	adc	\r4, \r0
	mov	\r1, (3+\res_offset)*8(\res)
	mov	\r0, (4+\res_offset)*8(\res)
	mulx	(5+\ap_offset)*8(\ap), \r3, \r4
	mulx	(6+\ap_offset)*8(\ap), \r2, \r1
	adc	\r3, \r5
	adc	\r4, \r2
	mov	\r5, (5+\res_offset)*8(\res)
	mov	\r2, (6+\res_offset)*8(\res)
	mulx	(7+\ap_offset)*8(\ap), \r3, \r4
	mulx	(8+\ap_offset)*8(\ap), \r0, \r5
	adc	\r3, \r1
	adc	\r4, \r0
	mov	\r1, (7+\res_offset)*8(\res)
	mov	\r0, (8+\res_offset)*8(\res)
	adc	$0, \r5
.endm

.macro	m10_str res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5
	mulx	(0+\ap_offset)*8(\ap), \r0, \r3
	mulx	(1+\ap_offset)*8(\ap), \r1, \r5
	mulx	(2+\ap_offset)*8(\ap), \r2, \r4
	add	\r3, \r1
	adc	\r5, \r2
	mov	\r0, (0+\res_offset)*8(\res)
	mov	\r1, (1+\res_offset)*8(\res)
	mov	\r2, (2+\res_offset)*8(\res)
	mulx	(3+\ap_offset)*8(\ap), \r3, \r5
	mulx	(4+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r4
	adc	\r5, \r0
	mov	\r4, (3+\res_offset)*8(\res)
	mov	\r0, (4+\res_offset)*8(\res)
	mulx	(5+\ap_offset)*8(\ap), \r3, \r5
	mulx	(6+\ap_offset)*8(\ap), \r2, \r4
	adc	\r3, \r1
	adc	\r5, \r2
	mov	\r1, (5+\res_offset)*8(\res)
	mov	\r2, (6+\res_offset)*8(\res)
	mulx	(7+\ap_offset)*8(\ap), \r3, \r5
	mulx	(8+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r4
	adc	\r5, \r0
	mov	\r4, (7+\res_offset)*8(\res)
	mov	\r0, (8+\res_offset)*8(\res)
	mulx	(9+\ap_offset)*8(\ap), \r3, \r5
	adc	\r3, \r1
	adc	$0, \r5
	mov	\r1, (9+\res_offset)*8(\res)
.endm

.macro	m11_str res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5
	mulx	(0+\ap_offset)*8(\ap), \r0, \r3
	mulx	(1+\ap_offset)*8(\ap), \r1, \r4
	mulx	(2+\ap_offset)*8(\ap), \r2, \r5
	add	\r3, \r1
	adc	\r4, \r2
	mov	\r0, (0+\res_offset)*8(\res)
	mov	\r1, (1+\res_offset)*8(\res)
	mov	\r2, (2+\res_offset)*8(\res)
	mulx	(3+\ap_offset)*8(\ap), \r3, \r4
	mulx	(4+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r5
	adc	\r4, \r0
	mov	\r5, (3+\res_offset)*8(\res)
	mov	\r0, (4+\res_offset)*8(\res)
	mulx	(5+\ap_offset)*8(\ap), \r3, \r4
	mulx	(6+\ap_offset)*8(\ap), \r2, \r5
	adc	\r3, \r1
	adc	\r4, \r2
	mov	\r1, (5+\res_offset)*8(\res)
	mov	\r2, (6+\res_offset)*8(\res)
	mulx	(7+\ap_offset)*8(\ap), \r3, \r4
	mulx	(8+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r5
	adc	\r4, \r0
	mov	\r5, (7+\res_offset)*8(\res)
	mov	\r0, (8+\res_offset)*8(\res)
	mulx	(9+\ap_offset)*8(\ap), \r3, \r4
	mulx	(10+\ap_offset)*8(\ap), \r2, \r5
	adc	\r3, \r1
	adc	\r4, \r2
	mov	\r1, (9+\res_offset)*8(\res)
	mov	\r2, (10+\res_offset)*8(\res)
	adc	$0, \r5
.endm

.macro	m12_str res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5
	mulx	(0+\ap_offset)*8(\ap), \r0, \r3
	mulx	(1+\ap_offset)*8(\ap), \r1, \r5
	mulx	(2+\ap_offset)*8(\ap), \r2, \r4
	add	\r3, \r1
	adc	\r5, \r2
	mov	\r0, (0+\res_offset)*8(\res)
	mov	\r1, (1+\res_offset)*8(\res)
	mov	\r2, (2+\res_offset)*8(\res)
	mulx	(3+\ap_offset)*8(\ap), \r3, \r5
	mulx	(4+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r4
	adc	\r5, \r0
	mov	\r4, (3+\res_offset)*8(\res)
	mov	\r0, (4+\res_offset)*8(\res)
	mulx	(5+\ap_offset)*8(\ap), \r3, \r5
	mulx	(6+\ap_offset)*8(\ap), \r2, \r4
	adc	\r3, \r1
	adc	\r5, \r2
	mov	\r1, (5+\res_offset)*8(\res)
	mov	\r2, (6+\res_offset)*8(\res)
	mulx	(7+\ap_offset)*8(\ap), \r3, \r5
	mulx	(8+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r4
	adc	\r5, \r0
	mov	\r4, (7+\res_offset)*8(\res)
	mov	\r0, (8+\res_offset)*8(\res)
	mulx	(9+\ap_offset)*8(\ap), \r3, \r5
	mulx	(10+\ap_offset)*8(\ap), \r2, \r4
	adc	\r3, \r1
	adc	\r5, \r2
	mov	\r1, (9+\res_offset)*8(\res)
	mov	\r2, (10+\res_offset)*8(\res)
	mulx	(11+\ap_offset)*8(\ap), \r3, \r5
	adc	\r3, \r4
	adc	$0, \r5
	mov	\r4, (11+\res_offset)*8(\res)
.endm

.macro	m13_str res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5
	mulx	(0+\ap_offset)*8(\ap), \r0, \r3
	mulx	(1+\ap_offset)*8(\ap), \r5, \r4
	mulx	(2+\ap_offset)*8(\ap), \r2, \r1
	add	\r3, \r5
	adc	\r4, \r2
	mov	\r0, (0+\res_offset)*8(\res)
	mov	\r5, (1+\res_offset)*8(\res)
	mov	\r2, (2+\res_offset)*8(\res)
	mulx	(3+\ap_offset)*8(\ap), \r3, \r4
	mulx	(4+\ap_offset)*8(\ap), \r0, \r5
	adc	\r3, \r1
	adc	\r4, \r0
	mov	\r1, (3+\res_offset)*8(\res)
	mov	\r0, (4+\res_offset)*8(\res)
	mulx	(5+\ap_offset)*8(\ap), \r3, \r4
	mulx	(6+\ap_offset)*8(\ap), \r2, \r1
	adc	\r3, \r5
	adc	\r4, \r2
	mov	\r5, (5+\res_offset)*8(\res)
	mov	\r2, (6+\res_offset)*8(\res)
	mulx	(7+\ap_offset)*8(\ap), \r3, \r4
	mulx	(8+\ap_offset)*8(\ap), \r0, \r5
	adc	\r3, \r1
	adc	\r4, \r0
	mov	\r1, (7+\res_offset)*8(\res)
	mov	\r0, (8+\res_offset)*8(\res)
	mulx	(9+\ap_offset)*8(\ap), \r3, \r4
	mulx	(10+\ap_offset)*8(\ap), \r2, \r1
	adc	\r3, \r5
	adc	\r4, \r2
	mov	\r5, (9+\res_offset)*8(\res)
	mov	\r2, (10+\res_offset)*8(\res)
	mulx	(11+\ap_offset)*8(\ap), \r3, \r4
	mulx	(12+\ap_offset)*8(\ap), \r0, \r5
	adc	\r3, \r1
	adc	\r4, \r0
	mov	\r1, (11+\res_offset)*8(\res)
	mov	\r0, (12+\res_offset)*8(\res)
	adc	$0, \r5
.endm

.macro	m14_str res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5
	mulx	(0+\ap_offset)*8(\ap), \r0, \r3
	mulx	(1+\ap_offset)*8(\ap), \r1, \r5
	mulx	(2+\ap_offset)*8(\ap), \r2, \r4
	add	\r3, \r1
	adc	\r5, \r2
	mov	\r0, (0+\res_offset)*8(\res)
	mov	\r1, (1+\res_offset)*8(\res)
	mov	\r2, (2+\res_offset)*8(\res)
	mulx	(3+\ap_offset)*8(\ap), \r3, \r5
	mulx	(4+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r4
	adc	\r5, \r0
	mov	\r4, (3+\res_offset)*8(\res)
	mov	\r0, (4+\res_offset)*8(\res)
	mulx	(5+\ap_offset)*8(\ap), \r3, \r5
	mulx	(6+\ap_offset)*8(\ap), \r2, \r4
	adc	\r3, \r1
	adc	\r5, \r2
	mov	\r1, (5+\res_offset)*8(\res)
	mov	\r2, (6+\res_offset)*8(\res)
	mulx	(7+\ap_offset)*8(\ap), \r3, \r5
	mulx	(8+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r4
	adc	\r5, \r0
	mov	\r4, (7+\res_offset)*8(\res)
	mov	\r0, (8+\res_offset)*8(\res)
	mulx	(9+\ap_offset)*8(\ap), \r3, \r5
	mulx	(10+\ap_offset)*8(\ap), \r2, \r4
	adc	\r3, \r1
	adc	\r5, \r2
	mov	\r1, (9+\res_offset)*8(\res)
	mov	\r2, (10+\res_offset)*8(\res)
	mulx	(11+\ap_offset)*8(\ap), \r3, \r5
	mulx	(12+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r4
	adc	\r5, \r0
	mov	\r4, (11+\res_offset)*8(\res)
	mov	\r0, (12+\res_offset)*8(\res)
	mulx	(13+\ap_offset)*8(\ap), \r3, \r5
	adc	\r3, \r1
	adc	$0, \r5
	mov	\r1, (13+\res_offset)*8(\res)
.endm

.macro	m15_str res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5
	mulx	(0+\ap_offset)*8(\ap), \r0, \r3
	mulx	(1+\ap_offset)*8(\ap), \r1, \r4
	mulx	(2+\ap_offset)*8(\ap), \r2, \r5
	add	\r3, \r1
	adc	\r4, \r2
	mov	\r0, (0+\res_offset)*8(\res)
	mov	\r1, (1+\res_offset)*8(\res)
	mov	\r2, (2+\res_offset)*8(\res)
	mulx	(3+\ap_offset)*8(\ap), \r3, \r4
	mulx	(4+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r5
	adc	\r4, \r0
	mov	\r5, (3+\res_offset)*8(\res)
	mov	\r0, (4+\res_offset)*8(\res)
	mulx	(5+\ap_offset)*8(\ap), \r3, \r4
	mulx	(6+\ap_offset)*8(\ap), \r2, \r5
	adc	\r3, \r1
	adc	\r4, \r2
	mov	\r1, (5+\res_offset)*8(\res)
	mov	\r2, (6+\res_offset)*8(\res)
	mulx	(7+\ap_offset)*8(\ap), \r3, \r4
	mulx	(8+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r5
	adc	\r4, \r0
	mov	\r5, (7+\res_offset)*8(\res)
	mov	\r0, (8+\res_offset)*8(\res)
	mulx	(9+\ap_offset)*8(\ap), \r3, \r4
	mulx	(10+\ap_offset)*8(\ap), \r2, \r5
	adc	\r3, \r1
	adc	\r4, \r2
	mov	\r1, (9+\res_offset)*8(\res)
	mov	\r2, (10+\res_offset)*8(\res)
	mulx	(11+\ap_offset)*8(\ap), \r3, \r4
	mulx	(12+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r5
	adc	\r4, \r0
	mov	\r5, (11+\res_offset)*8(\res)
	mov	\r0, (12+\res_offset)*8(\res)
	mulx	(13+\ap_offset)*8(\ap), \r3, \r4
	mulx	(14+\ap_offset)*8(\ap), \r2, \r5
	adc	\r3, \r1
	adc	\r4, \r2
	mov	\r1, (13+\res_offset)*8(\res)
	mov	\r2, (14+\res_offset)*8(\res)
	adc	$0, \r5
.endm

.macro	m16_str res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5
	mulx	(0+\ap_offset)*8(\ap), \r0, \r3
	mulx	(1+\ap_offset)*8(\ap), \r1, \r5
	mulx	(2+\ap_offset)*8(\ap), \r2, \r4
	add	\r3, \r1
	adc	\r5, \r2
	mov	\r0, (0+\res_offset)*8(\res)
	mov	\r1, (1+\res_offset)*8(\res)
	mov	\r2, (2+\res_offset)*8(\res)
	mulx	(3+\ap_offset)*8(\ap), \r3, \r5
	mulx	(4+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r4
	adc	\r5, \r0
	mov	\r4, (3+\res_offset)*8(\res)
	mov	\r0, (4+\res_offset)*8(\res)
	mulx	(5+\ap_offset)*8(\ap), \r3, \r5
	mulx	(6+\ap_offset)*8(\ap), \r2, \r4
	adc	\r3, \r1
	adc	\r5, \r2
	mov	\r1, (5+\res_offset)*8(\res)
	mov	\r2, (6+\res_offset)*8(\res)
	mulx	(7+\ap_offset)*8(\ap), \r3, \r5
	mulx	(8+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r4
	adc	\r5, \r0
	mov	\r4, (7+\res_offset)*8(\res)
	mov	\r0, (8+\res_offset)*8(\res)
	mulx	(9+\ap_offset)*8(\ap), \r3, \r5
	mulx	(10+\ap_offset)*8(\ap), \r2, \r4
	adc	\r3, \r1
	adc	\r5, \r2
	mov	\r1, (9+\res_offset)*8(\res)
	mov	\r2, (10+\res_offset)*8(\res)
	mulx	(11+\ap_offset)*8(\ap), \r3, \r5
	mulx	(12+\ap_offset)*8(\ap), \r0, \r1
	adc	\r3, \r4
	adc	\r5, \r0
	mov	\r4, (11+\res_offset)*8(\res)
	mov	\r0, (12+\res_offset)*8(\res)
	mulx	(13+\ap_offset)*8(\ap), \r3, \r5
	mulx	(14+\ap_offset)*8(\ap), \r2, \r4
	adc	\r3, \r1
	adc	\r5, \r2
	mov	\r1, (13+\res_offset)*8(\res)
	mov	\r2, (14+\res_offset)*8(\res)
	mulx	(15+\ap_offset)*8(\ap), \r3, \r5
	adc	\r3, \r4
	adc	$0, \r5
	mov	\r4, (15+\res_offset)*8(\res)
.endm

.macro	m3 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	mov	\scr1, \res_offset*8(\res)
	adcx	\scr2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	adcx	\scr1, \r1
	adcx	\zero, \r2
.endm

.macro	m4 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	mov	\scr1, \res_offset*8(\res)
	adcx	\scr2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	mulx	(3+\ap_offset)*8(\ap), \scr2, \r3
	adcx	\scr1, \r1
	adcx	\scr2, \r2
	adcx	\zero, \r3
.endm

.macro	m5 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	mov	\scr1, \res_offset*8(\res)
	adcx	\scr2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	mulx	(3+\ap_offset)*8(\ap), \scr2, \r3
	adcx	\scr1, \r1
	adcx	\scr2, \r2
	mulx	(4+\ap_offset)*8(\ap), \scr1, \r4
	adcx	\scr1, \r3
	adcx	\zero, \r4
.endm

.macro	m6 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	mov	\scr1, \res_offset*8(\res)
	adcx	\scr2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	mulx	(3+\ap_offset)*8(\ap), \scr2, \r3
	adcx	\scr1, \r1
	adcx	\scr2, \r2
	mulx	(4+\ap_offset)*8(\ap), \scr1, \r4
	mulx	(5+\ap_offset)*8(\ap), \scr2, \r5
	adcx	\scr1, \r3
	adcx	\scr2, \r4
	adcx	\zero, \r5
.endm

.macro	m7 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5, r6, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	mov	\scr1, \res_offset*8(\res)
	adcx	\scr2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	mulx	(3+\ap_offset)*8(\ap), \scr2, \r3
	adcx	\scr1, \r1
	adcx	\scr2, \r2
	mulx	(4+\ap_offset)*8(\ap), \scr1, \r4
	mulx	(5+\ap_offset)*8(\ap), \scr2, \r5
	adcx	\scr1, \r3
	adcx	\scr2, \r4
	mulx	(6+\ap_offset)*8(\ap), \scr1, \r6
	adcx	\scr1, \r5
	adcx	\zero, \r6
.endm

.macro	m8 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5, r6, r7, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	mov	\scr1, \res_offset*8(\res)
	adcx	\scr2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	mulx	(3+\ap_offset)*8(\ap), \scr2, \r3
	adcx	\scr1, \r1
	adcx	\scr2, \r2
	mulx	(4+\ap_offset)*8(\ap), \scr1, \r4
	mulx	(5+\ap_offset)*8(\ap), \scr2, \r5
	adcx	\scr1, \r3
	adcx	\scr2, \r4
	mulx	(6+\ap_offset)*8(\ap), \scr1, \r6
	mulx	(7+\ap_offset)*8(\ap), \scr2, \r7
	adcx	\scr1, \r5
	adcx	\scr2, \r6
	adcx	\zero, \r7
.endm

.macro	m3_chain res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, ip1, ip2, r0, r1, r2, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	adcx	\ip1, \scr1
	mov	\scr1, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	adcx	\scr2, \r0
	adox	\ip2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	adcx	\scr1, \r1
	adcx	\zero, \r2
.endm

.macro	m4_chain res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, ip1, ip2, r0, r1, r2, r3, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	adcx	\ip1, \scr1
	mov	\scr1, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	adcx	\scr2, \r0
	adox	\ip2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	mulx	(3+\ap_offset)*8(\ap), \scr2, \r3
	adcx	\scr1, \r1
	adcx	\scr2, \r2
	adcx	\zero, \r3
.endm

.macro	m2_chain res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, ip1, ip2, r0, r1, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	adcx	\ip1, \scr1
	mov	\scr1, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	adcx	\scr2, \r0
	adox	\ip2, \r0
	adcx	\zero, \r1
.endm

.macro	m5_chain res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, ip1, ip2, r0, r1, r2, r3, r4, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	adcx	\ip1, \scr1
	mov	\scr1, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	adcx	\scr2, \r0
	adox	\ip2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	mulx	(3+\ap_offset)*8(\ap), \scr2, \r3
	adcx	\scr1, \r1
	adcx	\scr2, \r2
	mulx	(4+\ap_offset)*8(\ap), \scr1, \r4
	adcx	\scr1, \r3
	adcx	\zero, \r4
.endm

.macro	m6_chain res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, ip1, ip2, r0, r1, r2, r3, r4, r5, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	adcx	\ip1, \scr1
	mov	\scr1, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	adcx	\scr2, \r0
	adox	\ip2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	mulx	(3+\ap_offset)*8(\ap), \scr2, \r3
	adcx	\scr1, \r1
	adcx	\scr2, \r2
	mulx	(4+\ap_offset)*8(\ap), \scr1, \r4
	mulx	(5+\ap_offset)*8(\ap), \scr2, \r5
	adcx	\scr1, \r3
	adcx	\scr2, \r4
	adcx	\zero, \r5
.endm

.macro	am2 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, scr, zero
	mulx	(0+\ap_offset)*8(\ap), \r2, \scr
	adcx	\r2, \r0
	mov	\r0, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \r0, \r2
	adcx	\scr, \r1
	adox	\r0, \r1
	adcx	\zero, \r2
	adox	\zero, \r2
.endm

.macro	am3 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, scr, zero
	mulx	(0+\ap_offset)*8(\ap), \scr, \r3
	adcx	\scr, \r0
	mov	\r0, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r3, \r1
	adox	\r0, \r1
	mulx	(2+\ap_offset)*8(\ap), \r0, \r3
	adcx	\scr, \r2
	adox	\r0, \r2
	adcx	\zero, \r3
	adox	\zero, \r3
.endm

.macro	am4 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, scr, zero
	mulx	(0+\ap_offset)*8(\ap), \r4, \scr
	adcx	\r4, \r0
	mov	\r0, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \r0, \r4
	adcx	\scr, \r1
	adox	\r0, \r1
	mulx	(2+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r4, \r2
	adox	\r0, \r2
	mulx	(3+\ap_offset)*8(\ap), \r0, \r4
	adcx	\scr, \r3
	adox	\r0, \r3
	adcx	\zero, \r4
	adox	\zero, \r4
.endm

.macro	am5 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5, scr, zero
	mulx	(0+\ap_offset)*8(\ap), \scr, \r5
	adcx	\scr, \r0
	mov	\r0, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r5, \r1
	adox	\r0, \r1
	mulx	(2+\ap_offset)*8(\ap), \r0, \r5
	adcx	\scr, \r2
	adox	\r0, \r2
	mulx	(3+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r5, \r3
	adox	\r0, \r3
	mulx	(4+\ap_offset)*8(\ap), \r0, \r5
	adcx	\scr, \r4
	adox	\r0, \r4
	adcx	\zero, \r5
	adox	\zero, \r5
.endm

.macro	am6 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5, r6, scr, zero
	mulx	(0+\ap_offset)*8(\ap), \r6, \scr
	adcx	\r6, \r0
	mov	\r0, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \r0, \r6
	adcx	\scr, \r1
	adox	\r0, \r1
	mulx	(2+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r6, \r2
	adox	\r0, \r2
	mulx	(3+\ap_offset)*8(\ap), \r0, \r6
	adcx	\scr, \r3
	adox	\r0, \r3
	mulx	(4+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r6, \r4
	adox	\r0, \r4
	mulx	(5+\ap_offset)*8(\ap), \r0, \r6
	adcx	\scr, \r5
	adox	\r0, \r5
	adcx	\zero, \r6
	adox	\zero, \r6
.endm

.macro	am7 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5, r6, r7, scr, zero
	mulx	(0+\ap_offset)*8(\ap), \scr, \r7
	adcx	\scr, \r0
	mov	\r0, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r7, \r1
	adox	\r0, \r1
	mulx	(2+\ap_offset)*8(\ap), \r0, \r7
	adcx	\scr, \r2
	adox	\r0, \r2
	mulx	(3+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r7, \r3
	adox	\r0, \r3
	mulx	(4+\ap_offset)*8(\ap), \r0, \r7
	adcx	\scr, \r4
	adox	\r0, \r4
	mulx	(5+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r7, \r5
	adox	\r0, \r5
	mulx	(6+\ap_offset)*8(\ap), \r0, \r7
	adcx	\scr, \r6
	adox	\r0, \r6
	adcx	\zero, \r7
	adox	\zero, \r7
.endm

.macro	am8 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5, r6, r7, r8, scr, zero
	mulx	(0+\ap_offset)*8(\ap), \r8, \scr
	adcx	\r8, \r0
	mov	\r0, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \r0, \r8
	adcx	\scr, \r1
	adox	\r0, \r1
	mulx	(2+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r8, \r2
	adox	\r0, \r2
	mulx	(3+\ap_offset)*8(\ap), \r0, \r8
	adcx	\scr, \r3
	adox	\r0, \r3
	mulx	(4+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r8, \r4
	adox	\r0, \r4
	mulx	(5+\ap_offset)*8(\ap), \r0, \r8
	adcx	\scr, \r5
	adox	\r0, \r5
	mulx	(6+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r8, \r6
	adox	\r0, \r6
	mulx	(7+\ap_offset)*8(\ap), \r0, \r8
	adcx	\scr, \r7
	adox	\r0, \r7
	adcx	\zero, \r8
	adox	\zero, \r8
.endm

	ALIGN(16)
PROLOGUE(flint_mpn_mul_1_1)
	mov	(%rdx), %rdx
	mulx	0*8(%rsi), %rcx, %rax
	mov	%rcx, 0*8(%rdi)
	mov	%rax, 1*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_2_1)
	mov	0*8(%rdx), %rdx
	xor	%r10d, %r10d
	mulx	0*8(%rsi), %rcx, %r8
	mulx	1*8(%rsi), %r9, %rax
	adcx	%r9, %r8
	adcx	%r10, %rax
	mov	%rcx, 0*8(%rdi)
	mov	%r8, 1*8(%rdi)
	mov	%rax, 2*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_2_2)
	mov	1*8(%rdx), %r10
	mov	0*8(%rdx), %rdx
	xor	%r11d, %r11d
	mulx	0*8(%rsi), %rcx, %r8
	mulx	1*8(%rsi), %rax, %r9
	adcx	%rax, %r8
	mov	%rcx, 0*8(%rdi)
	mov	%r10, %rdx
	mulx	0*8(%rsi), %rcx, %rax
	adox	%rcx, %r8
	adcx	%rax, %r9
	mov	%r8, 1*8(%rdi)
	mulx	1*8(%rsi), %rcx, %rax
	adox	%rcx, %r9
	adcx	%r11, %rax
	adox	%r11, %rax
	mov	%r9, 2*8(%rdi)
	mov	%rax, 3*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_3_1)
	mov	0*8(%rdx), %rdx
	m3_str	%rdi, 0, %rsi, 0, %rcx, %r8, %r9, %r10, %r11, %rax
	mov	%rax, 3*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_3_2)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx

	xor	%ebx, %ebx

	m3	%rdi, 0, %rsi, 0, %r10, %r8, %r9, %rax, %r11, %rbx

	mov	1*8(%rcx), %rdx
	am3	%rdi, 1, %rsi, 0, %r10, %r8, %r9, %rax, %r11, %rbx

	mov	%r8, 2*8(%rdi)
	mov	%r9, 3*8(%rdi)
	mov	%rax, 4*8(%rdi)

	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_3_3)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx

	xor	%ebx, %ebx

	m3	%rdi, 0, %rsi, 0, %rax, %r8, %r9, %r10, %r11, %rbx

	mov	1*8(%rcx), %rdx
	am3	%rdi, 1, %rsi, 0, %rax, %r8, %r9, %r10, %r11, %rbx
	mov	2*8(%rcx), %rdx
	am3	%rdi, 2, %rsi, 0, %r8, %r9, %r10, %rax, %r11, %rbx

	mov	%r9, 3*8(%rdi)
	mov	%r10, 4*8(%rdi)
	mov	%rax, 5*8(%rdi)

	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_4_1)
	mov	0*8(%rdx), %rdx
	m4_str	%rdi, 0, %rsi, 0, %rcx, %r8, %r9, %r10, %r11, %rax
	mov	%rax, 4*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_4_2)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp

	xor	%ebp, %ebp

	m4	%rdi, 0, %rsi, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp

	mov	1*8(%rcx), %rdx
	am4	%rdi, 1, %rsi, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp

	mov	%r8, 2*8(%rdi)
	mov	%r9, 3*8(%rdi)
	mov	%r10, 4*8(%rdi)
	mov	%rax, 5*8(%rdi)

	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_4_3)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp

	xor	%ebp, %ebp

	m4	%rdi, 0, %rsi, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp

	mov	1*8(%rcx), %rdx
	am4	%rdi, 1, %rsi, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp
	mov	2*8(%rcx), %rdx
	am4	%rdi, 2, %rsi, 0, %r8, %r9, %r10, %r11, %rax, %rbx, %rbp

	mov	%r9, 3*8(%rdi)
	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mov	%rax, 6*8(%rdi)

	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_4_4)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp

	xor	%ebp, %ebp

	m4	%rdi, 0, %rsi, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp

	mov	1*8(%rcx), %rdx
	am4	%rdi, 1, %rsi, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp
	mov	2*8(%rcx), %rdx
	am4	%rdi, 2, %rsi, 0, %rax, %r9, %r10, %r11, %r8, %rbx, %rbp
	mov	3*8(%rcx), %rdx
	am4	%rdi, 3, %rsi, 0, %r9, %r10, %r11, %r8, %rax, %rbx, %rbp

	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mov	%r8, 6*8(%rdi)
	mov	%rax, 7*8(%rdi)

	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_5_1)
	mov	0*8(%rdx), %rdx
	m5_str	%rdi, 0, %rsi, 0, %rcx, %r8, %r9, %r10, %r11, %rax
	mov	%rax, 5*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_5_2)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12

	xor	%r12d, %r12d

	m5	%rdi, 0, %rsi, 0, %rbx, %r8, %r9, %r10, %r11, %rax, %rbp, %r12

	mov	1*8(%rcx), %rdx
	am5	%rdi, 1, %rsi, 0, %rbx, %r8, %r9, %r10, %r11, %rax, %rbp, %r12

	mov	%r8, 2*8(%rdi)
	mov	%r9, 3*8(%rdi)
	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mov	%rax, 6*8(%rdi)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_5_3)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12

	xor	%r12d, %r12d

	m5	%rdi, 0, %rsi, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12

	mov	1*8(%rcx), %rdx
	am5	%rdi, 1, %rsi, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12
	mov	2*8(%rcx), %rdx
	am5	%rdi, 2, %rsi, 0, %r8, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12

	mov	%r9, 3*8(%rdi)
	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%rax, 7*8(%rdi)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_5_4)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12

	xor	%r12d, %r12d

	m5	%rdi, 0, %rsi, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12

	mov	1*8(%rcx), %rdx
	am5	%rdi, 1, %rsi, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12
	mov	2*8(%rcx), %rdx
	am5	%rdi, 2, %rsi, 0, %rax, %r9, %r10, %r11, %rbx, %r8, %rbp, %r12
	mov	3*8(%rcx), %rdx
	am5	%rdi, 3, %rsi, 0, %r9, %r10, %r11, %rbx, %r8, %rax, %rbp, %r12

	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%r8, 7*8(%rdi)
	mov	%rax, 8*8(%rdi)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_5_5)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12

	xor	%r12d, %r12d

	m5	%rdi, 0, %rsi, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12

	mov	1*8(%rcx), %rdx
	am5	%rdi, 1, %rsi, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12
	mov	2*8(%rcx), %rdx
	am5	%rdi, 2, %rsi, 0, %r8, %rax, %r10, %r11, %rbx, %r9, %rbp, %r12
	mov	3*8(%rcx), %rdx
	am5	%rdi, 3, %rsi, 0, %rax, %r10, %r11, %rbx, %r9, %r8, %rbp, %r12
	mov	4*8(%rcx), %rdx
	am5	%rdi, 4, %rsi, 0, %r10, %r11, %rbx, %r9, %r8, %rax, %rbp, %r12

	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%r9, 7*8(%rdi)
	mov	%r8, 8*8(%rdi)
	mov	%rax, 9*8(%rdi)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_6_1)
	mov	0*8(%rdx), %rdx
	m6_str	%rdi, 0, %rsi, 0, %rcx, %r8, %r9, %r10, %r11, %rax
	mov	%rax, 6*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_6_2)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	xor	%r13d, %r13d

	m6	%rdi, 0, %rsi, 0, %rbp, %r8, %r9, %r10, %r11, %rbx, %rax, %r12, %r13

	mov	1*8(%rcx), %rdx
	am6	%rdi, 1, %rsi, 0, %rbp, %r8, %r9, %r10, %r11, %rbx, %rax, %r12, %r13

	mov	%r8, 2*8(%rdi)
	mov	%r9, 3*8(%rdi)
	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%rax, 7*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_6_3)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	xor	%r13d, %r13d

	m6	%rdi, 0, %rsi, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13

	mov	1*8(%rcx), %rdx
	am6	%rdi, 1, %rsi, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13
	mov	2*8(%rcx), %rdx
	am6	%rdi, 2, %rsi, 0, %r8, %r9, %r10, %r11, %rbx, %rbp, %rax, %r12, %r13

	mov	%r9, 3*8(%rdi)
	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mov	%rax, 8*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_6_4)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	xor	%r13d, %r13d

	m6	%rdi, 0, %rsi, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13

	mov	1*8(%rcx), %rdx
	am6	%rdi, 1, %rsi, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13
	mov	2*8(%rcx), %rdx
	am6	%rdi, 2, %rsi, 0, %rax, %r9, %r10, %r11, %rbx, %rbp, %r8, %r12, %r13
	mov	3*8(%rcx), %rdx
	am6	%rdi, 3, %rsi, 0, %r9, %r10, %r11, %rbx, %rbp, %r8, %rax, %r12, %r13

	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mov	%r8, 8*8(%rdi)
	mov	%rax, 9*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_6_5)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	xor	%r13d, %r13d

	m6	%rdi, 0, %rsi, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13

	mov	1*8(%rcx), %rdx
	am6	%rdi, 1, %rsi, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13
	mov	2*8(%rcx), %rdx
	am6	%rdi, 2, %rsi, 0, %r8, %rax, %r10, %r11, %rbx, %rbp, %r9, %r12, %r13
	mov	3*8(%rcx), %rdx
	am6	%rdi, 3, %rsi, 0, %rax, %r10, %r11, %rbx, %rbp, %r9, %r8, %r12, %r13
	mov	4*8(%rcx), %rdx
	am6	%rdi, 4, %rsi, 0, %r10, %r11, %rbx, %rbp, %r9, %r8, %rax, %r12, %r13

	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mov	%r9, 8*8(%rdi)
	mov	%r8, 9*8(%rdi)
	mov	%rax, 10*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_6_6)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	xor	%r13d, %r13d

	m6	%rdi, 0, %rsi, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13

	mov	1*8(%rcx), %rdx
	am6	%rdi, 1, %rsi, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13
	mov	2*8(%rcx), %rdx
	am6	%rdi, 2, %rsi, 0, %r8, %r9, %rax, %r11, %rbx, %rbp, %r10, %r12, %r13
	mov	3*8(%rcx), %rdx
	am6	%rdi, 3, %rsi, 0, %r9, %rax, %r11, %rbx, %rbp, %r10, %r8, %r12, %r13
	mov	4*8(%rcx), %rdx
	am6	%rdi, 4, %rsi, 0, %rax, %r11, %rbx, %rbp, %r10, %r8, %r9, %r12, %r13
	mov	5*8(%rcx), %rdx
	am6	%rdi, 5, %rsi, 0, %r11, %rbx, %rbp, %r10, %r8, %r9, %rax, %r12, %r13

	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mov	%r10, 8*8(%rdi)
	mov	%r8, 9*8(%rdi)
	mov	%r9, 10*8(%rdi)
	mov	%rax, 11*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_7_1)
	mov	0*8(%rdx), %rdx
	m7_str	%rdi, 0, %rsi, 0, %rcx, %r8, %r9, %r10, %r11, %rax
	mov	%rax, 7*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_7_2)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	%r14d, %r14d

	m7	%rdi, 0, %rsi, 0, %r12, %r8, %r9, %r10, %r11, %rbx, %rbp, %rax, %r13, %r14

	mov	1*8(%rcx), %rdx
	am7	%rdi, 1, %rsi, 0, %r12, %r8, %r9, %r10, %r11, %rbx, %rbp, %rax, %r13, %r14

	mov	%r8, 2*8(%rdi)
	mov	%r9, 3*8(%rdi)
	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mov	%rax, 8*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_7_3)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	%r14d, %r14d

	m7	%rdi, 0, %rsi, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14

	mov	1*8(%rcx), %rdx
	am7	%rdi, 1, %rsi, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14
	mov	2*8(%rcx), %rdx
	am7	%rdi, 2, %rsi, 0, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %rax, %r13, %r14

	mov	%r9, 3*8(%rdi)
	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mov	%r12, 8*8(%rdi)
	mov	%rax, 9*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_7_4)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	%r14d, %r14d

	m7	%rdi, 0, %rsi, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14

	mov	1*8(%rcx), %rdx
	am7	%rdi, 1, %rsi, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14
	mov	2*8(%rcx), %rdx
	am7	%rdi, 2, %rsi, 0, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r8, %r13, %r14
	mov	3*8(%rcx), %rdx
	am7	%rdi, 3, %rsi, 0, %r9, %r10, %r11, %rbx, %rbp, %r12, %r8, %rax, %r13, %r14

	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mov	%r12, 8*8(%rdi)
	mov	%r8, 9*8(%rdi)
	mov	%rax, 10*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_7_5)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	%r14d, %r14d

	m7	%rdi, 0, %rsi, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14

	mov	1*8(%rcx), %rdx
	am7	%rdi, 1, %rsi, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14
	mov	2*8(%rcx), %rdx
	am7	%rdi, 2, %rsi, 0, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r9, %r13, %r14
	mov	3*8(%rcx), %rdx
	am7	%rdi, 3, %rsi, 0, %rax, %r10, %r11, %rbx, %rbp, %r12, %r9, %r8, %r13, %r14
	mov	4*8(%rcx), %rdx
	am7	%rdi, 4, %rsi, 0, %r10, %r11, %rbx, %rbp, %r12, %r9, %r8, %rax, %r13, %r14

	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mov	%r12, 8*8(%rdi)
	mov	%r9, 9*8(%rdi)
	mov	%r8, 10*8(%rdi)
	mov	%rax, 11*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_7_6)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	%r14d, %r14d

	m7	%rdi, 0, %rsi, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r14

	mov	1*8(%rcx), %rdx
	am7	%rdi, 1, %rsi, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r14
	mov	2*8(%rcx), %rdx
	am7	%rdi, 2, %rsi, 0, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r10, %r13, %r14
	mov	3*8(%rcx), %rdx
	am7	%rdi, 3, %rsi, 0, %r9, %rax, %r11, %rbx, %rbp, %r12, %r10, %r8, %r13, %r14
	mov	4*8(%rcx), %rdx
	am7	%rdi, 4, %rsi, 0, %rax, %r11, %rbx, %rbp, %r12, %r10, %r8, %r9, %r13, %r14
	mov	5*8(%rcx), %rdx
	am7	%rdi, 5, %rsi, 0, %r11, %rbx, %rbp, %r12, %r10, %r8, %r9, %rax, %r13, %r14

	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mov	%r12, 8*8(%rdi)
	mov	%r10, 9*8(%rdi)
	mov	%r8, 10*8(%rdi)
	mov	%r9, 11*8(%rdi)
	mov	%rax, 12*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_7_7)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	%r14d, %r14d

	m7	%rdi, 0, %rsi, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r14

	mov	1*8(%rcx), %rdx
	am7	%rdi, 1, %rsi, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r14
	mov	2*8(%rcx), %rdx
	am7	%rdi, 2, %rsi, 0, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r11, %r13, %r14
	mov	3*8(%rcx), %rdx
	am7	%rdi, 3, %rsi, 0, %r9, %r10, %rax, %rbx, %rbp, %r12, %r11, %r8, %r13, %r14
	mov	4*8(%rcx), %rdx
	am7	%rdi, 4, %rsi, 0, %r10, %rax, %rbx, %rbp, %r12, %r11, %r8, %r9, %r13, %r14
	mov	5*8(%rcx), %rdx
	am7	%rdi, 5, %rsi, 0, %rax, %rbx, %rbp, %r12, %r11, %r8, %r9, %r10, %r13, %r14
	mov	6*8(%rcx), %rdx
	am7	%rdi, 6, %rsi, 0, %rbx, %rbp, %r12, %r11, %r8, %r9, %r10, %rax, %r13, %r14

	mov	%rbp, 7*8(%rdi)
	mov	%r12, 8*8(%rdi)
	mov	%r11, 9*8(%rdi)
	mov	%r8, 10*8(%rdi)
	mov	%r9, 11*8(%rdi)
	mov	%r10, 12*8(%rdi)
	mov	%rax, 13*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_8_1)
	mov	0*8(%rdx), %rdx
	m8_str	%rdi, 0, %rsi, 0, %rcx, %r8, %r9, %r10, %r11, %rax
	mov	%rax, 8*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_8_2)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	xor	%r15d, %r15d

	m8	%rdi, 0, %rsi, 0, %r13, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %rax, %r14, %r15

	mov	1*8(%rcx), %rdx
	am8	%rdi, 1, %rsi, 0, %r13, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %rax, %r14, %r15

	mov	%r8, 2*8(%rdi)
	mov	%r9, 3*8(%rdi)
	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mov	%r12, 8*8(%rdi)
	mov	%rax, 9*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_8_3)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	xor	%r15d, %r15d

	m8	%rdi, 0, %rsi, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15

	mov	1*8(%rcx), %rdx
	am8	%rdi, 1, %rsi, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15
	mov	2*8(%rcx), %rdx
	am8	%rdi, 2, %rsi, 0, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %rax, %r14, %r15

	mov	%r9, 3*8(%rdi)
	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mov	%r12, 8*8(%rdi)
	mov	%r13, 9*8(%rdi)
	mov	%rax, 10*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_8_4)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	xor	%r15d, %r15d

	m8	%rdi, 0, %rsi, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15

	mov	1*8(%rcx), %rdx
	am8	%rdi, 1, %rsi, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15
	mov	2*8(%rcx), %rdx
	am8	%rdi, 2, %rsi, 0, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r8, %r14, %r15
	mov	3*8(%rcx), %rdx
	am8	%rdi, 3, %rsi, 0, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r8, %rax, %r14, %r15

	mov	%r10, 4*8(%rdi)
	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mov	%r12, 8*8(%rdi)
	mov	%r13, 9*8(%rdi)
	mov	%r8, 10*8(%rdi)
	mov	%rax, 11*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_8_5)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	xor	%r15d, %r15d

	m8	%rdi, 0, %rsi, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15

	mov	1*8(%rcx), %rdx
	am8	%rdi, 1, %rsi, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15
	mov	2*8(%rcx), %rdx
	am8	%rdi, 2, %rsi, 0, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r9, %r14, %r15
	mov	3*8(%rcx), %rdx
	am8	%rdi, 3, %rsi, 0, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r9, %r8, %r14, %r15
	mov	4*8(%rcx), %rdx
	am8	%rdi, 4, %rsi, 0, %r10, %r11, %rbx, %rbp, %r12, %r13, %r9, %r8, %rax, %r14, %r15

	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mov	%r12, 8*8(%rdi)
	mov	%r13, 9*8(%rdi)
	mov	%r9, 10*8(%rdi)
	mov	%r8, 11*8(%rdi)
	mov	%rax, 12*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_8_6)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	xor	%r15d, %r15d

	m8	%rdi, 0, %rsi, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15

	mov	1*8(%rcx), %rdx
	am8	%rdi, 1, %rsi, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15
	mov	2*8(%rcx), %rdx
	am8	%rdi, 2, %rsi, 0, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r10, %r14, %r15
	mov	3*8(%rcx), %rdx
	am8	%rdi, 3, %rsi, 0, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r10, %r8, %r14, %r15
	mov	4*8(%rcx), %rdx
	am8	%rdi, 4, %rsi, 0, %rax, %r11, %rbx, %rbp, %r12, %r13, %r10, %r8, %r9, %r14, %r15
	mov	5*8(%rcx), %rdx
	am8	%rdi, 5, %rsi, 0, %r11, %rbx, %rbp, %r12, %r13, %r10, %r8, %r9, %rax, %r14, %r15

	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mov	%r12, 8*8(%rdi)
	mov	%r13, 9*8(%rdi)
	mov	%r10, 10*8(%rdi)
	mov	%r8, 11*8(%rdi)
	mov	%r9, 12*8(%rdi)
	mov	%rax, 13*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_8_7)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	xor	%r15d, %r15d

	m8	%rdi, 0, %rsi, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r14, %r15

	mov	1*8(%rcx), %rdx
	am8	%rdi, 1, %rsi, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r14, %r15
	mov	2*8(%rcx), %rdx
	am8	%rdi, 2, %rsi, 0, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r11, %r14, %r15
	mov	3*8(%rcx), %rdx
	am8	%rdi, 3, %rsi, 0, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r11, %r8, %r14, %r15
	mov	4*8(%rcx), %rdx
	am8	%rdi, 4, %rsi, 0, %r10, %rax, %rbx, %rbp, %r12, %r13, %r11, %r8, %r9, %r14, %r15
	mov	5*8(%rcx), %rdx
	am8	%rdi, 5, %rsi, 0, %rax, %rbx, %rbp, %r12, %r13, %r11, %r8, %r9, %r10, %r14, %r15
	mov	6*8(%rcx), %rdx
	am8	%rdi, 6, %rsi, 0, %rbx, %rbp, %r12, %r13, %r11, %r8, %r9, %r10, %rax, %r14, %r15

	mov	%rbp, 7*8(%rdi)
	mov	%r12, 8*8(%rdi)
	mov	%r13, 9*8(%rdi)
	mov	%r11, 10*8(%rdi)
	mov	%r8, 11*8(%rdi)
	mov	%r9, 12*8(%rdi)
	mov	%r10, 13*8(%rdi)
	mov	%rax, 14*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_8_8)
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	xor	%r15d, %r15d

	m8	%rdi, 0, %rsi, 0, %rbx, %r8, %r9, %r10, %r11, %rax, %rbp, %r12, %r13, %r14, %r15

	mov	1*8(%rcx), %rdx
	am8	%rdi, 1, %rsi, 0, %rbx, %r8, %r9, %r10, %r11, %rax, %rbp, %r12, %r13, %r14, %r15
	mov	2*8(%rcx), %rdx
	am8	%rdi, 2, %rsi, 0, %r8, %r9, %r10, %r11, %rax, %rbp, %r12, %r13, %rbx, %r14, %r15
	mov	3*8(%rcx), %rdx
	am8	%rdi, 3, %rsi, 0, %r9, %r10, %r11, %rax, %rbp, %r12, %r13, %rbx, %r8, %r14, %r15
	mov	4*8(%rcx), %rdx
	am8	%rdi, 4, %rsi, 0, %r10, %r11, %rax, %rbp, %r12, %r13, %rbx, %r8, %r9, %r14, %r15
	mov	5*8(%rcx), %rdx
	am8	%rdi, 5, %rsi, 0, %r11, %rax, %rbp, %r12, %r13, %rbx, %r8, %r9, %r10, %r14, %r15
	mov	6*8(%rcx), %rdx
	am8	%rdi, 6, %rsi, 0, %rax, %rbp, %r12, %r13, %rbx, %r8, %r9, %r10, %r11, %r14, %r15
	mov	7*8(%rcx), %rdx
	am8	%rdi, 7, %rsi, 0, %rbp, %r12, %r13, %rbx, %r8, %r9, %r10, %r11, %rax, %r14, %r15

	mov	%r12, 8*8(%rdi)
	mov	%r13, 9*8(%rdi)
	mov	%rbx, 10*8(%rdi)
	mov	%r8, 11*8(%rdi)
	mov	%r9, 12*8(%rdi)
	mov	%r10, 13*8(%rdi)
	mov	%r11, 14*8(%rdi)
	mov	%rax, 15*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_9_1)
	mov	0*8(%rdx), %rdx
	m9_str	%rdi, 0, %rsi, 0, %rcx, %r8, %r9, %r10, %r11, %rax
	mov	%rax, 9*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_9_2)
	mov	0*8(%rdx), %rcx
	mov	1*8(%rdx), %r8
	push	%rbx
	push	%rbp

	mov	%rcx, %rdx
	xor	%ebp, %ebp
	m3	%rdi, 0, %rsi, 0, %r9, %r10, %r11, %rax, %rbx, %rbp
	mov	%r8, %rdx
	am3	%rdi, 1, %rsi, 0, %r9, %r10, %r11, %rax, %rbx, %rbp
	mov	%r10, 2*8(%rdi)

	mov	%rcx, %rdx
	m3_chain	%rdi, 3, %rsi, 3, %r11, %rax, %r9, %r10, %rax, %rbx, %r11, %rbp
	mov	%r8, %rdx
	am3	%rdi, 4, %rsi, 3, %r9, %r10, %rax, %r11, %rbx, %rbp
	mov	%r10, 5*8(%rdi)

	mov	%rcx, %rdx
	m3_chain	%rdi, 6, %rsi, 6, %rax, %r11, %r9, %r10, %r11, %rbx, %rax, %rbp
	mov	%r8, %rdx
	am3	%rdi, 7, %rsi, 6, %r9, %r10, %r11, %rax, %rbx, %rbp
	mov	%r10, 8*8(%rdi)

	mov	%r11, 9*8(%rdi)
	mov	%rax, 10*8(%rdi)
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_9_3)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx

	xor	%ebx, %ebx

	m3	%rdi, 0, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx

	mov	1*8(%rsi), %rdx
	am3	%rdi, 1, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx
	mov	2*8(%rsi), %rdx
	am3	%rdi, 2, %rcx, 0, %r8, %rax, %r10, %r9, %r11, %rbx
	mov	3*8(%rsi), %rdx
	am3	%rdi, 3, %rcx, 0, %rax, %r10, %r9, %r8, %r11, %rbx
	mov	4*8(%rsi), %rdx
	am3	%rdi, 4, %rcx, 0, %r10, %r9, %r8, %rax, %r11, %rbx
	mov	5*8(%rsi), %rdx
	am3	%rdi, 5, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx
	mov	6*8(%rsi), %rdx
	am3	%rdi, 6, %rcx, 0, %r8, %rax, %r10, %r9, %r11, %rbx
	mov	7*8(%rsi), %rdx
	am3	%rdi, 7, %rcx, 0, %rax, %r10, %r9, %r8, %r11, %rbx
	mov	8*8(%rsi), %rdx
	am3	%rdi, 8, %rcx, 0, %r10, %r9, %r8, %rax, %r11, %rbx

	mov	%r9, 9*8(%rdi)
	mov	%r8, 10*8(%rdi)
	mov	%rax, 11*8(%rdi)

	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_9_4)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp

	xor	%ebp, %ebp

	m4	%rdi, 0, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp

	mov	1*8(%rsi), %rdx
	am4	%rdi, 1, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp
	mov	2*8(%rsi), %rdx
	am4	%rdi, 2, %rcx, 0, %rax, %r9, %r10, %r11, %r8, %rbx, %rbp
	mov	3*8(%rsi), %rdx
	am4	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %r8, %rax, %rbx, %rbp
	mov	4*8(%rsi), %rdx
	am4	%rdi, 4, %rcx, 0, %r10, %r11, %r8, %rax, %r9, %rbx, %rbp
	mov	5*8(%rsi), %rdx
	am4	%rdi, 5, %rcx, 0, %r11, %r8, %rax, %r9, %r10, %rbx, %rbp
	mov	6*8(%rsi), %rdx
	am4	%rdi, 6, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp
	mov	7*8(%rsi), %rdx
	am4	%rdi, 7, %rcx, 0, %rax, %r9, %r10, %r11, %r8, %rbx, %rbp
	mov	8*8(%rsi), %rdx
	am4	%rdi, 8, %rcx, 0, %r9, %r10, %r11, %r8, %rax, %rbx, %rbp

	mov	%r10, 9*8(%rdi)
	mov	%r11, 10*8(%rdi)
	mov	%r8, 11*8(%rdi)
	mov	%rax, 12*8(%rdi)

	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_9_5)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12

	xor	%r12d, %r12d

	m5	%rdi, 0, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12

	mov	1*8(%rsi), %rdx
	am5	%rdi, 1, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12
	mov	2*8(%rsi), %rdx
	am5	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12
	mov	3*8(%rsi), %rdx
	am5	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rax, %r8, %rbp, %r12
	mov	4*8(%rsi), %rdx
	am5	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rax, %r8, %r9, %rbp, %r12
	mov	5*8(%rsi), %rdx
	am5	%rdi, 5, %rcx, 0, %r11, %rbx, %rax, %r8, %r9, %r10, %rbp, %r12
	mov	6*8(%rsi), %rdx
	am5	%rdi, 6, %rcx, 0, %rbx, %rax, %r8, %r9, %r10, %r11, %rbp, %r12
	mov	7*8(%rsi), %rdx
	am5	%rdi, 7, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12
	mov	8*8(%rsi), %rdx
	am5	%rdi, 8, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12

	mov	%r9, 9*8(%rdi)
	mov	%r10, 10*8(%rdi)
	mov	%r11, 11*8(%rdi)
	mov	%rbx, 12*8(%rdi)
	mov	%rax, 13*8(%rdi)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_9_6)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	xor	%r13d, %r13d

	m6	%rdi, 0, %rcx, 0, %rbp, %r8, %r9, %r10, %r11, %rbx, %rax, %r12, %r13

	mov	1*8(%rsi), %rdx
	am6	%rdi, 1, %rcx, 0, %rbp, %r8, %r9, %r10, %r11, %rbx, %rax, %r12, %r13
	mov	2*8(%rsi), %rdx
	am6	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12, %r13
	mov	3*8(%rsi), %rdx
	am6	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rax, %rbp, %r8, %r12, %r13
	mov	4*8(%rsi), %rdx
	am6	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rax, %rbp, %r8, %r9, %r12, %r13
	mov	5*8(%rsi), %rdx
	am6	%rdi, 5, %rcx, 0, %r11, %rbx, %rax, %rbp, %r8, %r9, %r10, %r12, %r13
	mov	6*8(%rsi), %rdx
	am6	%rdi, 6, %rcx, 0, %rbx, %rax, %rbp, %r8, %r9, %r10, %r11, %r12, %r13
	mov	7*8(%rsi), %rdx
	am6	%rdi, 7, %rcx, 0, %rax, %rbp, %r8, %r9, %r10, %r11, %rbx, %r12, %r13
	mov	8*8(%rsi), %rdx
	am6	%rdi, 8, %rcx, 0, %rbp, %r8, %r9, %r10, %r11, %rbx, %rax, %r12, %r13

	mov	%r8, 9*8(%rdi)
	mov	%r9, 10*8(%rdi)
	mov	%r10, 11*8(%rdi)
	mov	%r11, 12*8(%rdi)
	mov	%rbx, 13*8(%rdi)
	mov	%rax, 14*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_9_7)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	%r14d, %r14d

	m7	%rdi, 0, %rcx, 0, %rbp, %r8, %r9, %r10, %r11, %rbx, %rax, %r12, %r13, %r14

	mov	1*8(%rsi), %rdx
	am7	%rdi, 1, %rcx, 0, %rbp, %r8, %r9, %r10, %r11, %rbx, %rax, %r12, %r13, %r14
	mov	2*8(%rsi), %rdx
	am7	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rax, %r12, %rbp, %r13, %r14
	mov	3*8(%rsi), %rdx
	am7	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rax, %r12, %rbp, %r8, %r13, %r14
	mov	4*8(%rsi), %rdx
	am7	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rax, %r12, %rbp, %r8, %r9, %r13, %r14
	mov	5*8(%rsi), %rdx
	am7	%rdi, 5, %rcx, 0, %r11, %rbx, %rax, %r12, %rbp, %r8, %r9, %r10, %r13, %r14
	mov	6*8(%rsi), %rdx
	am7	%rdi, 6, %rcx, 0, %rbx, %rax, %r12, %rbp, %r8, %r9, %r10, %r11, %r13, %r14
	mov	7*8(%rsi), %rdx
	am7	%rdi, 7, %rcx, 0, %rax, %r12, %rbp, %r8, %r9, %r10, %r11, %rbx, %r13, %r14
	mov	8*8(%rsi), %rdx
	am7	%rdi, 8, %rcx, 0, %r12, %rbp, %r8, %r9, %r10, %r11, %rbx, %rax, %r13, %r14

	mov	%rbp, 9*8(%rdi)
	mov	%r8, 10*8(%rdi)
	mov	%r9, 11*8(%rdi)
	mov	%r10, 12*8(%rdi)
	mov	%r11, 13*8(%rdi)
	mov	%rbx, 14*8(%rdi)
	mov	%rax, 15*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_9_8)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	xor	%r15d, %r15d

	m8	%rdi, 0, %rcx, 0, %rbp, %r8, %r9, %r10, %r11, %rbx, %rax, %r12, %r13, %r14, %r15

	mov	1*8(%rsi), %rdx
	am8	%rdi, 1, %rcx, 0, %rbp, %r8, %r9, %r10, %r11, %rbx, %rax, %r12, %r13, %r14, %r15
	mov	2*8(%rsi), %rdx
	am8	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rax, %r12, %r13, %rbp, %r14, %r15
	mov	3*8(%rsi), %rdx
	am8	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rax, %r12, %r13, %rbp, %r8, %r14, %r15
	mov	4*8(%rsi), %rdx
	am8	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rax, %r12, %r13, %rbp, %r8, %r9, %r14, %r15
	mov	5*8(%rsi), %rdx
	am8	%rdi, 5, %rcx, 0, %r11, %rbx, %rax, %r12, %r13, %rbp, %r8, %r9, %r10, %r14, %r15
	mov	6*8(%rsi), %rdx
	am8	%rdi, 6, %rcx, 0, %rbx, %rax, %r12, %r13, %rbp, %r8, %r9, %r10, %r11, %r14, %r15
	mov	7*8(%rsi), %rdx
	am8	%rdi, 7, %rcx, 0, %rax, %r12, %r13, %rbp, %r8, %r9, %r10, %r11, %rbx, %r14, %r15
	mov	8*8(%rsi), %rdx
	am8	%rdi, 8, %rcx, 0, %r12, %r13, %rbp, %r8, %r9, %r10, %r11, %rbx, %rax, %r14, %r15

	mov	%r13, 9*8(%rdi)
	mov	%rbp, 10*8(%rdi)
	mov	%r8, 11*8(%rdi)
	mov	%r9, 12*8(%rdi)
	mov	%r10, 13*8(%rdi)
	mov	%r11, 14*8(%rdi)
	mov	%rbx, 15*8(%rdi)
	mov	%rax, 16*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_10_1)
	mov	0*8(%rdx), %rdx
	m10_str	%rdi, 0, %rsi, 0, %rcx, %r8, %r9, %r10, %r11, %rax
	mov	%rax, 10*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_10_2)
	mov	0*8(%rdx), %rcx
	mov	1*8(%rdx), %r8
	push	%rbx
	push	%rbp
	push	%r12

	mov	%rcx, %rdx
	xor	%r12d, %r12d
	m4	%rdi, 0, %rsi, 0, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12
	mov	%r8, %rdx
	am4	%rdi, 1, %rsi, 0, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12
	mov	%r10, 2*8(%rdi)
	mov	%r11, 3*8(%rdi)

	mov	%rcx, %rdx
	m4_chain	%rdi, 4, %rsi, 4, %rbx, %rax, %r9, %r10, %r11, %rax, %rbp, %rbx, %r12
	mov	%r8, %rdx
	am4	%rdi, 5, %rsi, 4, %r9, %r10, %r11, %rax, %rbx, %rbp, %r12
	mov	%r10, 6*8(%rdi)
	mov	%r11, 7*8(%rdi)

	mov	%rcx, %rdx
	m2_chain	%rdi, 8, %rsi, 8, %rax, %rbx, %r9, %r10, %rbp, %rax, %r12
	mov	%r8, %rdx
	am2	%rdi, 9, %rsi, 8, %r9, %r10, %rax, %rbp, %r12
	mov	%r10, 10*8(%rdi)
	mov	%rax, 11*8(%rdi)
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_10_3)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx

	xor	%ebx, %ebx

	m3	%rdi, 0, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx

	mov	1*8(%rsi), %rdx
	am3	%rdi, 1, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx
	mov	2*8(%rsi), %rdx
	am3	%rdi, 2, %rcx, 0, %r8, %r9, %rax, %r10, %r11, %rbx
	mov	3*8(%rsi), %rdx
	am3	%rdi, 3, %rcx, 0, %r9, %rax, %r10, %r8, %r11, %rbx
	mov	4*8(%rsi), %rdx
	am3	%rdi, 4, %rcx, 0, %rax, %r10, %r8, %r9, %r11, %rbx
	mov	5*8(%rsi), %rdx
	am3	%rdi, 5, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx
	mov	6*8(%rsi), %rdx
	am3	%rdi, 6, %rcx, 0, %r8, %r9, %rax, %r10, %r11, %rbx
	mov	7*8(%rsi), %rdx
	am3	%rdi, 7, %rcx, 0, %r9, %rax, %r10, %r8, %r11, %rbx
	mov	8*8(%rsi), %rdx
	am3	%rdi, 8, %rcx, 0, %rax, %r10, %r8, %r9, %r11, %rbx
	mov	9*8(%rsi), %rdx
	am3	%rdi, 9, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx

	mov	%r8, 10*8(%rdi)
	mov	%r9, 11*8(%rdi)
	mov	%rax, 12*8(%rdi)

	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_10_4)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp

	xor	%ebp, %ebp

	m4	%rdi, 0, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp

	mov	1*8(%rsi), %rdx
	am4	%rdi, 1, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp
	mov	2*8(%rsi), %rdx
	am4	%rdi, 2, %rcx, 0, %r8, %rax, %r10, %r11, %r9, %rbx, %rbp
	mov	3*8(%rsi), %rdx
	am4	%rdi, 3, %rcx, 0, %rax, %r10, %r11, %r9, %r8, %rbx, %rbp
	mov	4*8(%rsi), %rdx
	am4	%rdi, 4, %rcx, 0, %r10, %r11, %r9, %r8, %rax, %rbx, %rbp
	mov	5*8(%rsi), %rdx
	am4	%rdi, 5, %rcx, 0, %r11, %r9, %r8, %rax, %r10, %rbx, %rbp
	mov	6*8(%rsi), %rdx
	am4	%rdi, 6, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp
	mov	7*8(%rsi), %rdx
	am4	%rdi, 7, %rcx, 0, %r8, %rax, %r10, %r11, %r9, %rbx, %rbp
	mov	8*8(%rsi), %rdx
	am4	%rdi, 8, %rcx, 0, %rax, %r10, %r11, %r9, %r8, %rbx, %rbp
	mov	9*8(%rsi), %rdx
	am4	%rdi, 9, %rcx, 0, %r10, %r11, %r9, %r8, %rax, %rbx, %rbp

	mov	%r11, 10*8(%rdi)
	mov	%r9, 11*8(%rdi)
	mov	%r8, 12*8(%rdi)
	mov	%rax, 13*8(%rdi)

	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_10_5)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12

	xor	%r12d, %r12d

	m5	%rdi, 0, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12

	mov	1*8(%rsi), %rdx
	am5	%rdi, 1, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12
	mov	2*8(%rsi), %rdx
	am5	%rdi, 2, %rcx, 0, %rax, %r9, %r10, %r11, %rbx, %r8, %rbp, %r12
	mov	3*8(%rsi), %rdx
	am5	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %r8, %rax, %rbp, %r12
	mov	4*8(%rsi), %rdx
	am5	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %r8, %rax, %r9, %rbp, %r12
	mov	5*8(%rsi), %rdx
	am5	%rdi, 5, %rcx, 0, %r11, %rbx, %r8, %rax, %r9, %r10, %rbp, %r12
	mov	6*8(%rsi), %rdx
	am5	%rdi, 6, %rcx, 0, %rbx, %r8, %rax, %r9, %r10, %r11, %rbp, %r12
	mov	7*8(%rsi), %rdx
	am5	%rdi, 7, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12
	mov	8*8(%rsi), %rdx
	am5	%rdi, 8, %rcx, 0, %rax, %r9, %r10, %r11, %rbx, %r8, %rbp, %r12
	mov	9*8(%rsi), %rdx
	am5	%rdi, 9, %rcx, 0, %r9, %r10, %r11, %rbx, %r8, %rax, %rbp, %r12

	mov	%r10, 10*8(%rdi)
	mov	%r11, 11*8(%rdi)
	mov	%rbx, 12*8(%rdi)
	mov	%r8, 13*8(%rdi)
	mov	%rax, 14*8(%rdi)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_10_6)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	xor	%r13d, %r13d

	m6	%rdi, 0, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13

	mov	1*8(%rsi), %rdx
	am6	%rdi, 1, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13
	mov	2*8(%rsi), %rdx
	am6	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rbp, %rax, %r12, %r13
	mov	3*8(%rsi), %rdx
	am6	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rbp, %rax, %r8, %r12, %r13
	mov	4*8(%rsi), %rdx
	am6	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rbp, %rax, %r8, %r9, %r12, %r13
	mov	5*8(%rsi), %rdx
	am6	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %rax, %r8, %r9, %r10, %r12, %r13
	mov	6*8(%rsi), %rdx
	am6	%rdi, 6, %rcx, 0, %rbx, %rbp, %rax, %r8, %r9, %r10, %r11, %r12, %r13
	mov	7*8(%rsi), %rdx
	am6	%rdi, 7, %rcx, 0, %rbp, %rax, %r8, %r9, %r10, %r11, %rbx, %r12, %r13
	mov	8*8(%rsi), %rdx
	am6	%rdi, 8, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13
	mov	9*8(%rsi), %rdx
	am6	%rdi, 9, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rbp, %rax, %r12, %r13

	mov	%r9, 10*8(%rdi)
	mov	%r10, 11*8(%rdi)
	mov	%r11, 12*8(%rdi)
	mov	%rbx, 13*8(%rdi)
	mov	%rbp, 14*8(%rdi)
	mov	%rax, 15*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_10_7)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	%r14d, %r14d

	m7	%rdi, 0, %rcx, 0, %r12, %r8, %r9, %r10, %r11, %rbx, %rbp, %rax, %r13, %r14

	mov	1*8(%rsi), %rdx
	am7	%rdi, 1, %rcx, 0, %r12, %r8, %r9, %r10, %r11, %rbx, %rbp, %rax, %r13, %r14
	mov	2*8(%rsi), %rdx
	am7	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rbp, %rax, %r12, %r13, %r14
	mov	3*8(%rsi), %rdx
	am7	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rbp, %rax, %r12, %r8, %r13, %r14
	mov	4*8(%rsi), %rdx
	am7	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rbp, %rax, %r12, %r8, %r9, %r13, %r14
	mov	5*8(%rsi), %rdx
	am7	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %rax, %r12, %r8, %r9, %r10, %r13, %r14
	mov	6*8(%rsi), %rdx
	am7	%rdi, 6, %rcx, 0, %rbx, %rbp, %rax, %r12, %r8, %r9, %r10, %r11, %r13, %r14
	mov	7*8(%rsi), %rdx
	am7	%rdi, 7, %rcx, 0, %rbp, %rax, %r12, %r8, %r9, %r10, %r11, %rbx, %r13, %r14
	mov	8*8(%rsi), %rdx
	am7	%rdi, 8, %rcx, 0, %rax, %r12, %r8, %r9, %r10, %r11, %rbx, %rbp, %r13, %r14
	mov	9*8(%rsi), %rdx
	am7	%rdi, 9, %rcx, 0, %r12, %r8, %r9, %r10, %r11, %rbx, %rbp, %rax, %r13, %r14

	mov	%r8, 10*8(%rdi)
	mov	%r9, 11*8(%rdi)
	mov	%r10, 12*8(%rdi)
	mov	%r11, 13*8(%rdi)
	mov	%rbx, 14*8(%rdi)
	mov	%rbp, 15*8(%rdi)
	mov	%rax, 16*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_10_8)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	xor	%r15d, %r15d

	m8	%rdi, 0, %rcx, 0, %r12, %r8, %r9, %r10, %r11, %rbx, %rbp, %rax, %r13, %r14, %r15

	mov	1*8(%rsi), %rdx
	am8	%rdi, 1, %rcx, 0, %r12, %r8, %r9, %r10, %r11, %rbx, %rbp, %rax, %r13, %r14, %r15
	mov	2*8(%rsi), %rdx
	am8	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rbp, %rax, %r13, %r12, %r14, %r15
	mov	3*8(%rsi), %rdx
	am8	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rbp, %rax, %r13, %r12, %r8, %r14, %r15
	mov	4*8(%rsi), %rdx
	am8	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rbp, %rax, %r13, %r12, %r8, %r9, %r14, %r15
	mov	5*8(%rsi), %rdx
	am8	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %rax, %r13, %r12, %r8, %r9, %r10, %r14, %r15
	mov	6*8(%rsi), %rdx
	am8	%rdi, 6, %rcx, 0, %rbx, %rbp, %rax, %r13, %r12, %r8, %r9, %r10, %r11, %r14, %r15
	mov	7*8(%rsi), %rdx
	am8	%rdi, 7, %rcx, 0, %rbp, %rax, %r13, %r12, %r8, %r9, %r10, %r11, %rbx, %r14, %r15
	mov	8*8(%rsi), %rdx
	am8	%rdi, 8, %rcx, 0, %rax, %r13, %r12, %r8, %r9, %r10, %r11, %rbx, %rbp, %r14, %r15
	mov	9*8(%rsi), %rdx
	am8	%rdi, 9, %rcx, 0, %r13, %r12, %r8, %r9, %r10, %r11, %rbx, %rbp, %rax, %r14, %r15

	mov	%r12, 10*8(%rdi)
	mov	%r8, 11*8(%rdi)
	mov	%r9, 12*8(%rdi)
	mov	%r10, 13*8(%rdi)
	mov	%r11, 14*8(%rdi)
	mov	%rbx, 15*8(%rdi)
	mov	%rbp, 16*8(%rdi)
	mov	%rax, 17*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_11_1)
	mov	0*8(%rdx), %rdx
	m11_str	%rdi, 0, %rsi, 0, %rcx, %r8, %r9, %r10, %r11, %rax
	mov	%rax, 11*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_11_2)
	mov	0*8(%rdx), %rcx
	mov	1*8(%rdx), %r8
	push	%rbx
	push	%rbp
	push	%r12

	mov	%rcx, %rdx
	xor	%r12d, %r12d
	m4	%rdi, 0, %rsi, 0, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12
	mov	%r8, %rdx
	am4	%rdi, 1, %rsi, 0, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12
	mov	%r10, 2*8(%rdi)
	mov	%r11, 3*8(%rdi)

	mov	%rcx, %rdx
	m4_chain	%rdi, 4, %rsi, 4, %rbx, %rax, %r9, %r10, %r11, %rax, %rbp, %rbx, %r12
	mov	%r8, %rdx
	am4	%rdi, 5, %rsi, 4, %r9, %r10, %r11, %rax, %rbx, %rbp, %r12
	mov	%r10, 6*8(%rdi)
	mov	%r11, 7*8(%rdi)

	mov	%rcx, %rdx
	m3_chain	%rdi, 8, %rsi, 8, %rax, %rbx, %r9, %r10, %r11, %rbp, %rax, %r12
	mov	%r8, %rdx
	am3	%rdi, 9, %rsi, 8, %r9, %r10, %r11, %rax, %rbp, %r12
	mov	%r10, 10*8(%rdi)
	mov	%r11, 11*8(%rdi)
	mov	%rax, 12*8(%rdi)
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_11_3)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx

	xor	%ebx, %ebx

	m3	%rdi, 0, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx

	mov	1*8(%rsi), %rdx
	am3	%rdi, 1, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx
	mov	2*8(%rsi), %rdx
	am3	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %rax, %r11, %rbx
	mov	3*8(%rsi), %rdx
	am3	%rdi, 3, %rcx, 0, %r9, %r10, %rax, %r8, %r11, %rbx
	mov	4*8(%rsi), %rdx
	am3	%rdi, 4, %rcx, 0, %r10, %rax, %r8, %r9, %r11, %rbx
	mov	5*8(%rsi), %rdx
	am3	%rdi, 5, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx
	mov	6*8(%rsi), %rdx
	am3	%rdi, 6, %rcx, 0, %r8, %r9, %r10, %rax, %r11, %rbx
	mov	7*8(%rsi), %rdx
	am3	%rdi, 7, %rcx, 0, %r9, %r10, %rax, %r8, %r11, %rbx
	mov	8*8(%rsi), %rdx
	am3	%rdi, 8, %rcx, 0, %r10, %rax, %r8, %r9, %r11, %rbx
	mov	9*8(%rsi), %rdx
	am3	%rdi, 9, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx
	mov	10*8(%rsi), %rdx
	am3	%rdi, 10, %rcx, 0, %r8, %r9, %r10, %rax, %r11, %rbx

	mov	%r9, 11*8(%rdi)
	mov	%r10, 12*8(%rdi)
	mov	%rax, 13*8(%rdi)

	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_11_4)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp

	xor	%ebp, %ebp

	m4	%rdi, 0, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp

	mov	1*8(%rsi), %rdx
	am4	%rdi, 1, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp
	mov	2*8(%rsi), %rdx
	am4	%rdi, 2, %rcx, 0, %r8, %r9, %rax, %r11, %r10, %rbx, %rbp
	mov	3*8(%rsi), %rdx
	am4	%rdi, 3, %rcx, 0, %r9, %rax, %r11, %r10, %r8, %rbx, %rbp
	mov	4*8(%rsi), %rdx
	am4	%rdi, 4, %rcx, 0, %rax, %r11, %r10, %r8, %r9, %rbx, %rbp
	mov	5*8(%rsi), %rdx
	am4	%rdi, 5, %rcx, 0, %r11, %r10, %r8, %r9, %rax, %rbx, %rbp
	mov	6*8(%rsi), %rdx
	am4	%rdi, 6, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp
	mov	7*8(%rsi), %rdx
	am4	%rdi, 7, %rcx, 0, %r8, %r9, %rax, %r11, %r10, %rbx, %rbp
	mov	8*8(%rsi), %rdx
	am4	%rdi, 8, %rcx, 0, %r9, %rax, %r11, %r10, %r8, %rbx, %rbp
	mov	9*8(%rsi), %rdx
	am4	%rdi, 9, %rcx, 0, %rax, %r11, %r10, %r8, %r9, %rbx, %rbp
	mov	10*8(%rsi), %rdx
	am4	%rdi, 10, %rcx, 0, %r11, %r10, %r8, %r9, %rax, %rbx, %rbp

	mov	%r10, 11*8(%rdi)
	mov	%r8, 12*8(%rdi)
	mov	%r9, 13*8(%rdi)
	mov	%rax, 14*8(%rdi)

	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_11_5)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12

	xor	%r12d, %r12d

	m5	%rdi, 0, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12

	mov	1*8(%rsi), %rdx
	am5	%rdi, 1, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12
	mov	2*8(%rsi), %rdx
	am5	%rdi, 2, %rcx, 0, %r8, %rax, %r10, %r11, %rbx, %r9, %rbp, %r12
	mov	3*8(%rsi), %rdx
	am5	%rdi, 3, %rcx, 0, %rax, %r10, %r11, %rbx, %r9, %r8, %rbp, %r12
	mov	4*8(%rsi), %rdx
	am5	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %r9, %r8, %rax, %rbp, %r12
	mov	5*8(%rsi), %rdx
	am5	%rdi, 5, %rcx, 0, %r11, %rbx, %r9, %r8, %rax, %r10, %rbp, %r12
	mov	6*8(%rsi), %rdx
	am5	%rdi, 6, %rcx, 0, %rbx, %r9, %r8, %rax, %r10, %r11, %rbp, %r12
	mov	7*8(%rsi), %rdx
	am5	%rdi, 7, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12
	mov	8*8(%rsi), %rdx
	am5	%rdi, 8, %rcx, 0, %r8, %rax, %r10, %r11, %rbx, %r9, %rbp, %r12
	mov	9*8(%rsi), %rdx
	am5	%rdi, 9, %rcx, 0, %rax, %r10, %r11, %rbx, %r9, %r8, %rbp, %r12
	mov	10*8(%rsi), %rdx
	am5	%rdi, 10, %rcx, 0, %r10, %r11, %rbx, %r9, %r8, %rax, %rbp, %r12

	mov	%r11, 11*8(%rdi)
	mov	%rbx, 12*8(%rdi)
	mov	%r9, 13*8(%rdi)
	mov	%r8, 14*8(%rdi)
	mov	%rax, 15*8(%rdi)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_11_6)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	xor	%r13d, %r13d

	m6	%rdi, 0, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13

	mov	1*8(%rsi), %rdx
	am6	%rdi, 1, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13
	mov	2*8(%rsi), %rdx
	am6	%rdi, 2, %rcx, 0, %rax, %r9, %r10, %r11, %rbx, %rbp, %r8, %r12, %r13
	mov	3*8(%rsi), %rdx
	am6	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rbp, %r8, %rax, %r12, %r13
	mov	4*8(%rsi), %rdx
	am6	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rbp, %r8, %rax, %r9, %r12, %r13
	mov	5*8(%rsi), %rdx
	am6	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %r8, %rax, %r9, %r10, %r12, %r13
	mov	6*8(%rsi), %rdx
	am6	%rdi, 6, %rcx, 0, %rbx, %rbp, %r8, %rax, %r9, %r10, %r11, %r12, %r13
	mov	7*8(%rsi), %rdx
	am6	%rdi, 7, %rcx, 0, %rbp, %r8, %rax, %r9, %r10, %r11, %rbx, %r12, %r13
	mov	8*8(%rsi), %rdx
	am6	%rdi, 8, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13
	mov	9*8(%rsi), %rdx
	am6	%rdi, 9, %rcx, 0, %rax, %r9, %r10, %r11, %rbx, %rbp, %r8, %r12, %r13
	mov	10*8(%rsi), %rdx
	am6	%rdi, 10, %rcx, 0, %r9, %r10, %r11, %rbx, %rbp, %r8, %rax, %r12, %r13

	mov	%r10, 11*8(%rdi)
	mov	%r11, 12*8(%rdi)
	mov	%rbx, 13*8(%rdi)
	mov	%rbp, 14*8(%rdi)
	mov	%r8, 15*8(%rdi)
	mov	%rax, 16*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_11_7)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	%r14d, %r14d

	m7	%rdi, 0, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14

	mov	1*8(%rsi), %rdx
	am7	%rdi, 1, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14
	mov	2*8(%rsi), %rdx
	am7	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %rax, %r13, %r14
	mov	3*8(%rsi), %rdx
	am7	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rbp, %r12, %rax, %r8, %r13, %r14
	mov	4*8(%rsi), %rdx
	am7	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rbp, %r12, %rax, %r8, %r9, %r13, %r14
	mov	5*8(%rsi), %rdx
	am7	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %r12, %rax, %r8, %r9, %r10, %r13, %r14
	mov	6*8(%rsi), %rdx
	am7	%rdi, 6, %rcx, 0, %rbx, %rbp, %r12, %rax, %r8, %r9, %r10, %r11, %r13, %r14
	mov	7*8(%rsi), %rdx
	am7	%rdi, 7, %rcx, 0, %rbp, %r12, %rax, %r8, %r9, %r10, %r11, %rbx, %r13, %r14
	mov	8*8(%rsi), %rdx
	am7	%rdi, 8, %rcx, 0, %r12, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r13, %r14
	mov	9*8(%rsi), %rdx
	am7	%rdi, 9, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14
	mov	10*8(%rsi), %rdx
	am7	%rdi, 10, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %rax, %r13, %r14

	mov	%r9, 11*8(%rdi)
	mov	%r10, 12*8(%rdi)
	mov	%r11, 13*8(%rdi)
	mov	%rbx, 14*8(%rdi)
	mov	%rbp, 15*8(%rdi)
	mov	%r12, 16*8(%rdi)
	mov	%rax, 17*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_11_8)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	xor	%r15d, %r15d

	m8	%rdi, 0, %rcx, 0, %r13, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %rax, %r14, %r15

	mov	1*8(%rsi), %rdx
	am8	%rdi, 1, %rcx, 0, %r13, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %rax, %r14, %r15
	mov	2*8(%rsi), %rdx
	am8	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %rax, %r13, %r14, %r15
	mov	3*8(%rsi), %rdx
	am8	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rbp, %r12, %rax, %r13, %r8, %r14, %r15
	mov	4*8(%rsi), %rdx
	am8	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rbp, %r12, %rax, %r13, %r8, %r9, %r14, %r15
	mov	5*8(%rsi), %rdx
	am8	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %r12, %rax, %r13, %r8, %r9, %r10, %r14, %r15
	mov	6*8(%rsi), %rdx
	am8	%rdi, 6, %rcx, 0, %rbx, %rbp, %r12, %rax, %r13, %r8, %r9, %r10, %r11, %r14, %r15
	mov	7*8(%rsi), %rdx
	am8	%rdi, 7, %rcx, 0, %rbp, %r12, %rax, %r13, %r8, %r9, %r10, %r11, %rbx, %r14, %r15
	mov	8*8(%rsi), %rdx
	am8	%rdi, 8, %rcx, 0, %r12, %rax, %r13, %r8, %r9, %r10, %r11, %rbx, %rbp, %r14, %r15
	mov	9*8(%rsi), %rdx
	am8	%rdi, 9, %rcx, 0, %rax, %r13, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r14, %r15
	mov	10*8(%rsi), %rdx
	am8	%rdi, 10, %rcx, 0, %r13, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %rax, %r14, %r15

	mov	%r8, 11*8(%rdi)
	mov	%r9, 12*8(%rdi)
	mov	%r10, 13*8(%rdi)
	mov	%r11, 14*8(%rdi)
	mov	%rbx, 15*8(%rdi)
	mov	%rbp, 16*8(%rdi)
	mov	%r12, 17*8(%rdi)
	mov	%rax, 18*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_12_1)
	mov	0*8(%rdx), %rdx
	m12_str	%rdi, 0, %rsi, 0, %rcx, %r8, %r9, %r10, %r11, %rax
	mov	%rax, 12*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_12_2)
	mov	0*8(%rdx), %rcx
	mov	1*8(%rdx), %r8
	push	%rbx
	push	%rbp
	push	%r12

	mov	%rcx, %rdx
	xor	%r12d, %r12d
	m4	%rdi, 0, %rsi, 0, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12
	mov	%r8, %rdx
	am4	%rdi, 1, %rsi, 0, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12
	mov	%r10, 2*8(%rdi)
	mov	%r11, 3*8(%rdi)

	mov	%rcx, %rdx
	m4_chain	%rdi, 4, %rsi, 4, %rbx, %rax, %r9, %r10, %r11, %rax, %rbp, %rbx, %r12
	mov	%r8, %rdx
	am4	%rdi, 5, %rsi, 4, %r9, %r10, %r11, %rax, %rbx, %rbp, %r12
	mov	%r10, 6*8(%rdi)
	mov	%r11, 7*8(%rdi)

	mov	%rcx, %rdx
	m4_chain	%rdi, 8, %rsi, 8, %rax, %rbx, %r9, %r10, %r11, %rbx, %rbp, %rax, %r12
	mov	%r8, %rdx
	am4	%rdi, 9, %rsi, 8, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12
	mov	%r10, 10*8(%rdi)
	mov	%r11, 11*8(%rdi)

	mov	%rbx, 12*8(%rdi)
	mov	%rax, 13*8(%rdi)
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_12_3)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx

	xor	%ebx, %ebx

	m3	%rdi, 0, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx

	mov	1*8(%rsi), %rdx
	am3	%rdi, 1, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx
	mov	2*8(%rsi), %rdx
	am3	%rdi, 2, %rcx, 0, %rax, %r9, %r10, %r8, %r11, %rbx
	mov	3*8(%rsi), %rdx
	am3	%rdi, 3, %rcx, 0, %r9, %r10, %r8, %rax, %r11, %rbx
	mov	4*8(%rsi), %rdx
	am3	%rdi, 4, %rcx, 0, %r10, %r8, %rax, %r9, %r11, %rbx
	mov	5*8(%rsi), %rdx
	am3	%rdi, 5, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx
	mov	6*8(%rsi), %rdx
	am3	%rdi, 6, %rcx, 0, %rax, %r9, %r10, %r8, %r11, %rbx
	mov	7*8(%rsi), %rdx
	am3	%rdi, 7, %rcx, 0, %r9, %r10, %r8, %rax, %r11, %rbx
	mov	8*8(%rsi), %rdx
	am3	%rdi, 8, %rcx, 0, %r10, %r8, %rax, %r9, %r11, %rbx
	mov	9*8(%rsi), %rdx
	am3	%rdi, 9, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx
	mov	10*8(%rsi), %rdx
	am3	%rdi, 10, %rcx, 0, %rax, %r9, %r10, %r8, %r11, %rbx
	mov	11*8(%rsi), %rdx
	am3	%rdi, 11, %rcx, 0, %r9, %r10, %r8, %rax, %r11, %rbx

	mov	%r10, 12*8(%rdi)
	mov	%r8, 13*8(%rdi)
	mov	%rax, 14*8(%rdi)

	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_12_4)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp

	xor	%ebp, %ebp

	m4	%rdi, 0, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp

	mov	1*8(%rsi), %rdx
	am4	%rdi, 1, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp
	mov	2*8(%rsi), %rdx
	am4	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %rax, %r11, %rbx, %rbp
	mov	3*8(%rsi), %rdx
	am4	%rdi, 3, %rcx, 0, %r9, %r10, %rax, %r11, %r8, %rbx, %rbp
	mov	4*8(%rsi), %rdx
	am4	%rdi, 4, %rcx, 0, %r10, %rax, %r11, %r8, %r9, %rbx, %rbp
	mov	5*8(%rsi), %rdx
	am4	%rdi, 5, %rcx, 0, %rax, %r11, %r8, %r9, %r10, %rbx, %rbp
	mov	6*8(%rsi), %rdx
	am4	%rdi, 6, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp
	mov	7*8(%rsi), %rdx
	am4	%rdi, 7, %rcx, 0, %r8, %r9, %r10, %rax, %r11, %rbx, %rbp
	mov	8*8(%rsi), %rdx
	am4	%rdi, 8, %rcx, 0, %r9, %r10, %rax, %r11, %r8, %rbx, %rbp
	mov	9*8(%rsi), %rdx
	am4	%rdi, 9, %rcx, 0, %r10, %rax, %r11, %r8, %r9, %rbx, %rbp
	mov	10*8(%rsi), %rdx
	am4	%rdi, 10, %rcx, 0, %rax, %r11, %r8, %r9, %r10, %rbx, %rbp
	mov	11*8(%rsi), %rdx
	am4	%rdi, 11, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp

	mov	%r8, 12*8(%rdi)
	mov	%r9, 13*8(%rdi)
	mov	%r10, 14*8(%rdi)
	mov	%rax, 15*8(%rdi)

	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_12_5)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12

	xor	%r12d, %r12d

	m5	%rdi, 0, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12

	mov	1*8(%rsi), %rdx
	am5	%rdi, 1, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12
	mov	2*8(%rsi), %rdx
	am5	%rdi, 2, %rcx, 0, %r8, %r9, %rax, %r11, %rbx, %r10, %rbp, %r12
	mov	3*8(%rsi), %rdx
	am5	%rdi, 3, %rcx, 0, %r9, %rax, %r11, %rbx, %r10, %r8, %rbp, %r12
	mov	4*8(%rsi), %rdx
	am5	%rdi, 4, %rcx, 0, %rax, %r11, %rbx, %r10, %r8, %r9, %rbp, %r12
	mov	5*8(%rsi), %rdx
	am5	%rdi, 5, %rcx, 0, %r11, %rbx, %r10, %r8, %r9, %rax, %rbp, %r12
	mov	6*8(%rsi), %rdx
	am5	%rdi, 6, %rcx, 0, %rbx, %r10, %r8, %r9, %rax, %r11, %rbp, %r12
	mov	7*8(%rsi), %rdx
	am5	%rdi, 7, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12
	mov	8*8(%rsi), %rdx
	am5	%rdi, 8, %rcx, 0, %r8, %r9, %rax, %r11, %rbx, %r10, %rbp, %r12
	mov	9*8(%rsi), %rdx
	am5	%rdi, 9, %rcx, 0, %r9, %rax, %r11, %rbx, %r10, %r8, %rbp, %r12
	mov	10*8(%rsi), %rdx
	am5	%rdi, 10, %rcx, 0, %rax, %r11, %rbx, %r10, %r8, %r9, %rbp, %r12
	mov	11*8(%rsi), %rdx
	am5	%rdi, 11, %rcx, 0, %r11, %rbx, %r10, %r8, %r9, %rax, %rbp, %r12

	mov	%rbx, 12*8(%rdi)
	mov	%r10, 13*8(%rdi)
	mov	%r8, 14*8(%rdi)
	mov	%r9, 15*8(%rdi)
	mov	%rax, 16*8(%rdi)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_12_6)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	xor	%r13d, %r13d

	m6	%rdi, 0, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13

	mov	1*8(%rsi), %rdx
	am6	%rdi, 1, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13
	mov	2*8(%rsi), %rdx
	am6	%rdi, 2, %rcx, 0, %r8, %rax, %r10, %r11, %rbx, %rbp, %r9, %r12, %r13
	mov	3*8(%rsi), %rdx
	am6	%rdi, 3, %rcx, 0, %rax, %r10, %r11, %rbx, %rbp, %r9, %r8, %r12, %r13
	mov	4*8(%rsi), %rdx
	am6	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rbp, %r9, %r8, %rax, %r12, %r13
	mov	5*8(%rsi), %rdx
	am6	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %r9, %r8, %rax, %r10, %r12, %r13
	mov	6*8(%rsi), %rdx
	am6	%rdi, 6, %rcx, 0, %rbx, %rbp, %r9, %r8, %rax, %r10, %r11, %r12, %r13
	mov	7*8(%rsi), %rdx
	am6	%rdi, 7, %rcx, 0, %rbp, %r9, %r8, %rax, %r10, %r11, %rbx, %r12, %r13
	mov	8*8(%rsi), %rdx
	am6	%rdi, 8, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13
	mov	9*8(%rsi), %rdx
	am6	%rdi, 9, %rcx, 0, %r8, %rax, %r10, %r11, %rbx, %rbp, %r9, %r12, %r13
	mov	10*8(%rsi), %rdx
	am6	%rdi, 10, %rcx, 0, %rax, %r10, %r11, %rbx, %rbp, %r9, %r8, %r12, %r13
	mov	11*8(%rsi), %rdx
	am6	%rdi, 11, %rcx, 0, %r10, %r11, %rbx, %rbp, %r9, %r8, %rax, %r12, %r13

	mov	%r11, 12*8(%rdi)
	mov	%rbx, 13*8(%rdi)
	mov	%rbp, 14*8(%rdi)
	mov	%r9, 15*8(%rdi)
	mov	%r8, 16*8(%rdi)
	mov	%rax, 17*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_12_7)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	%r14d, %r14d

	m7	%rdi, 0, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14

	mov	1*8(%rsi), %rdx
	am7	%rdi, 1, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14
	mov	2*8(%rsi), %rdx
	am7	%rdi, 2, %rcx, 0, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r8, %r13, %r14
	mov	3*8(%rsi), %rdx
	am7	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rbp, %r12, %r8, %rax, %r13, %r14
	mov	4*8(%rsi), %rdx
	am7	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rbp, %r12, %r8, %rax, %r9, %r13, %r14
	mov	5*8(%rsi), %rdx
	am7	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %r12, %r8, %rax, %r9, %r10, %r13, %r14
	mov	6*8(%rsi), %rdx
	am7	%rdi, 6, %rcx, 0, %rbx, %rbp, %r12, %r8, %rax, %r9, %r10, %r11, %r13, %r14
	mov	7*8(%rsi), %rdx
	am7	%rdi, 7, %rcx, 0, %rbp, %r12, %r8, %rax, %r9, %r10, %r11, %rbx, %r13, %r14
	mov	8*8(%rsi), %rdx
	am7	%rdi, 8, %rcx, 0, %r12, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r13, %r14
	mov	9*8(%rsi), %rdx
	am7	%rdi, 9, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14
	mov	10*8(%rsi), %rdx
	am7	%rdi, 10, %rcx, 0, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r8, %r13, %r14
	mov	11*8(%rsi), %rdx
	am7	%rdi, 11, %rcx, 0, %r9, %r10, %r11, %rbx, %rbp, %r12, %r8, %rax, %r13, %r14

	mov	%r10, 12*8(%rdi)
	mov	%r11, 13*8(%rdi)
	mov	%rbx, 14*8(%rdi)
	mov	%rbp, 15*8(%rdi)
	mov	%r12, 16*8(%rdi)
	mov	%r8, 17*8(%rdi)
	mov	%rax, 18*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_12_8)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	xor	%r15d, %r15d

	m8	%rdi, 0, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15

	mov	1*8(%rsi), %rdx
	am8	%rdi, 1, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15
	mov	2*8(%rsi), %rdx
	am8	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %rax, %r14, %r15
	mov	3*8(%rsi), %rdx
	am8	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %rax, %r8, %r14, %r15
	mov	4*8(%rsi), %rdx
	am8	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rbp, %r12, %r13, %rax, %r8, %r9, %r14, %r15
	mov	5*8(%rsi), %rdx
	am8	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %r12, %r13, %rax, %r8, %r9, %r10, %r14, %r15
	mov	6*8(%rsi), %rdx
	am8	%rdi, 6, %rcx, 0, %rbx, %rbp, %r12, %r13, %rax, %r8, %r9, %r10, %r11, %r14, %r15
	mov	7*8(%rsi), %rdx
	am8	%rdi, 7, %rcx, 0, %rbp, %r12, %r13, %rax, %r8, %r9, %r10, %r11, %rbx, %r14, %r15
	mov	8*8(%rsi), %rdx
	am8	%rdi, 8, %rcx, 0, %r12, %r13, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r14, %r15
	mov	9*8(%rsi), %rdx
	am8	%rdi, 9, %rcx, 0, %r13, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r14, %r15
	mov	10*8(%rsi), %rdx
	am8	%rdi, 10, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15
	mov	11*8(%rsi), %rdx
	am8	%rdi, 11, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %rax, %r14, %r15

	mov	%r9, 12*8(%rdi)
	mov	%r10, 13*8(%rdi)
	mov	%r11, 14*8(%rdi)
	mov	%rbx, 15*8(%rdi)
	mov	%rbp, 16*8(%rdi)
	mov	%r12, 17*8(%rdi)
	mov	%r13, 18*8(%rdi)
	mov	%rax, 19*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_13_1)
	mov	0*8(%rdx), %rdx
	m13_str	%rdi, 0, %rsi, 0, %rcx, %r8, %r9, %r10, %r11, %rax
	mov	%rax, 13*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_13_2)
	mov	0*8(%rdx), %rcx
	mov	1*8(%rdx), %r8
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	mov	%rcx, %rdx
	xor	%r13d, %r13d
	m5	%rdi, 0, %rsi, 0, %r9, %r10, %r11, %rbx, %rbp, %rax, %r12, %r13
	mov	%r8, %rdx
	am5	%rdi, 1, %rsi, 0, %r9, %r10, %r11, %rbx, %rbp, %rax, %r12, %r13
	mov	%r10, 2*8(%rdi)
	mov	%r11, 3*8(%rdi)
	mov	%rbx, 4*8(%rdi)

	mov	%rcx, %rdx
	m5_chain	%rdi, 5, %rsi, 5, %rbp, %rax, %r9, %r10, %r11, %rbx, %rax, %r12, %rbp, %r13
	mov	%r8, %rdx
	am5	%rdi, 6, %rsi, 5, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12, %r13
	mov	%r10, 7*8(%rdi)
	mov	%r11, 8*8(%rdi)
	mov	%rbx, 9*8(%rdi)

	mov	%rcx, %rdx
	m3_chain	%rdi, 10, %rsi, 10, %rax, %rbp, %r9, %r10, %r11, %r12, %rax, %r13
	mov	%r8, %rdx
	am3	%rdi, 11, %rsi, 10, %r9, %r10, %r11, %rax, %r12, %r13
	mov	%r10, 12*8(%rdi)
	mov	%r11, 13*8(%rdi)
	mov	%rax, 14*8(%rdi)
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_13_3)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx

	xor	%ebx, %ebx

	m3	%rdi, 0, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx

	mov	1*8(%rsi), %rdx
	am3	%rdi, 1, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx
	mov	2*8(%rsi), %rdx
	am3	%rdi, 2, %rcx, 0, %r8, %rax, %r10, %r9, %r11, %rbx
	mov	3*8(%rsi), %rdx
	am3	%rdi, 3, %rcx, 0, %rax, %r10, %r9, %r8, %r11, %rbx
	mov	4*8(%rsi), %rdx
	am3	%rdi, 4, %rcx, 0, %r10, %r9, %r8, %rax, %r11, %rbx
	mov	5*8(%rsi), %rdx
	am3	%rdi, 5, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx
	mov	6*8(%rsi), %rdx
	am3	%rdi, 6, %rcx, 0, %r8, %rax, %r10, %r9, %r11, %rbx
	mov	7*8(%rsi), %rdx
	am3	%rdi, 7, %rcx, 0, %rax, %r10, %r9, %r8, %r11, %rbx
	mov	8*8(%rsi), %rdx
	am3	%rdi, 8, %rcx, 0, %r10, %r9, %r8, %rax, %r11, %rbx
	mov	9*8(%rsi), %rdx
	am3	%rdi, 9, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx
	mov	10*8(%rsi), %rdx
	am3	%rdi, 10, %rcx, 0, %r8, %rax, %r10, %r9, %r11, %rbx
	mov	11*8(%rsi), %rdx
	am3	%rdi, 11, %rcx, 0, %rax, %r10, %r9, %r8, %r11, %rbx
	mov	12*8(%rsi), %rdx
	am3	%rdi, 12, %rcx, 0, %r10, %r9, %r8, %rax, %r11, %rbx

	mov	%r9, 13*8(%rdi)
	mov	%r8, 14*8(%rdi)
	mov	%rax, 15*8(%rdi)

	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_13_4)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp

	xor	%ebp, %ebp

	m4	%rdi, 0, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp

	mov	1*8(%rsi), %rdx
	am4	%rdi, 1, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp
	mov	2*8(%rsi), %rdx
	am4	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rax, %rbx, %rbp
	mov	3*8(%rsi), %rdx
	am4	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rax, %r8, %rbx, %rbp
	mov	4*8(%rsi), %rdx
	am4	%rdi, 4, %rcx, 0, %r10, %r11, %rax, %r8, %r9, %rbx, %rbp
	mov	5*8(%rsi), %rdx
	am4	%rdi, 5, %rcx, 0, %r11, %rax, %r8, %r9, %r10, %rbx, %rbp
	mov	6*8(%rsi), %rdx
	am4	%rdi, 6, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp
	mov	7*8(%rsi), %rdx
	am4	%rdi, 7, %rcx, 0, %r8, %r9, %r10, %r11, %rax, %rbx, %rbp
	mov	8*8(%rsi), %rdx
	am4	%rdi, 8, %rcx, 0, %r9, %r10, %r11, %rax, %r8, %rbx, %rbp
	mov	9*8(%rsi), %rdx
	am4	%rdi, 9, %rcx, 0, %r10, %r11, %rax, %r8, %r9, %rbx, %rbp
	mov	10*8(%rsi), %rdx
	am4	%rdi, 10, %rcx, 0, %r11, %rax, %r8, %r9, %r10, %rbx, %rbp
	mov	11*8(%rsi), %rdx
	am4	%rdi, 11, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp
	mov	12*8(%rsi), %rdx
	am4	%rdi, 12, %rcx, 0, %r8, %r9, %r10, %r11, %rax, %rbx, %rbp

	mov	%r9, 13*8(%rdi)
	mov	%r10, 14*8(%rdi)
	mov	%r11, 15*8(%rdi)
	mov	%rax, 16*8(%rdi)

	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_13_5)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12

	xor	%r12d, %r12d

	m5	%rdi, 0, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12

	mov	1*8(%rsi), %rdx
	am5	%rdi, 1, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12
	mov	2*8(%rsi), %rdx
	am5	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %rax, %rbx, %r11, %rbp, %r12
	mov	3*8(%rsi), %rdx
	am5	%rdi, 3, %rcx, 0, %r9, %r10, %rax, %rbx, %r11, %r8, %rbp, %r12
	mov	4*8(%rsi), %rdx
	am5	%rdi, 4, %rcx, 0, %r10, %rax, %rbx, %r11, %r8, %r9, %rbp, %r12
	mov	5*8(%rsi), %rdx
	am5	%rdi, 5, %rcx, 0, %rax, %rbx, %r11, %r8, %r9, %r10, %rbp, %r12
	mov	6*8(%rsi), %rdx
	am5	%rdi, 6, %rcx, 0, %rbx, %r11, %r8, %r9, %r10, %rax, %rbp, %r12
	mov	7*8(%rsi), %rdx
	am5	%rdi, 7, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12
	mov	8*8(%rsi), %rdx
	am5	%rdi, 8, %rcx, 0, %r8, %r9, %r10, %rax, %rbx, %r11, %rbp, %r12
	mov	9*8(%rsi), %rdx
	am5	%rdi, 9, %rcx, 0, %r9, %r10, %rax, %rbx, %r11, %r8, %rbp, %r12
	mov	10*8(%rsi), %rdx
	am5	%rdi, 10, %rcx, 0, %r10, %rax, %rbx, %r11, %r8, %r9, %rbp, %r12
	mov	11*8(%rsi), %rdx
	am5	%rdi, 11, %rcx, 0, %rax, %rbx, %r11, %r8, %r9, %r10, %rbp, %r12
	mov	12*8(%rsi), %rdx
	am5	%rdi, 12, %rcx, 0, %rbx, %r11, %r8, %r9, %r10, %rax, %rbp, %r12

	mov	%r11, 13*8(%rdi)
	mov	%r8, 14*8(%rdi)
	mov	%r9, 15*8(%rdi)
	mov	%r10, 16*8(%rdi)
	mov	%rax, 17*8(%rdi)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_13_6)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	xor	%r13d, %r13d

	m6	%rdi, 0, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13

	mov	1*8(%rsi), %rdx
	am6	%rdi, 1, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13
	mov	2*8(%rsi), %rdx
	am6	%rdi, 2, %rcx, 0, %r8, %r9, %rax, %r11, %rbx, %rbp, %r10, %r12, %r13
	mov	3*8(%rsi), %rdx
	am6	%rdi, 3, %rcx, 0, %r9, %rax, %r11, %rbx, %rbp, %r10, %r8, %r12, %r13
	mov	4*8(%rsi), %rdx
	am6	%rdi, 4, %rcx, 0, %rax, %r11, %rbx, %rbp, %r10, %r8, %r9, %r12, %r13
	mov	5*8(%rsi), %rdx
	am6	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %r10, %r8, %r9, %rax, %r12, %r13
	mov	6*8(%rsi), %rdx
	am6	%rdi, 6, %rcx, 0, %rbx, %rbp, %r10, %r8, %r9, %rax, %r11, %r12, %r13
	mov	7*8(%rsi), %rdx
	am6	%rdi, 7, %rcx, 0, %rbp, %r10, %r8, %r9, %rax, %r11, %rbx, %r12, %r13
	mov	8*8(%rsi), %rdx
	am6	%rdi, 8, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13
	mov	9*8(%rsi), %rdx
	am6	%rdi, 9, %rcx, 0, %r8, %r9, %rax, %r11, %rbx, %rbp, %r10, %r12, %r13
	mov	10*8(%rsi), %rdx
	am6	%rdi, 10, %rcx, 0, %r9, %rax, %r11, %rbx, %rbp, %r10, %r8, %r12, %r13
	mov	11*8(%rsi), %rdx
	am6	%rdi, 11, %rcx, 0, %rax, %r11, %rbx, %rbp, %r10, %r8, %r9, %r12, %r13
	mov	12*8(%rsi), %rdx
	am6	%rdi, 12, %rcx, 0, %r11, %rbx, %rbp, %r10, %r8, %r9, %rax, %r12, %r13

	mov	%rbx, 13*8(%rdi)
	mov	%rbp, 14*8(%rdi)
	mov	%r10, 15*8(%rdi)
	mov	%r8, 16*8(%rdi)
	mov	%r9, 17*8(%rdi)
	mov	%rax, 18*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_13_7)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	%r14d, %r14d

	m7	%rdi, 0, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14

	mov	1*8(%rsi), %rdx
	am7	%rdi, 1, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14
	mov	2*8(%rsi), %rdx
	am7	%rdi, 2, %rcx, 0, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r9, %r13, %r14
	mov	3*8(%rsi), %rdx
	am7	%rdi, 3, %rcx, 0, %rax, %r10, %r11, %rbx, %rbp, %r12, %r9, %r8, %r13, %r14
	mov	4*8(%rsi), %rdx
	am7	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rbp, %r12, %r9, %r8, %rax, %r13, %r14
	mov	5*8(%rsi), %rdx
	am7	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %r12, %r9, %r8, %rax, %r10, %r13, %r14
	mov	6*8(%rsi), %rdx
	am7	%rdi, 6, %rcx, 0, %rbx, %rbp, %r12, %r9, %r8, %rax, %r10, %r11, %r13, %r14
	mov	7*8(%rsi), %rdx
	am7	%rdi, 7, %rcx, 0, %rbp, %r12, %r9, %r8, %rax, %r10, %r11, %rbx, %r13, %r14
	mov	8*8(%rsi), %rdx
	am7	%rdi, 8, %rcx, 0, %r12, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r13, %r14
	mov	9*8(%rsi), %rdx
	am7	%rdi, 9, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14
	mov	10*8(%rsi), %rdx
	am7	%rdi, 10, %rcx, 0, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r9, %r13, %r14
	mov	11*8(%rsi), %rdx
	am7	%rdi, 11, %rcx, 0, %rax, %r10, %r11, %rbx, %rbp, %r12, %r9, %r8, %r13, %r14
	mov	12*8(%rsi), %rdx
	am7	%rdi, 12, %rcx, 0, %r10, %r11, %rbx, %rbp, %r12, %r9, %r8, %rax, %r13, %r14

	mov	%r11, 13*8(%rdi)
	mov	%rbx, 14*8(%rdi)
	mov	%rbp, 15*8(%rdi)
	mov	%r12, 16*8(%rdi)
	mov	%r9, 17*8(%rdi)
	mov	%r8, 18*8(%rdi)
	mov	%rax, 19*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_13_8)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	xor	%r15d, %r15d

	m8	%rdi, 0, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15

	mov	1*8(%rsi), %rdx
	am8	%rdi, 1, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15
	mov	2*8(%rsi), %rdx
	am8	%rdi, 2, %rcx, 0, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r8, %r14, %r15
	mov	3*8(%rsi), %rdx
	am8	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r8, %rax, %r14, %r15
	mov	4*8(%rsi), %rdx
	am8	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rbp, %r12, %r13, %r8, %rax, %r9, %r14, %r15
	mov	5*8(%rsi), %rdx
	am8	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %r12, %r13, %r8, %rax, %r9, %r10, %r14, %r15
	mov	6*8(%rsi), %rdx
	am8	%rdi, 6, %rcx, 0, %rbx, %rbp, %r12, %r13, %r8, %rax, %r9, %r10, %r11, %r14, %r15
	mov	7*8(%rsi), %rdx
	am8	%rdi, 7, %rcx, 0, %rbp, %r12, %r13, %r8, %rax, %r9, %r10, %r11, %rbx, %r14, %r15
	mov	8*8(%rsi), %rdx
	am8	%rdi, 8, %rcx, 0, %r12, %r13, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r14, %r15
	mov	9*8(%rsi), %rdx
	am8	%rdi, 9, %rcx, 0, %r13, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r14, %r15
	mov	10*8(%rsi), %rdx
	am8	%rdi, 10, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15
	mov	11*8(%rsi), %rdx
	am8	%rdi, 11, %rcx, 0, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r8, %r14, %r15
	mov	12*8(%rsi), %rdx
	am8	%rdi, 12, %rcx, 0, %r9, %r10, %r11, %rbx, %rbp, %r12, %r13, %r8, %rax, %r14, %r15

	mov	%r10, 13*8(%rdi)
	mov	%r11, 14*8(%rdi)
	mov	%rbx, 15*8(%rdi)
	mov	%rbp, 16*8(%rdi)
	mov	%r12, 17*8(%rdi)
	mov	%r13, 18*8(%rdi)
	mov	%r8, 19*8(%rdi)
	mov	%rax, 20*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_14_1)
	mov	0*8(%rdx), %rdx
	m14_str	%rdi, 0, %rsi, 0, %rcx, %r8, %r9, %r10, %r11, %rax
	mov	%rax, 14*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_14_2)
	mov	0*8(%rdx), %rcx
	mov	1*8(%rdx), %r8
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	mov	%rcx, %rdx
	xor	%r13d, %r13d
	m5	%rdi, 0, %rsi, 0, %r9, %r10, %r11, %rbx, %rbp, %rax, %r12, %r13
	mov	%r8, %rdx
	am5	%rdi, 1, %rsi, 0, %r9, %r10, %r11, %rbx, %rbp, %rax, %r12, %r13
	mov	%r10, 2*8(%rdi)
	mov	%r11, 3*8(%rdi)
	mov	%rbx, 4*8(%rdi)

	mov	%rcx, %rdx
	m5_chain	%rdi, 5, %rsi, 5, %rbp, %rax, %r9, %r10, %r11, %rbx, %rax, %r12, %rbp, %r13
	mov	%r8, %rdx
	am5	%rdi, 6, %rsi, 5, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12, %r13
	mov	%r10, 7*8(%rdi)
	mov	%r11, 8*8(%rdi)
	mov	%rbx, 9*8(%rdi)

	mov	%rcx, %rdx
	m4_chain	%rdi, 10, %rsi, 10, %rax, %rbp, %r9, %r10, %r11, %rbx, %r12, %rax, %r13
	mov	%r8, %rdx
	am4	%rdi, 11, %rsi, 10, %r9, %r10, %r11, %rbx, %rax, %r12, %r13
	mov	%r10, 12*8(%rdi)
	mov	%r11, 13*8(%rdi)
	mov	%rbx, 14*8(%rdi)
	mov	%rax, 15*8(%rdi)
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_14_3)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx

	xor	%ebx, %ebx

	m3	%rdi, 0, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx

	mov	1*8(%rsi), %rdx
	am3	%rdi, 1, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx
	mov	2*8(%rsi), %rdx
	am3	%rdi, 2, %rcx, 0, %r8, %r9, %rax, %r10, %r11, %rbx
	mov	3*8(%rsi), %rdx
	am3	%rdi, 3, %rcx, 0, %r9, %rax, %r10, %r8, %r11, %rbx
	mov	4*8(%rsi), %rdx
	am3	%rdi, 4, %rcx, 0, %rax, %r10, %r8, %r9, %r11, %rbx
	mov	5*8(%rsi), %rdx
	am3	%rdi, 5, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx
	mov	6*8(%rsi), %rdx
	am3	%rdi, 6, %rcx, 0, %r8, %r9, %rax, %r10, %r11, %rbx
	mov	7*8(%rsi), %rdx
	am3	%rdi, 7, %rcx, 0, %r9, %rax, %r10, %r8, %r11, %rbx
	mov	8*8(%rsi), %rdx
	am3	%rdi, 8, %rcx, 0, %rax, %r10, %r8, %r9, %r11, %rbx
	mov	9*8(%rsi), %rdx
	am3	%rdi, 9, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx
	mov	10*8(%rsi), %rdx
	am3	%rdi, 10, %rcx, 0, %r8, %r9, %rax, %r10, %r11, %rbx
	mov	11*8(%rsi), %rdx
	am3	%rdi, 11, %rcx, 0, %r9, %rax, %r10, %r8, %r11, %rbx
	mov	12*8(%rsi), %rdx
	am3	%rdi, 12, %rcx, 0, %rax, %r10, %r8, %r9, %r11, %rbx
	mov	13*8(%rsi), %rdx
	am3	%rdi, 13, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx

	mov	%r8, 14*8(%rdi)
	mov	%r9, 15*8(%rdi)
	mov	%rax, 16*8(%rdi)

	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_14_4)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp

	xor	%ebp, %ebp

	m4	%rdi, 0, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp

	mov	1*8(%rsi), %rdx
	am4	%rdi, 1, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp
	mov	2*8(%rsi), %rdx
	am4	%rdi, 2, %rcx, 0, %rax, %r9, %r10, %r11, %r8, %rbx, %rbp
	mov	3*8(%rsi), %rdx
	am4	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %r8, %rax, %rbx, %rbp
	mov	4*8(%rsi), %rdx
	am4	%rdi, 4, %rcx, 0, %r10, %r11, %r8, %rax, %r9, %rbx, %rbp
	mov	5*8(%rsi), %rdx
	am4	%rdi, 5, %rcx, 0, %r11, %r8, %rax, %r9, %r10, %rbx, %rbp
	mov	6*8(%rsi), %rdx
	am4	%rdi, 6, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp
	mov	7*8(%rsi), %rdx
	am4	%rdi, 7, %rcx, 0, %rax, %r9, %r10, %r11, %r8, %rbx, %rbp
	mov	8*8(%rsi), %rdx
	am4	%rdi, 8, %rcx, 0, %r9, %r10, %r11, %r8, %rax, %rbx, %rbp
	mov	9*8(%rsi), %rdx
	am4	%rdi, 9, %rcx, 0, %r10, %r11, %r8, %rax, %r9, %rbx, %rbp
	mov	10*8(%rsi), %rdx
	am4	%rdi, 10, %rcx, 0, %r11, %r8, %rax, %r9, %r10, %rbx, %rbp
	mov	11*8(%rsi), %rdx
	am4	%rdi, 11, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp
	mov	12*8(%rsi), %rdx
	am4	%rdi, 12, %rcx, 0, %rax, %r9, %r10, %r11, %r8, %rbx, %rbp
	mov	13*8(%rsi), %rdx
	am4	%rdi, 13, %rcx, 0, %r9, %r10, %r11, %r8, %rax, %rbx, %rbp

	mov	%r10, 14*8(%rdi)
	mov	%r11, 15*8(%rdi)
	mov	%r8, 16*8(%rdi)
	mov	%rax, 17*8(%rdi)

	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_14_5)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12

	xor	%r12d, %r12d

	m5	%rdi, 0, %rcx, 0, %rbx, %r8, %r9, %r10, %r11, %rax, %rbp, %r12

	mov	1*8(%rsi), %rdx
	am5	%rdi, 1, %rcx, 0, %rbx, %r8, %r9, %r10, %r11, %rax, %rbp, %r12
	mov	2*8(%rsi), %rdx
	am5	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rax, %rbx, %rbp, %r12
	mov	3*8(%rsi), %rdx
	am5	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rax, %rbx, %r8, %rbp, %r12
	mov	4*8(%rsi), %rdx
	am5	%rdi, 4, %rcx, 0, %r10, %r11, %rax, %rbx, %r8, %r9, %rbp, %r12
	mov	5*8(%rsi), %rdx
	am5	%rdi, 5, %rcx, 0, %r11, %rax, %rbx, %r8, %r9, %r10, %rbp, %r12
	mov	6*8(%rsi), %rdx
	am5	%rdi, 6, %rcx, 0, %rax, %rbx, %r8, %r9, %r10, %r11, %rbp, %r12
	mov	7*8(%rsi), %rdx
	am5	%rdi, 7, %rcx, 0, %rbx, %r8, %r9, %r10, %r11, %rax, %rbp, %r12
	mov	8*8(%rsi), %rdx
	am5	%rdi, 8, %rcx, 0, %r8, %r9, %r10, %r11, %rax, %rbx, %rbp, %r12
	mov	9*8(%rsi), %rdx
	am5	%rdi, 9, %rcx, 0, %r9, %r10, %r11, %rax, %rbx, %r8, %rbp, %r12
	mov	10*8(%rsi), %rdx
	am5	%rdi, 10, %rcx, 0, %r10, %r11, %rax, %rbx, %r8, %r9, %rbp, %r12
	mov	11*8(%rsi), %rdx
	am5	%rdi, 11, %rcx, 0, %r11, %rax, %rbx, %r8, %r9, %r10, %rbp, %r12
	mov	12*8(%rsi), %rdx
	am5	%rdi, 12, %rcx, 0, %rax, %rbx, %r8, %r9, %r10, %r11, %rbp, %r12
	mov	13*8(%rsi), %rdx
	am5	%rdi, 13, %rcx, 0, %rbx, %r8, %r9, %r10, %r11, %rax, %rbp, %r12

	mov	%r8, 14*8(%rdi)
	mov	%r9, 15*8(%rdi)
	mov	%r10, 16*8(%rdi)
	mov	%r11, 17*8(%rdi)
	mov	%rax, 18*8(%rdi)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_14_6)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	xor	%r13d, %r13d

	m6	%rdi, 0, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13

	mov	1*8(%rsi), %rdx
	am6	%rdi, 1, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13
	mov	2*8(%rsi), %rdx
	am6	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %rax, %rbx, %rbp, %r11, %r12, %r13
	mov	3*8(%rsi), %rdx
	am6	%rdi, 3, %rcx, 0, %r9, %r10, %rax, %rbx, %rbp, %r11, %r8, %r12, %r13
	mov	4*8(%rsi), %rdx
	am6	%rdi, 4, %rcx, 0, %r10, %rax, %rbx, %rbp, %r11, %r8, %r9, %r12, %r13
	mov	5*8(%rsi), %rdx
	am6	%rdi, 5, %rcx, 0, %rax, %rbx, %rbp, %r11, %r8, %r9, %r10, %r12, %r13
	mov	6*8(%rsi), %rdx
	am6	%rdi, 6, %rcx, 0, %rbx, %rbp, %r11, %r8, %r9, %r10, %rax, %r12, %r13
	mov	7*8(%rsi), %rdx
	am6	%rdi, 7, %rcx, 0, %rbp, %r11, %r8, %r9, %r10, %rax, %rbx, %r12, %r13
	mov	8*8(%rsi), %rdx
	am6	%rdi, 8, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13
	mov	9*8(%rsi), %rdx
	am6	%rdi, 9, %rcx, 0, %r8, %r9, %r10, %rax, %rbx, %rbp, %r11, %r12, %r13
	mov	10*8(%rsi), %rdx
	am6	%rdi, 10, %rcx, 0, %r9, %r10, %rax, %rbx, %rbp, %r11, %r8, %r12, %r13
	mov	11*8(%rsi), %rdx
	am6	%rdi, 11, %rcx, 0, %r10, %rax, %rbx, %rbp, %r11, %r8, %r9, %r12, %r13
	mov	12*8(%rsi), %rdx
	am6	%rdi, 12, %rcx, 0, %rax, %rbx, %rbp, %r11, %r8, %r9, %r10, %r12, %r13
	mov	13*8(%rsi), %rdx
	am6	%rdi, 13, %rcx, 0, %rbx, %rbp, %r11, %r8, %r9, %r10, %rax, %r12, %r13

	mov	%rbp, 14*8(%rdi)
	mov	%r11, 15*8(%rdi)
	mov	%r8, 16*8(%rdi)
	mov	%r9, 17*8(%rdi)
	mov	%r10, 18*8(%rdi)
	mov	%rax, 19*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_14_7)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	%r14d, %r14d

	m7	%rdi, 0, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r14

	mov	1*8(%rsi), %rdx
	am7	%rdi, 1, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r14
	mov	2*8(%rsi), %rdx
	am7	%rdi, 2, %rcx, 0, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r10, %r13, %r14
	mov	3*8(%rsi), %rdx
	am7	%rdi, 3, %rcx, 0, %r9, %rax, %r11, %rbx, %rbp, %r12, %r10, %r8, %r13, %r14
	mov	4*8(%rsi), %rdx
	am7	%rdi, 4, %rcx, 0, %rax, %r11, %rbx, %rbp, %r12, %r10, %r8, %r9, %r13, %r14
	mov	5*8(%rsi), %rdx
	am7	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %r12, %r10, %r8, %r9, %rax, %r13, %r14
	mov	6*8(%rsi), %rdx
	am7	%rdi, 6, %rcx, 0, %rbx, %rbp, %r12, %r10, %r8, %r9, %rax, %r11, %r13, %r14
	mov	7*8(%rsi), %rdx
	am7	%rdi, 7, %rcx, 0, %rbp, %r12, %r10, %r8, %r9, %rax, %r11, %rbx, %r13, %r14
	mov	8*8(%rsi), %rdx
	am7	%rdi, 8, %rcx, 0, %r12, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r13, %r14
	mov	9*8(%rsi), %rdx
	am7	%rdi, 9, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r14
	mov	10*8(%rsi), %rdx
	am7	%rdi, 10, %rcx, 0, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r10, %r13, %r14
	mov	11*8(%rsi), %rdx
	am7	%rdi, 11, %rcx, 0, %r9, %rax, %r11, %rbx, %rbp, %r12, %r10, %r8, %r13, %r14
	mov	12*8(%rsi), %rdx
	am7	%rdi, 12, %rcx, 0, %rax, %r11, %rbx, %rbp, %r12, %r10, %r8, %r9, %r13, %r14
	mov	13*8(%rsi), %rdx
	am7	%rdi, 13, %rcx, 0, %r11, %rbx, %rbp, %r12, %r10, %r8, %r9, %rax, %r13, %r14

	mov	%rbx, 14*8(%rdi)
	mov	%rbp, 15*8(%rdi)
	mov	%r12, 16*8(%rdi)
	mov	%r10, 17*8(%rdi)
	mov	%r8, 18*8(%rdi)
	mov	%r9, 19*8(%rdi)
	mov	%rax, 20*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_14_8)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	xor	%r15d, %r15d

	m8	%rdi, 0, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15

	mov	1*8(%rsi), %rdx
	am8	%rdi, 1, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15
	mov	2*8(%rsi), %rdx
	am8	%rdi, 2, %rcx, 0, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r9, %r14, %r15
	mov	3*8(%rsi), %rdx
	am8	%rdi, 3, %rcx, 0, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r9, %r8, %r14, %r15
	mov	4*8(%rsi), %rdx
	am8	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rbp, %r12, %r13, %r9, %r8, %rax, %r14, %r15
	mov	5*8(%rsi), %rdx
	am8	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %r12, %r13, %r9, %r8, %rax, %r10, %r14, %r15
	mov	6*8(%rsi), %rdx
	am8	%rdi, 6, %rcx, 0, %rbx, %rbp, %r12, %r13, %r9, %r8, %rax, %r10, %r11, %r14, %r15
	mov	7*8(%rsi), %rdx
	am8	%rdi, 7, %rcx, 0, %rbp, %r12, %r13, %r9, %r8, %rax, %r10, %r11, %rbx, %r14, %r15
	mov	8*8(%rsi), %rdx
	am8	%rdi, 8, %rcx, 0, %r12, %r13, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r14, %r15
	mov	9*8(%rsi), %rdx
	am8	%rdi, 9, %rcx, 0, %r13, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r14, %r15
	mov	10*8(%rsi), %rdx
	am8	%rdi, 10, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15
	mov	11*8(%rsi), %rdx
	am8	%rdi, 11, %rcx, 0, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r9, %r14, %r15
	mov	12*8(%rsi), %rdx
	am8	%rdi, 12, %rcx, 0, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13, %r9, %r8, %r14, %r15
	mov	13*8(%rsi), %rdx
	am8	%rdi, 13, %rcx, 0, %r10, %r11, %rbx, %rbp, %r12, %r13, %r9, %r8, %rax, %r14, %r15

	mov	%r11, 14*8(%rdi)
	mov	%rbx, 15*8(%rdi)
	mov	%rbp, 16*8(%rdi)
	mov	%r12, 17*8(%rdi)
	mov	%r13, 18*8(%rdi)
	mov	%r9, 19*8(%rdi)
	mov	%r8, 20*8(%rdi)
	mov	%rax, 21*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_15_1)
	mov	0*8(%rdx), %rdx
	m15_str	%rdi, 0, %rsi, 0, %rcx, %r8, %r9, %r10, %r11, %rax
	mov	%rax, 15*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_15_2)
	mov	0*8(%rdx), %rcx
	mov	1*8(%rdx), %r8
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	mov	%rcx, %rdx
	xor	%r13d, %r13d
	m5	%rdi, 0, %rsi, 0, %r9, %r10, %r11, %rbx, %rbp, %rax, %r12, %r13
	mov	%r8, %rdx
	am5	%rdi, 1, %rsi, 0, %r9, %r10, %r11, %rbx, %rbp, %rax, %r12, %r13
	mov	%r10, 2*8(%rdi)
	mov	%r11, 3*8(%rdi)
	mov	%rbx, 4*8(%rdi)

	mov	%rcx, %rdx
	m5_chain	%rdi, 5, %rsi, 5, %rbp, %rax, %r9, %r10, %r11, %rbx, %rax, %r12, %rbp, %r13
	mov	%r8, %rdx
	am5	%rdi, 6, %rsi, 5, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12, %r13
	mov	%r10, 7*8(%rdi)
	mov	%r11, 8*8(%rdi)
	mov	%rbx, 9*8(%rdi)

	mov	%rcx, %rdx
	m5_chain	%rdi, 10, %rsi, 10, %rax, %rbp, %r9, %r10, %r11, %rbx, %rbp, %r12, %rax, %r13
	mov	%r8, %rdx
	am5	%rdi, 11, %rsi, 10, %r9, %r10, %r11, %rbx, %rbp, %rax, %r12, %r13
	mov	%r10, 12*8(%rdi)
	mov	%r11, 13*8(%rdi)
	mov	%rbx, 14*8(%rdi)

	mov	%rbp, 15*8(%rdi)
	mov	%rax, 16*8(%rdi)
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_15_3)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx

	xor	%ebx, %ebx

	m3	%rdi, 0, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx

	mov	1*8(%rsi), %rdx
	am3	%rdi, 1, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx
	mov	2*8(%rsi), %rdx
	am3	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %rax, %r11, %rbx
	mov	3*8(%rsi), %rdx
	am3	%rdi, 3, %rcx, 0, %r9, %r10, %rax, %r8, %r11, %rbx
	mov	4*8(%rsi), %rdx
	am3	%rdi, 4, %rcx, 0, %r10, %rax, %r8, %r9, %r11, %rbx
	mov	5*8(%rsi), %rdx
	am3	%rdi, 5, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx
	mov	6*8(%rsi), %rdx
	am3	%rdi, 6, %rcx, 0, %r8, %r9, %r10, %rax, %r11, %rbx
	mov	7*8(%rsi), %rdx
	am3	%rdi, 7, %rcx, 0, %r9, %r10, %rax, %r8, %r11, %rbx
	mov	8*8(%rsi), %rdx
	am3	%rdi, 8, %rcx, 0, %r10, %rax, %r8, %r9, %r11, %rbx
	mov	9*8(%rsi), %rdx
	am3	%rdi, 9, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx
	mov	10*8(%rsi), %rdx
	am3	%rdi, 10, %rcx, 0, %r8, %r9, %r10, %rax, %r11, %rbx
	mov	11*8(%rsi), %rdx
	am3	%rdi, 11, %rcx, 0, %r9, %r10, %rax, %r8, %r11, %rbx
	mov	12*8(%rsi), %rdx
	am3	%rdi, 12, %rcx, 0, %r10, %rax, %r8, %r9, %r11, %rbx
	mov	13*8(%rsi), %rdx
	am3	%rdi, 13, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx
	mov	14*8(%rsi), %rdx
	am3	%rdi, 14, %rcx, 0, %r8, %r9, %r10, %rax, %r11, %rbx

	mov	%r9, 15*8(%rdi)
	mov	%r10, 16*8(%rdi)
	mov	%rax, 17*8(%rdi)

	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_15_4)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp

	xor	%ebp, %ebp

	m4	%rdi, 0, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp

	mov	1*8(%rsi), %rdx
	am4	%rdi, 1, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp
	mov	2*8(%rsi), %rdx
	am4	%rdi, 2, %rcx, 0, %r8, %rax, %r10, %r11, %r9, %rbx, %rbp
	mov	3*8(%rsi), %rdx
	am4	%rdi, 3, %rcx, 0, %rax, %r10, %r11, %r9, %r8, %rbx, %rbp
	mov	4*8(%rsi), %rdx
	am4	%rdi, 4, %rcx, 0, %r10, %r11, %r9, %r8, %rax, %rbx, %rbp
	mov	5*8(%rsi), %rdx
	am4	%rdi, 5, %rcx, 0, %r11, %r9, %r8, %rax, %r10, %rbx, %rbp
	mov	6*8(%rsi), %rdx
	am4	%rdi, 6, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp
	mov	7*8(%rsi), %rdx
	am4	%rdi, 7, %rcx, 0, %r8, %rax, %r10, %r11, %r9, %rbx, %rbp
	mov	8*8(%rsi), %rdx
	am4	%rdi, 8, %rcx, 0, %rax, %r10, %r11, %r9, %r8, %rbx, %rbp
	mov	9*8(%rsi), %rdx
	am4	%rdi, 9, %rcx, 0, %r10, %r11, %r9, %r8, %rax, %rbx, %rbp
	mov	10*8(%rsi), %rdx
	am4	%rdi, 10, %rcx, 0, %r11, %r9, %r8, %rax, %r10, %rbx, %rbp
	mov	11*8(%rsi), %rdx
	am4	%rdi, 11, %rcx, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp
	mov	12*8(%rsi), %rdx
	am4	%rdi, 12, %rcx, 0, %r8, %rax, %r10, %r11, %r9, %rbx, %rbp
	mov	13*8(%rsi), %rdx
	am4	%rdi, 13, %rcx, 0, %rax, %r10, %r11, %r9, %r8, %rbx, %rbp
	mov	14*8(%rsi), %rdx
	am4	%rdi, 14, %rcx, 0, %r10, %r11, %r9, %r8, %rax, %rbx, %rbp

	mov	%r11, 15*8(%rdi)
	mov	%r9, 16*8(%rdi)
	mov	%r8, 17*8(%rdi)
	mov	%rax, 18*8(%rdi)

	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_15_5)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12

	xor	%r12d, %r12d

	m5	%rdi, 0, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12

	mov	1*8(%rsi), %rdx
	am5	%rdi, 1, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12
	mov	2*8(%rsi), %rdx
	am5	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12
	mov	3*8(%rsi), %rdx
	am5	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rax, %r8, %rbp, %r12
	mov	4*8(%rsi), %rdx
	am5	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rax, %r8, %r9, %rbp, %r12
	mov	5*8(%rsi), %rdx
	am5	%rdi, 5, %rcx, 0, %r11, %rbx, %rax, %r8, %r9, %r10, %rbp, %r12
	mov	6*8(%rsi), %rdx
	am5	%rdi, 6, %rcx, 0, %rbx, %rax, %r8, %r9, %r10, %r11, %rbp, %r12
	mov	7*8(%rsi), %rdx
	am5	%rdi, 7, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12
	mov	8*8(%rsi), %rdx
	am5	%rdi, 8, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12
	mov	9*8(%rsi), %rdx
	am5	%rdi, 9, %rcx, 0, %r9, %r10, %r11, %rbx, %rax, %r8, %rbp, %r12
	mov	10*8(%rsi), %rdx
	am5	%rdi, 10, %rcx, 0, %r10, %r11, %rbx, %rax, %r8, %r9, %rbp, %r12
	mov	11*8(%rsi), %rdx
	am5	%rdi, 11, %rcx, 0, %r11, %rbx, %rax, %r8, %r9, %r10, %rbp, %r12
	mov	12*8(%rsi), %rdx
	am5	%rdi, 12, %rcx, 0, %rbx, %rax, %r8, %r9, %r10, %r11, %rbp, %r12
	mov	13*8(%rsi), %rdx
	am5	%rdi, 13, %rcx, 0, %rax, %r8, %r9, %r10, %r11, %rbx, %rbp, %r12
	mov	14*8(%rsi), %rdx
	am5	%rdi, 14, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12

	mov	%r9, 15*8(%rdi)
	mov	%r10, 16*8(%rdi)
	mov	%r11, 17*8(%rdi)
	mov	%rbx, 18*8(%rdi)
	mov	%rax, 19*8(%rdi)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_15_6)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	xor	%r13d, %r13d

	m6	%rdi, 0, %rcx, 0, %rbx, %r8, %r9, %r10, %r11, %rax, %rbp, %r12, %r13

	mov	1*8(%rsi), %rdx
	am6	%rdi, 1, %rcx, 0, %rbx, %r8, %r9, %r10, %r11, %rax, %rbp, %r12, %r13
	mov	2*8(%rsi), %rdx
	am6	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rax, %rbp, %rbx, %r12, %r13
	mov	3*8(%rsi), %rdx
	am6	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rax, %rbp, %rbx, %r8, %r12, %r13
	mov	4*8(%rsi), %rdx
	am6	%rdi, 4, %rcx, 0, %r10, %r11, %rax, %rbp, %rbx, %r8, %r9, %r12, %r13
	mov	5*8(%rsi), %rdx
	am6	%rdi, 5, %rcx, 0, %r11, %rax, %rbp, %rbx, %r8, %r9, %r10, %r12, %r13
	mov	6*8(%rsi), %rdx
	am6	%rdi, 6, %rcx, 0, %rax, %rbp, %rbx, %r8, %r9, %r10, %r11, %r12, %r13
	mov	7*8(%rsi), %rdx
	am6	%rdi, 7, %rcx, 0, %rbp, %rbx, %r8, %r9, %r10, %r11, %rax, %r12, %r13
	mov	8*8(%rsi), %rdx
	am6	%rdi, 8, %rcx, 0, %rbx, %r8, %r9, %r10, %r11, %rax, %rbp, %r12, %r13
	mov	9*8(%rsi), %rdx
	am6	%rdi, 9, %rcx, 0, %r8, %r9, %r10, %r11, %rax, %rbp, %rbx, %r12, %r13
	mov	10*8(%rsi), %rdx
	am6	%rdi, 10, %rcx, 0, %r9, %r10, %r11, %rax, %rbp, %rbx, %r8, %r12, %r13
	mov	11*8(%rsi), %rdx
	am6	%rdi, 11, %rcx, 0, %r10, %r11, %rax, %rbp, %rbx, %r8, %r9, %r12, %r13
	mov	12*8(%rsi), %rdx
	am6	%rdi, 12, %rcx, 0, %r11, %rax, %rbp, %rbx, %r8, %r9, %r10, %r12, %r13
	mov	13*8(%rsi), %rdx
	am6	%rdi, 13, %rcx, 0, %rax, %rbp, %rbx, %r8, %r9, %r10, %r11, %r12, %r13
	mov	14*8(%rsi), %rdx
	am6	%rdi, 14, %rcx, 0, %rbp, %rbx, %r8, %r9, %r10, %r11, %rax, %r12, %r13

	mov	%rbx, 15*8(%rdi)
	mov	%r8, 16*8(%rdi)
	mov	%r9, 17*8(%rdi)
	mov	%r10, 18*8(%rdi)
	mov	%r11, 19*8(%rdi)
	mov	%rax, 20*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_15_7)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	%r14d, %r14d

	m7	%rdi, 0, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r14

	mov	1*8(%rsi), %rdx
	am7	%rdi, 1, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r14
	mov	2*8(%rsi), %rdx
	am7	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r11, %r13, %r14
	mov	3*8(%rsi), %rdx
	am7	%rdi, 3, %rcx, 0, %r9, %r10, %rax, %rbx, %rbp, %r12, %r11, %r8, %r13, %r14
	mov	4*8(%rsi), %rdx
	am7	%rdi, 4, %rcx, 0, %r10, %rax, %rbx, %rbp, %r12, %r11, %r8, %r9, %r13, %r14
	mov	5*8(%rsi), %rdx
	am7	%rdi, 5, %rcx, 0, %rax, %rbx, %rbp, %r12, %r11, %r8, %r9, %r10, %r13, %r14
	mov	6*8(%rsi), %rdx
	am7	%rdi, 6, %rcx, 0, %rbx, %rbp, %r12, %r11, %r8, %r9, %r10, %rax, %r13, %r14
	mov	7*8(%rsi), %rdx
	am7	%rdi, 7, %rcx, 0, %rbp, %r12, %r11, %r8, %r9, %r10, %rax, %rbx, %r13, %r14
	mov	8*8(%rsi), %rdx
	am7	%rdi, 8, %rcx, 0, %r12, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r13, %r14
	mov	9*8(%rsi), %rdx
	am7	%rdi, 9, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r14
	mov	10*8(%rsi), %rdx
	am7	%rdi, 10, %rcx, 0, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r11, %r13, %r14
	mov	11*8(%rsi), %rdx
	am7	%rdi, 11, %rcx, 0, %r9, %r10, %rax, %rbx, %rbp, %r12, %r11, %r8, %r13, %r14
	mov	12*8(%rsi), %rdx
	am7	%rdi, 12, %rcx, 0, %r10, %rax, %rbx, %rbp, %r12, %r11, %r8, %r9, %r13, %r14
	mov	13*8(%rsi), %rdx
	am7	%rdi, 13, %rcx, 0, %rax, %rbx, %rbp, %r12, %r11, %r8, %r9, %r10, %r13, %r14
	mov	14*8(%rsi), %rdx
	am7	%rdi, 14, %rcx, 0, %rbx, %rbp, %r12, %r11, %r8, %r9, %r10, %rax, %r13, %r14

	mov	%rbp, 15*8(%rdi)
	mov	%r12, 16*8(%rdi)
	mov	%r11, 17*8(%rdi)
	mov	%r8, 18*8(%rdi)
	mov	%r9, 19*8(%rdi)
	mov	%r10, 20*8(%rdi)
	mov	%rax, 21*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_15_8)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	xor	%r15d, %r15d

	m8	%rdi, 0, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15

	mov	1*8(%rsi), %rdx
	am8	%rdi, 1, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15
	mov	2*8(%rsi), %rdx
	am8	%rdi, 2, %rcx, 0, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r10, %r14, %r15
	mov	3*8(%rsi), %rdx
	am8	%rdi, 3, %rcx, 0, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r10, %r8, %r14, %r15
	mov	4*8(%rsi), %rdx
	am8	%rdi, 4, %rcx, 0, %rax, %r11, %rbx, %rbp, %r12, %r13, %r10, %r8, %r9, %r14, %r15
	mov	5*8(%rsi), %rdx
	am8	%rdi, 5, %rcx, 0, %r11, %rbx, %rbp, %r12, %r13, %r10, %r8, %r9, %rax, %r14, %r15
	mov	6*8(%rsi), %rdx
	am8	%rdi, 6, %rcx, 0, %rbx, %rbp, %r12, %r13, %r10, %r8, %r9, %rax, %r11, %r14, %r15
	mov	7*8(%rsi), %rdx
	am8	%rdi, 7, %rcx, 0, %rbp, %r12, %r13, %r10, %r8, %r9, %rax, %r11, %rbx, %r14, %r15
	mov	8*8(%rsi), %rdx
	am8	%rdi, 8, %rcx, 0, %r12, %r13, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r14, %r15
	mov	9*8(%rsi), %rdx
	am8	%rdi, 9, %rcx, 0, %r13, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r14, %r15
	mov	10*8(%rsi), %rdx
	am8	%rdi, 10, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r14, %r15
	mov	11*8(%rsi), %rdx
	am8	%rdi, 11, %rcx, 0, %r8, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r10, %r14, %r15
	mov	12*8(%rsi), %rdx
	am8	%rdi, 12, %rcx, 0, %r9, %rax, %r11, %rbx, %rbp, %r12, %r13, %r10, %r8, %r14, %r15
	mov	13*8(%rsi), %rdx
	am8	%rdi, 13, %rcx, 0, %rax, %r11, %rbx, %rbp, %r12, %r13, %r10, %r8, %r9, %r14, %r15
	mov	14*8(%rsi), %rdx
	am8	%rdi, 14, %rcx, 0, %r11, %rbx, %rbp, %r12, %r13, %r10, %r8, %r9, %rax, %r14, %r15

	mov	%rbx, 15*8(%rdi)
	mov	%rbp, 16*8(%rdi)
	mov	%r12, 17*8(%rdi)
	mov	%r13, 18*8(%rdi)
	mov	%r10, 19*8(%rdi)
	mov	%r8, 20*8(%rdi)
	mov	%r9, 21*8(%rdi)
	mov	%rax, 22*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_16_1)
	mov	0*8(%rdx), %rdx
	m16_str	%rdi, 0, %rsi, 0, %rcx, %r8, %r9, %r10, %r11, %rax
	mov	%rax, 16*8(%rdi)

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_16_2)
	mov	0*8(%rdx), %rcx
	mov	1*8(%rdx), %r8
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	mov	%rcx, %rdx
	xor	%r14d, %r14d
	m6	%rdi, 0, %rsi, 0, %r9, %r10, %r11, %rbx, %rbp, %r12, %rax, %r13, %r14
	mov	%r8, %rdx
	am6	%rdi, 1, %rsi, 0, %r9, %r10, %r11, %rbx, %rbp, %r12, %rax, %r13, %r14
	mov	%r10, 2*8(%rdi)
	mov	%r11, 3*8(%rdi)
	mov	%rbx, 4*8(%rdi)
	mov	%rbp, 5*8(%rdi)

	mov	%rcx, %rdx
	m6_chain	%rdi, 6, %rsi, 6, %r12, %rax, %r9, %r10, %r11, %rbx, %rbp, %rax, %r13, %r12, %r14
	mov	%r8, %rdx
	am6	%rdi, 7, %rsi, 6, %r9, %r10, %r11, %rbx, %rbp, %rax, %r12, %r13, %r14
	mov	%r10, 8*8(%rdi)
	mov	%r11, 9*8(%rdi)
	mov	%rbx, 10*8(%rdi)
	mov	%rbp, 11*8(%rdi)

	mov	%rcx, %rdx
	m4_chain	%rdi, 12, %rsi, 12, %rax, %r12, %r9, %r10, %r11, %rbx, %r13, %rax, %r14
	mov	%r8, %rdx
	am4	%rdi, 13, %rsi, 12, %r9, %r10, %r11, %rbx, %rax, %r13, %r14
	mov	%r10, 14*8(%rdi)
	mov	%r11, 15*8(%rdi)
	mov	%rbx, 16*8(%rdi)
	mov	%rax, 17*8(%rdi)
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_16_3)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx

	xor	%ebx, %ebx

	m3	%rdi, 0, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx

	mov	1*8(%rsi), %rdx
	am3	%rdi, 1, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx
	mov	2*8(%rsi), %rdx
	am3	%rdi, 2, %rcx, 0, %rax, %r9, %r10, %r8, %r11, %rbx
	mov	3*8(%rsi), %rdx
	am3	%rdi, 3, %rcx, 0, %r9, %r10, %r8, %rax, %r11, %rbx
	mov	4*8(%rsi), %rdx
	am3	%rdi, 4, %rcx, 0, %r10, %r8, %rax, %r9, %r11, %rbx
	mov	5*8(%rsi), %rdx
	am3	%rdi, 5, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx
	mov	6*8(%rsi), %rdx
	am3	%rdi, 6, %rcx, 0, %rax, %r9, %r10, %r8, %r11, %rbx
	mov	7*8(%rsi), %rdx
	am3	%rdi, 7, %rcx, 0, %r9, %r10, %r8, %rax, %r11, %rbx
	mov	8*8(%rsi), %rdx
	am3	%rdi, 8, %rcx, 0, %r10, %r8, %rax, %r9, %r11, %rbx
	mov	9*8(%rsi), %rdx
	am3	%rdi, 9, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx
	mov	10*8(%rsi), %rdx
	am3	%rdi, 10, %rcx, 0, %rax, %r9, %r10, %r8, %r11, %rbx
	mov	11*8(%rsi), %rdx
	am3	%rdi, 11, %rcx, 0, %r9, %r10, %r8, %rax, %r11, %rbx
	mov	12*8(%rsi), %rdx
	am3	%rdi, 12, %rcx, 0, %r10, %r8, %rax, %r9, %r11, %rbx
	mov	13*8(%rsi), %rdx
	am3	%rdi, 13, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx
	mov	14*8(%rsi), %rdx
	am3	%rdi, 14, %rcx, 0, %rax, %r9, %r10, %r8, %r11, %rbx
	mov	15*8(%rsi), %rdx
	am3	%rdi, 15, %rcx, 0, %r9, %r10, %r8, %rax, %r11, %rbx

	mov	%r10, 16*8(%rdi)
	mov	%r8, 17*8(%rdi)
	mov	%rax, 18*8(%rdi)

	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_16_4)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp

	xor	%ebp, %ebp

	m4	%rdi, 0, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp

	mov	1*8(%rsi), %rdx
	am4	%rdi, 1, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp
	mov	2*8(%rsi), %rdx
	am4	%rdi, 2, %rcx, 0, %r8, %r9, %rax, %r11, %r10, %rbx, %rbp
	mov	3*8(%rsi), %rdx
	am4	%rdi, 3, %rcx, 0, %r9, %rax, %r11, %r10, %r8, %rbx, %rbp
	mov	4*8(%rsi), %rdx
	am4	%rdi, 4, %rcx, 0, %rax, %r11, %r10, %r8, %r9, %rbx, %rbp
	mov	5*8(%rsi), %rdx
	am4	%rdi, 5, %rcx, 0, %r11, %r10, %r8, %r9, %rax, %rbx, %rbp
	mov	6*8(%rsi), %rdx
	am4	%rdi, 6, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp
	mov	7*8(%rsi), %rdx
	am4	%rdi, 7, %rcx, 0, %r8, %r9, %rax, %r11, %r10, %rbx, %rbp
	mov	8*8(%rsi), %rdx
	am4	%rdi, 8, %rcx, 0, %r9, %rax, %r11, %r10, %r8, %rbx, %rbp
	mov	9*8(%rsi), %rdx
	am4	%rdi, 9, %rcx, 0, %rax, %r11, %r10, %r8, %r9, %rbx, %rbp
	mov	10*8(%rsi), %rdx
	am4	%rdi, 10, %rcx, 0, %r11, %r10, %r8, %r9, %rax, %rbx, %rbp
	mov	11*8(%rsi), %rdx
	am4	%rdi, 11, %rcx, 0, %r10, %r8, %r9, %rax, %r11, %rbx, %rbp
	mov	12*8(%rsi), %rdx
	am4	%rdi, 12, %rcx, 0, %r8, %r9, %rax, %r11, %r10, %rbx, %rbp
	mov	13*8(%rsi), %rdx
	am4	%rdi, 13, %rcx, 0, %r9, %rax, %r11, %r10, %r8, %rbx, %rbp
	mov	14*8(%rsi), %rdx
	am4	%rdi, 14, %rcx, 0, %rax, %r11, %r10, %r8, %r9, %rbx, %rbp
	mov	15*8(%rsi), %rdx
	am4	%rdi, 15, %rcx, 0, %r11, %r10, %r8, %r9, %rax, %rbx, %rbp

	mov	%r10, 16*8(%rdi)
	mov	%r8, 17*8(%rdi)
	mov	%r9, 18*8(%rdi)
	mov	%rax, 19*8(%rdi)

	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_16_5)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12

	xor	%r12d, %r12d

	m5	%rdi, 0, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12

	mov	1*8(%rsi), %rdx
	am5	%rdi, 1, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12
	mov	2*8(%rsi), %rdx
	am5	%rdi, 2, %rcx, 0, %rax, %r9, %r10, %r11, %rbx, %r8, %rbp, %r12
	mov	3*8(%rsi), %rdx
	am5	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %r8, %rax, %rbp, %r12
	mov	4*8(%rsi), %rdx
	am5	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %r8, %rax, %r9, %rbp, %r12
	mov	5*8(%rsi), %rdx
	am5	%rdi, 5, %rcx, 0, %r11, %rbx, %r8, %rax, %r9, %r10, %rbp, %r12
	mov	6*8(%rsi), %rdx
	am5	%rdi, 6, %rcx, 0, %rbx, %r8, %rax, %r9, %r10, %r11, %rbp, %r12
	mov	7*8(%rsi), %rdx
	am5	%rdi, 7, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12
	mov	8*8(%rsi), %rdx
	am5	%rdi, 8, %rcx, 0, %rax, %r9, %r10, %r11, %rbx, %r8, %rbp, %r12
	mov	9*8(%rsi), %rdx
	am5	%rdi, 9, %rcx, 0, %r9, %r10, %r11, %rbx, %r8, %rax, %rbp, %r12
	mov	10*8(%rsi), %rdx
	am5	%rdi, 10, %rcx, 0, %r10, %r11, %rbx, %r8, %rax, %r9, %rbp, %r12
	mov	11*8(%rsi), %rdx
	am5	%rdi, 11, %rcx, 0, %r11, %rbx, %r8, %rax, %r9, %r10, %rbp, %r12
	mov	12*8(%rsi), %rdx
	am5	%rdi, 12, %rcx, 0, %rbx, %r8, %rax, %r9, %r10, %r11, %rbp, %r12
	mov	13*8(%rsi), %rdx
	am5	%rdi, 13, %rcx, 0, %r8, %rax, %r9, %r10, %r11, %rbx, %rbp, %r12
	mov	14*8(%rsi), %rdx
	am5	%rdi, 14, %rcx, 0, %rax, %r9, %r10, %r11, %rbx, %r8, %rbp, %r12
	mov	15*8(%rsi), %rdx
	am5	%rdi, 15, %rcx, 0, %r9, %r10, %r11, %rbx, %r8, %rax, %rbp, %r12

	mov	%r10, 16*8(%rdi)
	mov	%r11, 17*8(%rdi)
	mov	%rbx, 18*8(%rdi)
	mov	%r8, 19*8(%rdi)
	mov	%rax, 20*8(%rdi)

	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_16_6)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	xor	%r13d, %r13d

	m6	%rdi, 0, %rcx, 0, %rbp, %r8, %r9, %r10, %r11, %rbx, %rax, %r12, %r13

	mov	1*8(%rsi), %rdx
	am6	%rdi, 1, %rcx, 0, %rbp, %r8, %r9, %r10, %r11, %rbx, %rax, %r12, %r13
	mov	2*8(%rsi), %rdx
	am6	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12, %r13
	mov	3*8(%rsi), %rdx
	am6	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rbx, %rax, %rbp, %r8, %r12, %r13
	mov	4*8(%rsi), %rdx
	am6	%rdi, 4, %rcx, 0, %r10, %r11, %rbx, %rax, %rbp, %r8, %r9, %r12, %r13
	mov	5*8(%rsi), %rdx
	am6	%rdi, 5, %rcx, 0, %r11, %rbx, %rax, %rbp, %r8, %r9, %r10, %r12, %r13
	mov	6*8(%rsi), %rdx
	am6	%rdi, 6, %rcx, 0, %rbx, %rax, %rbp, %r8, %r9, %r10, %r11, %r12, %r13
	mov	7*8(%rsi), %rdx
	am6	%rdi, 7, %rcx, 0, %rax, %rbp, %r8, %r9, %r10, %r11, %rbx, %r12, %r13
	mov	8*8(%rsi), %rdx
	am6	%rdi, 8, %rcx, 0, %rbp, %r8, %r9, %r10, %r11, %rbx, %rax, %r12, %r13
	mov	9*8(%rsi), %rdx
	am6	%rdi, 9, %rcx, 0, %r8, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12, %r13
	mov	10*8(%rsi), %rdx
	am6	%rdi, 10, %rcx, 0, %r9, %r10, %r11, %rbx, %rax, %rbp, %r8, %r12, %r13
	mov	11*8(%rsi), %rdx
	am6	%rdi, 11, %rcx, 0, %r10, %r11, %rbx, %rax, %rbp, %r8, %r9, %r12, %r13
	mov	12*8(%rsi), %rdx
	am6	%rdi, 12, %rcx, 0, %r11, %rbx, %rax, %rbp, %r8, %r9, %r10, %r12, %r13
	mov	13*8(%rsi), %rdx
	am6	%rdi, 13, %rcx, 0, %rbx, %rax, %rbp, %r8, %r9, %r10, %r11, %r12, %r13
	mov	14*8(%rsi), %rdx
	am6	%rdi, 14, %rcx, 0, %rax, %rbp, %r8, %r9, %r10, %r11, %rbx, %r12, %r13
	mov	15*8(%rsi), %rdx
	am6	%rdi, 15, %rcx, 0, %rbp, %r8, %r9, %r10, %r11, %rbx, %rax, %r12, %r13

	mov	%r8, 16*8(%rdi)
	mov	%r9, 17*8(%rdi)
	mov	%r10, 18*8(%rdi)
	mov	%r11, 19*8(%rdi)
	mov	%rbx, 20*8(%rdi)
	mov	%rax, 21*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_16_7)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14

	xor	%r14d, %r14d

	m7	%rdi, 0, %rcx, 0, %rbx, %r8, %r9, %r10, %r11, %rax, %rbp, %r12, %r13, %r14

	mov	1*8(%rsi), %rdx
	am7	%rdi, 1, %rcx, 0, %rbx, %r8, %r9, %r10, %r11, %rax, %rbp, %r12, %r13, %r14
	mov	2*8(%rsi), %rdx
	am7	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %r11, %rax, %rbp, %r12, %rbx, %r13, %r14
	mov	3*8(%rsi), %rdx
	am7	%rdi, 3, %rcx, 0, %r9, %r10, %r11, %rax, %rbp, %r12, %rbx, %r8, %r13, %r14
	mov	4*8(%rsi), %rdx
	am7	%rdi, 4, %rcx, 0, %r10, %r11, %rax, %rbp, %r12, %rbx, %r8, %r9, %r13, %r14
	mov	5*8(%rsi), %rdx
	am7	%rdi, 5, %rcx, 0, %r11, %rax, %rbp, %r12, %rbx, %r8, %r9, %r10, %r13, %r14
	mov	6*8(%rsi), %rdx
	am7	%rdi, 6, %rcx, 0, %rax, %rbp, %r12, %rbx, %r8, %r9, %r10, %r11, %r13, %r14
	mov	7*8(%rsi), %rdx
	am7	%rdi, 7, %rcx, 0, %rbp, %r12, %rbx, %r8, %r9, %r10, %r11, %rax, %r13, %r14
	mov	8*8(%rsi), %rdx
	am7	%rdi, 8, %rcx, 0, %r12, %rbx, %r8, %r9, %r10, %r11, %rax, %rbp, %r13, %r14
	mov	9*8(%rsi), %rdx
	am7	%rdi, 9, %rcx, 0, %rbx, %r8, %r9, %r10, %r11, %rax, %rbp, %r12, %r13, %r14
	mov	10*8(%rsi), %rdx
	am7	%rdi, 10, %rcx, 0, %r8, %r9, %r10, %r11, %rax, %rbp, %r12, %rbx, %r13, %r14
	mov	11*8(%rsi), %rdx
	am7	%rdi, 11, %rcx, 0, %r9, %r10, %r11, %rax, %rbp, %r12, %rbx, %r8, %r13, %r14
	mov	12*8(%rsi), %rdx
	am7	%rdi, 12, %rcx, 0, %r10, %r11, %rax, %rbp, %r12, %rbx, %r8, %r9, %r13, %r14
	mov	13*8(%rsi), %rdx
	am7	%rdi, 13, %rcx, 0, %r11, %rax, %rbp, %r12, %rbx, %r8, %r9, %r10, %r13, %r14
	mov	14*8(%rsi), %rdx
	am7	%rdi, 14, %rcx, 0, %rax, %rbp, %r12, %rbx, %r8, %r9, %r10, %r11, %r13, %r14
	mov	15*8(%rsi), %rdx
	am7	%rdi, 15, %rcx, 0, %rbp, %r12, %rbx, %r8, %r9, %r10, %r11, %rax, %r13, %r14

	mov	%r12, 16*8(%rdi)
	mov	%rbx, 17*8(%rdi)
	mov	%r8, 18*8(%rdi)
	mov	%r9, 19*8(%rdi)
	mov	%r10, 20*8(%rdi)
	mov	%r11, 21*8(%rdi)
	mov	%rax, 22*8(%rdi)

	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()

	ALIGN(16)
PROLOGUE(flint_mpn_mul_16_8)
	mov	%rdx, %rcx
	mov	0*8(%rsi), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13
	push	%r14
	push	%r15

	xor	%r15d, %r15d

	m8	%rdi, 0, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r14, %r15

	mov	1*8(%rsi), %rdx
	am8	%rdi, 1, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r14, %r15
	mov	2*8(%rsi), %rdx
	am8	%rdi, 2, %rcx, 0, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r11, %r14, %r15
	mov	3*8(%rsi), %rdx
	am8	%rdi, 3, %rcx, 0, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r11, %r8, %r14, %r15
	mov	4*8(%rsi), %rdx
	am8	%rdi, 4, %rcx, 0, %r10, %rax, %rbx, %rbp, %r12, %r13, %r11, %r8, %r9, %r14, %r15
	mov	5*8(%rsi), %rdx
	am8	%rdi, 5, %rcx, 0, %rax, %rbx, %rbp, %r12, %r13, %r11, %r8, %r9, %r10, %r14, %r15
	mov	6*8(%rsi), %rdx
	am8	%rdi, 6, %rcx, 0, %rbx, %rbp, %r12, %r13, %r11, %r8, %r9, %r10, %rax, %r14, %r15
	mov	7*8(%rsi), %rdx
	am8	%rdi, 7, %rcx, 0, %rbp, %r12, %r13, %r11, %r8, %r9, %r10, %rax, %rbx, %r14, %r15
	mov	8*8(%rsi), %rdx
	am8	%rdi, 8, %rcx, 0, %r12, %r13, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r14, %r15
	mov	9*8(%rsi), %rdx
	am8	%rdi, 9, %rcx, 0, %r13, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r14, %r15
	mov	10*8(%rsi), %rdx
	am8	%rdi, 10, %rcx, 0, %r11, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r14, %r15
	mov	11*8(%rsi), %rdx
	am8	%rdi, 11, %rcx, 0, %r8, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r11, %r14, %r15
	mov	12*8(%rsi), %rdx
	am8	%rdi, 12, %rcx, 0, %r9, %r10, %rax, %rbx, %rbp, %r12, %r13, %r11, %r8, %r14, %r15
	mov	13*8(%rsi), %rdx
	am8	%rdi, 13, %rcx, 0, %r10, %rax, %rbx, %rbp, %r12, %r13, %r11, %r8, %r9, %r14, %r15
	mov	14*8(%rsi), %rdx
	am8	%rdi, 14, %rcx, 0, %rax, %rbx, %rbp, %r12, %r13, %r11, %r8, %r9, %r10, %r14, %r15
	mov	15*8(%rsi), %rdx
	am8	%rdi, 15, %rcx, 0, %rbx, %rbp, %r12, %r13, %r11, %r8, %r9, %r10, %rax, %r14, %r15

	mov	%rbp, 16*8(%rdi)
	mov	%r12, 17*8(%rdi)
	mov	%r13, 18*8(%rdi)
	mov	%r11, 19*8(%rdi)
	mov	%r8, 20*8(%rdi)
	mov	%r9, 21*8(%rdi)
	mov	%r10, 22*8(%rdi)
	mov	%rax, 23*8(%rdi)

	pop	%r15
	pop	%r14
	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
EPILOGUE()
