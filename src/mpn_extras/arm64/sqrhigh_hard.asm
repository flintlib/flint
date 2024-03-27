dnl
dnl Copyright (C) 2024 Albin Ahlbäck
dnl
dnl This file is part of FLINT.
dnl
dnl FLINT is free software: you can redistribute it and/or modify it under
dnl the terms of the GNU Lesser General Public License (LGPL) as published
dnl by the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.  See <https://www.gnu.org/licenses/>.
dnl

include(`config.m4')

define(`rp',  `x0')
define(`ap',  `x1')

define(`r0',  `x2')
define(`r1',  `x3')
define(`r2',  `x4')
define(`r3',  `x5')
define(`r4',  `x6')
define(`r5',  `x7')
define(`r6',  `x8')
define(`r7',  `x9')
define(`r8',  `x10')
define(`r9',  `x11')
define(`r10', `x12')
define(`r11', `x13')
define(`r12', `x14')
define(`r13', `x15')
define(`r14', `x16')
define(`r15', `x17')

define(`r16', `ap')	dnl Beware!

dnl NOTE: The following has to be pushed and popped
dnl NOTE: At least on Apple, we need to preserve both x18 and x29 at all times!
define(`r17', `x19')
define(`r18', `x20')
define(`r19', `x21')
define(`r20', `x22')
define(`r21', `x23')
define(`r22', `x24')
define(`r23', `x25')
define(`r24', `x26')
define(`r25', `x27')
define(`r26', `x28')
define(`r27', `x30')

define(`res', `rp') dnl NOTE: Synonymous with rp!

PROLOGUE(flint_mpn_sqrhigh_1)
	ldr	r0, [ap]
	umulh	r1, r0, r0
	mul	r0, r0, r0
	str	r1, [rp]
	mov	res, r0
	ret
EPILOGUE()

PROLOGUE(flint_mpn_sqrhigh_2)
	ldp	r0, r1, [ap]
	umulh	r2, r0, r0
	mul	r3, r0, r1
	umulh	r4, r0, r1
	mul	r5, r1, r1
	umulh	r6, r1, r1
	adds	r2, r2, r3
	adcs	r5, r5, r4
	cinc	r6, r6, cs
	adds	r2, r2, r3
	adcs	r5, r5, r4
	cinc	r6, r6, cs
	stp	r5, r6, [rp]
	mov	res, r2
	ret
EPILOGUE()

C   0 1 2
C 1   h x
C 2   d x
C 3     d
PROLOGUE(flint_mpn_sqrhigh_3)
	ldp	r0, r1, [ap]
	ldr	r2, [ap, #2*8]
	umulh	r3, r0, r1
	mul	r4, r0, r2
	umulh	r0, r0, r2
	mul	r5, r1, r2
	umulh	r6, r1, r2
	mul	r7, r1, r1
	adds	r3, r3, r4
	umulh	r1, r1, r1
	mul	r8, r2, r2
	adcs	r0, r0, r5
	umulh	r2, r2, r2
	cinc	r6, r6, cs
	adds	r3, r3, r3
	adcs	r0, r0, r0
	adcs	r6, r6, r6
	cset	r4, cs
	adds	r7, r7, r3
	adcs	r1, r1, r0
	adcs	r8, r8, r6
	adc	r2, r2, r4
	stp	r1, r8, [rp]
	str	r2, [rp, #2*8]
	mov	res, r7
	ret
EPILOGUE()

C   0 1 2 3
C 0     h x
C 1   e x x
C 2     d x
C 3       d
PROLOGUE(flint_mpn_sqrhigh_4)
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]
	umulh	r4, r0, r2
	mul	r5, r0, r3
	umulh	r6, r0, r3
	umulh	r7, r1, r2
	umulh	r8, r1, r3
	mul	r9, r1, r2
	adds	r4, r4, r5
	mul	r10, r1, r3
	mul	r11, r2, r3
	adcs	r6, r6, r7
	umulh	r12, r2, r3
	umulh	r1, r1, r1
	cinc	r8, r8, cs
	mul	r13, r2, r2
	umulh	r2, r2, r2
	adds	r4, r4, r9
	mul	r14, r3, r3
	umulh	r3, r3, r3
	adcs	r6, r6, r10
	adcs	r8, r8, r11
	cinc	r12, r12, cs
	adds	r4, r4, r4
	adcs	r6, r6, r6
	adcs	r8, r8, r8
	adcs	r12, r12, r12
	cset	r0, cs
	adds	r4, r4, r1
	adcs	r6, r6, r13
	adcs	r8, r8, r2
	adcs	r12, r12, r14
	adc	r0, r0, r3
	stp	r6, r8, [rp]
	stp	r12, r0, [rp, #2*8]
	mov	res, r4

	ret
EPILOGUE()

C   0 1 2 3 4
C 0       h x
C 1     h x x
C 2     d x x
C 3       d x
C 4         d
PROLOGUE(flint_mpn_sqrhigh_5)
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]
	ldr	r4, [ap, #4*8]

	umulh	r5, r0, r3
	mul	r6, r0, r4
	umulh	r7, r0, r4
	umulh	r8, r1, r3
	umulh	r9, r1, r4
	umulh	r10, r1, r2
	adds	r5, r5, r6
	mul	r11, r1, r4
	umulh	r12, r2, r3
	adcs	r7, r7, r8
	umulh	r13, r2, r4
	mul	r14, r1, r3
	cinc	r9, r9, cs
	mul	r15, r2, r3
	mul	r16, r2, r4
	adds	r5, r5, r10
	adcs	r7, r7, r11
	mul	r6, r3, r4
	umulh	r8, r3, r4
	adcs	r9, r9, r12
	cinc	r13, r13, cs
	mul	r0, r2, r2
	umulh	r2, r2, r2
	adds	r5, r5, r14
	adcs	r7, r7, r15
	mul	r1, r3, r3
	umulh	r3, r3, r3
	adcs	r9, r9, r16
	adcs	r13, r13, r6
	mul	r10, r4, r4
	umulh	r4, r4, r4
	cinc	r8, r8, cs
	adds	r5, r5, r5
	adcs	r7, r7, r7
	adcs	r9, r9, r9
	adcs	r13, r13, r13
	adcs	r8, r8, r8
	cset	r11, cs
	adds	r5, r5, r0
	adcs	r7, r7, r2
	adcs	r9, r9, r1
	adcs	r13, r13, r3
	stp	r7, r9, [rp]
	adcs	r8, r8, r10
	stp	r13, r8, [rp, #2*8]
	adc	r11, r11, r4
	str	r11, [rp, #4*8]
	mov	res, r5

	ret
EPILOGUE()

C   0 1 2 3 4 5
C 0         h x
C 1       h x x
C 2     e x x x
C 3       d x x
C 4         d x
C 5           d
PROLOGUE(flint_mpn_sqrhigh_6)
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]
	ldp	r4, r5, [ap, #4*8]
	umulh	r6, r1, r3
	mul	r7, r2, r3
	umulh	r8, r2, r3
	umulh	r9, r1, r4
	umulh	r10, r2, r4
	umulh	r11, r0, r4
	adds	r6, r6, r7
	mul	r12, r2, r4
	mul	r13, r3, r4
	adcs	r8, r8, r9
	umulh	r14, r3, r4
	mul	r15, r1, r4
	cinc	r10, r10, cs
	umulh	r16, r0, r5
	adds	r6, r6, r11
	umulh	r7, r1, r5
	adcs	r8, r8, r12
	umulh	r9, r2, r5
	adcs	r10, r10, r13
	umulh	r11, r3, r5
	cinc	r14, r14, cs
	mul	r0, r0, r5
	adds	r6, r6, r15
	mul	r1, r1, r5
	adcs	r8, r8, r16
	mul	r12, r2, r5
	adcs	r10, r10, r7
	mul	r13, r3, r5
	adcs	r14, r14, r9
	mul	r15, r4, r5
	cinc	r11, r11, cs
	adds	r6, r6, r0
	umulh	r16, r4, r5
	adcs	r8, r8, r1
	umulh	r2, r2, r2
	adcs	r10, r10, r12
	mul	r7, r3, r3
	adcs	r14, r14, r13
	umulh	r3, r3, r3
	adcs	r11, r11, r15
	mul	r9, r4, r4
	cinc	r16, r16, cs
	umulh	r4, r4, r4
	adds	r6, r6, r6
	adcs	r8, r8, r8
	mul	r0, r5, r5
	adcs	r10, r10, r10
	adcs	r14, r14, r14
	umulh	r5, r5, r5
	adcs	r11, r11, r11
	adcs	r16, r16, r16
	cset	r1, cs
	adds	r2, r2, r6
	adcs	r7, r7, r8
	adcs	r3, r3, r10
	stp	r7, r3, [rp]
	adcs	r9, r9, r14
	adcs	r4, r4, r11
	stp	r9, r4, [rp, #2*8]
	adcs	r0, r0, r16
	adc	r5, r5, r1
	stp	r0, r5, [rp, #4*8]
	mov	res, r2

	ret
EPILOGUE()

C   0 1 2 3 4 5 6
C 0           h x
C 1         h x x
C 2       h x x x
C 3       d x x x
C 4         d x x
C 5           d x
C 6             d
PROLOGUE(flint_mpn_sqrhigh_7)
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]
	ldp	r4, r5, [ap, #4*8]
	ldr	r6, [ap, #6*8]

	umulh	r7, r2, r3

	umulh	r8, r1, r4
	umulh	r9, r2, r4
	mul	r10, r3, r4
	umulh	r11, r3, r4
	mul	r12, r2, r4

	adds	r7, r7, r8
	umulh	r13, r1, r5
	umulh	r14, r2, r5
	adcs	r9, r9, r10
	umulh	r15, r3, r5
	cinc	r11, r11, cs
	umulh	r16, r0, r5
	adds	r7, r7, r12
	mul	r8, r2, r5
	adcs	r9, r9, r13
	mul	r10, r3, r5
	adcs	r11, r11, r14
	mul	r12, r4, r5
	cinc	r15, r15, cs
	umulh	r13, r4, r5

	adds	r7, r7, r16
	mul	r14, r1, r5

	adcs	r9, r9, r8
	umulh	r16, r0, r6
	adcs	r11, r11, r10
	umulh	r8, r1, r6
	adcs	r15, r15, r12
	umulh	r10, r2, r6
	cinc	r13, r13, cs
	umulh	r12, r3, r6
	adds	r7, r7, r14
	umulh	r14, r4, r6

	adcs	r9, r9, r16
	mul	r0, r0, r6
	adcs	r11, r11, r8
	mul	r1, r1, r6
	adcs	r15, r15, r10
	mul	r2, r2, r6
	adcs	r13, r13, r12
	mul	r16, r3, r6
	cinc	r14, r14, cs
	mul	r8, r4, r6
	adds	r7, r7, r0
	mul	r10, r5, r6
	adcs	r9, r9, r1
	umulh	r12, r5, r6

	mul	r0, r3, r3
	adcs	r11, r11, r2
	umulh	r3, r3, r3
	adcs	r15, r15, r16
	mul	r1, r4, r4
	adcs	r13, r13, r8
	umulh	r4, r4, r4
	adcs	r14, r14, r10
	mul	r2, r5, r5
	cinc	r12, r12, cs
	umulh	r5, r5, r5
	adds	r7, r7, r7
	adcs	r9, r9, r9
	mul	r16, r6, r6
	adcs	r11, r11, r11
	adcs	r15, r15, r15
	umulh	r6, r6, r6
	adcs	r13, r13, r13
	adcs	r14, r14, r14
	adcs	r12, r12, r12
	cset	r8, cs

	adds	r0, r0, r7
	adcs	r3, r3, r9
	adcs	r1, r1, r11
	stp	r3, r1, [rp]
	adcs	r4, r4, r15
	adcs	r2, r2, r13
	stp	r4, r2, [rp, #2*8]
	adcs	r5, r5, r14
	adcs	r16, r16, r12
	stp	r5, r16, [rp, #4*8]
	adc	r6, r6, r8
	str	r6, [rp, #6*8]
	mov	res, r0

	ret
EPILOGUE()

C   0 1 2 3 4 5 6 7
C 0             h x
C 1           h x x
C 2         h x x x
C 3       e x x x x
C 4         d x x x
C 5           d x x
C 6             d x
C 7               d
PROLOGUE(flint_mpn_sqrhigh_8)
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]
	ldp	r4, r5, [ap, #4*8]
	ldp	r6, r7, [ap, #6*8]

	umulh	r8, r2, r4
	mul	r9, r3, r4
	umulh	r10, r3, r4

	umulh	r11, r2, r5
	umulh	r12, r3, r5
	umulh	r13, r1, r5
	adds	r8, r8, r9
	mul	r14, r3, r5
	mul	r15, r4, r5
	adcs	r10, r10, r11
	umulh	r16, r4, r5
	cinc	r12, r12, cs
	mul	r9, r2, r5

	adds	r8, r8, r13
	umulh	r11, r1, r6
	adcs	r10, r10, r14
	umulh	r13, r2, r6
	adcs	r12, r12, r15
	umulh	r14, r3, r6
	cinc	r16, r16, cs
	umulh	r15, r4, r6
	adds	r8, r8, r9
	mul	r9, r1, r6
	adcs	r10, r10, r11
	mul	r11, r2, r6
	adcs	r12, r12, r13
	mul	r13, r3, r6
	adcs	r16, r16, r14
	mul	r14, r4, r6
	cinc	r15, r15, cs
	adds	r8, r8, r9
	mul	r9, r5, r6
	adcs	r10, r10, r11
	umulh	r11, r5, r6
	adcs	r12, r12, r13
	umulh	r13, r0, r6

	adcs	r16, r16, r14
	umulh	r14, r0, r7
	adcs	r15, r15, r9
	umulh	r9, r1, r7
	cinc	r11, r11, cs
	mul	r0, r0, r7	C Reduce latency
	adds	r8, r8, r13
	umulh	r13, r2, r7
	adcs	r10, r10, r14
	umulh	r14, r3, r7
	adcs	r12, r12, r9
	umulh	r9, r4, r7
	adcs	r16, r16, r13
	umulh	r13, r5, r7
	mul	r1, r1, r7
	adcs	r15, r15, r14
	mul	r2, r2, r7
	adcs	r11, r11, r9
	mul	r14, r3, r7
	cinc	r13, r13, cs
	mul	r9, r4, r7
	adds	r8, r8, r0
	mul	r0, r5, r7
	adcs	r10, r10, r1
	mul	r1, r6, r7
	adcs	r12, r12, r2
	umulh	r2, r6, r7

	adcs	r16, r16, r14
	umulh	r3, r3, r3
	adcs	r15, r15, r9
	mul	r14, r4, r4
	adcs	r11, r11, r0
	umulh	r4, r4, r4
	adcs	r13, r13, r1
	mul	r9, r5, r5
	cinc	r2, r2, cs
	umulh	r5, r5, r5

	adds	r8, r8, r8
	mul	r0, r6, r6
	adcs	r10, r10, r10
	umulh	r6, r6, r6
	adcs	r12, r12, r12
	adcs	r16, r16, r16
	adcs	r15, r15, r15
	adcs	r11, r11, r11
	adcs	r13, r13, r13
	adcs	r2, r2, r2
	cset	r1, cs

	adds	r3, r3, r8
	mul	r8, r7, r7
	adcs	r14, r14, r10
	umulh	r7, r7, r7
	adcs	r4, r4, r12
	stp	r14, r4, [rp]
	adcs	r9, r9, r16
	adcs	r5, r5, r15
	stp	r9, r5, [rp, #2*8]
	adcs	r0, r0, r11
	adcs	r6, r6, r13
	stp	r0, r6, [rp, #4*8]
	adcs	r8, r8, r2
	adc	r7, r7, r1
	stp	r8, r7, [rp, #6*8]
	mov	res, r3

	ret
EPILOGUE()
