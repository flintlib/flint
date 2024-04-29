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

define(`zr',  `xzr')

define(`res', `rp')	dnl NOTE: Synonymous with rp!

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

PROLOGUE(flint_mpn_sqr_1)
	ldr	r0, [ap]

	mul	r1, r0, r0
	umulh	r2, r0, r0

	stp	r1, r2, [rp]
	mov	res, r2

	ret
EPILOGUE()

PROLOGUE(flint_mpn_sqr_2)
	ldp	r0, r1, [ap]

	mul	r2, r0, r1
	umulh	r3, r0, r1
	C 2, 3

	mul	r4, r0, r0
	umulh	r0, r0, r0
	mul	r5, r1, r1
	umulh	r1, r1, r1
	C 4, 0, 5, 1

	adds	r2, r2, r2
	adcs	r3, r3, r3
	cset	r6, cs
	C 2, 3, 6

	adds	r0, r0, r2
	adcs	r5, r5, r3
	adc	r1, r1, r6
	C 4, 0, 5, 1

	stp	r4, r0, [rp]
	stp	r5, r1, [rp, #2*8]
	mov	res, r1

	ret
EPILOGUE()

PROLOGUE(flint_mpn_sqr_3)
	ldp	r0, r1, [ap]
	ldr	r2, [ap, #2*8]

	mul	r3, r0, r1
	umulh	r4, r0, r1
	mul	r5, r0, r2
	umulh	r6, r0, r2
	mul	r7, r1, r2
	umulh	r8, r1, r2
	C 3, (4, 5), (6, 7), 8

	mul	r9, r0, r0
	umulh	r0, r0, r0
	mul	r10, r1, r1
	umulh	r1, r1, r1

	adds	r4, r4, r5
	adcs	r6, r6, r7
	cinc	r8, r8, cs
	C 3, 4, 6, 8

	mul	r11, r2, r2
	umulh	r2, r2, r2
	C 9, 0, 10, 1, 11, 2

	adds	r3, r3, r3
	adcs	r4, r4, r4
	adcs	r6, r6, r6
	adcs	r8, r8, r8
	cset	r12, cs
	C 3, 4, 6, 8, 12

	C    3,  4, 6,  8, 12
	C 9, 0, 10, 1, 11,  2
	adds	r0, r0, r3

	adcs	r10, r10, r4
	adcs	r1, r1, r6
	adcs	r11, r11, r8

	stp	r9, r0, [rp]

	adc	r2, r2, r12

	stp	r10, r1, [rp, #2*8]
	stp	r11, r2, [rp, #4*8]
	mov	res, r2

	ret
EPILOGUE()

PROLOGUE(flint_mpn_sqr_4)
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]

	mul	r4, r0, r1
	umulh	r5, r0, r1
	mul	r6, r0, r2
	umulh	r7, r0, r2
	mul	r8, r0, r3
	umulh	r9, r0, r3
	umulh	r10, r1, r2
	umulh	r11, r1, r3
	C 4, (5, 6), (7, 8), (9, 10), 11

	mul	r12, r1, r2
	mul	r13, r1, r3

	adds	r5, r5, r6

	mul	r14, r2, r3
	umulh	r15, r2, r3

	adcs	r7, r7, r8

	mul	r6, r0, r0
	umulh	r0, r0, r0

	adcs	r9, r9, r10

	mul	r8, r1, r1
	umulh	r1, r1, r1

	cinc	r11, r11, cs

	adds	r7, r7, r12
	adcs	r9, r9, r13

	mul	r10, r2, r2
	umulh	r2, r2, r2

	adcs	r11, r11, r14
	cinc	r15, r15, cs
	C 4, 5, 7, 9, 11, 15

	adds	r4, r4, r4

	mul	r12, r3, r3
	umulh	r3, r3, r3
	C 6, 0, 8, 1, 10, 2, 12, 3

	adcs	r5, r5, r5
	adcs	r7, r7, r7
	adcs	r9, r9, r9
	adcs	r11, r11, r11
	adcs	r15, r15, r15
	cset	r13, cs
	C 4, 5, 7, 9, 11, 15, 13

	C    4, 5, 7,  9, 11, 15, 13
	C 6, 0, 8, 1, 10,  2, 12,  3
	adds	r0, r0, r4
	adcs	r8, r8, r5
	adcs	r1, r1, r7
	stp	r6, r0, [rp]
	adcs	r10, r10, r9
	adcs	r2, r2, r11
	adcs	r12, r12, r15
	stp	r8, r1, [rp, #2*8]
	adc	r3, r3, r13
	stp	r10, r2, [rp, #4*8]
	stp	r12, r3, [rp, #6*8]
	mov	res, r3

	ret
EPILOGUE()

PROLOGUE(flint_mpn_sqr_5)
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]
	ldr	r4, [ap, #4*8]

	stp	r17, r18, [sp, #-2*8]!

	mul	r5, r0, r1
	umulh	r6, r0, r1
	mul	r7, r0, r2
	umulh	r8, r0, r2
	mul	r9, r0, r3
	umulh	r10, r0, r3
	mul	r11, r0, r4
	umulh	r12, r0, r4
	umulh	r13, r1, r3
	umulh	r14, r1, r4

	adds	r6, r6, r7

	mul	r15, r1, r2
	umulh	r16, r1, r2

	adcs	r8, r8, r9

	mul	r17, r1, r4
	umulh	r18, r2, r3

	adcs	r10, r10, r11

	mul	r7, r1, r3
	mul	r9, r2, r3

	adcs	r12, r12, r13
	cinc	r14, r14, cs

	mul	r11, r2, r4
	umulh	r13, r2, r4

	adds	r8, r8, r15
	adcs	r10, r10, r16

	mul	r15, r3, r4
	umulh	r16, r3, r4

	adcs	r12, r12, r17
	adcs	r14, r14, r18
	cinc	r13, r13, cs
	C 5, 6, 8, (10, 7), (12, 9), (14, 11), (13, 15), 16
	C 17, 18

	mul	r17, r0, r0
	umulh	r0, r0, r0

	adds	r10, r10, r7
	adcs	r12, r12, r9

	mul	r18, r1, r1
	umulh	r1, r1, r1

	adcs	r14, r14, r11

	mul	r7, r2, r2
	umulh	r2, r2, r2

	adcs	r13, r13, r15
	cinc	r16, r16, cs
	C 5, 6, 8, 10, 12, 14, 13, 16
	C 17, 0, 18, 1, 7, 2
	C Free: 9, 11, 15

	adds	r5, r5, r5

	mul	r9, r3, r3
	umulh	r3, r3, r3

	adcs	r6, r6, r6
	adcs	r8, r8, r8
	adcs	r10, r10, r10

	mul	r11, r4, r4
	umulh	r4, r4, r4

	adcs	r12, r12, r12
	adcs	r14, r14, r14
	adcs	r13, r13, r13
	adcs	r16, r16, r16
	cset	r15, cs
	C     5,  6, 8, 10, 12, 14, 13, 16, 15
	C 17, 0, 18, 1,  7,  2,  9,  3, 11,  4

	adds	r0, r0, r5
	adcs	r6, r6, r18
	adcs	r1, r1, r8
	adcs	r7, r7, r10
	stp	r17, r0, [rp]
	adcs	r2, r2, r12
	adcs	r9, r9, r14
	adcs	r3, r3, r13
	stp	r6, r1, [rp, #2*8]
	adcs	r11, r11, r16
	adc	r4, r4, r15
	stp	r7, r2, [rp, #4*8]
	ldp	r17, r18, [sp], #2*8
	stp	r9, r3, [rp, #6*8]
	stp	r11, r4, [rp, #8*8]

	mov	res, r4

	ret
EPILOGUE()

PROLOGUE(flint_mpn_sqr_6)
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]
	ldp	r4, r5, [ap, #4*8]

	stp	r17, r18, [sp, #-4*8]!
	stp	r19, r20, [sp, #2*8]

	C mul	rX, r0, r1	First computed when doubling
	umulh	r6, r0, r1
	mul	r7, r0, r2
	umulh	r8, r0, r2
	mul	r9, r0, r3
	umulh	r10, r0, r3
	mul	r11, r0, r4
	umulh	r12, r0, r4
	mul	r13, r0, r5

	adds	r6, r6, r7

	umulh	r14, r0, r5
	C 6, (8 + cs, 9), (10, 11), (12, 13), 14

	mul	r15, r1, r2

	adcs	r8, r8, r9

	mul	r16, r1, r3
	mul	r17, r1, r4

	adcs	r10, r10, r11

	mul	r18, r1, r5
	umulh	r19, r1, r5

	adcs	r12, r12, r13

	umulh	r20, r1, r2
	umulh	r7, r1, r3

	cinc	r14, r14, cs

	umulh	r9, r1, r4
	C 6, (8, 15), (10, 16, 20), (12, 17, 7), (14, 18, 9), 19

	umulh	r11, r2, r4

	adds	r8, r8, r15
	adcs	r10, r10, r16
	adcs	r12, r12, r17

	umulh	r13, r2, r5
	mul	r15, r2, r3

	adcs	r14, r14, r18
	cinc	r19, r19, cs

	umulh	r16, r2, r3
	mul	r17, r2, r5

	adds	r10, r10, r20
	adcs	r12, r12, r7
	C 6, 8, 10, (12, 15), (14 + cs, 9, 16), (19, 11, 17), 13

	umulh	r18, r3, r4
	umulh	r20, r3, r5

	adcs	r14, r14, r9
	adcs	r19, r19, r11

	mul	r7, r2, r4
	mul	r9, r3, r4

	cinc	r13, r13, cs
	adds	r12, r12, r15

	mul	r11, r3, r5
	mul	r15, r4, r5

	adcs	r14, r14, r16
	adcs	r19, r19, r17

	umulh	r16, r4, r5
	mul	r17, r0, r1

	adcs	r13, r13, r18
	cinc	r20, r20, cs
	C 17, 6, 8, 10, 12, (14, 7), (19, 9), (13, 11), (20, 15), 16
	C Free: 18

	mul	r18, r0, r0
	umulh	r0, r0, r0

	adds	r14, r14, r7
	adcs	r19, r19, r9

	mul	r7, r1, r1
	umulh	r1, r1, r1

	adcs	r13, r13, r11
	adcs	r20, r20, r15

	mul	r9, r2, r2
	umulh	r2, r2, r2

	cinc	r16, r16, cs
	C 17, 6, 8, 10, 12, 14, 19, 13, 20, 16
	C 18, 0, 7, 1, 9, 2
	C Free: 11, 15

	adds	r17, r17, r17
	adcs	r6, r6, r6

	mul	r11, r3, r3
	umulh	r3, r3, r3

	adcs	r8, r8, r8
	adcs	r10, r10, r10
	adcs	r12, r12, r12
	adcs	r14, r14, r14
	lsr	r15, r16, #63
	adcs	r19, r19, r19
	adcs	r13, r13, r13
	adcs	r20, r20, r20
	adc	r16, r16, r16
	C     17, 6, 8, 10, 12, 14, 19, 13, 20, 16, 15
	C 18,  0, 7, 1,  9,  2, 11,  3
	C Free:

	adds	r0, r0, r17

	mul	r17, r4, r4
	umulh	r4, r4, r4

	adcs	r7, r7, r6
	adcs	r1, r1, r8
	adcs	r9, r9, r10

	stp	r18, r0, [rp]

	mul	r6, r5, r5
	umulh	r5, r5, r5

	adcs	r2, r2, r12
	adcs	r11, r11, r14
	adcs	r3, r3, r19

	stp	r7, r1, [rp, #2*8]

	adcs	r13, r13, r17
	adcs	r4, r4, r20

	stp	r9, r2, [rp, #4*8]

	ldp	r19, r20, [sp, #2*8]

	stp	r11, r3, [rp, #6*8]

	adcs	r6, r6, r16
	adcs	r5, r5, r15

	stp	r13, r4, [rp, #8*8]

	ldp	r17, r18, [sp], #4*8

	stp	r6, r5, [rp, #10*8]

	mov	res, r5

	ret
EPILOGUE()

PROLOGUE(flint_mpn_sqr_7)
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]
	ldp	r4, r5, [ap, #4*8]
	stp	r17, r18, [sp, #-6*8]!

	ldr	r6, [ap, #6*8]
	stp	r19, r20, [sp, #2*8]
	stp	r21, r22, [sp, #4*8]

	C mul	rX, r0, r1	First computed when doubling
	umulh	r7, r0, r1
	mul	r8, r0, r2
	umulh	r9, r0, r2
	mul	r10, r0, r3
	umulh	r11, r0, r3
	mul	r12, r0, r4
	umulh	r13, r0, r4
	mul	r14, r0, r5

	adds	r7, r7, r8

	umulh	r15, r0, r5
	mul	r16, r0, r6

	adcs	r9, r9, r10

	umulh	r17, r0, r6
	C 7, 9, (11 + cs, 12), (13, 14), (15, 16), 17

	mul	r18, r1, r2

	adcs	r11, r11, r12

	mul	r19, r1, r3
	mul	r20, r1, r4

	adcs	r13, r13, r14

	mul	r21, r1, r5
	mul	r22, r1, r6

	adcs	r15, r15, r16

	umulh	r8, r1, r6
	umulh	r10, r1, r2

	cinc	r17, r17, cs

	umulh	r12, r1, r3
	umulh	r14, r1, r4

	adds	r9, r9, r18

	umulh	r16, r1, r5
	mul	r18, r2, r3
	C 7, 9, (11 + cs, 19, 10), (13, 20, 12, 18), (15, 21, 14), (17, 22, 16), 8
	C Free:

	adcs	r11, r11, r19
	adcs	r13, r13, r20

	mul	r19, r2, r4
	mul	r20, r2, r5

	adcs	r15, r15, r21
	adcs	r17, r17, r22

	mul	r21, r2, r6
	umulh	r22, r2, r6

	cinc	r8, r8, cs

	adds	r11, r11, r10
	adcs	r13, r13, r12

	umulh	r10, r2, r3
	umulh	r12, r2, r4

	adcs	r15, r15, r14
	adcs	r17, r17, r16

	umulh	r14, r2, r5
	umulh	r16, r3, r5

	cinc	r8, r8, cs

	adds	r13, r13, r18
	adcs	r15, r15, r19
	C 7, 9, 11, 13, (15, 10), (17 + cs, 20, 12), (8, 21, 14), (22, 16)

	umulh	r18, r3, r6
	mul	r19, r3, r4

	adcs	r17, r17, r20
	adcs	r8, r8, r21

	umulh	r20, r3, r4
	mul	r21, r3, r6

	cinc	r22, r22, cs

	adds	r15, r15, r10
	adcs	r17, r17, r12

	umulh	r10, r4, r5
	umulh	r12, r4, r6

	adcs	r8, r8, r14
	adcs	r22, r22, r16

	mul	r14, r3, r5
	mul	r16, r4, r5

	cinc	r18, r18, cs

	adds	r17, r17, r19
	adcs	r8, r8, r20

	mul	r19, r4, r6
	mul	r20, r5, r6

	adcs	r22, r22, r21
	adcs	r18, r18, r10

	umulh	r21, r5, r6
	mul	r10, r0, r1

	cinc	r12, r12, cs
	C 10, 7, 9, 11, 13, 15, 17, (8, 14), (22, 16), (18, 19), (12, 20), 21

	adds	r8, r8, r14
	adcs	r22, r22, r16

	mul	r14, r0, r0
	umulh	r0, r0, r0

	adcs	r18, r18, r19
	adcs	r12, r12, r20
	cinc	r21, r21, cs

	mul	r16, r1, r1
	umulh	r1, r1, r1
	C 10, 7, 9, 11, 13, 15, 17, 8, 22, 18, 12, 21
	C 14, 0, 16, 1
	C Free: 19, 20

	adds	r10, r10, r10
	adcs	r7, r7, r7
	adcs	r9, r9, r9
	adcs	r11, r11, r11
	adcs	r13, r13, r13

	mul	r20, r2, r2
	umulh	r2, r2, r2

	adcs	r15, r15, r15
	adcs	r17, r17, r17
	adcs	r8, r8, r8
	lsr	r19, r21, #63
	adcs	r22, r22, r22
	adcs	r18, r18, r18
	adcs	r12, r12, r12
	adc	r21, r21, r21
	C    10,  7, 9, 11, 13, 15, 17, 8, 22, 18, 12, 21, 19
	C 14, 0, 16, 1, 20,  2
	C Free:

	adds	r0, r0, r10

	mul	r10, r3, r3
	umulh	r3, r3, r3

	adcs	r16, r16, r7
	adcs	r1, r1, r9

	stp	r14, r0, [rp]

	mul	r7, r4, r4
	umulh	r4, r4, r4

	adcs	r20, r20, r11
	adcs	r2, r2, r13

	stp	r16, r1, [rp, #2*8]

	mul	r9, r5, r5
	umulh	r5, r5, r5

	adcs	r10, r10, r15
	adcs	r3, r3, r17

	stp	r20, r2, [rp, #4*8]

	mul	r14, r6, r6
	umulh	r6, r6, r6

	adcs	r7, r7, r8
	adcs	r4, r4, r22

	stp	r10, r3, [rp, #6*8]

	adcs	r9, r9, r18
	adcs	r5, r5, r12
	adcs	r14, r14, r21
	adcs	r6, r6, r19

	stp	r7, r4, [rp, #8*8]
	ldp	r19, r20, [sp, #2*8]
	stp	r9, r5, [rp, #10*8]
	ldp	r21, r22, [sp, #4*8]
	stp	r14, r6, [rp, #12*8]
	ldp	r17, r18, [sp], #6*8

	mov	res, r6

	ret
EPILOGUE()

dnl From a{n - 4} forward, use the following sequence:
dnl
dnl	umulh	rX, r2, r4
dnl	umulh	rX, r2, r5
dnl	mul	rX, r2, r3
dnl	umulh	rX, r2, r3
dnl	mul	rX, r2, r5
dnl	umulh	rX, r3, r4
dnl	umulh	rX, r3, r5
dnl	mul	rX, r2, r4
dnl	mul	rX, r3, r4
dnl	mul	rX, r3, r5
dnl	mul	rX, r4, r5
dnl	umulh	rX, r4, r5
dnl	mul	rX, r0, r1

PROLOGUE(flint_mpn_sqr_8)
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]
	ldp	r4, r5, [ap, #4*8]
	stp	r17, r18, [sp, #-10*8]!
	stp	r19, r20, [sp, #2*8]
	ldp	r6, r7, [ap, #6*8]
	stp	r21, r22, [sp, #4*8]
	stp	r23, r24, [sp, #6*8]
	stp	r25, r26, [sp, #8*8]

	C mul	rX, r0, r1	Calculate this one later
	umulh	r8, r0, r1
	mul	r9, r0, r2
	umulh	r10, r0, r2
	mul	r11, r0, r3
	umulh	r12, r0, r3
	mul	r13, r0, r4
	umulh	r14, r0, r4
	mul	r15, r0, r5

	adds	r8, r8, r9

	umulh	r16, r0, r5
	mul	r17, r0, r6

	adcs	r10, r10, r11

	umulh	r18, r0, r6
	mul	r19, r0, r7

	adcs	r12, r12, r13

	umulh	r20, r0, r7
	mul	r21, r1, r2

	adcs	r14, r14, r15

	mul	r22, r1, r3
	mul	r23, r1, r4

	adcs	r16, r16, r17

	mul	r24, r1, r5
	mul	r25, r1, r6

	adcs	r18, r18, r19

	mul	r26, r1, r7
	umulh	r9, r1, r7

	cinc	r20, r20, cs
	C 8, (10, 21), (12, 22), (14, 23), (16, 24), (18, 25), (20, 26), 9
	C Free: 11, 13, 15, 17, 19

	umulh	r11, r1, r2
	umulh	r13, r1, r3

	adds	r10, r10, r21
	adcs	r12, r12, r22

	umulh	r15, r1, r4
	umulh	r17, r1, r5

	adcs	r14, r14, r23
	adcs	r16, r16, r24

	umulh	r19, r1, r6
	mul	r21, r2, r3

	adcs	r18, r18, r25
	adcs	r20, r20, r26

	mul	r22, r2, r4
	mul	r23, r2, r5

	cinc	r9, r9, cs
	C 8, 10, (12, 11), (14, 13, 21), (16, 15, 22), (18, 17, 23), (20, 19), 9
	C Free: 24, 25, 26

	adds	r12, r12, r11

	mul	r24, r2, r6
	mul	r25, r2, r7

	adcs	r14, r14, r13
	adcs	r16, r16, r15

	umulh	r26, r2, r7
	umulh	r11, r2, r3

	adcs	r18, r18, r17
	adcs	r20, r20, r19

	umulh	r13, r2, r4
	umulh	r15, r2, r5

	cinc	r9, r9, cs
	C 8, 10, 12, (14, 21), (16, 22, 11), (18, 23, 13), (20, 24, 15), (9, 25), 26
	C Free: 17, 19

	adds	r14, r14, r21

	umulh	r17, r2, r6
	mul	r19, r3, r4

	adcs	r16, r16, r22
	adcs	r18, r18, r23

	mul	r21, r3, r5
	mul	r22, r3, r6

	adcs	r20, r20, r24
	adcs	r9, r9, r25

	mul	r23, r3, r7
	umulh	r24, r3, r7

	cinc	r26, r26, cs
	C 8, 10, 12, 14, (16, 11), (18, 13, 19), (20, 15, 21), (9, 17, 22), (26, 23), 24
	C Free: 25

	adds	r16, r16, r11

	umulh	r25, r3, r4
	umulh	r11, r3, r5

	adcs	r18, r18, r13
	adcs	r20, r20, r15

	umulh	r13, r3, r6
	umulh	r15, r4, r6

	adcs	r9, r9, r17
	cinc	r26, r26, cs
	C 8, 10, 12, 14, 16, (18, 19), (20, 21, 25), (9, 22, 11), (26, 23, 13), (24, 15)
	C Free: 17

	adds	r18, r18, r19

	umulh	r17, r4, r7
	mul	r19, r4, r5

	adcs	r20, r20, r21
	adcs	r9, r9, r22

	umulh	r21, r4, r5
	mul	r22, r4, r7

	adcs	r26, r26, r23
	cinc	r24, r24, cs
	C 8, 10, 12, 14, 16, 18, (20, 25), (9, 11, 19), (26, 13, 21), (24, 15, 22), 17
	C Free: 23

	adds	r20, r20, r25

	umulh	r23, r5, r6
	umulh	r25, r5, r7

	adcs	r9, r9, r11
	adcs	r26, r26, r13

	mul	r11, r4, r6
	mul	r13, r5, r6

	adcs	r24, r24, r15
	cinc	r17, r17, cs
	C 8, 10, 12, 14, 16, 18, 20, (9, 19), (26, 21, 11), (24, 22, 13), (17, 23), 25
	C Free: 15

	adds	r9, r9, r19

	mul	r15, r5, r7
	mul	r19, r6, r7

	adcs	r26, r26, r21
	adcs	r24, r24, r22

	umulh	r21, r6, r7
	mul	r22, r0, r1

	adcs	r17, r17, r23
	cinc	r25, r25, cs

	mul	r23, r0, r0
	umulh	r0, r0, r0

	adds	r26, r26, r11
	adcs	r24, r24, r13

	mul	r11, r1, r1
	umulh	r1, r1, r1

	adcs	r17, r17, r15
	adcs	r25, r25, r19

	mul	r13, r2, r2
	umulh	r2, r2, r2

	cinc	r21, r21, cs
	C 22, 8, 10, 12, 14, 16, 18, 20, 9, 26, 24, 17, 25, 21
	C 23, 0, 11, 1, 13, 2
	C Free: 15, 19

	adds	r22, r22, r22
	adcs	r8, r8, r8
	adcs	r10, r10, r10
	adcs	r12, r12, r12
	adcs	r14, r14, r14
	adcs	r16, r16, r16

	mul	r19, r3, r3
	umulh	r3, r3, r3

	adcs	r18, r18, r18
	adcs	r20, r20, r20
	adcs	r9, r9, r9
	lsr	r15, r21, #63
	adcs	r26, r26, r26
	adcs	r24, r24, r24
	adcs	r17, r17, r17
	adcs	r25, r25, r25
	adc	r21, r21, r21
	C     22,  8, 10, 12, 14, 16, 18, 20, 9, 26, 24, 17, 25, 21, 15
	C 23,  0, 11,  1, 13,  2, 19,  3
	C Free:

	adds	r0, r0, r22

	mul	r22, r4, r4
	umulh	r4, r4, r4

	adcs	r11, r11, r8
	adcs	r1, r1, r10

	stp	r23, r0, [rp]

	mul	r8, r5, r5
	umulh	r5, r5, r5

	adcs	r13, r13, r12
	adcs	r2, r2, r14

	stp	r11, r1, [rp, #2*8]

	mul	r10, r6, r6
	umulh	r6, r6, r6

	adcs	r16, r16, r19
	adcs	r3, r3, r18

	stp	r13, r2, [rp, #4*8]

	mul	r12, r7, r7
	umulh	r7, r7, r7

	adcs	r22, r22, r20
	adcs	r4, r4, r9

	stp	r16, r3, [rp, #6*8]
	ldp	r19, r20, [sp, #2*8]

	adcs	r8, r8, r26
	adcs	r5, r5, r24

	stp	r22, r4, [rp, #8*8]
	ldp	r23, r24, [sp, #6*8]

	adcs	r10, r10, r17
	adcs	r6, r6, r25

	stp	r8, r5, [rp, #10*8]
	ldp	r25, r26, [sp, #8*8]

	adcs	r12, r12, r21
	adc	r7, r7, r15

	stp	r10, r6, [rp, #12*8]
	ldp	r21, r22, [sp, #4*8]
	stp	r12, r7, [rp, #14*8]
	ldp	r17, r18, [sp], #10*8

	mov	res, r7

	ret
EPILOGUE()

PROLOGUE(flint_mpn_sqr_9)
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]
	ldp	r4, r5, [ap, #4*8]
	stp	r17, r18, [sp, #-12*8]!
	stp	r19, r20, [sp, #2*8]
	ldp	r6, r7, [ap, #6*8]
	ldr	r8, [ap, #8*8]
	stp	r21, r22, [sp, #4*8]
	stp	r23, r24, [sp, #6*8]
	stp	r25, r26, [sp, #8*8]
	str	r27, [sp, #10*8]

	C mul	rX, r0, r1	Calculated later
	umulh	r9, r0, r1
	mul	r10, r0, r2
	umulh	r11, r0, r2
	mul	r12, r0, r3
	umulh	r13, r0, r3
	mul	r14, r0, r4
	umulh	r15, r0, r4
	mul	r16, r0, r5

	adds	r9, r9, r10

	umulh	r17, r0, r5
	mul	r18, r0, r6

	adcs	r11, r11, r12

	umulh	r19, r0, r6
	mul	r20, r0, r7

	adcs	r13, r13, r14

	umulh	r21, r0, r7
	mul	r22, r0, r8

	adcs	r15, r15, r16

	umulh	r23, r0, r8
	mul	r24, r1, r2

	adcs	r17, r17, r18

	mul	r25, r1, r3
	mul	r26, r1, r4

	adcs	r19, r19, r20

	mul	r27, r1, r5
	mul	r10, r1, r6

	adcs	r21, r21, r22

	mul	r12, r1, r7
	mul	r14, r1, r8

	cinc	r23, r23, cs
	C 9, (11, 24), (13, 25), (15, 26), (17, 27), (19, 10), (21, 12), (23, 14)
	C Free: 16, 18, 20, 22

	umulh	r16, r1, r8
	umulh	r18, r1, r2

	adds	r11, r11, r24
	adcs	r13, r13, r25

	umulh	r20, r1, r3
	umulh	r22, r1, r4

	adcs	r15, r15, r26
	adcs	r17, r17, r27

	umulh	r24, r1, r5
	umulh	r25, r1, r6

	adcs	r19, r19, r10
	adcs	r21, r21, r12

	umulh	r26, r1, r7
	mul	r27, r2, r3

	adcs	r23, r23, r14
	cinc	r16, r16, cs
	C 9, 11, (13, 18), (15, 20, 27), (17, 22), (19, 24), (21, 25), (23, 26), 16
	C Free: 10, 12, 14

	mul	r10, r2, r4
	mul	r12, r2, r5

	adds	r13, r13, r18
	adcs	r15, r15, r20

	mul	r14, r2, r6
	mul	r18, r2, r7

	adcs	r17, r17, r22
	adcs	r19, r19, r24

	mul	r20, r2, r8
	umulh	r22, r2, r8

	adcs	r21, r21, r25
	adcs	r23, r23, r26

	umulh	r24, r2, r3
	umulh	r25, r2, r4

	cinc	r16, r16, cs
	C 9, 11, 13, (15, 27), (17, 10, 24), (19, 12, 25), (21, 14), (23, 18), (16, 20), 22
	C Free: 26

	adds	r15, r15, r27

	umulh	r26, r2, r5
	umulh	r27, r2, r6

	adcs	r17, r17, r10
	adcs	r19, r19, r12

	umulh	r10, r2, r7
	mul	r12, r3, r4

	adcs	r21, r21, r14
	adcs	r23, r23, r18

	mul	r14, r3, r5
	mul	r18, r3, r6

	adcs	r16, r16, r20
	cinc	r22, r22, cs
	C 9, 11, 13, 15, (17, 24), (19, 25, 12), (21, 26, 14), (23, 27, 18), (16, 10), 22
	C Free: 20

	adds	r17, r17, r24

	mul	r20, r3, r7
	mul	r24, r3, r8

	adcs	r19, r19, r25
	adcs	r21, r21, r26

	umulh	r25, r3, r8
	umulh	r26, r3, r4

	adcs	r23, r23, r27
	adcs	r16, r16, r10

	umulh	r27, r3, r5
	umulh	r10, r3, r6

	cinc	r22, r22, cs
	C 9, 11, 13, 15, 17, (19, 12), (21, 14, 26), (23, 18, 27), (16, 20, 10), (22, 24), 25

	adds	r19, r19, r12
	adcs	r21, r21, r14

	umulh	r12, r3, r7
	mul	r14, r4, r5

	adcs	r23, r23, r18
	adcs	r16, r16, r20

	mul	r18, r4, r6
	mul	r20, r4, r7

	adcs	r22, r22, r24
	cinc	r25, r25, cs
	C 9, 11, 13, 15, 17, 19, (21, 26), (23, 27, 14), (16, 10, 18), (22, 12, 20), 25

	adds	r21, r21, r26

	mul	r24, r4, r8
	umulh	r26, r4, r8

	adcs	r23, r23, r27
	adcs	r16, r16, r10

	umulh	r27, r4, r5
	umulh	r10, r4, r6

	adcs	r22, r22, r12
	cinc	r25, r25, cs
	C 9, 11, 13, 15, 17, 19, 21, (23, 14), (16, 18, 27), (22, 20, 10), (25, 24), 26

	adds	r23, r23, r14

	umulh	r12, r4, r7
	umulh	r14, r5, r7

	adcs	r16, r16, r18
	adcs	r22, r22, r20

	umulh	r18, r5, r8
	mul	r20, r5, r6

	adcs	r25, r25, r24
	cinc	r26, r26, cs
	C 9, 11, 13, 15, 17, 19, 21, 23, (16, 27), (22, 10, 20), (25, 12), (26, 14), 18

	adds	r16, r16, r27

	umulh	r24, r5, r6
	mul	r27, r5, r8

	adcs	r22, r22, r10
	adcs	r25, r25, r12

	umulh	r10, r6, r7
	umulh	r12, r6, r8

	adcs	r26, r26, r14
	cinc	r18, r18, cs
	C 9, 11, 13, 15, 17, 19, 21, 23, 16, (22, 20), (25, 24), (26, 27), (18, 10), 12

	adds	r22, r22, r20

	mul	r14, r5, r7
	mul	r20, r6, r7

	adcs	r25, r25, r24
	adcs	r26, r26, r27

	mul	r24, r6, r8
	mul	r27, r7, r8

	adcs	r18, r18, r10
	cinc	r12, r12, cs

	umulh	r10, r7, r8

	adds	r25, r25, r14

	mul	r14, r0, r1

	adcs	r26, r26, r20

	mul	r20, r0, r0
	umulh	r0, r0, r0

	adcs	r18, r18, r24
	adcs	r12, r12, r27

	mul	r24, r1, r1
	umulh	r1, r1, r1

	cinc	r10, r10, cs
	C 14, 9, 11, 13, 15, 17, 19, 21, 23, 16, 22, 25, 26, 18, 12, 10
	C 20, 0, 24, 1
	C Free: 27

	adds	r14, r14, r14
	adcs	r9, r9, r9
	adcs	r11, r11, r11
	adcs	r13, r13, r13
	adcs	r15, r15, r15
	adcs	r17, r17, r17
	adcs	r19, r19, r19
	adcs	r21, r21, r21
	lsr	r27, r10, #63
	adcs	r23, r23, r23
	adcs	r16, r16, r16
	adcs	r22, r22, r22
	adcs	r25, r25, r25
	adcs	r26, r26, r26
	adcs	r18, r18, r18
	adcs	r12, r12, r12
	adc	r10, r10, r10
	C     14,  9, 11, 13, 15, 17, 19, 21, 23, 16, 22, 25, 26, 18, 12, 10, 27
	C 20,  0, 24,  1
	C Free:

	adds	r0, r0, r14

	mul	r14, r2, r2
	umulh	r2, r2, r2

	adcs	r24, r24, r9
	adcs	r1, r1, r11

	mul	r9, r3, r3
	umulh	r3, r3, r3

	stp	r20, r0, [rp]

	mul	r11, r4, r4
	umulh	r4, r4, r4

	stp	r24, r1, [rp, #2*8]

	adcs	r14, r14, r13
	adcs	r2, r2, r15

	mul	r0, r5, r5
	umulh	r5, r5, r5

	mul	r1, r6, r6
	umulh	r6, r6, r6

	adcs	r9, r9, r17
	adcs	r3, r3, r19

	stp	r14, r2, [rp, #4*8]

	mul	r13, r7, r7
	umulh	r7, r7, r7

	adcs	r11, r11, r21
	adcs	r4, r4, r23

	ldp	r19, r20, [sp, #2*8]
	stp	r9, r3, [rp, #6*8]

	mul	r15, r8, r8
	umulh	r8, r8, r8

	adcs	r0, r0, r16
	adcs	r5, r5, r22

	ldp	r21, r22, [sp, #4*8]

	adcs	r1, r1, r25
	adcs	r6, r6, r26

	stp	r11, r4, [rp, #8*8]
	ldp	r23, r24, [sp, #6*8]

	adcs	r13, r13, r18
	adcs	r7, r7, r12

	stp	r0, r5, [rp, #10*8]
	ldp	r25, r26, [sp, #8*8]

	adcs	r15, r15, r10
	adc	r8, r8, r27

	stp	r1, r6, [rp, #12*8]
	ldr	r27, [sp, #10*8]
	stp	r13, r7, [rp, #14*8]
	ldp	r17, r18, [sp], #12*8
	stp	r15, r8, [rp, #16*8]

	mov	res, r8

	ret
EPILOGUE()
