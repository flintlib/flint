#
#   Copyright (C) 2023 Albin Ahlbäck
#
#   This file is part of FLINT.
#
#   FLINT is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Lesser General Public License (LGPL) as published
#   by the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#

include(`config.m4')dnl
dnl
.text

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

.global	FUNC(flint_mpn_mul_16_2)
.p2align	4, 0x90
TYPE(flint_mpn_mul_16_2)

FUNC(flint_mpn_mul_16_2):
	.cfi_startproc
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
.flint_mpn_mul_16_2_end:
SIZE(flint_mpn_mul_16_2, .flint_mpn_mul_16_2_end)
.cfi_endproc
