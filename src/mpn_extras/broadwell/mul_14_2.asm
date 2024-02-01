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

.global	FUNC(flint_mpn_mul_14_2)
.p2align	4, 0x90
TYPE(flint_mpn_mul_14_2)

FUNC(flint_mpn_mul_14_2):
	.cfi_startproc
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
.flint_mpn_mul_14_2_end:
SIZE(flint_mpn_mul_14_2, .flint_mpn_mul_14_2_end)
.cfi_endproc
