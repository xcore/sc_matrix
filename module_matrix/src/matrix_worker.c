/*
 * Matrix op worker threads
 *
 * Copyright (C) 2011 Steve Kerrison <github@stevekerrison.com>
 *
 * This software is freely distributable under a derivative of the
 * University of Illinois/NCSA Open Source License posted in
 * LICENSE.txt and at <http://github.xcore.com/>
 */

#include "matrix_worker.h"
#include <stdlib.h>

/* If rand() isn't 32-bit let's fix it. */
#if RAND_MAX < 0xffff
#define myrand() (((rand() << 24) & 0xff000000) | ((rand() << 16) & 0x00ff000)) \
	((rand() << 8) & 0x0000ff00)) | (rand() & 0x000000ff))
#elif RAND_MAX < 0xffffffff
#define myrand() (((rand() << 16) & 0xffff0000) | (rand() & 0xffff))
#else
#define myrand() rand()
#endif

void matrix_mul_worker(int ptA, int ptDimA, int ptB, int ptDimB, int ptC,
	int ptOps, char nThreads, char offset)
{
	int *A = (int *)ptA, *B = (int *)ptB, *C = (int *)ptC,
		*ops = (int *)ptOps;
	short *dimA = (short *)ptDimA, *dimB = (short *)ptDimB;
	int p,r,c;
	ops += offset;
	*ops = 0;
	int rlim = dimA[0], clim = dimB[1], plim = dimA[1];
	for (r = 0; r < rlim; r++)
	{
		for (c = offset; c < clim; c += nThreads)
		{
			C[r * rlim + c] = 0;
			for (p = 0; p < plim; p += 1)
			{
				C[r * rlim + c] += A[r * plim + p] * B[p * clim + c];
				*ops += 1;
			}
		}
	}
	return;
}

void matrix_arr_worker(enum matrix_ops op, int ptA, int ptDimA, int ptB, int ptDimB, int ptC,
	int ptOps, char nThreads, char offset)
{
	int *A = (int *)ptA, *B = (int *)ptB, *C = (int *)ptC,
		*ops = (int *)ptOps;
	short *dimA = (short *)ptDimA, *dimB = (short *)ptDimB;
	int r,c;
	ops += offset;
	*ops = 0;
	int rlim = dimA[0], clim = dimB[1];
	for (r = 0; r < rlim; r++)
	{
		for (c = offset; c < clim; c += nThreads)
		{
			int res, a = A[r * rlim + c], b = B[r * rlim + c];
			switch (op)
			{
			case ADD:
				res = a + b;
				break;
			case SUB:
				res = a - b;
				break;
			case MUL:
				res = a * b;
				break;
			case DIV:
			case SDIV:
				res = a / b;
				break;
			case UDIV:
				res = (unsigned)a / (unsigned)b;
				break;
			case RAND: //Fall through to default
			default:
				break;	
			}
			C[r * rlim + c] = res;
			*ops += 1;
		}
	}
	return;
}

void matrix_sca_worker(enum matrix_ops op, int ptA, int ptDimA, int S, int ptC,
	int ptOps, char nThreads, char offset)
{
	int *A = (int *)ptA, *C = (int *)ptC,
		*ops = (int *)ptOps;
	short *dimA = (short *)ptDimA;
	int r,c;
	ops += offset;
	*ops = 0;
	int rlim = dimA[0], clim = dimA[1];
	for (r = 0; r < rlim; r++)
	{
		for (c = offset; c < clim; c += nThreads)
		{
			int a = A[r * rlim + c], res;
			switch (op)
			{
			case ADD:
				res = a + S;
				break;
			case SUB:
				res = a - S;
				break;
			case MUL:
				res = a * S;
				break;
			case DIV:
			case SDIV:
				res = a / S;
				break;
			case UDIV:
				res = (unsigned)a / (unsigned)S;
				break;
			case RAND:
				res = myrand();
				break;
			default:
				break;	
			}
			C[r * rlim + c] = res;
			*ops += 1;
		}
	}
	return;
}
