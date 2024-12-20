#include "ng_inner.h"

/*
 * Table format: first entry is kmax; it is followed by 2*kmax values
 * (formal table entries for -kmax to +kmax-1).
 */

/*
 * Type for a PRNG. When invoked, it writes len bytes into the buffer
 * pointed to by dst; the context parameter points to the PRNG state
 * structure. This is a wrapper around the externally provided RNG,
 * so that the external RNG is invoked only for producing chunks of
 * 512 bytes.
 */
typedef struct {
	uint8_t buf[512];
	size_t ptr;
	ntrugen_rng rng;
	void *rng_context;
} prng_buffer;

static inline void
prng_buffer_init(prng_buffer *pb, ntrugen_rng rng, void *rng_context)
{
	pb->ptr = sizeof pb->buf;
	pb->rng = rng;
	pb->rng_context = rng_context;
}

static inline uint16_t
prng_buffer_next_u16(prng_buffer *pb)
{
	if (pb->ptr > (sizeof pb->buf) - 2) {
		pb->rng(pb->rng_context, pb->buf, sizeof pb->buf);
		pb->ptr = 0;
	}
	unsigned x = pb->buf[pb->ptr ++];
	x |= (unsigned)pb->buf[pb->ptr ++] << 8;
	return x;
}



/* see ng_inner.h */
TARGET_AVX2
void
gauss_sample_poly(unsigned logn, int8_t *f,
	const uint16_t *tab, ntrugen_rng rng, void *rng_context)
{
	size_t n = (size_t)1 << logn;
	size_t kmax = tab[0];
	prng_buffer pb;
	prng_buffer_init(&pb, rng, rng_context);
	for (;;) {
		uint32_t parity = 0;
		for (size_t j = 0; j < n; j ++) {
			uint32_t v = -(uint32_t)kmax;
			uint32_t x = prng_buffer_next_u16(&pb);
			for (size_t k = 1; k <= (kmax << 1); k ++) {
				v += ((uint32_t)tab[k] - x) >> 31;
			}
			f[j] = (int8_t)*(int32_t *)&v;
			parity ^= v;
		}
		if ((parity & 1) != 0) {
			return;
		}
	}
}

/* see ng_inner.h */
void
gauss_sample_poly_reduced(unsigned logn, int8_t *f,
	const uint16_t *tab, ntrugen_rng rng, void *rng_context)
{
	size_t n = (size_t)1 << logn;
	int g = 1 << (8 - logn);
	size_t kmax = tab[0];
	prng_buffer pb;
	prng_buffer_init(&pb, rng, rng_context);
	for (;;) {
		uint32_t parity = 0;
		for (size_t j = 0; j < n;) {
			uint32_t v = -(uint32_t)kmax << (8 - logn);
			for (int i = 0; i < g; i ++) {
				uint32_t x = prng_buffer_next_u16(&pb);
				for (size_t k = 1; k <= (kmax << 1); k ++) {
					uint32_t z = tab[k];
					v += (z - x) >> 31;
				}
			}
			int32_t y = *(int32_t *)&v;
			if (y < -127 || y > +127) {
				continue;
			}
			f[j ++] = (int8_t)y;
			parity ^= v;
		}
		if ((parity & 1) != 0) {
			return;
		}
	}
}
