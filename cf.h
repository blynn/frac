// Requires gmp.h
//
// Opaque interface to continued fractions object.

#ifndef __CF_H__
#define __CF_H__

struct cf_s;
typedef struct cf_s *cf_t;

cf_t cf_new(void *(*func)(cf_t), void *data);
void cf_free(cf_t cf);

void cf_get(mpz_t z, cf_t cf);
void cf_put(cf_t cf, mpz_t z);

void cf_signal(cf_t cf);
int cf_wait(cf_t cf);

void *cf_data(cf_t cf);

// cf_converge.c:
// Compute convergents of a simple continued fraction.
// Outputs p then q on channel, where p/q is the last convergent computed.
cf_t cf_new_convergent(cf_t a);

#endif  // __CF_H__
