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

int cf_wait(cf_t cf);

void *cf_data(cf_t cf);

// cf_mobius.c:
// Compute convergents of a simple continued fraction x.
// Outputs p then q on channel, where p/q is the last convergent computed.
cf_t cf_new_convergent(cf_t x);

// Compute convergents of (a x + b)/(c x + d)
// where x is a regular continued fraction.
cf_t cf_new_mobius_convergent(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d);

// Compute convergents of (a x + b)/(c x + d)
// where x is a nonregular continued fraction.
cf_t cf_new_nonregular_mobius_convergent(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d);

cf_t cf_new_nonregular_to_cf(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d);

cf_t cf_new_mobius_to_decimal(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d);
cf_t cf_new_cf_to_decimal(cf_t x);
cf_t cf_new_nonregular_mobius_to_decimal(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d);

// cf_famous.c:
// e:
cf_t cf_new_e();
cf_t cf_new_pi();

#endif  // __CF_H__
