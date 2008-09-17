// Requires gmp.h
//
// Opaque interface to continued fractions object.

#ifndef __CF_H__
#define __CF_H__

struct cf_s;
typedef struct cf_s *cf_t;

cf_t cf_new(void *(*func)(cf_t), void *data);
static inline cf_t cf_new_const(void *(*func)(cf_t)) {
  return cf_new(func, NULL);
}
void cf_free(cf_t cf);

void cf_set_sign(cf_t cf, int sign);
int cf_sign(cf_t cf);
int cf_flip_sign(cf_t cf);
void cf_get(mpz_t z, cf_t cf);
void cf_put(cf_t cf, mpz_t z);
void cf_put_int(cf_t cf, int n);

int cf_wait(cf_t cf);

void *cf_data(cf_t cf);

void cf_signal(cf_t cf); // For tee.
void cf_wait_special(cf_t cf);

// From cf_tee.c:
//
void cf_tee(cf_t *out_array, cf_t in);

// From cf_mobius.c:
//
// Compute convergents of a simple continued fraction x.
// Outputs p then q on channel, where p/q is the last convergent computed.
cf_t cf_new_cf_convergent(cf_t x);
// Compute decimal representation of a simple continued fraction x.
// Outputs integer part first, then digits one at a time.
cf_t cf_new_cf_to_decimal(cf_t x);

// Compute convergents of (a x + b)/(c x + d)
// where x is a regular continued fraction.
cf_t cf_new_mobius_convergent(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d);
cf_t cf_new_mobius_to_decimal(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d);
cf_t cf_new_mobius_to_cf(cf_t x, mpz_t z[4]);

// Compute convergents of (a x + b)/(c x + d)
// where x is a nonregular continued fraction.
cf_t cf_new_nonregular_mobius_convergent(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d);
// Input: Mobius transformation and nonregular continued fraction.
// Output: Regular continued fraction. Assumes input fraction is well-behaved.
cf_t cf_new_nonregular_to_cf(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d);
// Does both of the above at once. Seems slow.
cf_t cf_new_nonregular_mobius_to_decimal(cf_t x, mpz_t a[4]);

// Well-known continued fraction expansions.
// cf_famous.c:
// e:
cf_t cf_new_sqrt2();
cf_t cf_new_e();
cf_t cf_new_pi();
cf_t cf_new_tan1();
cf_t cf_new_epow(mpz_t pow);
cf_t cf_new_tanh(mpz_t z);

// This won't work because my code cannot handle negative denominators,
// and also assumes the sequence of convergents alternatively overshoot
// and undershoots the target. The tan expansion leads to a sequence of
// strictly increasing convergents (for positive input).
cf_t cf_new_tan(mpz_t z);

// Gosper's method for computing bihomographic functions of continued fractions.
cf_t cf_new_bihom(cf_t x, cf_t y, mpz_t a[8]);

// From taylor.c:
cf_t cf_new_sin1();
cf_t cf_new_cos1();

// From newton.c:
// Use Newton's method to find solutions of:
//
//     a0 xy + a1 x + a2 y + a3
// y = ------------------------
//     a4 xy - a0 x + a5 y - a2
cf_t cf_new_newton(cf_t x, mpz_t a[6], mpz_t lower);
cf_t cf_new_sqrt(cf_t x);

#endif  // __CF_H__
