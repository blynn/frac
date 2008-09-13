#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cf.h"

// Converges extremely slowly. Not worth converting to decimal or
// continued fraction form.
static void *pi_six_sequence(cf_t cf) {
  mpz_t num, denom, t;
  mpz_init(num);
  mpz_init(denom);
  mpz_init(t);

  mpz_set_ui(denom, 3);
  cf_put(cf, denom);

  mpz_set_ui(num, 1);
  cf_put(cf, num);
  mpz_set_ui(t, 8);

  mpz_set_ui(denom, 6);
  while(cf_wait(cf)) {
    cf_put(cf, denom);
    mpz_add(num, num, t);
    cf_put(cf, num);
    mpz_add_ui(t, t, 8);
  }

  mpz_clear(num);
  mpz_clear(denom);
  mpz_clear(t);
  return NULL;
}

int main() {
  mpz_t z;
  mpz_t a, b, c, d;
  mpz_init(z);
  cf_t pi, conv;
  pi = cf_new(pi_six_sequence, NULL);
  mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d);
  mpz_set_ui(a, 1); mpz_set_ui(b, 0);
  mpz_set_ui(c, 0); mpz_set_ui(d, 1);

  conv = cf_new_nonregular_mobius_to_decimal(pi, a, b, c, d);
  for (int i = 1; i <= 5000; i++) {
    cf_get(z, conv);
    gmp_printf("%Zd\n",z);
  }
  printf("\n");
  cf_free(conv);
  cf_free(pi);
  mpz_clear(z);
  mpz_clear(d); mpz_clear(c); mpz_clear(b); mpz_clear(a);
  return 0;
}
