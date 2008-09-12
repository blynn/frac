#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cf.h"

static void *pi_arctan_sequence(cf_t cf) {
  mpz_t num, denom;
  mpz_init(num);
  mpz_init(denom);

  mpz_set_ui(denom, 1);
  cf_put(cf, denom);

  mpz_set_ui(num, 1);
  cf_put(cf, num);

  while(cf_wait(cf)) {
    mpz_add_ui(denom, denom, 2);
    cf_put(cf, denom);
    mpz_add(num, num, denom);
    cf_put(cf, num);
  }

  mpz_clear(num);
  mpz_clear(denom);
  return NULL;
}

int main() {
  mpz_t z;
  mpz_t a, b, c, d;
  mpz_init(z);
  cf_t pi, conv;
  pi = cf_new(pi_arctan_sequence, NULL);
  mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d);
  mpz_set_ui(a, 0); mpz_set_ui(b, 4);
  mpz_set_ui(c, 1); mpz_set_ui(d, 0);

  conv = cf_new_nonregular_to_cf(pi, a, b, c, d);
  cf_t dec = cf_new_cf_to_decimal(conv);
  cf_get(z, dec);
  for (int i = 1; i <= 5000; i++) {
    cf_get(z, dec);
    gmp_printf("%Zd",z);
    if (!(i%5)) printf(" ");
    if (!(i%50)) printf("\n");
  }
  printf("\n");
  cf_free(conv);
  cf_free(pi);
  mpz_clear(z);
  mpz_clear(d); mpz_clear(c); mpz_clear(b); mpz_clear(a);
  return 0;
}
