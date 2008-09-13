#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "cf.h"

int main() {
  mpz_t z;
  mpz_init(z);
  cf_t pi, conv;
  pi = cf_new_pi();

  conv = cf_new_cf_to_decimal(pi);
  for (int i = 0; i <= 5000; i++) {
    cf_get(z, conv);
    gmp_printf("%Zd", z);
    if (!(i % 5)) putchar(' ');
    if (!(i % 50)) putchar('\n');
  }
  cf_free(conv);
  cf_free(pi);
  mpz_clear(z);
  return 0;
}
