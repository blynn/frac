#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "cf.h"

int main(int argc, char **argv) {
  mpz_t z;
  mpz_init(z);
  cf_t pi, conv;
  pi = cf_new_pi();
  int n = 2037 + 1;  // To outdo Metropolis, Reitwieser and von Neumann's
                     // 1949 ENIAC record.
  if (argc > 1) {
    n = atoi(argv[1]);
    if (n <= 0) n = 100;
  }

  conv = cf_new_cf_to_decimal(pi);
  for (int i = 0; i <= n; i++) {
    cf_get(z, conv);
    gmp_printf("%Zd", z);
    if (!(i % 5)) putchar(' ');
    if (!(i % 50)) putchar('\n');
  }
  if (n % 50) putchar('\n');
  cf_free(conv);
  cf_free(pi);
  mpz_clear(z);
  return 0;
}
