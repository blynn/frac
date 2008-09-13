// Compute bihomographic functions of two continued fractions.
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cf.h"

struct bihom_data_s {
  cf_t x, y;
  mpz_t a[8];
};
typedef struct bihom_data_s bihom_data_t[1];
typedef struct bihom_data_s *bihom_data_ptr;

// In the 3D table, the four convergents are associated with the letter:
//   s  q
//   r  p
// e.g. the s convergent is s0/s1.
// I work left to right, and top to bottom, unlike Gosper.
struct pqrs_s {
  mpz_t p0, p1;
  mpz_t q0, q1;
  mpz_t r0, r1;
  mpz_t s0, s1;
};
typedef struct pqrs_s pqrs_t[1];

void pqrs_init(pqrs_t p) {
  mpz_init(p->p0); mpz_init(p->p1);
  mpz_init(p->q0); mpz_init(p->q1);
  mpz_init(p->r0); mpz_init(p->r1);
  mpz_init(p->s0); mpz_init(p->s1);
}

void pqrs_clear(pqrs_t p) {
  mpz_clear(p->p0); mpz_clear(p->p1);
  mpz_clear(p->q0); mpz_clear(p->q1);
  mpz_clear(p->r0); mpz_clear(p->r1);
  mpz_clear(p->s0); mpz_clear(p->s1);
}

void pqrs_set_coeff(pqrs_t p, mpz_t a[8]) {
  mpz_set(p->p0, a[0]); mpz_set(p->q0, a[1]);
  mpz_set(p->r0, a[2]); mpz_set(p->s0, a[3]);
  mpz_set(p->p1, a[4]); mpz_set(p->q1, a[5]);
  mpz_set(p->r1, a[6]); mpz_set(p->s1, a[7]);
}

void pqrs_print(pqrs_t p) {
  gmp_printf("%Zd/%Zd %Zd/%Zd\n", p->s0, p->s1, p->q0, p->q1);
  gmp_printf("%Zd/%Zd %Zd/%Zd\n", p->r0, p->r1, p->p0, p->p1);
}

static void *bihom(cf_t cf) {
  bihom_data_ptr bd = cf_data(cf);
  pqrs_t p;
  pqrs_init(p);
  pqrs_t qr;  // For quotient and remainders.
  pqrs_init(qr);
  pqrs_set_coeff(p, bd->a);
  cf_t x = bd->x;
  cf_t y = bd->y;
  mpz_t z, t0, t1;
  mpz_init(z);
  mpz_init(t0); mpz_init(t1);
  void move_down() {
    cf_get(z, y);
    mpz_mul(t0, z, p->r0);  mpz_mul(t1, z, p->r1);
    mpz_add(t0, t0, p->s0); mpz_add(t1, t1, p->s1);
    mpz_set(p->s0, p->r0);  mpz_set(p->s1, p->r1);
    mpz_set(p->r0, t0);     mpz_set(p->r1, t1);

    mpz_mul(t0, z, p->p0);  mpz_mul(t1, z, p->p1);
    mpz_add(t0, t0, p->q0); mpz_add(t1, t1, p->q1);
    mpz_set(p->q0, p->p0);  mpz_set(p->q1, p->p1);
    mpz_set(p->p0, t0);     mpz_set(p->p1, t1);
  }
  void move_right() {
    cf_get(z, x);
    mpz_mul(t0, z, p->q0);  mpz_mul(t1, z, p->q1);
    mpz_add(t0, t0, p->s0); mpz_add(t1, t1, p->s1);
    mpz_set(p->s0, p->q0);  mpz_set(p->s1, p->q1);
    mpz_set(p->q0, t0);     mpz_set(p->q1, t1);

    mpz_mul(t0, z, p->p0);  mpz_mul(t1, z, p->p1);
    mpz_add(t0, t0, p->r0); mpz_add(t1, t1, p->r1);
    mpz_set(p->r0, p->p0);  mpz_set(p->r1, p->p1);
    mpz_set(p->p0, t0);     mpz_set(p->p1, t1);
  }
  int recur() {
    if (!mpz_sgn(p->s1)) {
      move_right();
      move_down();
      return 0;
    }
    if (!mpz_sgn(p->p1)) {
      if (!mpz_sgn(p->q1)) {
	move_down();
      } else {
	move_right();
      }
      return 0;
    }
    if (!mpz_sgn(p->q1)) {
      move_down();
      return 0;
    }
    if (!mpz_sgn(p->r1)) {
      move_right();
      return 0;
    }
    mpz_fdiv_qr(qr->p0, qr->p1, p->p0, p->p1);
    // TODO: Optimize by removing other divisions.
    mpz_fdiv_qr(qr->q0, qr->q1, p->q0, p->q1);
    if (mpz_cmp(qr->p0, qr->q0)) {
      move_down();
      return 0;
    }
    mpz_fdiv_qr(qr->r0, qr->r1, p->r0, p->r1);
    if (mpz_cmp(qr->p0, qr->r0)) {
      move_right();
      return 0;
    }
    mpz_fdiv_qr(qr->s0, qr->s1, p->s0, p->s1);
    if (mpz_cmp(qr->q0, qr->q0)) {
      move_down();
      return 0;
    }
    cf_put(cf, qr->p0);
    mpz_set(p->p0, p->p1); mpz_set(p->p1, qr->p1);
    mpz_set(p->q0, p->q1); mpz_set(p->q1, qr->q1);
    mpz_set(p->r0, p->r1); mpz_set(p->r1, qr->r1);
    mpz_set(p->s0, p->s1); mpz_set(p->s1, qr->s1);
    return 1;
  }
  while(cf_wait(cf)) {
    while(!recur());
  }
  pqrs_clear(p);
  pqrs_clear(qr);
  mpz_clear(z);
  mpz_clear(t0); mpz_clear(t1);
  return NULL;
}

cf_t cf_new_bihom(cf_t x, cf_t y, mpz_t a[8]) {
  bihom_data_ptr p = malloc(sizeof(*p));
  p->x = x;
  p->y = y;
  for (int i = 0; i < 8; i++) {
    mpz_init(p->a[i]);
    mpz_set(p->a[i], a[i]);
  }
  return cf_new(bihom, p);
}
