#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "cf.h"

// Minds our p's and q's. The two last computed convergents,
// and room to compute the next one.
struct pqset_s {
  mpz_t pold, p, pnew;
  mpz_t qold, q, qnew;
};
typedef struct pqset_s *pqset_ptr;
typedef struct pqset_s pqset_t[1];

void pqset_init(pqset_t pq) {
  mpz_init(pq->pold);
  mpz_init(pq->qold);
  mpz_init(pq->p);
  mpz_init(pq->q);
  mpz_init(pq->pnew);
  mpz_init(pq->qnew);
}

void pqset_clear(pqset_t pq) {
  mpz_clear(pq->pold);
  mpz_clear(pq->qold);
  mpz_clear(pq->p);
  mpz_clear(pq->q);
  mpz_clear(pq->pnew);
  mpz_clear(pq->qnew);
}

// Compute the next convergent for regular continued fractions.
void pqset_regular_recur(pqset_t pq, mpz_t denom) {
  mpz_mul(pq->pnew, pq->p, denom);
  mpz_add(pq->pnew, pq->pnew, pq->pold);
  mpz_mul(pq->qnew, pq->q, denom);
  mpz_add(pq->qnew, pq->qnew, pq->qold);
}

// Compute the next convergent for nonregular continued fractions.
void pqset_nonregular_recur(pqset_t pq, mpz_t num, mpz_t denom) {
  mpz_mul(pq->pnew, pq->p, denom);
  mpz_addmul(pq->pnew, pq->pold, num);
  mpz_mul(pq->qnew, pq->q, denom);
  mpz_addmul(pq->qnew, pq->qold, num);
}

// Use pold, qold as temporary variables.
// Get rid of nontrivial GCD for {p, q, pnew, qnew}.
void pqset_remove_gcd(pqset_ptr pq) {
  mpz_gcd(pq->pold, pq->pnew, pq->qnew);
  mpz_gcd(pq->qold, pq->p, pq->q);
  mpz_gcd(pq->pold, pq->pold, pq->qold);
  if (mpz_cmp_ui(pq->pold, 1)) {
    mpz_divexact(pq->pnew, pq->pnew, pq->pold);
    mpz_divexact(pq->qnew, pq->qnew, pq->pold);
    mpz_divexact(pq->p, pq->p, pq->pold);
    mpz_divexact(pq->q, pq->q, pq->pold);
  }
}

void pqset_next(pqset_ptr pq) {
  mpz_set(pq->pold, pq->p); mpz_set(pq->p, pq->pnew);
  mpz_set(pq->qold, pq->q); mpz_set(pq->q, pq->qnew);
}

struct mobius_data_s {
  cf_t input;
  mpz_t a, b, c, d;
};
typedef struct mobius_data_s *mobius_data_ptr;

void pqset_set_mobius(pqset_t pq, mobius_data_ptr md) {
  mpz_set(pq->pold, md->b); mpz_set(pq->p, md->a);
  mpz_set(pq->qold, md->d); mpz_set(pq->q, md->c);
}

static void *mobius(cf_t cf) {
  mobius_data_ptr md = cf_data(cf);
  cf_t input = md->input;
  pqset_t pq;
  pqset_init(pq);
  pqset_set_mobius(pq, md);

  mpz_t denom;
  mpz_init(denom);
  while(cf_wait(cf)) {
    cf_get(denom, input);
    pqset_regular_recur(pq, denom);

    cf_put(cf, pq->pnew);
    cf_put(cf, pq->qnew);
    pqset_next(pq);
  }
  mpz_clear(denom);
  pqset_clear(pq);

  mpz_clear(md->a);
  mpz_clear(md->b);
  mpz_clear(md->c);
  mpz_clear(md->d);
  free(md);
  return NULL;
}
cf_t cf_new_convergent_mobius(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d) {
  mobius_data_ptr md = malloc(sizeof(*md));
  mpz_init(md->a);
  mpz_init(md->b);
  mpz_init(md->c);
  mpz_init(md->d);
  mpz_set(md->a, a);
  mpz_set(md->b, b);
  mpz_set(md->c, c);
  mpz_set(md->d, d);
  md->input = x;
  return cf_new(mobius, md);
}

cf_t cf_new_convergent(cf_t x) {
  mpz_t one, zero;
  mpz_init(one);
  mpz_init(zero);
  mpz_set_ui(one, 1);
  mpz_set_ui(zero, 0);
  return cf_new_convergent_mobius(x, one, zero, zero, one);
  mpz_clear(one);
  mpz_clear(zero);
}

static void *mobius_nonregular(cf_t cf) {
  mobius_data_ptr md = cf_data(cf);
  cf_t input = md->input;
  pqset_t pq;
  pqset_init(pq);
  pqset_set_mobius(pq, md);
  mpz_t num;
  mpz_init(num);
  mpz_t denom;
  mpz_init(denom);
  void recur() {
    pqset_nonregular_recur(pq, num, denom);
    pqset_remove_gcd(pq);

    cf_put(cf, pq->pnew);
    cf_put(cf, pq->qnew);
    pqset_next(pq);
  }
  mpz_set_ui(num, 1);
  cf_get(denom, input);
  recur();
  while(cf_wait(cf)) {
    cf_get(num, input);
    cf_get(denom, input);
    recur();
  }
  mpz_clear(denom);
  pqset_clear(pq);

  mpz_clear(md->a);
  mpz_clear(md->b);
  mpz_clear(md->c);
  mpz_clear(md->d);
  free(md);
  return NULL;
}
cf_t cf_new_convergent_nonregular(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d) {
  mobius_data_ptr md = malloc(sizeof(*md));
  mpz_init(md->a);
  mpz_init(md->b);
  mpz_init(md->c);
  mpz_init(md->d);
  mpz_set(md->a, a);
  mpz_set(md->b, b);
  mpz_set(md->c, c);
  mpz_set(md->d, d);
  md->input = x;
  return cf_new(mobius_nonregular, md);
}

static void *mobius_nonregular_throughput(cf_t cf) {
  mobius_data_ptr md = cf_data(cf);
  cf_t input = md->input;
  pqset_t pq;
  pqset_init(pq);
  pqset_set_mobius(pq, md);
  mpz_t num;
  mpz_init(num);
  mpz_t denom;
  mpz_init(denom);
  mpz_t t0, t1;
  mpz_init(t1); mpz_init(t0);
  int recur() {
    pqset_nonregular_recur(pq, num, denom);
    pqset_remove_gcd(pq);

    mpz_fdiv_qr(pq->pold, t0, pq->p, pq->q);
    // mpz_fdiv_qr(qold, t1, pnew, qnew);
    mpz_mul(pq->qold, pq->pold, pq->qnew);

    // TODO: What if one of these equals zero?
    if (mpz_cmp(pq->qold, pq->pnew) < 0) {
      mpz_add(pq->qold, pq->qold, pq->qnew);
      //if (!mpz_cmp(pold, qold)) {
      if (mpz_cmp(pq->qold, pq->pnew) > 0) {
	cf_put(cf, pq->pold);
	mpz_sub(t1, pq->qold, pq->pnew);
	mpz_sub(t1, pq->qnew, t1);
	mpz_set(pq->pold, pq->q);
	mpz_set(pq->qold, t0);
	mpz_set(pq->p, pq->qnew);
	mpz_set(pq->q, t1);
	return 1;
      }
    }
    pqset_next(pq);
    return 0;
  }
  mpz_set_ui(num, 1);
  cf_get(denom, input);
  recur();
  while(cf_wait(cf)) {
    do {
      cf_get(num, input);
      cf_get(denom, input);
    } while(!recur());
  }
  mpz_clear(denom);
  pqset_clear(pq);

  mpz_clear(md->a);
  mpz_clear(md->b);
  mpz_clear(md->c);
  mpz_clear(md->d);
  free(md);
  mpz_clear(t1); mpz_clear(t0);
  return NULL;
}
cf_t cf_new_nonregular_to_cf(cf_t x, mpz_t a, mpz_t b, mpz_t c, mpz_t d) {
  mobius_data_ptr md = malloc(sizeof(*md));
  mpz_init(md->a);
  mpz_init(md->b);
  mpz_init(md->c);
  mpz_init(md->d);
  mpz_set(md->a, a);
  mpz_set(md->b, b);
  mpz_set(md->c, c);
  mpz_set(md->d, d);
  md->input = x;
  return cf_new(mobius_nonregular_throughput, md);
}
