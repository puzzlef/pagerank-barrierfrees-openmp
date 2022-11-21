#pragma once
#include <vector>
#include <algorithm>
#include <omp.h>
#include "_main.hxx"
#include "transpose.hxx"
#include "dynamic.hxx"
#include "pagerank.hxx"
#include "pagerankSeq.hxx"
#include "pagerankOmp.hxx"

using std::vector;
using std::swap;




// PAGERANK-LOOP
// -------------

template <bool DEAD=false, bool SLEEP=false, bool HELP=false, bool ONE=false, class K, class T>
int pagerankMonolithicBarrierfreeOmpLoopU(vector<T>& a, vector<T>& r, vector<T>& c, const vector<T>& f, const vector<K>& vfrom, const vector<K>& efrom, const vector<K>& vdata, vector<PagerankThreadWork*>& works, K i, K n, K N, T p, T E, int L, int EF, float SP, int SD) {
  PRINTFI("pagerankMonolithicBarrierfreeOmpLoopU():\n");
  atomic<int> ls = 0;
  // 0. Reset thread works.
  for (int t=0; t<works.size(); ++t)
    works[t]->clear();
  // 1. Perform iterations.
  auto tstart = timeNow();
  PRINTFT("[%09.3f ms] parallel_out_begin\n", durationMilliseconds(tstart));
  #pragma omp parallel
  {
    int l = 0;
    vector<T> *rp = &r;
    vector<T> *ap = ONE? &r : &a;
    int t = omp_get_thread_num();
    PagerankThreadWork& me = *works[t];
    PRINTFT("[%09.3f ms; thread %02d] parallel_begin\n",   durationMilliseconds(tstart), t);
    PRINTFT("[%09.3f ms; thread %02d] iterations_begin\n", durationMilliseconds(tstart), t);
    while (l<L) {
      T c0 = DEAD? pagerankTeleportBarrierfreeOmp(*rp, vdata, N, p) : (1-p)/N;
      if (HELP) pagerankCalculateHelperBarrierfreeOmpW<SLEEP>(*ap, *rp, f, vfrom, efrom, i, n, c0, SP, SD, works, tstart);
      else      pagerankCalculateBarrierfreeOmpW<SLEEP>      (*ap, *rp, f, vfrom, efrom, i, n, c0, SP, SD, works, tstart);  // update ranks of vertices
      T el = pagerankErrorBarrierfreeOmp(works);  // compare previous and current ranks
      swap(ap, rp); ++l;                          // final ranks in (r)
      if (el<E) break;                            // check tolerance
    }
    if (l>ls) ls = l;
    float msparend = works[t]->blockedStartTime = durationMilliseconds(tstart);
    PRINTFT("[%09.3f ms; thread %02d] iterations_end\n", msparend, t);
    PRINTFT("[%09.3f ms; thread %02d] parallel_end\n",   msparend, t);
  }
  float msparoutend = durationMilliseconds(tstart);
  PRINTFT("[%09.3f ms] parallel_out_end\n", msparoutend);
  for (int t=0; t<works.size(); ++t)
    works[t]->blockedTime += msparoutend - works[t]->blockedStartTime;
  for (int t=0; t<works.size(); ++t)
    PRINTFI("[thread %02d] status {processed=%zu, stolen=%zu, slept=%zu, blocked=%.3f ms}\n", t, works[t]->processedCount, works[t]->stolenCount, works[t]->sleptCount, works[t]->blockedTime);
  if (!ONE && (ls & 1)==1) swap(a, r);
  return ls;
}




// PAGERANK (STATIC / INCREMENTAL)
// -------------------------------

// Find pagerank using multiple threads (pull, CSR).
// @param x  original graph
// @param xt transpose graph (with vertex-data=out-degree)
// @param q  initial ranks (optional)
// @param o  options {damping=0.85, tolerance=1e-6, maxIterations=500}
// @returns {ranks, iterations, time}
template <bool DEAD=false, bool SLEEP=false, bool HELP=false, bool ONE=false, class G, class H, class T=float>
PagerankResult<T> pagerankMonolithicBarrierfreeOmp(const G& x, const H& xt, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  using K = typename G::key_type;
  K    N  = xt.order();  if (N==0) return PagerankResult<T>::initial(xt, q);
  auto ks = pagerankVertices(x, xt, o, C);
  return pagerankOmp(xt, ks, K(0), N, pagerankMonolithicBarrierfreeOmpLoopU<DEAD, SLEEP, HELP, ONE, K, T>, q, o);
}

template <bool DEAD=false, bool SLEEP=false, bool HELP=false, bool ONE=false, class G, class T=float>
PagerankResult<T> pagerankMonolithicBarrierfreeOmp(const G& x, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  auto xt = transposeWithDegree(x);
  return pagerankMonolithicBarrierfreeOmp<DEAD, SLEEP, HELP, ONE>(x, xt, q, o, C);
}




// PAGERANK (DYNAMIC)
// ------------------

template <bool DEAD=false, bool SLEEP=false, bool HELP=false, bool ONE=false, class G, class H, class T=float>
PagerankResult<T> pagerankMonolithicBarrierfreeOmpDynamic(const G& x, const H& xt, const G& y, const H& yt, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  using K = typename G::key_type;
  K     N = yt.order();                                        if (N==0) return PagerankResult<T>::initial(yt, q);
  auto [ks, n] = pagerankDynamicVertices(x, xt, y, yt, o, C);  if (n==0) return PagerankResult<T>::initial(yt, q);
  return pagerankOmp(yt, ks, K(0), n, pagerankMonolithicBarrierfreeOmpLoopU<DEAD, SLEEP, HELP, ONE, K, T>, q, o);
}

template <bool DEAD=false, bool SLEEP=false, bool HELP=false, bool ONE=false, class G, class T=float>
PagerankResult<T> pagerankMonolithicBarrierfreeOmpDynamic(const G& x, const G& y, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}, const PagerankData<G> *C=nullptr) {
  auto xt = transposeWithDegree(x);
  auto yt = transposeWithDegree(y);
  return pagerankMonolithicBarrierfreeOmpDynamic<DEAD, SLEEP, HELP, ONE>(x, xt, y, yt, q, o, C);
}
