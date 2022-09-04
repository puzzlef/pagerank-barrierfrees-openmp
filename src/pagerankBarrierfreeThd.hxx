#pragma once
#include <thread>
#include <atomic>
#include <algorithm>
#include <vector>
#include "_main.hxx"
#include "transpose.hxx"
#include "dynamic.hxx"
#include "pagerank.hxx"
#include "pagerankSeq.hxx"
#include "pagerankOmp.hxx"
#include "pagerankMonolithicSeq.hxx"

using std::atomic;
using std::thread;
using std::vector;
using std::swap;
using std::min;




// PAGERANK-LOOP
// -------------

template <bool O, bool D, bool F, class T>
void pagerankBarrierfreeThdFuncU(int *l, vector<T> *a, vector<T> *r, vector<T> *c, const vector<T> *f, const vector<int> *vfrom, const vector<int> *efrom, const vector<int> *vdata, int i, int n, int N, T p, T E, int L, int EF) {
  *l = pagerankMonolithicSeqLoopU<O, D, T, F>(*a, *r, *c, *f, *vfrom, *efrom, *vdata, i, n, N, p, E, L, EF);
}

template <bool O, bool D, bool F, class T>
int pagerankBarrierfreeThdLoopU(vector<T>& a, vector<T>& r, vector<T>& c, const vector<T>& f, const vector<int>& vfrom, const vector<int>& efrom, const vector<int>& vdata, int i, int n, int N, T p, T E, int L, int EF) {
  float l = 0;
  if (!O) return 0;
  // Ordered approach
  int TS = omp_get_max_threads();
  int DN = ceilDiv(n, TS);
  vector<thread> threads;
  vector<int>    l1(TS);
  for (int t=0; t<TS; t++) {
    int i1 = i+t*DN, I1 = min(i1+DN, i+n), n1 = I1-i1;
    threads.emplace_back(pagerankBarrierfreeThdFuncU<O, D, F, T>, &l1[t], &a, &r, &c, &f, &vfrom, &efrom, &vdata, i1, n1, N, p, E, L, EF);
  }
  for (int t=0; t<TS; t++) {
    int i1 = i+t*DN, I1 = min(i1+DN, i+n), n1 = I1-i1;
    threads[t].join();
    l += l1[t] * n1/float(n);
  }
  return int(l + 0.5f);
}




// PAGERANK (STATIC / INCREMENTAL)
// -------------------------------

// Find pagerank using multiple threads (pull, CSR).
// @param x  original graph
// @param xt transpose graph (with vertex-data=out-degree)
// @param q  initial ranks (optional)
// @param o  options {damping=0.85, tolerance=1e-6, maxIterations=500}
// @returns {ranks, iterations, time}
template <bool O, bool D, bool F, class G, class H, class T=float>
PagerankResult<T> pagerankBarrierfreeThd(const G& x, const H& xt, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}) {
  int  N  = xt.order();  if (N==0) return PagerankResult<T>::initial(xt, q);
  auto ks = vertexKeys(xt);
  return pagerankOmp(xt, ks, 0, N, pagerankBarrierfreeThdLoopU<O, D, F, T>, q, o);
}

template <bool O, bool D, bool F, class G, class T=float>
PagerankResult<T> pagerankBarrierfreeThd(const G& x, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}) {
  auto xt = transposeWithDegree(x);
  return pagerankBarrierfreeThd<O, D, F>(x, xt, q, o);
}




// PAGERANK (DYNAMIC)
// ------------------

template <bool O, bool D, bool F, class G, class H, class T=float>
PagerankResult<T> pagerankBarrierfreeThdDynamic(const G& x, const H& xt, const G& y, const H& yt, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}) {
  int  N = yt.order();                             if (N==0) return PagerankResult<T>::initial(yt, q);
  auto [ks, n] = dynamicInVertices(x, xt, y, yt);  if (n==0) return PagerankResult<T>::initial(yt, q);
  return pagerankOmp(yt, ks, 0, n, pagerankBarrierfreeThdLoopU<O, D, F, T>, q, o);
}

template <bool O, bool D, bool F, class G, class T=float>
PagerankResult<T> pagerankBarrierfreeThdDynamic(const G& x, const G& y, const vector<T> *q=nullptr, const PagerankOptions<T>& o={}) {
  auto xt = transposeWithDegree(x);
  auto yt = transposeWithDegree(y);
  return pagerankBarrierfreeThdDynamic<O, D, F>(x, xt, y, yt, q, o);
}
