#include <string>
#include <vector>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include "src/main.hxx"

using namespace std;




// You can define datatype with -DTYPE=...
#ifndef TYPE
#define TYPE float
#endif
// You can define number of threads with -DMAX_THREADS=...
#ifndef MAX_THREADS
#define MAX_THREADS 12
#endif




template <class G, class H>
void runPagerank(const G& x, const H& xt, int repeat) {
  using T = TYPE;
  enum NormFunction { L0=0, L1=1, L2=2, Li=3 };
  vector<T> *init = nullptr;
  float damping   = 0.85;
  float tolerance = 1e-10;

  // Find pagerank using a single thread (unordered, no dead ends).
  auto a1 = pagerankMonolithicSeq<false, false>(x, xt, init, {repeat, Li, damping, tolerance});
  auto e1 = l1Norm(a1.ranks, a1.ranks);
  printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankSeqUnordered\n", a1.time, a1.iterations, e1);

  // Find pagerank using a single thread (ordered, no dead ends).
  auto a2 = pagerankMonolithicSeq<true, false>(x, xt, init, {repeat, Li, damping, tolerance});
  auto e2 = l1Norm(a2.ranks, a1.ranks);
  printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankSeqOrdered\n", a2.time, a2.iterations, e2);

  // Find pagerank accelerated with OpenMP (ordered, no dead ends).
  auto a3 = pagerankMonolithicOmp<true, false>(x, xt, init, {repeat, Li, damping, tolerance});
  auto e3 = l1Norm(a3.ranks, a1.ranks);
  printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankOmpOrdered\n", a3.time, a3.iterations, e3);

  // Find pagerank with barrier-free iterations accelerated with OpenMP (ordered, no dead ends).
  auto a4 = pagerankBarrierfreeOmp<true, false>(x, xt, init, {repeat, Li, damping, tolerance});
  auto e4 = l1Norm(a4.ranks, a1.ranks);
  printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreeOmpOrdered\n", a4.time, a4.iterations, e4);
}


int main(int argc, char **argv) {
  char *file = argv[1];
  int repeat = argc>2? stoi(argv[2]) : 5;
  printf("Loading graph %s ...\n", file);
  auto x  = readMtxOutDiGraph(file); println(x);
  auto fl = [](auto u) { return true; };
  selfLoopU(x, None(), fl); print(x); printf(" (selfLoopAllVertices)\n");
  auto xt = transposeWithDegree(x);  print(xt); printf(" (transposeWithDegree)\n");
  omp_set_num_threads(MAX_THREADS);
  printf("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  runPagerank(x, xt, repeat);
  printf("\n");
  return 0;
}
