#include <algorithm>
#include <string>
#include <vector>
#include <cstdio>
#include <fstream>
#include <iostream>
#include "src/main.hxx"

using namespace std;




// You can define datatype with -DTYPE=...
#ifndef TYPE
#define TYPE double
#endif
// You can define number of threads with -DMAX_THREADS=...
#ifndef MAX_THREADS
#define MAX_THREADS 12
#endif




template <class G, class H>
void runPagerank(const G& x, const H& xt, int repeat) {
  using V = TYPE;
  vector<V> *init = nullptr;

  for (int sleepDurationMs=1; sleepDurationMs<=1000; sleepDurationMs*=10) {
    for (float sleepProbability=0.0f; sleepProbability<1.01f; sleepProbability+=0.2f) {
      auto fv = [&](ThreadInfo *thread, auto v) {};
      // Find sequential pagerank for reference (synchronous, no dead ends).
      auto a0 = pagerankBasicSeq(xt, init, {1});
      // Find multi-threaded pagerank (synchronous, no dead ends).
      auto a1 = pagerankBasicOmp(xt, init, {repeat}, fv);
      auto e1 = l1NormOmp(a1.ranks, a0.ranks);
      printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankBasicOmp       {sleep_prob: %.1f, sleep_dur: %04d ms}\n", a1.time, a1.iterations, e1, sleepProbability, sleepDurationMs);
      // Find multi-threaded barrier-free pagerank (asynchronous, no dead ends).
      auto a2 = pagerankBarrierfreeOmp<true>(xt, init, {repeat}, fv);
      auto e2 = l1NormOmp(a2.ranks, a0.ranks);
      printf("[%09.3f ms; %03d iters.] [%.4e err.] pagerankBarrierfreeOmp {sleep_prob: %.1f, sleep_dur: %04d ms}\n", a2.time, a2.iterations, e2, sleepProbability, sleepDurationMs);
    }
  }
}


int main(int argc, char **argv) {
  char *file = argv[1];
  int repeat = argc>2? stoi(argv[2]) : 1;
  printf("Loading graph %s ...\n", file);
  OutDiGraph<int64_t> x;
  readMtxW(x, file); println(x);
  auto fl = [](auto u) { return true; };
  selfLoopU(x, None(), fl); print(x); printf(" (selfLoopAllVertices)\n");
  auto xt = transposeWithDegree(x); print(xt); printf(" (transposeWithDegree)\n");
  omp_set_num_threads(MAX_THREADS);
  printf("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  runPagerank(x, xt, repeat);
  printf("\n");
  return 0;
}
