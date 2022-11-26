#pragma once
#include <utility>
#include <random>
#include <atomic>
#include <vector>
#include "_main.hxx"

using std::random_device;
using std::default_random_engine;
using std::atomic;
using std::vector;
using std::move;




// PAGERANK-OPTIONS
// ----------------

enum NormFunction {
  L0_NORM = 0,
  L1_NORM = 1,
  L2_NORM = 2,
  LI_NORM = 3
};


template <class V>
struct PagerankOptions {
  int repeat;
  int toleranceNorm;
  V   tolerance;
  V   damping;
  int maxIterations;

  PagerankOptions(int repeat=1, int toleranceNorm=LI_NORM, V tolerance=1e-10, V damping=0.85, int maxIterations=500) :
  repeat(repeat), toleranceNorm(toleranceNorm), tolerance(tolerance), damping(damping), maxIterations(maxIterations) {}
};




// PAGERANK-RESULT
// ---------------

template <class V>
struct PagerankResult {
  vector<V> ranks;
  int   iterations;
  float time;

  PagerankResult() :
  ranks(), iterations(0), time(0) {}

  PagerankResult(vector<V>&& ranks, int iterations=0, float time=0) :
  ranks(ranks), iterations(iterations), time(time) {}

  PagerankResult(vector<V>& ranks, int iterations=0, float time=0) :
  ranks(move(ranks)), iterations(iterations), time(time) {}
};




// THREAD-INFO
// -----------

struct ThreadInfo {
  random_device dev;          // used for random sleeps
  default_random_engine rnd;  // used for random sleeps
  int id;                     // thread number
  int iteration;              // current iteration

  ThreadInfo(int id) :
  dev(), rnd(dev()), id(id), iteration(0) {}
  inline void clear() { iteration = 0; }
};
