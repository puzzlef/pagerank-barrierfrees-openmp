#pragma once
#include <utility>
#include <chrono>
#include "_debug.hxx"

using std::pair;
using std::chrono::microseconds;
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;




// PAIR
// ----

template <class K, class V>
struct PairFirst  { inline K& operator()(pair<K, V>& x) noexcept { return x.first; } };
template <class K, class V>
struct PairSecond { inline V& operator()(pair<K, V>& x) noexcept { return x.second; } };
template <class K, class V>
struct ConstPairFirst  { inline const K& operator()(const pair<K, V>& x) noexcept { return x.first; } };
template <class K, class V>
struct ConstPairSecond { inline const V& operator()(const pair<K, V>& x) noexcept { return x.second; } };
template <class K, class V>
struct PairFirstValue  { inline K operator()(const pair<K, V>& x) noexcept { return x.first; } };
template <class K, class V>
struct PairSecondValue { inline V operator()(const pair<K, V>& x) noexcept { return x.second; } };




// MEASURE-DURATION
// ----------------

inline auto timeNow() {
  return high_resolution_clock::now();
}
template <class T>
inline float durationMilliseconds(const T& start, const T& stop) {
  ASSERT(stop >= start);
  auto a = duration_cast<microseconds>(stop - start);
  return a.count()/1000.0f;
}
template <class T>
inline float durationMilliseconds(const T& start) {
  auto stop = timeNow();
  return durationMilliseconds(start, stop);
}


template <class F>
inline float measureDuration(F fn, int N=1) {
  ASSERT(N>0);
  auto start = high_resolution_clock::now();
  for (int i=0; i<N; i++)
    fn();
  auto stop = high_resolution_clock::now();
  return durationMilliseconds(start, stop)/N;
}


template <class F>
inline float measureDurationMarked(F fn, int N=1) {
  ASSERT(N>0);
  float duration = 0;
  for (int i=0; i<N; i++)
    fn([&](auto fm) { duration += measureDuration(fm); });
  return duration/N;
}




// RETRY
// -----

template <class F>
inline bool retry(F fn, int N=2) {
  ASSERT(N>0);
  for (int i=0; i<N; i++)
    if (fn()) return true;
  return false;
}




// MOVE
// ----
// Conditional move, otherwise value.

#define CMOVE(c, t, f) \
  ((c)? move(t) : (f))

#define CMOVE_VECTOR(t, f) \
  CMOVE(!(t).empty(), t, f)

#define CMOVE_GRAPH(t, f) \
  CMOVE((t).order()>0, t, f)
