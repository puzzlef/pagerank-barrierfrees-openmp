Effect of using different values of tolerance with barrier-free iterations in
[OpenMP]-based ordered [PageRank algorithm] for [link analysis].

**Unordered PageRank** is the *standard* approach of PageRank computation (as
described in the original paper by Larry Page et al. [(1)]), where *two*
*different rank vectors* are maintained; one representing the *current* ranks of
vertices, and the other representing the *previous* ranks. On the other hand,
**ordered PageRank** uses *a single rank vector*, representing the current ranks
of vertices [(2)]. This is similar to barrierless non-blocking implementations
of the PageRank algorithm by Hemalatha Eedi et al. [(3)]. As ranks are updated
in the same vector (with each iteration), the order in which vertices are
computed *affects* the final result (hence the modifier *ordered*). However,
as PageRank is an iteratively converging algorithm, results obtained with either
approach are *mostly the same*. **Barrier-free PageRank** is an *ordered*
*PageRank* where each thread processes a portion of vertices in the graph
independently, *without* waiting (with a barrier) for other threads to complete an
iteration. This minimizes unnecessary waits and allows each thread to be on a
*different iteration number* (which may or may not be beneficial for convergence)
[(3)].

In this experiment, we perform two different approaches of barrier-free
iterations of *OpenMP-based ordered PageRank*; one in which each thread detects
convergence by measuring the difference between the previous and the current
ranks of all the vertices (**full**), and the other in which the difference is
measured between the previous and current ranks of only the portion of vertices
being processed by each thread (**part**). We compare them with the same two
approaches of barrier-free iterations of *std::thread based ordered PageRank*
instead. This experiment is done while adjusting the tolerance `τ` from `10^-1`
to `10^-14` with three different tolerance functions: `L1-norm`, `L2-norm`, and
`L∞-norm`. We use a damping factor of `α = 0.85` and limit the maximum number of
iterations to `L = 500`. The error between the approaches is computed with
*L1-norm*. The *sequential unordered* approach is considered to be the *gold*
*standard* (wrt to which error is measured). *Dead ends* in the graph are
handled by adding self-loops to all vertices, such that there are no dead ends
(*loopall* approach [(4)]).

From the results, we observe that the **std::thread based approaches** have
performance **almost equivalent to that of the OpenMP-based approaches**. This
is not surprising, since OpenMP was essential used just to create threads (using
a `parallel for`).

All outputs are saved in a [gist] and a small part of the output is listed here.
Some [charts] are also included below, generated from [sheets]. The input data
used for this experiment is available from the [SuiteSparse Matrix Collection].
This experiment was done with guidance from [Prof. Kishore Kothapalli],
[Prof. Dip Sankar Banerjee], and [Prof. Sathya Peri].

<br>

```bash
$ g++ -std=c++17 -O3 -fopenmp main.cxx
$ ./a.out ~/data/web-Stanford.mtx
$ ./a.out ~/data/web-BerkStan.mtx
$ ...

# Loading graph /home/subhajit/data/web-Stanford.mtx ...
# order: 281903 size: 2312497 [directed] {}
# order: 281903 size: 2594400 [directed] {} (selfLoopAllVertices)
# order: 281903 size: 2594400 [directed] {} (transposeWithDegree)
# OMP_NUM_THREADS=12
# [00008.257 ms; 004 iters.] [7.0286e-02 err.] pagerankBarrierfreeFullOmp {tol_norm: L1, tolerance: 1e-01}
# [00003.034 ms; 001 iters.] [3.6772e-01 err.] pagerankBarrierfreePartOmp {tol_norm: L1, tolerance: 1e-01}
# [00008.753 ms; 004 iters.] [6.4310e-02 err.] pagerankBarrierfreeFullThd {tol_norm: L1, tolerance: 1e-01}
# [00003.192 ms; 001 iters.] [3.6264e-01 err.] pagerankBarrierfreePartThd {tol_norm: L1, tolerance: 1e-01}
# [00014.171 ms; 008 iters.] [2.1574e-02 err.] pagerankBarrierfreeFullOmp {tol_norm: L1, tolerance: 1e-02}
# [00005.684 ms; 003 iters.] [1.0533e-01 err.] pagerankBarrierfreePartOmp {tol_norm: L1, tolerance: 1e-02}
# [00015.681 ms; 008 iters.] [2.1268e-02 err.] pagerankBarrierfreeFullThd {tol_norm: L1, tolerance: 1e-02}
# [00006.014 ms; 003 iters.] [1.0038e-01 err.] pagerankBarrierfreePartThd {tol_norm: L1, tolerance: 1e-02}
# ...
# 00160.846 ms; 101 iters.] [1.9518e-07 err.] pagerankBarrierfreeFullOmp {tol_norm: Li, tolerance: 1e-14}
# [00114.781 ms; 094 iters.] [1.9623e-07 err.] pagerankBarrierfreePartOmp {tol_norm: Li, tolerance: 1e-14}
# [00155.972 ms; 101 iters.] [1.9418e-07 err.] pagerankBarrierfreeFullThd {tol_norm: Li, tolerance: 1e-14}
# [00119.789 ms; 094 iters.] [2.0086e-07 err.] pagerankBarrierfreePartThd {tol_norm: Li, tolerance: 1e-14}
#
# Loading graph /home/subhajit/data/web-BerkStan.mtx ...
# order: 685230 size: 7600595 [directed] {}
# order: 685230 size: 8285825 [directed] {} (selfLoopAllVertices)
# order: 685230 size: 8285825 [directed] {} (transposeWithDegree)
# OMP_NUM_THREADS=12
# [00012.362 ms; 004 iters.] [1.0036e-01 err.] pagerankBarrierfreeFullOmp {tol_norm: L1, tolerance: 1e-01}
# [00004.373 ms; 001 iters.] [3.3834e-01 err.] pagerankBarrierfreePartOmp {tol_norm: L1, tolerance: 1e-01}
# [00013.825 ms; 004 iters.] [1.0085e-01 err.] pagerankBarrierfreeFullThd {tol_norm: L1, tolerance: 1e-01}
# [00005.493 ms; 001 iters.] [3.4125e-01 err.] pagerankBarrierfreePartThd {tol_norm: L1, tolerance: 1e-01}
# [00022.049 ms; 009 iters.] [2.8209e-02 err.] pagerankBarrierfreeFullOmp {tol_norm: L1, tolerance: 1e-02}
# [00008.446 ms; 003 iters.] [1.0090e-01 err.] pagerankBarrierfreePartOmp {tol_norm: L1, tolerance: 1e-02}
# [00025.292 ms; 008 iters.] [2.7755e-02 err.] pagerankBarrierfreeFullThd {tol_norm: L1, tolerance: 1e-02}
# [00009.497 ms; 003 iters.] [1.0323e-01 err.] pagerankBarrierfreePartThd {tol_norm: L1, tolerance: 1e-02}
# ...
```

[![](https://i.imgur.com/kRQ3LmX.png)][sheetp]
[![](https://i.imgur.com/qnt5lsv.png)][sheetp]
[![](https://i.imgur.com/D6iZU0z.png)][sheetp]
[![](https://i.imgur.com/r9ktAoD.png)][sheetp]

[![](https://i.imgur.com/QFRJ9ix.png)][sheetp]
[![](https://i.imgur.com/lL2tQt6.png)][sheetp]
[![](https://i.imgur.com/0dnRETs.png)][sheetp]
[![](https://i.imgur.com/du7wOxH.png)][sheetp]

[![](https://i.imgur.com/XrIc7in.png)][sheetp]
[![](https://i.imgur.com/hfDNWjv.png)][sheetp]
[![](https://i.imgur.com/TTO2e7r.png)][sheetp]
[![](https://i.imgur.com/46KTsX5.png)][sheetp]

<br>
<br>


## References

- [An Efficient Practical Non-Blocking PageRank Algorithm for Large Scale Graphs; Hemalatha Eedi et al. (2021)](https://ieeexplore.ieee.org/document/9407114)
- [PageRank Algorithm, Mining massive Datasets (CS246), Stanford University](https://www.youtube.com/watch?v=ke9g8hB0MEo)
- [The PageRank Citation Ranking: Bringing Order to the Web; Larry Page et al. (1998)](https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427)
- [The University of Florida Sparse Matrix Collection; Timothy A. Davis et al. (2011)](https://doi.org/10.1145/2049662.2049663)
- [What's the difference between "static" and "dynamic" schedule in OpenMP?](https://stackoverflow.com/a/10852852/1413259)
- [OpenMP Dynamic vs Guided Scheduling](https://stackoverflow.com/a/43047074/1413259)

<br>
<br>


[![](https://i.imgur.com/cdNyVYm.jpg)](https://www.deviantart.com/jocelynminions/art/Sonic-Vs-Knuckles-Sonic-Movie2-906726255)<br>
[![DOI](https://zenodo.org/badge/532610283.svg)](https://zenodo.org/badge/latestdoi/532610283)


[(1)]: https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427
[(2)]: https://github.com/puzzlef/pagerank-ordered-vs-unordered
[(3)]: https://ieeexplore.ieee.org/document/9407114
[(4)]: https://gist.github.com/wolfram77/94c38b9cfbf0c855e5f42fa24a8602fc
[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[Prof. Sathya Peri]: https://people.iith.ac.in/sathya_p/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[OpenMP]: https://en.wikipedia.org/wiki/OpenMP
[PageRank algorithm]: https://en.wikipedia.org/wiki/PageRank
[link analysis]: https://en.wikipedia.org/wiki/Network_theory#Link_analysis
[gist]: https://gist.github.com/wolfram77/e9b2d0499fc67ecfb4f903f378e8734b
[charts]: https://imgur.com/a/hJ40WCV
[sheets]: https://docs.google.com/spreadsheets/d/13r2DIYgpE42x7fIU7KpQI7oNCoRg_wtrNZpOqw3CwQE/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vTVSeIDkmjmoU0WngTg1T9oYuT0uerdNGNf9eGyH5M1AiNejEd35TEIv_Dl8gEjM8MUliEoF0r53r-_/pubhtml
