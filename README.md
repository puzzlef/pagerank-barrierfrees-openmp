Comparison of barrier vs barrier-free iterations with fixed iterations in
[OpenMP]-based [PageRank algorithm] for [link analysis].

**Unordered PageRank** is the *standard* approach of PageRank computation (as
described in the original paper by Larry Page et al. [(1)]), where *two*
*different rank vectors* are maintained; one representing the *current* ranks of
vertices, and the other representing the *previous* ranks. On the other hand,
**ordered PageRank** uses *a single rank vector*, representing the current ranks
of vertices [(2)]. This is similar to barrierless non-blocking implementations
of the PageRank algorithm by Hemalatha Eedi et al. [(3)]. As ranks are updated
in the same vector (with each iteration), the order in which vertices are
processed *affects* the final result (hence the adjective *ordered*). However,
as PageRank is an iteratively converging algorithm, results obtained with either
approach are *mostly the same*. **Barrier-free PageRank** is an *ordered*
*PageRank* where each thread processes a subset of vertices in the graph
independently, *without* waiting (with a barrier) for other threads to complete an
iteration. This minimizes pointless waits and allows each thread to be on a
*different iteration number* (which may or may not be beneficial for convergence)
[(3)].

In this experiment, we fix the number of iterations of PageRank computation to
`5`, and compare the performance of two different approaches of barrier-free
iterations of *OpenMP-based ordered PageRank* (along with unordered and ordered
OpenMP-based approaches). With one barrier-free PageRank approach each thread
detects convergence by measuring the difference between the previous and the
current ranks of all the vertices (**full**). With the other approach,
difference is measured between the previous and current ranks of only the subset
of vertices being processed by each thread (**part**). We use the follwing
PageRank parameters: damping factor `α = 0.85`, tolerance `τ = 10^-6`, and limit
the maximum number of iterations to `L = 5.` The error between the approaches is
calculated with *L1-norm*. The *sequential unordered* approach is considered to
be the *gold standard* (wrt to which error is measured). *Dead ends* in the
graph are handled by always teleporting any vertex in the graph at random
(*teleport* approach [(4)]). The teleport contribution to all vertices is
calculated *once* (for all vertices) at the begining of each iteration.

From the results, we observe that **both the barrier-free PageRank** approaches
are **generally slower than OpenMP-based approaches** to complete the same
number of iterations, **but are faster on some graphs**. This indicates that
avoiding a barrier is beneficial on this subset of graphs, but not on the
others.

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
# [00005.631 ms; 005 iters.] [9.3543e-02 err.] pagerankOmpUnordered
# [00007.649 ms; 005 iters.] [5.1495e-02 err.] pagerankOmpOrdered
# [00010.924 ms; 005 iters.] [5.2899e-02 err.] pagerankBarrierfreeFullOmp
# [00008.421 ms; 005 iters.] [5.3969e-02 err.] pagerankBarrierfreePartOmp
#
# Loading graph /home/subhajit/data/web-BerkStan.mtx ...
# order: 685230 size: 7600595 [directed] {}
# order: 685230 size: 8285825 [directed] {} (selfLoopAllVertices)
# order: 685230 size: 8285825 [directed] {} (transposeWithDegree)
# OMP_NUM_THREADS=12
# [00009.800 ms; 005 iters.] [9.9160e-02 err.] pagerankOmpUnordered
# [00010.718 ms; 005 iters.] [6.3629e-02 err.] pagerankOmpOrdered
# [00016.638 ms; 005 iters.] [6.2224e-02 err.] pagerankBarrierfreeFullOmp
# [00010.984 ms; 005 iters.] [6.5731e-02 err.] pagerankBarrierfreePartOmp
#
# ...
```

[![](https://i.imgur.com/zcuooEK.png)][sheetp]
[![](https://i.imgur.com/W1k2kCS.png)][sheetp]

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


[![](https://i.imgur.com/6UYGStl.jpg)](https://www.youtube.com/watch?v=2k2orAjC4aw)<br>


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
[gist]: https://gist.github.com/wolfram77/eb498865f4650662bd562a6803031971
[charts]: https://imgur.com/a/5K6SfVe
[sheets]: https://docs.google.com/spreadsheets/d/1N-r6YU70b5nPqSibkYSfSE9GKopGI-xFJ8aju3tG0wA/edit?usp=sharing
[sheetp]: https://docs.google.com/spreadsheets/d/e/2PACX-1vQH1Jqm79qSz7GeQLlE2RU1iIYPXIiAvqEqxD9MwEk4_lb6go-j5iZREF90LVlGEeNg0UZfDKqtBQ06/pubhtml
