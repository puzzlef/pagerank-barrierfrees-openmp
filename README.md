Design of *OpenMP-based* **statically scheduled Barrier-free** *PageRank
algorithm* for link analysis.

**Unordered PageRank** is the *standard* approach of PageRank computation (as
described in the original paper by Larry Page et al. [(1)][page]), where *two*
*different rank vectors* are maintained; one representing the *current* ranks of
vertices, and the other representing the *previous* ranks. On the other hand,
**ordered PageRank** uses *a single rank vector*, representing the current ranks
of vertices [(2)][pagerank]. This is similar to barrierless non-blocking implementations
of the PageRank algorithm by Hemalatha Eedi et al. [(3)][eedi]. As ranks are updated
in the same vector (with each iteration), the order in which vertices are
processed *affects* the final result (hence the adjective *ordered*). However,
as PageRank is an iteratively converging algorithm, results obtained with either
approach are *mostly the same*. **Barrier-free PageRank** is an *ordered*
*PageRank* where each thread processes a subset of vertices in the graph
independently, *without* waiting (with a barrier) for other threads to complete
an iteration. This minimizes unnecessary waits and allows each thread to be on a
*different iteration number* (which may or may not be beneficial for
convergence) [(3)][eedi]. **Monolithic PageRank** is the standard PageRank computation
where vertices are grouped by *strongly connected components* (SCCs). This
improves *locality* of memory accesses and thus improves performance [(4)][sahu].
**Levelwise PageRank** is a decomposed form of PageRank computation, where each
*SCC* is processed by *topological order* in the *block-graph* (all components
in a *level* are processed together *until* convergence, after which we proceed
to the next level). This decomposition allows for distributed computation
*without per-iteration communication*. However, it does not work on a graph
which includes *dead ends* (vertices with no outgoing edges, also called
dangling nodes) [(5)][levelwise].

<br>


### Performing Fixed iterations

In this experiment ([fixed-iterations]), we fix the number of iterations of
PageRank computation to `5`, and compare the performance of two different
approaches of barrier-free iterations of *OpenMP-based ordered PageRank* (along
with unordered and ordered OpenMP-based approaches). With one barrier-free
PageRank approach each thread detects convergence by measuring the difference
between the previous and the current ranks of all the vertices (**full**). With
the other approach, difference is measured between the previous and current
ranks of only the subset of vertices being processed by each thread (**part**).
We use the follwing PageRank parameters: damping factor `α = 0.85`, tolerance
`τ = 10^-6`, and limit the maximum number of iterations to `L = 5.` The error
between the approaches is calculated with *L1-norm*. The *sequential unordered*
approach is considered to be the *gold standard* (wrt to which error is
measured). *Dead ends* in the graph are handled by always teleporting any vertex
in the graph at random (*teleport* approach [(1)][teleport]). The teleport
contribution to all vertices is calculated *once* (for all vertices) at the
begining of each iteration.

From the results, we observe that **both the barrier-free PageRank** approaches
are **generally slower than OpenMP-based approaches** to complete the same
number of iterations, **but are faster on some graphs**. This indicates that
avoiding a barrier is beneficial on this subset of graphs, but not on the
others.

[fixed-iterations]: https://github.com/puzzlef/pagerank-barrierfrees-openmp/tree/fixed-iterations

<br>


### Using std::thread

In this experiment ([using-threads]), we perform two different approaches of
barrier-free iterations of *OpenMP-based ordered PageRank*; one in which each
thread detects convergence by measuring the difference between the previous and
the current ranks of all the vertices (**full**), and the other in which the
difference is measured between the previous and current ranks of only the
portion of vertices being processed by each thread (**part**). We compare them
with the same two approaches of barrier-free iterations of *std::thread based*
*ordered PageRank* instead. This experiment is done while adjusting the tolerance
`τ` from `10^-1` to `10^-14` with three different tolerance functions:
`L1-norm`, `L2-norm`, and `L∞-norm`. We use a damping factor of `α = 0.85` and
limit the maximum number of iterations to `L = 500`. The error between the
approaches is computed with *L1-norm*. The *sequential unordered* approach is
considered to be the *gold* *standard* (wrt to which error is measured). *Dead*
*ends* in the graph are handled by adding self-loops to all vertices, such that
there are no dead ends (*loopall* approach [(1)][teleport]).

From the results, we observe that the **std::thread based approaches** have
performance **almost equivalent to that of the OpenMP-based approaches**. This
is not surprising, since OpenMP was essential used just to create threads (using
a `parallel for`).

[using-threads]: https://github.com/puzzlef/pagerank-barrierfrees-openmp/tree/using-threads

<br>


### Adjusting Tolerance

In this experiment ([adjust-tolerance]), we perform two diffrent approaches of
barrier-free iterations of *OpenMP-based ordered PageRank*; one in which each
thread detects convergence by measuring the difference between the previous and
the current ranks of all the vertices (**full**), and the other in which the
difference is measured between the previous and current ranks of only the subset
of vertices being processed by each thread (**part**). This is done while
adjusting the tolerance `τ` from `10^-1` to `10^-14` with three different
tolerance functions: `L1-norm`, `L2-norm`, and `L∞-norm`. We also compare it
with OpenMP-based unordered and ordered PageRank for the same tolerance and
tolerance function. We use a damping factor of `α = 0.85` and limit the maximum
number of iterations to `L = 500`. The error between the approaches is
calculated with *L1-norm*. The *sequential unordered* approach is considered to
be the *gold standard* (wrt to which error is measured). *Dead ends* in the
graph are handled by always teleporting any vertex in the graph at random
(*teleport* approach [(1)][teleport]). The teleport contribution to all vertices is
calculated *once* (for all vertices) at the begining of each iteration.

From the results, we observe that the **partial-error based barrier-free**
**PageRank** is **faster than full-error version**, but generally has **similar**
**performance as that of standard OpenMP-based PageRank (unordered, ordered)**,
although it **requires the fewest iterations to converge**. However the,
**partial barrier-free approach** also **generally yields greater error in**
**ranks** (with respect to sequential unordered approach).

[adjust-tolerance]: https://github.com/puzzlef/pagerank-barrierfrees-openmp/tree/adjust-tolerance

<br>


### Adjusting Tolerance with Monolithic technique

In this experiment ([adjust-tolerance-monolithic]), we perform two different
approaches of barrier-free iterations of *OpenMP-based ordered PageRank*; one in
which each thread detects convergence by measuring the difference between the
previous and the current ranks of all the vertices (**full**), and the other in
which the difference is measured between the previous and current ranks of only
the subset of vertices being processed by each thread (**part**). Both approahes
are performed in either the standard way (with vertices arranged in vertex-id
order) or with the monolithic approach (where vertices are grouped by SCCs).
This is done while adjusting the tolerance `τ` from `10^-1` to `10^-14` with
three different tolerance functions: `L1-norm`, `L2-norm`, and `L∞-norm`. We
also compare it with OpenMP-based unordered and ordered PageRank for the same
tolerance and tolerance function. We use a damping factor of `α = 0.85` and
limit the maximum number of iterations to `L = 500`. The error between the
approaches is calculated with *L1-norm*. The *sequential unordered* approach is
considered to be the *gold standard* (wrt to which error is measured). *Dead*
*ends* in the graph are handled by self loops to all the vertices (*loopall*
approach [(1)][teleport]).

From the results, we observe that **monolithic approaches are faster** than
default approaches in all the cases. We also observe that **monolithic**
**barrier-free** **approach with partial error measurement** is the **fastest** in
terms of *time*, and in terms of *iterations* for tolerance below `10^-10`.
However, *barrier-free approach with full error measurement* appears to take the
most time.

[adjust-tolerance-monolithic]: https://github.com/puzzlef/pagerank-barrierfrees-openmp/tree/adjust-tolerance-monolithic

<br>


### Adjusting Tolerance with Levelwise technique

In this experiment ([adjust-tolerance-levelwise]), we perform two different
approaches of barrier-free iterations of *OpenMP-based ordered PageRank*; one in
which each thread detects convergence by measuring the difference between the
previous and the current ranks of all the vertices (**full**), and the other in
which the difference is measured between the previous and current ranks of only
the subset of vertices being processed by each thread (**part**). Both approahes
are performed in either the *default* way (with vertices arranged in vertex-id
order), with the *monolithic* approach (where vertices are grouped by SCCs), or
with the *levelwise* approach (where vertices are processed until convergence in
groups of levels in the block-graph). This is done while adjusting the tolerance
`τ` from `10^-1` to `10^-14` with `L∞-norm` as the tolerance function. We also
compare it with OpenMP-based unordered and ordered PageRank for the same
tolerance and tolerance function. We use a damping factor of `α = 0.85` and
limit the maximum number of iterations to `L = 500`. The error between the
approaches is calculated with *L1-norm*. The *sequential* *unordered* approach
is considered to be the *gold standard* (wrt to which error is measured). *Dead*
*ends* in the graph are handled by self loops to all the vertices (*loopall*
approach [(1)][teleport]).

From the results, we observe that **monolithic approaches are faster** for
OpenMP-based unordered/unordered approaches, and **levelwise approaches are**
**faster** for barrier-free approaches with full/partial error measurement. While
we do see that **levelwise barrier-free approach with partial error**
**measurement** is the **fastest** both in terms of *time* and in the number of
*iterations* the difference is generally small. However when performing PageRank
computation in a distributed setting, we expect levelwise based approaches to
provide a good advantage, thanks to its reduced communication requirement.

[adjust-tolerance-levelwise]: https://github.com/puzzlef/pagerank-barrierfrees-openmp/tree/adjust-tolerance-levelwise

<br>


### Other experiments

- [compare-withbarrier](https://github.com/puzzlef/pagerank-barrierfrees-openmp/tree/compare-withbarrier)
- [with-random-sleep-uniform](https://github.com/puzzlef/pagerank-barrierfrees-openmp/tree/with-random-sleep-uniform)
- [with-random-sleep-nonuniform](https://github.com/puzzlef/pagerank-barrierfrees-openmp/tree/with-random-sleep-nonuniform)

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


[![](https://i.imgur.com/50yaKL7.jpg)](https://www.youtube.com/watch?v=g2tMcMQqSbA)<br>
[![ORG](https://img.shields.io/badge/org-puzzlef-green?logo=Org)](https://puzzlef.github.io)
[![DOI](https://zenodo.org/badge/530082702.svg)](https://zenodo.org/badge/latestdoi/530082702)


[Prof. Dip Sankar Banerjee]: https://sites.google.com/site/dipsankarban/
[Prof. Kishore Kothapalli]: https://faculty.iiit.ac.in/~kkishore/
[Prof. Sathya Peri]: https://people.iith.ac.in/sathya_p/
[SuiteSparse Matrix Collection]: https://sparse.tamu.edu
[page]: https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.38.5427
[eedi]: https://ieeexplore.ieee.org/document/9407114
[sahu]: https://ieeexplore.ieee.org/document/9835216
[levelwise]: https://gist.github.com/wolfram77/12e5a19ff081b2e3280d04331a9976ca
[openmp]: https://github.com/puzzlef/pagerank-openmp
[pagerank]: https://github.com/puzzlef/pagerank-ordered-vs-unordered
[teleport]: https://gist.github.com/wolfram77/94c38b9cfbf0c855e5f42fa24a8602fc
[OpenMP]: https://en.wikipedia.org/wiki/OpenMP
[PageRank algorithm]: https://en.wikipedia.org/wiki/PageRank
[link analysis]: https://en.wikipedia.org/wiki/Network_theory#Link_analysis
