# A+ (Plus Group) Implementation: Two Approaches Compared

## Overview

Both the `aplus` and `aplus3` branches implement the same fisheries "plus group" (A+)
accumulator — the oldest age class where survivors and newly-graduated individuals pool
together rather than continuing to age into new bins — but they take structurally different
approaches to scheduling and synchronizing A+ within the existing MPI cohort-parallel driver
(`main_cohort.cpp`). Both have been verified numerically equivalent: with A+ enabled they
produce identical checksums to each other, and with `-no-aplus` both exactly reproduce the
pre-A+ baseline, satisfying the conservation check `no-aplus checksum == (A+-on normal
checksum) + (A+-on A+ checksum)` to the last printed digit, on the real `skj_Fat.xml` config
(99/100 age groups, 480 months, 10 ranks).

## `aplus`: graph-based scheduling

A+ is folded directly into the same `SeapodymCohortDependencyAnalyzer` task graph that already
governs every ordinary cohort. Each calendar step of A+ becomes its own one-step task, appended
after the normal task IDs, with two explicit dependency edges: on A+'s own previous step, and on
the specific cohort step that "graduates" into it that month. Newborn cohorts get a third kind
of edge added automatically — a dependency on the matching A+ step, so a spawning cohort can
never be dispatched before A+'s density for that month is available. All of this is dispatched
by the same `TaskStepManager`/`TaskStepWorker` pair used for everything else; A+ tasks are just
another branch inside the existing `taskFunction`. No dedicated MPI rank, no extra communicator
split, no custom message protocol — correctness follows directly from the dependency graph's
existing "don't dispatch a task until every dependency is in the `completed` set" logic, a
mechanism that was already in use and already exercised long before A+ existed.

## `aplus3`: dedicated worker with ping-pong handoff

A+ is removed from the dependency graph entirely and instead runs as its own independent loop
on a dedicated MPI rank (`runAPlusWorker`). Cohorts that reach the oldest normal age class
("feeders") hand their graduating density directly to that rank via a blocking send, embedding
the payload in the message itself; A+ acknowledges once it has actually processed that calendar
step, and only after that ACK does the feeder tell the manager it's done — which is what
prevents anything downstream from running ahead of A+. Because feeders can complete out of
order, A+ carries a small reorder buffer that holds early arrivals and drains them strictly by
calendar step. This is a real custom protocol with its own state machine, and finding it correct
took real debugging effort — we found and fixed a genuine out-of-order feeding bug here that the
graph-based design was never exposed to, because it doesn't have a hand-rolled synchronization
protocol to get wrong in the first place.

## Pros and cons

| | `aplus` (graph) | `aplus3` (dedicated worker) |
|---|---|---|
| Minimum ranks | 2 | 3 |
| Extra machinery | None — reuses existing graph/dispatch | Dedicated worker loop, message tags, reorder buffer, ACK gating |
| Ordering guarantee | Declarative (graph edges) | Operational (blocking ACK + buffer) |
| A+ compute isolation | Shares dispatch queue with all cohorts | Own rank, never waits on farm backlog |
| Rank utilization | Every rank does real cohort work | A+'s rank is mostly idle (480 single-step advances vs. 50+ cohorts/rank elsewhere) |
| Debugging surface | Small — one well-tested mechanism | Larger — a protocol we had to design and validate ourselves |
| Verified correctness | Yes | Yes |

`aplus`'s case is really about simplicity and efficiency: it adds no new synchronization
primitive, needs no extra rank, and every rank stays busy with real cohort work. Its downside is
that A+'s own advancement rides in the same queue as everything else, so if the farm is ever
fully saturated, A+ could in principle wait behind unrelated backlog — though the graph
guarantees it's never *wrong*, just potentially not the very next thing scheduled.

`aplus3`'s case is about isolation: A+ gets guaranteed, uncontended throughput on its own rank,
decoupled from farm scheduling. The cost is a dedicated rank that, based on the timing captured
during testing, is doing meaningfully less work than any other rank (a single one-step-per-month
accumulator vs. dozens of full multi-step cohort trajectories), plus a bespoke protocol that's
more code to maintain and was the source of the one genuinely subtle bug found in this project
(the feeder-ordering gap).

## Recommendation

Converge on `aplus`. For the workload as it exists today — a single accumulator doing one step
of ordinary adult dynamics per calendar month — there's no evidence A+ is compute-bound enough
to need its own rank, and the graph-based design gets the same correctness guarantee with less
code, fewer ranks, and a synchronization mechanism that's already been proven out across the
rest of the scheduler rather than a new one built and debugged specifically for this feature.

`aplus3` isn't wrong, and it's a legitimate design if A+'s per-step cost ever grows to actually
compete for farm time — a much heavier per-step computation, or multiple accumulator bins each
needing real throughput — but that's not the situation today, and paying for a dedicated rank
plus a bespoke protocol ahead of that need isn't buying anything. Worth keeping `aplus3` around
as a documented alternative rather than deleting it, since the two are now proven equivalent and
it remains a reasonable fallback if requirements change, but it shouldn't be the branch new work
builds on by default.
