# The physics of better and better N-body methods

A qualitative picture, with no equations, of what the three algorithms in this project are doing and how they differ. The central object throughout is the gravitational potential: the field that all the masses together create in space. A body moves because of the slope of that field where it sits.

For the details of each method, see `README.md` (the project overview) and `src/fmm/PHYSICS.md` (the Fast Multipole Method in depth).

## The problem

Every body adds its own dimple to the potential, and the dimples add. To move the bodies, you need the potential's slope at each body's position, which means the combined dimples of all the others. With N bodies that is N bodies each asking about N − 1 others, so the work grows with the square of N. Every faster method is a way of not asking about everyone.

## The idea they all rely on: distance washes out detail

A distant group of masses cannot be told apart from a single lump with the same total mass. Look closer and you start to see that the group is lopsided, then elongated, then more complicated. Each finer feature fades with distance faster than the one before. So the farther away you are, the shorter the list of features you need to reproduce the group's field well.

Everything below is a different way of using that fact.

## Direct summation: ask everyone

Each body asks every other body. This is exact (up to the time integrator), and costs the full square. It is the reference the other methods are checked against.

## Barnes-Hut: each body consults the summaries

The bodies are organized into nested groups. For each body, a group that is far enough away, relative to its own size, is replaced by a single lump at its center of mass, and a group that is too close is opened up into its parts.

- **What is summarized:** each group's outgoing description, at the lowest order.
- **Who uses it:** each body, separately. Bodies that are next to each other repeat nearly the same consultations.
- **The accuracy knob:** how far away a group must be, compared with its size, to count as a lump.
- **Cost:** each body consults a number of groups that grows slowly with N, since distant groups are large and few. That gives roughly N log N.

## Fast Multipole Method: regions share an ambient field

FMM keeps the group summaries, but adds the field's description from the receiving end. Over a small region, the field from far-away masses is smooth: a uniform pull plus a gentle tidal stretch, plus finer variations. A region builds that ambient field once, and every body in it just reads it at its own position.

- **What is summarized:** each group's outgoing description, and each region's ambient field, both with several terms.
- **Who uses it:** a region builds the field once, and all of its bodies read it.
- **The accuracy knob:** the distance criterion, and how many terms are kept.
- **Cost:** the work per body doesn't grow with N, so the total is roughly proportional to N.

## Side by side

| | Direct | Barnes-Hut | FMM |
|---|---|---|---|
| picture | every body asks every body | each body asks each far group | each region builds a field, bodies read it |
| what the far world is to a body | everyone, individually | a list of lumps to consult | one smooth field |
| reuse between neighbouring bodies | none | none | all of it |
| approximation | none | far groups reduced to a lump | far-world fields cut off at a finest detail |
| scaling | N² | N log N | N |

## An Earth and Moon picture

Imagine asking what the Moon does to every water molecule in the ocean:

- **Direct:** every molecule asks every rock in the Moon.
- **Barnes-Hut:** each molecule treats the Moon as one lump and works out its pull itself.
- **FMM:** you work out the Moon's field over the whole Earth once, as a uniform pull plus a tidal stretch, and every molecule reads off its local value.

## What they have in common

All the faster methods are approximate in the same sort of way: a distant group's detailed structure is thrown away. The approximation is controlled, and tightening the knob recovers the exact answer in the limit. None of them averages anything beyond that. Neither Barnes-Hut nor FMM conserves momentum exactly, because the pull a body feels from a group is not exactly the pull it exerts back.

## When each wins

Direct summation is simplest and is the right choice for small N. Barnes-Hut has little overhead and works well at moderate accuracy and moderate N. FMM's savings come from sharing work, so it has more bookkeeping per region. It pays off at large N, or when high accuracy is needed, since adding terms is cheaper than tightening the distance criterion. Exactly where the methods cross over depends on the implementation, and the benchmarks in `README.md` are the place to measure it.
