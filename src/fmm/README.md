# Fast Multipole Method (FMM)

Direct summation costs O(N²), Barnes-Hut O(N log N), FMM O(N).
This is an overview of the idea, using a school as an analogy.

## The analogy

You want to know how loud the school is at your desk.

- **Slow way (direct):** ask every kid in the school. Everyone does this.
- **Barnes-Hut:** each classroom's teacher writes one note summarizing the class (the cell's total mass and center of mass). A kid reads the notes of far classrooms instead of asking each kid, and only asks individuals nearby. Every kid still does their own rounds, so the same notes get read over and over.
- **FMM:** the *classroom* reads the far notes once and writes a single "what it sounds like in this room" note. Every kid just reads that one note. The repeated work is done once per room, not once per kid.

| school | simulation |
|---|---|
| classroom / grade / school | tree cell at a given level |
| teacher's note about the kids | the cell's *summary* (multipole expansion) |
| "what it sounds like in this room" | the cell's *local field* (local expansion) |
| kid | body |
| sitting right next to you | near field, summed directly |

## Why a short note is enough

- **Summary (multipole):** seen from far away, a cell's mass distribution can be described by a few numbers (total mass, then dipole, quadrupole, ...). Barnes-Hut keeps only the first one, the mass at the center of mass.
- **Local field:** seen from inside a cell, the pull of far-away mass varies smoothly. A few numbers (value, slope, curvature, ...) describe it, and every body in the cell reads its acceleration from the same few numbers.
- Both are truncated at order `p`. This is the accuracy knob.

The tool for both is the expansion of the kernel

    1/|r − r'| = Σₙ (r'ⁿ / rⁿ⁺¹) Pₙ(cos γ),   |r'| < |r|

which separates the source position `r'` from the observation point `r`. The addition theorem rewrites `Pₙ` as products of spherical harmonics.

## The steps

Standard names are in brackets.

| step | what it does | school |
|---|---|---|
| summarize a leaf (P2M) | bodies in a leaf → its summary | a teacher asks each kid and writes one note |
| summarize upward (M2M) | children's summaries → parent's summary | a grade head combines class notes; nobody re-asks the kids |
| far note → local field (M2L) | a well-separated cell's summary → local field at a target cell | your teacher works out what that class sounds like *here* and adds it to the room's note |
| push field down (L2L) | parent's local field → children's local fields | the principal's note is handed down to grades, then to rooms |
| read off acceleration (L2P) | leaf's local field → acceleration at each body | a kid works out what the room's note means at their own desk |
| nearby directly (P2P) | body-to-body sum with close leaves | a kid listens to the neighbours themselves |

The only approximation is truncating the series at order `p`. Everything else is exact bookkeeping.

## Order of work

1. **Build the tree.**
2. **Upward pass:** summarize each leaf (P2M), then summarize upward to the root (M2M). All summaries must be finished before the next stage.
3. **Traversal:** walk pairs of cells. If two cells are well separated, turn the source's summary into a local field at the target (M2L). If they are close and both are leaves, sum bodies directly (P2P). Otherwise split the larger cell and recurse.
4. **Downward pass:** push local fields down to the leaves (L2L), then read off an acceleration for each body (L2P).
5. **Apply** the accelerations (as in `direct` and Barnes-Hut).

## What a body's acceleration is made of

For a body in leaf `B`:

    a(x) = −∇φ_near(x) − ∇φ_far(x)

- `φ_near`: bodies in `B` and its touching neighbours, summed directly.
- `φ_far`: evaluated from `B`'s local field. The local field of a cell `A` is

      L(A) = Σ_{C ∈ W(A)} M2L(C) + L2L(L(parent of A))

  where `W(A)` is the set of cells well separated from `A` whose parents were *not* well separated from `A`'s parent. This avoids counting any source twice: everything far from the parent has already come down from above.

## Why it's O(N)

Each cell interacts with a bounded number of partners, and the cost of each interaction depends only on `p`, not on `N`. The number of cells is proportional to `N`. Each body's own work (summarize leaf, read off acceleration, nearby bodies) is a constant.

## Design decisions

- Spherical-harmonic expansions.
- A separate `FmmTree`, reusing the ideas of `System` and `Tree` from Barnes-Hut. The Barnes-Hut code is unchanged.
- Dual tree traversal (no precomputed interaction lists).
- Leaves may hold several bodies, unlike Barnes-Hut.
- Not yet decided: the order `p`, the leaf size, the exact separation criterion.

## Checking correctness

Compare against `direct`. The error should shrink as `p` grows.
