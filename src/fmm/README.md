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
| teacher's note about the kids | the cell's *mass distribution* (multipole moments) |
| "what it sounds like in this room" | the cell's *external potential* (local expansion) |
| kid | body |
| sitting right next to you | near field, summed directly |

## Why a short note is enough

- **Mass distribution (multipole moments):** seen from far away, a cell's mass can be described by a few numbers (total mass, then dipole, quadrupole, ...). They describe what is *here*, packaged so that cells *there* can use it. Barnes-Hut keeps only the first number, the mass at the center of mass.
- **External potential (local expansion):** the potential felt inside a cell that is caused by everything outside its *near zone*. The near zone is the cell and its touching neighbours, whose bodies are summed directly. Seen from inside the cell, the potential from farther away varies smoothly, so a few numbers (value, slope, curvature, ...) describe it. It describes what is *there*, packaged for use *here*, and every body in the cell reads its acceleration from the same few numbers.
- Both are truncated at order `p`. This is the accuracy knob.

The tool for both is the expansion of the kernel

```text
1/|r − r'| = Σₙ (r'ⁿ / rⁿ⁺¹) Pₙ(cos γ),   |r'| < |r|
```

which separates the source position `r'` from the observation point `r`. The addition theorem rewrites `Pₙ` as products of spherical harmonics.

## Names in the code

Both are stored as a `SolidHarmonic`, a list of coefficients indexed by `(n, m)`.

| `FmmTree` field | standard name | stored as | valid |
|---|---|---|---|
| `mass_distribution` | multipole expansion | `SolidHarmonic::Irregular` | outside the cell |
| `external_potential` | local expansion | `SolidHarmonic::Regular` | inside the cell |

The potential at a body is split exactly, like an internal and an external potential of a subsystem:

    potential = external potential (series, truncated at p) + nearby bodies (summed directly)

The only approximation is the truncation of the external part. Nothing is averaged or neglected.

## The steps

Standard names are in brackets.

| step | what it does | school |
|---|---|---|
| summarize a leaf (P2M) | bodies in a leaf → its mass distribution | a teacher asks each kid and writes one note |
| summarize upward (M2M) | children's mass distributions → parent's | a grade head combines class notes; nobody re-asks the kids |
| far note → external potential (M2L) | a well-separated cell's mass distribution → external potential at a target cell | your teacher works out what that class sounds like *here* and adds it to the room's note |
| push potential down (L2L) | parent's external potential → children's external potentials | the principal's note is handed down to grades, then to rooms |
| read off acceleration (L2P) | leaf's external potential → acceleration at each body | a kid works out what the room's note means at their own desk |
| nearby directly (P2P) | body-to-body sum with close leaves | a kid listens to the neighbours themselves |

The only approximation is truncating the series at order `p`. Everything else is exact bookkeeping.

## Order of work

1. **Build the tree.**
2. **Accumulate the mass distributions (the upward pass):** summarize each leaf (P2M), then summarize upward to the root (M2M). All mass distributions must be finished before the next stage.
3. **Traversal:** walk pairs of cells. If two cells are well separated, turn the source's mass distribution into external potential at the target (M2L). If they are close and both are leaves, sum bodies directly (P2P). Otherwise split the larger cell and recurse.
4. **Downward pass:** push external potentials down to the leaves (L2L), then read off an acceleration for each body (L2P).
5. **Apply** the accelerations (as in `direct` and Barnes-Hut).

## What a body's acceleration is made of

For a body in leaf `B`:

```text
a(x) = −∇φ_near(x) − ∇φ_external(x)
```

- `φ_near`: bodies in `B` and its touching neighbours, summed directly.
- `φ_external`: evaluated from `B`'s external potential. The external potential of a cell `A` is

  ```text
  Φ(A) = Σ_{C ∈ W(A)} M2L(mass distribution of C) + L2L(Φ(parent of A))
  ```

  where `W(A)` is the set of cells well separated from `A` whose parents were *not* well separated from `A`'s parent. This avoids counting any source twice: everything far from the parent has already come down from above. The near zone grows with the size of the cell, so each level adds only the contributions that *became* external at that level.

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
