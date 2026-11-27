# better-n-better-n-body

A showcase of better and better algorithms for simulating gravitational N-body systems, written in Rust. Each algorithm trades a little accuracy for a lot of speed, and the code is meant to be read in order, from the simplest to the most sophisticated.

## Algorithms

| Algorithm      | Cost per step | Status  |
| -------------- | ------------- | ------- |
| Direct         | O(N²)         | Done    |
| Barnes-Hut     | O(N log N)    | Done    |
| Fast Multipole | O(N)          | Done    |

### Direct
Sums the pull of every body on every other body. It is exact (up to the integrator) and serves as the reference the other methods are checked against. See `src/direct.rs`.

### Barnes-Hut
Builds an octree over the bodies each step. A distant cell is treated as a single point mass at its center of mass, so the force on a body needs only about log N interactions instead of N. See `src/barnes_hut.rs` and `src/barnes_hut/tree.rs`.

The accuracy knob is the opening angle θ (called `resolution` in the code). A cell is treated as one point mass when `cell width < θ × distance to its center of mass`. At θ = 0 no cell is ever approximated, so the result matches direct summation. Larger θ is faster and less accurate.

### Fast Multipole Method
Barnes-Hut still makes every body consult the far cells one by one, and bodies that sit next to each other repeat nearly the same consultations. The FMM does that work once per cell. Each cell builds the smooth field that all the distant masses create over it, and every body in the cell reads its acceleration from that field. Near neighbours are still summed directly. The work per body stops growing with N, so the cost is O(N).

Each time step has four stages:

1. Build the tree (a cell is split while it holds more than 64 bodies).
2. Up: describe every cell's mass from afar (P2M at the leaves, M2M up the tree).
3. Across: one traversal of the tree. Cells that are far enough apart convert their description into a field over the other cell (M2L), and near leaves sum their bodies directly (P2P).
4. Down: each cell passes its field to its children (L2L), and each body reads its acceleration from its leaf's field (L2P).

The accuracy knob is the expansion order `p` (the `order` argument of `fmm::simulate`). A higher order is slower and more accurate. Two constants in `src/fmm/tree.rs` also matter: `MAX_LEAF_BODIES` (64) and `THETA` (0.5, how far apart two cells must be, compared with their size, for the far-field conversion to be used).

See `src/fmm.rs` and `src/fmm/`. For the idea in everyday terms see `src/fmm/README.md`, and for the physics in words, with no equations, see `src/fmm/PHYSICS.md` and `PHYSICS.md`.

## Running

```sh
cargo run --release   # times direct, Barnes-Hut and the FMM
cargo test            # unit tests, FMM and Barnes-Hut accuracy against direct summation
```

`main` simulates 100,000 random bodies for 10 steps (`dt` = 0.1) and prints the wall-clock time for direct summation, for Barnes-Hut at θ = 0.4, 0.6, 0.8 and 1.0, and for the FMM at order 2. The direct run takes several minutes at this size. Use `--release`: debug builds are far too slow to say anything meaningful about performance. The parameters (body count, time step, total time, order) are constants in `main` in `src/main.rs`.

## Benchmarks

What `main` prints: 100,000 uniform random bodies, 10 steps (`dt` = 0.1, total time 1.0), release build, a single run on one machine.

| Method              | Time  | Speedup vs. direct | Error    |
| ------------------- | ----- | ------------------ | -------- |
| Direct              | 456 s | 1.0×               | 0        |
| Barnes-Hut, θ = 0.4 | 61 s  | 7.5×               | 1.9×10⁻³ |
| Barnes-Hut, θ = 0.6 | 25 s  | 18×                | 5.8×10⁻³ |
| Barnes-Hut, θ = 0.8 | 13 s  | 35×                | 1.2×10⁻² |
| Barnes-Hut, θ = 1.0 | 8 s   | 57×                | 2.3×10⁻² |
| FMM, order 2        | 39 s  | 12×                | 6.0×10⁻³ |

`main` does not report the error. It is the relative root-mean-square difference from direct summation for one force evaluation, measured separately on a 20,000-body system, because direct summation at 100,000 bodies is slow. These are single runs, and the bodies move during the run, so the tree changes from step to step. Expect repeated runs to differ by about 10% or more.

Time and accuracy together, for **one force evaluation** of 100,000 uniform random bodies (the first step, with no time stepping), release build, one machine. Repeated single evaluations varied by about ±10%, and the per-step times in the run above are up to about 40% higher than these.

| Method              | Time    | Error   |
| ------------------- | ------- | ------- |
| Direct              | 41.6 s  | 0       |
| Barnes-Hut, θ = 0.3 | 10.6 s  | 9.7×10⁻⁴ |
| Barnes-Hut, θ = 0.4 | 5.3 s   | 1.9×10⁻³ |
| Barnes-Hut, θ = 0.5 | 3.0 s   | 3.7×10⁻³ |
| Barnes-Hut, θ = 0.6 | 1.8 s   | 5.8×10⁻³ |
| FMM, order 1        | 1.8 s   | 4.5×10⁻² |
| FMM, order 2        | 3.4 s   | 6.0×10⁻³ |
| FMM, order 3        | 6.8 s   | 1.1×10⁻³ |
| FMM, order 4        | 14.2 s  | 1.9×10⁻⁴ |

Compare methods at equal error, not at equal settings. At moderate accuracy (about 5×10⁻³) Barnes-Hut is faster: θ = 0.6 in 1.8 s against 3.4 s for order 2 here, and 25 s against 39 s in the full run above. The FMM at order 2 takes about as long as Barnes-Hut at θ = 0.5 but is less accurate. At finer accuracy the FMM is ahead (order 3 in 6.8 s against Barnes-Hut at θ = 0.3 in 10.6 s, for about the same error), but `main` only runs order 2, so the full run does not show it. The FMM has not been tuned yet; see the roadmap.

## Roadmap

Planned work on the Fast Multipole Method, in the order it is intended to be done. Each item is tracked as a GitHub issue. At order 2 the FMM's time splits roughly as 38% for the M2L conversions, 26% for the near field (P2P), and 22% for building the harmonics tables, so the first items go after the tables.

| # | Item | Why | Notes |
| - | ---- | --- | ----- |
| [#1](https://github.com/NonDairyNeutrino/better_n_better_n_body/issues/1) | Fuse and scale the M2L table without `powi` or `isqrt` | about 13% of the time is spent in `powi` and integer square roots when building each table | estimated 12–14% faster; to be measured |
| [#2](https://github.com/NonDairyNeutrino/better_n_better_n_body/issues/2) | Precompute the recurrence coefficients in `harmonics` | two square roots and a division per entry depend only on `(n, m)` | after 1 |
| [#3](https://github.com/NonDairyNeutrino/better_n_better_n_body/issues/3) | Cache M2L tables by level and relative offset | about a thousand offsets recur for 3.8 million conversions | after 1 and 2; floating-point keys need rounding |
| [#4](https://github.com/NonDairyNeutrino/better_n_better_n_body/issues/4) | Evaluate the gradient directly in L2P | O(p²) per body instead of an O(p⁴) shift followed by reading degree 1 | independent of the others; derivative factors must be checked numerically first |
| [#5](https://github.com/NonDairyNeutrino/better_n_better_n_body/issues/5) | Upgrade M2L to O(p²) per interaction (plane-wave expansions) | diagonal translation operators; an FFT is used for the azimuthal part | needs the standard directional interaction lists; likely pays off only from order 6 or so |
| [#6](https://github.com/NonDairyNeutrino/better_n_better_n_body/issues/6) | Speed up P2P | `powi`, a self-check branch per pair, and each pair done twice | Newton's third law, a direct `1/r³`, contiguous arrays |
| [#7](https://github.com/NonDairyNeutrino/better_n_better_n_body/issues/7) | Reduce the number of M2L conversions per cell | about a thousand per cell now, against at most 189 for the standard interaction list | trades accuracy; the end-to-end error bounds guard it |
| [#8](https://github.com/NonDairyNeutrino/better_n_better_n_body/issues/8) | Reuse buffers, then an arena / struct-of-arrays layout | allocation, cache behaviour, parallelism | a large rewrite, deliberately later |
| [#9](https://github.com/NonDairyNeutrino/better_n_better_n_body/issues/9) | Re-tune the leaf limit and θ; report error next to time in the benchmark | the best leaf size moves whenever the near or far field gets cheaper | the benchmark should compare at equal error |
| [#10](https://github.com/NonDairyNeutrino/better_n_better_n_body/issues/10) | Document the FMM in this README | status, FMM section, input contract and layout are done; the benchmark section reports timings and a separately measured error | finish once the benchmark itself reports error (9) |

## Conventions and assumptions

- Units are nondimensionalized: G = 1.
- The integrator is symplectic Euler, and the tree is rebuilt every time step.
- Barnes-Hut and the FMM do not conserve momentum exactly, because the pull a body feels from a group is not exactly the pull it exerts back.
- The input contract of Barnes-Hut and the FMM is assumed, not checked:
  - all masses are positive;
  - no two bodies share a position (coincident bodies are unphysical);
  - no body lies exactly on a cell boundary;
  - for Barnes-Hut, bodies are at least 1e-10 apart (closer ones merge into one point mass).
- The FMM has no minimum cell size. A group of more than 64 bodies at exactly one position can never be separated, so building the tree recurses until the stack overflows. This is outside the contract above.

## Layout

```
src/
├── main.rs            benchmark driver
├── shared.rs          Vec3, gravity, RNG, test helpers
├── shared/
│   ├── body.rs        Body (mass, position, velocity)
│   └── system.rs      System: a region of space and the bodies in it
├── direct.rs          direct summation
├── barnes_hut.rs      Barnes-Hut time stepping
├── barnes_hut/
│   └── tree.rs        octree construction and acceleration queries
├── fmm.rs             FMM time stepping
└── fmm/
    ├── README.md      the idea, using a school analogy
    ├── PHYSICS.md     the physics in words, with no equations
    ├── complex.rs     complex numbers
    ├── harmonic.rs    solid harmonics and the expansions built from them
    ├── tree.rs        the tree, the traversal, the downward pass, the accelerations
    ├── operators.rs   re-exports the operators below
    └── operators/     one file per operator: p2m, m2m, m2l, l2l, l2p, p2p
```

`PHYSICS.md` in the repository root compares the three methods in words.

## License

GPL-3.0. See [LICENSE](LICENSE).
