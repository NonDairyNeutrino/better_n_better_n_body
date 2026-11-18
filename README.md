# better-n-better-n-body

A showcase of better and better algorithms for simulating gravitational N-body systems, written in Rust. Each algorithm trades a little accuracy for a lot of speed, and the code is meant to be read in order, from the simplest to the most sophisticated.

## Algorithms

| Algorithm      | Cost per step | Status  |
| -------------- | ------------- | ------- |
| Direct         | O(N²)         | Done    |
| Barnes-Hut     | O(N log N)    | Done    |
| Fast Multipole | O(N)          | Planned |

### Direct
Sums the pull of every body on every other body. It is exact (up to the integrator) and serves as the reference the other methods are checked against. See `src/direct.rs`.

### Barnes-Hut
Builds an octree over the bodies each step. A distant cell is treated as a single point mass at its center of mass, so the force on a body needs only about log N interactions instead of N. See `src/barnes_hut.rs` and `src/barnes_hut/tree.rs`.

The accuracy knob is the opening angle θ (called `resolution` in the code). A cell is treated as one point mass when `cell width < θ × distance to its center of mass`. At θ = 0 no cell is ever approximated, so the result matches direct summation. Larger θ is faster and less accurate.

## Running

```sh
cargo run --release   # times direct vs. Barnes-Hut at several θ values
cargo test            # includes Barnes-Hut (θ = 0) vs. direct checks
```

`main` simulates 10,000 random bodies and prints the wall-clock time for direct summation and for Barnes-Hut at θ = 0.4, 0.6, 0.8 and 1.0. Use `--release`: debug builds are far too slow to say anything meaningful about performance. The parameters (body count, time step, total time) are constants at the top of `main` in `src/main.rs`.

## Benchmarks

10,000 bodies, 100 steps (`dt` = 0.01, total time 1.0), release build, single run on one machine:

| Method              | Time    | Speedup vs. direct |
| ------------------- | ------- | ------------------ |
| Direct              | 39.67 s | 1.0×               |
| Barnes-Hut, θ = 0.4 | 23.74 s | 1.7×               |
| Barnes-Hut, θ = 0.6 | 10.78 s | 3.7×               |
| Barnes-Hut, θ = 0.8 | 6.45 s  | 6.1×               |
| Barnes-Hut, θ = 1.0 | 4.73 s  | 8.4×               |

These are timings only; accuracy versus θ is not measured here.

## Conventions and assumptions

- Units are nondimensionalized: G = 1.
- The integrator is symplectic Euler, and the tree is rebuilt every time step.
- The Barnes-Hut input contract is assumed, not checked:
  - all masses are positive;
  - no two bodies share a position;
  - bodies are at least 1e-10 apart (closer ones merge into one point mass);
  - no body lies exactly on a cell boundary.

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
└── barnes_hut/
    └── tree.rs        octree construction and acceleration queries
```

## License

GPL-3.0. See [LICENSE](LICENSE).
