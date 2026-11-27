# Stage 1: compile. Edition 2024 needs Rust >= 1.85, which rust:1 satisfies.
FROM rust:1-bookworm AS build
WORKDIR /app
COPY Cargo.toml Cargo.lock ./
# main.rs embeds these as crate docs via include_str!, so the build needs them.
COPY README.md PHYSICS.md ./
COPY src ./src
RUN cargo build --release --locked

# Stage 2: run. The project has no dependencies, so the binary only needs libc.
FROM debian:bookworm-slim
COPY --from=build /app/target/release/better-n-better-n-body /usr/local/bin/better-n-better-n-body
ENTRYPOINT ["better-n-better-n-body"]
