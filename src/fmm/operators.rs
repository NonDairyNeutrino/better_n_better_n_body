//! The six FMM operators, which move information between bodies and cells.
//!
//! "Multipole" is a cell's *mass distribution* seen from far away; "local" is the
//! *field* that the rest of the universe makes inside a cell. Both are
//! [`SolidHarmonic`](crate::fmm::harmonic::SolidHarmonic)s about a center point.
//!
//! | Operator                                | Name | Takes                         | Gives                          |
//! | --------------------------------------- | ---- | ----------------------------- | ------------------------------ |
//! | [`particle_to_multipole`]               | P2M  | bodies in a leaf              | mass distribution of the leaf  |
//! | [`multipole_to_multipole`]              | M2M  | mass distribution of a child  | the same, about the parent     |
//! | [`multipole_to_local`]                  | M2L  | mass distribution of a source | field over a far target cell   |
//! | [`local_to_local`]                      | L2L  | field of a parent             | the same field, about a child  |
//! | [`local_to_particle`]                   | L2P  | field of a leaf               | acceleration of each body      |
//! | [`particle_to_particle`]                | P2P  | bodies of a near leaf         | their direct pull on targets   |
//!
//! Operators that move an expansion to a new center take a `shift` vector, and
//! the direction matters. For M2M, `shift` is the old center minus the new one;
//! for L2L it is the new center minus the old one; for M2L it is the target
//! center minus the source center. Every operator works at one fixed order.

pub mod l2l;
pub mod l2p;
pub mod m2l;
pub mod m2m;
pub mod p2m;
pub mod p2p;

pub use l2l::local_to_local;
pub use l2p::local_to_particle;
pub use m2l::multipole_to_local;
pub use m2m::multipole_to_multipole;
pub use p2m::particle_to_multipole;
pub use p2p::particle_to_particle;
