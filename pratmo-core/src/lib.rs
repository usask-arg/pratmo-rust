// The port intentionally preserves several Fortran-shaped loops, slice copies,
// and subroutine signatures so indices remain auditable against the reference
// source. These style lints are therefore disabled at the crate boundary;
// correctness and general compiler warnings remain part of the lint gate.
#![allow(
    clippy::clone_on_copy,
    clippy::manual_memcpy,
    clippy::needless_range_loop,
    clippy::too_many_arguments
)]

pub mod api;
pub mod chemistry;
pub mod clno3;
pub mod constants;
pub mod ctm;
pub mod diurnal;
pub mod heterogeneous;
pub mod init;
pub mod jvalue;
pub mod output;
pub mod path;
pub mod reader;
pub mod solver;
pub mod state;
pub mod tracers;
