//! This crate is for loading and serializing pymol's PSE filetype
//!
//! it is currently a work in progress.
//!
#[macro_use]
extern crate num_derive;

pub mod molviewspec;
pub mod pymolparsing;

pub use crate::pymolparsing::psedata::PSEData;
