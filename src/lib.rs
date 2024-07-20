//! # PyMOL PSE Parser
//!
//! A Rust crate for working with PyMOL PSE (PyMOL Session) files.
//!
//! ## Features
//!
//! - Load and parse PSE files
//! - Serialize PSE data
//! - Access molecular structures and visualization settings
//!
//! ## Usage
//!
//! ```no_run
//! use pseutils::psedata::PSEData;
//! let psedata = PSEData::load("path/to/file.pse");
//! // Work with the loaded PSE data
//! psedata.to_disk_full("my_output_directory");
//! ```
//!
//! ## Modules
//!
//! - `molviewspec`: Handles molecular viewing specifications
//! - `pymolparsing`: Core parsing functionality for PSE files
//!
pub mod molviewspec;
pub mod pymolparsing;

pub use crate::pymolparsing::psedata::PSEData;
