//! Error types for the whole crate.

use std::path::PathBuf;

/// Crate-wide result alias.
pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error(transparent)]
    Io(#[from] std::io::Error),

    /// A malformed line in an input file. `path`/`line` are filled in when known.
    #[error("{}{}", .location.as_ref().map(|l| format!("{l}: ")).unwrap_or_default(), message)]
    Parse {
        message: String,
        location: Option<String>,
    },

    #[error("unsupported input format for {path}: {hint}")]
    UnsupportedFormat { path: PathBuf, hint: String },

    #[error("no {thing} was loaded (checked {checked} input file(s))")]
    EmptyInput { thing: &'static str, checked: usize },

    /// Referenced something by name that the input does not contain.
    #[error("{kind} {name:?} not found {context}")]
    NotFound {
        kind: &'static str,
        name: String,
        context: String,
    },

    #[error(transparent)]
    Tree(#[from] tree::TreeError),
}

impl Error {
    pub(crate) fn parse(message: impl Into<String>) -> Self {
        Error::Parse {
            message: message.into(),
            location: None,
        }
    }
}

/// Newick/tree errors live in [`crate::tree`] and are wrapped by [`Error`].
pub mod tree {
    #[derive(Debug, thiserror::Error)]
    pub enum TreeError {
        #[error("newick parse error at offset {offset}: {message}\n  input: {snippet}")]
        Parse {
            offset: usize,
            message: String,
            snippet: String,
        },

        #[error("tree has no tips")]
        Empty,

        #[error("duplicate tip name {0:?}")]
        DuplicateTip(String),

        #[error("outgroup {0:?} is not a tip of the tree")]
        UnknownOutgroup(String),

        #[error("cannot reroot: {0}")]
        Unrootable(String),
    }
}
