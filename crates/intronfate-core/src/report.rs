//! Non-fatal quality diagnostics.
//!
//! Annotation files in the wild are rarely clean. Rather than failing the whole run on the first
//! odd transcript, ingestion records structured diagnostics and the CLI reports them at the end
//! together with the offending rows, so the projection tables stay usable for a 20k-gene genome.

use std::fmt;

/// How loud a diagnostic is.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub enum Level {
    Info,
    Warn,
    Error,
}

impl fmt::Display for Level {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Level::Info => "info",
            Level::Warn => "warn",
            Level::Error => "error",
        })
    }
}

/// One diagnostic: a kind so callers can filter, a subject so rows can be traced back.
#[derive(Clone, Debug)]
pub struct Diagnostic {
    pub level: Level,
    /// Stable machine-readable code, e.g. `frame_mismatch`, `no_protein`.
    pub kind: &'static str,
    /// Transcript / protein / file the message is about.
    pub subject: String,
    pub message: String,
}

/// Collects diagnostics, keeping full counts but capping retained rows so a pathological input
/// cannot exhaust memory.
#[derive(Clone, Debug, Default)]
pub struct Diagnostics {
    counts: Vec<(&'static str, u64)>,
    retained: Vec<Diagnostic>,
    cap: usize,
    pub errors: u64,
    pub warnings: u64,
    pub infos: u64,
}

impl Diagnostics {
    pub fn new(cap: usize) -> Self {
        Diagnostics {
            counts: Vec::new(),
            retained: Vec::new(),
            cap: cap.max(1),
            ..Default::default()
        }
    }

    fn bump_count(&mut self, kind: &'static str) {
        if let Some(e) = self.counts.iter_mut().find(|(k, _)| *k == kind) {
            e.1 += 1;
        } else {
            self.counts.push((kind, 1));
        }
    }

    pub fn push(
        &mut self,
        level: Level,
        kind: &'static str,
        subject: impl Into<String>,
        message: impl Into<String>,
    ) {
        let subject = subject.into();
        let message = message.into();
        match level {
            Level::Info => self.infos += 1,
            Level::Warn => self.warnings += 1,
            Level::Error => self.errors += 1,
        }
        self.bump_count(kind);
        if level >= Level::Warn && self.retained.len() < self.cap {
            self.retained.push(Diagnostic {
                level,
                kind,
                subject,
                message,
            });
        }
    }

    pub fn warn(
        &mut self,
        kind: &'static str,
        subject: impl Into<String>,
        message: impl Into<String>,
    ) {
        self.push(Level::Warn, kind, subject, message);
    }

    pub fn error(
        &mut self,
        kind: &'static str,
        subject: impl Into<String>,
        message: impl Into<String>,
    ) {
        self.push(Level::Error, kind, subject, message);
    }

    pub fn info(
        &mut self,
        kind: &'static str,
        subject: impl Into<String>,
        message: impl Into<String>,
    ) {
        self.push(Level::Info, kind, subject, message);
    }

    pub fn retained(&self) -> &[Diagnostic] {
        &self.retained
    }

    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    /// `(kind, count)` sorted by descending count then name — this is the run's quality summary.
    pub fn tally(&self) -> Vec<(&'static str, u64)> {
        let mut c = self.counts.clone();
        c.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));
        c
    }

    pub fn count_of(&self, kind: &str) -> u64 {
        self.counts
            .iter()
            .find(|(k, _)| *k == kind)
            .map(|(_, n)| *n)
            .unwrap_or(0)
    }

    pub fn merge(&mut self, other: &Diagnostics) {
        for (k, n) in &other.counts {
            for _ in 0..*n {
                self.bump_count(k);
            }
        }
        self.errors += other.errors;
        self.warnings += other.warnings;
        self.infos += other.infos;
        for d in &other.retained {
            if self.retained.len() < self.cap {
                self.retained.push(d.clone());
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_and_cap() {
        let mut d = Diagnostics::new(2);
        for i in 0..10 {
            d.warn("no_protein", format!("t{i}"), "missing");
        }
        d.error("frame_mismatch", "tX", "bad");
        assert_eq!(d.warnings, 10);
        assert_eq!(d.errors, 1);
        assert_eq!(d.count_of("no_protein"), 10);
        assert_eq!(d.retained().len(), 2, "retained rows must respect the cap");
        assert_eq!(d.tally()[0], ("no_protein", 10), "tally sorted by count");
    }
}
