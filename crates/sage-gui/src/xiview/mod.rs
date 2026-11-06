//! Local xiVIEW: the xiVIEW web app (github.com/Rappsilber-Laboratory/xiview)
//! run on this computer with XNova's results.
//!
//! xiVIEW is a browser app that loads its data from the PRIDE crosslinking
//! API (github.com/Rappsilber-Laboratory/crosslinking-api). [`server`] is a
//! small HTTP server on 127.0.0.1 that serves the xiVIEW files and answers
//! that API from an XNova results folder ([`data`]), with annotated spectra
//! in the form xiVIEW's spectrum viewer (xiSPEC) reads ([`annotate`]).

pub mod annotate;
pub mod data;
pub mod server;
