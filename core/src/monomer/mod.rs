use std::{fmt::Display, str::FromStr};

use itertools::Itertools;

use chrom::Chromosome;
use mon_type::MonomerType;
use sf::SF;
use status::Status;

mod chrom;
mod mon_type;
mod parse;
mod sf;
mod status;
mod token;

/// An alpha-satellite higher-order repeat monomer.
///
/// ```
/// use rs_asat_hor::Monomer;
///
/// let mon = Monomer::new("S1C16H1L.2");
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Monomer {
    pub monomers: Vec<u8>,
    pub suprachromosomal_family: Vec<SF>,
    pub chromosomes: Vec<Chromosome>,
    pub monomer_type: MonomerType,
    pub monomer_type_desc: Option<String>,
    pub status: Option<Status>,
}

impl Monomer {
    pub fn new(s: &str) -> eyre::Result<Self> {
        Monomer::from_str(s)
    }
}

impl Display for Monomer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let status = match self.status {
            Some(Status::Live) => "L",
            Some(Status::Divergent) => "d",
            None => "",
        };
        let monomers = self.monomers.iter().join("/");
        let monomer_type_desc = self
            .monomer_type_desc
            .as_ref()
            .map(|desc| format!("-{desc}"))
            .unwrap_or_default();
        let chromosomes = self.chromosomes.iter().join("/");
        let sfs = self.suprachromosomal_family.iter().join("/");
        write!(
            f,
            "S{}C{}{:?}{monomer_type_desc}{status}.{monomers}",
            sfs, chromosomes, self.monomer_type,
        )
    }
}
