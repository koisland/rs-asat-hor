use std::{fmt::Display, str::FromStr};

use itertools::Itertools;

use super::{chrom::Chromosome, mon_type::MonomerType, ord::Strand, sf::SF, status::Status};

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
    pub strand: Option<Strand>,
}

impl Monomer {
    /// Construct a new [`Monomer`] from a given string.
    ///
    /// ```
    /// use rs_asat_hor::Monomer;
    ///
    /// let mon = Monomer::new("S1C16H1L.2");
    /// assert!(mon.is_ok());
    /// ```
    pub fn new(s: &str) -> eyre::Result<Self> {
        Monomer::from_str(s)
    }

    /// Add [`Strand`] information. Affects chimeric monomer comparison with other [`Monomer`]s.
    /// * If omitted, the default ordering (assumed `+`) is retained.
    /// * Does not alter [`Monomer::monomers`].
    /// ```
    /// use rs_asat_hor::{Monomer, Strand};
    ///
    /// let mon1 = Monomer::new("S1C1/5/19H1L.5").unwrap();
    /// let mon2 = Monomer::new("S1C1/5/19H1L.4/6").unwrap();
    /// let mut mon2_inv = mon2.clone().with_strand(Strand::Minus);
    /// // If (+): 5 > 4/6
    /// assert!(mon1 > mon2);
    /// // If (-): 5 < 6/4
    /// assert!(mon1 < mon2_inv);
    /// ```
    pub fn with_strand(mut self, strand: Strand) -> Self {
        self.strand = Some(strand);
        self
    }

    /// Check if this [`Monomer`] is hybrid/chimeric and contains multiple numbers.
    ///
    /// ```
    /// use rs_asat_hor::Monomer;
    ///
    /// let mon1 = Monomer::new("S1C1/5/19H1L.5").unwrap();
    /// let mon2 = Monomer::new("S1C1/5/19H1L.4/6").unwrap();
    /// assert!(!mon1.is_chimeric());
    /// assert!(mon2.is_chimeric());
    /// ```
    pub fn is_chimeric(&self) -> bool {
        self.monomers.len() > 1
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
