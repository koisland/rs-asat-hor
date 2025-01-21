use std::{
    fmt::Display,
    ops::{Deref, Range},
    str::FromStr,
};

use itertools::Itertools;

use crate::{monomer::Monomer, monomers_to_hor};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MonomerUnit {
    Range(Range<u8>),
    Single(u8),
    Chimera(Vec<u8>),
}

/// An alpha-satellite higher-order repeat composed of one or more [`Monomer`]s.
/// ```
/// use rs_asat_hor::HOR;
///
/// // A 6-monomer chr3, SF01/1 HOR.
/// let hor = HOR::new("S01/1C3H1L.11-6").unwrap();
/// assert_eq!(hor.len(), 6)
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HOR {
    pub(crate) monomer_structure: Vec<MonomerUnit>,
    pub(crate) monomers: Vec<Monomer>,
}

impl HOR {
    /// Generate a new [`HOR`] from an input string.
    ///
    /// ```
    /// use rs_asat_hor::HOR;
    ///
    /// let hor = HOR::new("S01/1C3H1L.11-6");
    /// assert!(hor.is_ok());
    /// ```
    pub fn new(s: &str) -> eyre::Result<Self> {
        HOR::from_str(s)
    }

    /// Get the number of monomers in a [`HOR`].
    ///
    /// ```
    /// use rs_asat_hor::HOR;
    ///
    /// let hor = HOR::new("S01/1C3H1L.11-6").unwrap();
    /// assert_eq!(hor.n_monomers(), 6)
    /// ```
    pub fn n_monomers(&self) -> usize {
        self.monomers.len()
    }

    /// Get all [`Monomer`]s within this [`HOR`].
    ///
    /// ```
    /// use rs_asat_hor::{HOR, Monomer};
    ///
    /// let hor = HOR::new("S01/1C3H1L.11-9").unwrap();
    /// assert_eq!(
    ///     hor.monomers(),
    ///     &[
    ///         Monomer::new("S01/1C3H1L.11").unwrap(),
    ///         Monomer::new("S01/1C3H1L.10").unwrap(),
    ///         Monomer::new("S01/1C3H1L.9").unwrap(),
    ///     ]
    /// )
    /// ```
    pub fn monomers(&self) -> &[Monomer] {
        &self.monomers[..]
    }

    /// Get the [`MonomerUnit`]s within this [`HOR`].
    ///
    /// ```
    /// use rs_asat_hor::{HOR, MonomerUnit};
    ///
    /// let hor = HOR::new("S01/1C3H1L.7-9_11").unwrap();
    /// assert_eq!(
    ///     hor.monomer_units(),
    ///     &[MonomerUnit::Range(7..9), MonomerUnit::Single(11)]
    /// )
    /// ```
    pub fn monomer_units(&self) -> &[MonomerUnit] {
        &self.monomer_structure[..]
    }

    /// Generate [`HOR`]s from monomers.
    /// * Convenience function for [`crate::monomers_to_hor`].
    ///
    /// ```
    /// use rs_asat_hor::{HOR, Monomer};
    ///
    /// let mons = [
    ///     Monomer::new("S1C1/5/19H1L.1").unwrap(),
    ///     Monomer::new("S1C1/5/19H1L.2").unwrap(),
    ///     Monomer::new("S1C1/5/19H1L.3").unwrap(),
    /// ];
    /// let hor = HOR::from_monomers(&mons).unwrap();
    /// assert_eq!(
    ///     format!("{}", hor[0]),
    ///     "S1C1/5/19H1L.1-3"
    /// )
    /// ```
    pub fn from_monomers(monomers: &[Monomer]) -> eyre::Result<Vec<Self>> {
        monomers_to_hor(monomers.iter(), None)
    }

    /// Generate the reversed version of this [`HOR`].
    ///
    /// ```
    /// use rs_asat_hor::HOR;
    ///
    /// let hor = HOR::new("S01/1C3H1L.11-6").unwrap();
    /// let rev_hor = hor.reversed();
    ///
    /// assert_eq!(format!("{rev_hor}"), "S01/1C3H1L.6-11");
    /// ```
    pub fn reversed(&self) -> Self {
        let new_monomer_structure = self
            .monomer_structure
            .iter()
            .rev()
            .map(|m| match m {
                MonomerUnit::Range(range) => MonomerUnit::Range(range.end..range.start),
                MonomerUnit::Chimera(monomers) => {
                    MonomerUnit::Chimera(monomers.iter().rev().cloned().collect())
                }
                MonomerUnit::Single(_) => m.clone(),
            })
            .collect_vec();
        let new_monomers = self
            .monomers
            .iter()
            .cloned()
            .rev()
            .map(|mut mon| {
                // Swap chimeric monomer if present.
                mon.monomers.reverse();
                mon
            })
            .collect_vec();
        Self {
            monomer_structure: new_monomer_structure,
            monomers: new_monomers,
        }
    }
}

// https://stackoverflow.com/a/70547964
impl IntoIterator for HOR {
    type Item = Monomer;
    type IntoIter = <Vec<Self::Item> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        self.monomers.into_iter()
    }
}

impl Deref for HOR {
    type Target = [Monomer];

    fn deref(&self) -> &[Monomer] {
        &self.monomers[..]
    }
}

impl Display for HOR {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Some(mon_1) = self.monomers.first() else {
            write!(f, "")?;
            return Ok(());
        };
        let mon_elems_str = format!("{mon_1}");

        let Some((mon_elems_slice, _)) = mon_elems_str.split_once('.') else {
            unreachable!("Safe. Should always have . at this point.")
        };

        // Write monomer information.
        write!(f, "{mon_elems_slice}.")?;

        for (i, mon_order) in self.monomer_structure.iter().enumerate() {
            match mon_order {
                MonomerUnit::Range(range) => {
                    write!(f, "{}-{}", range.start, range.end)?;
                }
                MonomerUnit::Single(mon) => {
                    write!(f, "{mon}")?;
                }
                MonomerUnit::Chimera(monomers) => {
                    write!(f, "{}", monomers.iter().join("/"))?;
                }
            }
            if i != self.monomer_structure.len() - 1 {
                write!(f, "_")?;
            }
        }
        Ok(())
    }
}
