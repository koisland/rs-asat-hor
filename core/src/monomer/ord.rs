use super::Monomer;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

impl Monomer {
    /// Get right-most mon based on [`Monomer::strand`].
    /// * If not chimeric, return the only number.
    ///
    /// ```
    /// use rs_asat_hor::{Monomer, Strand};
    ///
    /// let mon = Monomer::new("S1C1/5/19H1L.4/6").unwrap();
    /// assert_eq!(mon.right_most_num(), Some(&6));
    /// assert_eq!(mon.with_strand(Strand::Minus).right_most_num(), Some(&4))
    /// ```
    pub fn right_most_num(&self) -> Option<&u8> {
        match self.strand {
            // If (+)
            // 6/4*
            Some(Strand::Plus) => self.monomers.last(),
            // If (-) based on alignment, the order should be inverted.
            // *6/4 -> 4/6*
            Some(Strand::Minus) => self.monomers.first(),
            None => self.monomers.last(),
        }
    }

    /// Get left-most mon based on [`Monomer::strand`].
    /// * If not chimeric, return the only number.
    ///
    /// ```
    /// use rs_asat_hor::{Monomer, Strand};
    ///
    /// let mon = Monomer::new("S1C1/5/19H1L.4/6").unwrap();
    /// assert_eq!(mon.left_most_num(), Some(&4));
    /// assert_eq!(mon.with_strand(Strand::Minus).left_most_num(), Some(&6))
    /// ```
    pub fn left_most_num(&self) -> Option<&u8> {
        match self.strand {
            // If (+)
            // *6/4
            Some(Strand::Plus) => self.monomers.first(),
            // If (-)
            // 6/4* -> *4/6
            Some(Strand::Minus) => self.monomers.last(),
            None => self.monomers.first(),
        }
    }
}

impl PartialOrd for Monomer {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let last_mon = self.right_most_num();
        let first_mon = other.left_most_num();
        let (Some(last_mon), Some(first_mon)) = (last_mon, first_mon) else {
            return None;
        };
        Some(last_mon.cmp(first_mon))
    }
}
