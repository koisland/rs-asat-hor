use eyre::bail;

use crate::{
    as_hor::{hor_monomer_structure_to_monomers, MonomerUnit},
    Monomer, Strand, HOR,
};

fn get_hor_num(start_mon: Option<&Monomer>, current_num: &u8) -> eyre::Result<MonomerUnit> {
    let Some(Some(start_num)) = start_mon.map(|mon| mon.right_most_num()) else {
        bail!("Start ({start_mon:?}) monomer number not found.")
    };
    // Add previous range or single monomer.
    // If same number, is single monomer.
    if start_num == current_num {
        Ok(MonomerUnit::Single(*start_num))
    } else {
        Ok(MonomerUnit::Range(*start_num..*current_num))
    }
}

/// Convert a sequence of [`Monomer`]s into a [`HOR`].
/// * This assumes that the input sequence has been chunked by strand, gaps, and chrom name.
///
/// ```
/// use rs_asat_hor::{Monomer, Strand, HOR, monomers_to_hor};
///
/// let monomers = [
///     Monomer::new("S1C1/5/19H1L.1").unwrap(),
///     Monomer::new("S1C1/5/19H1L.2").unwrap(),
///     Monomer::new("S1C1/5/19H1L.3").unwrap(),
/// ];
/// let hor = monomers_to_hor(monomers.iter(), Some(Strand::Plus)).unwrap();
/// assert_eq!(
///     hor[0],
///     HOR::new("S1C1/5/19H1L.1-3").unwrap(),
/// )
/// ```
pub fn monomers_to_hor<'a, M>(monomers: M, enforce_strand: Option<Strand>) -> eyre::Result<Vec<HOR>>
where
    M: Iterator<Item = &'a Monomer>,
    M: ExactSizeIterator,
{
    let mut hors = Vec::new();
    if monomers.len() <= 1 {
        return Ok(hors);
    }

    let mut monomers_iter = monomers.into_iter().peekable();
    // Create a base monomer to clone.
    let mut monomer_base = monomers_iter.peek().cloned().unwrap().clone();
    monomer_base.monomers.clear();

    // Bookkeeping vars.
    // Keep track of start and store units. Clear when add new HOR.
    let mut start_mon: Option<&Monomer> = None;
    let mut hor_units: Vec<MonomerUnit> = Vec::new();

    while let Some(mon_1) = monomers_iter.next() {
        // Set start monomer.
        let mon_1_chimeric = mon_1.is_chimeric();
        if start_mon.is_none() {
            start_mon = Some(mon_1)
        }
        let mon_1_num = mon_1.right_most_num().unwrap();

        // Check mon_2.
        let Some(mon_2) = monomers_iter.peek() else {
            // Add remainder once hit end of iterator.
            let final_hor_unit = if mon_1_chimeric {
                MonomerUnit::Chimera(mon_1.monomers.to_vec())
            } else {
                get_hor_num(start_mon, mon_1_num)?
            };
            hor_units.push(final_hor_unit);
            break;
        };
        let mon_2_chimeric = mon_2.is_chimeric();
        let mon_2_num = mon_2.left_most_num().unwrap();

        // Check if want to enforce strand.
        // If not, get strand based on mon_2.
        // If no strand provided, assume forward ort.
        let strand = enforce_strand.unwrap_or_else(|| mon_2.strand.unwrap_or(Strand::Plus));

        // > > x >
        // 5 6 - 1
        let is_gap = mon_1_num.abs_diff(*mon_2_num) > 1;
        let is_broken = match strand {
            // > x > >
            // 6 - 5 6
            Strand::Plus => mon_1 > mon_2,
            // < x < <
            // 5 - 6 5
            Strand::Minus => mon_1 < mon_2,
        };

        if is_gap || is_broken {
            // Case 1: Gap in range.
            // > >  > x  >
            // 1 2 *3 - *6
            // Case 2: Broken range.
            // > >  > x  >
            // 4 5 *6 - *1
            let hor_unit = if mon_1_chimeric {
                MonomerUnit::Chimera(mon_1.monomers.to_vec())
            } else {
                get_hor_num(start_mon, mon_1_num)?
            };
            hor_units.push(hor_unit);
            // Create new HOR and add it.
            let hor = HOR {
                monomer_structure: hor_units.clone(),
                monomers: hor_monomer_structure_to_monomers(hor_units.iter(), &monomer_base),
            };
            hors.push(hor);

            // Clear HOR and start tracking new mon.
            hor_units.clear();
            start_mon.take();
        } else if mon_1_chimeric {
            // Case: Chimeric at start of a HOR.
            //  >  > >
            // 6/4 5 6
            // This should be unreachable in correct data.
            // TODO: Check if anything breaks.
        } else if mon_2_chimeric {
            // Case: Chimeric monomer. End of range.
            hor_units.push(get_hor_num(start_mon, mon_1_num)?);

            // Add chimeric monomer.
            hor_units.push(MonomerUnit::Chimera(mon_2.monomers.to_vec()));

            // Clear start mon.
            // And consumer chimeric mon.
            start_mon.take();
            monomers_iter.next();
        }
    }
    // Add final HOR.
    let monomers = hor_monomer_structure_to_monomers(hor_units.iter(), &monomer_base);
    let hor = HOR {
        monomer_structure: hor_units,
        monomers,
    };
    hors.push(hor);
    Ok(hors)
}

#[cfg(test)]
mod test {
    use itertools::Itertools;

    use crate::{Monomer, HOR};

    use super::{monomers_to_hor, Strand};

    fn assert_hors_equal(
        hors: impl Iterator<Item = HOR>,
        exp_hors: impl Iterator<Item = HOR>,
        string_only: bool,
    ) {
        for (res, exp) in hors.into_iter().zip(exp_hors) {
            if string_only {
                assert_eq!(format!("{res}"), format!("{exp}"))
            } else {
                assert_eq!(res, exp)
            }
        }
    }

    #[test]
    fn test_stv_rev_break() {
        let mons = [
            Monomer::new("S2C13/21H1L.11").unwrap(),
            Monomer::new("S2C13/21H1L.10").unwrap(),
            Monomer::new("S2C13/21H1L.9").unwrap(),
            Monomer::new("S2C13/21H1L.2").unwrap(),
            Monomer::new("S2C13/21H1L.1").unwrap(),
        ];
        assert_hors_equal(
            monomers_to_hor(mons.iter(), Some(Strand::Minus))
                .unwrap()
                .into_iter(),
            [
                HOR::new("S2C13/21H1L.11-9").unwrap(),
                HOR::new("S2C13/21H1L.2-1").unwrap(),
            ]
            .into_iter(),
            false,
        );
    }

    #[test]
    fn test_stv_strand_break() {
        let mons = [
            Monomer::new("S1C1/5/19H1L.1").unwrap(),
            Monomer::new("S1C1/5/19H1L.2").unwrap(),
            Monomer::new("S1C1/5/19H1L.3").unwrap(),
            Monomer::new("S1C1/5/19H1L.4").unwrap(),
            Monomer::new("S1C1/5/19H1L.5").unwrap(),
        ];
        // Going in minus strand breaks order for each step.
        assert_hors_equal(
            monomers_to_hor(mons.iter(), Some(Strand::Minus))
                .unwrap()
                .into_iter(),
            [
                HOR::new("S1C1/5/19H1L.1").unwrap(),
                HOR::new("S1C1/5/19H1L.2").unwrap(),
                HOR::new("S1C1/5/19H1L.3").unwrap(),
                HOR::new("S1C1/5/19H1L.4").unwrap(),
                HOR::new("S1C1/5/19H1L.5").unwrap(),
            ]
            .into_iter(),
            false,
        );
    }

    #[test]
    fn test_stv_break() {
        let mons = [
            Monomer::new("S1C1/5/19H1L.1").unwrap(),
            Monomer::new("S1C1/5/19H1L.2").unwrap(),
            Monomer::new("S1C1/5/19H1L.6").unwrap(),
            Monomer::new("S1C1/5/19H1L.4").unwrap(),
            Monomer::new("S1C1/5/19H1L.5").unwrap(),
        ];
        assert_hors_equal(
            monomers_to_hor(mons.iter(), Some(Strand::Plus))
                .unwrap()
                .into_iter(),
            [
                HOR::new("S1C1/5/19H1L.1-2").unwrap(),
                HOR::new("S1C1/5/19H1L.6").unwrap(),
                HOR::new("S1C1/5/19H1L.4-5").unwrap(),
            ]
            .into_iter(),
            false,
        );
    }

    #[test]
    fn test_stv_no_break() {
        let mons = [
            Monomer::new("S1C1/5/19H1L.1").unwrap(),
            Monomer::new("S1C1/5/19H1L.2").unwrap(),
            Monomer::new("S1C1/5/19H1L.3").unwrap(),
            Monomer::new("S1C1/5/19H1L.4").unwrap(),
            Monomer::new("S1C1/5/19H1L.5").unwrap(),
        ];
        assert_hors_equal(
            monomers_to_hor(mons.iter(), Some(Strand::Plus))
                .unwrap()
                .into_iter(),
            [HOR::new("S1C1/5/19H1L.1-5").unwrap()].into_iter(),
            false,
        );
    }

    #[test]
    fn test_stv_chimera_no_break() {
        let mons = [
            Monomer::new("S1C1/5/19H1L.1").unwrap(),
            Monomer::new("S1C1/5/19H1L.2").unwrap(),
            Monomer::new("S1C1/5/19H1L.3").unwrap(),
            Monomer::new("S1C1/5/19H1L.4").unwrap(),
            Monomer::new("S1C1/5/19H1L.5").unwrap(),
            Monomer::new("S1C1/5/19H1L.6/4").unwrap(),
            Monomer::new("S1C1/5/19H1L.5").unwrap(),
            Monomer::new("S1C1/5/19H1L.6").unwrap(),
        ];
        assert_hors_equal(
            monomers_to_hor(mons.iter(), Some(Strand::Plus))
                .unwrap()
                .into_iter(),
            [HOR::new("S1C1/5/19H1L.1-5_6/4_5-6").unwrap()].into_iter(),
            false,
        );
    }

    #[test]
    fn test_stv_chimera_without_strand() {
        let mons = [
            Monomer::new("S1C1/5/19H1L.5").unwrap(),
            Monomer::new("S1C1/5/19H1L.6/4").unwrap(),
            Monomer::new("S1C1/5/19H1L.5").unwrap(),
            Monomer::new("S1C1/5/19H1L.6/4").unwrap(),
        ];
        assert_hors_equal(
            monomers_to_hor(mons.iter(), Some(Strand::Minus))
                .unwrap()
                .into_iter(),
            [
                HOR::new("S1C1/5/19H1L.5").unwrap(),
                HOR::new("S1C1/5/19H1L.6/4").unwrap(),
                HOR::new("S1C1/5/19H1L.5").unwrap(),
                HOR::new("S1C1/5/19H1L.6/4").unwrap(),
            ]
            .into_iter(),
            false,
        );
    }
    #[test]
    fn test_stv_chimera_with_strand() {
        let mons = [
            Monomer::new("S1C1/5/19H1L.5").unwrap(),
            Monomer::new("S1C1/5/19H1L.6/4").unwrap(),
            Monomer::new("S1C1/5/19H1L.5").unwrap(),
            Monomer::new("S1C1/5/19H1L.6/4").unwrap(),
        ];
        let rev_mons: Vec<Monomer> = mons
            .into_iter()
            .map(|mon| mon.with_strand(Strand::Minus))
            .collect_vec();
        assert_hors_equal(
            monomers_to_hor(rev_mons.iter(), Some(Strand::Minus))
                .unwrap()
                .into_iter(),
            [HOR::new("S1C1/5/19H1L.5_6/4_5_6/4").unwrap()].into_iter(),
            // Because HOR is read-only. We cannot remove monomer strand information.
            // We only compare the formatted string instead.
            true,
        );
    }
}
