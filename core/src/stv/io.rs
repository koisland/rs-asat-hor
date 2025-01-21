use std::{collections::HashMap, io::BufRead, path::Path, str::FromStr};

use itertools::Itertools;

use crate::{Monomer, Strand, HOR};

use super::monomers_to_hor;

/// A `BED9` HOR monomer record.
/// ```no_run
///
/// let record = ("chr1", 1, 170, "S1C1/5/19H1L.6", 100.0, "+", 1, 170, "0,0,0");
/// ```
pub type MonomerRecord<'a> = (&'a str, u64, u64, &'a str, f32, &'a str, u64, u64, &'a str);

/// An `BED4` HOR structural variation record.
/// ```no_run
/// use rs_asat_hor::HOR;
///
/// let record = (String::from("chr1"), 1, 1020, HOR::new("S1C1/5/19H1L.1-6"));
/// ```
pub type StvRecord = (String, u64, u64, HOR);

/// Read a `BED9` file of [`MonomerRecord`]s and convert them to [`StvRecord`]s.
///
/// # Args
/// * `bedfile`
///     * Path to `BED9` file.
/// * `fn_filter`
///     * Function to filter records if `true`.
///     * A noop can be achieved with `|_| true`
///
/// # Returns
/// * Valid [`StvRecord`]s
///
/// # Examples
/// Filter monomers that have identity greater than `85.0`.
/// ```
/// use rs_asat_hor::{read_from_monomer_bed, MonomerRecord};
///
/// let records = read_from_monomer_bed(
///     "test/mons.bed",
///     |rec: MonomerRecord| rec.4 < 85.0
/// );
/// assert!(records.is_ok())
/// ```
pub fn read_from_monomer_bed<F>(
    bedfile: impl AsRef<Path>,
    fn_filter: F,
) -> eyre::Result<Vec<StvRecord>>
where
    F: Fn(MonomerRecord) -> bool,
{
    let file = std::fs::File::open(bedfile).unwrap();
    let fh = std::io::BufReader::new(file);
    let mut records: Vec<StvRecord> = vec![];

    let mut chr_mons: HashMap<String, Vec<(u64, u64, Monomer)>> = HashMap::new();
    for line in fh.lines() {
        let line = line?;
        let Some((chrom, st, end, name, score, ort, tst, tend, rgb)) =
            line.trim().split('\t').collect_tuple()
        else {
            continue;
        };
        let (st, end, score, tst, tend) = (
            st.parse::<u64>()?,
            end.parse::<u64>()?,
            score.parse::<f32>()?,
            tst.parse::<u64>()?,
            tend.parse::<u64>()?,
        );

        // Allow filter function.
        if fn_filter((chrom, st, end, name, score, ort, tst, tend, rgb)) {
            continue;
        }

        // Add strand.
        let strand = Strand::from_str(ort)?;
        if let Ok(mon) = Monomer::new(name).map(|m| m.with_strand(strand)) {
            if let Some(mons) = chr_mons.get_mut(chrom) {
                mons.push((st, end, mon));
            } else {
                chr_mons
                    .entry(chrom.to_owned())
                    .or_insert(vec![(st, end, mon)]);
            }
        } else {
            log::error!("Cannot convert monomer ({name}) at {chrom}:{st}-{end}. Skipping.");
        }
    }
    for (chrom, mons) in chr_mons.iter() {
        // Convert monomers in chromosome to HOR.
        // We don't enforce strand here or chunk to avoid breaking HORs.
        let hors = monomers_to_hor(mons.iter().map(|m| &m.2), None)?;

        // Keep track of monomer index positions with cumulative sum of indices.
        // ex.
        //    mon: 1 2 3 7 8
        //    hor: 0 0 0 1 1
        //    idx: 0 1 2 3 4
        // res.
        //    [0, 3, 5]
        let mut idxs_mon = vec![0; hors.len() + 1];

        for (i, idx_mon) in hors
            .iter()
            .map(|h| h.n_monomers())
            .enumerate()
            // Offset by 1 for starting position 0.
            .map(|(i, m)| (i + 1, m))
        {
            // Safe as always i < idxs_mon.
            let idx_mon_offset = idxs_mon.get(i - 1).unwrap();
            idxs_mon[i] = idx_mon + idx_mon_offset
        }

        // Convert to idx intervals.
        // ex.  [0, 3, 5]
        // res. (0, 3), (3, 5)
        for ((st, end), hor) in idxs_mon
            .into_iter()
            .tuple_windows::<(usize, usize)>()
            .zip(hors.into_iter())
        {
            let Some(mons) = mons.get(st..end) else {
                continue;
            };
            // Find min and max coordinates of HOR.
            let mut min_st = u64::MAX;
            let mut max_end = 0;
            for (st, end, _) in mons {
                min_st = std::cmp::min(min_st, *st);
                max_end = std::cmp::max(max_end, *end);
            }
            assert!(
                min_st != u64::MAX,
                "Logic error with indexing with {chrom}:{st}-{end} and {hor}. Report on GitHub issue tracker."
            );
            records.push((chrom.to_string(), min_st, max_end, hor));
        }
    }
    Ok(records)
}

#[cfg(test)]
mod test {
    use crate::{read_from_monomer_bed, Monomer, Strand, HOR};

    #[test]
    fn test_read_mon_bed() {
        let records = read_from_monomer_bed("test/mons.bed", |_| false).unwrap();

        let mut monomers: Vec<Monomer> = (1..12)
            .map(|i| {
                Monomer::new(&format!("S2C15H1L.{i}"))
                    .unwrap()
                    .with_strand(Strand::Minus)
            })
            .rev()
            .collect();
        monomers.extend(monomers.clone());
        let exp_hor = HOR::from_monomers(&monomers).unwrap();

        assert_eq!(
            records,
            vec![
                (
                    String::from("chm1_chr15:3977696-8919402"),
                    2732060,
                    2733936,
                    exp_hor[0].clone()
                ),
                (
                    String::from("chm1_chr15:3977696-8919402"),
                    2734617,
                    2736493,
                    exp_hor[1].clone()
                ),
            ]
        );
    }

    #[test]
    fn test_read_mon_bed_filter() {
        // Filter less than 99% identity.
        let records = read_from_monomer_bed("test/mons.bed", |row| row.4 < 99.0).unwrap();

        let mut monomers: Vec<Monomer> = (4..12)
            .map(|i| {
                Monomer::new(&format!("S2C15H1L.{i}"))
                    .unwrap()
                    .with_strand(Strand::Minus)
            })
            .rev()
            .collect();
        monomers.extend(monomers.clone());
        monomers.push(
            Monomer::new("S2C15H1L.1")
                .unwrap()
                .with_strand(Strand::Minus),
        );
        let exp_hor = HOR::from_monomers(&monomers).unwrap();

        assert_eq!(
            records,
            vec![
                (
                    String::from("chm1_chr15:3977696-8919402"),
                    2732060,
                    2733427,
                    exp_hor[0].clone()
                ),
                (
                    String::from("chm1_chr15:3977696-8919402"),
                    2734617,
                    2735984,
                    exp_hor[1].clone()
                ),
                (
                    String::from("chm1_chr15:3977696-8919402"),
                    2736323,
                    2736493,
                    exp_hor[2].clone()
                ),
            ]
        );
    }
}
