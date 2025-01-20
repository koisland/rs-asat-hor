use std::io::BufRead;

use eyre::bail;
use itertools::Itertools;

use crate::{hor::MonomerNumber, Monomer, HOR};

#[derive(Debug, Clone)]
pub enum Strand {
    Plus,
    Minus,
}

/// Convert a sequence of [`Monomer`]s to [`HOR`]s.
/// This assumes that the input has been chunked by strand, gaps, and chrom name.
fn monomer_to_hor(monomers: &[Monomer], strand: Strand) -> eyre::Result<Vec<HOR>> {
    let mut hors = vec![];

    let Some(first_monomer) = monomers.first() else {
        return Ok(hors);
    };
    if monomers.len() == 1 {
        let monomer = first_monomer.clone();
        let monomer_structure = if monomer.monomers.len() > 1 {
            MonomerNumber::Single(monomer.monomers[0])
        } else {
            MonomerNumber::Chimera(monomer.monomers.to_vec())
        };
        hors.push(HOR {
            monomer_structure: vec![monomer_structure],
            monomers: vec![monomer],
        });
        return Ok(hors);
    }

    let fn_brkpt_detect = match strand {
        Strand::Plus => |a: u8, b: u8| a < b,
        Strand::Minus => |a: u8, b: u8| a > b,
    };

    let mut monomers_iter = monomers.iter().peekable();
    let mut monomer_units: Vec<Monomer> = vec![];
    let mut start_idx = 0;
    while let Some(mon_1) = monomers_iter.next() {
        let Some(mon_2) = monomers_iter.next() else {
            // Add mon_1
            break;
        };

        let (Some(mon_1_id), Some(mon_2_id)) = (mon_1.monomers.last(), mon_2.monomers.first())
        else {
            bail!("One or more monomers lacks an ID. ({mon_1}, {mon_2})")
        };
    }

    Ok(hors)
}

#[test]
fn test_read_mon_bed() {
    let file = std::fs::File::open("test/mon.bed").unwrap();
    let fh = std::io::BufReader::new(file);

    let mut mons = vec![];
    for line in fh.lines() {
        let line = line.unwrap();
        let Some((chrom, st, end, name, score, ort, tst, tend, rgb)) =
            line.trim().split('\t').collect_tuple()
        else {
            continue;
        };

        if let Ok(mon) = Monomer::new(name) {
            mons.push(mon);
        } else {
            dbg!(name);
        }
    }
}

#[test]
fn test_read_hor_bed() {
    let file = std::fs::File::open("test/stv_all.bed").unwrap();
    let fh = std::io::BufReader::new(file);

    let mut hors = vec![];
    for line in fh.lines() {
        let line = line.unwrap();
        let Some((chrom, st, end, name, score, ort, tst, tend, rgb)) =
            line.trim().split('\t').collect_tuple()
        else {
            continue;
        };

        if let Ok(hor) = HOR::new(name) {
            hors.push(hor);
        } else {
            dbg!(name);
        }
    }
}

#[cfg(test)]
mod test {
    use crate::HOR;

    #[test]
    fn test() {
        let hor = HOR::new("S1C10H1L.1-5_6/2/4").unwrap();
        println!("{hor}")
    }
}
