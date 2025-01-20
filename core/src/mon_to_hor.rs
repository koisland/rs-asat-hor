use std::io::BufRead;

use itertools::Itertools;

use crate::{Monomer, HOR};

fn monomer_to_hor(monomers: &[Monomer]) -> Vec<HOR> {
    let hors = vec![];

    for monomer in monomers {}

    hors
}

#[test]
fn test_read_bed() {
    let file = std::fs::File::open("test/stv_all.bed").unwrap();
    let fh = std::io::BufReader::new(file);

    let mut monomers = vec![];
    for line in fh.lines() {
        let line = line.unwrap();
        let Some((chrom, st, end, name, score, ort, tst, tend, rgb)) =
            line.trim().split('\t').collect_tuple()
        else {
            continue;
        };

        if let Ok(monomer) = HOR::new(name) {
            monomers.push(monomer);
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
