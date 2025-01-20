use std::str::FromStr;

use eyre::bail;
use itertools::Itertools;

use crate::Monomer;

use super::{
    hor::{MonomerUnit, HOR},
    token::Token,
};

#[inline]
fn chars2num(chars: impl Iterator<Item = char>) -> eyre::Result<u8> {
    Ok(chars.into_iter().join("").parse::<u8>()?)
}

#[inline]
// https://stackoverflow.com/a/69298721
fn n_digits(num: u8) -> u32 {
    num.checked_ilog10().unwrap_or(0) + 1
}

pub fn hor_monomer_structure_to_monomers<'a>(
    monomers: impl Iterator<Item = &'a MonomerUnit>,
    monomer_base: &Monomer,
) -> Vec<Monomer> {
    let mut new_monomers = vec![];
    let fn_get_new_mon = |m| {
        let mut final_mon = monomer_base.clone();
        final_mon.monomers.push(m);
        final_mon
    };

    for mon in monomers.into_iter() {
        match mon {
            MonomerUnit::Range(range) => {
                if range.end < range.start {
                    let new_range = range.end.saturating_sub(1)..range.start + 1;
                    // First reverse to make range iterable.
                    // Second reverse to restore order.
                    new_monomers.extend(new_range.rev().map(fn_get_new_mon))
                } else {
                    new_monomers.extend(range.clone().map(fn_get_new_mon))
                }
            }
            MonomerUnit::Single(m) => {
                let mut final_mon = monomer_base.clone();
                final_mon.monomers.push(*m);
                new_monomers.push(final_mon);
            }
            MonomerUnit::Chimera(mons) => {
                let mut final_mon = monomer_base.clone();
                final_mon.monomers.extend(mons);
                new_monomers.push(final_mon);
            }
        }
    }
    new_monomers
}

impl FromStr for HOR {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let Some((mon_info, mons)) = s.split('.').collect_tuple::<(&str, &str)>() else {
            bail!("Invalid HOR, {s}. HOR requires monomer info and monomers delimited by '.'")
        };
        let monomers = extract_monomer_order(mons, mon_info)?;
        // Start with base template.
        let mut monomer_base = Monomer::new(&format!("{mon_info}.1"))?;
        monomer_base.monomers.clear();

        let new_monomers = hor_monomer_structure_to_monomers(monomers.iter(), &monomer_base);
        Ok(HOR {
            monomer_structure: monomers,
            monomers: new_monomers,
        })
    }
}

fn extract_monomer_order(mons: &str, mon_info: &str) -> eyre::Result<Vec<MonomerUnit>> {
    let mut ranges = vec![];

    let tokens = &mons.chars().chunk_by(|c| Token::from(*c));
    let mut tokens_iter = tokens.into_iter().peekable();

    let mon_info_len = mon_info.len().try_into()?;
    let mut curr_pos: u32 = mon_info_len;
    while let Some((token, values)) = tokens_iter.next() {
        // Must start with number.
        if token == Token::Number {
            let start_num = chars2num(values.into_iter())?;
            // 45
            curr_pos += n_digits(start_num);

            // Edge-case of 1-monomer.
            if tokens_iter.peek().is_none() {
                ranges.push(MonomerUnit::Single(start_num));
                break;
            }
            let Some((next_token, _)) = tokens_iter.next_if(|(tk, _)| {
                matches!(tk, Token::Chimera | Token::Hyphen | Token::Underscore)
            }) else {
                bail!(
                    "Invalid token ('{}') following number {start_num}, at position {curr_pos}.",
                    tokens_iter
                        .next()
                        .map(|mut t| t.1.join(""))
                        .unwrap_or_default()
                )
            };
            curr_pos += 1;

            match next_token {
                // Case 1: 3/10
                // Chimeric monomers
                Token::Chimera => {
                    let mut chimeric_monomers = vec![start_num];
                    while let Some((chimera_token, chimera_token_vals)) =
                        tokens_iter.next_if(|(tk, _)| matches!(tk, Token::Number | Token::Chimera))
                    {
                        if chimera_token == Token::Chimera {
                            curr_pos += 1;
                            continue;
                        }
                        let num = chars2num(chimera_token_vals.into_iter())?;
                        curr_pos += n_digits(num);
                        chimeric_monomers.push(num);
                    }
                    ranges.push(MonomerUnit::Chimera(chimeric_monomers));
                }
                // Case 2: 1-2
                // Range of monomers.
                Token::Hyphen => {
                    let Some((_, end_num_vals)) =
                        tokens_iter.next_if(|(tk, _)| *tk == Token::Number)
                    else {
                        bail!(
                            "Unexpected token ('{}') at pos {curr_pos}. Expect number after '-'.",
                            tokens_iter
                                .next()
                                .map(|mut t| t.1.join(""))
                                .unwrap_or_default()
                        )
                    };
                    let end_num = chars2num(end_num_vals.into_iter())?;
                    curr_pos += n_digits(end_num);
                    ranges.push(MonomerUnit::Range(start_num..end_num + 1));
                }
                // Case 3: 1_
                // Start of monomer sequence.
                Token::Underscore => {
                    curr_pos += 1;
                    ranges.push(MonomerUnit::Single(start_num));
                }
                _ => unreachable!(),
            }
        } else if token == Token::Underscore && curr_pos != mon_info_len {
            // Do nothing if break in monomer sequence.
            // But don't allow at start.
            curr_pos += 1;
            continue;
        } else {
            bail!(
                "Invalid token ('{}') at {curr_pos}",
                values.into_iter().join("")
            );
        }
    }
    Ok(ranges)
}

#[cfg(test)]
mod test {
    use itertools::Itertools;

    use crate::as_hor::hor::HOR;

    #[test]
    fn test_one_mon_stv() {
        const HOR_SINGLE: &str = "S01/1C3H1L.11";
        let res = HOR::new(HOR_SINGLE).unwrap();
        assert_eq!(format!("{res}"), HOR_SINGLE);
    }

    #[test]
    fn test_invalid_start_stv() {
        const HOR_INV_ST: &str = "S1C10H1L._6/2/4";
        assert!(HOR::new(HOR_INV_ST).is_err());
    }

    #[test]
    fn test_simple_stv() {
        const HOR_SIMPLE: &str = "S01/1C3H1L.11-6";
        let res = HOR::new(HOR_SIMPLE).unwrap();
        assert_eq!(format!("{res}"), HOR_SIMPLE);
    }

    #[test]
    fn test_split_mon_stv() {
        const HOR_SPLIT: &str = "S2C16H2-A.4_7-8";
        let res = HOR::new(HOR_SPLIT).unwrap();
        assert_eq!(format!("{res}"), HOR_SPLIT);
    }

    #[test]
    fn test_long_stv() {
        const HOR: &str = "S2C4H1L.5-14_8-9_3-14_8-9_3-14_8-14_8-9_3-14_8-9_3-14_8-9_3-14_8-9_3-14_8-10_4-14_8-9_3-14_8-14_8-9_3-14_8-9_3-14_8-9_3-14_8-9_3-14_8-9_3-19";
        let res = HOR::new(HOR).unwrap();
        assert_eq!(format!("{res}"), HOR);
    }

    #[test]
    fn test_chim_stv() {
        const HOR_CHIM: &str = "S4CYH1L.46-35_32/34_31/32_31-26_15-1";
        let res = HOR::new(HOR_CHIM).unwrap();
        assert_eq!(format!("{res}"), HOR_CHIM);
    }

    #[test]
    fn test_chim_multi_stv() {
        const HOR_CHIM_MULT: &str = "S1C10H1L.1-5_6/2/4";
        let res = HOR::new(HOR_CHIM_MULT).unwrap();
        assert_eq!(format!("{res}"), HOR_CHIM_MULT);
    }

    #[test]
    fn test_reverse_valid_stv() {
        const HOR_SIMPLE: &str = "S01/1C3H1L.11-6";
        let res = HOR::new(HOR_SIMPLE).unwrap();
        let rev_res = res.reversed();
        assert_eq!(format!("{rev_res}"), "S01/1C3H1L.6-11")
    }

    #[test]
    fn test_reverse_split_mon_stv() {
        const HOR_SPLIT: &str = "S2C16H2-A.4_7-8";
        let res = HOR::new(HOR_SPLIT).unwrap();
        let rev_res = res.reversed();
        assert_eq!(format!("{rev_res}"), "S2C16H2-A.8-7_4")
    }

    #[test]
    fn test_reverse_chim_stv() {
        const HOR_CHIM: &str = "S4CYH1L.46-35_32/34_31/32_31-26_15-1";
        let res = HOR::new(HOR_CHIM).unwrap();
        let rev_res = res.reversed();
        assert_eq!(format!("{rev_res}"), "S4CYH1L.1-15_26-31_32/31_34/32_35-46");
    }

    #[test]
    fn test_iter_hor_stv() {
        const HOR_SIMPLE: &str = "S01/1C3H1L.11-6";
        let res = HOR::new(HOR_SIMPLE).unwrap();
        assert_eq!(
            res.iter().map(|m| format!("{m}")).collect_vec(),
            [
                "S01/1C3H1L.11",
                "S01/1C3H1L.10",
                "S01/1C3H1L.9",
                "S01/1C3H1L.8",
                "S01/1C3H1L.7",
                "S01/1C3H1L.6",
            ]
        )
    }
}
