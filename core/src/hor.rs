use std::{
    fmt::Display,
    ops::{Deref, Range},
    str::FromStr,
};

use eyre::bail;
use itertools::Itertools;

use crate::monomer::Monomer;

#[derive(Debug, Clone, PartialEq, Eq)]
enum Token {
    Number,
    Underscore,
    Hyphen,
    Chimera,
    Other(char),
}

impl From<char> for Token {
    fn from(value: char) -> Self {
        match value {
            '0'..='9' => Token::Number,
            '_' => Token::Underscore,
            '-' => Token::Hyphen,
            '/' => Token::Chimera,
            _ => Token::Other(value),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MonomerNumber {
    Range(Range<u8>),
    Single(u8),
    Chimera(Vec<u8>),
}

#[inline]
fn chars2num(chars: impl Iterator<Item = char>) -> eyre::Result<u8> {
    Ok(chars.into_iter().join("").parse::<u8>()?)
}

#[inline]
// https://stackoverflow.com/a/69298721
fn n_digits(num: u8) -> u32 {
    num.checked_ilog10().unwrap_or(0) + 1
}

fn extract_monomer_order(mons: &str, mon_info: &str) -> eyre::Result<Vec<MonomerNumber>> {
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
                ranges.push(MonomerNumber::Single(start_num));
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
                    ranges.push(MonomerNumber::Chimera(chimeric_monomers));
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
                    ranges.push(MonomerNumber::Range(start_num..end_num + 1));
                }
                // Case 3: 1_
                // Start of monomer sequence.
                Token::Underscore => {
                    curr_pos += 1;
                    ranges.push(MonomerNumber::Single(start_num));
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

/// An alpha-satellite higher-order repeat composed of [`Monomer`]s.
#[derive(Debug, Clone)]
pub struct HOR {
    monomer_structure: Vec<MonomerNumber>,
    monomers: Vec<Monomer>,
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

        let fn_get_new_mon = |m| {
            let mut final_mon = monomer_base.clone();
            final_mon.monomers.push(m);
            final_mon
        };
        let mut new_monomers = vec![];
        for mon in monomers.iter() {
            match mon {
                MonomerNumber::Range(range) => {
                    if range.end < range.start {
                        let new_range = range.end.saturating_sub(1)..range.start + 1;
                        // First reverse to make range iterable.
                        // Second reverse to restore order.
                        new_monomers.extend(new_range.rev().map(fn_get_new_mon))
                    } else {
                        new_monomers.extend(range.clone().map(fn_get_new_mon))
                    }
                }
                MonomerNumber::Single(m) => {
                    let mut final_mon = monomer_base.clone();
                    final_mon.monomers.push(*m);
                    new_monomers.push(final_mon);
                }
                MonomerNumber::Chimera(mons) => {
                    let mut final_mon = monomer_base.clone();
                    final_mon.monomers.extend(mons);
                    new_monomers.push(final_mon);
                }
            }
        }
        Ok(HOR {
            monomer_structure: monomers,
            monomers: new_monomers,
        })
    }
}

impl HOR {
    /// Generate a new [`HOR`] from an input string.
    ///
    /// ```
    /// use rs_asat_hor::HOR;
    ///
    /// let hor = HOR::new("S01/1C3H1L.11-6").unwrap();
    /// ```
    pub fn new(s: &str) -> eyre::Result<Self> {
        HOR::from_str(s)
    }

    /// Geneate a new [`HOR`] from monomers.
    pub fn from_monomers(_monomers: &[Monomer]) -> Self {
        // impl from iter.
        todo!()
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
                MonomerNumber::Range(range) => {
                    MonomerNumber::Range(range.end.saturating_sub(1)..range.start + 1)
                }
                MonomerNumber::Chimera(monomers) => {
                    MonomerNumber::Chimera(monomers.iter().rev().cloned().collect())
                }
                MonomerNumber::Single(_) => m.clone(),
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
                MonomerNumber::Range(range) => {
                    write!(f, "{}-{}", range.start, range.end.saturating_sub(1))?;
                }
                MonomerNumber::Single(mon) => {
                    write!(f, "{mon}")?;
                }
                MonomerNumber::Chimera(monomers) => {
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

#[cfg(test)]
mod test {
    use itertools::Itertools;

    use crate::hor::HOR;

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
