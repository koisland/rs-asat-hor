use std::{fmt::Display, ops::Range, str::FromStr};

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
    Chimera(u8, u8),
}

#[inline]
fn chars2num(chars: impl Iterator<Item = char>) -> eyre::Result<u8> {
    Ok(chars.into_iter().join("").parse::<u8>()?)
}

fn extract_monomer_order(s: &str) -> eyre::Result<Vec<MonomerNumber>> {
    let mut start_token: Option<Token> = None;
    let mut start_num: Option<u8> = None;
    let mut prev_token: Option<Token> = None;
    let mut ranges = vec![];

    for (token, grps) in &s.chars().chunk_by(|c| Token::from(*c)) {
        // 4_7-8
        // 46-35_32/34_31/32_31-26_15-1
        match (&start_token, &prev_token, &token) {
            (None, None, Token::Number) => {
                start_token = Some(token);
                start_num = Some(chars2num(grps.into_iter())?);
            }
            // Skip if new state.
            // _
            (None, None, Token::Underscore) => {
                continue;
            }
            (None, None, _) => {
                bail!("Invalid starting token, {token:?}.");
            }
            (None, Some(_), _) => {
                unreachable!("{start_token:?},{prev_token:?},{token:?}")
            }
            // Start is single mon.
            // 1_
            (Some(Token::Number), None, Token::Underscore) => {
                // Add 1-length range.
                ranges.push(MonomerNumber::Single(start_num.unwrap()));
                // Reset start.
                start_token.take();
                start_num.take();
            }
            // Continuation of starting mon.
            // 1-
            // Start is chimeric mon.
            // 1/
            (Some(Token::Number), None, _) => prev_token = Some(token),
            (Some(_), None, _) => {
                unreachable!("{start_token:?},{prev_token:?},{token:?}")
            }
            // Range
            // 1-12
            (Some(Token::Number), Some(Token::Hyphen), Token::Number) => {
                let start_num_end = chars2num(grps.into_iter())?;
                ranges.push(MonomerNumber::Range(start_num.unwrap()..start_num_end + 1));
                start_token.take();
                start_num.take();
                prev_token.take();
            }
            // Chimera
            // 3/10
            (Some(Token::Number), Some(Token::Chimera), Token::Number) => {
                let start_num_second = chars2num(grps.into_iter())?;
                ranges.push(MonomerNumber::Chimera(start_num.unwrap(), start_num_second));
                start_token.take();
                start_num.take();
                prev_token.take();
            }
            (Some(Token::Number), Some(_), Token::Other(_)) => {
                bail!("Invalid pattern. {s:?}")
            }
            (Some(_), Some(_), _) => {
                unreachable!("{start_token:?},{prev_token:?},{token:?}")
            }
        }
    }
    Ok(ranges)
}

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
        let monomers = extract_monomer_order(mons)?;
        let monomer_base = Monomer::new(&format!("{mon_info}.1"))?;
        let fn_get_new_mon = |m| {
            let mut final_mon = monomer_base.clone();
            final_mon.monomer_1 = m;
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
                        new_monomers.extend(new_range.rev().map(fn_get_new_mon).rev())
                    } else {
                        new_monomers.extend(range.clone().map(fn_get_new_mon))
                    }
                }
                MonomerNumber::Single(m) => {
                    let mut final_mon = monomer_base.clone();
                    final_mon.monomer_1 = *m;
                    new_monomers.push(final_mon);
                }
                MonomerNumber::Chimera(m1, m2) => {
                    let mut final_mon = monomer_base.clone();
                    final_mon.monomer_1 = *m1;
                    final_mon.monomer_2 = Some(*m2);
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
    /// let hor = HOR::new("S01/1C3H1L.11-6").unwrap();
    /// ```
    pub fn new(s: &str) -> eyre::Result<Self> {
        HOR::from_str(s)
    }

    /// Generate the reversed version of this [`HOR`].
    ///
    /// ```
    /// let hor = HOR::new("S01/1C3H1L.11-6").unwrap();
    /// let rev_hor = hor.reversed();
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
                MonomerNumber::Chimera(m1, m2) => MonomerNumber::Chimera(*m2, *m1),
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
                if let Some(m2) = mon.monomer_2.as_mut() {
                    std::mem::swap(m2, &mut mon.monomer_1);
                }
                mon
            })
            .collect_vec();
        Self {
            monomer_structure: new_monomer_structure,
            monomers: new_monomers,
        }
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
                MonomerNumber::Chimera(m1, m2) => {
                    write!(f, "{m1}/{m2}")?;
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

    use crate::hor::HOR;

    #[test]
    fn test_valid_stv() {
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
        assert_eq!(format!("{rev_res}"), "S4CH1L.1-15_26-31_32/31_34/32_35-46");
    }
}
