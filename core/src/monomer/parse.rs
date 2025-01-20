use eyre::{bail, ContextCompat};
use std::str::FromStr;

use itertools::Itertools;

use super::{
    chrom::Chromosome, mon_type::MonomerType, sf::SF, status::Status, token::Token, Monomer,
};

impl FromStr for Monomer {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut monomers: Vec<u8> = vec![];
        let mut suprachromosomal_family: Vec<SF> = Vec::with_capacity(2);
        let mut chromosomes: Vec<Chromosome> = vec![];
        let mut monomer_type: Option<MonomerType> = None;
        let mut monomer_type_desc: Option<String> = None;
        let mut status: Option<Status> = None;

        // TODO: Bookkeeping for position.
        // Create peekable iterator.
        let tokens = &s.chars().chunk_by(|c| Token::from(*c));
        let mut tokens_iter = tokens.into_iter().peekable();

        while let Some((token, values)) = tokens_iter.next() {
            match token {
                Token::SF => {
                    while let Some((tk, sf_values)) =
                        tokens_iter.next_if(|(tk, _)| *tk == Token::Number || *tk == Token::Chimera)
                    {
                        // Skip / in 1/01
                        if tk == Token::Chimera {
                            continue;
                        }
                        let sf_str = String::from_iter(sf_values);
                        suprachromosomal_family.push(SF::from_str(&sf_str)?);
                    }
                }
                Token::Chrom => {
                    while let Some((tk, chr_values)) = tokens_iter.next_if(|(tk, _)| {
                        let tk_is_num = *tk == Token::Number;
                        let tk_is_alpha = std::mem::discriminant(tk)
                            == std::mem::discriminant(&Token::Value('a'));
                        let tk_is_delim = *tk == Token::Chimera;
                        tk_is_alpha || tk_is_num || tk_is_delim
                    }) {
                        // Skip / in cases like 1/5/19
                        if tk == Token::Chimera {
                            continue;
                        }
                        let chrom_str = String::from_iter(chr_values);
                        chromosomes.push(Chromosome::from_str(&chrom_str)?);
                    }
                }
                Token::Monomer => {
                    let Some((_, mon_values)) = tokens_iter.next_if(|(tk, _)| *tk == Token::Number)
                    else {
                        bail!("No numeric value after '.'");
                    };
                    let mon1_str = String::from_iter(mon_values);
                    monomers.push(u8::from_str(&mon1_str)?);

                    let Some(_) = tokens_iter.next_if(|(tk, _)| *tk == Token::Chimera) else {
                        continue;
                    };
                    let Some((_, mon_2_values)) =
                        tokens_iter.next_if(|(tk, _)| *tk == Token::Number)
                    else {
                        bail!(
                            "Unexpected token, {:?}, after chimeric monomer delimiter.",
                            tokens_iter
                                .next()
                                .map(|(_, tk_vals)| tk_vals.into_iter().join(","))
                        );
                    };
                    let mon2_str = String::from_iter(mon_2_values);
                    monomers.push(u8::from_str(&mon2_str)?);
                }
                Token::Live | Token::Divergent => {
                    status = Some(Status::try_from(token.char())?);
                }
                Token::MType => {
                    // Take num
                    let Some((_, mtype_vals)) = tokens_iter.next_if(|(tk, _)| *tk == Token::Number)
                    else {
                        bail!(
                            "Unexpected token, {:?}, after chimeric monomer delimiter.",
                            tokens_iter
                                .next()
                                .map(|(_, tk_vals)| tk_vals.into_iter().join(","))
                        )
                    };
                    let mut mtype = String::from_iter(mtype_vals);
                    mtype.insert(0, 'H');
                    monomer_type = Some(MonomerType::from_str(&mtype)?);

                    // Hyphen found. Is commented.
                    let Some(_) = tokens_iter.next_if(|(tk, _)| *tk == Token::Hyphen) else {
                        continue;
                    };
                    // For mtype comment.
                    // Consume values until non-value.
                    let mut monomer_desc = String::new();
                    while let Some((desc_token, _)) = tokens_iter.next_if(|(tk, _)| {
                        std::mem::discriminant(tk) == std::mem::discriminant(&Token::Value('a')) ||
                        // Edge case since C can be chrom or a comment.
                        tk.char() == 'C'
                    }) {
                        monomer_desc.push(desc_token.char());
                    }
                    monomer_type_desc = (!monomer_desc.is_empty()).then_some(monomer_desc);
                }
                Token::Hyphen | Token::Number | Token::Chimera => {
                    bail!(
                        "Invalid monomer str, {s}. Unconsumed token, {}.",
                        values.into_iter().join(",")
                    )
                }
                Token::Value(v) => {
                    bail!("Invalid monomer str, {s}. Unknown character, {v}.")
                }
            }
        }

        Ok(Monomer {
            monomers,
            suprachromosomal_family,
            chromosomes,
            monomer_type: monomer_type
                .with_context(|| format!("Invalid monomer, {s}. Monomer type is required."))?,
            monomer_type_desc,
            status,
        })
    }
}

#[cfg(test)]
mod test {
    use crate::monomer::{chrom::Chromosome, sf::SF, Monomer, MonomerType, Status};

    #[test]
    fn test_invalid_mon() {
        const MON_NO_REQ_ATTR: &str = "S1";
        // Doesn't start with an attribute
        const MON_NO_START_ATTR: &str = "1C16H1.2";

        assert!(Monomer::new(MON_NO_REQ_ATTR).is_err());
        assert!(Monomer::new(MON_NO_START_ATTR).is_err());
    }

    #[test]
    fn test_nonnumber_chrom_mon() {
        const MON: &str = "S4CYH1L.46";
        assert_eq!(
            Monomer {
                monomers: vec![46],
                suprachromosomal_family: vec![SF::SF4],
                chromosomes: vec![Chromosome::CY],
                monomer_type: MonomerType::H1,
                monomer_type_desc: None,
                status: Some(Status::Live),
            },
            Monomer::new(MON).unwrap()
        )
    }

    #[test]
    fn test_live_mon() {
        const MON_LIVE: &str = "S1C16H1L.2";
        assert_eq!(
            Monomer {
                monomers: vec![2],
                suprachromosomal_family: vec![SF::SF1],
                chromosomes: vec![Chromosome::C16],
                monomer_type: MonomerType::H1,
                monomer_type_desc: None,
                status: Some(Status::Live),
            },
            Monomer::new(MON_LIVE).unwrap()
        )
    }

    #[test]
    fn test_nonlive_mon() {
        const MON_NON_LIVE: &str = "S4C20H7.11";
        assert_eq!(
            Monomer {
                monomers: vec![11],
                suprachromosomal_family: vec![SF::SF4],
                chromosomes: vec![Chromosome::C20],
                monomer_type: MonomerType::H7,
                monomer_type_desc: None,
                status: None,
            },
            Monomer::new(MON_NON_LIVE).unwrap()
        )
    }

    #[test]
    fn test_div_mon() {
        const MON_DIV: &str = "S5C1H6d.1";
        assert_eq!(
            Monomer {
                monomers: vec![1],
                suprachromosomal_family: vec![SF::SF5],
                chromosomes: vec![Chromosome::C1],
                monomer_type: MonomerType::H6,
                monomer_type_desc: None,
                status: Some(Status::Divergent),
            },
            Monomer::new(MON_DIV).unwrap()
        )
    }

    #[test]
    fn test_mon_chimeric() {
        const MON_CHIMERIC: &str = "S2C2H1L.3/1";
        assert_eq!(
            Monomer {
                monomers: vec![3, 1],
                suprachromosomal_family: vec![SF::SF2],
                chromosomes: vec![Chromosome::C2],
                monomer_type: MonomerType::H1,
                monomer_type_desc: None,
                status: Some(Status::Live),
            },
            Monomer::new(MON_CHIMERIC).unwrap()
        );
    }

    #[test]
    fn test_hyphen_mon_type() {
        const MON_HYPHEN_1: &str = "S3C1H2-B.4";
        const MON_HYPHEN_2: &str = "S2C2H2-C.6";
        assert_eq!(
            Monomer {
                monomers: vec![4],
                suprachromosomal_family: vec![SF::SF3],
                chromosomes: vec![Chromosome::C1],
                monomer_type: MonomerType::H2,
                monomer_type_desc: Some(String::from("B")),
                status: None,
            },
            Monomer::new(MON_HYPHEN_1).unwrap()
        );
        assert_eq!(
            Monomer {
                monomers: vec![6],
                suprachromosomal_family: vec![SF::SF2],
                chromosomes: vec![Chromosome::C2],
                monomer_type: MonomerType::H2,
                monomer_type_desc: Some(String::from("C")),
                status: None,
            },
            Monomer::new(MON_HYPHEN_2).unwrap()
        );
    }

    #[test]
    fn test_ambig_mon() {
        const MON_AMBIG: &str = "S1C1/5/19H1L.6/4";
        assert_eq!(
            Monomer {
                monomers: vec![6, 4],
                suprachromosomal_family: vec![SF::SF1],
                chromosomes: vec![Chromosome::C1, Chromosome::C5, Chromosome::C19],
                monomer_type: MonomerType::H1,
                monomer_type_desc: None,
                status: Some(Status::Live),
            },
            Monomer::new(MON_AMBIG).unwrap()
        );
    }

    #[test]
    fn test_multiple_sfs() {
        const MON_SFS: &str = "S01/1C3H1L.17";
        assert_eq!(
            Monomer {
                monomers: vec![17],
                suprachromosomal_family: vec![SF::SF01, SF::SF1],
                chromosomes: vec![Chromosome::C3],
                monomer_type: MonomerType::H1,
                monomer_type_desc: None,
                status: Some(Status::Live),
            },
            Monomer::new(MON_SFS).unwrap()
        )
    }
}
