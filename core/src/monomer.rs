use std::fmt::Display;
use std::str::FromStr;

use eyre::{bail, ContextCompat};
use itertools::Itertools;

use crate::chrom::Chromosome;
use crate::sf::SF;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Status {
    Live,
    Divergent,
}

impl FromStr for Status {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "L" | "live" => Status::Live,
            "d" | "divergent" => Status::Divergent,
            _ => bail!("Invalid status, {s}."),
        })
    }
}
impl Display for Status {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Status::Live => write!(f, "Live"),
            Status::Divergent => write!(f, "Divergent"),
        }
    }
}

impl TryFrom<char> for Status {
    type Error = eyre::Error;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        Ok(match value {
            'L' => Status::Live,
            'd' => Status::Divergent,
            _ => bail!("Invalid status, {value}."),
        })
    }
}

/// An alpha-satellite higher-order repeat monomer.
///
/// ```
/// use rs_asat_hor::Monomer;
///
/// let mon = Monomer::new("S1C16H1L.2");
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Monomer {
    pub monomer_1: u8,
    pub monomer_2: Option<u8>,
    pub suprachromosomal_family: Vec<SF>,
    pub chromosomes: Vec<Chromosome>,
    pub monomer_type: MonomerType,
    pub monomer_type_desc: Option<char>,
    pub status: Option<Status>,
}

impl Monomer {
    pub fn new(s: &str) -> eyre::Result<Self> {
        Monomer::from_str(s)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum Token {
    SF,
    Chrom,
    Monomer,
    Divergent,
    Live,
    MType,
    Number,
    Chimera,
    Hyphen,
    Value(char),
}

impl Token {
    fn char(&self) -> char {
        match self {
            Token::SF => 'S',
            Token::Chrom => 'C',
            Token::Number => 'n',
            Token::Divergent => 'd',
            Token::Live => 'L',
            Token::MType => 'H',
            Token::Monomer => '.',
            Token::Chimera => '/',
            Token::Hyphen => '-',
            Token::Value(value) => *value,
        }
    }
}

impl From<char> for Token {
    fn from(value: char) -> Self {
        match value {
            'S' => Token::SF,
            'C' => Token::Chrom,
            'H' => Token::MType,
            '0'..='9' => Token::Number,
            'd' => Token::Divergent,
            'L' => Token::Live,
            '.' => Token::Monomer,
            '/' => Token::Chimera,
            '-' => Token::Hyphen,
            _ => Token::Value(value),
        }
    }
}

impl FromStr for Monomer {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut monomer_1: Option<u8> = None;
        let mut monomer_2: Option<u8> = None;
        let mut suprachromosomal_family: Vec<SF> = Vec::with_capacity(2);
        let mut chromosomes: Vec<Chromosome> = vec![];
        let mut monomer_type: Option<MonomerType> = None;
        let mut monomer_type_desc: Option<char> = None;
        let mut status: Option<Status> = None;

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
                    monomer_1 = Some(u8::from_str(&mon1_str)?);

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
                    monomer_2 = Some(u8::from_str(&mon2_str)?);
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
                    let Some((desc_token, _)) = tokens_iter.next_if(|(tk, _)| {
                        std::mem::discriminant(tk) == std::mem::discriminant(&Token::Value('a')) ||
                        // Edge case since C can be chrom or a comment.
                        tk.char() == 'C'
                    }) else {
                        bail!(
                            "Unexpected token, {:?}, after monomer type hyphen delimiter.",
                            tokens_iter
                                .next()
                                .map(|(_, tk_vals)| tk_vals.into_iter().join(","))
                        )
                    };
                    monomer_type_desc = Some(desc_token.char())
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
            monomer_1: monomer_1.with_context(|| {
                format!("Invalid monomer, {s}. At least one monomer is required.")
            })?,
            monomer_2,
            suprachromosomal_family,
            chromosomes,
            monomer_type: monomer_type
                .with_context(|| format!("Invalid monomer, {s}. Monomer type is required."))?,
            monomer_type_desc,
            status,
        })
    }
}

impl Display for Monomer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let status = match self.status {
            Some(Status::Live) => "L",
            Some(Status::Divergent) => "d",
            None => "",
        };
        let monomer_2 = self
            .monomer_2
            .map(|m2| format!("/{m2}"))
            .unwrap_or_default();
        let monomer_type_desc = self
            .monomer_type_desc
            .map(|desc| format!("-{desc}"))
            .unwrap_or_default();
        let chromosomes = self.chromosomes.iter().join("/");
        let sfs = self.suprachromosomal_family.iter().join("/");
        write!(
            f,
            "S{}C{}{:?}{monomer_type_desc}{status}.{}{monomer_2}",
            sfs, chromosomes, self.monomer_type, self.monomer_1
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MonomerType {
    H1,
    H2,
    H3,
    H4,
    H5,
    H6,
    H7,
    H8,
}

impl FromStr for MonomerType {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "H1" => MonomerType::H1,
            "H2" => MonomerType::H2,
            "H3" => MonomerType::H3,
            "H4" => MonomerType::H4,
            "H5" => MonomerType::H5,
            "H6" => MonomerType::H6,
            "H7" => MonomerType::H7,
            "H8" => MonomerType::H8,
            _ => bail!("Unknown monomer type, {s}."),
        })
    }
}

impl Display for MonomerType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                MonomerType::H1 => "H1",
                MonomerType::H2 => "H2",
                MonomerType::H3 => "H3",
                MonomerType::H4 => "H4",
                MonomerType::H5 => "H5",
                MonomerType::H6 => "H6",
                MonomerType::H7 => "H7",
                MonomerType::H8 => "H8",
            }
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AncestralMonomerTypes {
    W1,
    Ca,
    La,
    Ba,
    Ja,
    Na,
    Fa,
    Oa,
    J1,
    R1,
    W3,
    Aa,
    M1,
    R2,
    Ea,
    Ia,
    W5,
    Qa,
    Ga,
    Ta,
    D2,
    W2,
    D1,
    Ka,
    Ha,
    Pa,
    J2,
    W4,
}

impl FromStr for AncestralMonomerTypes {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "W1" => AncestralMonomerTypes::W1,
            "Ca" => AncestralMonomerTypes::Ca,
            "La" => AncestralMonomerTypes::La,
            "Ba" => AncestralMonomerTypes::Ba,
            "Ja" => AncestralMonomerTypes::Ja,
            "Na" => AncestralMonomerTypes::Na,
            "Fa" => AncestralMonomerTypes::Fa,
            "Oa" => AncestralMonomerTypes::Oa,
            "J1" => AncestralMonomerTypes::J1,
            "R1" => AncestralMonomerTypes::R1,
            "W3" => AncestralMonomerTypes::W3,
            "Aa" => AncestralMonomerTypes::Aa,
            "M1" => AncestralMonomerTypes::M1,
            "R2" => AncestralMonomerTypes::R2,
            "Ea" => AncestralMonomerTypes::Ea,
            "Ia" => AncestralMonomerTypes::Ia,
            "W5" => AncestralMonomerTypes::W5,
            "Qa" => AncestralMonomerTypes::Qa,
            "Ga" => AncestralMonomerTypes::Ga,
            "Ta" => AncestralMonomerTypes::Ta,
            "D2" => AncestralMonomerTypes::D2,
            "W2" => AncestralMonomerTypes::W2,
            "D1" => AncestralMonomerTypes::D1,
            "Ka" => AncestralMonomerTypes::Ka,
            "Ha" => AncestralMonomerTypes::Ha,
            "Pa" => AncestralMonomerTypes::Pa,
            "J2" => AncestralMonomerTypes::J2,
            "W4" => AncestralMonomerTypes::W4,
            _ => bail!("Unknown monomer type, {s}."),
        })
    }
}

#[cfg(test)]
mod test {
    use crate::{
        chrom::Chromosome,
        monomer::{Monomer, MonomerType, Status},
        sf::SF,
    };

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
                monomer_1: 46,
                monomer_2: None,
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
                monomer_1: 2,
                monomer_2: None,
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
                monomer_1: 11,
                monomer_2: None,
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
                monomer_1: 1,
                monomer_2: None,
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
                monomer_1: 3,
                monomer_2: Some(1),
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
                monomer_1: 4,
                monomer_2: None,
                suprachromosomal_family: vec![SF::SF3],
                chromosomes: vec![Chromosome::C1],
                monomer_type: MonomerType::H2,
                monomer_type_desc: Some('B'),
                status: None,
            },
            Monomer::new(MON_HYPHEN_1).unwrap()
        );
        assert_eq!(
            Monomer {
                monomer_1: 6,
                monomer_2: None,
                suprachromosomal_family: vec![SF::SF2],
                chromosomes: vec![Chromosome::C2],
                monomer_type: MonomerType::H2,
                monomer_type_desc: Some('C'),
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
                monomer_1: 6,
                monomer_2: Some(4,),
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
                monomer_1: 17,
                monomer_2: None,
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
