use std::{fmt::Display, str::FromStr};

use eyre::bail;

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
    H9,
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
            "H9" => MonomerType::H9,
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
                MonomerType::H9 => "H9",
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
