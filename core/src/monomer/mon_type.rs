use std::{fmt::Display, str::FromStr};

use eyre::bail;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum MonomerHOR {
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

impl FromStr for MonomerHOR {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "H1" => MonomerHOR::H1,
            "H2" => MonomerHOR::H2,
            "H3" => MonomerHOR::H3,
            "H4" => MonomerHOR::H4,
            "H5" => MonomerHOR::H5,
            "H6" => MonomerHOR::H6,
            "H7" => MonomerHOR::H7,
            "H8" => MonomerHOR::H8,
            "H9" => MonomerHOR::H9,
            _ => bail!("Unknown monomer type, {s}."),
        })
    }
}

impl Display for MonomerHOR {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                MonomerHOR::H1 => "H1",
                MonomerHOR::H2 => "H2",
                MonomerHOR::H3 => "H3",
                MonomerHOR::H4 => "H4",
                MonomerHOR::H5 => "H5",
                MonomerHOR::H6 => "H6",
                MonomerHOR::H7 => "H7",
                MonomerHOR::H8 => "H8",
                MonomerHOR::H9 => "H9",
            }
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AncestralMonomer {
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

impl FromStr for AncestralMonomer {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "W1" => AncestralMonomer::W1,
            "Ca" => AncestralMonomer::Ca,
            "La" => AncestralMonomer::La,
            "Ba" => AncestralMonomer::Ba,
            "Ja" => AncestralMonomer::Ja,
            "Na" => AncestralMonomer::Na,
            "Fa" => AncestralMonomer::Fa,
            "Oa" => AncestralMonomer::Oa,
            "J1" => AncestralMonomer::J1,
            "R1" => AncestralMonomer::R1,
            "W3" => AncestralMonomer::W3,
            "Aa" => AncestralMonomer::Aa,
            "M1" => AncestralMonomer::M1,
            "R2" => AncestralMonomer::R2,
            "Ea" => AncestralMonomer::Ea,
            "Ia" => AncestralMonomer::Ia,
            "W5" => AncestralMonomer::W5,
            "Qa" => AncestralMonomer::Qa,
            "Ga" => AncestralMonomer::Ga,
            "Ta" => AncestralMonomer::Ta,
            "D2" => AncestralMonomer::D2,
            "W2" => AncestralMonomer::W2,
            "D1" => AncestralMonomer::D1,
            "Ka" => AncestralMonomer::Ka,
            "Ha" => AncestralMonomer::Ha,
            "Pa" => AncestralMonomer::Pa,
            "J2" => AncestralMonomer::J2,
            "W4" => AncestralMonomer::W4,
            _ => bail!("Unknown monomer type, {s}."),
        })
    }
}
