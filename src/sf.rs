use std::{fmt::Display, str::FromStr};

use eyre::bail;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SF {
    SF01,
    SF02,
    SF1,
    SF2,
    SF3,
    SF4,
    SF5,
}

impl FromStr for SF {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "01" | "SF01" => SF::SF01,
            "02" | "SF02" => SF::SF02,
            "1" | "SF1" => SF::SF1,
            "2" | "SF2" => SF::SF2,
            "3" | "SF3" => SF::SF3,
            "4" | "SF4" => SF::SF4,
            "5" | "SF5" => SF::SF5,
            _ => bail!("Invalid SF class, {s}"),
        })
    }
}

impl Display for SF {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                SF::SF01 => "01",
                SF::SF02 => "02",
                SF::SF1 => "1",
                SF::SF2 => "2",
                SF::SF3 => "3",
                SF::SF4 => "4",
                SF::SF5 => "5",
            }
        )
    }
}
