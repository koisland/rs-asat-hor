use std::{fmt::Display, str::FromStr};

use eyre::bail;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Chromosome {
    C1,
    C2,
    C3,
    C4,
    C5,
    C6,
    C7,
    C8,
    C9,
    C10,
    C11,
    C12,
    C13,
    C14,
    C15,
    C16,
    C17,
    C18,
    C19,
    C20,
    C21,
    C22,
    CX,
    CY,
}

impl Display for Chromosome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Chromosome::C1 => "1",
                Chromosome::C2 => "2",
                Chromosome::C3 => "3",
                Chromosome::C4 => "4",
                Chromosome::C5 => "5",
                Chromosome::C6 => "6",
                Chromosome::C7 => "7",
                Chromosome::C8 => "8",
                Chromosome::C9 => "9",
                Chromosome::C10 => "10",
                Chromosome::C11 => "11",
                Chromosome::C12 => "12",
                Chromosome::C13 => "13",
                Chromosome::C14 => "14",
                Chromosome::C15 => "15",
                Chromosome::C16 => "16",
                Chromosome::C17 => "17",
                Chromosome::C18 => "18",
                Chromosome::C19 => "19",
                Chromosome::C20 => "20",
                Chromosome::C21 => "21",
                Chromosome::C22 => "22",
                Chromosome::CX => "X",
                Chromosome::CY => "Y",
            }
        )
    }
}

impl FromStr for Chromosome {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "1" | "chr1" => Chromosome::C1,
            "2" | "chr2" => Chromosome::C2,
            "3" | "chr3" => Chromosome::C3,
            "4" | "chr4" => Chromosome::C4,
            "5" | "chr5" => Chromosome::C5,
            "6" | "chr6" => Chromosome::C6,
            "7" | "chr7" => Chromosome::C7,
            "8" | "chr8" => Chromosome::C8,
            "9" | "chr9" => Chromosome::C9,
            "10" | "chr10" => Chromosome::C10,
            "11" | "chr11" => Chromosome::C11,
            "12" | "chr12" => Chromosome::C12,
            "13" | "chr13" => Chromosome::C13,
            "14" | "chr14" => Chromosome::C14,
            "15" | "chr15" => Chromosome::C15,
            "16" | "chr16" => Chromosome::C16,
            "17" | "chr17" => Chromosome::C17,
            "18" | "chr18" => Chromosome::C18,
            "19" | "chr19" => Chromosome::C19,
            "20" | "chr20" => Chromosome::C20,
            "21" | "chr21" => Chromosome::C21,
            "22" | "chr22" => Chromosome::C22,
            "X" | "chrX" => Chromosome::CX,
            "Y" | "chrY" => Chromosome::CY,
            _ => bail!("Invalid chromosome, {s}."),
        })
    }
}
