use std::{fmt::Display, str::FromStr};

use eyre::bail;

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
