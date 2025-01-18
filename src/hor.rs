use std::str::FromStr;

use crate::monomer::Monomer;

#[derive(Debug, Clone)]
pub struct HOR {
    monomers: Vec<Monomer>,
}

impl FromStr for HOR {
    type Err = eyre::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        todo!()
    }
}

impl HOR {
    pub fn new(s: &str) -> eyre::Result<Self> {
        HOR::from_str(s)
    }
}
