mod chrom;
mod hor;
mod monomer;
mod sf;

use hor::HOR;

fn main() -> eyre::Result<()> {
    const HOR: &str = "S2C16H2-A.4_7-8";
    HOR::new(HOR)?;

    Ok(())
}
