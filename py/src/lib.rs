pub use pyo3::prelude::*;

mod hor;
mod monomer;

use hor::PyHOR;
use monomer::PyMonomer;

#[pymodule]
fn asat_hor(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMonomer>()?;
    m.add_class::<PyHOR>()?;
    Ok(())
}
