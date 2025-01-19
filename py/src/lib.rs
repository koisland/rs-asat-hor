use pyo3::{exceptions::PyValueError, prelude::*};

use rs_asat_hor::Monomer;

#[pyclass(name = "Monomer")]
struct PyMonomer(Monomer);

#[pymethods]
impl PyMonomer {
    #[new]
    fn new(value: &str) -> PyResult<Self> {
        Monomer::new(value)
            .map_err(|err| PyValueError::new_err(err.to_string()))
            .map(PyMonomer)
    }

    #[getter]
    fn monomers(&self) -> PyResult<(u8, Option<u8>)> {
        Ok((self.0.monomer_1, self.0.monomer_2))
    }

    #[getter]
    fn sfs(&self) -> PyResult<Vec<String>> {
        Ok(self
            .0
            .suprachromosomal_family
            .iter()
            .map(|sf| sf.to_string())
            .collect())
    }

    #[getter]
    fn chromosomes(&self) -> PyResult<Vec<String>> {
        Ok(self
            .0
            .chromosomes
            .iter()
            .map(|chr| chr.to_string())
            .collect())
    }

    #[getter]
    fn monomer_type(&self) -> PyResult<String> {
        let mut mon_type = self.0.monomer_type.to_string();
        if let Some(mon_type_desc) = self.0.monomer_type_desc {
            mon_type.push('-');
            mon_type.push(mon_type_desc);
        }
        Ok(mon_type)
    }

    #[getter]
    fn status(&self) -> PyResult<Option<String>> {
        Ok(self.0.status.as_ref().map(|status| status.to_string()))
    }
}

/// A Python module implemented in Rust.
#[pymodule]
fn asat_hor(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyMonomer>()?;
    Ok(())
}
