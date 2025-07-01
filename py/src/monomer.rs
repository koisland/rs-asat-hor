use std::ops::Deref;

use pyo3::{exceptions::PyValueError, prelude::*};

use rs_asat_hor::Monomer;

#[pyclass(name = "Monomer")]
/// A Python wrapper class for [`Monomer`]
pub(crate) struct PyMonomer(pub(crate) Monomer);

#[pymethods]
impl PyMonomer {
    #[new]
    fn new(value: &str) -> PyResult<Self> {
        Monomer::new(value)
            .map_err(|err| PyValueError::new_err(err.to_string()))
            .map(PyMonomer)
    }

    fn __str__(slf: PyRef<'_, Self>) -> String {
        format!("{}", slf.0)
    }

    #[getter]
    fn monomers(&self) -> PyResult<Vec<u8>> {
        Ok(self.0.monomers.to_vec())
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
    fn hor(&self) -> PyResult<String> {
        let mut hor = self.0.hor.to_string();
        if let Some(hor_desc) = self.0.hor_desc.as_ref() {
            hor.push('-');
            hor.push_str(hor_desc.deref());
        }
        Ok(hor)
    }

    #[getter]
    fn status(&self) -> PyResult<Option<String>> {
        Ok(self.0.status.as_ref().map(|status| status.to_string()))
    }
}
