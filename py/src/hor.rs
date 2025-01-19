use pyo3::{exceptions::PyValueError, prelude::*};

use rs_asat_hor::HOR;

use crate::monomer::PyMonomer;

#[pyclass]
struct PyHORIterator {
    iter: std::vec::IntoIter<PyMonomer>,
}

#[pymethods]
impl PyHORIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<PyMonomer> {
        slf.iter.next()
    }
}

#[pyclass(name = "HOR")]
/// A Python wrapper class for [`HOR`]
pub(crate) struct PyHOR(HOR);

#[pymethods]
impl PyHOR {
    #[new]
    fn new(value: &str) -> PyResult<Self> {
        HOR::new(value)
            .map_err(|err| PyValueError::new_err(err.to_string()))
            .map(PyHOR)
    }

    fn __str__(slf: PyRef<'_, Self>) -> String {
        format!("{}", slf.0)
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyHORIterator {
        // Expensive.
        PyHORIterator {
            iter: slf
                .0
                .iter()
                .map(|monomer| PyMonomer(monomer.clone()))
                .collect::<Vec<_>>()
                .into_iter(),
        }
    }

    fn reversed(slf: PyRef<'_, Self>) -> Self {
        Self(slf.0.reversed())
    }
}
