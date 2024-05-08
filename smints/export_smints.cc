#include <pybind11/pybind11.h>
#include <sointegrals>
namespace py = pybind11;


void export_smints(py::module& m) {
    m.def("smints_benbee", &smints_benbee, "Dummy");
}