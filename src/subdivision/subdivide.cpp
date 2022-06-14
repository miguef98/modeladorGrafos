#include "subdivide.h"

namespace py = pybind11;

PYBIND11_MODULE(MeshGen, m) {
    py::class_<CMesh>(m, "CMesh")
        .def(py::init< const py::list& , const py::list& > ())
        .def("subdivide", &CMesh::subdivide, py::arg("step") = 1, py::return_value_policy::reference )
        .def("getVertices", &CMesh::getVertices)
        .def("getCaras", &CMesh::getCaras);
}