#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "subdivide.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector< std::vector< double > >);
PYBIND11_MAKE_OPAQUE(std::vector< std::vector< int > >);

PYBIND11_MODULE(MeshGen, m) {
    py::class_<CMesh>(m, "CMesh")
        .def(py::init< const std::vector< std::vector<double> > &, const vector< vector<int> >& > ())
        .def("subdivide", &CMesh::subdivide, py::arg("step") = 1, py::return_value_policy::reference )
        .def("getVertices", &CMesh::getVertices)
        .def("getCaras", &CMesh::getCaras);
}