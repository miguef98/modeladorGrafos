#include <pybind11/pybind11.h>
#include "subdivide.h"

namespace py = pybind11;

PYBIND11_MODULE(MeshGen, m) {
    py::class_<CMesh>(m, "CMesh")
        .def(py::init<const std::vector< std::vector<float> > &, const vector< vector<unsigned int> >& >())
        .def("subdivide", &CMesh::subdivide, py::arg("step") = 1)
        .def("getVertices", &CMesh::getVertices)
        .def("getCaras", &CMesh::getCaras);
}