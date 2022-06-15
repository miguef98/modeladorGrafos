#ifndef SUBDIVIDE_H
#define SUBDIVIDE_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/subdivision_method_3.h>
#include <pybind11/pybind11.h>
#include <vector>
#include <fstream>

typedef CGAL::Simple_cartesian<double>         Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3>    PolygonMesh;
typedef PolygonMesh::Vertex_index vIndex;
typedef PolygonMesh::Face_index fIndex;

namespace py = pybind11;
using namespace CGAL;
namespace params = CGAL::parameters;

class CMesh{
public:
    CMesh( const py::list& vertices, const py::list& caras ) : mesh() {
        auto indices = std::vector<vIndex>();
        for (auto vertice : vertices){
            py::list vert = vertice.cast<py::list>();
            auto punto = Kernel::Point_3( vert[0].cast<float>(), vert[1].cast<float>(), vert[2].cast<float>() );
            indices.push_back( mesh.add_vertex( punto ) );
        }

        for(py::handle cara : caras){
            py::list caraList = cara.cast<py::list>();
            fIndex f = mesh.add_face( indices[ caraList[0].cast<int>() ], indices[ caraList[1].cast<int>() ], indices[ caraList[2].cast<int>() ], indices[ caraList[3].cast<int>() ]);
        }
    }

    void subdivide( int step = 1 ) {
        Subdivision_method_3::CatmullClark_subdivision(mesh, params::number_of_iterations(step));
    }

    py::list getVertices() const{
        auto res = py::list();
        for( vIndex vi : mesh.vertices() ){
            auto vertice = mesh.point(vi);
            py::list verticeList = py::list();
            verticeList.append( vertice.x() ); verticeList.append( vertice.y() ); verticeList.append( vertice.z() );
            res.append( verticeList );
        }
        return res;
    }

    py::list getCaras() const{
        auto res = py::list();
        for(fIndex fd : mesh.faces()){
            CGAL::Vertex_around_face_circulator<PolygonMesh> vcirc(mesh.halfedge(fd), mesh), done(vcirc);
            auto cara = py::list();
            do{
                cara.append( int(*vcirc) );
                vcirc++;
            }while (vcirc != done);
            res.append(cara);
        }
        return res;
    }

    void exportar( const std::string& path ) const {
        std::ofstream out(path);
        out << mesh;
    }

private:
    PolygonMesh mesh;
};

#endif