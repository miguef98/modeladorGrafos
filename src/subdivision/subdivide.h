#ifndef SUBDIVIDE_H
#define SUBDIVIDE_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/subdivision_method_3.h>

#include <vector>

typedef CGAL::Simple_cartesian<double>         Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3>    PolygonMesh;
typedef PolygonMesh::Vertex_index vIndex;
typedef PolygonMesh::Face_index fIndex;

using namespace std;
using namespace CGAL;
namespace params = CGAL::parameters;

class CMesh{
public:
    CMesh( const vector< vector<double> >& vertices, const vector< vector<int> >& caras ) : mesh() {
        auto indices = vector< vIndex >();
        for( auto vertice = vertices.begin() ; vertice != vertices.end() ; vertice++ ){
            indices.push_back( mesh.add_vertex( Kernel::Point_3( (*vertice)[0], (*vertice)[1], (*vertice)[2] ) ) );
        }
        for( auto cara = caras.begin() ; cara != caras.end() ; cara++ ){
            fIndex f = mesh.add_face( indices[(*cara)[0]], indices[(*cara)[1]], indices[(*cara)[2]], indices[(*cara)[3]]);
        }
    }

    void subdivide( int step = 1 ) {
        Subdivision_method_3::CatmullClark_subdivision(mesh, params::number_of_iterations(step));
    }

    vector< vector<double> > getVertices() const{
        auto res = vector< vector<double> >();
        for( vIndex vi : mesh.vertices() ){
            auto vertice = mesh.point(vi);
            res.push_back( { vertice.x(), vertice.y(), vertice.z() } );
        }

        return res;
    }

    vector< vector<int> > getCaras() const{
        auto res = vector< vector<int> >();
        for(fIndex fd : mesh.faces()){
            CGAL::Vertex_around_face_circulator<PolygonMesh> vcirc(mesh.halfedge(fd), mesh), done(vcirc);
            auto cara = vector<int>();
            do{
                cara.push_back( int(*vcirc) );
                vcirc++;
            }while (vcirc != done);
            res.push_back(cara);
        }
        return res;
    }

private:
    PolygonMesh mesh;
};

#endif