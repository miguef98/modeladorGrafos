import numpy as np
from itertools import permutations
import src.mesh_gen.MeshGen as mg

class Cuadrado:
    ''' 
        Uso esta clase para representar los cuadrados que voy a ir creando a lo largo
        del centerline.
    '''
    def __init__( self, posicion, normal, upVector ):
        self.posicion = posicion
        self.normal = normal
        self.upVector = upVector
        self.sideVector = self.normal.cross( self.upVector )

        self.vertices = [
            self.posicion + self.upVector,
            self.posicion + self.sideVector,
            self.posicion - self.upVector,
            self.posicion - self.sideVector ]

class MeshGrafo:
    '''
        En esta clase almaceno los datos de la malla en si, es decir vertices y caras.
    '''
    def __init__( self, G ):
        self.cuadradoNodo = {}
        self.caras = []
        self.G = G
        self.meshEnC = None

    def agregarCuadrado( self, nodo, normal, upVector ):
        self.cuadradoNodo[nodo] = Cuadrado( self.G.posicionNodo(nodo) , normal, upVector )
        return self
    
    def agregarCuadradoOrientado( self, nodoFrom, nodo, indicesVertices=None ):
        if indicesVertices is None:
            return self.agregarCuadradoOrientadoPorNodo(nodoFrom, nodo )
        else:
            return self.agregarCuadradoOrientadoPorVertices( nodoFrom, nodo, indicesVertices)

    def agregarCuadradoOrientadoPorNodo( self, nodoFrom, nodo ):
        ''' 
            Crea cuadrado de nodo suponiendo que se encuentra en un segmento, realizando
            la cuestion de calcular el plano normal sumando direcciones y proyectando el upVector a partir de un
            nodoFrom anterior en el segmento
        '''
        
        normalNodo = self.calcularNormalNodo( nodoFrom, nodo )
        upVectorNodo = self.getUpVectorNodo(nodoFrom).projectToPlane( normalNodo ).setSize(self.G.radioNodo(nodo) * np.sqrt(2))

        return self.agregarCuadrado( nodo, normalNodo, upVectorNodo )

    def agregarCuadradoOrientadoPorVertices( self, nodoFrom, nodo, indicesVertices ):
        ''' 
            Crea cuadrado de nodo suponiendo que se encuentra en un segmento, realizando
            la cuestion de calcular el plano normal sumando direcciones tomando en cuenta el nodoFrom 
            y proyectando el upVector a partir de un grupo de vertices
        '''

        normalNodo = self.calcularNormalNodo( nodoFrom, nodo )
        centroDeMasaVertices = np.sum( [ self.vertice(i) for i in indicesVertices ] ) / 4
        upVectorDeVertices = self.vertice( indicesVertices[0] ) - centroDeMasaVertices
        upVectorDeNodo = upVectorDeVertices.projectToPlane( normalNodo ).setSize(self.G.radioNodo(nodo) * np.sqrt(2))

        self.agregarCuadrado( nodo, normalNodo, upVectorDeNodo )

    def agregarTapaANodo( self, nodo ):
        self.caras.append( [ self.indiceVertice(nodo, i) for i in range(4) ] )

    def getUpVectorNodo( self, nodo ):
        return self.cuadradoNodo[nodo].upVector

    def getNormalNodo( self, nodo ):
        return self.cuadradoNodo[nodo].normal
    
    def calcularNormalNodo( self, nodoFrom, nodo ):
        if self.G.gradoNodo(nodo) == 1:
            normalNodo = self.G.direccion( nodoFrom, nodo )
        elif self.G.gradoNodo( nodo ) == 2:
            nodoTo = [ vecino for vecino in self.G.vecinos( nodo ) if vecino != nodoFrom ][0]
            normalNodo = ( self.G.direccion( nodoFrom, nodo ) + self.G.direccion( nodo, nodoTo ) ).normalizar()
        else:
            normalNodo = self.G.planoPromedioJoint( nodoFrom, nodo )  
        
        return normalNodo

    def calcularCaraCuadranteEntreNodos( self, nodoFrom, nodoTo, cuadrante ):
        '''
            Calculo cara del cuadrante indicado por dos nodos
            suponiendo que estan orientados con un upVector proyectado.
        '''
        return [ 
            self.indiceVertice( nodoTo, cuadrante ),
            self.indiceVertice( nodoFrom, cuadrante ),
            self.indiceVertice( nodoFrom, (cuadrante + 1) % 4 ),
            self.indiceVertice( nodoTo, (cuadrante + 1) % 4 )
        ]

    def calcularCaraCuadranteEntreNodoYVertices( self, nodo, indicesVertices, cuadrante ):
        '''
            Dado un nodo de vertices (v1,v2,v3,v4) , y 4 vertices (w1,w2,w3,w4) con los que conectarse 
            busco las conexiones (vi, wj) tales que la sumatoria de || vi - wj || sea minima para todas las 
            posibles permutaciones de los w. 
            Luego, calculo los indices que pertenecen a cada cara.
        '''
        verticesNodo = np.array( [ verticeNodo.toNumpy() for verticeNodo in self.cuadradoNodo[nodo].vertices ] )
        verticesAConectar = np.array( [ self.vertice(i).toNumpy() for i in indicesVertices ] )

        permutacionesVertices = [ verticesAConectar[ list(permutacion) ] for permutacion in permutations([0,1,2,3]) ]
        permutacionesIndices = [ np.array( indicesVertices )[ list(permutacion) ] for permutacion in permutations([0,1,2,3]) ]
        normasDeFrobenius = [ np.exp( np.linalg.norm( (np.full_like( permutacionesVertices, verticesNodo) - permutacionesVertices)[i] , 'fro') ) / np.exp( self.G.radioNodo(nodo) * 2) for i in range(len(permutacionesVertices))]

        verticesOrdenOptimo = permutacionesIndices[ np.argmin( normasDeFrobenius ) ]
        indicesVerticesCara = [ 
            self.indiceVertice( nodo, cuadrante ),
            verticesOrdenOptimo[ cuadrante ],
            self.indiceVertice( nodo, (cuadrante + 1) % 4),
            verticesOrdenOptimo[ (cuadrante + 1) % 4 ]
        ]

        centroMasaCara = np.sum( [ self.vertice(i) for i in indicesVerticesCara ] ) / 4
        verticesDesdeCentroMasa = [ self.vertice(i) - centroMasaCara for i in indicesVerticesCara ]
        angulos = [0] + [ verticesDesdeCentroMasa[0].angleTo( vertice, verticesDesdeCentroMasa[0].cross(verticesDesdeCentroMasa[1]).normalizar() ) for vertice in verticesDesdeCentroMasa[1:] ]
        
        ordenAngulos = np.argsort( angulos )

        return [
            indicesVerticesCara[ ordenAngulos[0] ],
            indicesVerticesCara[ ordenAngulos[1] ],
            indicesVerticesCara[ ordenAngulos[2] ],
            indicesVerticesCara[ ordenAngulos[3] ]
        ]




    def tileTrivially( self, nodoFrom, nodoTo ):
        if not nodoTo in self.cuadradoNodo:
            self.agregarCuadradoOrientado( nodoFrom, nodoTo )
            for cuadrante in range(4):
                self.caras.append( self.calcularCaraCuadranteEntreNodos( nodoFrom, nodoTo, cuadrante ) )

        else:
            indicesVerticesNodoTo = [ self.indiceVertice(nodoTo, i) for i in range(4) ]
            for cuadrante in range(4):
                self.caras.append( list(reversed(self.calcularCaraCuadranteEntreNodoYVertices( nodoFrom, indicesVerticesNodoTo, cuadrante))))

        self.G.setearAristaProcesada( nodoFrom, nodoTo )

    def tileJoint( self, nodoFrom, nodoJoint, nodosTo, indicesVertices, cola ):
        if not nodoJoint in self.cuadradoNodo:
            self.agregarCuadradoOrientado( nodoFrom, nodoJoint, indicesVertices )

        cola = cola.union( set(nodosTo) )
        if nodoJoint in cola:
            cola.remove( nodoJoint )

        nodosPorCuadrante = [ [], [], [], [] ]
        [ nodosPorCuadrante[ self.cuadrante( nodoJoint, nodo ) ].append(nodo) for nodo in nodosTo ]
        
        for cuadrante, nodosCuadrante in enumerate(nodosPorCuadrante):
            if len( nodosCuadrante ) == 0:
                if indicesVertices is None:
                    self.caras.append( self.calcularCaraCuadranteEntreNodos( nodoFrom, nodoJoint, cuadrante ) )
                else:
                    self.caras.append( self.calcularCaraCuadranteEntreNodoYVertices( nodoJoint, indicesVertices, cuadrante) )
            else:
                vecinoMasCercano = self.G.nodoMasCercano( nodoJoint, nodosCuadrante)

                if indicesVertices is None:
                    indicesVerticesConexionConMasCercano = self.calcularCaraCuadranteEntreNodos( nodoFrom, nodoJoint, cuadrante)
                else:
                    indicesVerticesConexionConMasCercano = self.calcularCaraCuadranteEntreNodoYVertices( nodoJoint, indicesVertices, cuadrante)

                if not vecinoMasCercano in self.cuadradoNodo:
                    self.agregarCuadradoOrientado( nodoJoint, vecinoMasCercano, indicesVerticesConexionConMasCercano )
                
                nodosCuadrante.remove( vecinoMasCercano )
                cola2 = self.tileJoint( nodoJoint, vecinoMasCercano, nodosCuadrante, indicesVerticesConexionConMasCercano, cola)

                [ self.G.setearAristaProcesada( nodoJoint, i ) for i in nodosCuadrante ]

        self.G.setearAristaProcesada( nodoFrom, nodoJoint )
        return cola

    def cuadrante( self, nodoFrom, nodoTo ):
        '''
            Devuelvo el numero de cuadrante al que pertenece el nodoTo del nodoFrom
        '''
        direccionBifurcacion = self.G.direccion( nodoFrom, nodoTo )
        angulo = self.getUpVectorNodo( nodoFrom ).angleTo( direccionBifurcacion, self.getNormalNodo( nodoFrom ) )
        return np.floor(( 4*angulo / (2*np.pi) )).astype(np.uint8)

    def vertice( self, indice ):
        '''
            Devuelve el vertice que corresponde al indice.
            Cada nodo tiene 4 vertices.
            *** Supongo que los nodos vienen indexados desde el 0
            y los vertices desde el 1... es decir el nodo 0 tiene los vertices 1,2,3,4 ;
             el nodo 1 los 5,6,7,8 ,etc etc ***
        '''
        if indice == 0 or indice > len( self.G.nodos() ) * 4:
            raise ValueError( "El indice " + str(indice) + " esta fuera de rango. El rango posible es (1, " + str(len(self.G.nodos()) * 4) + ")" )
        return self.cuadradoNodo[ int( (indice - 1) / 4 ) ].vertices[ int( (indice - 1) % 4 ) ]

    @staticmethod
    def indiceVertice( nodo, nroVertice ):
        return ( nodo * 4 ) + ( nroVertice + 1 )

    def subdivide( self, step ):
        if self.meshEnC is None:
            self.meshEnC = mg.CMesh( self.getVertices(), self.getCaras() )
        
        self.meshEnC.subdivide( step )

    def getVertices( self ):
        if self.meshEnC is None:
            return [ self.vertice(i).toNumpy() for i in range(1, len(self.G.nodos()) * 4 + 1)]
        else:
            return self.meshEnC.getVertices()
    
    def getCaras( self ):
        if self.meshEnC is None:
            return list( np.array( self.caras ) - np.ones_like( self.caras ) )
        else:
            return self.meshEnC.getCaras()

    def exportar( self, path ):
        if self.meshEnC is None:
            self.meshEnC = mg.CMesh( self.getVertices(), self.getCaras() )
        
        self.meshEnC.exportar( path )

