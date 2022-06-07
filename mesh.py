import numpy as np
from itertools import permutations

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

    def agregarCuadrado( self, nodo, normal, upVector ):
        self.cuadradoNodo[nodo] = Cuadrado( self.G.posicionNodo(nodo) , normal, upVector )
        return self
    
    def agregarCuadradoOrientadoPorNodo( self, nodoFrom, nodo ):
        ''' 
            Crea cuadrado de nodo suponiendo que se encuentra en un segmento, realizando
            la cuestion de calcular el plano normal sumando direcciones y proyectando el upVector a partir de un
            nodoFrom anterior en el segmento
            Se espera -- deg(nodo) = 2 --
        '''
        if self.G.gradoNodo( nodo ) > 2:
            raise ValueError("Se espera nodo de grado 2 o 1, y el nodo " + str(nodo) + " es de grado " + str(self.G.gradoNodo(nodo)) )

        elif self.G.gradoNodo( nodo ) == 2:
            nodoTo = [ vecino for vecino in self.G.vecinos( nodo ) if vecino != nodoFrom ][0]
            normalNodo = ( self.G.direccion( nodoFrom, nodo ) + self.G.direccion( nodo, nodoTo ) ).normalizar()
        else:
            normalNodo = self.G.direccion( nodoFrom, nodo )

        upVectorNodo = self.cuadradoNodo[nodoFrom].upVector.projectToPlane( normalNodo ).normalizar() * (self.G.radioNodo(nodo) * np.sqrt(2))

        return self.agregarCuadrado( nodo, normalNodo, upVectorNodo )

    def agregarCuadradoOrientadoPorVertices( self, nodoFrom, nodo, indicesVertices ):
        ''' 
            Crea cuadrado de nodo suponiendo que se encuentra en un segmento, realizando
            la cuestion de calcular el plano normal sumando direcciones tomando en cuenta el nodoFrom 
            y proyectando el upVector a partir de un grupo de vertices
            Se espera -- deg(nodo) = 2 --
        '''
        if self.G.gradoNodo( nodo ) > 2:
            raise ValueError("Se espera nodo de grado 2 o 1, y el nodo " + str(nodo) + " es de grado " + str(self.G.gradoNodo(nodo)) )

        elif self.G.gradoNodo( nodo ) == 2:
            nodoTo = [ vecino for vecino in self.G.vecinos( nodo ) if vecino != nodoFrom ][0]
            normalNodo = ( self.G.direccion( nodoFrom, nodo ) + self.G.direccion( nodo, nodoTo ) ).normalizar()
        else:
            normalNodo = self.G.direccion( nodoFrom, nodo )

        centroDeMasaVertices = np.sum( [ self.vertice(i) for i in indicesVertices ] ) / 4
        upVectorDeVertices = self.vertice( indicesVertices[0] ) - centroDeMasaVertices
        upVectorDeNodo = upVectorDeVertices.projectToPlane( normalNodo ).normalizar() * self.G.radioNodo(nodo) * np.sqrt(2)

        self.agregarCuadrado( nodo, normalNodo, upVectorDeNodo )

    def getUpVectorNodo( self, nodo ):
        return self.cuadradoNodo[nodo].upVector

    def getNormalNodo( self, nodo ):
        return self.cuadradoNodo[nodo].normal

    def calcularCaraCuadranteEntreNodos( self, nodoFrom, nodoTo, cuadrante ):
        '''
            Calculo cara del cuadrante indicado por dos nodos
            suponiendo que estan orientados con un upVector proyectado.
        '''
        return [ 
            self.indiceVertice( nodoFrom, cuadrante ),
            self.indiceVertice( nodoTo, cuadrante ),
            self.indiceVertice( nodoTo, (cuadrante + 1) % 4 ),
            self.indiceVertice( nodoFrom, (cuadrante + 1) % 4 )
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
        normasDeFrobenius = [ np.linalg.norm( (np.full_like( permutacionesVertices, verticesNodo) - permutacionesVertices)[i] , 'fro') for i in range(len(permutacionesVertices))]

        verticesOrdenOptimo = permutacionesIndices[ np.argmin( normasDeFrobenius ) ]
        return [ 
            self.indiceVertice( nodo, cuadrante ),
            verticesOrdenOptimo[ cuadrante ],
            self.indiceVertice( nodo, (cuadrante + 1) % 4),
            verticesOrdenOptimo[ (cuadrante + 1) % 4 ]
        ]

    def tileTrivially( self, nodoFrom, nodoTo ):
        if not nodoTo in self.cuadradoNodo:
            self.agregarCuadradoOrientadoPorNodo( nodoFrom, nodoTo )

        # TO DO: else => tengo que ver que hago si ya esta...
        # porque no puedo asumir que tienen el upVector orientado
        # por lo pronto puedo usar agregarCuadranteEntreNodoYVertices
        
        for cuadrante in range(4):
            self.caras.append( self.calcularCaraCuadranteEntreNodos( nodoFrom, nodoTo, cuadrante ) )

        self.G.setearAristaProcesada( nodoFrom, nodoTo )

    def tileJoint( self, nodoFrom, nodoJoint, nodosTo, indicesVertices, cola ):
        if not nodoJoint in self.cuadradoNodo:
            if indicesVertices is None:
                self.agregarCuadradoOrientadoPorNodo( nodoFrom, nodoJoint )
            else:
                self.agregarCuadradoOrientadoPorVertices( nodoFrom, nodoJoint, indicesVertices )

        # si tiene 0 o 1 vecinos a conectar => tengo que agregarlo a la cola
        if len(nodosTo) < 2:            
            cola.append( nodoJoint )

        nodosPorCuadrante = [ [], [], [], [] ]
        [ nodosPorCuadrante[ self.cuadrante( nodoJoint, nodo ) ] for nodo in nodosTo ]
        
        for cuadrante, nodosCuadrante in enumerate(nodosPorCuadrante):
            if len( nodosCuadrante ) == 0:
                if indicesVertices is None:
                    self.caras.append( self.calcularCaraCuadranteEntreNodos( nodoFrom, nodoJoint, cuadrante ) )
                else:
                    self.caras.append( self.calcularCaraCuadranteEntreNodoYVertices( nodoJoint, indicesVertices, cuadrante ) )
            else:
                vecinoMasCercano = self.G.vecinoMasCercano( nodoJoint, nodosCuadrante)

                if indicesVertices is None:
                    indicesVerticesProximaIt = self.calcularCaraCuadranteEntreNodos( nodoFrom, nodoJoint, cuadrante)
                    #normal = vecinoMasCercano.direccion(0)
                else:
                    indicesVerticesProximaIt = self.calcularCaraCuadranteEntreNodoYVertices( nodoJoint, indicesVertices, cuadrante )
                    #normal = vecinoMasCercano.direccion(0) + mesh.vertices[indiceVerticesCuadranteFrom[0]].planoFormado( *np.array(mesh.vertices)[indiceVerticesCuadranteFrom[1:3]] )


                # como mi nuevo joint va a pasar a ser el vecinoMasCercano, tengo que asegurarme que su otro 
                # vecino (se de antemano que vecinoMasCercano es de grado 2) no sea un joint
                vecinoDelMasCercano = self.G.vecinosDistintos( vecinoMasCercano, [ nodoJoint ] )[0]
                if self.G.gradoNodo( vecinoDelMasCercano ) > 2:
                    self.G.generarNodoIntermedio( vecinoMasCercano, vecinoDelMasCercano )

                nodosCuadrante.remove(vecinoMasCercano)
                self.G.pasarVecinosANodo( nodoJoint, vecinoMasCercano, nodosCuadrante )

                self.tileJoint( nodoJoint, vecinoMasCercano, self.G.vecinosDistintos( vecinoMasCercano, [nodoJoint] ), indicesVerticesProximaIt, cola )

        self.G.setearAristaProcesada( nodoFrom, nodoJoint )        

    def cuadrante( self, nodoFrom, nodoTo ):
        '''
            Devuelvo el numero de cuadrante al que pertenece el nodoTo del nodoFrom
        '''
        direccionBifurcacion = self.G.direccion( nodoFrom, nodoTo )
        angulo = self.getUpVectorNodo( nodoFrom ).angleTo( direccionBifurcacion, self.getNormalNodo( nodoFrom ) )
        return [ np.floor(( 4*angulo / (2*np.pi) )).astype(np.uint8) ]

    def vertice( self, indice ):
        '''
            Devuelve el vertice que corresponde al indice.
            Cada nodo tiene 4 vertices.
            *** Supongo que los nodos vienen indexados desde el 0
            y los vertices desde el 1... es decir el nodo 0 tiene los vertices 1,2,3,4 ;
             el nodo 1 los 5,6,7,8 ,etc etc ***
        '''
        return self.cuadradoNodo[ int( (indice - 1) / 4 ) ].vertices[ int( (indice - 1) % 4 ) ]

    @staticmethod
    def indiceVertice( nodo, nroVertice ):
        return ( nodo * 4 ) + ( nroVertice )

    def getVertices( self ):
        return np.array(
            [ self.vertice(i).toNumpy() for i in range(1, len(self.G.nodos()) * 4 + 1)]
        )
    
    def getCaras( self ):
        return np.array( self.caras )

