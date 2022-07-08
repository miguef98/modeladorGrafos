from src.mesh_gen.mesh import MeshGrafo
import numpy as np
import networkx as nx
from src.mesh_gen.vec3 import Vec3, Interpolada

def centroMasa( x, y ):
    return (x + y) / 2

class GrafoCentros:
    '''
        Clase grafo de centerline
    '''
    def __init__( self, grafo ):
        self.G = nx.convert_node_labels_to_integers( grafo )
        if nx.number_connected_components( grafo ) != 1:
            raise ValueError( "El grafo tiene mas de 1 componente conexa" )

        self.mesh = MeshGrafo( self )

        self.maxNombre = self.cantNodos( ) - 1

    def tile( self ):
        
        # tengo que elegir algun nodo de grado 1 donde comenzar a recorrer
        # el grafo
        nodoInicial = self.elegirNodoGrado( 1 )

        cola = set()
        self.procesarTile( nodoInicial, cola )

    def procesarTile( self, nodo, cola ):
        if self.gradoNodo( nodo ) == 1:
            vecinoAProcesar = self.iesimoVecino(nodo, 0)
            normalCuadrado = self.direccion(nodo, vecinoAProcesar )
            self.mesh.agregarCuadrado( nodo, normalCuadrado, Vec3.random().projectToPlane(normalCuadrado).normalizar() )       
            self.mesh.tileTrivially( nodo, vecinoAProcesar )
            self.procesarTile( vecinoAProcesar, cola )
        
        else:
            vecinosAProcesar = self.vecinosAProcesar(nodo)
            cantVecinosAProcesar = len(vecinosAProcesar)
            if cantVecinosAProcesar == 0:
                self.mesh.agregarTapaANodo( nodo )
                return self.finCamino( cola )

            elif cantVecinosAProcesar == 1:
                vecinoAProcesar = vecinosAProcesar[0]
                if self.gradoNodo( vecinoAProcesar ) == 1:
                    self.mesh.tileTrivially( nodo, vecinoAProcesar )
                    return self.finCamino( cola )
                elif self.gradoNodo(vecinoAProcesar) == 2:
                    self.mesh.tileTrivially( nodo, vecinoAProcesar )
                    self.procesarTile( vecinoAProcesar, cola )
                else:
                    self.generarIntermediosAJoint( nodo, vecinoAProcesar )
                    vecinosFwd, vecinosBwd = self.clasificarVecinosFwdBwd( nodo, vecinoAProcesar )
                    cola = self.mesh.tileJoint( nodo, vecinoAProcesar, vecinosBwd, None, cola )

                    self.procesarTile( vecinoAProcesar, cola )
            
            else: # osea tengo mas de un vecino a procesar... significa que tengo que unirlos forward!
                vecinoFwdMasCercano = self.nodoMasCercano( nodo, vecinosAProcesar )

                vecinosAProcesar.remove( vecinoFwdMasCercano )
                self.pasarVecinosANodo( nodo, vecinoFwdMasCercano, vecinosAProcesar)
                self.procesarTile( nodo, cola )

    def finCamino(self, cola):
        if len(cola) == 0:
            return # listo termine !
        else:
            proximo = cola.pop()
            if self.gradoNodo( proximo ) != 1:
                self.procesarTile( proximo , cola )
            else:
                return

    def generarNodoIntermedio( self, nodoFrom, nodoTo ):
        self.G.remove_edge( nodoFrom, nodoTo )
        # --- ATENCION ---
        # que nombre tiene el nodo?? puedo estar seguro que len(self.G.nodes) va a andar bien?
        # todo depende de como lo arme al grafo... Lo unico que se es que GrafoCentros no va a borrar nodos
        # asi que por esta parte deberia estar bien (aunque sea por ahora)
        nodoIntermedio = len( self.G.nodes )
        posicionNodoIntermedio = (self.posicionNodo(nodoFrom) + self.posicionNodo( nodoTo )) / 2
        radioNodoIntermedio = (self.radioNodo( nodoFrom ) + self.radioNodo(nodoTo)) / 2
        self.G.add_node( nodoIntermedio, posicion=posicionNodoIntermedio, radio=radioNodoIntermedio )

        self.G.add_edges_from( [(nodoFrom, nodoIntermedio), (nodoIntermedio, nodoTo) ] )
        self.setearAristaNoProcesada( nodoFrom, nodoIntermedio )
        self.setearAristaNoProcesada( nodoIntermedio, nodoTo )

    def generarIntermediosAJoint( self, nodoFrom, nodoJoint ):
        '''
            Este metodo no existe en el algoritmo original. La idea es agregar nodos
            intermedios en el caso de que haya dos nodos joints que son vecinos.
            Si no hago esto tengo un problema porque podria pasar que una joint ya tenga alguna
            de las bifurcaciones tileadas como si fuera de grado 2.
        ''' 
        for vecino in list(self.vecinos( nodoJoint )):
            if self.gradoNodo(vecino) > 2 and vecino != nodoFrom:
                self.generarNodoIntermedio( nodoJoint, vecino )

    def planoPromedioJoint( self, nodoFrom, nodoJoint ):
        nIn = self.direccion( nodoFrom, nodoJoint )
        esPositiva = lambda n: 1 if n.dot(nIn) > 0 else 0

        normales = [ self.direccion( nodoJoint, nodoBifurcacion ) for nodoBifurcacion in self.vecinos( nodoJoint ) if nodoBifurcacion != nodoFrom ]
        return np.sum( [ n_i * esPositiva( n_i ) for n_i in normales ] + [ nIn ] ).normalizar()
    
    def clasificarVecinosFwdBwd( self, nodoFrom, nodoJoint ):
        fwd, bwd = [], []
        nAvg = self.planoPromedioJoint( nodoFrom, nodoJoint )
        for vecino in self.vecinos( nodoJoint ):
            if vecino != nodoFrom:
                grupo = fwd if self.direccion( nodoJoint, vecino ).dot( nAvg ) > 0 else bwd
                grupo.append( vecino )
        
        return fwd, bwd

    def pasarVecinosANodo( self, nodoOriginal, nodoActual, vecinos ):
        for vecino in vecinos:
            self.G.remove_edge( nodoOriginal, vecino )
            self.G.add_edge( nodoActual, vecino )
            self.setearAristaNoProcesada( nodoActual, vecino )

    def direccion( self, nodoFrom, nodoTo ):
        return self.posicionNodo( nodoFrom ).dirTo( self.posicionNodo(nodoTo) )

    def setearAristaProcesada( self, nodoFrom, nodoTo ):
        nx.set_edge_attributes(self.G, {(nodoFrom, nodoTo) : {'procesada':True}})
    
    def setearAristaNoProcesada( self, nodoFrom, nodoTo ):
        nx.set_edge_attributes(self.G, {(nodoFrom, nodoTo) : {'procesada':False}})

    def vecinos( self, nodo ):
        return ( list(arista)[1] for arista in self.G.edges(nodo) )
    
    def vecinosDistintos( self, nodo, nodosDist ):
        return [ vecino for vecino in self.vecinos(nodo) if not vecino in nodosDist ]

    def iesimoVecino( self, nodo, i ):
        return list(self.vecinos( nodo ))[i]

    def vecinosAProcesar( self, nodo ):
        return [ vecino for vecino in self.vecinos(nodo) if not self.aristaFueProcesada(nodo, vecino )]

    def aristaFueProcesada(self, nodoFrom, nodoTo):
        return self.G.get_edge_data( nodoFrom, nodoTo )['procesada']

    def nodoMasCercano( self, nodo, listaNodos ):
        return listaNodos[ np.argmin([ self.posicionNodo(nodo).distTo( self.posicionNodo(otro) ) for otro in listaNodos ] ) ] 

    def vecinoMasCercano( self, nodo ):
        return self.nodoMasCercano( nodo, list( self.vecinos(nodo) ))
    
    def cantNodos( self ):
        return len(self.nodos())

    def nodos( self ):
        return self.G.nodes

    def posicionNodo( self, nodo ):
        # cuando arme el nodo pongo la posicion como un Vec3
        posicion = nx.get_node_attributes( self.G, 'posicion' )[nodo]
        if not isinstance(posicion, Vec3):
            nx.set_node_attributes( self.G, {nodo: Vec3(*posicion)}, 'posicion')
            posicion = nx.get_node_attributes( self.G, 'posicion' )[nodo]

        return posicion

    def radioNodo( self, nodo ):
        return nx.get_node_attributes( self.G, 'radio' )[nodo]

    def gradoNodo( self, nodo ):
        return self.G.degree( nodo )

    def getNuevoNombreNodo( self ):
        self.maxNombre += 1
        return self.maxNombre

    def getVertices( self ):
        return np.array( self.mesh.getVertices() )

    def getCaras( self ):
        return np.array( self.mesh.getCaras() )

    def subdivide( self, step = 1):
        self.mesh.subdivide( step )
        return self

    def elegirNodoGrado( self, grado ):
        for nodo in self.nodos():
            if self.gradoNodo(nodo) == grado:
                return nodo

    def exportar( self, path="result.off" ):
        self.mesh.exportar( path )

    def crearNodo( self, posicion, radio ):
        nombre = self.getNuevoNombreNodo()
        self.G.add_node( nombre, posicion=posicion, radio=radio)
        return nombre
    
    def eliminarNodo( self, nodo ):
        if self.gradoNodo( nodo ) == 2:
            vecinos = list(self.vecinos(nodo))
            self.G.add_edge( *vecinos )
            self.G.remove_node( nodo )
            self.setearAristaProcesada( *vecinos )
        else:
            vecinoMasCercano = self.vecinoMasCercano( nodo )
            for vecino in self.vecinosDistintos(nodo, [vecinoMasCercano]):
                self.G.add_edge( vecinoMasCercano, vecino )
                self.setearAristaNoProcesada(vecinoMasCercano, vecino)  

    def obtenerRamasDesdeNodo( self, nodoInicial, nodoProcedencia=None ):
        '''
            Devuelvo los nodos de una rama, partiendo de un nodo inicial, que presunpongo de grado 1 o n > 2.
        '''
        ramas = []
        nodoPrevio = nodoInicial
        for nodoActual in self.vecinos(nodoPrevio):
            if not nodoProcedencia is None and nodoActual == nodoProcedencia:
                continue
            
            nodosRama = [ nodoPrevio ]

            while self.gradoNodo(nodoActual) == 2:
                nodosRama.append(nodoActual)
                nodoProximo = self.vecinosDistintos( nodoActual, [ nodoPrevio ] )[0]
                nodoPrevio = nodoActual
                nodoActual = nodoProximo

            nodosRama.append(nodoActual)

            ramas.append(nodosRama)
            nodoPrevio = nodoInicial

        return ramas

    def resamplear( self, puntosPorRama ):
        self.resamplearGrafo( self.elegirNodoGrado( 1 ), None, {}, puntosPorRama )
        self.G = nx.convert_node_labels_to_integers( self.G )

    def resamplearGrafo( self, nodoInicial, nodoProcedente, nodosJointVisitados, puntosPorRama ):
        ramas = self.obtenerRamasDesdeNodo( nodoInicial, nodoProcedente )

        for rama in ramas:
            nodosNuevos = self.resamplearRama( rama, puntosPorRama )

            proxNodoInicial = rama[-1]

            if not proxNodoInicial in nodosJointVisitados and self.gradoNodo(proxNodoInicial) != 1:
                nodosJointVisitados[proxNodoInicial] = True
                self.resamplearGrafo( proxNodoInicial, nodosNuevos[-1], nodosJointVisitados, puntosPorRama )

    def resamplearRama( self, listaNodos, N ):
        '''
            Para resamplear una rama (es decir un camino donde salvo el nodo inicial y el final todos son de grado 2),
            interpolo los puntos con una curva usando scipy y luego tomo N-2 puntos de esa curva con un espaciado uniforme.
            No anda para N = 2.
        '''

        # duplico bordes para no perderlos al hacer Catmull-Rom
        posicionesNodos = [ self.posicionNodo( nodo ) for nodo in listaNodos ]
        posicionesNodos = [posicionesNodos[0]] + posicionesNodos + [posicionesNodos[-1]]

        radioNodos = [ self.radioNodo( nodo ) for i, nodo in enumerate(listaNodos) ]
        radioNodos = [ radioNodos[0] ] + radioNodos + [ radioNodos[-1] ]

        curvaPosicionesInterpolada = Interpolada( posicionesNodos )
        curvaRadiosInterpolada = Interpolada( radioNodos )

        nuevasPosiciones = [ curvaPosicionesInterpolada.evaluar( i / (N - 1) ) for i in range(N ) ]
        nuevosRadios = [ curvaRadiosInterpolada.evaluar( i / (N - 1) ) for i in range(N ) ]

        # elimino todos los nodos entre puntas
        for nodo in listaNodos[1:-1]:
            self.G.remove_node( nodo )

        nodosNuevos = []
        ultimoNodo = listaNodos[0]
        for posicion, radio in zip( nuevasPosiciones[1:-1], nuevosRadios[1:-1] ) :
            
            nuevoNodo = self.crearNodo( posicion, radio )
            nodosNuevos.append( nuevoNodo )
            
            self.G.add_edge( ultimoNodo, nuevoNodo )
            self.setearAristaNoProcesada( ultimoNodo, nuevoNodo )

            ultimoNodo = nuevoNodo
        
        self.G.add_edge( ultimoNodo, listaNodos[-1] )
        self.setearAristaNoProcesada( ultimoNodo, listaNodos[-1] )

        return nodosNuevos

    def resamplearRamaBidireccional( self, listaNodos, *, alpha=0.1, beta=0.1 ):
        curvaturas = curvatura( [ self.posicion(nodo) for nodo in listaNodos ] )
        g = lambda i: alpha * self.radioNodo(listaNodos[i]) / (1 * beta * curvaturas[i])

        indice0 = 1
        indice1 = listaNodos[-2]

        while True:
            





def curvatura( puntos ):
    curvaturas = []

    for i in range(1, len(puntos) - 1):
        p1 = puntos[i+1] - puntos[i]
        p2 = puntos[i-1] - puntos[i]

        angulo = p1.angleTo( p2, p1.cross(p2).normalizar() )

        curvaturas.append( angulo / (p1.norm2() + p2.norm2()))

    return curvaturas