from mesh import MeshGrafo
import numpy as np
import networkx as nx
from vec3 import Vec3

def centroMasa( x, y ):
    return (x + y) / 2

class GrafoCentros:
    '''
        Clase grafo de centerline
    '''
    def __init__( self, grafo ):
        self.G = grafo
        if nx.number_connected_components( grafo ) != 1:
            raise ValueError( "El grafo tiene mas de 1 componente conexa" )

        self.mesh = MeshGrafo( self )

    def tile( self ):
        
        # tengo que elegir algun nodo de grado 1 donde comenzar a recorrer
        # el grafo
        nodoInicial = None
        for nodo in self.nodos():
            if self.gradoNodo(nodo) == 1:
                nodoInicial = nodo
                break

        cola = []
        self.procesarTile( nodoInicial, cola )
    
    def procesarTile( self, nodo, cola ):

        if self.gradoNodo( nodo ) == 1:
            vecino = self.iesimoVecino( nodo, 0)
            if self.aristaFueProcesada(nodo, vecino):
                return self.finCamino(cola)
            else:
                normalCuadrado = self.direccion(nodo, vecino)
                self.mesh.agregarCuadrado( nodo, normalCuadrado, Vec3.random().projectToPlane(normalCuadrado).normalizar() ) 
                self.mesh.tileTrivially( nodo, vecino )
                self.procesarTile( vecino, cola )
        
        elif self.gradoNodo( nodo ) == 2:
            aristasProcesadas = [ self.aristaFueProcesada(nodo, vecino) for vecino in self.vecinos( nodo ) ]
            if all( aristasProcesadas ):
                return self.finCamino(cola)
            else:
                vecinoAProcesar = self.iesimoVecino( nodo, 0 if aristasProcesadas[1] else 1 )

                if self.gradoNodo( vecinoAProcesar ) > 2:
                    self.generarIntermediosAJoint( vecinoAProcesar )
                    vecinosFwd, vecinosBwd = self.clasificarVecinosFwdBwd( nodo, vecinoAProcesar )
                    self.mesh.tileJoint( nodo, vecinoAProcesar, vecinosBwd, None, cola )
                    
                    if len( vecinosFwd ) > 0:
                        vecinoFwdMasCercano = self.nodoMasCercano( vecinoAProcesar, vecinosFwd )

                        # como mi nuevo joint va a pasar a ser el vecinoMasCercano, tengo que asegurarme que su otro 
                        # vecino (se de antemano que vecinoMasCercano es de grado 2) no sea un joint
                        vecinoDeFwdMasCercano = self.vecinosDistintos( vecinoFwdMasCercano, [ vecinoAProcesar ] )[0]
                        if self.gradoNodo( vecinoDeFwdMasCercano ) > 2:
                            self.generarNodoIntermedio( vecinoFwdMasCercano, vecinoDeFwdMasCercano )

                        vecinosFwd.remove( vecinoFwdMasCercano )
                        self.pasarVecinosANodo( vecinoAProcesar, vecinoFwdMasCercano, vecinosFwd)
                        self.mesh.tileJoint( vecinoAProcesar, vecinoFwdMasCercano, vecinosFwd, None, cola )
                else:
                    self.mesh.tileTrivially( nodo, vecinoAProcesar )
                    self.procesarTile( vecinoAProcesar, cola )
        else:
            raise Exception( "Nunca deberia pasar esto... " )

    def finCamino(self, cola):
        if len(cola) == 0:
            return # listo termine !
        else:
            self.procesarTile( cola.pop() , cola )

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
        nx.set_edge_attributes( self.G, {(nodoFrom, nodoIntermedio):{'procesado': False}, (nodoIntermedio, nodoTo):{'procesado': False}})

    def generarIntermediosAJoint( self, nodoJoint ):
        '''
            Este metodo no existe en el algoritmo original. La idea es agregar nodos
            intermedios en el caso de que haya dos nodos joints que son vecinos.
            Si no hago esto tengo un problema porque podria pasar que una joint ya tenga alguna
            de las bifurcaciones tileadas como si fuera de grado 2.
        ''' 
        for vecino in self.vecinos( nodoJoint ):
            if self.gradoNodo(vecino) > 2:
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
            nx.set_edge_attributes( self.G, {(nodoActual, vecino) : {'procesada' : False}})

    def direccion( self, nodoFrom, nodoTo ):
        return self.posicionNodo( nodoFrom ).dirTo( self.posicionNodo(nodoTo) )
    
    def setearAristaProcesada( self, nodoFrom, nodoTo ):
        nx.set_edge_attributes(self.G, {(nodoFrom, nodoTo) : {'procesada':True}})

    def vecinos( self, nodo ):
        return ( list(arista)[1] for arista in self.G.edges(nodo) )
    
    def vecinosDistintos( self, nodo, nodosDist ):
        return [ vecino for vecino in self.vecinos(nodo) if not vecino in nodosDist ]

    def iesimoVecino( self, nodo, i ):
        return list(self.vecinos( nodo ))[i]

    def aristaFueProcesada(self, nodoFrom, nodoTo):
        return self.G.get_edge_data( nodoFrom, nodoTo )['procesada']

    def nodoMasCercano( self, nodo, listaNodos ):
        return listaNodos[ np.argmin([ self.posicionNodo(nodo).distTo( self.posicionNodo(otro) ) for otro in listaNodos ] ) ] 

    def vecinoMasCercano( self, nodo ):
        return self.nodoMasCercano( nodo, list( self.vecinos(nodo) ))
    
    def nodos( self ):
        return self.G.nodes

    def posicionNodo( self, nodo ):
        # cuando arme el nodo pongo la posicion como un Vec3
        return nx.get_node_attributes( self.G, 'posicion' )[nodo]

    def radioNodo( self, nodo ):
        return nx.get_node_attributes( self.G, 'radio' )[nodo]

    def gradoNodo( self, nodo ):
        return self.G.degree( nodo )

        
def prepararAristas( grafo ):
    for arista in grafo.edges():
        nx.set_edge_attributes( grafo, {arista : {'procesada':False}})

if __name__ == '__main__':
    G = nx.Graph()
    G.add_nodes_from( [(0, {'posicion': Vec3(0,0,0), 'radio': 1.0}), (1, {'posicion': Vec3(1,0,0), 'radio': 1.3}), (2, {'posicion': Vec3(2,0,0), 'radio': 1.5}) ] )
    G.add_edges_from( [(0,1), (1,2)])
    prepararAristas( G )

    GC = GrafoCentros( G )
    GC.tile()