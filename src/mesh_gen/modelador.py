from src.mesh_gen.mesh import MeshGrafo
import numpy as np
import networkx as nx
from scipy.optimize import minimize
from scipy import integrate
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
    
    def crearArista( self, nodoOrigen, nodoFin ):
        self.G.add_edge( nodoOrigen, nodoFin )
        self.setearAristaNoProcesada( nodoOrigen, nodoFin )
    
    def eliminarNodo( self, nodo ):
        self.G.remove_node( nodo )

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

    def grafoDeRamas( self ):
        grafo = nx.Graph()

        for nodo in self.nodos():
            if self.gradoNodo( nodo ) == 1 or self.gradoNodo( nodo ) > 2:
                grafo.add_node( nodo )

        for nodo in grafo.nodes:
            ramas = self.obtenerRamasDesdeNodo( nodo )
            for rama in ramas:
                if not (nodo, rama[-1]) in grafo.edges:
                    grafo.add_edge(nodo, rama[-1], rama=rama )

        return grafo

    def resamplear( self, alpha=0.1, beta=0.1, w=0.01 ):
        grafoRamas = self.grafoDeRamas( )
        diccionarioRamas = nx.get_edge_attributes(grafoRamas, 'rama')

        for edge in grafoRamas.edges:
            self.resamplearRama( diccionarioRamas[edge], alpha, beta, w )

        self.G = nx.convert_node_labels_to_integers( self.G )

    def resamplearRama( self, listaNodos, alpha, beta, w ):
        
        posicionesNodos = [ self.posicionNodo( nodo ) for nodo in listaNodos ]
        curvaPosicionesInterpolada = Interpolada(  posicionesNodos ).reparametrizar( lambda x : np.clip(x, 0, 1))

        radioNodos = [ self.radioNodo( nodo ) for nodo in listaNodos ]
        radiosInterpolados = Interpolada( radioNodos ).reparametrizar( lambda x : np.clip(x, 0, 1))
        
        curvaturaNodos = curvatura( posicionesNodos )
        curvaturaInterpolada = Interpolada( curvaturaNodos  ).reparametrizar( lambda x : np.clip(x, 0, 1))

        termino = lambda j : radiosInterpolados[j] / (1 + beta * curvaturaInterpolada[j])
        
        def costo( ts ):
            
            sigmoide = lambda x : 1 / (1 + np.exp(-10*x))
            f = lambda x, y: np.exp( -x + y )
            h = lambda x, y: -x + f(0,y)
            g = lambda x, y: h(x,y) * sigmoide(-x) + f(x,y) * sigmoide(x)

            def CostoPrincipal( xs ):
                return np.sum( 
                    [ g( xs[0], 1 / len(xs) )]+
                    [ g( xs[i+1] - xs[i], 1 / len(xs) ) for i in range(0, len(xs)-1)] +
                    [ g( 1 - xs[-1], 1 / len(xs) )] )

            def CostoSecundario( xs ):
                return np.sum( 
                    [ g( xs[0], alpha * (termino(0) + termino(xs[0])) )]+
                    [ g( xs[i+1] - xs[i], alpha* (termino(xs[i]) + termino(xs[i+1])) ) for i in range(0, len(xs)-1)] +
                    [ g( 1 - xs[-1], alpha*(termino(xs[1]) + termino(1) ) )] )

            
            return CostoPrincipal( ts ) + w * CostoSecundario( ts )
        
        cantPuntos = np.min( [self.estimadorCantPuntos( termino, alpha), int(curvaPosicionesInterpolada.longitudDeArco() // (0.5 * np.max(radioNodos) ))] )
        paso = 1 / cantPuntos
        ts = np.linspace(0 + paso, 1 - paso, cantPuntos)
        parametros = minimize( costo, ts )        

        self.actualizarRama( listaNodos, curvaPosicionesInterpolada.evaluarLista(parametros.x), radiosInterpolados.evaluarLista(parametros.x) )

    @staticmethod
    def curvaInterpoladaConBordes( puntos, bordeIzq, bordeDer, cantPuntos ):
        radioNodos = np.concatenate( [ bordeIzq, puntos, bordeDer  ] )
        primerIndice = ( 1 / len(radioNodos) ) * cantPuntos
        ultimoIndice = ( 1 / len(radioNodos) ) * ( cantPuntos + len(puntos) )
        return Interpolada( radioNodos ).reparametrizar( lambda x : (ultimoIndice - primerIndice) * x + primerIndice ), ( -primerIndice / (ultimoIndice - primerIndice) + 0.01, (1-primerIndice) / (ultimoIndice - primerIndice) - 0.01)
    
    @staticmethod
    def estimadorCantPuntos( h, alpha, grado=5 ):
        integral = integrate.quad(h, 0, 1)
        ak = integrate.newton_cotes( grado )[0]
        return np.max( [2, int( ((1 + alpha * (h(0) - h(1))) * np.max(ak)) / (2 * alpha * integral[0]) )])

    def obtenerCurvas( self ):
        return self.obtenerCurvasGrafo( self.elegirNodoGrado( 1 ), None, {} )

    def obtenerCurvasGrafo( self, nodoInicial, nodoProcedente, nodosJointVisitados ):
        ramas = self.obtenerRamasDesdeNodo( nodoInicial, nodoProcedente )
    
        curvasPosiciones = []
        curvasRadios = []
        curvasCurvaturas = []

        for rama in ramas:
            radioNodos = [ self.radioNodo( nodo ) for nodo in rama ]
            bordeIzq = np.linspace( np.power(10, 1), radioNodos[0], len(rama) )
            bordeDer = np.linspace( radioNodos[-1], np.power(10, 1), len(rama) )

            radiosInterpolados, limitesRadios = self.curvaInterpoladaConBordes( radioNodos, bordeIzq, bordeDer, len(rama) )

            posicionesNodos = [ self.posicionNodo( nodo ) for nodo in rama ]
            curvaturaNodos = curvatura( posicionesNodos )
            bordeIzq = list(reversed( np.exp(np.linspace(0, -10, len(rama)) * 1)* curvaturaNodos[0] + 0.5 ) )
            bordeDer = np.exp( np.linspace(0, -10, len(rama)) * 1 ) * curvaturaNodos[-1] + 0.5
            curvaturaInterpolada, limitesCurvatura = self.curvaInterpoladaConBordes( curvaturaNodos, bordeIzq, bordeDer, len(rama) )
            
            posicionesNodos = [ self.posicionNodo( nodo ) for nodo in rama ]
            curvaPosicionesInterpolada = Interpolada(  [posicionesNodos[0]] + posicionesNodos + [posicionesNodos[-1]] )
            
            curvasPosiciones.append(curvaPosicionesInterpolada)
            curvasRadios.append(radiosInterpolados)
            curvasCurvaturas.append(curvaturaInterpolada)
            

            ultimoNodoNuevo = rama[-2]
            
            proxNodoInicial = rama[-1]

            if not proxNodoInicial in nodosJointVisitados and self.gradoNodo(proxNodoInicial) != 1:
                nodosJointVisitados[proxNodoInicial] = True
                cp, cr, cc = self.obtenerCurvasGrafo( proxNodoInicial, ultimoNodoNuevo, nodosJointVisitados )

                curvasPosiciones += cp
                curvasRadios += cr
                curvasCurvaturas += cc

        return curvasPosiciones, curvasRadios, curvasCurvaturas


    def actualizarRama( self, nodosARemplazar, nuevasPosiciones, nuevosRadios ):
        [ self.eliminarNodo( nodo ) for nodo in nodosARemplazar[1:-1] ] # elimino los nodos menos los de las puntas

        ultimoNodo = nodosARemplazar[0]
        for posicion, radio in zip( nuevasPosiciones, nuevosRadios):
            nodoNuevo = self.crearNodo( posicion, radio )
            self.crearArista( ultimoNodo, nodoNuevo )
            ultimoNodo = nodoNuevo

        self.crearArista( ultimoNodo, nodosARemplazar[-1] )

        return ultimoNodo


def curvatura( puntos ):
    curvaturas = []

    for i in range(1, len(puntos) - 1):
        punto = puntos[i]
        j = i - 1
        puntoAnterior = puntos[j]
        while j > 0 and puntoAnterior.isClose( punto, rtol=1e-3, atol=1e-3 ):
            j -= 1
            puntoAnterior = puntos[j]

        j = i + 1
        puntoSiguiente = puntos[j]
        while j < len(puntos) and puntoSiguiente.isClose( punto, rtol=1e-3, atol=1e-3 ):
            j += 1
            puntoSiguiente = puntos[j]

        p1 = puntoSiguiente - puntos[i]
        p2 = puntoAnterior - puntos[i]

        angulo = p1.angleTo( p2, p1.cross(p2).normalizar() )

        curvaturas.append( 1 / (angulo / (p1.norm2() + p2.norm2())))

    return [0] + curvaturas + [0]