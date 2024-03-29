import numpy as np
from dataclasses import dataclass

@dataclass
class Vec3:
    x: float
    y: float
    z: float

    def __add__(self, other ):
        return Vec3( self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__( self, other ):
        return self + (-1 * other)
    
    def __mul__( self, other ):
        return Vec3( self.x * other, self.y * other, self.z * other)

    __rmul__ = __mul__

    def __floordiv__( self, other ):
        return Vec3( self.x // other, self.y // other, self.z // other )
    
    def __truediv__( self, other ):
        return self * (1 / other)

    def __len__( self ):
        return 3
            
    def dot( self, other ):
        return self.x * other.x + self.y * other.y + self.z * other.z 
    
    def cross( self, other ):
        return Vec3(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x
        )
    
    def normalizar( self ):
        norm2 = self.norm2()
        if np.isclose( norm2, 0 ):
            raise ValueError(" El vector nulo no se puede normalizar ")

        self.x /= norm2
        self.y /= norm2
        self.z /= norm2

        return self

    def dirTo( self, other ):
        return (other - self).normalizar()

    def distTo(self, other, norm='2' ):
        if norm == '2':
            return (self - other).norm2()
        elif norm == 'sq2':
            return (self - other).sqNorm2()
        else:
            raise ValueError('Norma invalida')

    def angleTo( self, other, normal ): # asumo normal normalizada
        proy = other.projectToPlane( normal )

        cruz = self.cross(proy)
        cotang = normal.cross( self )

        #uso max coordenada para minimizar error num (se que no es 0 por norma = 1)
        coord = normal.argmaxAbsCoord()
        
        t1 = cruz.getCoord(coord) / (self.norm2() * proy.norm2() * normal.getCoord(coord))
        if np.isclose( t1, 1 ):
            t1 = 1
        elif np.isclose( t1, -1):
            t1 = -1
        elif t1 > 1 or t1 < -1:
            raise Exception( "Algun error numerico..." )

        anguloSelf = np.arcsin( t1 )

        t2 = (cotang.cross(proy)).getCoord(coord) / (cotang.norm2() * proy.norm2() * normal.getCoord(coord))
        if np.isclose( t2, 1 ):
            t2 = 1
        elif np.isclose( t2, -1):
            t2 = -1
        elif t2 > 1 or t2 < -1:
            raise Exception( "Algun error numerico..." )

        anguloCot = np.arcsin( t2 )

        esPositivo = lambda t : t > 0

        if esPositivo( anguloSelf ) and not esPositivo( anguloCot ):
            return anguloSelf
        elif esPositivo( anguloSelf ) and esPositivo( anguloCot ):
            return np.pi / 2 + anguloCot
        elif not esPositivo( anguloSelf ) and not esPositivo( anguloCot ):
            return 2 * np.pi + anguloSelf
        else:
            return 2*np.pi + (-np.pi - anguloSelf)


    def getCoord( self, coord ):
        return getattr( self, coord )

    def argmaxCoord( self ):
        coord = np.argmax( self.toNumpy( ))
        return { 0: 'x', 1: 'y', 2: 'z'}[coord]
    
    def argmaxAbsCoord( self ):
        coord = np.argmax( np.abs(self.toNumpy() ))
        return { 0: 'x', 1: 'y', 2: 'z'}[coord]

    def maxCoord( self ):
        return np.max( [ self.x, self.y, self.z ] )

    def maxAbsCoord( self ):
        return np.max( np.abs( [ self.x, self.y, self.z ] ))

    def norm2( self ):
        return np.sqrt( self.sqNorm2() )

    def sqNorm2( self ):
        return self.dot( self )

    def projectToVector( self, other ):
        return other * ( self.dot(other) / other.sqNorm2() )

    def projectToPlane( self, normal ):
        return self - self.projectToVector( normal )

    def planoFormado( self, v1, v2 ):
        return (( v1 - self ).cross( v2 - self )).normalizar()

    def isClose( self, other, rtol=1e-05, atol=1e-08 ):
        return np.allclose(self.toNumpy(), other.toNumpy() , rtol=rtol, atol=atol)

    def toList( self ):
        return [ self.x, self.y, self.z ]

    def toNumpy( self ):
        return np.array( self.toList() )
    
    def setSize( self, size ):
        return self.normalizar() * size

    @classmethod
    def fromCsv( cls, filaCsv ):
        return cls( filaCsv.loc['X'], filaCsv.loc['Y'], filaCsv.loc['Z'] )

    @classmethod
    def random( cls ):
        # devuelve vector normalizado random
        return Vec3( *np.random.uniform(-1, 1, 3) ).normalizar()

    @classmethod
    def fromList( cls, lista ):
        return [ cls(*i) for i in lista ]


class Interpolada:
    def __init__( self, puntos ):
        if len(puntos) < 2:
            raise ValueError("No se puede interpolar menos de 2 puntos")
        
        self.puntos = [puntos[0]] + puntos + [puntos[-1]]
        self.parametrizacion = lambda x : x

        try:
            self.dimension = len(puntos[0])
        except:
            self.dimension = 1
        
    def __getitem__(self, t):
        return self.evaluar(t)
    
    def evaluar( self, t ):
        '''
            Recibe t entre [0,1].
        '''

        t = self.parametrizacion(t)

        if t < 0 or t > 1:
            raise ValueError("Se espera t en rango [0,1]. Pero se provio t=" + str(t) +".")

        cantCurvas = (len(self.puntos) - 3)
        
        if np.isclose(t, 1):
          return self.evaluarCurva( cantCurvas, 1 )

        indicePunto = np.floor( t * cantCurvas ).astype(np.uint32) + 1
        
        return self.evaluarCurva(indicePunto, ( t - (indicePunto - 1) / cantCurvas ) * cantCurvas )

    def evaluarCurva( self, indice, t ):

        def spline_4p( t, p_1, p0, p1, p2 ):

            return (
                t*((2-t)*t - 1)   * p_1
                + (t*t*(3*t - 5) + 2) * p0
                + t*((4 - 3*t)*t + 1) * p1
                + (t-1)*t*t         * p2 ) / 2

        return spline_4p(t, self.puntos[indice - 1], self.puntos[indice], self.puntos[indice + 1], self.puntos[indice + 2])

    def evaluarLista( self , ts ):
        '''
            ts array-like.
        '''

        return ( self.evaluar(t) for t in ts )

    def longitudDeArco( self, *, eps=0.01, tInicial=0, tFinal=1 ):
        if tFinal - tInicial <= eps:
            return self.evaluar(tInicial).distTo( self.evaluar(tFinal))

        longitud = 0
        ultimoValor = self.evaluar(tInicial)
        for step in np.arange(eps, tFinal + eps, eps):
            nuevoValor = self.evaluar( step )
            longitud += ultimoValor.distTo( nuevoValor )
            ultimoValor = nuevoValor

        return longitud

    def puntosADistancia( self, distancia, *, eps=0.01, tInicial=0, tFinal=1 ):
        '''
            Calculo puntos espaciados por distancia, recorriendo la curva desde tInicial a tFinal con paso epsilon.
        '''

        indices = [ tInicial ]

        tActual = tInicial + eps
        while tActual < tFinal:
            if self.longitudDeArco( tInicial=indices[-1], tFinal=tActual ) >= distancia :
                indices.append(tActual)
            
            tActual += eps

        return indices

    def reparametrizar( self, funcion ):
        self.parametrizacion = funcion
        return self

    def gradiente( self, t0=0, t1=1, eps=0.001, normalizado=False ):
        muestras = np.array( [ muestra.toNumpy() for muestra in self.evaluarLista( np.arange(t0, t1, eps ))])
        ds = np.array([ np.gradient( muestras.T[i], eps ) for i in range(self.dimension) ])

        if self.dimension == 3:
            if normalizado:
                return Interpolada( Vec3.fromList( ds.T ) ), Interpolada( list( map( lambda x : x.normalizar(),  Vec3.fromList( ds.T ) )) )
            else:
                return Interpolada( Vec3.fromList( ds.T ) )
        else:
            return Interpolada( ds.T )

    def curvatura( self ):
        '''
            k = norm( dT / ds ) = norm( dT / dt ) / norm( dC / dt )
        '''

        dC_dt, Tp = self.gradiente(normalizado=True)
        dTp_dt = Tp.gradiente()

        return Interpolada( [ dTp_dt[t].norm2() / dC_dt[t].norm2() for t in np.linspace(0, 1, 100)] )
    



