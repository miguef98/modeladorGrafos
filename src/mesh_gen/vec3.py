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
        maxCoord = self.maxCoord()
        if np.isclose( self.x , maxCoord ):
            return 'x'
        elif np.isclose(self.y, maxCoord):
            return 'y'
        else:
            return 'z'
    
    def argmaxAbsCoord( self ):
        maxCoord = self.maxCoord()
        if np.isclose( np.abs(self.x) , maxCoord ):
            return 'x'
        elif np.isclose( np.abs(self.y) , maxCoord ):
            return 'y'
        else:
            return 'z'

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


class Interpolada:
    def __init__( self, puntos ):
        '''
            Supongo puntos es lista de Vec3
        '''
        self.puntos = puntos
    
    def evaluar( self, t ):
        '''
            Recibe t entre [0,1).
        '''
        if t < 0 or t > 1:
            raise ValueError("Se espera t en rango [0,1)")

        cantCurvas = (len(self.puntos) - 3)
        
        if np.isclose(t, 1):
          return self.evaluarCurva( cantCurvas, 1 )

        indicePunto = np.floor( t * cantCurvas ).astype(np.uint32) + 1
        
        return self.evaluarCurva(indicePunto, ( t - (indicePunto - 1) / cantCurvas ) * cantCurvas )

    def evaluarCurva( self, indice, t ):
        if t < 0 or t > 1:
            raise ValueError("Se espera t en rango [0,1)")

        def spline_4p( t, p_1, p0, p1, p2 ):

            return (
                t*((2-t)*t - 1)   * p_1
                + (t*t*(3*t - 5) + 2) * p0
                + t*((4 - 3*t)*t + 1) * p1
                + (t-1)*t*t         * p2 ) / 2

        return spline_4p(t, self.puntos[indice - 1], self.puntos[indice], self.puntos[indice + 1], self.puntos[indice + 2])

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




