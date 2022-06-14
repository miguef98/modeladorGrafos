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