import numpy

luminosity = 200 * 10**9 # 1/barn
sigma_ee_phi = 3*10**-6  # barn
BR_phi_KK = 0.34
BR_KL_2pi = 0.002
BR_KS_2pi = 0.69

epsL = 0.5
epsS = 0.7
eps = epsS*epsL
acceptance = 0.25       #in the fiducial volume

p_K0 = 110.6            #MeV/c
tau_S = 0.8954e-10      # seconds
tau_L = 5.116e-8        # seconds
l_S = 0.6               #cm
l_L = 340               #cm

zlim = 120              #cm
xmin = 35               #cm
ymin = 35               #cm
xmax = 150              #cm
ymax = 150              #cm

m_K0 = 497.611          #MeV
m_pi = 139.57039        #MeV
m_e = 0.51099895000     #MeV
m_mu = 105.6583745      #MeV


def getE(m, px,py,pz):
    return numpy.sqrt(m*m + px*px + py*py + pz*pz)

def sumvec(x1,y1,z1,x2,y2,z2):
    return x1+x2, y1+y2, z1+z2

def dotprod(x1,y1,z1,x2,y2,z2):
    return x1*x2 + y1*y2 + z1*z2

def norm(x,y,z):
    return numpy.sqrt(x*x + y*y + z*z)

def invmass(E, x1, y1, z1):
    sum_squared = dotprod(x1, y1, z1, x1, y1, z1)
    return numpy.sqrt(E**2 - sum_squared)

def getComponents(mag, x, y, z):
    normalization = norm(x,y,z)**2
    versor = [x/normalization, y/normalization, z/normalization]
    theta = numpy.arctan2(numpy.sqrt(versor[0]**2 + versor[1]**2),versor[2])
    phi = numpy.arctan2(versor[1],versor[0])
    pz = mag*numpy.cos(theta)
    pt = mag*numpy.sin(theta)
    px = pt*numpy.cos(phi)
    py = pt*numpy.sin(phi)
    return px,py,pz

def getTheta(x,y,z):
    return numpy.arctan2(numpy.sqrt(x**2 + y**2), z)

def getPhi(x, y, z):
    return numpy.arctan2(y,x)
