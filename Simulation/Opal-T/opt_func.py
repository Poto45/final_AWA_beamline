import numpy as np
import scipy
from scipy import linalg as la
from sklearn.neighbors import KernelDensity

# Prepare Input
def getMesh(N, type='t'):
    '''
    compute Nx, Ny, Nz for IMPACT-T input. 
    Number of macro-particles should be less than five in each cell.
    Use 't' to focus on transverse plane, 'l' to focus on longitudinal plane. 
    '''
    i = 2
    if type=='t':
        while 2**(3*i-1)<N/5:i+=1
        return 2**i,2**i,2**(i-1)
    elif type=='l':
        while 2**(3*i-2)<N/5:i+=1
        return 2**(i-1),2**(i-1),2**i
    else:
        raise ValueError('Wrong type. Use t or l')
        
def parametersRestore(x, rng):
    '''
    x is a 1d array of 0~1 (deap individual)
    rng is the ranges of parameters
    return the random sample drawn from the rng(the real ranges)
    '''
    a = np.diff(rng).flatten()
    b = np.min(rng, axis=1)
    return a * x + b


# Fitness function and output

def rcim(stp,col):
    '''
    Read Columns from Impact-t output. Col starts at 1.  
    '''
    return np.loadtxt('fort.' + str(stp), usecols=col - 1)

def rcas(root,ext,stp,col):
    '''
    Read Columns from ASTRA output. Col starts at 1.
    assumes filenames are root.ext.stp
    ext=Xemit, Yemit, ...
    stp is the run number
    '''
    return np.loadtxt(root+'.' + ext +'.' + str(stp).zfill(3), usecols=col - 1)


def linearRamp(z, bins=100):
    '''
    output: skewness, linearity
    '''
    (zy, bin_edges) = np.histogram(z, bins)
    zx = .5*(bin_edges[1:]+bin_edges[:-1])
    b = zy > 3        # remove bins with fewer counts
    zx = zx[b]
    zy = zy[b]
    A = np.vstack([zx, np.ones(len(zx))]).T
    m, c = np.linalg.lstsq(A, zy, rcond=None)[0]
    zy_fit = m * zx + c
    d = scipy.spatial.distance.euclidean(zy,zy_fit) / bins
    return scipy.stats.skew(z), d
    
def powerTriangle(mu, bounds, N):
    """return powerlaw distribution with mu, between bounds(a,b)"""
    assert isinstance(bounds, tuple)
    assert bounds[0] < bounds[1]
    a = bounds[0]
    b = bounds[1]
    T = -powerlaw.rvs(mu+1,scale=b-a,size=N)
    T = T - np.max(T) + a
    return T

def calWass(t):
    t2 = powerTriangle(2, (np.min(t), np.max(t)), len(t))
    return wasserstein_distance(t,t2) / (np.max(t) - np.min(t))
       
def magBeam(fortout):
    disti = np.loadtxt(fortout)
    # Calculate beam matrix

    # sigma6=np.zeros((6,6))
    # for i in range(0,6):
    #   for j in range(0,6):
    #      sigma6[i,j]=np.mean(dist[:,i]*dist[:,j])

    betagamma = np.sqrt(np.mean(disti[:, 1]**2+disti[:, 3]**2+disti[:, 5]**2))
    gamma = np.sqrt(1.+betagamma**2)
    beta = np.sqrt(1.-1/gamma**2)

    dist = disti.copy()

    # convert from impact format
    dist[:, 1] = disti[:, 1]/disti[:, 5]
    dist[:, 3] = disti[:, 3]/disti[:, 5]

    sigma4 = np.zeros((4, 4))
    for i in range(0, 4):
        for j in range(0, 4):
            sigma4[i, j] = np.mean(dist[:, i]*dist[:, j])

    # Correlation matrix  C = <YX>.<XX>^-1
    YX = np.array([[sigma4[2, 0], sigma4[2, 1]], [sigma4[3, 0], sigma4[3, 1]]])
    XY = np.array([[sigma4[0, 2], sigma4[0, 3]], [sigma4[1, 2], sigma4[1, 3]]])
    XX = np.array([[sigma4[0, 0], sigma4[0, 1]], [sigma4[1, 0], sigma4[1, 1]]])
    YY = np.array([[sigma4[2, 2], sigma4[2, 3]], [sigma4[3, 2], sigma4[3, 3]]])

    C = np.dot(YX, np.linalg.inv(XX))
    # print("Correlation matrix determinant: " + str(np.linalg.det(C)))
    # print("Sigma4 matrix determinant "+str(np.linalg.det(sigma4)))

    # Eigenemittance calculation
    null2x2 = np.array([[0, 0], [0, 0]])
    J = np.array([[0, 1], [-1, 0]])
    J4 = np.vstack((np.hstack((J, null2x2)), np.hstack((null2x2, J))))

    emittances = beta*gamma*np.abs(np.imag(np.linalg.eigvals(np.dot(sigma4, J4))))
    # print(emittances)

    # print("emittances")
    # print(beta*gamma*np.sqrt(np.abs(np.linalg.det(XX))))
    # print(beta*gamma*np.sqrt(np.abs(np.linalg.det(YY))))

    ## Reminder: emittance4d = Sigma4 matrix determinant**0.25*betagamma 
    ## OR emittance4d = sqrt(emittances)
    emittance4d = np.sqrt(emittances[0]*emittances[2])
    ratio = np.divide(emittances[0],emittances[2])
    if ratio < 1:
        ratio = 1. / ratio 
        
    return emittance4d, ratio    
    
    
def getTwiss():
    # load x,y,z emit, filename is fort.31
    sigma = np.dtype({'names': ['t', 'z',
                                'x2', 'xpx', 'xy',  'xpy',
                                'px2', 'pxy', 'pxpy',
                                'y2',  'ypy',
                                'py2'],
                      'formats': [np.double, np.double,
                                  np.double, np.double, np.double, np.double,
                                  np.double, np.double, np.double,
                                  np.double, np.double,
                                  np.double]})
    S = np.loadtxt('fort.31', dtype=sigma,skiprows=1)
    z, rel_gamma,rel_beta = np.loadtxt('fort.18', usecols=(1,2,4), skiprows=1,unpack=True)
    # Px in IMPACT-T is actually Px/MC. This is to convert Px to x'
    S['xpx'] = S['xpx'] / (rel_gamma * rel_beta)
    S['xpy'] = S['xpy'] / (rel_gamma * rel_beta)
    S['px2'] = S['px2'] / (rel_gamma * rel_beta)**2
    S['pxy'] = S['pxy'] / (rel_gamma * rel_beta)
    S['pxpy'] = S['pxpy'] / (rel_gamma * rel_beta)**2
    S['ypy'] = S['ypy'] / (rel_gamma * rel_beta)
    S['py2'] = S['py2'] / (rel_gamma * rel_beta)**2
    
    twiss_beta_x = S['x2'] / np.sqrt(S['x2'] * S['px2'] - S['xpx'] * S['xpx'])
    twiss_alpha_x = - S['xpx'] / np.sqrt(S['x2'] * S['px2'] - S['xpx'] * S['xpx'])
    twiss_beta_y = S['y2'] / np.sqrt(S['y2'] * S['py2'] - S['ypy'] * S['ypy'])
    twiss_alpha_y = - S['ypy'] / np.sqrt(S['y2'] * S['py2'] - S['ypy'] * S['ypy'])
    
    return z, twiss_alpha_x, twiss_beta_x, twiss_alpha_y, twiss_beta_y
    
def LoadImpactTSigma(filename):

    # load x,y,z emit
    # Note, the sigma here is normalized
    sigma = np.dtype({'names': ['t', 'z',
                                'x2', 'xpx', 'xy',  'xpy',
                                'px2', 'pxy', 'pxpy',
                                'y2',  'ypy',
                                'py2'],
                      'formats': [np.double, np.double,
                                  np.double, np.double, np.double, np.double,
                                  np.double, np.double, np.double,
                                  np.double, np.double,
                                  np.double]})
    
    S = np.loadtxt(filename, dtype=sigma,skiprows=1)

    S4 = np.zeros((4, 4))
    J4 = np.zeros((4, 4))

    J4 = [[0., 1., 0., 0.],
          [-1., 0., 0., 0.],
          [0., 0., 0., 1.],
          [0., 0., -1., 0.]]
    enx = np.zeros((len(S['z'])))
    eny = np.zeros((len(S['z'])))
    emit_4d = np.zeros((len(S['z'])))
    curlyL = np.abs(S['xpy'] - S['pxy'])/2
    
    for jj in range(len(S['z'])):
        S4 = [[S['x2'][jj],   S['xpx'][jj],  S['xy'][jj],  S['xpy'][jj]],
              [S['xpx'][jj],  S['px2'][jj],  S['pxy'][jj], S['pxpy'][jj]],
              [S['xy'][jj],   S['pxy'][jj],  S['y2'][jj],  S['ypy'][jj]],
              [S['xpy'][jj],  S['pxpy'][jj], S['ypy'][jj], S['py2'][jj]]]

        vals, vecs = la.eig(np.dot(J4, S4))

        enx[jj] = np.abs(np.imag(vals[0]))
        eny[jj] = np.abs(np.imag(vals[2]))
        emit_4d[jj] = np.sqrt(enx[jj]*eny[jj])

    return S, enx, eny, emit_4d, curlyL
    
def compare_triangle(data):
    """calculate the Wasserstein distance between the input distribution and triangle"""
    # rescale data to 0 ~ 1
    data = ( data - data.min() ) / ( data.max() - data.min() )
    # padding
    dz = 1/2000
    z = np.arange(data.min()-data.std(), data.max()+data.std(), dz)
    # KDE
    bandwidth = (4 * data.std()**5/3 /data.size)**(1/5)
    kde = KernelDensity(bandwidth=bandwidth, kernel='gaussian', rtol=1E-8, atol=1E-8)
    kde.fit(data[:, None])
    logprob = kde.score_samples(z[:, None])
    data_pdf = np.exp(logprob)
    data_cdf = np.cumsum(data_pdf) * dz
    # ideal distribution
    triangle = lambda a, x: ((x >= 0) & (x <=1)) * (a + 1) * (1-x)**a
    
    tri_pdf_1 = triangle(1, z)
    tri_pdf_2 = triangle(2, z)
    
    tri_cdf_1  = np.cumsum(tri_pdf_1)  * dz
    tri_cdf_2  = np.cumsum(tri_pdf_2)  * dz
    
    w1 = np.sum(np.multiply(np.abs(data_cdf - tri_cdf_1), dz))
    w2 = np.sum(np.multiply(np.abs(data_cdf - tri_cdf_2), dz))
    w = (w1 > w2) * w2 + (w1 < w2) * w1
    tri_pdf = (w1 > w2) * tri_pdf_2 + (w1 < w2) * tri_pdf_1

    return w
