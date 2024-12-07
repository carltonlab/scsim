import numpy as np
from numpy.core.umath_tests import inner1d
from numpy.matlib import repmat

def solver(kPos, kAnchor, link0, link1, w0, cycles=1000, precision=0.001, dampening=0.1, debug=False):
    """
    kPos       : vector array - knot position
    kAnchor    : float array  - knot's anchor state, 0 = moves freely, 1 = anchored (not moving)
    link0      : int array    - array of links connecting each knot. each index corresponds to a knot
    link1      : int array    - array of links connecting each knot. each index corresponds to a knot
    w0         : float array  - initial link length
    cycles     : int          - eval stops when n cycles reached
    precision  : float        - eval stops when highest applied force is below this value
    dampening  : float        - keeps system stable during each iteration
    """

    kPos        = np.asarray(kPos)
    pos         = np.array(kPos) # copy of kPos
    kAnchor     = 1-np.clip(np.asarray(kAnchor).astype(float),0,1)[:,None]
    link0       = np.asarray(link0).astype(int)
    link1       = np.asarray(link1).astype(int)
    w0          = np.asarray(w0).astype(float)

    F = np.zeros(pos.shape)
    i = 0

    for i in xrange(cycles):

        # Init force applied per knot
        F = np.zeros(pos.shape)

        # Calculate forces
        AB = pos[link1] - pos[link0] # get link vectors between knots
        w1 = np.sqrt(inner1d(AB,AB)) # get link current lengths
        AB/=w1[:,None] # normalize link vectors
        f = (w1 - w0) # calculate force
        f = f[:,None] * AB # calculate force vector

        # Apply force vectors on each knot
        np.add.at(F, link0, f)      # F[link0] += f*AB
        np.subtract.at(F, link1, f) # F[link1] -= f*AB

        # Update point positions       
        pos += F * dampening * kAnchor

        # If the maximum force applied is below our precision criteria, exit
        if np.amax(F) < precision:
            break

    # Debug info
    if debug:
        print 'Iterations: %s'%i
        print 'Max Force:  %s'%np.amax(F)

    return pos


# Create the chromosome knots and links 
kPos = np.loadtxt('minit.txt')
# as of 2021-04-02 (what?) it's 25 nodes x 4 chromatids

# Define the links connecting each knots
link0 = np.loadtxt('l0.txt').astype(int)
link1 = np.loadtxt('l1.txt').astype(int)
AB    = kPos[link0]-kPos[link1]
w0    = np.sqrt(inner1d(AB,AB)) # this is a square grid, each link's initial length will be 0.25
#w0 = [7.0 for i in range(146)]
#w0=np.array(w0)
#w0[0:96]-=3.0

# Set the anchor states
kAnchor = np.zeros(len(kPos)) # All knots will be free floating
#kAnchor[12] = 1 # Middle knot will be anchored

#print np.allclose(kPos,solver(kPos, kAnchor, link0, link1, w0, debug=True))
# Returns True
# Iterations: 0
# Max Force:  0.0

# Move the center knot up a little
#kPos[12] = np.array([0,0.3,0])

# eval the system
new = solver(kPos, kAnchor, link0, link1, w0, debug=True, cycles=10000) # positions will have moved
print(new)

fin=np.loadtxt('mfinal.txt')
for tryit in range(300):
    finforce=fin-new
#    lf=np.linalg.norm(finforce,axis=1) # get magnitude
#    llf=len(lf)
#    lf=np.reshape(lf,(llf,1))
#    lf=repmat(lf,1,3)
#    finforce/=lf
    finforce/=100 # differently, ok lets try
    kPos+=finforce

    new = solver(kPos, kAnchor, link0, link1, w0, debug=True, cycles=10000) # positions will have moved
    print(new)
