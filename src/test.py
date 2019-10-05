from numpy import *
from geometry import segment_intersect, Line, Point
def perp( a ) :
    b = empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b

# line segment a given by endpoints a1, a2
# line segment b given by endpoints b1, b2
# return 
def seg_intersect(a1,a2, b1,b2) :
    da = a2-a1
    db = b2-b1
    dp = a1-b1
    dap = perp(da)
    denom = dot( dap, db)
    num = dot( dap, dp )
    return (num / denom.astype(float))*db + b1

p1 = Point( [0.0, 0.0] )
p2 = Point( [1.0, 1.0] )

p3 = Point( [1.0, 0.0] )
p4 = Point( [0.0, 1.0] )

L1 = Line([p1,p2])
L2 = Line([p3,p4])

print segment_intersect(L1, L2) 
