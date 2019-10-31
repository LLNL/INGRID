import Ingrid
import numpy as np
import yaml

class v():

    def __init__(self):
        self.yaml = yaml.read('../test_params/params3.yaml')
        self.xpt = np.array([self.yaml['grid_params']['rxpt'], self.yaml['grid_params']['zxpt']])
        self.mag = np.array([self.yaml['grid_params']['rmagx'], self.yaml['grid_params']['zmagx']])
        self.grid = Ingrid.Ingrid(self.yaml)
        self.grid.setup()
        self.grid.OMFIT_read_psi()
        self.grid.calc_psinorm()
        
        self.psi = self.grid.psi_norm
        self.eps = 0.001
    
    def vr(self,x,y):
        return (self.v(x + self.eps, y) - self.v(x -self.eps, y))/(2 * self.eps)
    
    def vz(self,x,y):
        return (self.v(x, y + self.eps) - self.v(x,y -self.eps))/(2 * self.eps)
    def vrz(self,x,y):
        numer = self.v(x + self.eps, y +self.eps) \
                - self.v(x + self.eps, y - self.eps) \
                - self.v(x - self.eps, y + self.eps) \
                + self.v(x - self.eps, y - self.eps)
        denom = 4 * self.eps ** 2
        return numer / denom
    def vrr(self,x,y):
        return (self.v(x + self.eps, y) - 2 * self.v(x,y) + self.v(x -self.eps, y))/(self.eps ** 2)
    def vzz(self,x,y):
        return (self.v(x, y+self.eps) - 2 * self.v(x,y) + self.v(x, y -self.eps))/(self.eps ** 2)
    def v(self,x,y):
        return self.psi.get_psi(x,y)
    def lap(self,x,y):
        return np.array([self.vrr(x,y), self.vzz(x,y)])
    
    def check_derivs(self,x,y):
        print('Psi values: {}'.format(self.v(x,y)))
        print('Finite Difference fr: {}'.format(self.vr(x,y)))
        print('Ingrid vr: {}'.format(self.psi.get_psi(x,y, tag = 'vr')))
        print('Finite Difference fz: {}'.format(self.vz(x,y)))
        print('Ingrid vz: {}'.format(self.psi.get_psi(x,y, tag = 'vz')))
        print('Finite Difference frz: {}'.format(self.vrz(x,y)))
        print('Ingrid vrz: {}'.format(self.psi.get_psi(x,y, tag = 'vrz')))
        print('Finite Difference frr: {}'.format(self.vrr(x,y)))
        print('Ingrid vrr: {}'.format(self.psi.get_psi(x,y, tag = 'vrr')))
        print('Finite Difference fzz: {}'.format(self.vzz(x,y)))
        print('Ingrid vzz: {}'.format(self.psi.get_psi(x,y, tag = 'vzz')))

    def get_theta(self, x, y):
        res = (2 * self.psi.get_psi(x,y,tag='vrz') /(self.psi.get_psi(x,y,tag='vrr') - self.psi.get_psi(x,y,tag='vzz')))
        theta = np.zeros(4)

        for i in range(4):
            print(res)
            theta[i-1] = np.array([np.arctan(res) * 1/2]) + np.pi/2*(i-1)
            print(np.tan(2*theta[-1]))
        return theta


psi = v()
print('Checking primary x-point:\n')
psi.check_derivs(psi.xpt[0], psi.xpt[1])

print('Checking magnetic axis:\n')
psi.check_derivs(psi.mag[0], psi.mag[1])

print('Get theta:\n')
theta = psi.get_theta(psi.xpt[0], psi.xpt[1])
print(theta)

for i in range(len(theta)):
    res = -psi.psi.get_psi(psi.xpt[0], psi.xpt[1],tag='vrr') * np.cos(2*theta[i]) \
          +psi.psi.get_psi(psi.xpt[0], psi.xpt[1],tag='vzz') * np.cos(2*theta[i]) \
          -2 *2 *psi.psi.get_psi(psi.xpt[0], psi.xpt[1],tag='vrz') * np.sin(2*theta[i])
    print(np.sign(res))
for i in range(len(theta)):
    print((-1/2 * psi.psi.get_psi(psi.xpt[0],psi.xpt[1], tag='vrr') * np.sin(2*theta[i])) \
            + (1/2 * psi.psi.get_psi(psi.xpt[0],psi.xpt[1],tag='vzz')*np.sin(2*theta[i])) \
            + (psi.psi.get_psi(psi.xpt[0],psi.xpt[1],tag='vrz') * np.cos(2*theta[i])))

for i in range(len(theta)):
    print((-1/2 * psi.psi.get_psi(psi.xpt[0],psi.xpt[1], tag='vrr') * np.sin(2*theta[i])) \
            + (1/2 * psi.psi.get_psi(psi.xpt[0],psi.xpt[1],tag='vzz')*np.sin(2*theta[i])) \
            + (psi.psi.get_psi(psi.xpt[0],psi.xpt[1],tag='vrz') * np.cos(2*theta[i])))

