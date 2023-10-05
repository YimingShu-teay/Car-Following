#动力学模型
class KinematicModel:
    def __init__(self,x,v,dt,f0,f1,f2,m):
        self.x = x
        self.v = v
        self.dt = dt
        self.f0 = f0
        self.f1 = f1
        self.f2 = f2
        self.m = m

    def update_state(self,a_x):
        self.x = self.x + self.v*self.dt
        self.v = self.v + (a_x/self.m - self.v*(self.f1+2*self.f2*self.v)/self.m)*self.dt

    def get_state(self):
        return self.x, self.v
    
    
