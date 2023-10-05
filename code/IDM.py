#单辆lear car IDM
class IDM:
    def __init__(self,Leader_V,Leader_X,Leader_Y,Leader_a,A_max,dt,V0):
        
        self.Leader_V = Leader_V
        self.Leader_X = Leader_X
        self.Leader_Y = Leader_Y
        self.Leader_a = Leader_a
        self.dt = dt
        self.V0 = V0
        self.A_max = A_max
    
    def update_state(self):
        self.Leader_X = self.Leader_X + self.Leader_V*self.dt
        self.Leader_V = self.Leader_V + self.Leader_a*self.dt   
        self.Leader_a = self.A_max*(1-round((self.Leader_V/self.V0),1)**4)

    def get_state(self):
        return self.Leader_X, self.Leader_V, self.Leader_a
    
    #定义lead car预测序列
    def get_x_lead_list(self,H):
        x_lead_list = [self.Leader_X]
        for item in range(H):
            x_mid = self.Leader_X + self.Leader_V*self.dt
            x_lead_list.append(x_mid)
        return x_lead_list