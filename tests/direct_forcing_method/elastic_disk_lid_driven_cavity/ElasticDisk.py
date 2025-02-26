from fenics import *
from local_parameters import dt

class ElasticDisk:
    eta = 0.05/dt
    kappa = 0.01/dt/dt
    def __init__(self, W, dt, nu_s, rho):
        self.W = W
        self.xi = Function(W)
        self.xi_ = Function(W)
        self.xi_n = Function(W)
        self.xi.interpolate( Expression(("x[0]", "x[1]"), degree=2))
        self.xi_.interpolate( Expression(("x[0]", "x[1]"), degree=2))
        self.xi_n.interpolate( Expression(("x[0]", "x[1]"), degree=2))
        self.disp_n = Function(W) # 需要初始化一下
        self.dt = dt
        self.Sigma_part_1 = Function(W)
        self.Sigma_part_2 = Function(W)
        Sigma = as_tensor([[self.Sigma_part_1[0], self.Sigma_part_1[1]],[self.Sigma_part_2[0],self.Sigma_part_2[1]]])
                
        N = FacetNormal(W.mesh())

        # Define variables
        dxi = TrialFunction(W)
        v = TestFunction(W)
        
        k = Constant(dt)
        rho = Constant(rho)
        nu_s = Constant(nu_s)
        
        # 计算固体方程
        F = grad(self.xi_)
        C = F.T*F
        P1 = nu_s/det(F) * ( F - 0.5*tr(C)*inv(F).T)
        nv = 0.499
        K_s = 2.0*nu_s*(1.0+nv)/3/(1-2*nv)
        J = det(F)
        P2 = K_s*J*ln(J)*inv(F.T)
        P = P1 + P2
        j = J*sqrt(inner(inv(F).T* N,inv(F).T* N))
        # F = rho/k/k*inner((self.xi_-2*self.xi + self.xi_n),v)*dx + inner(P,grad(v))*dx - inner(det(F)*Sigma*inv(F.T)*N, v)*ds
        F = rho/k/k*inner((self.xi_-2*self.xi + self.xi_n),v)*dx + inner(P,grad(v))*dx - inner(j*Sigma*inv(F.T)*N, v)*ds
        J = derivative(F, self.xi_, dxi)

        self.F = F
        self.J = J
        

    # # 更新上一步的位移和速度
    # def update_displacement(self, disp):
    #     self.xi.assign(disp)
        
    # def update_displacement_before(self, disp):
    #     self.xi_n.assign(disp)
    
    def solve(self):
        # 求解出固体的位移self.u和速度self.v和拉格朗日乘子self.p
        solve(
            self.F == 0,
            self.xi_,
            J=self.J,
            solver_parameters={
                "newton_solver": {
                    "linear_solver": "mumps",
                    "absolute_tolerance": 1.0e-7,
                    "relative_tolerance": 1e-7,
                }
            },
        )
        self.xi_n.assign(self.xi)
        self.xi.assign(self.xi_)
    
    def pently_force(self, disp):
        print("PENTLY FORCE")
        F = Function(self.W)
        print(len(F.vector()))
        # kappa*(xi-disp) + eta/dt*(xi-xi_n-disp+disp_n)
        for i in range(len(F.vector())):
            F.vector()[i] = self.kappa*(self.xi.vector()[i]-disp.vector()[i]) 
            + self.eta/self.dt*(self.xi.vector()[i]-self.xi_n.vector()[i]-disp.vector()[i]+self.disp_n.vector()[i])  
        self.disp_n.assign(disp)
        return F
