from dolfin import *
from mshr import *
import numpy as np 
from scipy.optimize import fsolve


def energy(u, E_hyp, rho):
    return 0.5*rho*assemble(inner(u, u)*dx)  + E_hyp

def modified_energy(u, E_hyp, rho, delta):
    tmp = energy(u, E_hyp, rho) + delta # + 弹性势能
    return tmp

def SAV_X1(un, u1, u2, E_hyp, dt, mu, delta, rho):
    return 2.0/dt*modified_energy(un, E_hyp, rho, delta) + mu*assemble(inner(grad(u2), grad(u2))*dx)

def SAV_X2(un, u1, u2, E_hyp, dt, mu, delta, rho, qn):
    return - 2.0*qn/dt*np.sqrt(modified_energy(un, E_hyp, rho, delta)) + 2.0*mu*assemble(inner(grad(u1), grad(u2))*dx)

def SAV_X3(un, u1, u2, E_hyp, dt, mu, delta, rho):
    return mu*assemble(inner(grad(u1), grad(u1))*dx)

def calculate_SAV(u0, u1, u2, E, dt, nu, delta, qn):
    X1 = SAV_X1(u0, u1, u2, E, dt, nu, delta, 1.0)
    X2 = SAV_X2(u0, u1, u2, E, dt, nu, delta, 1.0, qn)
    X3 = SAV_X3(u0, u1, u2, E, dt, nu, delta, 1.0)
    # 定义目标函数
    def func(x):
        return X1*x*x + X2*x + X3
    # 使用 fsolve 求解
    # initial_guess = 1.0  # 初始猜测值
    # solution = fsolve(func, initial_guess)
    # print(f"解是: {solution[0]}")
    # return solution[0]
    # 直接求解
    S = (-X2 + np.sqrt(X2*X2-4.0*X1*X3))/(2*X1)
    if np.isnan(S):
        S = -X2/2.0/X1
    return S





# from dolfin import *
# import numpy as np

def calculate_SAV_2(h, dt, mu, En, qn, qnm1, unp1, un, unm1, u1, u2, p1, p2, delta, alpha, n, rho):
    # h = 0.1
    # dt = 0.1
    # mu = 0.1
    # En = 0.1
    # qn = 0.1
    # qnm1 = 0.1

    # unp1 = 0.1 # velocity
    # un = 0.1
    # unm1 = 0.1
    # u1 = 0.1
    # u2 = 0.1
    # p1 = 0.1
    # p2 = 0.1
    # delta = 0.1
    # alpha = 0.1
    # n = 0.1 # normal vector
    # rho = 0.1

    _0 = 3.0/dt*(En + delta)
    _1 = (mu+3.0*alpha*h)*assemble(inner(grad(u2), grad(u2))*dx)
    _2 = 1.0/dt*(-4.0*qn + qnm1)*np.sqrt(En + delta)
    _3 = 2.0*mu*assemble(inner(grad(u1), grad(u2))*dx)
    _4 = alpha*h*assemble(inner(3.0*grad(u1)-4.0*grad(un)+grad(unm1), grad(u2))*dx)
    _5 = 3.0*alpha*h*assemble(inner(grad(u2), grad(u1))*dx)
    _6 = mu*assemble(inner(dot(n, grad(u2)),unp1)*ds)
    _7 = assemble(inner(n, u2)*p2*ds) 
    _8 = 3.0*alpha*h*assemble(inner(dot(n, grad(u2)),unp1)*ds)
    _9 = mu*assemble(inner(grad(u1), grad(u1))*dx)
    _10 = alpha*h*assemble(inner(3.0*grad(u1)-4.0*grad(un)+grad(unm1), grad(u1))*dx)
    _11 = mu*assemble(inner(dot(n, grad(u1)),unp1)*ds)
    _12 = alpha*h*assemble(inner(dot(n, 3.0*grad(u1) - 4.0*grad(un) + grad(unm1)), unp1)*ds)
    _13 = assemble(inner(n, unp1)*p1*ds)
    _14 = 0.5*rho*assemble(inner(n, unp1)*inner(unp1,unp1)*ds)

    A = _0 + _1
    B = _2 + _3 + _4 + _5 - _6 + _7 - _8 
    C = _9 + _10 - _11 - _12 + _13 + _14

    S = (-B+np.sqrt(B*B-4.0*A*C))/2.0/A
    return S



class TaylorHoodSolverBDF2_1:
    def __init__(self, u0_, u0, p0, dt, nu, stab=False, alpha = 0.1):
        # Reconstruct element space
        mesh = u0.function_space().mesh()
        element1 = u0.function_space()._ufl_element
        element2 = p0.function_space()._ufl_element
        TH = element1*element2
        W = FunctionSpace(mesh, TH)

        # Define variables
        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)
        k = Constant(dt)
        self.w_ = Function(W)
        self.un , self.pn = Function(W).split(True)
        self.un_, self.pn_ = Function(W).split(True)
        
        # dt approximation
        du_dt = (3.0*u-4.0*self.un+self.un_)/(2.0*k)
        u_hat = 2.0*self.un - self.un_
        u_grad_u = grad(u_hat)*u_hat

        # Define variational problem
        F = inner(du_dt, v)*dx + nu * inner(grad(u), grad(v))*dx - div(v)*p * dx # - inner(f, v)*dx
        if stab:
            F += alpha * inner(grad(u-self.un), grad(v))*dx  # 稳定项
        F += q*div(u)*dx
        a = lhs(F)
        self.W = W
        self.A = assemble(a)
        self.L = rhs(F)
        self.un.assign(u0_)

    def update(self, un, pn):
        self.un_.assign(self.un) # u_{n-1} = u_n
        self.un.assign(un)       # u_n = u_{n+1}
        self.pn.assign(pn)       # p_n = p_{n+1}

    def solve(self, bcu, bcp):
        b = assemble(self.L)
        [bc.apply(self.A, b) for bc in bcu]
        [bc.apply(self.A, b) for bc in bcp]
        solve(self.A, self.w_.vector(), b)
        return self.w_.split(True)



class TaylorHoodSolverBDF2_2:
    def __init__(self, u0_, u0, p0, f, dt, nu, stab=False, alpha = 0.1):
        # Reconstruct element space
        mesh = u0.function_space().mesh()
        element1 = u0.function_space()._ufl_element
        element2 = p0.function_space()._ufl_element
        TH = element1*element2
        W = FunctionSpace(mesh, TH)

        # Define variables
        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)
        k = Constant(dt)
        self.w_ = Function(W)
        self.un , self.pn = Function(W).split(True)
        self.un_, self.pn_ = Function(W).split(True)
        
        # dt approximation
        du_dt = (3.0*u)/(2.0*k)
        u_hat = 2.0*self.un - self.un_
        u_grad_u = grad(u_hat)*u_hat

        # Define variational problem
        F = inner(du_dt, v)*dx + nu * inner(grad(u), grad(v))*dx + div(v)*p * dx - inner(f, v)*dx + inner(u_grad_u, v)*dx
        if stab:
            F += alpha * inner(grad(u), grad(v))*dx  # 稳定项
        F += q*div(u)*dx
        a = lhs(F)
        self.W = W
        self.A = assemble(a)
        self.L = rhs(F)
        self.un.assign(u0_)

    def update(self, un, pn):
        self.un_.assign(self.un) # u_{n-1} = u_n
        self.un.assign(un)       # u_n = u_{n+1}
        self.pn.assign(pn)       # p_n = p_{n+1}

    def solve(self, bcu, bcp):
        b = assemble(self.L)
        [bc.apply(self.A, b) for bc in bcu]
        [bc.apply(self.A, b) for bc in bcp]
        solve(self.A, self.w_.vector(), b)
        return self.w_.split(True)


