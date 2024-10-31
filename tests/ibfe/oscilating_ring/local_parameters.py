# Define time parameters
T = 1.0
dt = 1/8000
num_steps = int(T/dt)

# Define fluid parameters
rho = 1.0
nu  = 0.001

# Define stablization parameters
alpha = 10.0*dt
stab  = False
delta = 0.1
SAV   = 1.0

# 32 40
# 64 80
# 128 160
n_mesh_fluid = 32
n_mesh_solid = 40
nu_s = 1.0/0.0625 
order_velocity = 2
order_pressure = 1
order_displacement = 1