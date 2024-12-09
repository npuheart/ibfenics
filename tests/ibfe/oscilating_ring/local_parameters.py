# Define time parameters
T = 1.2
dt = 1 / 800
num_steps = int(T / dt)

# Define fluid parameters
rho = 1.0
nu = 0.001

# Define stablization parameters
alpha = 10.0 * dt
conv = True
stab = True
delta = 0.1
SAV = 1.0

# 32 40
# 64 80
# 128 160
n_mesh_fluid = 32
n_mesh_solid = 40
h = 1.0 / n_mesh_fluid
nu_s = 1.0 * 1.0 / 0.0625
order_velocity = 2
order_pressure = 1
order_displacement = 1
