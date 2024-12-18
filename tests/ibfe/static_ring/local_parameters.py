# Define time parameters
T = 0.000001
num_steps = 10
dt = T/num_steps

# Define fluid parameters
rho = 1.0
nu = 0.01

# Define solid parameters


# Define stablization parameters
alpha = 1.0 * dt
conv = True
stab = False
delta = 0.1
SAV = 1.0

# 32 40
# 64 80
# 128 160
n_mesh_fluid = 32
n_mesh_solid = 40
nu_s = 0.2
order_velocity = 2
order_pressure = 1
order_displacement = 1


h = 1.0 / n_mesh_fluid