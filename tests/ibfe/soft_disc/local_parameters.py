# Define time parameters
T = 12
dt = 1 / 40
num_steps = int(T / dt)

# Define fluid parameters
rho = 1.0
nu = 0.01

# Define solid parameters


# Define stablization parameters
alpha = 10 * dt
conv = True
stab = True
delta = 1000.0
SAV = 1.0

# 32 40
# 64 80
# 128 160
n_mesh_fluid = 32
n_mesh_solid = 32
nu_s = 10
order_velocity = 2
order_pressure = 1
order_displacement = 1


h = 1.0 / n_mesh_fluid