# Define time parameters
T = 12
dt = 1 / 400
num_steps = int(T / dt)

# Define fluid parameters
rho = 1.0
nu = 0.01

# Define solid parameters


# Define stablization parameters
alpha = 10.0 * dt
conv = True
stab = False
delta = 0.1
SAV = 1.0

# 32 40
# 64 80
# 128 160
n_mesh_fluid = 40
n_mesh_solid = 32
nu_s = 0.1
order_velocity = 2
order_pressure = 1
order_displacement = 1
