# # Standard units
# units = UserIU(g=0.001, cm=0.01)

# Define time parameters
T = 12.0
dt = 1 / 8000
num_steps = int(T / dt)

# Define fluid parameters
rho = 1.0
nu = 0.001

# Define solid parameters
beta_s = 1e6
kappa_stab = 0.0
G_s = 32

# Define stablization parameters
alpha = 10.0 * dt
conv = True
stab = False
delta = 0.1
SAV = 1.0

# 32 40
# 64 80
# 128 160
n_mesh_fluid = 32
n_mesh_solid = 40
order_velocity = 2
order_pressure = 1
order_displacement = 1
