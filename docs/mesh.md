# 网格

## 1. 左心室网格

### 1.1. 网格生成

左心室网格由[fenicsx-ldrb]()生成，需要进入docker容器执行代码。

`sudo docker run --rm -v $PWD:/home/shared -w /home/shared -it ghcr.io/finsberg/ldrb:latest bash`

```python
import dolfin as df
import ldrb
import cardiac_geometries
geometry = cardiac_geometries.create_lv_ellipsoid()
mesh = geometry.mesh
ffun = geometry.ffun

with df.XDMFFile(mesh.mpi_comm(), "mesh.xdmf") as xdmf:
    xdmf.write(mesh)


df.File("boundaries.xml") << ffun

markers = geometry.markers
markers = {
    "base": geometry.markers["BASE"][0],
    "epi": geometry.markers["EPI"][0],
    "lv": geometry.markers["ENDO"][0],
}
angles = dict(
    alpha_endo_lv=60,  # Fiber angle on the endocardium
    alpha_epi_lv=-60,  # Fiber angle on the epicardium
    beta_endo_lv=0,  # Sheet angle on the endocardium
    beta_epi_lv=0,  # Sheet angle on the epicardium
)
fiber_space = "Lagrange_2"
fiber, sheet, sheet_normal = ldrb.dolfin_ldrb(
    mesh=mesh,
    fiber_space=fiber_space,
    ffun=ffun,
    markers=markers,
    **angles,
)

file1 = df.XDMFFile(mesh.mpi_comm(), "microstructure.xdmf")
file1.parameters['rewrite_function_mesh'] = False
file1.parameters["functions_share_mesh"] = True
file1.parameters["flush_output"] = True

t = 0.0
file1.write_checkpoint(fiber,  "fiber", t, df.XDMFFile.Encoding.HDF5, True)
file1.write_checkpoint(sheet,  "sheet", t, df.XDMFFile.Encoding.HDF5, True)  
file1.write_checkpoint(sheet_normal, "normal", t, df.XDMFFile.Encoding.HDF5, True)  
file1.close()
```

### 1.2. 读取网格、边界和纤维
```python
mesh_path = "/kokkos_2/mesh-benchmark/B009-ldrb/"
order=2

def read_fibers(mesh, data_file, order, index = 0):
    W = VectorFunctionSpace(mesh, 'P', order)
    f0 = Function(W)
    s0 = Function(W)
    n0 = Function(W)
    f_in = XDMFFile(mesh.mpi_comm(), data_file)
    f_in.read_checkpoint(f0, "fiber", 0)
    f_in.read_checkpoint(s0, "sheet", 0)
    f_in.read_checkpoint(n0, "normal", 0)
    f_in.close()
    return f0,s0,n0

def read_mesh_general(mesh_path):
    mesh_file = XDMFFile(mesh_path)
    mesh = Mesh()
    mesh_file.read(mesh)
    mesh_file.close()
    return mesh

mesh = read_mesh_general(mesh_path+"mesh.xdmf")
bdry = MeshFunction("size_t", mesh, mesh_path + "boundaries.xml")
f0,s0,n0 = read_fibers(mesh, mesh_path+"microstructure.xdmf", order, 0):
```