import icepack
import firedrake
import rasterio
import numpy as np
import tqdm


mesh = firedrake.Mesh("domain.msh")

Q = firedrake.FunctionSpace(mesh, family="CG", degree=3)
V = firedrake.VectorFunctionSpace(mesh, family="CG", degree=3)

W = firedrake.VectorFunctionSpace(mesh, family="CG", degree=3)
x = firedrake.interpolate(mesh.coordinates, W)
meshx = x.dat.data_ro[:, 0]
meshy = x.dat.data_ro[:, 1]

b = firedrake.Function(Q)
with rasterio.open("./output_SRTMGL3_32159.tif") as src:
    b.dat.data[:] = np.array([pnt[0] for pnt in src.sample(zip(meshx, meshy))])

firedrake.VTKFile("b.pvd").write(b)

h0 = firedrake.Function(Q).interpolate(firedrake.Constant(0.0))
s0 = firedrake.Function(Q).interpolate(b + h0)

model = icepack.models.ShallowIce()

solver = icepack.solvers.FlowSolver(model)

T = firedrake.Constant(273.15 - 1)
A = icepack.rate_factor(T)

u0 = firedrake.Function(V)
h = h0.copy(deepcopy=True)
u = solver.diagnostic_solve(
    velocity=u0,
    thickness=h,
    surface=s0,
    fluidity=A,
)


def mass_balance(s, max_a=0.5, da_ds=0.5 / 1000, ela=300.0):
    return firedrake.min_value((s - ela) * da_ds, max_a)


ela = 2800.0
max_a = 5
da_ds = 0.5 / 100

a = mass_balance(s0, ela=ela, max_a=max_a, da_ds=da_ds)

firedrake.VTKFile("a.pvd").write(firedrake.interpolate(a, Q))

dt = 0.1
num_timesteps = 4000

dh_max = np.zeros(num_timesteps) * np.nan
a = firedrake.Function(Q)

hfile = firedrake.VTKFile("h.pvd")
hfile.write(h0, time=0)

ufile = firedrake.VTKFile("u.pvd")
ufile.write(u0, time=0)

for step in tqdm.trange(num_timesteps):
    h = solver.prognostic_solve(
        dt,
        thickness=h,
        accumulation=a,
        velocity=u,
    )

    h.interpolate(firedrake.max_value(h, 0))
    s = firedrake.Function(Q).interpolate(h + b)
    u = solver.diagnostic_solve(
        velocity=u,
        thickness=h,
        surface=s,
        fluidity=A,
    )

    a.interpolate(mass_balance(s, ela=ela, max_a=max_a, da_ds=da_ds))

    if not step % 10:
        hfile.write(h, time=dt * (step + 1))
        ufile.write(u, time=dt * (step + 1))
