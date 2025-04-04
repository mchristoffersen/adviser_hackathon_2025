import firedrake as df

W = df.VectorFunctionSpace(mesh, family="CG", degree=2)
x = df.interpolate(mesh.coordinates, W)
meshx = x.dat.data_ro[:,0]
meshy = x.dat.data_ro[:,1]

b = df.function(Q)
with rio.open("./data/frenchjoe_srtm_utm12n.tif") as src:
            b.dat.data[:] = np.array([pnt[0] for pnt in src.sample(zip(meshx, meshy))])
