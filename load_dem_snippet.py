v_dg = df.VectorFunctionSpace(self.model.mesh,self.model.E_thk)
X = df.interpolate(self.model.mesh.coordinates,v_dg)
meshx = X.dat.data_ro[:,0]*1000
meshy = X.dat.data_ro[:,1]*1000
with rio.open("./data/frenchjoe_srtm_utm12n.tif") as src:
            self.model.B.dat.data[:] = np.array([pnt[0] for pnt in src.sample(zip(meshx, meshy))])/1000
