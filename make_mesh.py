import fiona
import numpy as np
import gmsh

shp = fiona.open("frenchjoe.gpkg", mode="r")
coords = np.array(list(shp.values())[1]["geometry"]["coordinates"]).squeeze()
shp.close()

# Resolution
lc = 200

# Trim last point (is identical to first)
coords = coords[:-1, :]

gmsh.initialize()

gmsh.model.add("domain")

points = []

for x, y in coords:
    points.append(gmsh.model.geo.addPoint(x, y, 0, lc))


lines = []
for i in range(len(points)):
    lines.append(gmsh.model.geo.addLine(points[i - 1], points[i]))

cl = gmsh.model.geo.addCurveLoop(lines)

ps = gmsh.model.geo.addPlaneSurface([cl], 1)

gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(2)

gmsh.write("domain.msh")

gmsh.finalize()
