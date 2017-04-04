import numpy as np
import mshr
import dolfin as df

# mesh = mshr.Sphere(df.Point(0, 0, 0), 3)
mesh = mshr.Sphere(df.Point(0, 0, 0), 3)
mesh = mshr.generate_mesh(mesh, 8)

df.plot(mesh, interactive=True)

VS = df.VectorFunctionSpace(mesh, 'CG', 1, 3)
# VF = df.Function(VS)

# Hedgehog-Vortex:
out = df.Expression(["sin(pi * sqrt(x[0] * x[0] + x[1] * x[1])) * cos(atan2(x[1], x[0]))",
                     "sin(pi * sqrt(x[0] * x[0] + x[1] * x[1])) * sin(atan2(x[1], x[0]))",
                     "-cos(pi * sqrt(x[0] * x[0] + x[1] * x[1]))"])
VF = df.interpolate(out, VS)

# VF.vector().set_local(u)

df.plot(VF, interactive=True)

_out = df.File('sphere.pvd')
_out << VF
