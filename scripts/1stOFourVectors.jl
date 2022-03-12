using Revise 
using GeodesicTracer
using CarterBoyerLindquist
using StaticArrays
using Plots
gr()

m = CarterMethodBL(M=1.0, a=1.0) #chosen metric
u = @SVector [0.0, 4.0, deg2rad(90.0), 0.0] #initial position (time, radius, theta, phi)
# v = @SVector [0.0, -1.0, 0.0, sin(deg2rad(0.0))] #initial velocity (time, -1/+1= in/out, theta, phi)
v = map_impact_parameters(m, u, 0.0, 1.0) #initial velocity with better ability to adjust parameters
sols = tracegeodesics(m, u, v, (0.0, 100.0);
    abstol=1e-9,
    reltol=1e-9
)

w = CarterBoyerLindquist.carter_velocity(u, m.E, m.M, m.a, sols.prob.p)
print(w)

plot(sols, vars=(4,2), proj=:polar, range=(0.0, 10.0))

