include("temperature_render.jl")

<<<<<<< HEAD
angles = 75.0:1.0:85.0
anim = @animate for angle in angles
    fangle = Float64(angle)
    temperature_render(obs_angle=fangle, tolerance = 1e-4)
end
gif(anim, "hello.gif", fps=10)
=======
# angles = 75.0:1.0:85.0
# angles = [-179.0:1.0:-1.0;1.0:1.0:179.0]
masses = 1.0:5.0:200
n_frames = length(masses)
anim = @animate for (i, mass) in enumerate(masses)
    print("Frame: $i / $n_frames complete")
    # fangle = Float64(angle)
    fmass = Float64(mass)
    temperature_render(mass=mass, tolerance=1e-10, resolution=1080, size_multiplier=4, dtmax=2.5)
end
gif(anim, "massgif.gif", fps=20)
>>>>>>> refs/remotes/origin/main
