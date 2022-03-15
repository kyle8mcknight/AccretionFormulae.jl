using GeodesicRendering
using AccretionGeometry
using CarterBoyerLindquist
using AccretionFormulae
using StaticArrays
<<<<<<< HEAD
=======
using Printf
>>>>>>> refs/remotes/origin/main

using Plots
gr()

function temperature(m, sol, max_time; kwargs...)
    g =  AccretionFormulae.redshift(m, sol, max_time; kwargs...)

    u = sol.u[end]
<<<<<<< HEAD

    M = m.M*1.99e30
    # M = 5.0
    a = m.a
    a_star = a/M
    temperature = AccretionFormulae.observed_temperature(u[2], a_star, M, g)
end

function temperature_render(mass=1.0, spin=0.998;obs_angle=85.0, disc_angle=90.0, 
                            tolerance=1e-8, dtmax=1000.0, 
                            size_multiplier::Int64=1, resolution=400, fov=3.0
                            )
                            
    m = CarterMethodBL(M=mass, a=spin)
=======
    M = m.M
    a_star = m.a
    temperature = AccretionFormulae.observed_temperature(u[2], a_star, M, g)
    
end

function radius(m, sol, max_time; kwargs...)
    u = sol.u[end]
    radius = u[2]
end

function temperature_render(;mass=1, spin=0.998, obs_angle=85.0, disc_angle=90.0, 
                            tolerance=1e-8, dtmax=1000.0, 
                            size_multiplier::Int64=1, resolution=400, fov=3.0,
                            η=0.1, η_phys=0.1, edd_ratio=0.1, edd_ratio_phys=0.1
                            )
                            
    m = CarterMethodBL(M=1.0, a=spin)
    M = m.M
>>>>>>> refs/remotes/origin/main

    # observer position
    u = [0.0, 1000.0, deg2rad(obs_angle), 0.0]
    R_isco = AccretionFormulae.r_isco(m.a, m.M)

<<<<<<< HEAD
    # new method using the CarterBoyerLindquist RMS gave very different values
    # and wasn't working:
    # R_isco = CarterBoyerLindquist.rms(m.M, m.a)

    # disc has inner radius 10, outer radius 50, perpendicular to the spin axis
    d = GeometricThinDisc(R_isco, 50.0, deg2rad(disc_angle))


    # create and compose the ValueFunction
    redshift_vf = (
        # calculate redshift
        ValueFunction(temperature)
        # filter only pixels with r > 9, i.e. on the disc
        ∘ FilterValueFunction((m, sol, max_time; kwargs...) -> sol.u[end][2] > R_isco, NaN)
        # filter only pixels that terminated early (i.e., intersected with something)
=======
    # disc
    d = GeometricThinDisc(R_isco+1, 50.0, deg2rad(disc_angle))


    # create and compose the ValueFunctions
    temperature_vf = (
        ValueFunction(temperature)
        ∘ FilterValueFunction((m, sol, max_time; kwargs...) -> sol.u[end][2] > R_isco, NaN)
        ∘ ConstValueFunctions.filter_early_term
    )

    radius_vf = (
        ValueFunction(radius)
        ∘ FilterValueFunction((m, sol, max_time; kwargs...) -> sol.u[end][2] > R_isco, NaN)
>>>>>>> refs/remotes/origin/main
        ∘ ConstValueFunctions.filter_early_term
    )

    # do the render
<<<<<<< HEAD
    img = @time rendergeodesics(
=======
    temperature_img = @time rendergeodesics(
>>>>>>> refs/remotes/origin/main
        m, u, 2000.0, 
        d, 
        fov_factor=fov*size_multiplier, abstol=tolerance, reltol=tolerance,
        image_width = 350*size_multiplier,
        image_height = 250*size_multiplier,
<<<<<<< HEAD
        vf = redshift_vf,
        dtmax = dtmax
    )

    # plot
    scale = 1e4
    new_img = reverse(img, dims=1)
    new_img ./= scale

    heatmap(new_img, aspect_ratio=1.0, size=(resolution*3/2, resolution))
    title!("Temperature $scale, $(m.M)")
end

temperature_render(obs_angle=85.0, size_multiplier=2)
=======
        vf = temperature_vf,
        dtmax = dtmax
    )

    radius_img = @time rendergeodesics(
        m, u, 2000.0, 
        d, 
        fov_factor=fov*size_multiplier, abstol=tolerance, reltol=tolerance,
        image_width = 350*size_multiplier,
        image_height = 250*size_multiplier,
        vf = radius_vf,
        dtmax = dtmax
    )

    # correcting for physical mass
    M_phys = mass*1.99e30
    r_isco_phys = AccretionFormulae.r_isco(M_phys, 0.998)
    r_g_phys = M_phys

    numerators =   AccretionFormulae.mass_scale_fraction.(M_phys, η_phys, edd_ratio_phys, r_isco_phys, r_g_phys, radius_img*M_phys)
    denominators = AccretionFormulae.mass_scale_fraction.(M, η, edd_ratio, R_isco, M, radius_img)

    fractions = numerators./denominators
    temperature_img .*= fractions

    # # autoscale
    # scale = maximum(filter(!isnan,temperature_img))
    # scale = floor(log(10, scale))
    # scale = 10^scale
    
    # scaling image
    scale = 1e7
    scalestr = @sprintf "%.E" scale
    new_img = reverse(temperature_img, dims=1)
    new_img ./= scale

    heatmap(new_img, aspect_ratio=1.0, size=(resolution*3/2, resolution), 
    clim=(0,3)
    )
    # contour(new_img, aspect_ratio=1.0, size=(resolution*3/2, resolution), clim=(0,3))
    title!("Temperature Scale = $scalestr, Mass = $mass M_☼")
end

temperature_render(obs_angle=85.0, mass=1)
>>>>>>> refs/remotes/origin/main
