using GeodesicRendering
using AccretionGeometry
using CarterBoyerLindquist
using AccretionFormulae
using StaticArrays
using Printf
<<<<<<< HEAD
using LaTeXStrings
=======
>>>>>>> 9dfda8835958162ae46e0e575201462dd749fd60

using Plots
gr()

<<<<<<< HEAD
function mass_scale_fraction(M, η, edd_ratio, r_isco, r_g, r)
    mass_scale_fraction = ((M * η)^(-1) * (edd_ratio * (1-√(r_isco/r))) * (r_g/r)^3)^(1/4)
end

function temperature(m, sol, max_time; kwargs...)
    g =  AccretionFormulae.redshift(m, sol, max_time; kwargs...)

    u = sol.u[end]
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
=======
function temperature(m, gp, max_time; kwargs...)
    g = AccretionFormulae.redshift(m, gp, max_time; kwargs...)

    M = m.M
    a_star = m.a
    temperature = AccretionFormulae.observed_temperature(gp.u[2], a_star, M, g)

end

function radius(m, gp, max_time; kwargs...)
    gp.u[2]
end

function temperature_render(;
    mass = 1,
    spin = 0.998,
    obs_angle = 85.0,
    disc_angle = 90.0,
    tolerance = 1e-8,
    dtmax = 1000.0,
    size_multiplier::Int64 = 1,
    resolution = 400,
    fov = 3.0,
    η = 0.1,
    η_phys = 0.1,
    edd_ratio = 0.1,
    edd_ratio_phys = 0.1
)

    m = CarterMethodBL(M = 1.0, a = spin)
>>>>>>> 9dfda8835958162ae46e0e575201462dd749fd60
    M = m.M

    # observer position
    u = [0.0, 1000.0, deg2rad(obs_angle), 0.0]
    R_isco = AccretionFormulae.r_isco(m.a, m.M)

    # disc
<<<<<<< HEAD
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
        ∘ ConstValueFunctions.filter_early_term
    )

    # do the render
    temperature_img = @time rendergeodesics(
        m, u, 2000.0, 
        d, 
        fov_factor=fov*size_multiplier, abstol=tolerance, reltol=tolerance,
        image_width = 350*size_multiplier,
        image_height = 250*size_multiplier,
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

    numerators = mass_scale_fraction.(M_phys, η_phys, edd_ratio_phys, r_isco_phys, r_g_phys, radius_img*M_phys)
    denominators = mass_scale_fraction.(M, η, edd_ratio, R_isco, M, radius_img)

    fractions = numerators./denominators
=======
    d = GeometricThinDisc(R_isco + 1, 50.0, deg2rad(disc_angle))


    # cache the render
    cache = @time prerendergeodesics(
        m,
        u,
        2000.0,
        d,
        fov_factor = fov * size_multiplier,
        abstol = tolerance,
        reltol = tolerance,
        image_width = 350 * size_multiplier,
        image_height = 350 * size_multiplier,
        dtmax = dtmax,
    )

    # create and compose the PointFunctions
    temperature_vf = (
        PointFunction(temperature) ∘
        FilterPointFunction((m, gp, max_time; kwargs...) -> gp.u[2] > R_isco, NaN) ∘
        ConstPointFunctions.filter_early_term
    )

    radius_vf = (
        PointFunction(radius) ∘
        FilterPointFunction((m, gp, max_time; kwargs...) -> gp.u[2] > R_isco, NaN) ∘
        ConstPointFunctions.filter_early_term
    )

    radius_img = GeodesicRendering.apply(radius_vf, cache)
    temperature_img = GeodesicRendering.apply(temperature_vf, cache)

    # correcting for physical mass
    M_phys = mass * 1.99e30
    r_isco_phys = AccretionFormulae.r_isco(M_phys, 0.998)
    r_g_phys = M_phys

    numerators =
        AccretionFormulae.mass_scale_fraction.(
            M_phys,
            η_phys,
            edd_ratio_phys,
            r_isco_phys,
            r_g_phys,
            radius_img * M_phys,
        )
    denominators =
        AccretionFormulae.mass_scale_fraction.(M, η, edd_ratio, R_isco, M, radius_img)

    fractions = numerators ./ denominators
>>>>>>> 9dfda8835958162ae46e0e575201462dd749fd60
    temperature_img .*= fractions

    # # autoscale
    # scale = maximum(filter(!isnan,temperature_img))
    # scale = floor(log(10, scale))
    # scale = 10^scale
<<<<<<< HEAD
    
    # scaling image
    scale = 1e6
    scalestr = @sprintf "%.E" scale
    new_img = reverse(temperature_img, dims=1)
    new_img ./= scale

    heatmap(new_img, aspect_ratio=1.0, size=(resolution*3/2, resolution), 
    clim=(0,10)
=======

    # scaling image
    scale = 1e7
    scalestr = @sprintf "%.E" scale
    exponent =  Int32(floor(log(10, scale)))
    new_img = reverse(temperature_img, dims = 1)
    new_img ./= scale

    hmap = heatmap(
        new_img,
        aspect_ratio = 1.0,
        size = (resolution * 3 / 2, resolution),
        clim = (0, 3),
>>>>>>> 9dfda8835958162ae46e0e575201462dd749fd60
    )
    # contour(new_img, aspect_ratio=1.0, size=(resolution*3/2, resolution), clim=(0,3))
    # title!("Temperature Scale = $scalestr, Mass = $mass M_☼, Obs Angle = $obs_angle")
    # title = "Temperature Scale = $scalestr K, Mass = $mass M_☼, Obs Angle = $obs_angle"
    title = "Temperature Scale = \$10^{$exponent}\$ K, Mass = $mass \$\\mathrm{M}_{\\odot}\$, Obs Angle = $obs_angle\$^{\\circ}\$"
    # "")
    return hmap, cache, title
end

<<<<<<< HEAD
# hmap, cache, title = temperature_render(obs_angle = 62.5, mass = 10, resolution=1080)
# # hmap, cache, title = temperature_render(
# #                                         mass=10,
# #                                         spin=0.998,
# #                                         obs_angle=0.01,
# #                                         tolerance=1e-12,
# #                                         size_multiplier=6,
# #                                         dtmax=1,
# #                                         resolution=2000
# #                                         )
# title!(title)
# display(hmap)
=======
<<<<<<< HEAD
temperature_render(obs_angle=85.0, mass=100)
=======
temperature_render(obs_angle = 85.0, mass = 1)
>>>>>>> 9dfda8835958162ae46e0e575201462dd749fd60
>>>>>>> c815900e2b495adc3d734f141dbbeaf9c5d28fc5
