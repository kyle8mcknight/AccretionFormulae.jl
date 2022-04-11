include("emissivity_deviations.jl")

"""
From Reynolds (2020) (Eq2-4)
Calculates the radius of the inner most stable circular orbit of the a black hole,
with spin a_star, and mass M.
"""
function r_isco(a_star, M)
    a = a_star
    r_g = M
    z1 = 1 + ∛(1 - a^2) * (∛(1 + a) + ∛(1 - a))
    z2 = √(3 * a^2 + z1^2)
    if a >= 0
        r_isco = (3 + z2 - √((3 - z1) * (3 + z1 + 2 * z2))) * r_g
    elseif a < 0
        r_isco = (3 + z2 + √((3 - z1) * (3 + z1 + 2 * z2))) * r_g
    end
end

r_vals, f_vals = f_arbitrary()

function mdot(M)
    # need something here to convert M into SI units, as M_☼ is in SI, althought we 
    # aren't currently clear about what the units M currently is in
    # would look something like:
    # M_kg = M*c^2
    # L_edd = 3e4*L_☼*(M_kg/M_☼)

    L_edd = 3e4 * L_☼ * (M / M_☼)
    Mdot = -L_edd / (c^2 * η)
    mdot = -0.1 * Mdot
end

function diss(mdot, r)
    index = findall(x->x==r, r_vals)
    diss = ((c^6) / (G^2)) * mdot * f_vals[index] / (4 * π * r)
end

function temperature(r, M)
    m_dot = mdot(M)
    temperature = (diss(m_dot, r)[1] / σ_SB)^(1 / 4)
end

"""
From Fanton et al. (1997) (Eq78)
"""
function observed_temperature(r, g)
    T = temperature(r, a_star, M)
    observed_temperature = g * T
end

function mass_scale_fraction(M, η, edd_ratio, r_isco, r_g, r)
    mass_scale_fraction =
        ((M * η)^(-1) * (edd_ratio * (1 - √(r_isco / r))) * (r_g / r)^3)^(1 / 4)
end

# constants
G = 6.67e-11
c = 3e8
L_☼ = 3.8e26
M_☼ = 1.99e30
σ_SB = 5.67e-8
η = 0.1

export observed_temperature, r_isco, temperature, mass_scale_fraction

temperature_vals = temperature.(r_vals, 1.99e31)

plt = plot(r_vals, temperature_vals)
display(plt)