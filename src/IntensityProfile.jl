using AccretionFormulae

# constants
h = 6.63e-34
c = 3.00e8
f_col = 1.7
Υ = 1
k_b = 1.38e-23
M_☼ = 1.99e30

function intensity(r, ν_e, a_star, M) # \nu not v
    temperature_eff = AccretionFormulae.temperature(r, a_star, M)
    temperature_col = f_col*temperature_eff

    intensity = (2*h*ν_e^3)/(c^2) * 1/(f_col^4) * Υ/(exp((h*ν_e)/(k_b*temperature_col)) - 1)
end

M = 10*M_☼
a_star = 0.8
r = 6e30

λ_vals = LinRange(0.3e-9, 3e-9, 1000)
ν_vals = c/λ_vals

intensity_vals = intensity.(r, ν_vals, a_star, M)

plot(λ_vals, intensity_vals)

