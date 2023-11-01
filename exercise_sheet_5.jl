using Unitful

w = 150u"mm"
t_cap = 12u"mm"
I = w * t_cap^3 / 12
A = t_cap * w
m = 800u"kg"
v = 59.16u"m/s"
kg = 0.773
c_lα = 0.95 * 2π
ρ_0 = 1.225u"kg/m^3"
w_g = 15.2u"m/s"
S = 13u"m^2"
g = 9.81u"m/s^2"
n_z = 1 + kg * ρ_0 * w_g * S * c_lα / 2 / m / g * v
L = m * g * n_z .|> u"N"
a = 4.5u"m"

m_wing = 100u"kg"

Q_wing = L/2 - m_wing * g * n_z
M_wing = a/2 * Q

σ(H) = M_wing .* H ./ 2 ./ (I .+ A ./ 4 .* (H .- t_cap) .^ 2)
τ(t_web, H) = Q_wing ./ (2 .* t_web .* (H .- 2t_cap))

σ_cfrp = 500u"MPa"
τ_max = 100u"MPa"

RF_bend = σ_cfrp / σ(191u"mm") .|> u"m/m"
RF_shear = τ_max / τ(1u"mm", 191u"mm") .|> u"m/m"

println(L)
println(Q_wing)
println(M_wing)
println(RF_bend)
println(RF_shear)
