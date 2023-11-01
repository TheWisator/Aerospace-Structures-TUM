using Plots, Roots

gr()

G = 10^7
K_a = 0.7
A_ant = 70^2 * π / 4
a = 0.2
P_emit = 450_000

distance_range = (10000:100:100_000_000) .* 1e3
size_range = 10:10:1000



P_r(d_ast, r_ast) = P_emit .* G .* A_ant .* K_a .* a .* π ./ 4 .* d_ast .^ 2 ./ ((4 .* π .* r_ast .^ 2).^2)

R_(d_ast) = (P_emit .* G .* A_ant .* K_a .* a .* π ./ 4 .* d_ast .^ 2 ./ ((4 .* π).^2) ./ 1e-18) .^(1/4)

P_1 = P_r(800, distance_range)
R1 = R_(size_range)
fun1(x) = exp.(-x) .* sin.(100x) .+ x .* exp.(-0.1x)
xrange = range(0.001,100,length=10000)


#p1 = plot(xrange, fun1(xrange), xlabel="x", ylabel = "f(x)")
r = 6178000 + 800000

Lrange = 0:0.1:10
F(L) = 0.001 .- 3.986e14 .* (50 ./ ((r .- 0.923 .* L).^2 ) .- 50 ./ ((r .- 0.923 .* L) .* r ) .- 600 ./ ((r .+ 0.076 .* L).^2 ) .+ 600 ./ ((r .+ 0.076 .* L) .* r ))
F2(L) = 3.986e14 .* (50 ./ ((r .- 0.923 .* L).^2 ) .- 50 ./ ((r .- 0.923 .* L) .* r ) .- 600 ./ ((r .+ 0.076 .* L).^2 ) .+ 600 ./ ((r .+ 0.076 .* L) .* r ))

print(find_zero(F, 0))
p1 = plot(Lrange, F(Lrange), xlabel="L [m]", ylabel="delta F [N]", label="Function")
vline!(p1, [9.29], label="Root")
