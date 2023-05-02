using Plots
gr()

## definition of important constants
k1 = 50
k2 = 75

## TASK 1 and TASK 2
# define Master Stiffness Matrix
K(p1, p2) = [   k1 .* cos.(p1).^2 -k1 .* sin.(p1) .* cos.(p1) -k1 .* cos.(p1).^2 k1 .* sin.(p1) .* cos.(p1) 0 0;
                -k1 .* sin.(p1) .* cos.(p1) -k1 .* sin.(p1).^2 k1 .* sin.(p1) .* cos.(p1) -k1 .* sin.(p1).^2 0 0;
                k1 .* cos.(p1).^2 k1 .* sin.(p1) .* cos.(p1) k1 .* cos.(p1).^2 .+ k2 .* cos.(p2).^2 -k1 .* sin.(p1) .* cos.(p1) + k2 .* sin.(p2) .* cos.(p2) -k2 .* cos.(p2).^2 -k2.* sin.(p2) .* cos.(p2);
                k1 .* sin.(p1) .* cos.(p1) -k1 .* sin.(p1).^2 -k1 .* sin.(p1) .* cos.(p1) .+ k2 .* sin.(p2) .* cos.(p2) k1 .* sin.(p1).^2 .+ k2 .* sin.(p2).^2 -k2 .* sin.(p2) .* cos.(p2) -k2 .* sin.(p2).^2;
                0 0 -k2 .* cos.(p2).^2 -k2 .* sin.(p2) .* cos.(p2) k2 .* cos.(p2).^2 k2 .* sin.(p2) .* cos.(p2);
                0 0 -k2 .* sin.(p2) .* cos.(p2) -k2 .* sin.(p2).^2 k2 .* sin.(p2) .* cos.(p2) k2 .* sin.(p2).^2]


# define External Forces vector
F(fx, fy) = [   0;
                0;
                fx;
                fy;
                0;
                0]

# define general Displacement vector
u(u3, u4)   =  [0;
                0;
                u3;
                u4;
                0;
                0]



## TASK 3
# solve stiffness equation
u_solve = inv(K(π/4,π/4)[3:4, 3:4]) * F(1,-0.8)[3:4]

# solve for Forces
f = K(π/4,π/4) * u(u_solve[1], u_solve[2])

println("Displacement Vector:")
println(u(u_solve[1], u_solve[2]))
println("Force Vector:")
println(f)

