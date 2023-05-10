using Plots, LinearAlgebra
plotly()



## definition of important constants
k1 = 50
k2 = 75

## TASK 1 and TASK 2
# define Master Stiffness Matrix
K(p1, p2) = [   k1 .* cos.(p1).^2 -k1 .* sin.(p1) .* cos.(p1) -k1 .* cos.(p1).^2 k1 .* sin.(p1) .* cos.(p1) 0 0;
                -k1 .* sin.(p1) .* cos.(p1) k1 .* sin.(p1).^2 k1 .* sin.(p1) .* cos.(p1) -k1 .* sin.(p1).^2 0 0;
                -k1 .* cos.(p1).^2 k1 .* sin.(p1) .* cos.(p1) k1 .* cos.(p1).^2 .+ k2 .* cos.(p2).^2 -k1 .* sin.(p1) .* cos.(p1) + k2 .* sin.(p2) .* cos.(p2) -k2 .* cos.(p2).^2 -k2.* sin.(p2) .* cos.(p2);
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

function solve_task_three(ϕ1, ϕ2, Fx, Fy)
    # 6x6 stiffness matrix according to sheet
    K_mat = K(ϕ1, ϕ2)
    # 6x1 displacement vector
    F_vec = F(Fx, Fy)
    # solve stiffness equation for two free displacements
    u_free = inv(K_mat[3:4, 3:4]) * F_vec[3:4]
    # calculate reaction Forces
    F_r = K_mat * u(u_free[1], u_free[2])

    return (u_free, F_r)
end

## TASK 3
solved = solve_task_three(π/4,π/4,1,-0.8)

println("Displacement Vector:")
println(u(solved[1][1], solved[1][2]))
println("Force Vector:")
println(solved[2])

## TASK 4
ϕ_series = range(0.1, deg2rad(89), length=100)

# result vactors
u_x = []
u_y = []
detK = []
eigenvalK = []
eigenvecK = []

for ϕ in ϕ_series
    solved_cache = solve_task_three(ϕ, ϕ, 1, -0.8)
    K_matrix = K(ϕ, ϕ)
    append!(u_y, solved_cache[1][2])
    append!(detK, det(K_matrix))
    append!(eigenvalK, [round.(eigvals(K_matrix), digits=2)])
    append!(eigenvecK, [round.(eigvecs(K_matrix), digits=2)])
end

p1 = plot(rad2deg.(ϕ_series), sqrt.(u_y.^2 .+ u_x.^2))

xlabel!(p1, "ϕ [degrees]")
ylabel!(p1, "u_y")
#=println("Eigenvalues of K")
display(eigenvalK)
println("Eigenvectors of K")
display(eigenvecK)=#

display(p1)

## ADD-ON: Alternative way to find K
function solve_task_alternatively(ϕ1, ϕ2, Fx, Fy)
    # rotation matrix
    T(α) = [    cos(α) -sin(α) 0 0;
                sin(α) cos(α) 0 0;
                0 0 cos(α) -sin(α);
                0 0 sin(α) cos(α)]

    ## Element 1
    K_1 = k1 * [1 0 -1 0;
                0 0 0 0;
                -1 0 1 0;
                0 0 0 0;] 
    T1 = T(ϕ1)
    # 4x4 in global sytem
    K_1e = transpose(T1) * K_1 * T1
    # extend to 6x6 in global system
    K_1e = [    K_1e[1,1] K_1e[1,2] K_1e[1,3] K_1e[1,4] 0 0;
                K_1e[2,1] K_1e[2,2] K_1e[2,3] K_1e[2,4] 0 0;
                K_1e[3,1] K_1e[3,2] K_1e[3,3] K_1e[3,4] 0 0;
                K_1e[4,1] K_1e[4,2] K_1e[4,3] K_1e[4,4] 0 0;
                0 0 0 0 0 0;
                0 0 0 0 0 0] 

    ## Element 2
    K_2 = k2 * [1 0 -1 0;
                0 0 0 0;
                -1 0 1 0;
                0 0 0 0;]  
    T2 = T(-ϕ2)
    #4x4 in global system
    K_2e = transpose(T2) * K_2 * T2
    # extend to 6x6 in global system
    K_2e = [    0 0 0 0 0 0;
                0 0 0 0 0 0;
                0 0 K_2e[1,1] K_2e[1,2] K_2e[1,3] K_2e[1,4];
                0 0 K_2e[2,1] K_2e[2,2] K_2e[2,3] K_2e[2,4];
                0 0 K_2e[3,1] K_2e[3,2] K_2e[3,3] K_2e[3,4];
                0 0 K_2e[4,1] K_2e[4,2] K_2e[4,3] K_2e[4,4]]

    ## Global Evaluation
    # global 6x6 stiffness matrix
    K_ges = K_1e + K_2e
    # 6x1 displacement vector
    F_vec = F(Fx, Fy)
    # solve stiffness equation for two free displacements
    u_free = inv(K_ges[3:4, 3:4]) * F_vec[3:4]
    # calculate reaction Forces
    F_r = K_ges * u(u_free[1], u_free[2])

    return (u_free, F_r)
end

solved = solve_task_alternatively(π/4,π/4,1,-0.8)

println("Displacement Vector alternatively:")
println(u(solved[1][1], solved[1][2]))
println("Force Vector alternatively:")
println(solved[2])