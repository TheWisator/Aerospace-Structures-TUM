using Plots, LinearAlgebra, Unitful
#plotly()

## Material Definition of Ti6Al4V
E_ti = 114000
ν_ti = 0.33
R_p02_ti = 830

## Definition Case 1: Orbiter / B747
A_c1 = 621.1
D_c1 = 28.1
L_1 = 1796.6
k_c1 = (E_ti * A_c1 / L_1)
Fx_c1 = 90300
Fy_c1 = 154800
ϕ_c1 = 52.8 * π / 180

## Definition Case 2: Orbiter / B747
A_c2 = 2410.9
D_c2 = 55.4
L_2 = 1796.6
k_c2 = (E_ti * A_c2 / L_2)
Fx_c2 = 404800
Fy_c2 = 529300
ϕ_c2 = 52.8 * π / 180




## Case 1
# define Master Stiffness Matrix
K(p1, p2, k1, k2) = ustrip.([   k1 .* cos.(p1).^2 -k1 .* sin.(p1) .* cos.(p1) -k1 .* cos.(p1).^2 k1 .* sin.(p1) .* cos.(p1) 0 0;
                -k1 .* sin.(p1) .* cos.(p1) k1 .* sin.(p1).^2 k1 .* sin.(p1) .* cos.(p1) -k1 .* sin.(p1).^2 0 0;
                -k1 .* cos.(p1).^2 k1 .* sin.(p1) .* cos.(p1) k1 .* cos.(p1).^2 .+ k2 .* cos.(p2).^2 -k1 .* sin.(p1) .* cos.(p1) + k2 .* sin.(p2) .* cos.(p2) -k2 .* cos.(p2).^2 -k2.* sin.(p2) .* cos.(p2);
                k1 .* sin.(p1) .* cos.(p1) -k1 .* sin.(p1).^2 -k1 .* sin.(p1) .* cos.(p1) .+ k2 .* sin.(p2) .* cos.(p2) k1 .* sin.(p1).^2 .+ k2 .* sin.(p2).^2 -k2 .* sin.(p2) .* cos.(p2) -k2 .* sin.(p2).^2;
                0 0 -k2 .* cos.(p2).^2 -k2 .* sin.(p2) .* cos.(p2) k2 .* cos.(p2).^2 k2 .* sin.(p2) .* cos.(p2);
                0 0 -k2 .* sin.(p2) .* cos.(p2) -k2 .* sin.(p2).^2 k2 .* sin.(p2) .* cos.(p2) k2 .* sin.(p2).^2]) .* unit(k_c1)


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

function solve(ϕ1, ϕ2, Fx, Fy, k1, k2, A1, A2)
    # 6x6 stiffness matrix according to sheet
    K_mat = K(ϕ1, ϕ2, k1, k2)
    # 6x1 displacement vector
    F_vec = F(Fx, Fy)
    # solve stiffness equation for two free displacements
    u_free = inv(K_mat[3:4, 3:4]) * F_vec[3:4]
    # calculate reaction Forces
    F_r = K_mat * u(u_free[1], u_free[2])
    # rotation matrix
    T(α) = [    cos(α) -sin(α);
                sin(α) cos(α)]
    
    # stresses in element 1
    u1 = T((-ϕ1)) * u_free
    σ1 = u1[1] * k1 / A1

    # stresses in element 2
    u2 = T((ϕ2)) * u_free
    σ2 = u2[1] * k2 / A2

    σ_r = [σ1, σ2]

    return (u_free, F_r, σ_r)
end

function iterate(ϕ1, ϕ2, Fx, Fy, E)
    A_iter = 100
    k_iter = E * A_iter / L_1
    σ_req = R_p02_ti / 3
    solved_iter = solve(ϕ1, ϕ2, Fx, Fy, k_iter, k_iter, A_iter, A_iter)
    σ_cache = maximum(solved_iter[3])
    residual = σ_req - σ_cache



    while abs(residual) > 0.001
        A_iter = A_iter * (1 - residual / σ_req)
        k_iter = E * A_iter / L_1
        solved_iter = solve(ϕ1, ϕ2, Fx, Fy, k_iter, k_iter, A_iter, A_iter)
        σ_cache = maximum(solved_iter[3])
        residual = σ_req - σ_cache
    end
    return A_iter
end



## Case 1
solved = solve(ϕ_c1,ϕ_c1, Fx_c1, Fy_c1, k_c1, k_c1, A_c1, A_c1)

println("CASE 1: B747")
println("Solved Area Case 1 [mm²]:")
println(iterate(ϕ_c1,ϕ_c1, Fx_c1, Fy_c1, E_ti))
println("Displacement Vector:")
println(u(solved[1][1], solved[1][2]))
println("Force Vector:")
println(solved[2])
println("Stresses")
println(solved[3])

## Case 2
solved = solve(ϕ_c1,ϕ_c1, Fx_c2, Fy_c2, k_c2, k_c2, A_c2, A_c2)

println("CASE 2: Tank")
println("Solved Area Case 2 [mm²]:")
println(iterate(ϕ_c2, ϕ_c2, Fx_c2, Fy_c2, E_ti))
println("Displacement Vector:")
println(u(solved[1][1], solved[1][2]))
println("Force Vector:")
println(solved[2])
println("Stresses")
println(solved[3])


#= TASK 4
ϕ_series = range(0.1, deg2rad(89), length=100)

# result vactors
u_y = []
detK = []
eigenvalK = []
eigenvecK = []

for ϕ in ϕ_series
    solved_cache = solve_task_three(ϕ, ϕ, Fx_c1, Fy_c1)
    K_matrix = K(ϕ, ϕ)
    append!(u_y, solved_cache[1][2])
    append!(detK, det(K_matrix))
    append!(eigenvalK, [round.(eigvals(K_matrix), digits=2)])
    append!(eigenvecK, [round.(eigvecs(K_matrix), digits=2)])
end

p1 = plot(rad2deg.(ϕ_series), u_y)

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
=#