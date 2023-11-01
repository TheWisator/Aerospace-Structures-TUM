mutable struct Body
    # struct for a body with initial conditions
    # position in m
    posx::Float64
    posy::Float64
    posz::Float64
    
    # velocity in m/s
    velx::Float64
    vely::Float64
    velz::Float64
    
    # acceleration in m/s²
    accx::Float64
    accy::Float64
    accz::Float64
    
    # mass in kg
    mass::Float64
    # density in kg/m³
    density::Float64
end

const G0 = 6.6743e-11               # gravitational constant in m³/kg/s²

# definition of bodies
sun = Body(0,0,0,0,0,0,0,0,0,1.98847e30,1410)
earth = Body(146_900_000_000,0,0,0,sqrt(G0 * sun.mass / 146_900_000_000),0,0,0,0,5.9722e24,5541)

bodies = [sun, earth]               # vector holding all bodies
const N = length(bodies)            # number of bodies


function assemble_matrices(bodies)
    positions =  zeros(N,3)
    velocities = zeros(N,3)
    mass = zeros(N)

    for i in 1:1:length(bodies)
        positions[i, :] = [bodies[i].posx, bodies[i].posy, bodies[i].posz]
        velocities[i, :] = [bodies[i].velx, bodies[i].vely, bodies[i].velz]
        mass[i] = bodies[i].mass
    end

    return (positions, velocities, mass)
end

function calculate_acceleration(positions, masses, N)
    accelerations = zeros(N,3)
    for i in 1:1:N
        for j in 1:1:N
            r0 = positions[i,:]
        end
    end
end

print(assemble_matrices(bodies))
