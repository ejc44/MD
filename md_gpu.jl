using CuArrays


# Generate strengths of interactions between particle types
# Very simple but user could replace with anything they wanted

function gen_interaction(num_part_types)
    interaction_params = zeros(Float64, num_part_types,num_part_types)
    for i=1:num_part_types
        for j = 1:num_part_types
            if (i==j) # Self-interaction is randomly repulsive
                interaction_params[i,j] = -1*rand(Float64)
            elseif (i<j) # Others randomly attractive
                val = rand(Float64)
                interaction_params[i,j] = val
                interaction_params[j,i] = val
            end
        end        
    end
    
    return interaction_params
end


# Global constants 
# If convert to script, some should be user inputs

const number_of_steps = 10000; # Number of steps to execute in simulation
const dim = 3; # Dimensions of simulation
const box_size = 10; # Size of one side of box
const finite_box = true; # Whether box if finite or just where particles are initially placed
const periodic = true; # Whether simulation is periodic
const part_num = 100; # Number of particles in simulation
const dt = 0.01; # Time step
const num_part_types = 2; # Number of types of particles
const interaction_params = gen_interaction(num_part_types); # Interations parameters for types of particles
const mass = rand(Float64); # Masses of types of particles
const save_interval = 10; # How often save position
const save_data = false; # Whether to save data to file so can restart


# Initialize position, velocity, acceleration, and particle type
# Currently start with zero velocity and acceleration
# Currently start with random positions and randomly assigned particle types

function initialize(part_num, dim, box_size, num_part_types)
    
    pos = cu(box_size*rand(Float64, part_num, dim)) # Initialized to be randomly placed within a box
    
    vel = cu(zeros(Float64, part_num, dim)) # Initialized to zero
    acc = cu(zeros(Float64, part_num, dim)) # Initialized to zero
    
    part_types = cu(rand(1:num_part_types, part_num)) # Randomly assign type of each particle
        
    return pos, vel, acc, part_types
end


# Update position, velocity, and acceleration using Velocity Verlet Algorithm
# Currently system size is infinite (positions not constrained)

function step_update(pos, vel, acc, force, mass, dt)
    f_m  = cu(similar(force))
    pos .= pos .+ vel .* dt .+ acc .* 0.5 .* (dt^2)
    f_m .= force ./ mass
    vel .= vel .+ acc .* f_m .* 0.5 .* dt
    acc .= f_m  
    
    return pos, vel, acc
end


# Find force on each particle
# 1/r^2 interactions currently
# Again, very simple but user can replace with anything they want

function find_force(part_num, dim, pos, vel, acc, part_types, interaction_params, mass)
    
    force = cu(zeros(part_num, dim))
    
    # Loop is parallelizable
    for i = 1:part_num # For every particle
        for k = 1:part_num # Contribution from every other particle
            if (i != k) # No self-interaction
                for j = 1:dim # For each dimension
                    int_strength = interaction_params[part_types[i],part_types[k]] # Strength of interaction between particles
                    if (pos[i,j] > pos[k,j])
                       int_strength = -1*int_strength # Reverses direction of force if positions flopped
                    end

                    # Find distance between particles
                    dist = 0.0
                    if (periodic) # If periodic, check whether a periodic distance is shorter
                        dist = abs(pos[i,j] - pos[k,j])
                        if (pos[i,j] < pos[k,j])
                            new_dist = abs(box_size + pos[i,j] - pos[k,j])
                            if new_dist > dist
                                dist = new_dist
                                int_strength = -1*int_strength # Reverses direction of force
                            end
                        else
                            new_dist = abs(box_size + pos[k,j] - pos[i,j])
                            if new_dist > dist
                                dist = new_dist
                                int_strength = -1*int_strength # Reverses direction of force
                            end
                        end
                        else # Otherwise, just regular distance
                        dist = abs(pos[i,j] - pos[k,j])
                    end
                    
                    force[i,j] += int_strength / dist^2 # 1/r^2 interaction
                end
            end
        end
    end
    
    return force
end


# Write each variable to its own output file in current diretory
# May need to add more imfo if want to restart

function write_output(part_num, dim, part_types, interaction_params, mass, dt, saved_positions, saved_velocities, saved_accelerations, saved_forces) 
    writedlm("part_num.txt", part_num)
    writedlm("dim.txt", dim)
    writedlm("part_types.txt", part_types)
    writedlm("interaction_params.txt", interaction_params)
    writedlm("mass.txt", mass)
    writedlm("dt.txt", dt)
    writedlm("saved_positions.txt", saved_positions)
    writedlm("saved_velocities.txt", saved_velocities)
    writedlm("saved_accelerations.txt", saved_accelerations)
    writedlm("saved_forces.txt", saved_forces)
end


# Main - where program executes
# Convert to main(args) if convert to script
# This is the part of the program that is timed

function main()
    
    # Set up matrices to save positions, velocities, and accelerations into
    num_pos = floor(Int, number_of_steps/save_interval)+1 # Number of positions to save
    saved_positions = zeros(num_pos, part_num, dim) # Save positions with specified frequency
    if (save_data) # Only save vel, acc, and force if plan on saving to file
        saved_velocities = zeros(num_pos, part_num, dim) # Save velocities with specified frequency
        saved_accelerations = zeros(num_pos, part_num, dim) # Save accelerations with specified frequency
        saved_forces = zeros(num_pos, part_num, dim) # Save forces with specified frequency
    end
    save_index = 1 # Keep track of index so can save positions with specified frequency

    pos, vel, acc, part_types = initialize(part_num, dim, box_size, num_part_types) # Initialize
    force = find_force(part_num, dim, pos, vel, acc, part_types, interaction_params, mass) # Find forces on particles
    saved_positions[save_index,:,:] = pos # Save position
    if (save_data) # Only save vel, acc, and force if plan on saving to file
        saved_velocities[save_index,:,:] = vel # Save velocity
        saved_accelerations[save_index,:,:] = acc # Save accelerations
        saved_forces[save_index,:,:] = force # Save forces
    end
    save_index += 1 # Increment index

    for i = 1:number_of_steps
        pos, vel, acc = step_update(pos, vel, acc, force, mass, dt) # Update
        force = find_force(part_num, dim, pos, vel, acc, part_types, interaction_params, mass) # Find new forces

        if (i % save_interval == 0) # If are saving this time step
            saved_positions[save_index,:,:] = pos # Save position
            if (save_data) # Only save vel, acc, and force if plan on saving to file
                saved_velocities[save_index,:,:] = vel # Save velocity
                saved_accelerations[save_index,:,:] = acc # Save accelerations
                saved_forces[save_index,:,:] = force # Save forces
            end
            save_index += 1 # Increment index
        end
    end
    
    # Don't save data in benchmarking version
    #=if (save_data) # If save data
        write_output(part_num, dim, part_types, interaction_params, mass, dt, saved_positions, saved_velocities, saved_accelerations, saved_forces)
    end=#
    
end


using BenchmarkTools

println(@benchmark main())
