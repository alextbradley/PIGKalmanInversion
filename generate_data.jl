include("observe_sinusoid.jl")

using LinearAlgebra, Distributions, JLD2, Random

function main()

    # Paths
    output_dir = joinpath(ARGS[1])
    if !isdir(output_dir)
        mkdir(output_dir)
    end
    data_path = joinpath(output_dir, ARGS[2])

    if isfile(data_path)
        @warn data_path * " already exists"
        exit(0)
    end

    rng_seed = parse(Int64, ARGS[3])
    rng_model = Random.MersenneTwister(rng_seed)

    # the true parameters
    # theta_true = Dict("initial_thickness" => 500.0*ones(4) )

    #randomness

    # create noise
    Γ = 1.0 * I
    #dim_output = length(parameter_to_data_map(theta_true,data_path)) #just to get size here
    dim_output = 7 #must match size of output
    noise_dist = MvNormal(zeros(dim_output), Γ)
    
    observed_volume = [8.0803, 8.0514, 8.0226, 7.9937, 7.9649, 7.9360,  7.9071]; #observed ice volume in 10^14 m^3 at times 1990, 1995, 2000, 2005, 2010, 2015, 2020]. NB: currently based on linear retreat with retreat initiated in the 1940s 

    # evaluate map with noise to create data
    # y = parameter_to_data_map(theta_true,data_path) .+ rand(rng_model, noise_dist) 
    y = observed_volume .+ rand(rng_model, noise_dist)

    # save
    @save data_path y Γ rng_model

end

main()
