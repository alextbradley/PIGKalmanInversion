include("shared.jl") # packages

function main()

    output_dir = ARGS[1]
    data_path = joinpath(output_dir, ARGS[2])
    eki_path = joinpath(output_dir, ARGS[3])

    if isfile(eki_path)
        @warn eki_path * " already exists"
        exit(0)
    end

    toml_path = ARGS[4]
    N_ensemble = parse(Int64, ARGS[5])
    rng_seed = parse(Int64, ARGS[6])
    rng_ekp = Random.MersenneTwister(rng_seed)

    # We construct the prior from file. Parameters as follows:
    # Physical parameters
    # weertman_c_prefactor: premultiplier of the weertman c field (weertmanC = weertman_c_prefactor .* observationally_constrained_weertmanC)
    # dimensionless_ungrounded_weertmanC_exponent: value of the N(0,1) distribution (weertman C in ungrounded regions ~ 10^(\mu + \sigma^2 * dimensionless_ungrounded_weermanC))
    # glen_a_ref_prefactor: premultiplier of the glen a field (glenA = glen_a_ref_prefactor .* observationally_constrained_glenA)
    # melt_rate_prefactor: dimensionless premultiplier of the melt rate field
    #
    # Forcing parameters
    # per_century_trend: per century shoaling of the pycnocline
    # bump_amplitude: amplitude of the 1940s bump 
    # bump_duration: duraction of the 1940s anomaly
    #
    param_dict = TOML.parsefile(toml_path)
    #names = ["weertman_c_prefactor", "glen_a_ref_prefactor", "bump_amplitude", "melt_rate_prefactor", "per_century_trend"]
    names = ["weertman_c_prefactor", "dimensionless_ungrounded_weertmanC_exponent", "glen_a_ref_prefactor",  "melt_rate_prefactor", "per_century_trend","bump_amplitude","bump_duration"]

    
    prior_vec = [get_parameter_distribution(param_dict, n) for n in names]
    prior = combine_distributions(prior_vec)

    # we load the data, noise, from file
    @load data_path y Γ

    # initialize ensemble Kalman inversion
    initial_ensemble = EKP.construct_initial_ensemble(rng_ekp, prior, N_ensemble)
    #eki = EKP.EnsembleKalmanProcess(initial_ensemble, y, Γ, Inversion(); rng = rng_ekp, localization_method = SEC(3)) #apply some localization
    eki = EKP.EnsembleKalmanProcess(initial_ensemble, y, Γ, Inversion(); rng = rng_ekp)

    # save the parameter ensemble and EKP
    save_parameter_ensemble(
        get_u_final(eki), # constraints applied when saving
        prior,
        param_dict,
        output_dir,
        "parameters.toml",
        0, # We consider the initial ensemble to be the 0th iteration
    )

    #save new state
    @save eki_path eki param_dict prior


end

main()
