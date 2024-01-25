#script to manually initialize the EKP, used for debugging. Currently based on the idealized model!
using LinearAlgebra
using Random
using Distributions
using EnsembleKalmanProcesses
using EnsembleKalmanProcesses.ParameterDistributions
using JLD2
using EnsembleKalmanProcesses.Localizers
const EKP = EnsembleKalmanProcesses

y = [13.77, 13.37, 12.97, 12.52] #observed SLR
Γ = 1.0 * I #noise level

#setup priors
prior_weertman = constrained_gaussian("weertman_c_prefactor", 1.0, 0.2, -Inf, Inf)
prior_glen_a = constrained_gaussian("glen_a_prefactor", 1.0, 0.2, -Inf, Inf)
prior_bump_amplitude = constrained_gaussian("glen_a_prefactor", 200.0, 100.0, -Inf, Inf)
prior_melt_prefactor= constrained_gaussian("melt_rate_prefactor", 8.0, 2.5, -Inf, Inf)
prior_per_century_trend = constrained_gaussian("per_century_trend", 100.0, 100.0, -Inf, Inf)
prior = combine_distributions([prior_weertman, prior_glen_a,prior_bump_amplitude,  prior_melt_prefactor, prior_per_century_trend]) #order must match that in initialize_EKP.jl

N_ensemble = 10 #number of ensemble members
 

rng_seed = 87653 #random seed, hard coded?
rng = Random.MersenneTwister(rng_seed)
initial_ensemble = EKP.construct_initial_ensemble(rng, prior, N_ensemble)

ensemble_kalman_process = EKP.EnsembleKalmanProcess(initial_ensemble, y, Γ, Inversion(); rng = rng)

params_i = get_ϕ_final(prior, ensemble_kalman_process) #initial parameters used in simulations

#get outputs from first simulations
outputs = zeros(4,10)
for i = 1:9
       fname = "/data/hpcdata/users/aleey/projects/AttribEKF/MyKalmanEnsembling/KalmanEnsembling/ensemble/output/iteration_000/member_00" *  string(i) * "/output.jld2"; dd = load(fname); print(dd); outputs[:,i] = dd["model_output"];
end
for i = 10
       fname = "/data/hpcdata/users/aleey/projects/AttribEKF/MyKalmanEnsembling/KalmanEnsembling/ensemble/output/iteration_000/member_0" *  string(i) * "/output.jld2"; dd = load(fname); print(dd); outputs[:,i] = dd["model_output"];
end

G_ens = outputs #values from the simulation
EKP.update_ensemble!(ensemble_kalman_process, G_ens)
get_ϕ_final(prior, ensemble_kalman_process) #new values of parameters
