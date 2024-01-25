using WAVI, Interpolations, JLD2

# The function parameter_to_data_map takes in the ensemble parameters and returns the observation associated with these parameters. In this case, it's the sea level rise at the end of the simulation
function parameter_to_data_map(ensemble_parameters, member_path)
    driver(ensemble_parameters, member_path)

    #specify the paths of the output files. 
    # This is hard coded based on a timestep of 0.05, and simulation start of 1700. 
    # 1990: timestep 5800
    # 1995: timestep 5900
    # 2000: timestep 6000
    # 2005: timestep 6100
    # 2010: timestep 6200
    # 2015: timestep 6300
    # 2020: timestep 6400
    paths = ["outfile0000005800.mat", "outfile0000005900.mat", "outfile0000006000.mat", "outfile0000006100.mat", "outfile0000006200.mat", "outfile0000006300.mat", "outfile0000006400.mat"]; #these must match output times
    ice_mass = zeros(1,length(paths)) 
    dx = 2.0e3
    dy = 2.0e3
    
    for i = 1:length(paths)
       data = load(member_path*"/"*paths[i])
       ice_mass[i] = sum(sum(data["h"])) .* dx .* dy/1e14
    end 
 
    return ice_mass
end

function driver(ensemble_parameters, member_path)

#
#Grid and boundary conditions
#
nx = 398
ny = 250
nσ = 12
x0 = -1789000.0
y0 = -397000.0
dx = 2000.0
dy = 2000.0
h_mask=trues(nx,ny)

#
# mask and boundary conditions
#
h_mask=Array{Float64}(undef,nx,ny);
read!(joinpath(dirname(@__FILE__),"input_data/Inverse_2km_h_mask_clip_BedmachineV3_FULL_stripe_fix.bin"),h_mask)
h_mask.=ntoh.(h_mask)

u_iszero=Array{Float64}(undef,nx+1,ny);
read!(joinpath(dirname(@__FILE__),"input_data/Inverse_2km_uiszero_clip_BedmachineV3_FULL_stripe_fix.bin"),u_iszero)
u_iszero.=ntoh.(u_iszero)

v_iszero=Array{Float64}(undef,nx,ny+1);
read!(joinpath(dirname(@__FILE__),"input_data/Inverse_2km_viszero_clip_BedmachineV3_FULL_stripe_fix.bin"),v_iszero)
v_iszero.=ntoh.(v_iszero)

sigma_grid=Array{Float64}(undef,nσ);
read!(joinpath(dirname(@__FILE__),"input_data/Inverse_2km_sigma_grid_BedmachineV3_FULL_stripe_fix.bin"),sigma_grid)
sigma_grid.=ntoh.(sigma_grid)

#
# make grid
#
grid = Grid(nx = nx,
            ny = ny,
            nσ = nσ,
            x0 = x0,
            y0 = y0,
            dx = dx,
            dy = dy,
            h_mask = h_mask,
            u_iszero = u_iszero,
            v_iszero = v_iszero,
            σ = sigma_grid)

#
#Bed
#
bathy=Array{Float64}(undef,nx,ny);
read!(joinpath(dirname(@__FILE__), "input_data/Inverse_2km_bed_clip_noNan_BedmachineV3_FULL_stripe_fix.bin"),bathy)
bathy.=ntoh.(bathy)

#
# initial conditions
#
h=Array{Float64}(undef,nx,ny);
read!(joinpath(dirname(@__FILE__), "input_data/steady_thickness_nomelt.bin"),h)
h.=ntoh.(h)

viscosity=Array{Float64}(undef,nx,ny,nσ);
read!(joinpath(dirname(@__FILE__),"input_data/Inverse_2km_viscosity3D_clip_noNan_BedmachineV3_FULL_stripe_fix.bin"),viscosity)
viscosity.=ntoh.(viscosity)

temp=Array{Float64}(undef,nx,ny,nσ);
read!(joinpath(dirname(@__FILE__),"input_data/Inverse_2km_3Dtemp_clip_noNan_BedmachineV3_FULL_stripe_fix.bin"),temp)
temp.=ntoh.(temp)

damage=Array{Float64}(undef,nx,ny,nσ);
read!(joinpath(dirname(@__FILE__),"input_data/Inverse_2km_damage3D_clip_noNan_BedmachineV3_FULL_stripe_fix.bin"),damage)
damage.=ntoh.(damage)


dhdt=Array{Float64}(undef,nx,ny);
read!(joinpath(dirname(@__FILE__),"input_data/Inverse_2km_dhdt_clip_noNan_BedmachineV3_FULL_stripe_fix.bin"),dhdt)
dhdt.=ntoh.(dhdt)

initial_conditions = InitialConditions(initial_thickness = h,
                                        initial_viscosity = viscosity,
                                        initial_temperature = temp,
                                        initial_damage = damage)


#
#solver parameters
#
maxiter_picard = 6
#parallel_spec = SharedMemorySpec(ngridsx=2,ngridsy=1,overlap=1,damping=0.0,niterations=1)
parallel_spec = BasicParallelSpec()
tol_picard = 1.0e-4
solver_params = SolverParams(maxiter_picard = maxiter_picard,tol_picard = tol_picard)

#
#Fixed physical parameters
#
accumulation_rate=Array{Float64}(undef,nx,ny);
read!(joinpath(dirname(@__FILE__),"input_data/Inverse_2km_accumulation_clip_noNan_BedmachineV3_FULL_stripe_fix.bin"),accumulation_rate)
accumulation_rate.=ntoh.(accumulation_rate)

#
# variable physical parameters
#
weertman_c=Array{Float64}(undef,nx,ny);
read!(joinpath(dirname(@__FILE__),"input_data/Inverse_2km_WeertmanC_clip_adjusted_noNan_BedmachineV3_FULL_stripe_fix.bin"),weertman_c)
weertman_c.=ntoh.(weertman_c)
weertman_c_prefactor = ensemble_parameters["weertman_c_prefactor"]
weertman_c_prefactor = max(weertman_c_prefactor, 1e-3) #avoid negative weertman c
weertman_c = weertman_c .* weertman_c_prefactor

sec_per_year = 365.25*24*60^2
glen_a_ref= 4.9e-16 *sec_per_year * 1.0e-9 #standard value used in WAVI 
glen_a_ref_prefactor = ensemble_parameters["glen_a_ref_prefactor"]

glen_a_ref = 10^(glen_a_ref_prefactor) .* glen_a_ref


# adjust the weertman c in floating areas
dimensionless_ungrounded_weertmanC_exponent = ensemble_parameters["dimensionless_ungrounded_weertmanC_exponent"]
mu_wc = 4.3
sigma_wc = 0.5
weertman_c[weertman_c .== 10000] .= 10 .^ (mu_wc + sigma_wc * dimensionless_ungrounded_weertmanC_exponent);

params = Params(accumulation_rate = accumulation_rate,
                                  glen_a_ref = glen_a_ref,
                                  weertman_c = weertman_c)



#
# Melt rate model
#

# variable components
bump_amplitude      = ensemble_parameters["bump_amplitude"]
bump_duration       = ensemble_parameters["bump_duration"]
bump_duration       = max(bump_duration, 1) #must be positive
melt_rate_prefactor = ensemble_parameters["melt_rate_prefactor"]
melt_rate_prefactor = max(melt_rate_prefactor, 1e-3) #must be positive
per_century_trend   = ensemble_parameters["per_century_trend"]
random_seed         = 3

# fixed comonents
end_time = 320. # end of the simulation is 2020
bump_width = bump_duration/2 #standard deviation of the duration
bump_time = 245. #start simulation corresponds to 1700
trend_onset = 260. #1960
pc_max = -400.0
pc_min = -600.0 #pyclocline center limits within internal variability only
pw     = 400.0
rf_threshold = 2.0 #random forcing threshold for pc max/min

idealized_anthro_melt_rate = IdealizedAnthroMeltRate(bump_amplitude = bump_amplitude,
                    bump_width = bump_width,
                    bump_time = bump_time,
                    per_century_trend = per_century_trend,
                    trend_onset = trend_onset,
                    pc_max = pc_max,
                    pc_min = pc_min,
                    M = 2.6*melt_rate_prefactor, #melt_rate_prefactor gives the correct basal flux from PIG
                    random_seed = random_seed,
                    rf_threshold = rf_threshold,
                    pw = pw)

#
#
#make the model
#
model = Model(grid = grid,
              bed_elevation = bathy,
              params = params,
              solver_params = solver_params,
              initial_conditions = initial_conditions,
              parallel_spec = parallel_spec,
              melt_rate = idealized_anthro_melt_rate);


#
#timestepping parameters
#
niter0 = 0
dt = 0.05
#end_time = 1000.
chkpt_freq = 1000.0
pchkpt_freq = 1000.0
timestepping_params = TimesteppingParams(niter0 = niter0,
                                        dt = dt,
                                        end_time = end_time,
                                        chkpt_freq = chkpt_freq,
                                        pchkpt_freq = pchkpt_freq)

#
#output parameters
#
outputs = (h = model.fields.gh.h,
           u = model.fields.gh.u,
           v = model.fields.gh.v,
           b = model.fields.gh.b,
           s = model.fields.gh.s,
           a = model.fields.gh.accumulation,
           grfrac = model.fields.gh.grounded_fraction,
           m = model.fields.gh.basal_melt)

output_freq = 1.0
output_params = OutputParams(outputs = outputs,
                            output_freq = output_freq,
                            output_format = "mat",
                            zip_format = "nc",
                            output_start = true,
		                    output_path = member_path)

#
# assemble the simulation
#
simulation = Simulation(model = model,
                        timestepping_params = timestepping_params,
                        output_params = output_params)

#
#perform the simulation
#
run_simulation!(simulation)

return nothing

end
    
