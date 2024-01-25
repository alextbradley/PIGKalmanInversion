# convert the model outputs and parameters files to mat format
using Pkg
Pkg.activate(".")
using TOML
using MAT
using JLD2


function main()

fpath = "./ensemble/output/";

members = 1:10
iterations = 0:19

for member in members
for iteration in iterations
	#output folder
	folder = fpath * "iteration_" * lpad(string(iteration),3,"0") * "/member_" *  lpad(string(member),3,"0");
	
	output_path = folder * "/output.jld2";

	#load it
	output_dict = load(output_path)
	model_output = output_dict["model_output"];
	
	#save in mat format
	matwrite(folder * "/output.mat", Dict(
	"model_output" => model_output))

	#repeat for parameters
	toml_path = folder * "/parameters.toml"
	param_dict = TOML.parsefile(toml_path)

	 matwrite(folder* "/parameters.mat", Dict(
        "parameters" => param_dict))

end
end

return nothing
end

main()