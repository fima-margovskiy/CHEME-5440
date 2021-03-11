function calculate_transcriptional_control_array(t::Float64,x::Array{Float64,1}, problem::Dict{String,Any})::Float64

    # initialize -
    u_variable = 0.0
    
    # alias elements of the species vector -
    mRNA = x[1]
    G = x[2]
    σ70 = x[3]

    # get stuff from the problem dictionary -
    E1 = problem["E1"]
    E2 = problem["E2"]
    R = problem["ideal_gas_constant_R"]
    T_K = problem["temperature_K"]
    KD = problem["inducer_dissociation_constant"]
    n = problem["inducer_cooperativity_parameter"]

    # TODO: write u-varible function here 
    dG1 = 41.94
    dG2 = -27.67
    W1 = exp(-dG1/ (R * T_K))
    W2 = exp(-dG2/ (R * T_K))
    Fi = (σ70^n)/(((KD*1000)^n)+(σ70^n))
    u_variable = (W1 + (W2 * Fi))/(1 + W1 + (W2 * Fi))

    # return -
    return u_variable
end