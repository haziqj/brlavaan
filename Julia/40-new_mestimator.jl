using DataFrames, Random, Distributions, LinearAlgebra, Cubature, Combinatorics
using MEstimation, Optim

function gen_data_bin(n = 1000; seed = nothing)
    if seed !== nothing
        Random.seed!(seed)
    end

    # Set up the loadings and covariance matrices
    Lambda = [0.80, 0.70, 0.47, 0.38, 0.34]
    neta = length(Lambda)  # q
    nitems = length(Lambda)  # p
    Psi = 1.0
    Theta = Diagonal(1 .- (Lambda .* Lambda') .* Psi)

    tau = [-1.43, -0.55, -0.13, -0.72, -1.13]

    # Generate the data
    eta = rand(MvNormal(zeros(neta), Psi), n)
    epsilon = rand(MvNormal(zeros(nitems), Theta), n)
    ystar = eta' * Lambda .+ epsilon'

    # Define y as binary data (0 and 1)
    y = DataFrame(Int.((ystar .> tau')), :auto)

    return y
end

# Function to create the summary table
function create_summary_table(data, remove_zeroes = true)
    n, p = size(data)
    summary = DataFrame(i = Int[], j = Int[], pair = String[], pattern = String[], count = Int[])
    for pair in combinations(1:p, 2)
        i, j = pair
        for a in 0:1, b in 0:1
            pattern_count = sum((data[!, i] .== a) .& (data[!, j] .== b))

            if !(remove_zeroes && pattern_count == 0)
                push!(summary, (i, j, join(string.(pair), "-"), string(a, b), pattern_count))
            end
        end
    end

    return summary
end




# Function to calculate model pairwise probabilities
function calc_model_pairwise_prob(i, j, pattern, Vy, tau)
    Vy_small = Vy[[i, j], [i, j]]
    
    function f(x)
        pdf(MvNormal(zeros(2), Vy_small), x)  # x is an array with two elements [x1, x2]
    end

    # Determine bounds for integration based on pattern
    lower_x = pattern[1] == '1' ? tau[i] : -1e2
    upper_x = pattern[1] == '1' ? 1e2 : tau[i]
    lower_y = pattern[2] == '1' ? tau[j] : -1e2
    upper_y = pattern[2] == '1' ? 1e2 : tau[j]

    # Numerically integrate the PDF over the specified region
    prob, err = hcubature(f, [lower_x, lower_y], [upper_x, upper_y])

    return prob
end

# Function to compute the pairwise log-likelihood
function pl_fn(theta, data)
    nitems = ncol(data)
    lambdas = theta[1:nitems]
    tau = theta[(nitems + 1):end]
    
    print(lambdas, "\n", tau, "\n")

    # Build the variance-covariance matrix Vy
    Lambda = reshape(lambdas, nitems, 1)
    Psi = 1.0
    Theta = zeros(nitems, nitems)
    LPLt = Lambda * Psi * Lambda'
    for i in 1:nitems
        Theta[i, i] = 1 - LPLt[i, i]
    end
    Vy = LPLt + Theta 

    summary_table = create_summary_table(data, true)  

    pl = 0.0
    for row in eachrow(summary_table)
        prob = calc_model_pairwise_prob(row.i, row.j, row.pattern, Vy, tau)
        pl += row.count * log(prob)
    end

    return pl
end

# Example usage:
dat = gen_data_bin(1000)
theta0 = [0.80, 0.70, 0.47, 0.38, 0.34, -1.43, -0.55, -0.13, -0.72, -1.13]
pl_value = pl_fn(theta0, dat)
println("Pairwise likelihood value: ", pl_value)

# Function to compute the pairwise log-likelihood for each row of the summary table
function pl_one(theta, summary_table, i)
    nitems = maximum(summary_table.j)
    lambdas = theta[1:nitems]
    tau = theta[(nitems + 1):end]

    # Build the variance-covariance matrix Vy
    Lambda = reshape(lambdas, nitems, 1)
    Psi = 1.0
    Theta = zeros(nitems, nitems)
    LPLt = Lambda * Psi * Lambda'
    for k in 1:nitems
        Theta[k, k] = 1 - LPLt[k, k]
    end
    Vy = LPLt + Theta

    pl = 0.0
    for row in eachrow(summary_table[i:i, :])
        prob = calc_model_pairwise_prob(row.i, row.j, row.pattern, Vy, tau)
        pl += row.count * log(prob)
    end

    return pl
end

function cfa_nobs(data)
    size(data)[1]
end

# Test function
theta = [0.80, 0.70, 0.47, 0.38, 0.34, -1.43, -0.55, -0.13, -0.72, -1.13]
my_data = create_summary_table(gen_data_bin(1000))
pl_one(theta, my_data, 1)

# M-estimation
cfa_template = objective_function_template(cfa_nobs, pl_one)
o1_ml = fit(cfa_template, my_data, theta)