using MEstimation
using Random
using Distributions
using Optim

struct logistic_data
    y::Vector
    x::Array{Float64}
    m::Vector
end

function logistic_nobs(data::logistic_data)
    nx = size(data.x)[1]
    ny = length(data.y)
    nm = length(data.m)
    if (nx != ny)
        error("number of rows in of x is not equal to the length of y")
    elseif (nx != nm)
        error("number of rows in of x is not equal to the length of m")
    elseif (ny != nm)
        error("length of y is not equal to the length of m")
    end
    nx
end

function logistic_loglik(theta::Vector,
    data::logistic_data,
    i::Int64)
eta = sum(data.x[i, :] .* theta)
mu = exp.(eta)./(1 .+ exp.(eta))
data.y[i] .* log.(mu) + (data.m[i] - data.y[i]) .* log.(1 .- mu)
end

# Simulate some logistic regression data
Random.seed!(123)
n = 100
m = 1
p = 10

x = Array{Float64}(undef, n, p);
x[:, 1] .= 1.0;

for j in 2:p
    x[:, j] .= rand(n);
end

true_betas = randn(p) * sqrt(p);
y = rand.(Binomial.(m, cdf.(Logistic(), x * true_betas)));
my_data = logistic_data(y, x, fill(m, n));

logistic_template = objective_function_template(logistic_nobs, logistic_loglik)
o1_ml = fit(logistic_template, my_data, true_betas)