# Week 2 NumMeth/PDE Work
using Distributions, Random, StatsBase

# Assuming this script is run in the project folder
# Parse data into data array
file = open("data5.txt")
dat = [parse(Float64,l) for l in eachline(file)]

"""
  estimate_gaussian_params(data)
This function estimates the gaussian parameters of a vector.
Returns a dictionary of results
"""
function estimate_gaussian_params(data::Vector)
  N = length(data)
  mean = (1/N)*sum(data)
  # Using the modified variance function
  variance = (1/(N-1))*(sum((data .- mean).^2) - (1/N)*(sum(data .- mean))^2)
  std = sqrt(variance)
  skew = (1/N)*sum(((data .- mean)./std).^3)
  kurt = ((1/N)*sum(((data .- mean)./std).^4))-3

  return Dict("Mean"=>mean,"Variance"=>variance,"Std. Deviation"=>std,"Skew"=>skew,"Kurtosis"=>kurt)
end

# Find Stats for dat
stats = estimate_gaussian_params(dat)

"""
  check_params(mean,std,length)
Given a mean and standard deviation and sample size, check our estimate_gaussian_params
routine.
Returns the gaussian params
"""
function check_params(mean::Number=0,std_dev::Number=1,length::Int64=200)
  # Generate random data
  dat = rand(Normal(mean, std_dev), length)
  stats = estimate_gaussian_params(dat)
end

# Check params using the defualt arguments
check_params()

# Check params using custom arguments
check_params(34,2,1000)

function epanechnikov(z::Number)
  if abs(t) < h
    return 1/(2*h)
  else
    return 0
  end
end

function gauss(t::Number,h::Number)
  return 1/(sqrt(pi*2)*h) * exp(-t^2/(2*h^2))
end

function rectangular(t::Number,h::Number)
  if abs(t) < h
    return 1/(2*h)
  else
    return 0
  end
end


function kernel_density_estimator(data::Vector,bw::Number,numPoints::Int64=200,kernel::Function=rectangular)
  return [(1/numPoints) * sum([kernel(x-xi,bw) for xi in data]) for x = range(minimum(data),
                                                                       stop=maximum(data),
                                                                       length=numPoints)]
end
