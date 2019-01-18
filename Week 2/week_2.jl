# Week 2 NumMeth/PDE Work
using StatsBase
using Distributions
using BenchmarkTools
# Assuming this script is run in the project folder
# Parse data into data array
file = open("data1.txt")
dat = [parse(Float64,l) for l in eachline(file)]

function estimate_gaussian_params(data::Vector)
  N = length(data)
  mean = (1/N)*sum(data)
  # Using the modified variance function
  variance = (1/(N-1))*(sum((data .- mean).^2) - (1/N)*(sum(data .- mean))^2)
  std = sqrt(variance)
  skew = (1/N)*sum(((data .- mean)./std).^3)
  kurt = ((1/N)*sum(((data .- mean)./std).^4))-3

  return Dict("Mean"=>mean,"Variance"=>variance,"Std. Deviation"=>std,"Skew"=>skew,"Kutosis"=>kurt)
end

# Find stats for dat
@benchmark stats = estimate_gaussian_params(dat)

function check_params(stats::Dict,data)
  
