#!/usr/bin/env julia
include("BridgeFinder_.jl")
# using BridgeFinder
if length(ARGS)<2
revcompFQ(ARGS[1],"",true)
else
revcompFQ(ARGS[1],ARGS[2],true)
end
