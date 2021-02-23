module AttemptAtQNM

include("SpectralSolver.jl")
include("NewtonSolver.jl")
include("SchwarszchildModes.jl")
include("Interface.jl")

using .Interface

# Write your package code here.

struct Potato
    Root::Float64
end

print(GetModes(2,2,2,2))


export Potato
export GetModes

end
