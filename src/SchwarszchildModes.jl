module Schwarzschild

include("SpectralSolver.jl")
include("NewtonSolver.jl")
using .NewtonSolver
using .SpectralSolver
using NLsolve
using CSV

function LoadLowSchwarzschildQnms()
    A = CSV.read(joinpath(@__DIR__, "SwData.csv"),NamedTuple)
    ωlist = Dict{NamedTuple,Complex{Float64}}()
    for i=1:length(A[:l])
        merge!(ωlist,Dict((l = A[:l][i], s = A[:s][i], n = A[:n][i]) => A[:ωreal][i]+im*A[:ωimag][i]))
    end
    ωlist
end

SchwMODESLOWWℓ = LoadLowSchwarzschildQnms();

function High_l_Guess(l,n,s)
    #Taken Directly from Eq: 17-22 of https://arxiv.org/pdf/0908.0329.pdf
    L = l+0.5
    N = n+0.5
    β = 1-s^2

    ω₋₁ = 1
    ω₀ᵢ = 1
    ω₁ = β/3 - (5/36)*N^2 - (115/432)
    ω₂ᵢ = β/9 + (235/3888)*N^2 - (1415/15552)
    ω₃ = -(1/27)*β^2 + β*(204*N^2 +211)/3888 + (854160*N^4 - 1664760*N^2 - 776939)/(40310784)
    ω₄ᵢ = (1/27)*β^2 + (1100*N^2 - 2719)/(46656) + (11273136*N^4 - 52753800*N^2 + 66480535)/(2902376448)

    ω = (ω₋₁*L + ω₁*(L^(-1)) + ω₃*(L^(-3))) - im*N*(ω₀ᵢ + ω₂ᵢ*(L^(-2)) - ω₄ᵢ*(L^(-4)))
    ω = (1/sqrt(27))*ω
end

function Low_l_Guess(l,n,s)
    SchwMODESLOWWℓ[(l=l,s=s,n=n)]
end

function High_n_Guess(l,n,s)
    print("Used higher n guess")
    # From H.-J. Blome, B. Mashhoon, Phys. Lett. A 110, 231 (1984)
    return (0.5 + l - im*(0.5 + n))/(3.0*sqrt(3))
end


function SchwarzschildGuess(l,n,s)
    if l >= 3
        return High_l_Guess(l,n,s)
    else
        return Low_l_Guess(l,n,s)
    end
end


function SchwarzschildModes(l,n,s)
    #if l < 3
    #   return SchwarzschildGuess(l,n,s),l*(l+1) - s*(s+1)
    #end
    ωguess = SchwarzschildGuess(l,n,s)
    #println(ωguess)
    Almsch = l*(l+1) - s*(s+1)
    P = ModeParameters((l=l,s=s,m=0,n=n,a=0.0,ω = ωguess,Alm = Almsch, Nmax = 10000, lmax = 10))
    #Newton!(P,0.0000001)

    function CFSch(P::ModeParameters)
        n = P.n; N = P.Nmax; ω = P.ω; Alm = P.Alm
        αₙsch(n) = (1+n)*(1+n-s-4*im*ω)
        βₙsch(n) = -1 - Alm - 2*n^2 - s + n*(-2 + im*16*ω) + im*6*ω + 24*ω^2 + ω*(im*2 + 8*ω)
        γₙsch(n) = (n-4*im*ω)*(n+s-4*im*ω)
        ContinuedFraction(αₙsch,βₙsch,γₙsch,n,N)
    end

    function f!(F,x)
        #ωᵢ = P.ω
        P.ω = x[1]+im*x[2]
        val = CF(P)
        #val = CFSch(P)
        #P.ω = ωᵢ
        F[1] = real(val)
        F[2] = imag(val)
    end

    NmaxTrials = [100,300,500,1000,5_000,10_000,50_000,100_000]
    ω₀ = P.ω; Alm₀ = P.Alm
    ωnew = P.ω + 1.0;
    iterations = 1;
    while (abs(ωnew-ω₀) > 1e-10) & (iterations < 8)
        P.Nmax = NmaxTrials[iterations]
        ω₀ = P.ω;
        sol = nlsolve(f!, [real(ωguess),imag(ωguess)],ftol = 1e-14, method = :newton)
        ωnews = sol.zero
        P.ω = ωnews[1] + im*ωnews[2]
        ωnew = P.ω;
        iterations += 1
        #println("ω diff = ",abs(ωnew-ω₀))
    end
    #println("Number of Iterations = ", iterations)
    P.ω, P.Alm
end


export SchwarzschildModes

#SchwarzschildModes(2,6,0) |> print


#0.30105345461250616 - 0.4782769832230698im
end
