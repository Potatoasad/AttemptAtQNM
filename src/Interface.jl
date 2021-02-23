module Interface

include("SpectralSolver.jl")
include("NewtonSolver.jl")
include("SchwarszchildModes.jl")

using .NewtonSolver
using .Schwarzschild
using .SpectralSolver

#P = ModeParameters((l=2,s=2,m=0,n=0,a=0.01,ω = 0.373672 - 0.0889623*im,Alm = 0, Nmax = 300, lmax = 10))

function GetModes(l,m,n,s; amax = 0.99, ϵ = 0.01, Nmax = 300, lmax = 10)
    ##Initialize
    as = [a for a in 0.0:ϵ:amax]
    ωs = zero(im*(as))
    Almss = zero(im*(as))

    #Get Schwarzschild Modes
    ωsch, Almsch = SchwarzschildModes(l,n,s)

    ϵs = 0.0000001 #Finite difference step size

    P = ModeParameters((l=l,s=s,m=m,n=n,a=0.0,ω = ωsch,Alm = Almsch, Nmax = Nmax, lmax = lmax))

    NmaxTrials = [100,300,500,1000,5_000,10_000,50_000]
    lmaxTrials = [10,12,15,15,20,25]

    for i =1:length(as)
        #print(i," ")
        P.a = as[i]
        ω₀ = P.ω; Alm₀ = P.Alm
        ωnew = P.ω + 1.0; Almnew = Alm₀ + 1.0;
        iterations = 1;
        while ((abs(ωnew-ω₀) > 1e-10) | (abs(Almnew-Alm₀) > 1e-10)) & (iterations < 6)
            P.Nmax = NmaxTrials[iterations]
            P.lmax = lmaxTrials[iterations]
            ω₀ = P.ω; Alm₀ = P.Alm;
            Newton!(P,ϵs)
            ωnew = P.ω; Almnew = P.Alm
            iterations += 1
            #println("ω diff = ",abs(ωnew-ω₀)," Alm diff =  ",abs(Almnew-Alm₀))
        end
        #println("Number of Iterations = ", iterations)
        ωs[i] = P.ω
        Almss[i] = P.Alm
    end

    final = [i for i in zip(as,ωs,Almss)]
end

Ms = GetModes(2,2,2,2)

#for i in Ms
#   println(i)
#end

export GetModes
end
