module NewtonSolver

using NLsolve
include("SpectralSolver.jl")
using .SpectralSolver

mutable struct ModeParameters
    l::Int64
    m::Int64
    n::Int64
    s::Int64
    a::Float64
    ω::Complex{Float64}
    Alm::Complex{Float64}
    Nmax::Int64
    lmax::Int64
end

ModeParameters(NT::NamedTuple) = ModeParameters(NT[:l],NT[:m],NT[:n],NT[:s],NT[:a],NT[:ω],NT[:Alm],NT[:Nmax],NT[:lmax])

#Parameter Transformation
function ParameterTransformations(l,m,n,s,a,ω,Alm,Nmax,lmax)
    M = 1
    r₊ = M + sqrt(M^2 - a^2)
    r₋ = M - sqrt(M^2 - a^2)

    σ₊ = (2*ω*M*r₊ - m*a)/(r₊-r₋)
    σ₋ = (2*ω*M*r₋ - m*a)/(r₊-r₋)

    ζ₊ = im*ω
    ξ₋ = (-s-(s+2*im*σ₊))/2
    η₋ = -s+im*σ₋

    ζ = ζ₊
    ξ = ξ₋
    η = η₋

    p = (r₊ - r₋)*(ζ/2)
    α = 1 + s + ξ + η - 2*ζ + s*im*ω/ζ
    γ = 1 + s + 2*η
    δ = 1 + s + 2*ξ
    σ = Alm + (a*ω)^2 - 8*ω^2 + p*(2*α + γ - δ) + (1 + s - (γ + δ)/2 )*(s + (γ + δ)/2 )

    D₀ = δ
    D₁ = 4*p - 2*α + γ - δ -2
    D₂ = 2*α - γ + 2
    D₃ = α*(4*p - δ) - σ
    D₄ = α*(α - γ + 1)

    return ((ζ,ξ,η),(p,α,γ,δ,σ),(D₀,D₁,D₂,D₃,D₄))
end

ParameterTransformations(P::ModeParameters) = ParameterTransformations(P.l,P.m,P.n,P.s,P.a,P.ω,P.Alm,P.Nmax,P.lmax)

function ContinuedFraction(α,β,γ,n::Int,N::Int; rN = zero(α(1)))
    num = β(0)
    for i = 1:n
        num = β(i) - α(i-1)*γ(i)/num
    end
    lhs = num

    num = β(N)+α(N)*rN
    for i = (N):-1:(n+1)
       num = β(i-1) - α(i-1)*γ(i)/num
    end
    rhs = num - β(n)

    lhs + rhs
end

function rNCoeffs(D0,D1,D2,D3,D4)
    u1 = sqrt(-D0 - D1 - D2)
    u2 = (1/2)*(-4 - 2*D0 - D1)
    u3 = (8*u1 + 16*D0*u1 + 8*(D0^2)*u1 + 12*D1*u1 +
     8*D0*D1*u1 + (D1^2)*u1 + 8*D2*u1 +
     4*D0*D2*u1 - 4*D3*u1 - 4*D4*u1)/
   (8*(D0 + D1 + D2))
    u4 = (1/2)*(4 + 4*D0 + 2*D0^2 + D1 + D0*D1 - D3)
    u1,u2,u3,u4
end

function CF(P)
    _,_,(D₀,D₁,D₂,D₃,D₄) = ParameterTransformations(P)
    n = P.n; N = P.Nmax;
    αₙ(n) = n^2 + (D₀+1)*n + D₀
    βₙ(n) = -2*n^2 + (D₁+2)*n +D₃
    γₙ(n) = n^2 + (D₂-3)*n  + D₄ - D₂ + 2
    u1,u2,u3,u4 = rNCoeffs(D₀,D₁,D₂,D₃,D₄)
    rN = 1 + u1*N^(-0.5) + u2*N^(-1) + u3*N^(1.5) + u4*N^(-2)
    ContinuedFraction(αₙ,βₙ,γₙ,n,N; rN = rN)
end


function FiniteDifferenceAtParameters(P,ϵs)
    y = CF(P)
    ωᵢ = P.ω
    Almᵢ = P.Alm

    P.ω = ωᵢ + ϵs+im*0
    P.Alm,_,_ = ComputeAₗₘ(P.s,P.m,P.a*(P.ω), P.Alm,P.lmax)
    yʼr = (CF(P) - y)/ϵs

    P.ω = ωᵢ + im*ϵs
    P.Alm,_,_ = ComputeAₗₘ(P.s,P.m,P.a*(P.ω), P.Alm,P.lmax)
    yʼi = (CF(P) - y)/ϵs

    P.ω = ωᵢ
    P.Alm = Almᵢ
    return yʼr,yʼi,y
end

function Newton!(P,ϵs)
    ω₀ = P.ω
    Alm₀ = P.Alm
    function f!(F,x)
        ωᵢ = P.ω
        Almᵢ = P.Alm
        P.ω = x[1]+im*x[2]
        P.Alm,_,_ = ComputeAₗₘ(P.s,P.m,P.a*(P.ω), P.Alm,P.lmax)
        val = CF(P)
        P.ω = ωᵢ
        P.Alm = Almᵢ
        F[1] = real(val)
        F[2] = imag(val)
    end

    function j!(J,x)
        ω = x[1]+im*x[2]
        P.ω = ω
        ydr,ydi,y = FiniteDifferenceAtParameters(P,ϵs)
        J[1,1] = real(ydr)
        J[1,2] = real(ydi)
        J[2,1] = imag(ydr)
        J[2,2] = imag(ydi)
    end

    function fj!(F,J,x)
        ωᵢ = P.ω
        Almᵢ = P.Alm
        P.ω = x[1]+im*x[2]
        P.Alm,_,_ = ComputeAₗₘ(P.s,P.m,P.a*(P.ω), P.Alm,P.lmax)
        val = CF(P)
        F[1] = real(val)
        F[2] = imag(val)

        P.ω = x[1]+im*x[2] + ϵs+im*0
        P.Alm,_,_ = ComputeAₗₘ(P.s,P.m,P.a*(P.ω), P.Alm,P.lmax)
        yʼr = (CF(P) - val)/ϵs

        P.ω = x[1]+im*x[2] + im*ϵs
        P.Alm,_,_ = ComputeAₗₘ(P.s,P.m,P.a*(P.ω), P.Alm,P.lmax)
        yʼi = (CF(P) - val)/ϵs

        P.ω = ωᵢ
        P.Alm = Almᵢ

        J[1,1] = real(yʼr )
        J[1,2] = real(yʼi)
        J[2,1] = imag(yʼr )
        J[2,2] = imag(yʼi)
    end

    val = CF(P)
    df = OnceDifferentiable(f!, j!, fj!,[real(ω₀),imag(ω₀)],[real(val),imag(val)])
    sol = nlsolve(df, [real(ω₀),imag(ω₀)],ftol = 1e-14, method = :newton)

    ωnew = sol.zero
    P.ω = ωnew[1] + im*ωnew[2]
    Almnew,_,_ = ComputeAₗₘ(P.s,P.m,P.a*(P.ω), P.Alm,P.lmax)
    P.Alm = Almnew
    P
end

export ParameterTransformations, ContinuedFraction, CF, FiniteDifferenceAtParameters, Newton!
export ModeParameters


end
