#Structure for MD Simulations
module MD_Structure

using Parameters
using Random

@with_kw struct Params
    n::Int64;
    ρ::Float64;
    T::Float64;
end

mutable struct MD
    p::Params;
    n::Int64;
    L::Float64;
    r::Matrix{Float64};
    v::Matrix{Float64};
    f::Matrix{Float64};
    kinetic_energy::Float64;
    potential_energy::Float64;
    total_energy::Float64;
    function MD()
        p_ = Params( n=1000, ρ=0.8, T=1.0 );#Default constructor: n = 1000, ρ = 0.8, T = 1.0.
        n_ = p_.n;
        ρ_ = p_.ρ;
        L_ = (n_/ρ_)^(1.0/3.0);
        r_ = Matrix{Float64}(undef,3,n_);#Column-major order
        v_ = Matrix{Float64}(undef,3,n_);
        f_ = Matrix{Float64}(undef,3,n_);
        this = new(p_,p_.n,L_,r_,v_,f_,0.0,0.0,0.0);
        Initalize!(this);
        return this;
    end
    function MD(params)
        p_ = params;#Call it with parameters defined elsewhere
        n_ = p_.n;
        ρ_ = p_.ρ;
        L_ = (n_/ρ_)^(1.0/3.0);
        r_ = Matrix{Float64}(undef,3,n_);
        v_ = Matrix{Float64}(undef,3,n_);
        f_ = Matrix{Float64}(undef,3,n_);
        this = new(p_,p_.n,L_,r_,v_,f_,0.0,0.0,0.0);
        Initalize!(this);
        return this;
    end
end

function Initalize!(md::MD)
    n = md.n;
    nside = round(n^(1.0/3.0));
    a = md.L/nside;
    T = md.p.T;
    xhat = [1, 0, 0];#Column vectors for cardinal directions
    yhat = [0, 1, 0];
    zhat = [0, 0, 1];
    i = 1;
    for nx in 0:(nside-1)
        for ny in 0:(nside-1)
            for nz in 0:(nside-1)
                md.r[:,i] =  a*( (nx+0.5)*xhat + (ny+0.5)*yhat + (nz+ 0.5)*zhat );#Build lattice
                i += 1;
            end
        end
    end
    md.v = randn(size(md.v));#Choose velocities from a normal distribution.
    for k in 1:3
        sumv = 0.0;
        for i in 1:n
            sumv += md.v[k,i];
        end
        md.v[k,:] = md.v[k,:] .- sumv/n;#Subtract off center-of-mass
    end
    Thermalize!(md);
end

function KineticEnergy(md::MD)
    return 24.0*sum(md.v.*md.v);
end

function Temperature(md::MD)
    K = KineticEnergy(md);
    T = 2.0*K/(3.0*(md.n-1));
    return T;
end

function Thermalize!(md::MD)
    Tcurrent = Temperature(md);
    σ2 = md.p.T/Tcurrent;
    md.v .*= sqrt(σ2);
end

function CalculateForces!(md::MD)
    md.f .= 0.0;
    n = md.n;
    L = md.L;
    Linv = 1.0/L;
    rij = zeros(3);
    fij = 0.0;
    vij = 0.0;
    d2 = 0.0;
    d2cut = 2.5*2.5;
    v = r -> 4*(r^(-12) - r^(-6));
    vcut = v(d2cut);
    function potential_and_force(d2)
        d2_inv = 1.0/d2;
        d6_inv = d2_inv*d2_inv*d2_inv;
        return 48.0*(d6_inv*d6_inv - 0.5*d6_inv)*d2_inv, (4.0*(d6_inv*d6_inv - d6_inv) - vcut);
    end
    md.potential_energy = 0.0;
    for i in 1:(n-1)
        ri = md.r[:,i];
        for j in (i+1):n
            rj = md.r[:,j];
            rij .= ri .- rj;
            rij .-= L*round.(rij*Linv);
            d2 = sum(rij.*rij);
            if d2 < 6.25
                fij,vij = potential_and_force(d2);
                md.f[:,i] += fij*rij;
                md.f[:,j] -= fij*rij;
                md.potential_energy += vij;
            end
        end
    end
end



end
