module LennardJones

n = 1000;
ρ = 0.8;
T = 1.5;
L = pow(n/ρ,1.0/3.0);

r::Matrix{Float64}(undef,3,n);
v::Matrix{Float64}(undef,3,n);
f::Matrix{Float64}(undef,3,n);

function Init()
    nside = round(pow(n,1.0/3.0));
    a = L/nside;
    xhat = [1, 0, 0];
    yhat = [0, 1, 0];
    zhat = [0, 0, 1];
    i=1;
    for nx in 0:(nside-1)
        for ny in 0:(nside-1)
            for nz in 0:(nside-1)
                    r[:,i] = xhat*nx + yhat*ny + zhat*nz;
                    i+=1;
                end
            end
        end
    end

# Write your package code here.

end
