using LennardJones

function CheckLattice(md::MD)
    nside = round( (md.n)^(1.0/3.0) );
    a = md.L/nside;
    n = md.n;
    L = md.L;
    Linv = 1.0/md.L;
    min_d2 = md.L;
    rij=zeros(3);
    for i in 1:(md.n-1)
        ri = md.r[:,i];
        for j in (i+1):md.n
            rj = md.r[:,j];
            rij = ri - rj;
            rij .-= L*round.(rij*Linv);
            d2 = sum(rij.*rij);
            min_d2 = (d2 < min_d2 ? d2 : min_d2);
        end
    end
    return sqrt(min_d2) ≈ a;
end
