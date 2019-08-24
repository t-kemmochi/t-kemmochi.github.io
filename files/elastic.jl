include("newton.jl")

# Gauss-Legendre formula
function coeflege(n)
    if n<=1
        error("n must be > 1")
    end
    a = 1.0*zeros(1:n)
    b = 1.0*zeros(1:n)
    b[1] = 2;
    for k = 2:n
        b[k] = 1/(4-1/(k-1)^2);
    end
return (a, b)
end

function gaulege(n)
    if n<=1
        error("n must be > 1")
    end
    (a, b) = coeflege(n);
    JacM = diagm(a) + diagm(sqrt(b[2:n]),1) +  diagm(sqrt(b[2:n]),-1);
    (x, w) = eig(JacM)
    w = w[1,:]'.^2*2;
    return (x, w)
end


# B-spline basis functions
function BSbasis(knot,p,i,zeta)
    m = length(knot)
    if i<1 || i>m-p-1
        error("The index is invalid!")
    end
    Npi = zeros(zeta)
    for k=1:length(zeta)
        z = zeta[k]
        if p==0
            if z >= knot[i] && z < knot[i+1]
                Npi[k] = 1
            end
        else
            tmp1 = BSbasis(knot,p-1,i,[z])
            tmp2 = BSbasis(knot,p-1,i+1,[z])
            Npi[k] = (z-knot[i])*tmp1[1]/(knot[i+p]-knot[i]) + (knot[i+p+1]-z)*tmp2[1]/(knot[i+p+1]-knot[i+1])
        end
    end
    return Npi
end

function pBSbasis(N,p,i,zeta)
    h = 1/N
    knot = -p:(N+p)
    knot = knot*h
    Npi = BSbasis(knot,p,i,zeta)
    return Npi
end

function dBSbasis(knot,p,i,zeta)
    tmp1 = BSbasis(knot,p-1,i,zeta)
    tmp2 = BSbasis(knot,p-1,i+1,zeta)
    dNpi = p*( tmp1/(knot[i+p] - knot[i]) - tmp2/(knot[i+p+1] - knot[i+1]) )
    return dNpi
end

function ddBSbasis(knot,p,i,zeta)
    tmp1 = BSbasis(knot,p-2,i,zeta)
    tmp2 = BSbasis(knot,p-2,i+1,zeta)
    tmp3 = BSbasis(knot,p-2,i+2,zeta)
    dNpi1 = (p-1)*( tmp1/(knot[i+p] - knot[i+1]) -tmp2/(knot[i+p] - knot[i+1]) )
    dNpi2 = (p-1)*( tmp2/(knot[i+p] - knot[i+1]) -tmp3/(knot[i+p+1] - knot[i+2]) )
    ddNpi = p*( dNpi1/(knot[i+p] - knot[i]) - dNpi2/(knot[i+p+1] - knot[i+1]) )
    return ddNpi
end


function coordBScurve(p,U)
    N = Int64(length(U)/2)
    Ux = U[1:N]
    Uy = U[N+1:end]
    zeta = linspace(0,1,1000+1)
    x = zeros(zeta)
    y = zeros(zeta)
    for i=1:N
        Nip = pBSbasis(N,p,i,zeta)
        for k=1:length(zeta)
            x[k] = x[k] + Ux[i]*Nip[k]
            y[k] = y[k] + Uy[i]*Nip[k]
        end
    end
    for i=(N+1):(N+p)
        Nip = pBSbasis(N,p,i,zeta)
        for k=1:length(zeta)
            x[k] = x[k] + Ux[i-N]*Nip[k]
            y[k] = y[k] + Uy[i-N]*Nip[k]
        end
    end
    return (x, y, [Ux;Ux[1]], [Uy;Uy[1]])
end



# initialize
function initial(p,U0list,epsilon)
    N = Int64(length(U0list)/2)
    if N<2*p
        error("Take N >= 2*p")
    end
    # sample points and weights for Gauss-Legendre formula
    M = 2*p
    (xi, w) = gaulege(M)
    # reference B-spline
    refknot = 0.0:2.0:(2*p+2)
    Np0 = zeros(M,p+1)  # 0th derivative
    Np1 = zeros(M,p+1)  # 1st derivative
    Np2 = zeros(M,p+1)  # 2nd derivative
    for n=1:p+1
        Np0[:,n] =   BSbasis(refknot,p,1,xi+2*(n-1)+1)
        Np1[:,n] =  dBSbasis(refknot,p,1,xi+2*(n-1)+1)
        Np2[:,n] = ddBSbasis(refknot,p,1,xi+2*(n-1)+1)
    end
    Ux = U0list[:,1]
    Uy = U0list[:,2]
    U = [Ux; Uy]
    E = Energy(p, N, epsilon, M, w, Np1, Np2, U)
    return (N, U, E, M, xi, w, Np0, Np1, Np2)
end
