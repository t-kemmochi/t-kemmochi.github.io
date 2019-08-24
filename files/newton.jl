function elastic_newton(p, N, deltaT, epsilon, Nitrmax, Ntol, M, xi, w, Np0, Np1, Np2, Uold)
    U = copy(Uold)
    Nitr = 0
    E = 0.0
    dU = Inf
    while Nitr<Nitrmax && dU>Ntol
        (difU, E) = Newton(p, N, deltaT, epsilon, M, w, Np0, Np1, Np2, U, Uold)
        for i=1:2*N
            U[i] = U[i] - difU[i]
        end
        Nitr = Nitr+1
        dU = norm(difU, Inf)
    end
    return (U, E, Nitr)
end


function len_pk(p, N, k, Np1, M, U)
    len = zeros(M)
    for m=1:M
        len1 = 0.0
        len2 = 0.0
        for l=(k-p):k
            U1l = U[mod(l-1,N)+1]
            U2l = U[mod(l-1,N)+1+N]
            len1 = len1 + U1l*Np1[m,k-l+1]
            len2 = len2 + U2l*Np1[m,k-l+1]
        end
        len[m] = sqrt( len1^2 + len2^2 )
    end
    return len
end

function val10_pk(p, N, k, Np0, M, U)
    val1 = zeros(M)
    for m=1:M
        for l=(k-p):k
            val1[m] = val1[m] + U[mod(l-1,N)+1]*Np0[m,k-l+1]
        end
    end
    return val1
end

function val20_pk(p, N, k, Np0, M, U)
    val2 = zeros(M)
    for m=1:M
        for l=(k-p):k
            val2[m] = val2[m] + U[mod(l-1,N)+1+N]*Np0[m,k-l+1]
        end
    end
    return val2
end

function val11_pk(p, N, k, Np1, M, U)
    val1 = zeros(M)
    for m=1:M
        for l=(k-p):k
            val1[m] = val1[m] + U[mod(l-1,N)+1]*Np1[m,k-l+1]
        end
    end
    return val1
end

function val21_pk(p, N, k, Np1, M, U)
    val2 = zeros(M)
    for m=1:M
        for l=(k-p):k
            val2[m] = val2[m] + U[mod(l-1,N)+1+N]*Np1[m,k-l+1]
        end
    end
    return val2
end

function val12_pk(p, N, k, Np2, M, U)
    val1 = zeros(M)
    for m=1:M
        for l=(k-p):k
            val1[m] = val1[m] + U[mod(l-1,N)+1]*Np2[m,k-l+1]
        end
    end
    return val1
end

function val22_pk(p, N, k, Np2, M, U)
    val2 = zeros(M)
    for m=1:M
        for l=(k-p):k
            val2[m] = val2[m] + U[mod(l-1,N)+1+N]*Np2[m,k-l+1]
        end
    end
    return val2
end

function det_pk(p, N, M, val11U, val21U, val12U, val22U)
    val = zeros(M)
    for m=1:M
        val[m] = val11U[m]*val22U[m] - val12U[m]*val21U[m]
    end
    return val
end


function Dk(p, N, k, M, w, Np0, lenUV, val10U, val10V, val20U, val20V)
    D1 = zeros(N+p)
    D2 = zeros(N+p)
    for i=(k-p):k
        for m=1:M
            common = 0.5*lenUV[m]*Np0[m,k-i+1]*w[m]
            D1[i] = D1[i] + (val10U[m] - val10V[m])*common
            D2[i] = D2[i] + (val20U[m] - val20V[m])*common
        end
    end
    D = storeVec(p, N, k, D1, D2)
    return D
end

function storeVec(p, N, k, V1, V2)
    vec = zeros(2*N)
    if k<=N
        for i=(k-p):k
            vec[i] = V1[i]
            vec[i+N] = V2[i]
        end
    else
        for i=(k-p):N
            vec[i] = V1[i]
            vec[i+N] = V2[i]
        end
        for i=1:(k-N)
            vec[i] = V1[i+N]
            vec[i+N] = V2[i+N]
        end
    end
    return vec
end

function storeMat(p, N, k, M11, M12, M21, M22)
    mat = sparse(zeros(2*N,2*N))
    if k<=N
        for i=(k-p):k
            for j=(k-p):k
                mat[i,j] = M11[i,j]
                mat[i,j+N] = M12[i,j]
                mat[i+N,j] = M21[i,j]
                mat[i+N,j+N] = M22[i,j]
            end
        end
    else
        for i=(k-p):N
            for j=(k-p):N
                mat[i,j] = M11[i,j]
                mat[i,j+N] = M12[i,j]
                mat[i+N,j] = M21[i,j]
                mat[i+N,j+N] = M22[i,j]
            end
            for j=1:(k-N)
                mat[i,j] = M11[i,j+N]
                mat[i,j+N] = M12[i,j+N]
                mat[i+N,j] = M21[i,j+N]
                mat[i+N,j+N] = M22[i,j+N]
            end
        end
        for i=1:(k-N)
            for j=(k-p):N
                mat[i,j] = M11[i+N,j]
                mat[i,j+N] = M12[i+N,j]
                mat[i+N,j] = M21[i+N,j]
                mat[i+N,j+N] = M22[i+N,j]
            end
            for j=1:(k-N)
                mat[i,  j]   = M11[i+N,j+N]
                mat[i,  j+N] = M12[i+N,j+N]
                mat[i+N,j]   = M21[i+N,j+N]
                mat[i+N,j+N] = M22[i+N,j+N]
            end
        end
    end
    return mat
end


function JDk(p, N, k, M, w, Np0, Np1, lenUV, val10U, val10V, val20U, val20V, val11U, val11V, val21U, val21V)
    JD11 = sparse(zeros(N+p,N+p))
    JD12 = sparse(zeros(N+p,N+p))
    JD21 = sparse(zeros(N+p,N+p))
    JD22 = sparse(zeros(N+p,N+p))
    tmp = sparse(zeros(N+p,N+p))
    for i=(k-p):k
        for j=(k-p):k
            for m=1:M
                L = lenUV[m]
                V11 = val11U[m] + val11V[m]
                V21 = val21U[m] + val21V[m]
                V10 = val10U[m] - val10V[m]
                V20 = val20U[m] - val20V[m]
                common = Np1[m,k-j+1]*Np0[m,k-i+1]*w[m]/(2*L)
                JD11[i,j] = JD11[i,j] + V11*V10*common
                JD12[i,j] = JD12[i,j] + V21*V10*common
                JD21[i,j] = JD21[i,j] + V11*V20*common
                JD22[i,j] = JD22[i,j] + V21*V20*common
                tmp[i,j] = tmp[i,j] + 0.5*L*Np0[m,k-j+1]*Np0[m,k-i+1]*w[m]
            end
            JD11[i,j] = JD11[i,j] + tmp[i,j]
            JD22[i,j] = JD22[i,j] + tmp[i,j]
        end
    end
    JD = storeMat(p, N, k, JD11, JD12, JD21, JD22)
    return JD
end

function Ak(p, N, k, M, w, Np1, lenU, lenV, val11U, val11V, val21U, val21V)
    A1 = zeros(N+p)
    A2 = zeros(N+p)
    for i=(k-p):k
        for m=1:M
            common = -1*Np1[m,k-i+1]*w[m]/(lenU[m] + lenV[m])
            A1[i] = A1[i] + (val11U[m] + val11V[m])*common
            A2[i] = A2[i] + (val21U[m] + val21V[m])*common
        end
    end
    A = storeVec(p, N, k, A1, A2)
    return A
end


function JAk(p, N, k, M, w, Np1, lenU, lenV, val11U, val11V, val21U, val21V)
    JA11 = sparse(zeros(N+p,N+p))
    JA12 = sparse(zeros(N+p,N+p))
    JA21 = sparse(zeros(N+p,N+p))
    JA22 = sparse(zeros(N+p,N+p))
    tmp = sparse(zeros(N+p,N+p))
    for i=(k-p):k
        for j=(k-p):k
            for m=1:M
                LU = lenU[m]
                LV = lenV[m]
                V11 = val11U[m] + val11V[m]
                V21 = val21U[m] + val21V[m]
                common = Np1[m,k-i+1]*Np1[m,k-j+1]*w[m]/(LU*(LU + LV)^2)
                JA11[i,j] = JA11[i,j] + V11*val11U[m]*common
                JA12[i,j] = JA12[i,j] + V11*val21U[m]*common
                JA21[i,j] = JA21[i,j] + V21*val11U[m]*common
                JA22[i,j] = JA22[i,j] + V21*val21U[m]*common
                tmp[i,j] = tmp[i,j] - Np1[m,k-i+1]*Np1[m,k-j+1]*w[m]/(LU + LV)
            end
            JA11[i,j] = JA11[i,j] + tmp[i,j]
            JA22[i,j] = JA22[i,j] + tmp[i,j]
        end
    end
    JA = storeMat(p, N, k, JA11, JA12, JA21, JA22)
    return JA
end

function K1_pk(lenU, detU, detV)
    K1 = (detU + detV) ./ lenU.^5
    return K1
end

function K2_pk(lenU, lenV, detV)
    tmp1 = -1*detV.^2 .* (lenU.^8 + lenU.^6 .* lenV.^2 + lenU.^4 .* lenV.^4 + lenU.^2 .* lenV.^6 + lenV.^8)
    tmp2 = lenU.^5 .* lenV.^5 .* (lenU.^5 + lenV.^5)
    K2 = tmp1 ./ tmp2
    return K2
end



function Bk(p, N, k, M, w, Np1, lenU, lenV, val11U, val11V, val21U, val21V, val12V, val22V, detU, detV)
    B1 = zeros(N+p)
    B2 = zeros(N+p)
    K1 = K1_pk(lenU, detU, detV)
    K2 = K2_pk(lenU, lenV, detV)
    for i=(k-p):k
        for m=1:M
            common = Np1[m,k-i+1]*w[m]
            tmp1 = -K1[m]*val22V[m] - K2[m]*(val11U[m] + val11V[m])
            B1[i] = B1[i] + tmp1*common
            tmp2 = K1[m]*val12V[m] - K2[m]*(val21U[m] + val21V[m])
            B2[i] = B2[i] + tmp2*common
        end
    end
    B = storeVec(p, N, k, B1, B2)
    return B
end


function k11_pk(lenU, val11U, val21U, val22U, detU, detV)
    k111 = val22U ./ lenU.^5
    k112 = -val21U ./ lenU.^5
    k113 = -5 * (detU + detV) .* val11U ./ lenU.^7
    return (k111, k112, k113)
end

function k12_pk(lenU, val11U, val12U, val21U, detU, detV)
    k121 = -val12U ./ lenU.^5
    k122 = val11U ./ lenU.^5
    k123 = -5 * (detU + detV) .* val21U ./ lenU.^7
    return (k121, k122, k123)
end

function k2l_pk_0(lenU, lenV, detV)
    tmp1n = 2*( 4*lenU.^6 + 3*lenU.^4 .* lenV.^2 + 2*lenU.^2 .* lenV.^4 + lenV.^6 )
    tmp1d = lenU.^5 .* lenV.^5 .* (lenU.^5 + lenV.^5)
    tmp2n1 = 5*(2*lenU.^5 + lenV.^5)
    tmp2n2 = lenU.^8 + lenU.^6 .* lenV.^2 + lenU.^4 .* lenV.^4 + lenU.^2 .* lenV.^6 + lenV.^8
    tmp2n = tmp2n1 .* tmp2n2
    tmp2d = lenU.^7 .* lenV.^5 .* (lenU.^5 + lenV.^5).^2
    k2l_0 = -1*detV.^2 .* ( tmp1n ./ tmp1d - tmp2n ./ tmp2d )
    return k2l_0
end


function JBk(p, N, k, M, w, Np1, Np2, lenU, lenV, val11U, val11V, val21U, val21V, val12U, val12V, val22U, val22V, detU, detV)
    JB11 = sparse(zeros(N+p,N+p))
    JB12 = sparse(zeros(N+p,N+p))
    JB21 = sparse(zeros(N+p,N+p))
    JB22 = sparse(zeros(N+p,N+p))
    tmp = sparse(zeros(N+p,N+p))

    (k111, k112, k113) = k11_pk(lenU, val11U, val21U, val22U, detU, detV)
    (k121, k122, k123) = k12_pk(lenU, val11U, val12U, val21U, detU, detV)
    k2l_0 = k2l_pk_0(lenU, lenV, detV)
    k21 = k2l_0 .* val11U
    k22 = k2l_0 .* val21U
    K2 = K2_pk(lenU, lenV, detV)

    for i=(k-p):k
        for j=(k-p):k
            for m=1:M
                common1 = Np1[m,k-i+1]*Np1[m,k-j+1]*w[m]
                common2 = Np1[m,k-i+1]*Np2[m,k-j+1]*w[m]
                b111 = -k111[m]*val22V[m] - k113[m]*val22V[m] - k21[m]*(val11U[m]+val11V[m])
                b112 = -k112[m]*val22V[m]
                b121 = -k121[m]*val22V[m] - k123[m]*val22V[m] - k22[m]*(val11U[m]+val11V[m])
                b122 = -k122[m]*val22V[m]
                b211 = k111[m]*val12V[m] + k113[m]*val12V[m] - k21[m]*(val21U[m]+val21V[m])
                b212 = k112[m]*val12V[m]
                b221 = k121[m]*val12V[m] + k123[m]*val12V[m] - k22[m]*(val21U[m]+val21V[m])
                b222 = k122[m]*val12V[m]
                JB11[i,j] = JB11[i,j] + b111*common1 + b112*common2
                JB12[i,j] = JB12[i,j] + b121*common1 + b122*common2
                JB21[i,j] = JB21[i,j] + b211*common1 + b212*common2
                JB22[i,j] = JB22[i,j] + b221*common1 + b222*common2
                tmp[i,j] = tmp[i,j] - K2[m]*common1
            end
            JB11[i,j] = JB11[i,j] + tmp[i,j]
            JB22[i,j] = JB22[i,j] + tmp[i,j]
        end
    end
    JB = storeMat(p, N, k, JB11, JB12, JB21, JB22)
    return JB
end

function Ck(p, N, k, M, w, Np2, lenU, val11U, val21U, detU, detV)
    C1 = zeros(N+p)
    C2 = zeros(N+p)
    K1 = K1_pk(lenU, detU, detV)
    for i=(k-p):k
        for m=1:M
            common = K1[m]*Np2[m,k-i+1]*w[m]
            C1[i] = C1[i] + val21U[m]*common
            C2[i] = C2[i] - val11U[m]*common
        end
    end
    C = storeVec(p, N, k, C1, C2)
    return C
end


function JCk(p, N, k, M, w, Np1, Np2, lenU, lenV, val11U, val21U, val12U, val22U, detU, detV)
    JC11 = sparse(zeros(N+p,N+p))
    JC12 = sparse(zeros(N+p,N+p))
    JC21 = sparse(zeros(N+p,N+p))
    JC22 = sparse(zeros(N+p,N+p))
    tmp = sparse(zeros(N+p,N+p))

    (k111, k112, k113) = k11_pk(lenU, val11U, val21U, val22U, detU, detV)
    (k121, k122, k123) = k12_pk(lenU, val11U, val12U, val21U, detU, detV)
    k2l_0 = k2l_pk_0(lenU, lenV, detV)
    k21 = k2l_0 .* val11U
    k22 = k2l_0 .* val21U
    K1 = K1_pk(lenU, detU, detV)

    for i=(k-p):k
        for j=(k-p):k
            for m=1:M
                common1 = Np2[m,k-i+1]*Np1[m,k-j+1]*w[m]
                common2 = Np2[m,k-i+1]*Np2[m,k-j+1]*w[m]
                c111 = (k111[m] + k113[m])*val21U[m]
                c112 = k112[m]*val21U[m]
                c121 = (k121[m] + k123[m])*val21U[m]
                c122 = k122[m]*val21U[m]
                c211 = -(k111[m] + k113[m])*val11U[m]
                c212 = -k112[m]*val11U[m]
                c221 = -(k121[m] + k123[m])*val11U[m]
                c222 = -k122[m]*val11U[m]
                JC11[i,j] = JC11[i,j] + c111*common1 + c112*common2
                JC12[i,j] = JC12[i,j] + c121*common1 + c122*common2
                JC21[i,j] = JC21[i,j] + c211*common1 + c212*common2
                JC22[i,j] = JC22[i,j] + c221*common1 + c222*common2
                tmp[i,j] = tmp[i,j] + K1[m]*common1
            end
            JC12[i,j] = JC12[i,j] + tmp[i,j]
            JC21[i,j] = JC21[i,j] - tmp[i,j]
        end
    end
    JC = storeMat(p, N, k, JC11, JC12, JC21, JC22)
    return JC
end



function Energyk(p, N, epsilon, k, M, w, lenU, detU)
    Ek = 0.0
    for m=1:M
        cvt2 = detU[m]^2/lenU[m]^5
        Ek = Ek + (epsilon^2*cvt2 + lenU[m])*w[m]
    end
    return Ek
end


function Energy(p, N, epsilon, M, w, Np1, Np2, U)
    E = 0.0
    for k=(p+1):(p+N)
        #
        lenU = len_pk(p, N, k, Np1, M, U)
        val11U = val11_pk(p, N, k, Np1, M, U)
        val21U = val21_pk(p, N, k, Np1, M, U)
        val12U = val12_pk(p, N, k, Np2, M, U)
        val22U = val22_pk(p, N, k, Np2, M, U)
        detU = det_pk(p, N, M, val11U, val21U, val12U, val22U)
        #
        Ek = Energyk(p, N, epsilon, k, M, w, lenU, detU)
        E = E + Ek
    end
    return E
end



function Newton(p, N, deltaT, epsilon, M, w, Np0, Np1, Np2, U, V)
    F = zeros(2*N)
    JF = sparse(zeros(2*N,2*N))
    E = 0.0
    for k=(p+1):(p+N)
        #
        lenUV = len_pk(p, N, k, Np1, M, U+V)
        lenU = len_pk(p, N, k, Np1, M, U)
        lenV = len_pk(p, N, k, Np1, M, V)
        val10U = val10_pk(p, N, k, Np0, M, U)
        val10V = val10_pk(p, N, k, Np0, M, V)
        val20U = val20_pk(p, N, k, Np0, M, U)
        val20V = val20_pk(p, N, k, Np0, M, V)
        val11U = val11_pk(p, N, k, Np1, M, U)
        val11V = val11_pk(p, N, k, Np1, M, V)
        val21U = val21_pk(p, N, k, Np1, M, U)
        val21V = val21_pk(p, N, k, Np1, M, V)
        val12U = val12_pk(p, N, k, Np2, M, U)
        val12V = val12_pk(p, N, k, Np2, M, V)
        val22U = val22_pk(p, N, k, Np2, M, U)
        val22V = val22_pk(p, N, k, Np2, M, V)
        detU = det_pk(p, N, M, val11U, val21U, val12U, val22U)
        detV = det_pk(p, N, M, val11V, val21V, val12V, val22V)
        #
        D = Dk(p, N, k, M, w, Np0, lenUV, val10U, val10V, val20U, val20V)
        A = Ak(p, N, k, M, w, Np1, lenU, lenV, val11U, val11V, val21U, val21V)
        B = Bk(p, N, k, M, w, Np1, lenU, lenV, val11U, val11V, val21U, val21V, val12V, val22V, detU, detV)
        C = Ck(p, N, k, M, w, Np2, lenU, val11U, val21U, detU, detV)
        Fk = D - deltaT*A - epsilon^2*deltaT*(B + C)
        for i=1:(2*N)
            F[i] = F[i] + Fk[i]
        end
        #
        JD = JDk(p, N, k, M, w, Np0, Np1, lenUV, val10U, val10V, val20U, val20V, val11U, val11V, val21U, val21V)
        JA = JAk(p, N, k, M, w, Np1, lenU, lenV, val11U, val11V, val21U, val21V)
        JB = JBk(p, N, k, M, w, Np1, Np2, lenU, lenV, val11U, val11V, val21U, val21V, val12U, val12V, val22U, val22V, detU, detV)
        JC = JCk(p, N, k, M, w, Np1, Np2, lenU, lenV, val11U, val21U, val12U, val22U, detU, detV)
        JFk = JD - deltaT*JA - epsilon^2*deltaT*(JB + JC)
        for i=1:(2*N)
            for j=1:(2*N)
                JF[i,j] = JF[i,j] + JFk[i,j]
            end
        end
        #
        Ek = Energyk(p, N, epsilon, k, M, w, lenU, detU)
        E = E + Ek
    end
    return (JF\F, E)
end
