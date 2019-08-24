using PyPlot

# Example 1
function sine(nmax)
    close("all")
    x = linspace(0,2*pi,100)
    y = sin(x)
    line, = plot(x,y)
    for n=1:nmax
        wait(Timer(0.01))
        line[:set_ydata](sin(x-n/20))
    end
end


# Example 2
function randomwalk(n,r)
    close("all")
    x = [0.0]
    y = [0.0]
    p, = plot(x,y,marker="o")
    for i=1:n
        wait(Timer(0.01))
        theta = 2*pi*rand()
        xdif = r*cos(theta)
        ydif = r*sin(theta)
        x = x + xdif
        y = y + ydif
        p[:set_data](x,y)
    end
end


# Excercise 1
# FDM for the heat equation without external forces
function heat(u0,N,lambda,T)
    # close graphic windows
    close("all")

    # discrete parameters
    h = 1/(N+1)
    tau = lambda*h^2

    # mesh
    x = linspace(...)

    # plot initial value
    u = u0(x[2:N+1])
    sol, = plot(x,[0;u;0])

    # matrices for iteration
    A = SymTridiagonal(...)
    H = eye(N) + lambda*A

    # number of iteration
    itrmax = floor(Int64,T/tau)

    # iteration
    for i=1:itrmax
        # wait 0.01 seconds
        wait(Timer(0.01))

        # update the solution
        u = H*u

        # update the graph
        ...
    end
end


# functions for initial values
function u01(x)
    y = 10*sin(pi*x).*exp(x)
end

function u02(x)
    y = 10*(x.>0.3) .* (x.<0.9)
end

function u03(x)
    y = min(x, 1.0-x)
end






# Excercise 2
include("elastic.jl")

# functions to plot B-spline curves
function plotBScurve(p,Ulist)
    close("all")
    Ux = Ulist[:,1]
    Uy = Ulist[:,2]
    (x,y,Px,Py) = coordBScurve(p,[Ux;Uy])
    plot(x,y, c="blue")
    ax = gca()
    ax[:set_aspect](1)
end

function plotBScurveCP(p,Ulist)
    close("all")
    Ux = Ulist[:,1]
    Uy = Ulist[:,2]
    (x,y,Px,Py) = coordBScurve(p,[Ux;Uy])
    plot(Px,Py, c="skyblue", marker=".", ms=8.0)
    plot(x,y, c="blue")
    ax = gca()
    ax[:set_aspect](1)
end

# list of control points

CP0 = [cos(2*pi*(1:6)/6) sin(2*pi*(1:6)/6)]
CP1 = float([2 0; 1 1; 0 0; -1 -1; -2 0; -1 1; 0 0; 1 -1])
CP2 = [1.0 0; 1.0 1.0; 0.4 0.5; 0 0.4; -0.6 1; -1 1.5; -0.9 -0.7; -0.3 -0.5; 0 0; 0.4 -0.6; 0.9 -0.5]
CP3 = float([0 0; 1 0; 1 2; 2 2; 2 0; 5 0; 5 3; 4 3; 4 1; 3 1; 3 3; 0 3])
CP4 = float([0 0; 2 0; 2 3; -2 3; -2 -2; 4 -2; 4 3; 3 3; 3 -1; -1 -1; -1 2; 1 2; 1 1; 0 1])
CP5 = float([0 0; 2 0; 2 3; -2 3; -2 -2; 4 -2; 4 5; -3 5; -3 4; 3 4; 3 -1; -1 -1; -1 2; 1 2; 1 1; 0 1])





# elastic flow
function elastic(p,deltaT,T,U0list)
    # close graphic windows
    close("all")

    # initialize
    epsilon = 1.0
    (N, U, E, M, xi, w, Np0, Np1, Np2) = initial(p,U0list,epsilon)

    # plot initial curve
    (x,y,Px,Py) = coordBScurve(p,U)
    line, = plot(x,y)
    ax = gca()
    ax[:set_aspect](1)

    # number of iteration
    itrmax = floor(Int64,T/deltaT)

    # condition to stop iteration
    itr = 0
    gradE = 1.0
    tol = 1.0e-6
    Nitr = 0
    Nitrmax = 10
    Ntol = 1.0e-4

    # iteration
    while itr<itrmax && gradE>tol && Nitr<Nitrmax
        # wait deltaT seconds
        wait(Timer(0.1*deltaT))

        # Newton method
        Uold = copy(U)
        Eold = E
        (U, E, Nitr) = elastic_newton(p, N, deltaT, epsilon, Nitrmax, Ntol, M, xi, w, Np0, Np1, Np2, Uold)

        # update data
        (x,y,Px,Py) = coordBScurve(p,U)
        line[:set_data](x,y)

        # condition to stop
        gradE = abs(E-Eold)/deltaT
        itr = itr + 1
    end

    # end comment
    if itr == itrmax
        println("Time is up.")
    elseif gradE <= tol
        println("This is the steady state.")
    else
        println("The Newton method did not converge.")
    end
end
