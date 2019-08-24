function sum1(a,b,n)
    x = a
    for i=1:n
        x = x + b
    end
    println(maxabs(x))
end

function sum2(a,b,n)
    x = a
    m = length(a)
    for i=1:n
        for j=1:m
            x[j] = x[j] + b[j]
        end
    end
    println(maxabs(x))
end

function test(m,n)
    a = rand(m)
    b = rand(m)
    @time sum1(a,b,n)
    @time sum2(a,b,n)
end


# test(10000,100000)
