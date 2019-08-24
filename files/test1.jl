function func1(x)
    y = sin(x) .* exp(-x) ./ (1+x.^2)
    return y
end

function func2(n)
    for i=1:n
        if i%3 == 0 || search(string(i), '3') > 0
            println("Aho")
        else
            println(i)
        end
    end
end

function func3(a,b)
    x1 = a .* b
    x2 = a'*b
    x3 = a*b'
    return (x1,x2,x3)
end
