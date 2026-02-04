
function read_dataframe(p, ft, m, ct)
    
    d = CSV.read(p, DataFrame,stringtype=String)

    # Check for mandatory headers
    for tm in m
        if ~any(tm .== names(d))
            error("The $(ft) file doesn't contain the mandatory header: " * tm)
        end
    end

    # Convert headers
    for  (im,tm) in enumerate(m)
        if typeof(d[:, tm]) != Vector{ct[im]}
            d[!, tm] = string.(d[:, tm])
        end
    end
    return d
end
