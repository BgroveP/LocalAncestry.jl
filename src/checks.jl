function assertFile(f::String, e::String)
    extensionsize = length(e)
    ## Does the file exist?
    if !isfile(f)
        throw("$f not found.")
    end

    ## Is it a .vcf file?
    if f[(end-extensionsize + 1):end] != e
        throw("File is not a $(extension) file.")
    end
end

function assertVector(x, name = "")
    if length(x) == 0
        throw("$(name) is empty.")
    end
end

function assertDataFrame(x, name = "")
        throw("Function not developed yet")
end

function assertPositiveInteger(x::Int64; errormessage = "")
    if x < 1
        throw(errormessage)
    end
end

function assertInRange(x::Int64, r::UnitRange{Int64}; errormessage = "")
    if !in(x, r)
        throw(errormessage)
    end
end
