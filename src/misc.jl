
function vvToV(x)
    # Initialize an empty Vector{Int64} to store the result
    flattened_vector = Int64[]

    # Iterate over each vector of ranges
    for range_vector in x
        # Iterate over each range in the current vector
        for range in range_vector
            # Append each integer in the range to the flattened vector
            append!(flattened_vector, collect(range))
        end
    end

    return flattened_vector
end


function haplotypeOrigins(i, o)
    x = DataFrames.DataFrame(individual = repeat(i, inner = 2))
    return DataFrames.leftjoin(x, o, on = "individual")[!, "population"]
end


rangeChange(rObject;firstInc=0,firstDec=0,lastInc=0,lastDec=0) = (first(rObject) + firstInc - firstDec):(last(rObject) + lastInc - lastDec)

function count_tokens_T_in_class(X, t)
    count = 0
      for x in eachrow(X)
          if all(x .== t)
              count += 1
          end
      end #Loop over training set
      return count
  end
  

  
function searchForward(est0_ind,prob0_ind,pos)
    countForward = 1
    while ismissing(est0_ind[pos+countForward]) 
        countForward += 1
        if pos+countForward > length(est0_ind)
            countForward = 0
            break
        end
    end
#    println("counted $countForward steps forward for pos $pos")
    forward = getindex.(collect(values.(Ref(prob0_ind))),pos+countForward)
    return forward,countForward
end

function searchBackwards(est0_ind,prob0_ind,pos)
    countBackwards = 1
    while ismissing(est0_ind[pos-countBackwards])
        countBackwards += 1
        if pos-countBackwards < 1
        countBackwards = 0
            break
        end
    end
#    println("counted $countBackwards steps backwards for pos $pos")
    backwards = getindex.(collect(values.(Ref(prob0_ind))),pos-countBackwards)
    return backwards,countBackwards
end


"""
    rangeFromString(x)

...
# Arguments
- `x::String`: The range that should be converted from string.
...

Converts a integer range from string to actual range.
Can only convert to ranges with increments of one.

# Examples
```julia-repl
julia> BoA.rangeFromString("1:2")
1:2

julia> BoA.rangeFromString("1:10")
1:10

julia> typeof(BoA.rangeFromString("1:10"))
UnitRange{Int64}
```
"""
function rangeFromString(x)
    startAndStop = parse.(Int, split(x, ":"))
    range = startAndStop[1]:startAndStop[2]
end 

