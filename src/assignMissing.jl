
function refineBoA(est0_ind::Vector{Union{Missing, String}},prob0_ind::OrderedDict{String, Vector{Union{Missing, Float64}}})
    newList = deepcopy(est0_ind)
    for (pos,est) in enumerate(est0_ind)
        if ismissing(est)
            current = getindex.(collect(values.(Ref(prob0_ind))),pos)
            if (pos>1) && (pos<length(est0_ind))
                backwards,stepBack = searchBackwards(est0_ind,prob0_ind,pos)
                forward,stepForward   = searchForward(est0_ind,prob0_ind,pos)
                newProb = map(*,backwards,current,forward)
            elseif pos==1
                forward,stepForward   = searchForward(est0_ind,prob0_ind,pos)
                newProb = map(*,current,forward)
            elseif pos==length(est0_ind)
                backwards,stepBack = searchBackwards(est0_ind,prob0_ind,pos)
                newProb = map(*,current,backwards)
            end
            newClass = argmax(newProb)
#            println("newClass for $pos is $newClass with $(newProb[newClass]./sum(newProb))")
            newList[pos] = collect(keys(prob0_ind))[newClass]
        end
    end
    return newList
end
