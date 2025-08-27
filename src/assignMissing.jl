


function hmm_assign!(postProb, postClass, h, columns, probRowdict)

    # Constants
    statechange = 0.00000001

    # Initialize
    initialAssignments = postClass[h, [first(columns), last(columns)]]
    nBlocks = length(columns)
    tmpmat = zeros(Float64, 2, nBlocks)
    forward = zeros(Float64, 2, nBlocks)
    backward = zeros(Float64, 2, nBlocks)

    for i in 1:nBlocks
        if i == 1
            tmpmat[:, i] = [1.00, 0.0]
        elseif i == nBlocks
            tmpmat[:, i] = [0.0, 1.00]
        else
            tmpmat[:, i] = postProb[probRowdict[h], columns[i]][[postClass[h, first(columns)], postClass[h, last(columns)]]]
        end
        tmpmat[:, i] = tmpmat[:, i] ./ sum(tmpmat[:, i])
    end

    # 
    for i in 1:nBlocks
        ib = nBlocks-i+1
        if i == 1
            forward[:, i] = tmpmat[:, 1]
            backward[:, ib] = tmpmat[:, ib]
        else
            for j in 1:2
                forward[j, i] = tmpmat[j, i] * (forward[j, i-1] * (1 - statechange) + forward[(2:-1:1)[j], i-1] * statechange)
                backward[j, ib] = tmpmat[j, ib] * (backward[j, ib+1] * (1 - statechange) + backward[(2:-1:1)[j], ib+1] * statechange)
            end
            forward[:,i] = forward[:,i] ./ sum(forward[:,i])
            backward[:,ib] = backward[:,ib] ./ sum(backward[:,ib])
        end
    end
    
    #Viterby assignment
    postClass[h, columns] = initialAssignments[[last(i)  for i in findmax.(eachcol(forward .* backward))]]

    return nothing
end


function assign_missing!(postProb, postClass, nPopulations, probRowdict, ploidity)
    nHaplotypes = size(postClass, 1)
    unassigned = postClass[1, :] .== 0
    assigned = .~unassigned
    pos = 1
    a1::Int = 1
    a2::Union{Int,Nothing} = 1

    for h in 1:nHaplotypes
        unassigned = postClass[h, :] .== 0
        assigned = .~unassigned
        pos = 1
        notdone = any(unassigned)
        while notdone
            u = findfirst(unassigned[pos:end]) + pos - 1
            a1 = u == pos ? findfirst(assigned[pos:end]) + pos - 1 : u - 1

            if a1 > u
                postClass[h, u:(a1-1)] .= postClass[h, a1]
                pos = a1
                notdone = any(unassigned[pos:end])
            else
                a2 = findfirst(assigned[u:end])
                if isnothing(a2)
                    postClass[h, u:end] .= postClass[h, a1]
                    notdone = false
                else
                    a2 = a2 + u - 1
                    if postClass[h, a1] == postClass[h, a2]
                        postClass[h, a1:a2] .= postClass[h, a1]
                    else
                        hmm_assign!(postProb, postClass, h, a1:a2, probRowdict)
                    end
                    pos = a2
                    notdone = any(unassigned[pos:end])
                end
            end
        end
    end
end
