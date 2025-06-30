
function assignCertain!(postClass, postProb, populations, LL, minProb)
    nPopulations = length(populations)
    nBlocks = length(LL)

    tmp = zeros(Union{Float64,Missing}, nPopulations, nBlocks)
    maxProb = 0.0
    maxLoc = Int(0)
    for k in keys(postClass)
        for (p, pop) in enumerate(populations)
            tmp[p, :] = postProb[k][pop]
        end

        for r in 1:nBlocks
            maxProb = maximum(tmp[:, r])
            if maxProb > minProb
                maxLoc = findfirst(maxProb .== tmp[:, r])
                postClass[k][r] = populations[maxLoc]
            end
        end
    end
end

function assignFirst!(postClass, postProb)
    for k in keys(postClass)
        prob0_ind = postProb[k]
        est0_ind = postClass[k]
        originalAssignments = deepcopy(est0_ind)

        for (pos, est) in enumerate(originalAssignments)
            if ismissing(est)
                current = getDictionaryCrosssection(prob0_ind, pos)

                if pos > 1
                    backwards, stepBack = searchBackwards(
                        originalAssignments, prob0_ind, pos
                    )
                else
                    stepBack = 0
                end

                if pos < length(originalAssignments)
                    forward, stepForward = searchForward(
                        originalAssignments, prob0_ind, pos
                    )
                else
                    stepForward = 0
                end

                # Use both
                if (stepForward > 0) && (stepBack > 0)
                    newProb = map(*, backwards, current, forward)
                elseif (stepForward > 0) && (stepBack == 0)
                    newProb = map(*, current, forward)
                elseif (stepForward == 0) && (stepBack > 0)
                    newProb = map(*, current, backwards)
                else
                    println(join(["Error for: ", k, pos], " "))
                end

                newClass = argmax(newProb)
                postClass[k][pos] = collect(keys(prob0_ind))[newClass]
            end
        end
    end
end

function assignHMM!(postClass, postProb)
    populations = string.(keys(postProb[first(keys(postProb))]))
    numberOfPopulations = length(populations)
    startProbs = zeros(Float64, numberOfPopulations)
    endProbs = zeros(Float64, numberOfPopulations)
    forwardAssign = 0
    backwardsAssign = 0
    transition_matrix = getTransitionMatrix(2, 0.999)

    for k in keys(postClass)
        prob0_ind = postProb[k]
        est0_ind = postClass[k]

        for (pos, est) in enumerate(est0_ind)
            if ismissing(est)

                # Search
                if pos > 1
                    backwards, stepBack = searchBackwards(est0_ind, prob0_ind, pos)
                else
                    stepBack = 0
                    backwards = zeros(Float64, numberOfPopulations)
                end

                if pos < length(est0_ind)
                    forward, stepForward = searchForward(est0_ind, prob0_ind, pos)
                else
                    stepForward = 0
                    forward = zeros(Float64, numberOfPopulations)
                end

                # Get probability vectors
                ## Forward
                if stepForward == 0
                    endProbs .= 1 / 2
                    forwardAssign = 0
                    endPos = length(est0_ind) + 1
                else
                    endProbs .= forward
                    endPos = pos + stepForward
                    forwardAssign = argmax(endProbs)
                end

                # Backwards
                if stepBack == 0
                    startProbs .= 1 / 2
                    backwardsAssign = 0
                    startPos = 0
                else
                    startProbs .= backwards
                    startPos = pos - stepBack
                    backwardsAssign = argmax(startProbs)
                end

                HMMpos = (startPos + 1):(endPos - 1)
                if (forwardAssign * backwardsAssign > 0) &&
                    (forwardAssign != backwardsAssign)

                    # Perform HMM
                    HMMpops = [argmax(startProbs); argmax(endProbs)]
                    HMMblocks = length(HMMpos)
                    HMMforward = zeros(Float64, 2, HMMblocks)
                    HMMbackwards = zeros(Float64, 2, HMMblocks)

                    for (i, j) in enumerate(HMMpos)

                        ## Forward
                        emission = getDictionaryCrosssection(prob0_ind, j)[HMMpops]
                        if i == 1
                            HMMforward[:, i] =
                                getDictionaryCrosssection(prob0_ind, j)[HMMpops] .*
                                startProbs[HMMpops]
                        else
                            for m in 1:2
                                HMMforward[m, i] =
                                    emission[m] * sum(
                                        HMMforward[s, i - 1] * transition_matrix[m, s] for
                                        s in 1:2
                                    )
                            end
                        end
                        HMMforward[:, i] .= HMMforward[:, i] ./ sum(HMMforward[:, i])

                        ## Backwards
                        ib = HMMblocks - i + 1
                        jb = HMMpos[ib]
                        emission = getDictionaryCrosssection(prob0_ind, jb)[HMMpops]
                        if ib == HMMblocks
                            HMMbackwards[:, ib] =
                                getDictionaryCrosssection(prob0_ind, jb)[HMMpops] .*
                                endProbs[HMMpops]
                        else
                            for m in 1:2
                                HMMbackwards[m, ib] =
                                    emission[m] * sum(
                                        HMMbackwards[s, ib + 1] * transition_matrix[m, s]
                                        for s in 1:2
                                    )
                            end
                        end
                        HMMbackwards[:, ib] .=
                            HMMbackwards[:, ib] ./ sum(HMMbackwards[:, ib])
                    end

                    viterbyHMM(HMMforward, HMMbackwards, HMMblocks, collect(1:2))

                    for b in 1:HMMblocks
                        postClass[k][HMMpos[b]] = populations[HMMpops[argmax(
                            HMMforward[:, b]
                        )]]
                    end
                elseif forwardAssign == 0
                    postClass[k][HMMpos] .= postClass[k][startPos]
                elseif backwardsAssign == 0
                    postClass[k][HMMpos] .= postClass[k][endPos]
                else
                    postClass[k][HMMpos] .= postClass[k][endPos]
                end
            end
        end
    end
end



function hmm_assign2!(postProb, postClass, h, columns, probRowdict)

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


function assign_missing2!(postProb, postClass, nPopulations, probRowdict, ploidity)
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
                        hmm_assign2!(postProb, postClass, h, a1:a2, probRowdict)
                    end
                    pos = a2
                    notdone = any(unassigned[pos:end])
                end
            end
        end
    end
end
