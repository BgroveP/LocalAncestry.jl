
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
                    backwards, stepBack = searchBackwards(originalAssignments, prob0_ind, pos)
                else
                    stepBack = 0
                end

                if pos < length(originalAssignments)
                    forward, stepForward = searchForward(originalAssignments, prob0_ind, pos)
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

                HMMpos = (startPos+1):(endPos-1)
                if (forwardAssign * backwardsAssign > 0) && (forwardAssign != backwardsAssign)

                    # Perform HMM
                    HMMpops = [argmax(startProbs); argmax(endProbs)]
                    HMMblocks = length(HMMpos)
                    HMMforward = zeros(Float64, 2, HMMblocks)
                    HMMbackwards = zeros(Float64, 2, HMMblocks)

                    for (i, j) in enumerate(HMMpos)

                        ## Forward
                        emission = getDictionaryCrosssection(prob0_ind, j)[HMMpops]
                        if i == 1
                            HMMforward[:, i] = getDictionaryCrosssection(prob0_ind, j)[HMMpops] .* startProbs[HMMpops]
                        else
                            for m in 1:2
                                HMMforward[m, i] = emission[m] * sum(HMMforward[s, i-1] * transition_matrix[m, s] for s in 1:2)
                            end
                        end
                        HMMforward[:, i] .= HMMforward[:, i] ./ sum(HMMforward[:, i])

                        ## Backwards
                        ib = HMMblocks - i + 1
                        jb = HMMpos[ib]
                        emission = getDictionaryCrosssection(prob0_ind, jb)[HMMpops]
                        if ib == HMMblocks
                            HMMbackwards[:, ib] = getDictionaryCrosssection(prob0_ind, jb)[HMMpops] .* endProbs[HMMpops]
                        else
                            for m in 1:2
                                HMMbackwards[m, ib] = emission[m] * sum(HMMbackwards[s, ib+1] * transition_matrix[m, s] for s in 1:2)
                            end
                        end
                        HMMbackwards[:, ib] .= HMMbackwards[:, ib] ./ sum(HMMbackwards[:, ib])

                    end

                    viterbyHMM(HMMforward, HMMbackwards, HMMblocks, collect(1:2))

                    for b in 1:HMMblocks
                        postClass[k][HMMpos[b]] = populations[HMMpops[argmax(HMMforward[:, b])]]
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
