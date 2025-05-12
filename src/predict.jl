
# Naive Bayes
function predictNaiveBayes!(postProb, targetData, targetIndividuals, ploidity, LL, populations, nHaplotypeBlocks, priorProb, priorLevel)
    ProbVec = OrderedDict(populations .=> [Vector{Union{Missing,Float64}}(missing, nHaplotypeBlocks) for i in 1:length(populations)])
    for (i, id) in enumerate(targetIndividuals)
        for h in 1:ploidity
            idname = id * "_hap" * string(h)
            for (r, reg) in enumerate(keys(LL))
                thisprior = returnPrior(priorProb, priorLevel, id, h, r)
                z = targetData[2*i+h-2, reg]
                if in(z, keys(LL[reg]))
                    tmpdict = naive_bayes_predict2(thisprior, LL[reg][z], populations)
                else
                    tmpdict = missingProbabilities(populations)
                end
                for p in populations
                    ProbVec[p][r] = tmpdict[p]
                end
            end
            postProb[idname] = deepcopy(ProbVec)
        end
    end
end

# Naive Bayes support function
function naive_bayes_predict2(logpr, ll, populations)
    post = merge(+, logpr, ll)
    sumExpLogProb = sum(exp.(values(post)))
    map!(x -> (exp(x) / sumExpLogProb), values(post))
    return post
end

function missingProbabilities(populations)
    return Dict{String,Union{Missing,Float64}}(zip(populations, repeat([1/length(populations)], length(populations))))
end

function returnPrior(prior, priorLevel, individual, haplotype, region)

    if priorLevel == "individual"
        x = prior[individual]
    elseif priorLevel == "haplotype"
        x = prior[individual*"_hap"*string(haplotype)]
    else
        throw(DomainError(priorLevel, "Expected 'individual', or 'haplotype'"))
    end
    return x
end

# Hidden markov model
function predictHiddenMarkov!(postProb, LL, populations, targetIndividuals, probStayState, targetData; ploidity=2)

    # Instantiate
    blocks = string.(keys(LL))
    numberOfPopulations = length(populations)
    numberOfBlocks = length(blocks)
    numberOfIndividuals = length(targetIndividuals)
    out = OrderedDict{String, OrderedDict{String, Vector{Union{Missing, Float64}}}}()
    forward = zeros(Union{Missing,Float64}, numberOfPopulations, numberOfBlocks)
    backwards = zeros(Union{Missing,Float64}, numberOfPopulations, numberOfBlocks)
    initialStateProbs = repeat([1 / numberOfPopulations], numberOfPopulations)

    # Transition matrix
    transition_matrix = getTransitionMatrix(numberOfPopulations, probStayState)
    emission_matrix = zeros(numberOfPopulations, 1)
    
    # forward
    for i in 1:numberOfIndividuals
        ind = targetIndividuals[i]
        for h in 1:ploidity
            thisRow = 2 * i + h - 2
            outname = ind * "_hap" * string(h)

            forwardsHMM(targetData, blocks, LL, thisRow, populations, emission_matrix, transition_matrix, initialStateProbs, forward)
            backwardsHMM(targetData, blocks, LL, thisRow, populations, emission_matrix, transition_matrix, initialStateProbs, backwards, numberOfBlocks)
            viterbyHMM(forward, backwards, numberOfBlocks, populations)
            
            for (p, pop) in enumerate(populations)
            postProb[outname][pop] = forward[p, :]
            end
        end
    end

    return out
end

# Support functions
function getTransitionMatrix(numberOfPopulations, probStayState)

    x = zeros(numberOfPopulations, numberOfPopulations)
    xo = (1 - probStayState) / (numberOfPopulations - 1)
    x = x .+ xo

    for i in 1:numberOfPopulations
        x[i, i] = probStayState
    end
    return x
end

function forwardsHMM(targetData, blocks, LL, thisRow, populations, emission_matrix, transition_matrix, initialStateProbs, forward)
    # forward
    for (b, r) in enumerate(blocks)
        r = rangeFromString(r)
        thisb = targetData[thisRow, r]
        # Emission matrix
        if in(thisb, keys(LL[r]))
            for (p, pop) in enumerate(populations)
                emission_matrix[p] = exp(LL[r][thisb][pop])
            end
        else
            for (p, _) in enumerate(populations)
                emission_matrix[p] = 0.5
            end
        end

        if b == 1
            forward[:, 1] = initialStateProbs .* emission_matrix
        else
            for (p1, pop1) in enumerate(populations)
                forward[p1, b] = emission_matrix[p1] * sum(forward[p2, b-1] * transition_matrix[p2, p1] for (p2, pop2) in enumerate(populations))
            end
        end
        forward[:, b] = forward[:, b] ./ sum(forward[:, b])
    end

end

function backwardsHMM(targetData, blocks, LL, thisRow, populations, emission_matrix, transition_matrix, initialStateProbs, backwards, numberOfBlocks)
            # Backwards
            for (b, r) in enumerate(reverse(blocks))
                r = rangeFromString(r)
                thisb = targetData[thisRow, r]

                # Emission matrix
                if in(thisb, keys(LL[r]))
                    for (p, pop) in enumerate(populations)
                        emission_matrix[p] = exp(LL[r][thisb][pop])
                    end
                else
                    for (p, _) in enumerate(populations)
                        emission_matrix[p] = 0.5
                    end
                end

                b = numberOfBlocks - b + 1
                if b == numberOfBlocks
                    backwards[:, b] = initialStateProbs .* emission_matrix
                else
                    for (p1, pop1) in enumerate(populations)
                        backwards[p1, b] = emission_matrix[p1] * sum(backwards[p2, b+1] * transition_matrix[p2, p1] for (p2, pop2) in enumerate(populations))
                    end
                end
                backwards[:, b] = backwards[:, b] ./ sum(backwards[:, b])
            end

end

function viterbyHMM(forward, backwards, numberOfBlocks, populations)
    for b in 1:numberOfBlocks
        forward[:, b] = forward[:, b] .* backwards[:, b] ./ sum(forward[p, b] .* backwards[p, b] for (p, pop) in enumerate(populations))
    end
end



