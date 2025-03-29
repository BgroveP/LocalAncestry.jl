function predInd(prior_z, data_z, LL, breeds, nHaplo, minProb)
    ProbVec = OrderedDict(breeds .=> [Vector{Union{Missing,Float64}}(missing, nHaplo) for i in 1:length(breeds)])
    ClassVec = Vector{Union{Missing,String}}(missing, nHaplo)
    #    println("size of class vec $(size(ClassVec))")
    for (r, reg) in enumerate(keys(LL))
        ###Since predicted pop data has ID at the first column I skip it, when taking z
        z = data_z[reg]
        if in(z, keys(LL[reg]))
            indZhat, indPredClass, classProb = naive_bayes_predict(prior_z, LL[reg][z], breeds, z) #combines postprobs, and postclass
        else
            indZhat = OrderedDict(zip(breeds, zeros(length(breeds))))
            indPredClass = missing
            classProb = 0.0 #####change for refinement step may be missing
        end
        for k in breeds
            ProbVec[k][r] = indZhat[k]
        end
        #        println("classProb $(classProb)")
        ######
        if classProb >= minProb
            ClassVec[r] = indPredClass
        else
            nothing
        end
        ######
    end
    #    println("class vec $(ClassVec)")
    return ProbVec, ClassVec
end

function naive_bayes_predict(logpr, ll, C, z)
    #post = Dict()
    post = merge(+, logpr, ll)
    #get probs
    sumExpLogProb = sum(exp.(values(post)))
    map!(x -> (exp(x) / sumExpLogProb), values(post))
    prob = maximum(values(post))
    pos = argmax(collect(values(post)))
    class = collect(keys(post))[pos]
    return post, class, prob
end



function predictHMM(LL, populations, targetIndividuals, probStayState, targetData, haplotypeLibrary)
    blocks = string.(keys(LL))
    numberOfPopulations = length(populations)
    numberOfBlocks = length(blocks)
    numberOfIndividuals = length(targetIndividuals)
    out = OrderedDict()
    forward = zeros(Union{Missing,Float64}, numberOfPopulations, numberOfBlocks)
    backwards = zeros(Union{Missing,Float64}, numberOfPopulations, numberOfBlocks)
    initialStateProbs = repeat([1 / numberOfPopulations], numberOfPopulations)

    # Transition matrix
    transition_matrix = zeros(numberOfPopulations, numberOfPopulations)
    transition_matrix_offdiagonal = (1 - probStayState) / (numberOfPopulations - 1)
    transition_matrix = transition_matrix .+ transition_matrix_offdiagonal
    emission_matrix = zeros(numberOfPopulations, 1)

    for i in 1:numberOfPopulations
        transition_matrix[i, i] = probStayState
    end

    for i in ProgressBars.ProgressBar(1:numberOfIndividuals)
        ind = targetIndividuals[i]
        for h in 1:2
            thisRow = 2 * i + h - 2
            # forward
            for (b, r) in enumerate(blocks)
                r = ARV.rangeFromString(r)
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

            # Backwards
            for (b, r) in enumerate(reverse(blocks))
                r = ARV.rangeFromString(r)
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

            for b in 1:numberOfBlocks
                forward[:, b] = forward[:, b] .* backwards[:, b] ./ sum(forward[p, b] .* backwards[p, b] for (p, pop) in enumerate(populations))
            end

            # Return
            outname = ind * "_hap" * string(h)
            out[outname] = OrderedDict(zip(populations, [forward[j, :] for j in 1:numberOfPopulations]))

        end
    end

    # predict
    postClass = OrderedDict(zip([i * "_hap" * string(h) for h in 1:ploidity for i in targetIndividuals],
        [Vector{Union{Missing,String}}(missing, length(haplotypeLibrary)) for l in 1:length(targetIndividuals) for h in 1:ploidity]))

    for (i, id) in enumerate(targetIndividuals)
        for h in 1:ploidity
            idname = id * "_hap" * string(h)

            for (j, k) in enumerate(postClass[idname])
                maxProb = maximum(out[idname][p][j] for p in populations)
                if maxProb < minProb
                    postClass[idname][j] = missing
                else
                    postClass[idname][j] = populations[findfirst([out[idname][p][j] == maxProb for p in populations])]
                end
            end
            # Assign missing
            postClass[idname] = ARV.refineBoA(postClass[idname], out[idname])
        end
    end
    return out, postClass
end
