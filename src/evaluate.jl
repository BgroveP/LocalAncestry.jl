# Evaluate
function evaluate(classes, map_path, true_origins, library, chromosome)
    individuals = unique(replace.(string.(keys(classes)), r"_hap.$" => ""))
    trueO, trueI = readTrue(true_origins, map_path, chromosome, individuals)
    popDict = Dict{String,Int}("holstein" => 1, "jersey" => 2, "reddairy" => 3)

    individuals_with_predictions = string.(keys(classes))
    individuals_with_true =
        string.(trueI) .* "_hap" .* repeat(string.(1:2); outer=Int(size(trueI, 1) / 2))

    individuals_with_both = intersect(individuals_with_true, individuals_with_predictions)

    # Print if some animals not present
    if length(individuals_with_both) != length(individuals_with_predictions)
        println(
            string(length(individuals_with_both) - length(individuals_with_predictions)) *
            " individuals with predictions not in reference!",
        )
    end

    indices_in_both = [
        findfirst(ind .== individuals_with_true) for ind in individuals_with_both
    ]
    out = zeros(Bool, size(individuals_with_both, 1), size(trueO, 2))
    for (i, ind) in enumerate(individuals_with_both)
        for (j, key) in enumerate(keys(library))
            out[i, key] = trueO[indices_in_both[i], key] .== popDict[classes[ind][j]]
        end
    end
    overall = mean(out)
    per_locus = mean.(eachcol(out))

    return overall, per_locus
end

# Evaluate
function evaluate2(
    classes, map_path, true_origins, library, probs, chromosome, haplotype=0, minProb=0.0
)
    individuals = unique(replace.(string.(keys(classes)), r"_hap.$" => ""))
    trueO, trueI = readTrue(true_origins, map_path, chromosome, individuals)
    popDict = Dict{String,Int}("holstein" => 1, "jersey" => 2, "reddairy" => 3)

    individuals_with_predictions = string.(keys(classes))
    individuals_with_true =
        string.(trueI) .* "_hap" .* repeat(string.(1:2); outer=Int(size(trueI, 1) / 2))

    individuals_with_both = intersect(individuals_with_true, individuals_with_predictions)

    # Print if some animals not present
    if length(individuals_with_both) != length(individuals_with_predictions)
        println(
            string(length(individuals_with_both) - length(individuals_with_predictions)) *
            " individuals with predictions not in reference!",
        )
    end

    if haplotype > 0
        individuals_with_both = individuals_with_both[findall(
            occursin("_hap$haplotype", i) for i in individuals_with_both
        )]
    end
    indices_in_both = [
        findfirst(ind .== individuals_with_true) for ind in individuals_with_both
    ]

    x = DataFrame(;
        block=string.(keys(library)),
        size=length.(keys(library)),
        accuracy=0.0,
        meanMaxPostProb=0.0,
        meanRefine=0.0,
        meanComplexBlock=0.0,
    )

    for i in string.(keys(popDict)), j in string.(keys(popDict))
        k = "O" * i * "P" * j
        x[:, k] .= 0.0
    end

    # Initialize 
    firsti = string.(keys(probs))[1]
    probVec = [probs[firsti][k][1] for k in keys(probs[firsti])]
    divisor = size(individuals_with_both, 1)
    pops = string.(keys(popDict))

    # Overall accuracy
    for (j, b) in enumerate(x.block)
        r = rangeFromString(b)
        Ovec = [sum(trueO[indices_in_both, r] .== popDict[k]) for k in pops]

        for (i, ind) in enumerate(individuals_with_both)

            # Overall accuracy
            x.accuracy[j] +=
                mean(trueO[indices_in_both[i], r] .== popDict[classes[ind][j]]) / divisor

            # 
            for (m, k) in enumerate(keys(probs[ind]))
                probVec[m] = probs[ind][k][j]
            end
            x.meanMaxPostProb[j] += maximum(probVec) / divisor

            # Rate of refinement
            x.meanRefine[j] += (maximum(probVec) <= minProb) / divisor

            # Crossover mid-block
            x.meanComplexBlock[j] +=
                (size(unique(trueO[indices_in_both[i], r]), 1) > 1) / divisor

            # Confusion matrix
            predPop = classes[ind][j]
            for m in r
                truePopIndice = findfirst(trueO[indices_in_both[i], m] .== values(popDict))
                truePop = pops[truePopIndice]

                tempColumn = "O" * truePop * "P" * predPop
                x[j, tempColumn] += 1 / Ovec[truePopIndice]
            end
        end
    end

    return x
end

# Evaluate
function evaluatePerLocus(classes, map_path, true_origins, library, chromosome, haplotype=0)
    individuals = unique(replace.(string.(keys(classes)), r"_hap.$" => ""))
    trueO, trueI = readTrue(true_origins, map_path, chromosome, individuals)
    popDict = Dict{String,Int}("holstein" => 1, "jersey" => 2, "reddairy" => 3)

    individuals_with_predictions = string.(keys(classes))
    individuals_with_true =
        string.(trueI) .* "_hap" .* repeat(string.(1:2); outer=Int(size(trueI, 1) / 2))

    individuals_with_both = intersect(individuals_with_true, individuals_with_predictions)

    # Print if some animals not present
    if length(individuals_with_both) != length(individuals_with_predictions)
        println(
            string(length(individuals_with_both) - length(individuals_with_predictions)) *
            " individuals with predictions not in reference!",
        )
    end

    if haplotype > 0
        individuals_with_both = individuals_with_both[findall(
            occursin("_hap$haplotype", i) for i in individuals_with_both
        )]
    end
    indices_in_both = [
        findfirst(ind .== individuals_with_true) for ind in individuals_with_both
    ]

    x = DataFrame(; locus=1:size(trueO, 2), correct=0.0, incorrect=0.0, unassigned=0.0)

    # Overall accuracy
    divisor = size(individuals_with_both, 1)
    for (j, r) in enumerate(keys(library))
        for (i, ind) in enumerate(individuals_with_both)
            for l in r
                x.correct[l] +=
                    sum(trueO[indices_in_both[i], l] == popDict[classes[ind][j]]) / divisor
                x.unassigned[l] += ismissing(classes[ind][j]) / divisor
            end
        end
    end
    x.incorrect .= 1 .- (x.correct .+ x.unassigned)
    return x
end


# Evaluate2
function evaluatePerLocus2(classes, individuals, map_path, true_origins, library, chromosome, populations, haplotype=0)
    trueO, _ = LocalAncestry.readTrue(true_origins, map_path, chromosome, individuals)
    intpops = ["holstein", "jersey", "reddairy"]

    nLoci = size(trueO, 2)
    nIndividuals = size(trueO, 1)
    out = zeros(Int, nLoci)
    skips = haplotype > 0 ? 2 : 1
    haplotype = max(1, haplotype)

    divisor = nIndividuals/skips
    for i in 1:skips:nIndividuals
        for (k, b) in enumerate(keys(library))
            for j in b
                if intpops[trueO[i, j]] == string.(populations)[classes[i, k]]
                    out[j] += 1
                end
            end
        end
    end

    tout = out ./ divisor
    return tout
end

