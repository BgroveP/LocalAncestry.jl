# Evaluate
function evaluate(classes, map_path, true_origins, chromosome)
    individuals = unique(replace.(string.(keys(classes)), r"_hap.$" => ""))
    trueO, trueI = AlleleOrigins.readTrue(true_origins, map_path, chromosome, individuals)
    popDict = Dict{String,Int}("holstein" => 1, "jersey" => 2, "reddairy" => 3)

    individuals_with_predictions = string.(keys(classes))
    individuals_with_true = string.(trueI) .* "_hap" .* repeat(string.(1:2), outer=Int(size(trueI, 1) / 2))

    individuals_with_both = intersect(individuals_with_true, individuals_with_predictions)

    # Print if some animals not present
    if length(individuals_with_both) != length(individuals_with_predictions)
        println(string(length(individuals_with_both) - length(individuals_with_predictions)) * " individuals with predictions not in reference!")
    end

    indices_in_both = [findfirst(ind .== individuals_with_true) for ind in individuals_with_both]
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