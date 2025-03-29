
function makePriors(C, id, method)
    !isa(id, Vector) ? id = [id] : nothing
    sameLogPrior = log.(ones(length(C)) .* (1.0 / length(C)))
    logPrior = Dict()
    if isempty(inputPrior) && isempty(id)
    elseif (isempty(inputPrior) && !isempty(id))
        for hvec in id
            logPrior[hvec] = Dict(zip(C, sameLogPrior)) 
        end
    elseif (!isempty(inputPrior) && !isempty(id))
        if isa(inputPrior, String)
            println("Reading priors from the file (and all regions!)")
            priorRead = CSV.read(inputPrior, Tables.matrix)
        elseif isa(inputPrior, Tables.matrix)
            priorRead = deepcopy(inputPrior)
        elseif isa(inputPrior, DataFrame)
            priorRead = Matrix(inputPrior)
        elseif isa(inputPrior, Vector)
            if isa(first(inputPrior), Dict)
                
            end
        else
            throw(DomainError(ch, "expects A,C,G, or T"))
        end

        for hvec in id
            thisH = findall(priorRead[:, 1] .== hvec) #returns a vector always by default
            if isempty(thisH)
                println("no prior for $hvec")
                logPrior[hvec] = Dict(zip(C, sameLogPrior)) ### LOG!!!!
            else
                logPrior[hvec] = Dict(zip(C, log.(priorRead[thisH[], 2:end]))) ### LOG!!!!
            end
        end

    else
        throw(DomainError(ch, "expects A,C,G, or T"))

    end
    #    println(logPrior)
    return logPrior
end

function priorsFlat(C, I, P = 2)
    !isa(I, Vector) ? I = [I] : nothing
    v = log.(ones(length(C)) .* (1.0 / length(C)))
    o = Dict{String, Dict{String, Float64}}()

    for i in I 
        for p in 1:P
        o[i * "_hap" * string(p)] = Dict(zip(C, v))
    end
end
    return o
end

function priorsCGR(referenceData, targetData, targetIndividuals, referenceOriginsVector, certainty="full")

    populations = unique(referenceOriginsVector)
    p = zeros(Float32, size(referenceData, 2), length(populations))

    for (j, pop) in enumerate(populations)
        rows = findall(pop .== referenceOriginsVector)
        p[:, j] = mean(referenceData[rows, :], dims=1)
    end

    # Initialize return object
    out = Dict{String, Dict{String, Float64}}()

    opt = NLopt.Opt(:LN_COBYLA, 3)
    lower_bounds!(opt, zeros(Float64, 3))
    NLopt.equality_constraint!(opt, (x, g) -> nloptConstraint(x, g), 1e-8)
    maxeval!(opt, 1000)
    priors = repeat([1 ./ length(populations)], length(populations))

    for (i, ind) in enumerate(targetIndividuals)
        for h in 1:2
            outname = ind * "_hap" * string(h)
            inindice = 2 * i + h - 2
            y = (targetData')[:, inindice]
            NLopt.min_objective!(opt, (x, g) -> nloptObjective(x, g, y, p))
            min_f, min_x, _ = NLopt.optimize(opt, priors)

            if certainty == "full"
                min_x = max.(log.(min_x), repeat([-10^10], length(populations)))
            elseif certainty == "auto"
                min_x = max.(log.(min_f .* priors .+ (1 - min_f) .* min_x), repeat([-10^10], length(populations)))
            elseif certainty == "autosqrt"
                scalar = sqrt(min_f)
                min_x = max.(log.(scalar .* priors .+ (1 - scalar) .* min_x), repeat([-10^10], length(populations)))
            else
                throw(DomainError(certainty, "expects 'full' or 'auto'"))
            end

            out[outname] = Dict(zip(populations, min_x))
        end
    end

    return out
end


