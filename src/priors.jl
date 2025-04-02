
# Archived function for constructing priors.
# It may be useful to keep some input flexibility from this function.
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

# Main function for flat priors
function priorsFlat(C, I)
    n = "individual"
    c = unique(C)
    !isa(I, Vector) ? I = [I] : nothing
    v = log.(ones(length(c)) .* (1.0 / length(c)))
    o = Dict{String, Dict{String,Float64}}()

    for i in I
                o[i] = Dict(zip(c, v))
    end
    return o, n
end

# Main function for constrained genomic regression (CGR)
function priorsCGR(referenceData, targetData, targetIndividuals, referenceOriginsVector, certainty)
    n = "haplotype"
    populations = getPopulations(referenceOriginsVector)
    p = alleleFrequencies(referenceData, referenceOriginsVector)

    # Initialize return object
    opt = optimizerCGR(populations)
    priors = repeat([1 ./ length(populations)], length(populations))

    out = performCGR(targetData, targetIndividuals, ploidity, opt, priors, certainty)

    return out, n
end

# Support function for CGR
function optimizerCGR(x::Vector{String})
    n = length(x)
    o = NLopt.Opt(:LN_COBYLA, n)
    lower_bounds!(o, zeros(Float64, n))
    NLopt.equality_constraint!(o, (x, g) -> constraintCGR(x, g), 1e-8)
    maxeval!(o, 1000)
    return o
end

function objectiveCGR(x::Vector, grad::Vector, y::Vector, X::Matrix)
    if length(grad) > 0
        grad[1] = 0
        grad[2] = 0.0
    end
    return mean((y - X * x).^2) 
end

function constraintCGR(x::Vector, grad::Vector)
    if length(grad) > 0
        grad[1] = 0.0
        grad[2] = 0.0
    end
    return sum(x) - 1
end

function certaintyScaleCGR(certainty, min_x, min_f, populations, priors)
    
    if certainty == "CGRfull"
        min_x = max.(log.(min_x), repeat([-10^10], length(populations)))
    elseif certainty == "CGR"
        min_x = max.(log.(min_f .* priors .+ (1 - min_f) .* min_x), repeat([-10^10], length(populations)))
    elseif certainty == "CGRsqrt"
        scalar = sqrt(min_f)
        min_x = max.(log.(scalar .* priors .+ (1 - scalar) .* min_x), repeat([-10^10], length(populations)))
    else
        throw(DomainError(certainty, "expects 'CGRfull' or 'CGR'"))
    end

    return min_x
end

function performCGR(targetData, targetIndividuals, ploidity, opt, priors, certainty)
    out = Dict{String,Dict{String,Float64}}()
    for (i, ind) in enumerate(targetIndividuals)
        for h in 1:ploidity
            outname = ind * "_hap" * string(h)
            inindice = 2 * i + h - 2
            y = (targetData')[:, inindice]
            NLopt.min_objective!(opt, (x, g) -> objectiveCGR(x, g, y, p))
            min_f, min_x, _ = NLopt.optimize(opt, priors)

            out[outname] = Dict(zip(populations, certaintyScaleCGR(certainty, min_x, min_f, populations, priors)))
        end
    end
    return out
end


