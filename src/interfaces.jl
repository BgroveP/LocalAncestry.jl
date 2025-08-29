"""
    get_local_ancestries(chromosome::Union{Int64,String}, 
                          referenceVCF::String, 
                          targetVCF::String, 
                          referenceAncestries::DataFrame; 
                          priorsMethod::String = "flat", 
                          minBlockSize::Int64 = 5, 
                          incrBlockSize::Int64 = 1, 
                          blockCrit::Float64 = 0.2, 
                          minNBCProb::Float64 = 0.95)
# Purpose
This function infers local ancestries. It is meant as a one-function interface to the entire inference process. 

# Arguments
- `chromosome::Union{Int64, String}`: The focal chromosome. Autosomal chromosomes can be denoted by their number as e.g.: 1, "1", or "chr1".
- `referenceVCF::String`: The relative path to .vcf file with phased genotypes of reference individuals.
- `targetVCF::String`: The relative path to .vcf file with phased genotypes of target individuals.
- `referenceAncestries::DataFrame`: Two-column (["individual", "ancestry"]) DataFrame with ancestries of reference individuals.
- `priorsMethod::String`: The method for calculating priors for the Naive Bayes Classification step (flat, CGR). We recommend flat priors for now.
- `minBlockSize::Int64`: The minimal size of haplotype blocks.
- `incrBlockSize::Int64`: The incremental size of haplotype blocks.
- `blockCrit::Float64`: The stopping criterion for building haplotype blocks. Smaller values provide larger haplotype blocks.
- `minNBCProb::Float64`: The lower threshold for posterior probabilities. Posterior probabilities above this threshold is assigned with the Naive Bayes Classification step, while those below the threshold will be assigned with the Hidden Markov step. 

# Returns
- `postProb::OrderedDict{String, Vector{OrderedDict{String, Float64}}}`: The posterior probabilities from the Naive Bayes step.
- `postClass::OrderedDict{String, Vector{String}}`: The assigned populations after the Hidden Markov model step.
- `haplotypeLibrary::OrderedDict{}`: The library of haplotype blocks.

"""
function get_local_ancestries(
    chromosome::Union{Int,AbstractString},
    referencepath::AbstractString,
    targetpath::AbstractString,
    ancestries::DataFrame;
    threshold::Float64=0.01,
    nbcprob::Float64=0.95,
)
    ## Read reference haplotype data
    println("Reading reference data")
    refdata, refloci, refind = LocalAncestry.readVCF(referencepath, chromosome)
    refancestries = LocalAncestry.haplotype_ancestries(refind, ancestries)

    # Get population information
    println("Getting population information")
    popDict = LocalAncestry.LocalAncestry.get_pop_dict(refancestries)

    # Get haplotype library
    println("Constructing haplotype blocks")
    library = LocalAncestry.get_haplotype_library(refdata, popDict, threshold)

    # Load target haplotype data
    println("Reading target data")
    targetdata, targetloci, targetind = LocalAncestry.readVCF(targetpath, chromosome)

    # Estimating Local ancestries
    println("Estimating local ancestries")
    return assign2(library, targetdata, targetind, nbcprob, popDict)
end

function assign2(library, targetdata, targetind, nbcprob, popDict)

    # Split into worker threads
    chunks = LocalAncestry.vecsplit(targetind, NCHUNKS)
    npopulations = length(keys(popDict))
    blocks = sort(UnitRange.(keys(library)))

    inddict = Dict{String,Int}(targetind .=> 1:length(targetind))
    ancestries = zeros(Int8, size(targetdata, 1), length(blocks))
    chunkends = PLOIDITY * cumsum(length.(chunks))
    chunkstarts = [1; [i + 1 for i in chunkends[1:(end-1)]]]

    # Locks
    writelock = ReentrantLock()

    # Do work
    Threads.@threads for c in 1:NCHUNKS

        # Internal initialization
        default_probabilities = repeat([1], npopulations) ./ npopulations
        individuals = chunks[c]
        probabilities = zeros(Float64, length(keys(popDict)), length(blocks))
        ancestry = zeros(Int8, length(blocks))
        internal_ancestries = zeros(Int8, PLOIDITY * length(individuals), length(blocks))
        for (i, ind) in enumerate(individuals)
            for h in 1:PLOIDITY
                ancestry .= 0

                haplotyperow = PLOIDITY * inddict[ind] - abs(h - 2)
                haplotype = @view targetdata[haplotyperow, :]

                for (ib, b) in enumerate(blocks)
                    probabilities[:, ib] = get(library[b], haplotype[b],default_probabilities)

                    if any(probabilities[:, ib] .>= nbcprob)
                        _, ancestry[ib] = findmax(probabilities[:, ib])
                    end
                end
                # Assign missing
                if any(ancestry .> 0)
                    assign_missing2!(probabilities, ancestry)
                else
                    println("Individual: ", ind, ", haplotype: ", h, " was not assigned")
                end
                ## Allocate to internal memory
                internal_ancestries[i, :] = ancestry
            end

        end
        # Allocate to common memory
        @lock writelock ancestries[chunkstarts[c]:chunkends[c], :] = internal_ancestries
    end

    return ancestries
end


function assign_missing2!(probabilities, ancestry)
    unassigned = ancestry .== 0
    assigned = .~unassigned
    pos = 1
    a1::Int = 1
    a2::Union{Int,Nothing} = 1
    notdone = any(unassigned)
    while notdone
        u = findfirst(unassigned[pos:end]) + pos - 1
        a1 = u == pos ? findfirst(assigned[pos:end]) + pos - 1 : u - 1
        if a1 > u
            ancestry[u:(a1-1)] .= ancestry[a1]
            pos = a1
            notdone = any(unassigned[pos:end])
        else
            a2 = findfirst(assigned[u:end])
            if isnothing(a2)
                ancestry[u:end] .= ancestry[a1]
                notdone = false
            else
                a2 = a2 + u - 1
                if ancestry[a1] == ancestry[a2]
                    ancestry[a1:a2] .= ancestry[a1]
                else
                    hmm_assign2!(probabilities, ancestry, a1:a2)
                end
                pos = a2
                notdone = any(unassigned[pos:end])
            end
        end
    end
end




function hmm_assign2!(probabilities, ancestry, columns)
    # Initialize
    initialAssignments = ancestry[[first(columns), last(columns)]]
    nBlocks = length(columns)
    tmpmat = zeros(Float64, 2, nBlocks)
    forward = zeros(Float64, 2, nBlocks)
    backward = zeros(Float64, 2, nBlocks)

    # Forward
    for i in 1:nBlocks
        if i == 1
            tmpmat[:, i] = [1.00, 0.0]
        elseif i == nBlocks
            tmpmat[:, i] = [0.0, 1.00]
        else
            tmpmat[:, i] = probabilities[initialAssignments, columns[i]]
        end
        tmpmat[:, i] = tmpmat[:, i] ./ sum(tmpmat[:, i])
    end

    # HMM
    for i in 1:nBlocks
        ib = nBlocks - i + 1
        if i == 1
            forward[:, i] = tmpmat[:, 1]
            backward[:, ib] = tmpmat[:, ib]
        else
            for j in 1:2
                forward[j, i] = tmpmat[j, i] * (forward[j, i-1] * (1 - HMM_STATECHANGE_PROB) + forward[(2:-1:1)[j], i-1] * HMM_STATECHANGE_PROB)
                backward[j, ib] = tmpmat[j, ib] * (backward[j, ib+1] * (1 - HMM_STATECHANGE_PROB) + backward[(2:-1:1)[j], ib+1] * HMM_STATECHANGE_PROB)
            end
            forward[:, i] = forward[:, i] ./ sum(forward[:, i])
            backward[:, ib] = backward[:, ib] ./ sum(backward[:, ib])
        end
    end

    #Viterby assignment
    ancestry[columns] = initialAssignments[[last(i) for i in findmax.(eachcol(forward .* backward))]]

    return nothing
end
