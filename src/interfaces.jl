"""

"""
function get_local_ancestries(
    referencepath::AbstractString,
    targetpath::AbstractString,
    ancestrypath::String;
    omitpath::String="",
    chromosome::Union{Int,AbstractString}="",
    threshold::Float64=0.66,
    maf::Float64=0.0001,
    printlevel::String="standard"
)
    # Print input
    LocalAncestry.printinput(chromosome, referencepath, targetpath, ancestrypath, threshold, printlevel, maf)

    # Checks
    println("\nReading the locus information")
    wloci = LocalAncestry.check_loci(referencepath, targetpath)
    LocalAncestry.QGIO._print_loci(wloci)

    # Read reference ancestries
    println("\nReading the reference ancestries")
    ancestry = LocalAncestry.QGIO.read_dataframe(ancestrypath, "ancestry", ["individual", "population"], [String, String])
    LocalAncestry.QGIO._print_ancestry(ancestry)

    # Read reference omits
    if omitpath != ""
        omits = QGIO.read_dataframe(omitpath, "omission", ["individual", "haplotype"], [String, Int])
    else
        omits = DataFrame(individual=String[], haplotype=Int[])
    end
    
    # Print sample summary
    println("\nNumber of haplotypes")
    LocalAncestry.QGIO._print_samples(LocalAncestry.QGIO.samples(referencepath), ancestry, omits)

    # Subset loci (This needs a speed-up)
    println("\nSubsetting loci")
    chromosome = chromosome == "" ? wloci.chromosome[1] : LocalAncestry.QGIO.convert_chromosome(chromosome, wloci)
    if maf > NEARZERO_FLOAT
        LocalAncestry.QGIO.allelefreq!(wloci, referencepath, omits = omits, ancestries = ancestry)
        wloci[:,"maftoolow"] = wloci.allelefreq .< maf
    else
        QGIO.allelefreq!(wloci, referencepath, omits = omits)
        wloci[:,"maftoolow"] .= false
    end

    # Subset based on chromosome
    wloci[:,"wrongchromosome"] = wloci.chromosome .!= chromosome

    # Delete omitted loci
    LocalAncestry.QGIO._print_locussubset(wloci, maf)
    deleteat!(wloci, findall(wloci.maftoolow .| wloci.wrongchromosome))

    # Calculate per-locus informativeness for assignment
    QGIO.inform_for_assign!(wloci, mode = "min")
    
    println("\nReading the reference haplotypes")
    refdata, refsamples = LocalAncestry.QGIO.haplotypes(referencepath, loci = wloci, ancestries = ancestry, omits = omits)
    
    # Set popdict
    println("\nMapping reference ancestries")
    pdf = combine(groupby(refsamples, "population"), :xrow .=> [minimum, maximum] .=> ["min", "max"])
    popDict = Dict(String.(pdf.population) .=> UnitRange.(pdf.min, pdf.max))

    # Get haplotype library
    println("\nGetting haplotype library")
    library = LocalAncestry.get_haplotype_library(refdata, popDict, threshold, wloci)
    println("-number of blocks: $(length(keys(library)))")

    # Estimating Local ancestries
    println("\nReading the target haplotypes")
    targetdata, targetsamples = QGIO.haplotypes(targetpath, loci = wloci)
    
    # Assign
    println("\nAssigning local ancestries")
    assignments =  LocalAncestry.assign(library, targetdata, targetsamples, popDict)
    pretty!(assignments, chromosome, wloci)
    return assignments[:,["individual", "chromosome", "haplotype", "basepairs", "ancestry"]]
end
