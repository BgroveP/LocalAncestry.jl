
function get_local_ancestries(
    referencepath::AbstractString,
    targetpath::AbstractString,
    ancestrypath::String;
    omitpath::String="",
    chromosome::Union{Int,AbstractString}="",
    threshold::Float64=0.13,
    nbcprob::Float64=0.99,
    maf::Float64=0.0001,
    printlevel::String="standard"
)
    # Print input
    LocalAncestry.printinput(chromosome, referencepath, targetpath, ancestrypath, threshold, nbcprob, printlevel, maf)

    # Read reference locus information
    println("\nReading the reference loci")
    refloci = LocalAncestry.QGIO.loci(referencepath)
    LocalAncestry.QGIO._print_loci(refloci)

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

    # Subset loci
    println("\nSubsetting loci")
    chromosome = chromosome == "" ? refloci.chromosome[1] : LocalAncestry.QGIO.convert_chromosome(chromosome, refloci)
    if maf > NEARZERO_FLOAT
        LocalAncestry.QGIO.allelefreq!(refloci, referencepath, omits = omits, ancestries = ancestry)
        refloci[:,"maftoolow"] = refloci.allelefreq .< maf
    else
        QGIO.allelefreq!(refloci, referencepath, omits = omits)
        refloci[:,"maftoolow"] .= false
    end

    # Subset based on chromosome
    refloci[:,"wrongchromosome"] = refloci.chromosome .!= chromosome

    # Delete omitted loci
    LocalAncestry.QGIO._print_locussubset(refloci, maf)
    deleteat!(refloci, findall(refloci.maftoolow .| refloci.wrongchromosome))

    # Calculate per-locus informativeness for assignment
    QGIO.inform_for_assign!(refloci, mode = "min")
    
    println("\nReading the reference haplotypes")
    refdata, refsamples = LocalAncestry.QGIO.haplotypes(referencepath, loci = refloci, ancestries = ancestry, omits = omits)
    
    # Set popdict
    println("\nMapping reference ancestries")
    pdf = combine(groupby(refsamples, "population"), :xrow .=> [minimum, maximum] .=> ["min", "max"])
    popDict = Dict(String.(pdf.population) .=> UnitRange.(pdf.min, pdf.max))

    # Get haplotype library
    println("\nGetting haplotype library")
    library = LocalAncestry.get_haplotype_library(refdata, popDict, threshold, refloci)

    # Estimating Local ancestries
    println("\nReading the target haplotypes")
    targetdata, targetsamples = QGIO.haplotypes(targetpath, loci = refloci)
    
    println("\nAssigning local ancestries")
    assignments =  LocalAncestry.assign(library, targetdata, targetsamples, nbcprob, popDict, printlevel)
    pretty!(assignments, chromosome, refloci)
    return assignments[:,["individual", "chromosome", "haplotype", "basepairs", "ancestry"]]
end
