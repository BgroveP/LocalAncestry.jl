"""
    localancestry(
        `referencepath`::AbstractString,
        `targetpath`::AbstractString,
        `ancestrypath`::String;
        `omitpath`::String="",
        `chromosome`::Union{Int,AbstractString}="",
        `threshold`::Float64=0.66,
        `maf`::Float64=0.0001,
        `printlevel`::String="standard"
        )

    Mandatory inputs:
    ``referencepath``::AbstractString: The path to a VCF file or compressed VCF file with phased genotypes for reference individuals.
    `targetpath`::AbstractString: The path to a VCF file or compressed VCF file with phased genotypes for target individuals.
    `ancestrypath`::String: The path to a delimited file with two named columns: individual and population. 
    `omitpath`::String="": The path to a delimited file with two named columns: individual and haplotype. An empty string denotes that no haplotypes should be omitted.
    `chromosome`::Union{Int,AbstractString}="": The focal chromosome. An empty string denotes that the first encountered chromosome in the reference VCF is the focal chromosome.
    `threshold`::Float64=0.66: The lower-limit for informativeness for assignment when building haplotype blocks.
    `maf`::Float64=0.0001: The lower limit for the average minor allele frequencies across ancestral population for inclusion into the analysis.
    `printlevel`::String="standard": Options are standard and debug. Debug prints more text.

    Returns:
    `x`::DataFrame: A DataFrame object with columns *individual* with String elements, *chromosome* with String elements, *haplotype* with Int elements, *basepairs* with UnitRange elements, and *ancestry* with String elements.
"""
function localancestry(
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
    println("   Haplotype blocks: $(length(keys(library)))")

    # Estimating Local ancestries
    println("\nReading the target haplotypes")
    targetdata, targetsamples = QGIO.haplotypes(targetpath, loci = wloci)
    
    # Assign
    println("\nAssigning local ancestries")
    assignments =  LocalAncestry.assign(library, targetdata, targetsamples, popDict)
    pretty!(assignments, chromosome, wloci)
    return assignments[:,["individual", "chromosome", "haplotype", "basepairs", "ancestry"]]
end
