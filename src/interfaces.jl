"""
    localancestry(referencepath, targetpath, ancestrypath;
                  omitpath="", chromosome="", threshold=0.66, maf=0.0001,
                  totcm=100.0, printlevel="standard")

Estimate local ancestries for the target individuals.

Reads phased reference and target genotypes in Variant Call Format, builds a
library of haplotype blocks from the reference individuals grouped by their
global ancestries, and infers local ancestries for the target individuals
using a hidden Markov model.

# Arguments
- `referencepath::AbstractString`: path to the reference VCF file.
- `targetpath::AbstractString`: path to the target VCF file.
- `ancestrypath::String`: path to a delimited file with the global ancestries of the reference individuals. Must contain the columns `individual` and `population`.
- `omitpath::String=""`: path to a delimited file that denotes the haplotypes to exclude from the reference panel. Must contain the columns `individual` and `haplotype`.
- `chromosome::Union{Int,AbstractString}=""`: the focal chromosome. Defaults to the first chromosome encountered in the reference VCF file.
- `threshold::Float64=0.66`: stopping criterion for building genome regions based on the informativeness for assignment statistic.
- `maf::Float64=0.0001`: lower limit for the across-population minor allele frequency for a locus to be included in the analysis.
- `totcm::Float64=100.0`: genetic length of the focal chromosome in centiMorgans, used for the state-transition probabilities of the hidden Markov model.
- `printlevel::String="standard"`: either `"standard"` or `"debug"`, where `"debug"` prints additional host information.

# Returns
- `x::DataFrame`: a DataFrame with the columns `individual` (String), `chromosome` (String), `haplotype` (Int), `basepairs` (UnitRange), and `ancestry` (String).
"""
function localancestry(
    referencepath::AbstractString,
    targetpath::AbstractString,
    ancestrypath::String;
    omitpath::String="",
    chromosome::Union{Int,AbstractString}="",
    threshold::Float64=0.66,
    maf::Float64=0.0001,
    totcm::Float64=100.0,
    printlevel::String="standard"
)
    # Print input
    printinput(chromosome, referencepath, targetpath, ancestrypath, threshold, printlevel, maf, totcm)

    # Checks
    println("\nReading the locus information")
    wloci = check_loci(referencepath, targetpath)
    QGIO._print_loci(wloci)

    # Read reference ancestries
    println("\nReading the reference ancestries")
    ancestry = QGIO.read_dataframe(ancestrypath, "ancestry", ["individual", "population"], [String, String])
    QGIO._print_ancestry(ancestry)

    # Read reference omits
    if omitpath != ""
        omits = QGIO.read_dataframe(omitpath, "omission", ["individual", "haplotype"], [String, Int])
    else
        omits = DataFrame(individual=String[], haplotype=Int[])
    end

    # Print sample summary
    println("\nNumber of haplotypes")
    QGIO._print_samples(QGIO.samples(referencepath), ancestry, omits)

    # Subset loci (This needs a speed-up)
    println("\nSubsetting loci")
    chromosome = chromosome == "" ? wloci.chromosome[1] : QGIO.convert_chromosome(chromosome, wloci)
    if maf > NEARZERO_FLOAT
        QGIO.allelefreq!(wloci, referencepath, omits = omits, ancestries = ancestry)
        wloci[:, "maftoolow"] = wloci.allelefreq .< maf
    else
        QGIO.allelefreq!(wloci, referencepath, omits = omits)
        wloci[:, "maftoolow"] .= false
    end

    # Subset based on chromosome
    wloci[:, "wrongchromosome"] = wloci.chromosome .!= chromosome

    # Delete omitted loci
    QGIO._print_locussubset(wloci, maf)
    deleteat!(wloci, findall(wloci.maftoolow .| wloci.wrongchromosome))

    # Calculate per-locus informativeness for assignment
    QGIO.inform_for_assign!(wloci, mode = "min")

    println("\nReading the reference haplotypes")
    refdata, refsamples = QGIO.haplotypes(referencepath, loci = wloci, ancestries = ancestry, omits = omits)

    # Set popdict
    println("\nMapping reference ancestries")
    pdf = combine(groupby(refsamples, "population"), :xrow .=> [minimum, maximum] .=> ["min", "max"])
    popDict = Dict(String.(pdf.population) .=> UnitRange.(pdf.min, pdf.max))

    # Get haplotype library
    println("\nGetting haplotype library")
    library = get_haplotype_library(refdata, popDict, threshold, wloci)
    println("   Genome regions: $(length(keys(library)))")

    # Estimating Local ancestries
    println("\nReading the target haplotypes")
    targetdata, targetsamples = QGIO.haplotypes(targetpath, loci = wloci)

    # Assign
    println("\nAssigning local ancestries")
    assignments = assign(library, targetdata, targetsamples, popDict; totcm = totcm)
    pretty!(assignments, chromosome, wloci)
    return assignments[:, ["individual", "chromosome", "haplotype", "basepairs", "ancestry"]]
end
