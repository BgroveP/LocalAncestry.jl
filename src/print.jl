chromosome = "chr22"
targetpath = "test/data/AdmixPop_human1_chr22_admixed_all.vcf.gz"
referencepath = "test/data/AdmixPop_human1_chr22_ancestral_all.vcf.gz"
ancestries = "test/data/AdmixPop_human1_ancestral.csv"
threshold = 0.01
nbcprob = 0.95
printlevel = "standard"
i = [referencepath, targetpath, ancestries, chromosome, threshold, nbcprob, printlevel]
function printlog(message, level, type)

end


function printinput(i)

    o = [
        "reference haplotypes",
        "target haplotypes",
        "reference ancestries",
        "focal chromosome",
        "informativeness for assignment threshold",
        "NBC probability",
        "printlevel",
        "number of threads"
    ]
    #
    ancestryi = findfirst(o .== "reference ancestry")
    if isa(i[ancestryi], DataFrame)
        i[ancestryi] = "Julia DataFrame"
    end

    indict = Dict{String,String}(o .=> [string.(i); string(Threads.nthreads())])
    indent = maximum(length.(keys(indict)))

    println("")
    println("LocalAncestry.jl (v1.0.0)")
    for oi in o
        printneat(oi, indict[oi], 1, indent + 3, ":")
    end

end



using ProgressMeter

function do_work(i)
    sleep(rand() * 0.1)  # Simulate work
    return i^2
end

function main()
    n = 10000
    results = Vector{Int}(undef, n)
    progress = Progress(n, dt=0.1, desc="Processing: ")
    lock = ReentrantLock()

    Threads.@threads for i in 1:n
        results[i] = do_work(i)
        @lock lock update!(progress, 1)  # Thread-safe update
    end
end

main()

function get_popdict(individuals, ancestries)

    # Checks
    ## All vcf individuals are unique
    if length(individuals) != length(unique(individuals))
        error("Some individuals are repeated in the vcf file")
    end

    ## All ancestry individuals are unique
    if length(ancestries.individual) != length(unique(ancestries.individual))
        error("Some individuals are repeated in the reference ancestries")
    end

    
    # Standardize missing values
    ## Join the two information sources
    o = outerjoin(ancestries, DataFrame(individual = individuals, invcf = 1:length(individuals)), on = "individual")
    
    ## Translate missing values
    o.population[ismissing.(o.population)] .= "unknown"
    o.invcf[ismissing.(o.invcf)] .= 0
    o.omit[ismissing.(o.omit)] .= 2

    # Combine to obtain summaries
    o2 = sort(combine(
        groupby(o, "population"),
        [
            "population" => (x -> sum(x .!= "unknown")) => "col1",
            "invcf" => (x -> sum(x .> 0)) => "col2",
            ["population", "invcf"] => ((x, y) -> sum((x .!= "unknown") .& (y .> 0))) => "overlap"
        ]
    ))
  
    # Print
    print_popdict(o2)

    # Get: VCF column => x row, VCF columns to use
    o3 = sort(repeat(o[(o.invcf .> 0), :], inner = PLOIDITY), "invcf")
    o3.vcfcol = 1:nrow(o3)
    o3.haplotype = repeat(1:PLOIDITY, outer = Int(nrow(o3)/PLOIDITY))

    out1 = (o3.haplotype .!= o3.omit) .& (o3.population .!= "unknown")
    o4 = sort(o3[out1,:], "population")
    o4.xrow = 1:nrow(o4)
    sort!(o4, "invcf")
    out2 = o4.xrow

    o2.start = (1 .+ cumsum([0;o2.overlap]))[1:nrow(o2)]
    o2.end = o2.start + o2.overlap .- 1
    o2 = o2[(o2.overlap .> 0) .& (1:nrow(o2) .< nrow(o2)),:]
    out3 = Dict(o2.population .=> UnitRange.(o2.start, o2.end))
    return out1, out2, out3
end

function print_popdict(o2)

     # Print
    println("")
    println("Number of individuals per population")
    append!(o2, DataFrame(population = "Any", col1 = sum(o2.col1), col2 = sum(o2.col2), overlap = sum(o2.overlap)))
    pretty_table(o2, header=["Reference population", "Has known ancestry", "In .vcf file", "Overlap"], hlines =[0,1,nrow(o2),nrow(o2)+1])
end
