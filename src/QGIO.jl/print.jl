function _print_ancestries(x) 
    tmp = sort(combine(
        groupby(x, "population"),
        [
            "population" => (x -> sum(x .!= "unknown")) => "col1",
            "invcf" => (x -> sum(x .> 0)) => "col2",
            "omit" => (x -> sum(x)) => "col3",
            ["population", "invcf", "omit"] => ((x, y, z) -> sum((y .> 0) .& (.~z))) => "used"
        ]
    ))
  
    println("")
    println("Number of haplotypes per population")
    append!(tmp, DataFrame(population = "Any", col1 = sum(tmp.col1), col2 = sum(tmp.col2), col3 = sum(tmp.col3), used = sum(tmp.used)))
    pretty_table(tmp, header=["Reference population", "Ancestries", "Haplotypes", "Omitted", "Included"], hlines =[0,1,nrow(tmp),nrow(tmp)+1])
end


function _print_loci(loci)

    t = combine(groupby(loci, "chromosome"),
    ["position" => length => "loci",
    "format" => (x -> sum(occursin.(x,"GT"))) => "haplotypes"])
    println("   Chromosomes        $(join(t.chromosome, ", "))")
    println("   Loci               $(join(t.loci, ", "))")
    println("   Loci w. GT field   $(join(t.haplotypes, ", "))")

    return nothing
end

function _print_ancestry(ancestry)

    t = combine(groupby(ancestry, "population"), :individual => length => :N)
    sort!(t)
    println("   Populations     $(join(t.population, ", "))")
    println("   N individuals   $(join(t.N, ", "))")

    return nothing
end

function _print_samples(samples, ancestry, omits)

    t = QGIO._mergesamples(samples, ancestry, omits, removeomits = false)
    s = combine(groupby(t, :population), 
    [
        :individual => length => :N,
        [:individual, :haplotype] => ((x,y) -> length(unique(zip(x,y)))) => :unique,
        :omit => sum => :omit,
        :omit => (x -> sum(.~x)) => :use
    ])

    sappend = DataFrame(last(s))
    sappend.population .= "Total"
    sappend[1,2:end] = sum.(eachcol(s[:,2:end]))
    append!(s, sappend)

    show(s, show_row_number=false, eltypes = false, summary = false)
    println("")
end


function _print_locussubset(loci, maf)
    
    s = combine(loci, 
    [
        :chromosome => length => "n",
        :wrongchromosome => (x -> sum(.~ x)) => "chrom_ok",
        :maftoolow => (x -> sum(.~ x)) => "maf_ok",
        [:maftoolow, :wrongchromosome] => ((x,y) -> sum(.~(x .| y))) => "use",
    ])

    println("   Loci on focal chromosome:   $(s.chrom_ok[1])/$(s.n[1])")
    if maf > NEARZERO_FLOAT; println("   Loci with MAF above limit:  $(s.maf_ok[1])/$(s.n[1])");end
    println("   Loci for later analysis:    $(s.use[1])/$(s.n[1])")

    return nothing
end
