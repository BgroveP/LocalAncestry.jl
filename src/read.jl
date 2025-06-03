function getVCFrows(file, c)
    # Initialize
    record = VCF.Record()
    loci::Int64 = 0
    foundFocalChromosome = false
    this_chromosome = [string(c), "chr" * string(c)]
    reader = VCF.Reader(open(file, "r"))

    while !eof(reader)
        try
            read!(reader, record)
        catch e
            break
        end
        recordChromosome = VCF.chrom(record)

        if any(recordChromosome .== this_chromosome)
            loci += 1
            if !foundFocalChromosome
                foundFocalChromosome = true
            end
        else
            if foundFocalChromosome
                break
            end
        end
    end

    close(reader)

    return loci
end

function getVCFdata(f, c, l, n)

    # Initialize
    x = zeros(Int8, 2 * n, l)
    this_chromosome = [string(c), "chr" * string(c)]
    foundFocalChromosome = false
    record = VCF.Record()
    xRows = sort([collect(3:3:(3 * n)); collect(1:3:(3 * n))])
    i = 0
    dict_record = Dict{UInt8,UInt8}(0x30 => 0, 0x31 => 1, 0x7c => 0)
    reader = VCF.Reader(open(f, "r"))

    while !eof(reader)
        try
            read!(reader, record)
        catch e
            break
        end
        recordChromosome = VCF.chrom(record)

        if any(recordChromosome .== this_chromosome)
            i += 1
            x[:, i] = Int8.(
                get.([dict_record], record.data[vvToV(record.genotype)][xRows], 2)
            )
            if !foundFocalChromosome
                foundFocalChromosome = true
            end
        else
            if foundFocalChromosome
                break
            end
        end
    end
    close(reader)
    return x
end

function getVCFindividuals(f)::Vector{String}
    r = VCF.Reader(open(f, "r"))
    x = r.header.sampleID
    close(r)
    return x
end

function readVCF(file::String, chromosome::Int64)

    # Checks
    assertFile(file, ".vcf")
    assertPositiveInteger(
        chromosome;
        errormessage="chromosome was $(chromosome), but should be an integer larger than zero",
    )

    # Get the number of loci from the focal chromosome
    individuals = getVCFindividuals(file)
    loci = getVCFrows(file, chromosome)
    x = getVCFdata(file, chromosome, loci, length(individuals))

    return x, individuals
end

function readTrue(file::String, map::String, chromosome::Int64, individuals::Vector{String})
    genomic_map = CSV.read(map, DataFrame)
    columns_from_file = findall(genomic_map.chromosome .== chromosome)
    out = zeros(Int32, 2 * length(individuals), length(columns_from_file))
    outInd = string.(zeros(Int8, 2 * length(individuals)))
    inviduals_strings = string.(individuals)
    open(file, "r") do reader
        iterator = 1
        for line in eachline(reader)
            # Convert the line to an integer
            splitted_string = split(line, "\t")
            numbers = parse.(Int, splitted_string[2:end])
            if any(splitted_string[1] .== inviduals_strings)
                out[iterator, :] = numbers[columns_from_file]
                outInd[iterator] = splitted_string[1]
                iterator = iterator + 1
            end
        end
    end

    return out, outInd
end

function readRfmix(path::String, mappath::String)
    if !isfile(path)
        error("Assignment file not found!")
    end

    if !isfile(mappath)
        error("Map file not found!")
    end

    # Read the data
    data = CSV.read(path, DataFrame; header=2)
    genomic_map = CSV.read(mappath, DataFrame)

    individuals = unique(replace.(names(data)[5:end], r":::.*" => ""))
    populations = unique(replace.(names(data)[5:end], r".*:::hap\d+:::" => ""))
    haplotypes = unique(
        replace.(replace.(names(data)[5:end], r".*:::hap" => ""), r":::.*" => "")
    )
    chromosome = data.chromosome[1]
    nMarkers = sum(genomic_map.chromosome .== chromosome)
    nBlocks = size(data, 1)
    blocks = Vector{UnitRange{Int64}}(undef, size(data, 1))
    for (i, j) in enumerate(data.genetic_marker_index)
        blockstart = j + 1
        blockend = i == nBlocks ? nMarkers : data.genetic_marker_index[i + 1]
        blocks[i] = blockstart:blockend
    end

    # Fill library
    library = OrderedDict{UnitRange{Int64},Vector{Matrix{Int8}}}()
    for i in 1:nBlocks
        library[blocks[i]] = [zeros(Int8, 1, 1)]
    end

    # Fill class
    class = OrderedDict{String,Vector{Union{Missing,String}}}()
    levels1 = [i * "_hap" * h for i in individuals for h in haplotypes]
    levels2 = [i * ":::hap" * h * ":::" for i in individuals for h in haplotypes]

    for (i, ind) in enumerate(levels1)
        probcolumns = levels2[i] .* populations
        class[ind] = [
            populations[argmax([data[j, c] for c in probcolumns])] for j in 1:nBlocks
        ]
    end

    # Fill probabilities
    probs = OrderedDict{String,OrderedDict{String,Vector{Union{Missing,Float64}}}}()
    for (i, ind) in enumerate(levels1)
        probs[ind] = OrderedDict{String,Vector{Union{Missing,Float64}}}()
        for p in populations
            probs[ind][p] = data[:, levels2[i] .* p]
        end
    end

    # Return
    return probs, class, library
end
