function readVCF(file, chromosome)

    # Initialize chromosome strings from input
    this_chromosome = string(chromosome)
    next_chromosome = string(chromosome + 1)

    # Open the file
    reader = VCF.Reader(open(file, "r"))

    # Get individuals from VCF file
    individuals = reader.header.sampleID
    numberOfIndividuals = length(individuals)

    # Get the number of loci from the focal chromosome
    ## Initialize variables
    record = VCF.Record()
    dict_record = Dict{UInt8,UInt8}(0x30 => 0, 0x31 => 1, 0x7c => 0)
    loci = 0
    foundLargerChromosome = false

    ## Loop over rows in the vcf file
    while !(eof(reader) | foundLargerChromosome)
        record = read!(reader, record)
        recordChromosome = VCF.chrom(record)

        if recordChromosome .== this_chromosome
            loci += 1
        end
        if recordChromosome .== next_chromosome
            foundLargerChromosome = true
        end
    end
    
    # Fill matrix
    ## Initialize variables
    x = zeros(Int8, 2*numberOfIndividuals, loci)
    xRows = sort([collect(3:3:(3*numberOfIndividuals)); collect(1:3:(3*numberOfIndividuals))])
    reader = VCF.Reader(open(file, "r"))
    i = 0
    foundLargerChromosome = false

    ## Loop over rows in the VCF file
    while !(eof(reader) | foundLargerChromosome)
        record = read!(reader, record)
        recordChromosome = VCF.chrom(record)

        if recordChromosome .== this_chromosome
            i += 1
            x[:, i] = Int8.(get.([dict_record], record.data[vvToV(record.genotype)][xRows], 2))
        end
        if recordChromosome .== next_chromosome
            foundLargerChromosome = true
        end
    end

    return x, individuals
end

function readTrue(file, map, chromosome, individuals)

    genomic_map = CSV.read(map, DataFrame)
    columns_from_file = findall(genomic_map.chromosome .== chromosome) 
    out = zeros(Int32, 2*length(individuals), length(columns_from_file))
    outInd = string.(zeros(Int8, 2*length(individuals)))
    inviduals_strings = string.(individuals)
    open(file, "r") do reader
        iterator = 1
        for line in eachline(reader)
            # Convert the line to an integer
            splitted_string = split(line, "\t")
            numbers = parse.(Int, splitted_string[2:end])
            if any(splitted_string[1] .== inviduals_strings)
                out[iterator, :] = numbers[columns_from_file]
                outInd[iterator] =  splitted_string[1]
                iterator = iterator + 1
            end
        end
    end

    return out, outInd
end

