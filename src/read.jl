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