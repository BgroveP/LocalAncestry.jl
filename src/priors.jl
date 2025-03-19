
function makePriors(C,id,inputPrior)
    #    println("length ID = $(length(id))")
        !isa(id,Vector) ? id = [id] : nothing 
        sameLogPrior = log.(ones(length(C)).*(1.0 / length(C)))
        logPrior = Dict()
        if isempty(inputPrior) && isempty(id)
    #        println("no prediction will be performed, only computing frequencies")
        elseif (isempty(inputPrior) && !isempty(id))
    #        println("same log prior ($(sameLogPrior)) for all breeds (and all regions!)")
            for hvec in id
               logPrior[hvec] = Dict(zip(C, sameLogPrior)) ### LOG!!!!
            end
        elseif (!isempty(inputPrior) && !isempty(id))
            if isa(inputPrior,String)
                println("Reading priors from the file (and all regions!)")
                priorRead = CSV.read(inputPrior,Tables.matrix)
            elseif isa(inputPrior,Tables.matrix)
                priorRead = deepcopy(inputPrior)
            elseif isa(inputPrior,DataFrame)
                priorRead = Matrix(inputPrior)
            else  throw(DomainError(ch, "expects A,C,G, or T"))
            end
    
            for hvec in id
                thisH = findall(priorRead[:,1] .== hvec) #returns a vector always by default
                if isempty(thisH)
                    println("no prior for $hvec")
                    logPrior[hvec] = Dict(zip(C, sameLogPrior)) ### LOG!!!!
                else
                    logPrior[hvec] = Dict(zip(C, log.(priorRead[thisH[],2:end]))) ### LOG!!!!
                end
            end
            
        else throw(DomainError(ch, "expects A,C,G, or T"))
        
        end
    #    println(logPrior)
        return logPrior
    end
    