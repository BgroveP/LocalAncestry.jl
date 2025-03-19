function predInd(prior_z,data_z,LL,breeds,nHaplo,minProb)
    ProbVec = OrderedDict(breeds .=> [Vector{Union{Missing,Float64}}(missing,nHaplo) for i in 1:length(breeds)])
    ClassVec = Vector{Union{Missing,String}}(missing,nHaplo)
#    println("size of class vec $(size(ClassVec))")
    for (r,reg) in enumerate(keys(LL))
        ###Since predicted pop data has ID at the first column I skip it, when taking z
        z     = data_z[reg]
        if in(z,keys(LL[reg]))
            indZhat,indPredClass,classProb = naive_bayes_predict(prior_z,LL[reg][z],breeds,z) #combines postprobs, and postclass
        else
            indZhat = OrderedDict(zip(breeds, zeros(length(breeds))))
            indPredClass = missing
            classProb = 0.0 #####change for refinement step may be missing
        end
        for k in breeds
            ProbVec[k][r] = indZhat[k]
        end
#        println("classProb $(classProb)")
        ######
        if classProb >= minProb
            ClassVec[r] = indPredClass
        else
            nothing
        end
        ######
    end
#    println("class vec $(ClassVec)")
    return ProbVec,ClassVec
end

function naive_bayes_predict(logpr, ll, C, z)
    #post = Dict()
        post = merge(+,logpr,ll)
        #get probs
        sumExpLogProb = sum(exp.(values(post)))
        map!(x->(exp(x)/sumExpLogProb),values(post))
        prob = maximum(values(post))
        pos = argmax(collect(values(post)))
        class = collect(keys(post))[pos]
        return post,class,prob
    end
    