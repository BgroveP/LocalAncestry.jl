
function evaluate(predla, truela)

	nhaplotypes = nrow(unique(predla[:, ["individual", "haplotype"]]))
	out = DataFrame(locus = 1:maximum(last.(truela.basepairs)), ncorrect = 0)
	reference = copy(truela)
	deleteat!(reference, .~in.(reference.individual, [predla.individual]))
	wlock = ReentrantLock()
	chunks = LocalAncestry.vecsplit(UnitRange(axes(reference, 1)), nthreads())

	#
	@threads for t in 1:nthreads()
        outt = copy(out)
        predt = copy(predla)
        reft = copy(reference)
		for i in chunks[t]
			@views tmp = predt[
				(predt.individual .== reft.individual[i]) .& (predt.haplotype .== reft.haplotype[i]) .& 
                (length.(intersect.(UnitRange.(first(reft.basepairs[i]), last(reft.basepairs[i])), predt.basepairs), ) .> 0) .& 
                (predt.ancestry .== reft.ancestry[i]),
				:,
			]
			for j in axes(tmp, 1)
				outt.ncorrect[intersect(reft.basepairs[i], tmp.basepairs[j])] .+= 1
			end
		end
        
		@lock wlock  out.ncorrect .+= outt.ncorrect
	end
	out.pcorrect = out.ncorrect/nhaplotypes
	return out
end
