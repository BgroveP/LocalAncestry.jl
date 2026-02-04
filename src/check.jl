
function check_loci(r, t)

    rloci = LocalAncestry.QGIO.loci(r)
    tloci = LocalAncestry.QGIO.loci(t)

    if ~(rloci[:,1:3] == tloci[:,1:3])
        println("Warning: the reference and target don't contain the same loci. Continuing with the overlap.")
    end
    
    return innerjoin(rloci, tloci[:,1:3], on = names(rloci)[1:3])
end

