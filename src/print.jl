
function printinput(chromosome, referencepath, targetpath, ancestrypath, threshold, nbcprob, printlevel, maf)

    println("")
    println("LocalAncestry.jl v1.0.0")
    println("Started at $(Dates.format(now(), "HH:MM:SS on u d, Y")  )")
    println("")
    println("Files")
    println("   Reference        $(referencepath)")
    println("   Target           $(targetpath)")
    println("   Ancestry         $(ancestrypath)")
    println("")
    println("Parameters")
    println("   Chromosome       $(chromosome)")
    println("   MAF threshold    $(maf)")
    println("   IA threshold     $(threshold)")
    println("   NBC threshold    $(nbcprob)")
    println("   N threads        $(nthreads())")
    
    if printlevel == "debug"
        println("")
        println("Debug")
        println("   Hostname         $(gethostname())")
        println("   CPU              $(_print_cpu())")
   end

end

function _print_cpu()
    if Sys.islinux() || Sys.isapple()
        return readchomp(pipeline(`cat /proc/cpuinfo`, `grep "model name"`))
    elseif Sys.iswindows()
        return replace(split(readchomp(pipeline(`wmic cpu get name`)), '\n')[2], '\r' => "")
    else
        return "Unsupported OS for CPU name detection."
    end
end


