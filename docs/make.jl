using Documenter
using LocalAncestry

makedocs(
    sitename = "LocalAncestry",
    format = Documenter.HTML(),
    modules = [LocalAncestry],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
