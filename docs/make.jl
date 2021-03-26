using Documenter, Bounds

makedocs(
    modules = [Bounds],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Arie Beresteanu",
    sitename = "Bounds.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/ArieBeresteanu/Bounds.jl.git",
    push_preview = true
)
