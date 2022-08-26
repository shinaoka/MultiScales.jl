using MultiScales
using Documenter

DocMeta.setdocmeta!(MultiScales, :DocTestSetup, :(using MultiScales); recursive=true)

makedocs(;
    modules=[MultiScales],
    authors="Hiroshi Shinaoka <h.shinaoka@gmail.com> and contributors",
    repo="https://github.com/shinaoka/MultiScales.jl/blob/{commit}{path}#{line}",
    sitename="MultiScales.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
