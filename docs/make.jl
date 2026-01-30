using TSConditioning
using Documenter

DocMeta.setdocmeta!(TSConditioning, :DocTestSetup, :(using TSConditioning); recursive=true)

makedocs(;
    modules=[TSConditioning],
    authors="Galen Lynch <galen@galenlynch.com>",
    sitename="TSConditioning.jl",
    format=Documenter.HTML(;
        canonical="https://galenlynch.github.io/TSConditioning.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/galenlynch/TSConditioning.jl",
    devbranch="main",
)
