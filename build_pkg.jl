Template(; 
    user="tomouellette",
    authors="Tom W Ouellette",
    julia=v"1.5.4",
    plugins=[
        License(; name="MIT"),
        Git(; manifest=true, ssh=true),
        GitHubActions(),
    ],
)

generate("CanEvolve.jl",t)