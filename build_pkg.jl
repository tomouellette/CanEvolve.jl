t = Template(; 
    user="tomouellette",
    authors="Tom W Ouellette",
    julia=v"1.5.4",
    plugins=[
        License(; name="MIT"),
        Git(; manifest=true, ssh=true),
        TravisCI(; linux = true, osx = true, x64 = true, x86 = true),
    ],
)

generate("CanEvolve.jl",t)