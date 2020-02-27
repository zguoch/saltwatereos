
figpath=source/_figures

function makelatexfigures()
{
    cd figures/latex
    # file tree of main program
    xelatex filetree_main.tex
    pdf2svg filetree_main.pdf filetree_main.svg
    mv filetree_main.svg ../../$figpath
    mv filetree_main.pdf ../../$figpath

    # mv out to the original path
    cd ../../
}

# latex figures
makelatexfigures