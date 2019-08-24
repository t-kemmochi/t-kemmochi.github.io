function compile()

    open("table.tex","w") do io
        println(io,"\\begin{longtable}{|l|l|l|} \\hline")
    end

    open("table.tex","a") do io
        tbl = sortrows(readcsv("pronunciation.csv"))
        (row,col) = size(tbl)

        println(io,"単語 & 発音 & 備考 \\\\ \\hline\\hline \\endhead")

        for i=1:row
            word = tbl[i,1]
            pron = tbl[i,2]
            note = tbl[i,3]
            println(io,"$word & \\textipa{$pron} & $note \\\\\\hline")
        end

        println(io,"\\end{longtable}")
    end

    run(`ptex2pdf -l pronunciation.tex`)
    run(`ptex2pdf -l pronunciation.tex`)
    run(`SumatraPDF pronunciation.pdf`)

    run(`rm *.toc`)
    run(`rm *.log`)
    run(`rm *.aux`)

end


function preview(input)
    filename = randstring()
    open("$filename.tex","a") do io
        println(io,"\\documentclass[12pt,a4j]{jsarticle}")
        println(io,"\\usepackage[T1]{fontenc}")
        println(io,"\\usepackage{lmodern}")
        println(io,"\\usepackage[T1]{tipa}")
        println(io,"\\pagestyle{empty}")
        println(io,"\\begin{document}")
        println(io,"\\Huge")
        println(io,"\\textipa{$input}")
        println(io,"\\end{document}")
    end
    run(`ptex2pdf -l $filename.tex`)
    run(`pdfcrop $filename.pdf $filename.pdf`)
    run(`SumatraPDF $filename.pdf`)
    run(`rm $(filename)*`)
end


function add(word,pron,note)
    open("pronunciation.csv","a") do io
        println(io,"$word,$pron,$note")
    end
end


function rewrite(word,note,pron="")
    tbl = readcsv("pronunciation.csv")
    rows = find(tbl[:,1] .== word)
    if length(rows) == 0
        error("There is no '$word' in the list!")
    elseif length(rows) > 1
        error("There is many '$word's in the list!")
    else
        row = rows[1]
    end
    if length(note)>0
        tbl[row,3] = note
    end
    if length(pron)>0
        tbl[row,2] = pron
    end
    writecsv("pronunciation.csv",tbl)
end


function git()
    run(`git add pronunciation.pdf`)
    run(`git commit pronunciation.pdf -m "pronunciation.pdf has been commited via Julia."`)
    # cd("../../")
    run(`git push`)
    # cd("./memo/english")
end


"""
# Commands for phonetic symbols:

ŋ -> N, θ -> T, ð -> D, ʃ -> S

ʒ -> Z, ʊ -> U, ɔ -> O, ə -> @

ʌ -> 2, æ -> \\ae, ɑ -> A, ɪ -> I

first accent -> \\', second accent -> \\`

To input \\, you should type \\ twice!!

# Functions

## add(word,pron,note)
Add word and others to the list.

## preview(input)
Preview the phonetic symbols.

## compile()
Compile the tex file.

## rewrite(word,note,pron="")
Rewire the note (and optionally symbols) for the 'word'.

## git()
Add, commit, and push the pdf file to git.

"""
function readme()
end


println("")
println("This is pronunciation.jl. You can type '?readme' for help.")
