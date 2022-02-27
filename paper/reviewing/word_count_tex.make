DingelNeiman_wordcount.tex: DingelNeiman.tex
	cat DingelNeiman.tex | \
	sed '/begin{equation}/,/end{equation}/d' | \
	sed '/begin{equation\*}/,/end{equation\*}/d' | \
	sed '/begin{align}/,/end{align}/d' | \
	sed '/begin{align\*}/,/end{align\*}/d' | \
	sed '/begin{figure}/,/end{figure}/d' | \
	sed '/begin{table}/,/end{table}/d' | \
	sed '/begin{abstract}/,/end{abstract}/d' | \
	grep -v 'maketitle' | \
	sed 's/^\\bibliography{/\\nobibliography{/' > $@
