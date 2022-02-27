hyperlinkscheck:
	diff --side-by-side --suppress-common-lines <(grep -o --no-filename 'hypertarget{[A-Za-z0-9:_]*}' ../sections/*.tex | sed 's/hypertarget//' | sort | uniq) <(grep -o --no-filename 'hyperlink{[A-Za-z0-9:_]*}' ../sections/*.tex | sed 's/hyperlink//' | sort | uniq)
