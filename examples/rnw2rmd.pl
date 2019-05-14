#!/usr/bin/perl

## Based on rnw2rmd.pl by Laurent Gatto
## https://gist.github.com/lgatto/d9d0e3afcc0a4417e5084e5ca46a4d9e

## The script ignores labels and references, un-numbered sections
## (section*), quotes and probably a couple of more. It won't deal
## with the pre-amble, bibliography and document tags either. Still
## useful, though.

## Now it removes \begin{*} and \end{*} environments, and basically
## every LaTeX command that starts with a \. It also removes all
## trailing white spaces.

## NOTE: Some things here can be a little specific, so take care when
## using in your own documents. Depending on the document you are trying
## to convert, some residual text may be found (but should be easy to
## remove by hand).

## Usage:
## rnw2rmd file.Rnw > file.Rmd

use warnings;

my $filename = $ARGV[0];

open (RNW, $filename);
while (<RNW>) {
    s|\\Robject\{(.+?)}|`$1`|g;
    s|\\Rcode\{(.+?)}|`$1`|g;
    s|\\Rclass\{(.+?)}|*$1*|g;
    s|\\Rfunction\{(.+?)}|`$1`|g;
    s|\\texttt\{(.+?)}|`$1`|g;
    s|\\textit\{(.+?)}|*$1*|g;
    s|\\textbf\{(.+?)}|**$1**|g;
    s|\\emph\{(.+?)}|*$1*|g;
    s|<<|```\{r |g;
    s|>>=|}|g;
    s|@|```|g;
    s|\\section\{(.+?)}|# $1|g;
    s|\\section<presentation>\*\{\}||g;
    s|\\frametitle\{(.+?)}|# $1|g;
    s|\\subsection\{(.+?)}|## $1|g;
    s|\\subsubsection\{(.+?)}|### $1|g;
    s|\\Biocexptpkg\{(.+?)}|`r Biocexptpkg("$1")`|g;
    s|\\Biocannopkg\{(.+?)}|`r Biocannopkg("$1")`|g;
    s|\\Biocpkg\{(.+?)}|`r Biocpkg("$1")`|g;
    s|\\cite\{(.+?)}|[\@$1]|g;
    s|\\ref\{(.+?)}|\\\@ref($1)|g;
    s|\\url\{(.+?)}|<$1>|g;
    s|\\ldots|\.\.\.|g;
    # s|\\label\{| \{#|g; ## only for sections
    s|\\label\{(.+?)}||g;
    s|\\item(.+?)|-$1|g;
    s|\\begin\{table\}\[\]||g;
    s|^\\begin\{tabular\}\{(.+?)}||g;
    s|\\hline||g;
    ## Remove \begin{*}[*]
    s|\\begin\{(.+?)}?\[*[A-Za-z]+\]||g;
    ## Remove \begin{*}
    s|\\begin\{(.+?)}||g;
    s|\\end\{(.+?)}||g;
    s|\\documentclass\{(.+?)}||g;
    s|\\usepackage\{(.+?)}||g;
    s|\\usepackage\[*[A-Za-z0-9]+\]\{(.+?)}||g;
    s|\\setlength\{(.+?)}\{(.+?)}||g;
    s|\\SweaveOpts\{(.+?)}||g;
    s|\\frame\{(.+?)}||g;
    s|\\setkeys\{Gin\}\{(.+?)}||g;
    s|\\centering||g;
    s|\\includegraphics\{(.+?)}||g;
    s|\\caption\{(.+?)}||g;
    ## Try to remove comments
    s|^(%+)*||g;
    ## Try to remove basically everything starting with \
    s|^\\.+||g;
    ## Remove all trailing whitespace before every line
    s|^(\s+)||g;
    ## Substitute lines with = (rather specific) by a comment
    s|^=+|<!-- -->\n|g;
    print $_;
}
close (RNW);
