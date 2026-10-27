# latexmk configuration for the NE 806 notes.
#
#   latexmk                 draft build (shows transcription markup) -> main.pdf
#   latexmk -r student.rc   student build (hides it) -> build/student/main.pdf
#   latexmk -c              remove auxiliary files
$pdf_mode = 1;
@default_files = ('main.tex');
$bibtex_use = 2;
$clean_ext = 'bbl run.xml';
