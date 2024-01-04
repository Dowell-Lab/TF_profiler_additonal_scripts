This was pulled from cistromeDB http://cistrome.org/db/#/ batch download on 12/1/22.
Please see cistrome db website for citation instructions.
Raw data is in human_factor.tar.gz, described in human_factor_full_QC.txt.
Using the qc metrics described in qc_explained.pdf I parsed out data in parsing_cistrome.ipynb.
The gist of the parsing is that cistromeDB provides 6 main scoring metrics, 3 mapping metrics and 3 peak metrics.
If there was a 0 in mapping or a 0 in peak the dataset was not kept. The max score is 3 meaning that all 3 scoring metrics were met in that category.
The bed files were renamed using the following set up:
<cistromDB id>_<tf>_<cell line>_<mapping score>_<peak score>
This was repeated for mouse on 12/2/22.
