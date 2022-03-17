#!/bin/bash

#----------
# Defaults
#----------

log_file=''
output_dir='.'
prefix='Mulberry_report'
quantitative='FALSE'
report_template=''
template='report_template.Rnw'
title=''


#-------
# Usage
#-------

display_usage() {
  echo -e "\nUsage: $0 -o 'output_dir' -p 'prefix' -t 'report_title'\n"
  echo -e "  -a  analysis results from mulberry."
  echo -e "      Default: ar_prediction.tsv"
  echo -e "  -l  Log file."
  echo -e "      Default: <prefix>.log"
  echo -e "  -m  message for inclusion in the results section."
  echo -e "      Default: ''"
  echo -e "  -o  output directory."
  echo -e "      Default: ."
  echo -e "  -t  Report title."
  echo -e "      Default: ''"
  echo -e "  -T  Report template (.Rnw)."
  echo -e "      Default: report_template.Rnw"
  echo -e "  -h  Displays this help message\n"
  echo -e "create_mulberry_report.sh"
  echo -e "This script is part of an AMR gene analysis pipeline and produces "
  echo -e "quantitative and graphical reports.\n"
  exit
}


#---------
# Options
#---------

while getopts ":a:l:m:o:p:t:T:h" opt; do
  case $opt in
    a) ar_prediction=$OPTARG;;
    l) log_file=$OPTARG;;
    m) message=$OPTARG;;
    o) output_dir=${OPTARG%/};;
    p) prefix=$OPTARG;;
    t) title=$OPTARG;;
    T) report_template=$OPTARG;;
    h) display_usage;;
   \?) #unrecognized option - show help
      echo -e \\n"Option -$OPTARG not allowed."
      display_usage;;
  esac
done

if [ -z $log_file ]; then
  log_file="${output_dir}/create_mulberry_report.log"
fi


#----------------------
# Create report
#----------------------

sweave_report="${output_dir}/${prefix}.Rnw"
tex_report="${output_dir}/${prefix}.tex"
cp $report_template $sweave_report

#module load R

echo "Converting Sweave to Tex" | tee $log_file

cd ${output_dir}

R -e "Sweave('$sweave_report')" --args  --ar_prediction="$ar_prediction" --message="$message"
exit_status=$?

if [ "$exit_status" -ne 0 ]; then
  echo "Sweave was unable to complete the conversion. Please check the log file." \
    | tee -a $log_file
  exit 1
fi

#----------------------------------------------------
# Insert title into tex report
#----------------------------------------------------
echo $title
title_clean=$(echo $title | sed 's/"//g')
echo $title_clean
sed -i['.bak'] "s/\"Title goes here\"/ ${title_clean}/" $tex_report || true

echo "Tex conversion successful." | tee -a $log_file

#---------------
# Create report
#---------------

echo "Converting Tex to Pdf" | tee -a $log_file
#module load latex
# Convert it twice to avoid problems
R CMD pdflatex -interaction=nonstopmode $tex_report 2>&1 >/dev/null
R CMD pdflatex -interaction=nonstopmode $tex_report >>$log_file 2>&1
exit_status=$?

if [ "$exit_status" -ne 0 ]; then
  echo "The tex to pdf conversion had a non-zero exit status. Please check the pdf file" \
  " and the log file." | tee -a $log_file
  exit 1
fi

echo "Report creation successful." | tee -a $log_file
