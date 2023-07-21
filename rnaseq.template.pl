#!/usr/bin/perl
use strict;
use warnings;

use CQS::ClassFactory;
use CQS::FileUtils;
use CQS::SystemUtils;
use CQS::ConfigUtils;
use CQS::PerformRNAseq;

my $def = {

  #define task name, this name will be used as prefix of a few result, such as read count table file name.
  task_name => "%prj_name",

  #email which will be used for notification if you run through cluster
  email => "hua-chang.chen\@vumc.org",
  emailType => "FAIL",

  #target dir which will be automatically created and used to save code and result, you need to change it for each project.
  target_dir         => "%pw_prj/result",

  #DEseq2 fold change, you can use either 1.5 or 2 or other option you want.
  DE_fold_change     => 1.5,

  #Since we have most of figure, we don't need multiqc anymore. But if you want, you can set it to 1.
  perform_multiqc    => 0,

  #We use webgestalt to do gene enrichment analysis using differential expressed genes.
  perform_webgestalt => 1,

  #We use GSEA for gene set enrichment analysis. It works for human genome only.
  perform_gsea       => 1,

  #If we need to trim the adapter from reads. Set to 0 if you don't find adapter in raw data.
  perform_cutadapt => 1,

  #If you find adapter through FastQC report, you need to specify adapter here.
  cutadapt_option  => "-O 1 -q 20 -a AGATCGGAAGAGCGTC -a AGATCGGAAGAGCACA -A AGATCGGAAGAGCACA -A AGATCGGAAGAGCGTC -a CTGTCTCTTATACACA -A CTGTCTCTTATACACA",

  #discard reads with length less than 30 after adapter trimming
  min_read_length  => 30,

  #Is the data pairend data or single end data
  pairend => 1,

  #Call variant using GATK pipeline
  perform_call_variants => 0,
  
  #source files, it's a hashmap with key (sample name) points to array of files. For single end data, the array should contains one file only.

  files => { %files

  },

  #group definition, group name points to array of sample name defined in files.
  groups => {
    "case"   => %case_group_name,
    %ctrl_str
  },

  #comparison definition, comparison name points to array of group name defined in groups.
  #for each comparison, only two group names allowed while the first group will be used as control.
  %compare_str
};

my $config = performRNASeq_gencode_hg19($def, 1);

1;
