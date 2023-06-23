#!/usr/bin/perl
use strict;
use warnings;
use threads;
use threads::shared;
use Data::Dumper;
use YAML::XS 'LoadFile';
use Cwd qw(abs_path);


# This directory where the files created during the CA process will be written and read.
our $working_dir = shift;
chomp $working_dir;

# This is the input BAM file of aligned reads
our $experement_file = shift;
chomp $experement_file;

# This is the input Library type of aligned reads
my $lib_type = "empty";
$lib_type = shift;
chomp $lib_type;

# This is the reference type of aligned reads
my $ref_file = "/biodb/human/gencode/v35/gene_annotations.gtf";
$ref_file = shift;
chomp $ref_file;

#print "$ref_file";

# This controls the number of parallel steps to take in each direction,
# it will actually run twice this number of processes because its that many steps
# in each direction.
my $num_threads = 3;
my $in_threads = shift;
$num_threads = $in_threads if $in_threads ne "";

#print "$num_threads";

# max iteration of coordinate ascent
my $max_iters = 0;
my $in_iters = shift;
$max_iters = $in_iters if $in_iters ne "";

#print "$max_iters";

my $scallop_path = "";
$scallop_path = shift;
chomp $scallop_path;

my $sub_sample = 1;
my $new_sub_sample = shift;
$sub_sample = $new_sub_sample if $new_sub_sample ne "";

# This is used if you already ran this once and changed something,
# this will be used to check to see if you don't need to rerun a parameter choice vector
my $old_working_dir = "";
$old_working_dir = shift;
chomp $old_working_dir;

#adapted with yaml config file
#get current path, read in yaml file
my $path = abs_path();
my $filename = $path .'/..'. ('/scallop.yml'); #change to config.yml later
my $yaml = LoadFile($filename);

#get paras
our $name = $yaml->{initial_option}->{name};
our $input_command = $yaml->{initial_option}->{input_option};
our $output_command = $yaml->{initial_option}->{output_option};
our $additional_command = $yaml->{initial_option}->{additional_option};

my $gffcompare_path = $yaml->{gffcompare}->{directory};
my $gffcompare_command = $yaml->{gffcompare}->{command};
my $gtfcuff_path = $yaml->{gtfcuff}->{directory};
my $gtfcuff_command = $yaml->{gtfcuff}->{command};

#get parameter values/step/type/prefix
our %parameter_values;
our %step_size;
our %type;
our %command;

# These are the default parameter values, the place to start the CA
# each of the elements is a tuneable parameter for Scallop.
my $parameter_bounds = $yaml->{parameter_bounds};
foreach my $i (@{$parameter_bounds}){
    while(my ($key, $value) = each %{$i}){
        # As CA progresses this hash also keeps track of the best parameter vector
        $parameter_values{$key} = $value->{default};
        # This is the initial step size for each parameter, this will be slowly decreased
        $step_size{$key} = $value->{step};
        # This keeps track of the parameter type when decreasing the step size and stop conditions
        $type{$key} = $value->{type};
        # This keeps track of the parameter sign ' -' or ' --'
        $command{$key} = $value->{command};
    }
}



# This subroutine runs Scallop using the parameter vector stored in %parameter_values
# with one parameter changed.
# If both arguments are "", it will run the full vector (i.e. the defaults on the first run).
# Returns the AUC of this new parameter vector.
# $check is used to make a first pass to see if Scallop needs to be run, when equal to "check"
#  returns 1 if if the AUC has already been computed, and 0 if Scallop would have been run.
sub run_with_one_change{
  my $param_to_change = shift;
  my $param_value = shift;
  my $check = "";
  $check = shift;

  # "--min_transcript_coverage 0" sets allows the AUC to work on the quality threshold output by Scallop
  my $command = "${scallop_path}$name$input_command$experement_file$additional_command";
  my $out_fname = "";
  #this only allow boolen case with 0,1
  return 0 if($param_value ne "" && (($type{$param_to_change} eq "bool" && $param_value > 1) || $param_value < 0));

  # ensures the floating point parameters don't get too complex
  if($type{$param_to_change} eq "float"){
    $param_value = sprintf("%.2f", $param_value);
  }

  # builds the command line and output file name from the parameter vector, including the parameter to change.
  for my $p (sort keys(%type)){
    if($p ne $param_to_change){
      if($type{$p} ne "cag"){
        $command .= "$command{$p}$p $parameter_values{$p} ";
      }else{
        #this only allow true/ false
        $command .= "$command{$p} ".(($parameter_values{$p}==1)?"true":"false");
      }
      $out_fname .= "_$parameter_values{$p}";
    }else{
        if($type{$p} ne "cag"){
          $command .= "$command{$p} $param_value ";
        }else{
          $command .= "$command{$p} ".(($param_value==1)?"true":"false");
        }
        $out_fname .= "_$param_value";
    }
  }
  if($sub_sample < 1){
    $command .= " --subsampling $sub_sample";
  }
  my $auc = "";
  print STDERR "Would run '$command'\n";
  #return 0;
  system("mkdir -p $working_dir");

  #TODO: leave this two conditions unchanged
  # if AUC is already computed return it rather than rerunning Scallop
  if(-e "$working_dir/$out_fname.auc" && `grep -c auc $working_dir/$out_fname.auc` > 0){
    $auc = `cat $working_dir/$out_fname.auc`;

  # if AUC is in the old directory, copy these files into the current directory and return the AUC
  }elsif($old_working_dir ne "" && -e "$old_working_dir/$out_fname.auc" && `grep -c auc $old_working_dir/$out_fname.auc` > 0){
    system("cp $old_working_dir/$out_fname* $working_dir/");
    $auc = `cat $working_dir/$out_fname.auc`;

  # Otherwise, something needs to be run.
  }else{
    # if you get to this point and you're just checking, you would run Scallop, so something needs to be done.
    return 0 if $check eq "check";

    # either extract the compressed GTF file, or run Scallop
    system("gunzip $working_dir/$out_fname.gtf.gz") if(-e "$working_dir/$out_fname.gtf.gz");
    system("$command -o $working_dir/$out_fname.gtf >/dev/null 2>&1 ") if(!(-e "$working_dir/$out_fname.gtf"));

    if((-e "$working_dir/$out_fname.gtf")){

      # make sure the output formatting from the aligner matches the reference
      #if(`grep -c "^chr" $working_dir/$out_fname.gtf` == 0){
      #  system("sed -i 's/^/chr/' $working_dir/$out_fname.gtf");
      #}

      # Compute AUC
      #TODO: change hard code
      if(`grep -c "^chr" $working_dir/$out_fname.gtf` > 0 && `grep -c "^chr" $ref_file` == 0){
        system("sed -i 's/^chr//' $working_dir/$out_fname.gtf");
      }
      system("$gffcompare_path$gffcompare_command$ref_file -o $working_dir/$out_fname $working_dir/$out_fname.gtf");
      if((-e "$working_dir/$out_fname.$out_fname.gtf.tmap")){
        # $num_transcripts should match the number of transcripts in the reference
        #if($ref_type eq "GRCh37"){
        #  my $num_transcripts = ;
        #}else{
        #  my $num_transcripts = ;
        #}
        #get number of reference transcript for future gftcuff
        my $get_transcript = `cat $ref_file | awk \'{print \$3}\' | grep -c transcript`;
        my $num_transcripts = `$get_transcript`;
        print 'number of transcripts:'.$num_transcripts;
        #my $num_transcripts = 229580;

        $auc = `$gtfcuff_path$gtfcuff_command$working_dir/$out_fname.$out_fname.gtf.tmap $num_transcripts | tee $working_dir/$out_fname.auc`;
        system("rm $working_dir/$out_fname.$out_fname.gtf.refmap $working_dir/$out_fname.loci $working_dir/$out_fname.annotated.gtf $working_dir/$out_fname.tracking");
        system("gzip $working_dir/$out_fname.gtf");
      }else{
        # the tmap file wasn't created, this is an error.
        return 0;
      }
    }else{
      # the GTF file wasn't created, this is an error.
      return 0;
    }
  }

  #if you got here and $check is set, then nothing needs to be done.
  return 1 if $check eq "check";
  chomp $auc;
  $auc =~ s/.*auc = //;
  return $auc;
}

# gets default AUC
my $cur_auc = run_with_one_change("","","");

# loops as long as you were able to decrease one of the step sizes
my $decreased_steps = 1;
while($decreased_steps == 1){
  $decreased_steps = 0;
  my $made_one_change = 1;

  # loops as long as you have made a change in the parameter vector
  # without a change in step size
  while($made_one_change == 1){
    $made_one_change = 0;
    foreach my $param(sort keys(%type)){
      print("$cur_auc\n");

      # loops as long as you can continue moving on this parameter (coordinate) and still increase AUC
      my $single_param_change = 1;
      while($single_param_change==1){
        $single_param_change = 0;

        my $iter_num = `ls $working_dir | grep .gtf.gz | wc -l`;
        print "$iter_num";
        if($max_iters > 0 && $iter_num > $max_iters){
            exit 0;
        }

        #no-thread results are those which nothing needs to be run, avoids spawning threads when nothing needs to be computed
        # threads are spawned anytime any of Scallop, gtfcuff, or gffcompare need to be run.
        print STDERR "Updating $param, type: ". $type{$param} ."\n";
        if($type{$param} ne "bool"){
          my @threads;
          my @threads_index;
          my @nothreads;
          my @nothreads_index;
          foreach my $t(1...$num_threads){
            if(run_with_one_change($param, $parameter_values{$param} + ($t * $step_size{$param}),"check") == 0){
              push @threads, threads->new (sub { return run_with_one_change($param, $parameter_values{$param} + ($t * $step_size{$param}),""); } );
              push @threads_index, $t;
            }
            else{
              push @nothreads, run_with_one_change($param, $parameter_values{$param} + ($t * $step_size{$param}),"");
              push @nothreads_index, $t;
            }
          }
          foreach my $t(1...$num_threads){
            if(run_with_one_change($param, $parameter_values{$param} - ($t * $step_size{$param}),"check") == 0){
              push @threads, threads->new (sub { run_with_one_change($param, $parameter_values{$param} - ($t * $step_size{$param}),""); } );
              push @threads_index, -1 * $t;
            }else{
              push @nothreads, run_with_one_change($param, $parameter_values{$param} - ($t * $step_size{$param}),"");
              push @nothreads_index, -1 * $t;
            }
          }

          # join threads and update parameter vector
          print STDERR "Num threads running: " . scalar(@threads) . "\n";
          my $max_change = 0;
          foreach my $t (0...scalar(@threads)-1){
          	my $in_auc = $threads[$t]->join();
            if($in_auc > $cur_auc){
              $cur_auc = $in_auc;
              $max_change = $threads_index[$t];
              $single_param_change = 1;
              $made_one_change = 1;
            }
          }

          # join the non-run processes, but keep the farthest increasing step
          foreach my $nt (0...scalar(@nothreads)-1){
            my $in_auc = $nothreads[$nt];
            if($in_auc > $cur_auc || ($in_auc == $cur_auc && abs($max_change) > abs($nothreads_index[$nt]))){
              $cur_auc = $in_auc;
              $max_change = $nothreads_index[$nt];
              $single_param_change = 1;
              $made_one_change = 1;
            }
          }
          $parameter_values{$param} += $max_change * $step_size{$param};

        #threads are not used for boolean parameters since only one will ever need to be run.
        }else{
          my $auc_plus = run_with_one_change($param, $parameter_values{$param} + $step_size{$param},"");
          if($auc_plus > $cur_auc){
            $cur_auc = $auc_plus;
            $parameter_values{$param} += $step_size{$param};
            $single_param_change = 1;
            $made_one_change = 1;
          }else{
            my $auc_minus = run_with_one_change($param, $parameter_values{$param} - $step_size{$param},"");
            if($auc_minus > $cur_auc){
              $cur_auc = $auc_minus;
              $parameter_values{$param} -= $step_size{$param};
              $single_param_change = 1;
              $made_one_change = 1;
            }
          }
        }
      }
    }
  }
  last;

  # decrease step sizes as long as you can without the steps being too small (0.01 and 1 for float and int/bool)
  foreach my $param(sort keys(%type)){
    if($type{$param} eq "int" ){
      my $temp = int($step_size{$param} * 0.75);
      $temp = ($temp < $step_size{$param} - 1)?$temp:$step_size{$param} - 1;
      if($temp > 0){
        $step_size{$param} = $temp;
        $decreased_steps = 1;
      }
    }

    if($type{$param} eq "float"){
      my $temp = sprintf "%.2f", $step_size{$param} * 0.75;
      $temp = ($temp < $step_size{$param} - 0.01)?$temp:$step_size{$param} - 0.01;
      if($temp > 0){
        $step_size{$param} = $temp;
        $decreased_steps = 1;
      }
    }
  }
}
