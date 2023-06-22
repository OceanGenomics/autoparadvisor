#!/usr/bin/perl
use strict;
use warnings;
use threads;
use threads::shared;
use Data::Dumper;
use YAML qw(LoadFile);
#use YAML::XS 'LoadFile';
use Cwd qw(abs_path);
use Data::Dumper;
use feature 'say';

#get current path, read in yaml file
my $path = abs_path();
my $filename = $path . ('/scallop.yml'); #change to config.yml later
my $yaml = LoadFile($filename);
my %yaml_data = %{ LoadFile($filename) };
#print %yaml_data;


#get parameter values/step/type
my @paranames;
our %parameter_values;
our %step_size;
our %type;

my $parameter_bounds = $yaml->{parameter_bounds};
foreach my $i (@{$parameter_bounds}){
    while(my ($key, $value) = each %{$i}){
        #add values/step/type
        $parameter_values{$key} = $value->{default};
        $type{$key} = $value->{type};
    }
    push @paranames, keys %{$i};
}


my $scallop_path = 'biodb/';
my $experement_file = 'input.gtf';
my $lib_type = 'empty';
my $command = "${scallop_path}scallop -i $experement_file --verbose 0 --min_transcript_coverage 0 --library_type $lib_type\n";
print $command;

#get paras
my $name = $yaml->{initial_command}->{name};
my $input_command = $yaml->{initial_command}->{input_command};
my $output_command = $yaml->{initial_command}->{output_command};
my $additional_command = $yaml->{initial_command}->{additional_command};

my $command_test = "${scallop_path}$name$input_command$experement_file$additional_command";
print $command_test;

#test from original file 
#ultimate goal is to reimplement these dicts 
# As CA progresses this hash also keeps track of the best parameter vector.
our %parameter_values_test = (
  "min_flank_length" => 3,
  "max_edit_distance" => 10,
  "min_bundle_gap" => 50,
  "min_num_hits_in_bundle" => 20,
  "min_mapping_quality" => 1,
  "min_splice_boundary_hits" => 1,
  "use_second_alignment" => 0, #true/false
  "uniquely_mapped_only" => 0, #true/false
  "min_subregion_gap" => 3,
  "min_subregion_overlap" => 1.5, #real value
  "min_subregion_length" => 15,
  "max_intron_contamination_coverage" => 2.0, #real value
  "min_transcript_length_base" => 150,
  "min_transcript_length_increase" => 50,
  "min_exon_length" => 20,
  "max_num_exons" => 1000,
  "max_dp_table_size" => 10000,
  "min_router_count" => 1
);

# This is the initial step size for each parameter, this will be slowly decreased
our %step_size_test = (
  "min_flank_length" => 10,
  "max_edit_distance" => 20,
  "min_bundle_gap" => 100,
  "min_num_hits_in_bundle" => 100,
  "min_mapping_quality" => 10,
  "min_splice_boundary_hits" => 10,
  "use_second_alignment" => 1, #true/false
  "uniquely_mapped_only" => 1, #true/false
  "min_subregion_gap" => 10,
  "min_subregion_overlap" => 10, #real value
  "min_subregion_length" => 20,
  "max_intron_contamination_coverage" => 10, #real value
  "min_transcript_length_base" => 500,
  "min_transcript_length_increase" => 100,
  "min_exon_length" => 100,
  "max_num_exons" => 10000,
  "max_dp_table_size" => 100000,
  "min_router_count" => 10
);

# This keeps track of the parameter type when decreasing the step size and stop conditions
our %type_test = (
  "min_flank_length" => "int",
  "max_edit_distance" => "int",
  "min_bundle_gap" => "int",
  "min_num_hits_in_bundle" => "int",
  "min_mapping_quality" => "int",
  "min_splice_boundary_hits" => "int",
  "use_second_alignment" => "bool", #true/false
  "uniquely_mapped_only" => "bool", #true/false
  "min_subregion_gap" => "int",
  "min_subregion_overlap" => "float", #real value
  "min_subregion_length" => "int",
  "max_intron_contamination_coverage" => "float", #real value
  "min_transcript_length_base" => "int",
  "min_transcript_length_increase" => "int",
  "min_exon_length" => "int",
  "max_num_exons" => "int",
  "max_dp_table_size" => "int",
  "min_router_count" => "int"
);

my $ref_file = 'C:/Users/Zhiwen Yan/Downloads/gencode.v35.annotation.gtf';
my $get_transcript = 'cat ' + $ref_file + ' | awk \'{print $3}\' | grep -c transcript';
my $num_transcripts = `$get_transcript`;

print 'number of transcripts:' + $num_transcripts;