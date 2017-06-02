#!/usr/bin/perl -w
#Author: Thomas Thiel (thiel@ipk-gatersleben.de)
#v 0.6
#
#use LWP::Simple;
#require 'dumpvar.pl';

# Modified by Junli Zhang (zhjl86@gmail.com) on 05/14/2016
# Now the script check all the Restriction enzymes in the RE file
# so the 3rd parameter can be any enzyme name
# Usage: ./SNP2CAPS.pl for_SNP2CAPS.fa REgcg.txt EcoRV,BamHI > CAPS_output.txt

#######################
##### DECLARATION #####
#######################

# INTERNATIONAL UNION OF BIOCHEMISTRY AND MOLECULAR BIOLOGY (IUBMB)
# Nomenclature for Incompletely Specified Bases in Nucleic Acid Sequences
# http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html

my %bas_amb = (
    'A' => 'A-*',
    'C' => 'C-*',
    'G' => 'G-*',
    'T' => 'T-*',
    'm' => '[ac]-*',    # A or C
    'r' => '[ag]-*',    # A or G
    'w' => '[at]-*',    # A or T
    's' => '[cg]-*',    # C or G
    'y' => '[ct]-*',    # C or T
    'k' => '[gt]-*',    # G or T
    'v' => '[acg]-*',   # A or C or G (not T)
    'h' => '[act]-*',   # A or C or T (not G)
    'd' => '[agt]-*',   # A or G or T (not C)
    'b' => '[cgt]-*',   # C or G or T (not A)
    'x' => '[acgtn]-*', # G or A or T or C - not defined
    'n' => '[acgtn]-*'  # G or A or T or C
     );

my %bas_amb_fuzzy = ('A'=>'[AN]-*','C'=>'[CN]-*','G'=>'[GN]-*','T'=>'[TN]-*',
    'm'=>'[acN]-*','r'=>'[agN]-*','w'=>'[atN]-*','s'=>'[cgN]-*','y'=>'[ctN]-*',  
    'k'=>'[gtN]-*','v'=>'[acgN]-*','h'=>'[actN]-*','d'=>'[agtN]-*','b'=>'[cgtN]-*',
    'n' => '[acgtn]-*');

my $rebase_file;
my %enzyme; # stores REBASE restriction enzyme data
my %rebase_customer; # stores REBASE customer info
my %marker; # stores sequence data to each marker
my @candidate; # stores CAPS marker candidates


#######################
##### HELP/SYNTAX #####
#######################

@ARGV == 3 || die (
"\n_______________________________________________________________________________

SNP2CAPS
========

A SNP and INDEL analysis tool for CAPS marker development

AUTHOR:  Thomas Thiel

SYNTAX:  (1) GUI version: SNP2CAPS.pl -gui

         (2) Command line version:

             SNP2CAPS.pl <msa_file> <rebase_db> <enzymes> > <output>
             
             <msa_file>  file containing the multiple alignments
             <rebase_db> file containing the REBASE database in GCG format
             <enzymes>   comma separated list of enzymes (e.g. EcoRV,BamHI)
             <output>    save the results in the output file (optional)

see http://pgrc.ipk-gatersleben.de/SNP2CAPS/ for updates

_______________________________________________________________________________\n\n");


###################
##### M A I N #####
###################


GetSequences2();
GetRebaseData('f');
#GetCandidates(split ',', $ARGV[2]);
GetCandidates(keys %enzyme); # use all the enzymes
Results();



#######################
##### SUBROUTINES #####
#######################

#>>> READ DNA sequences from file <<<#

sub GetSequences2
  {
  $file = $ARGV[0];
  %marker = ();
  open (IN, "<$file") || die ("Could not open $file\n");
  my $msa_type; # type of alignment file
  while (<IN>)
    {
    if (/^>[^_]+_[^_]+/) {$msa_type = 'S2C';last}
    elsif (/^MSF/)  {$msa_type = 'MSF';last}
    elsif (/^CLUSTAL W/)  {$msa_type = 'ALN';last}
    };
  my ($count_sequences, $count_marker);
  seek (IN, 0, 0);
  if ($msa_type eq 'S2C') # SNP2CAPS format
    {
    $/ = ">";
    while (<IN>)
      {
      my ($id, $seq);
      next unless (($id, $seq) = /(.*?)\n(.*)/s);
      my ($marker_name, $cultivar) = $id =~ /(.*?)_(.*)/;
      $seq =~ s/[\s>]//g; # remove whitespace
      $count_sequences++;
      $marker{$marker_name}{$cultivar} = [$seq];
      push @{$marker{$marker_name}{$cultivar}}, FindGaps(\$seq);
      };
    $count_marker = scalar keys %marker;
    }
  elsif ($msa_type eq 'MSF')
    {
    $/ = "\n";
    my @parts = split '//', join '', <IN>; # $parts[2] contains sequence alignment
    foreach (split "\n", $parts[2])
      {
      my ($cultivar, $seq);
      next unless (($cultivar, $seq) = /^(\S+)\s+(\S+)/);
      $marker{$file}{$cultivar}[0] .= $seq;
      };
    foreach my $cultivar (keys %{$marker{$file}})
      {
      push @{$marker{$file}{$cultivar}}, FindGaps(\$marker{$file}{$cultivar}[0]);
      };
    $count_marker = 1;
    $count_sequences = scalar keys %{$marker{$file}};
    }
  elsif ($msa_type eq 'ALN') # CLUSTAL W format
    {
    $/ = "\n";
    while (<IN>)
      {
      chomp;
      next if /^CLUSTAL W/;
      my ($cultivar, $seq);
      next unless (($cultivar, $seq) = /^(\S+)\s+(\S+)/);
      $marker{$file}{$cultivar}[0] .= $seq;
      };
    foreach my $cultivar (keys %{$marker{$file}})
      {
      push @{$marker{$file}{$cultivar}}, FindGaps(\$marker{$file}{$cultivar}[0]);
      };
    $count_marker = 1;
    $count_sequences = scalar keys %{$marker{$file}};
    }
  else {return};
  $/ = "\n";
  close (IN);
  if ($ARGV[0] eq "-gui")
    {
    $info -> insert("insert","\nLoaded data of $count_marker marker(s) and $count_sequences sequence(s).");
    $info -> yviewMoveto(1);
    };
  };

sub FindGaps
  {
  my ($seqref) = @_;
  my @gaps; # stores gap position and type
  while ($$seqref =~ /([-_\.]+)/g) # find out positions of gaps
    {
    my $gap_pos = pos ($$seqref) + 1 - length $1;
    push @gaps, $gap_pos, length $1;
    };
  return @gaps
  };

#>>> READ REBASE enzyme data <<<#

sub GetRebaseData
  { # Option 'w' -> obtain via WWW, 'f' -> Read from File
  my ($rebase_file, $rebase_version, $rebase_reldate, $count);
  if ($_[0] eq 'w') {while (!($rebase_file = get("http://rebase.neb.com/rebase/link_gcg"))) {}}
  else
    {
    $file = $ARGV[1];
    open (IN, "<$file");
    $rebase_file = join '',<IN>;
    close (IN)
    };
  %enzyme = ();
  # EXTRACT items from the REBASE data file #
  (my $rebase_file1, my $rebase_file2) = split /\.\.\n/,$rebase_file;  # ".." separates intro and enzyme list
  ($rebase_version) = $rebase_file1 =~ /\nREBASE version (\d+)/s;
  ($rebase_reldate) = $rebase_file1 =~ /\nRich Roberts\s+(.*?)\n/s;
  (%rebase_customer) = $rebase_file1 =~ /\n\s+([A-Z])\s+(.*?)(?=\n)/gs;

  while ($rebase_file2 =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)(?: \? !| !)\s+(\S+|)\s+>(\S+|)/gs)
    # $1->name  $2->cutting pos.  $3->recognition seq.  $4->overhang
    # $5->isoschizomers  $6->commercial sources
    {
    my $isoschizomer = 0; # isoschizomer flag (0->false, 1-true)
    my $recseq = $3;
    my $enzname = $1;
    {
      $recseq =~ s/[_']//g; # remove cut position markers
      $isoschizomer = 1 if $enzname =~ s/;//; # preceding ; marks isoschizomers
      $recseq =~ s/N+$//i;
    };   # curly brackets keep matched variables
    (my $recseq_regex, my $recseq_fuzzy) = ($recseq, $recseq);
    foreach (keys %bas_amb) {$recseq_regex =~ s/$_/$bas_amb{$_}/g}; # translate to regex like style
    foreach (keys %bas_amb_fuzzy) {$recseq_fuzzy =~ s/$_/$bas_amb_fuzzy{$_}/g}; # for enzymes which are listed twice, e.g. AloI
    if (exists $enzyme{$enzname}) {$enzname .= ' !'} else {$count++};
    $enzyme{$enzname} = [$recseq_regex, $recseq, $isoschizomer, length ($recseq), $2, $4, $5, $6, $recseq_fuzzy];
    };
  };


#>>> SEARCH for CAPS candidates <<<#

sub GetCandidates
  { # call by enzymes to be used for screening
  if (scalar keys %enzyme == 0) {$info->insert("insert","\nError: Load REBASE data base first"); $info->yviewMoveto(1); return};
  if (scalar keys %marker == 0) {$info->insert("insert","\nError: Load Sequence Alignment first"); $info->yviewMoveto(1); return};
  @candidate = ();
  foreach my $enzyme_name(@_)
    {
    my ($not_palindromic, $rc_recseq, $rc_recseq_regex, $rc_recseq_fuzzy);
    study $enzyme{$enzyme_name}[0];   # of any benefit ?
    if (! Palindromic($enzyme{$enzyme_name}[1]) && ! (exists $enzyme{$enzyme_name.' !'} || $enzyme_name =~ /!/)) # test if palindromic or of type as AloI
      {
      $not_palindromic = 1;
      $rc_recseq = reverse DNAComplement($enzyme{$enzyme_name}[1]);
      ($rc_recseq_regex, $rc_recseq_fuzzy) = ($rc_recseq, $rc_recseq);
      foreach (keys %bas_amb) {$rc_recseq_regex =~ s/$_/$bas_amb{$_}/g}; # translate to regex like style
      foreach (keys %bas_amb_fuzzy) {$rc_recseq_fuzzy =~ s/$_/$bas_amb_fuzzy{$_}/g}
      }
    foreach my $marker_name (keys %marker)
      {
      my %cut; # stores enzyme binding positions temporarily
      foreach my $cultivar (keys %{$marker{$marker_name}})
        {
        while ($marker{$marker_name}{$cultivar}[0] =~ /($enzyme{$enzyme_name}[0])/ig)
          {
          pos ($marker{$marker_name}{$cultivar}[0]) += 1 - length $1;
          my $start = pos ($marker{$marker_name}{$cultivar}[0]);
          $cut{$start}{$cultivar} = ['f',$1];
          };
        if ($not_palindromic) # rec. seq. not palindromic => check the rc as well
          {
          while ($marker{$marker_name}{$cultivar}[0] =~ /($rc_recseq_regex)/ig)
            {
            pos ($marker{$marker_name}{$cultivar}[0]) += 1 - length $1;
            my $start = pos ($marker{$marker_name}{$cultivar}[0]);
            $cut{$start}{$cultivar} = ['r',$1];
            };
          };
        };
      my %test; # test for CAPS candidates
      # (1) true CAPS candidate (alternative restriction due to SNP)
      # (2) Not sure: probably true CAPS candidate (N within DNA sequence, but SNP at different position)
      # (3) Not sure: probably false positive CAPS candidate at at least one site (N within DNA sequence)
      # cases (2) and (3) need further inspection by the experimenter
      foreach my $site (keys %cut)
        {
        my $difference = scalar (keys %{$marker{$marker_name}}) - scalar (keys %{$cut{$site}});
        if ($difference > 0) # test if all sequences were cut
          {
          my %intersection; # intersection between 2 hashes
          @intersection{ keys %{$marker{$marker_name}} } = ();
          delete @intersection{ grep (exists $intersection{$_}, keys %{$cut{$site}}) };
          foreach my $cultivar (keys %intersection) # for all that were not cut
            {
            my $test_seq = substr($marker{$marker_name}{$cultivar}[0], $site - 1);
            $test_seq =~ s/((?:[^-]-*){$enzyme{$enzyme_name}[3]}).*/$1/;
            $cut{$site}{$cultivar} = ['n',$test_seq];
            if ($test_seq !~ /N/ig) {$test{1} = ""} # see above
            elsif ($test_seq !~ /$enzyme{$enzyme_name}[8]/i)
              {
              if ($not_palindromic) {if ($test_seq !~ /$rc_recseq_fuzzy/i) {$test{2} = ""}}
              else {$test{2} = ""}
              }
            else {$test{3} = ""}
            }
          }
        };
      if (%test) # true if alternative restriction was found
        {
        my $test_id = Test(\%test);
        foreach my $site (sort {$a <=> $b} keys %cut)
          {
          foreach my $cultivar (keys %{$cut{$site}})
            {
            my $cutsite = 0; # calculate the cut site position
            if ($cut{$site}{$cultivar}[0] eq "f") {$cutsite = $site - PrecedingGaps($marker_name, $cultivar, $site) + $enzyme{$enzyme_name}[4] - 1}
            elsif ($cut{$site}{$cultivar}[0] eq "r") {$cutsite = $site - PrecedingGaps($marker_name, $cultivar, $site) + $enzyme{$enzyme_name}[3] - $enzyme{$enzyme_name}[4] - 1}; ## gaps within rec.seq???? TEST
            push @{$candidate[$test_id]{$marker_name}{$enzyme_name}{$cultivar}}, $cut{$site}{$cultivar}[0], $cutsite, $site, $cut{$site}{$cultivar}[1];
            };
          };
        $candidate[$test_id]{$marker_name}{$enzyme_name} = group($marker_name,\%{$candidate[$test_id]{$marker_name}{$enzyme_name}});
        }
      }
    };
  my (%not_converted, @converted);
  @not_converted{keys %marker} = ();
  for (my $i = 1; $i < 4; $i++)
    {
    delete @not_converted{keys %{$candidate[$i]}};
    push @converted, keys %{$candidate[$i]}
    };
  @{$candidate[0]} = sort keys %not_converted;
  @converted = sort @converted;
  };


#>>> SORTING ROUTINE <<<#

sub group
  { # two keys are appended during this routine containing sorting matrizes
    # according to the presence / absence of restriction sites (group_align)
    # and according to the resulting fragment sizes (group_native)
  my $marker = shift;
  my %group = %{$_[0]};
  my (%sort_align, %sort_native); # contain the sorting keys
  my @cultivar = sort keys %group;
  for (my $i = 0; $i < scalar @cultivar; $i++)
    {
    for (my $j = 0; $j < scalar @{$group{$cultivar[$i]}}; $j += 4)
      {
      $sort_align{$cultivar[$i]} .= $group{$cultivar[$i]}[$j + 3]; # rec.seq.
      $sort_native{$cultivar[$i]} .= $group{$cultivar[$i]}[$j + 1] # cut.pos.
      };
    $sort_native{$cultivar[$i]} .= length $marker{$marker}{$cultivar[$i]}[0]; #considers also total fragment length
    };
  $group{'group_align'}[0][0] = $cultivar[0]; # initialize sorting matrix with first element
  $group{'group_native'}[0][0] = $cultivar[0];
  for (my $i = 1; $i < scalar @cultivar; $i++) # check all cultivars
    {
    my $found_group_align = 0;
    my $found_group_native = 0;
    for (my $j = 0; $j < scalar @{$group{'group_align'}}; $j++) # check in all existing groups
      {
      if ($sort_align{$cultivar[$i]} eq $sort_align{$group{'group_align'}[$j][0]}) # compare with the first element of each group
        {
        push @{$group{'group_align'}[$j]}, $cultivar[$i];
        $found_group_align = 1;
        last
        }
      };
    for (my $j = 0; $j < scalar @{$group{'group_native'}}; $j++)
      {
      if ($sort_native{$cultivar[$i]} eq $sort_native{$group{'group_native'}[$j][0]})
        {
        push @{$group{'group_native'}[$j]}, $cultivar[$i];
        $found_group_native = 1;
        last
        }
      };
    if (not $found_group_align) {$group{'group_align'}[scalar @{$group{'group_align'}}][0] = $cultivar[$i]}; # form a new group
    if (not $found_group_native) {$group{'group_native'}[scalar @{$group{'group_native'}}][0] = $cultivar[$i]}
    };
  return \%group
  };

#>>> Get the number of preceding gaps <<<#

sub PrecedingGaps
  { # invoke with (1) marker name, (2) cultivar name, (3) position to be tested
  my $gaps = 0;
  my @array = @{$marker{$_[0]}{$_[1]}};
  my $position = $_[2];
  shift @array;
  while (@array)
    {
    my $pos = shift @array;
    my $length = shift @array;
    if ($pos <= $position) {$gaps += $length}
    else {last};
    };
  return $gaps
  };

#>>> Generate the complement sequence <<<#

sub DNAComplement
  { # invoke by DNA sequence
  my $seq = shift;
  $seq =~ tr/ACTGactgMRVHmrvhKYBDkybd/TGACtgacKYBDkybdMRVHmrvh/;
  return $seq
  };

#>>> Check if sequence is complementary <<<#

sub Palindromic
  { # invoke by DNA sequence
  my $seq = shift;
  return 1 if ($seq eq reverse DNAComplement($seq));
  return 0
  };

#>>> Does the quality grouping <<<#

sub Test
  {
  my @test = keys %{$_[0]};
  return $test[0] if scalar @test == 1;
  return 2;
  };

sub Results
  {
  my %info = (
    1 => 'Predicted CAPS candidates',
    2 => 'Predicted CAPS candidates (need further inspection: at least one N found within target)',
    3 => 'Probably false positive candidates (need further inspection: at least one N found within target)',
    0 => 'Following markers could not be converted into CAPS'
    );
  print "#Format:\n" ,join "\t", '#Marker', 'Enzyme', 'Total size', 'Restriction Sites', 'Expected Fragments', 'Members', "\n\n";

  for (my $i = 1; $i < 4; $i++) # iterate Test Groups
    {
    next if scalar keys %{$candidate[$i]} == 0;
    print  "$info{$i}\n";
    print '-' x length $info{$i}, "\n\n";
    foreach my $marker (sort keys %{$candidate[$i]})
      {
      foreach my $enzyme (sort keys %{$candidate[$i]{$marker}})
        {
        for (my $j = 0; $j < scalar @{$candidate[$i]{$marker}{$enzyme}{'group_native'}}; $j++)
          {
          my @cutpos; # store cut positions
          for (my $x = 1; $x < scalar @{$candidate[$i]{$marker}{$enzyme}{$candidate[$i]{$marker}{$enzyme}{'group_native'}[$j][0]}}; $x += 4)
            {
            my $cutpos = $candidate[$i]{$marker}{$enzyme}{$candidate[$i]{$marker}{$enzyme}{'group_native'}[$j][0]}[$x];
            push @cutpos, $cutpos if $cutpos > 0;
            };
          @cutpos = sort {$a <=> $b} @cutpos;
          my @fragmentsizes; # stores fragment sizes
          my $length_align = length $marker{$marker}{$candidate[$i]{$marker}{$enzyme}{'group_native'}[$j][0]}[0];
          my $lengthseq = $length_align - PrecedingGaps($marker, $candidate[$i]{$marker}{$enzyme}{'group_native'}[$j][0], $length_align);
          if (@cutpos)
            {
            push @fragmentsizes, $cutpos[0], $lengthseq - $cutpos[$#cutpos];
            for (my $x = 1; $x < scalar @cutpos; $x++) {push @fragmentsizes, $cutpos[$x] - $cutpos[$x - 1]};
            }
          else {$fragmentsizes[0] = $lengthseq};
          @fragmentsizes = sort {$b<=>$a} @fragmentsizes;
          if (@cutpos < 5) # if less than 5 cut
            {
            print join "\t", $marker, $enzyme, $lengthseq;
            print "\t", join ";", @cutpos;
            print "\t", join ";",@fragmentsizes;
            print "\t", join ";",@{$candidate[$i]{$marker}{$enzyme}{'group_native'}[$j]};
            print "\n";
          }
          };
        print "\n";
        }
      print "\n";
      }
    };
  print "$info{0}\n";
  print '-' x length $info{0};
  print "\n\n";
  print join "\n", @{$candidate[0]};
  };
