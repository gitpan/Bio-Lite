package Bio::Lite;
{
  $Bio::Lite::DIST = 'Bio-Lite';
}
# ABSTRACT: Bio::Lite is a perl module that aims to answer the same questions as Bio-perl, FASTER and using a SIMPLIFIED API.
$Bio::Lite::VERSION = '0.001';
use Carp;

use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(seqFileIterator reverseComplemente pairedEndSeqFileIterator gffFileIterator);
our @EXPORT_OK = qw();


sub reverseComplemente {
  my $seq = reverse shift;
  $seq =~ tr/ACGTacgt/TGCAtgca/;
  return $seq;
}



sub seqFileIterator {
  my ($file,$format) = @_;

  croak "Missing file in argument of seqFileIterator" if !defined $file;

  # Get file handle for $file
  my $fh = getReadingFileHandle($file);

  # Automatic file extension detection
  if(!defined $format) {
    if($file =~ /\.(fasta|fa)(\.|$)/) {
      $format = 'fasta';
    } elsif($file =~ /\.(fastq|fq)(\.|$)/) {
      $format = 'fastq';
    } else {
      croak "Undefined file extension";
    }
  } else {
    $format = lc $format;
  }

  # FASTA ITERATOR
  if ($format eq 'fasta') {
    # Read prev line for FASTA because we dont know the number
    # of line used for the sequence
    my $prev_line = <$fh>;
    chomp $prev_line;
    return sub {
      my ($name,$seq,$qual); 
      if(defined $prev_line) {
        ($name) = $prev_line =~ />(.*)$/;
        $prev_line = <$fh>;
        # Until we find a new sequence identifier ">", we
        # concatenate the lines corresponding to the sequence
        while(defined $prev_line && $prev_line !~ /^>/) {
          chomp $prev_line;
          $seq .= $prev_line;
          $prev_line = <$fh>;
        }
        return {name => $name, seq => $seq, qual => $qual};
      } else {
        return undef;
      }
    };
  # FASTQ ITERATOR
  } elsif ($format eq 'fastq') {
    return sub {
      my ($name,$seq,$qual); 
      ($name) = <$fh> =~ /@(.*)$/;
      if(defined $name) {
        $seq = <$fh>;
        chomp $seq;
        <$fh>; # skip second seq name (useless line)
        $qual = <$fh>;
        chomp $qual;
        return {name => $name, seq => $seq, qual => $qual};
      } else {
        return undef;
      }
    };
  } else {
    croak "Undefined file format";
  }
}


sub pairedEndSeqFileIterator {
  my($file1,$file2,$format) = @_;

  my $it1 = seqFileIterator($file1,$format);
  my $it2 = seqFileIterator($file2,$format);

  return sub {
    my $entry1 = $it1->();
    my $entry2 = $it2->();
    if(defined $entry1 && defined $entry2) {
      return { read1 => $entry1, read2 => $entry2 };
    } else {
      return undef;
    }
  };
}


sub gffFileIterator {
  my ($file,$format) = @_;

  croak "Missing arguments in gffFileIterator" if !defined $file || !defined $format;

  $format = lc $format;
  my $fh = getReadingFileHandle($file);

  my $attribute_split;

  if ($format eq 'gff3') {
    $attribute_split = sub {my $attr = shift; return $attr =~ /(\S+)=(.*)/;};
  } elsif ($format eq 'gtf' || $format eq 'gff2') {
    $attribute_split = sub {my $attr = shift; return $attr  =~ /(\S+)\s+"(.*)"/;};
  } else {
    croak "Undefined gff format";
  }

  return sub {
    my $line = <$fh>;
    if (defined $line) {
      my($chr,$source,$feature,$start,$end,$score,$strand,$frame,$attributes) = split("\t",$line);
      my @attributes_tab = split(";",$attributes);
      my %attributes_hash;
      foreach my $attr (@attributes_tab) {
        my ($k,$v) = $attribute_split->($attr);
        $attributes_hash{$k} = $v;
      }
      return { chr        => $chr,
        source     => $source,
        feature    => $feature, 
        start      => $start, 
        end        => $end, 
        score      => $score, 
        strand     => $strand,
        frame      => $frame,
        attributes => \%attributes_hash,
      };
    } else {
      return undef;
    }
  };
}


sub getReadingFileHandle {
  my $file = shift;
  my $fh;
  if($file =~ /\.gz$/) {
    open($fh,"gunzip -c $file |") or die ("Cannot open $file");
  } else {
    open($fh,"< $file") or die ("Cannot open $file");
  }
  return $fh;
}


1;

__END__

=pod

=encoding UTF-8

=head1 NAME

Bio::Lite - Bio::Lite is a perl module that aims to answer the same questions as Bio-perl, FASTER and using a SIMPLIFIED API.

=head1 VERSION

version 0.001

=head1 SYNOPSIS

Keep it simple, keep it fast.

Bio::Lite is a set of subroutines that aim to answer similar question as
Bio-perl distribution without the complexity. 

Bio::Lite is fast, simple and does not make use of
complexe data struture and object that lead to slow execution.

All methods can be imported with a single "use Bio::Lite"

Bio::Lite is also a lightweigth, single, module with no dependencies.

=head1 UTILS

=head2 reverseComplemente

Reverse complemente the (nucleotid) sequence in arguement.

=head1 PARSING

This are some tools that aim to read (bio) files like
- Sequence files : FASTA, FASTQ
- Annotation files : GFF3, GTF2, BED6, BED12, ...
- Alignement files : SAM, BAM
- 

=head2 seqFileIterator

Open Fasta, or Fastq files (can be gziped).
seqFileIterator has an automatic file extension detection but you can force it
using a second parameter with the format : 'fasta' or 'fastq'

  my $it = seqFileIterator('file.fastq','fastq');
  while(my $entry = $it->()) {
    print "Sequence name   : $entry->{name}
           Sequence        : $entry->{seq}
           Sequence quality: $entry->{qual}","\n";
  }

=head2 pairedEndSeqFileIterator

Open Paired-End Sequence file unsing seqFileIterator()

  my $it = pairedEndSeqFileIterator($file);
  while (my $entry = $it->()) {
    print "Read_1 : $entry->{read1}->{seq}
           Read_2 : $entry->{read2}->{seq}";
  }

=head2 gffFileIterator 

manage GFF3 and GTF2 file format

Return a hashref with the annotation parsed :

  { chr => ...,
    source => ...,
    feature => ...,
    start => ...,
    end => ...,
    strand ...,
    frame ...,
    attributes => { id => val, ...}
  }

=head1 FILES IO

=head2 getReadingFileHandle

Return a file handle for the file in argument.
Display errors if file cannot be oppenned and manage gzipped files (based on .gz file extension)

=head1 TODO

- add a seqFileIterator that manage paired-end files

=head1 AUTHOR

Jérôme Audoux <jerome.audoux@gmail.com>

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2014 by Jérôme Audoux.

This is free software, licensed under:

  The GNU General Public License, Version 3, June 2007

=cut
