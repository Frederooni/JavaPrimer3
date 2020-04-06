#!/usr/bin/perl -w

my $license = <<END;
#  Copyright (c) 2018
#  Whitehead Institute for Biomedical Research, Steve Rozen
#  (http://purl.com/STEVEROZEN/), and Helen Skaletsky
#  All rights reserved.
#
#        This file is part of primer3 software suite.
#
#        This software suite is is free software;
#        you can redistribute it and/or modify it under the terms
#        of the GNU General Public License as published by the Free
#        Software Foundation; either version 2 of the License, or (at
#        your option) any later version.
#
#        This software is distributed in the hope that it will be useful,
#        but WITHOUT ANY WARRANTY; without even the implied warranty of
#        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#        GNU General Public License for more details.
#
#        You should have received a copy of the GNU General Public License
#        along with this software (file gpl-2.0.txt in the source
#        distribution); if not, write to the Free Software
#        Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#  OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON A THEORY
#  OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
END

use strict;
use warnings;

my $file = "src/main/resources/thal.properties.new";
my $final = "src/main/resources/thal.properties";
open (FILE, ">$file") or die "Error creating $file";
print FILE $license;

my %files;
foreach (<primer3_config/*.d?>) {
  my $fname = $_;
  $fname =~ s|.*/||;
  $fname =~ s/\./_/g;
  $files{$fname} = $_;
}
foreach (sort keys %files) {
  my $key = $_;
  $key =~ s/\./_/g;
  print FILE "$key=";
  my $fname = $files{$key};
  warn "Adding $fname\n";

  open (LOADFILE, "$fname") or die "Error opening $fname";
  my $i = 0;
  while (<LOADFILE>) {
    my $line = $_;
    $line =~ s/\t/\\t/g;
    $line =~ s/\n/\\n/g;
    if ($line =~ /[^A-Za-z0-9\.\-\\]/) {
      die "  Unknown char: $line\n";
    }
    print FILE "$line";
  }
  print FILE "\n";
  close(LOADFILE);
}

close(FILE);
rename($file, $final);
