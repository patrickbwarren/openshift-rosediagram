#!/usr/bin/perl

use strict;
use Getopt::Long;

my $analyse3d = "/sandbox/survex/bin/analyse3d";
my $url_file = "url.txt";
my $downloaded_file = "downloaded.txt";
my $survey_file = "survey.3d";
my $svg_file = "diagram.svg";
my $png_file = "diagram.png";
my $domain = "http://www-patrickbwarren.rhcloud.com/";

my ($command, $stats, $downloaded, $ans);
my ($url, $img, $height, $width, @fields);
my $ns = 16;
my $title = "";
my $maxtime = 5;
my $color = 1;

GetOptions ("ns=i" => \$ns, "title=s" => \$title, 
	    "maxtime=i" => \$maxtime, "color=i" => \$color);

$url = `/usr/bin/awk 'NR == 1 { print \$1 }' $url_file`; chomp($url);

if (-e $downloaded_file) {
    $downloaded = `cat $downloaded_file`; chomp($downloaded);
} else {
    $downloaded = "";
}

if ($url eq $downloaded) {
    print "Using existing download...<br>\n";
} else {
    print "Downloading...<br>\n";
    system("curl --silent --max-time $maxtime '$url' > $survey_file");
    system("echo '$url' > $downloaded_file");
}

if (!-e $survey_file || -z $survey_file) {
    print "<tt>$url</tt><br><br>\n";
    print "Didn't result in a downloaded 3d file for some reason (may have timed out)<br>\n";
    system("rm -f $downloaded_file");
    exit 0;
}

$stats = `$analyse3d --nsector=$ns --title=\"$title\" --plot=$svg_file $survey_file`;

unless ($stats =~ /STATISTICS/) {
    print "<tt>$url</tt><br><br>\n";
    print "Didn't result in a valid 3d file for some reason (not parsed by survex file reader)<br>\n";
    system("rm -f $downloaded_file");
    exit 0;
}

unless ($color) {
    system("sed -i 's/yellow/none/; s/blue/black/; s/<circle/& stroke-dasharray=\"10 10\"/' $svg_file");
}

system("convert $svg_file $png_file");

print "<pre>\n";
print "$url\n";
print "number of sectors = $ns\n";
print "title = $title\n" if ($title ne "");
print $stats, "\n";
print "</pre><br>\n";

$ans = `identify $png_file`; chomp($ans);
@fields = split(/ /, $ans);
($width, $height) = split(/x/, @fields[2]);

printf "<img src=\"%s?%s\" width=\"%i\" height=\"%i\"><br><br>\n", $png_file, time(), $width, $height;

print "<button onclick=\"fetchSVG()\">Download SVG version</button>\n\n";

print "<script>\n";
print "  function fetchSVG() { window.open(\"" . $svg_file . "\"); } \n";
print "</script>\n";

# That's it!
