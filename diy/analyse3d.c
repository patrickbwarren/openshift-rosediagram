/* analyse3d.c */
/* Generate passage orientation statistics from a .3d file */
/* Based on dump3d.c written by Olly Betts
 * Modified to produce passage orientation statistics by Patrick Warren
 *
 * Copyright (C) 2001,2002,2006 Olly Betts, 2010 Patrick Warren
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */

/*
gcc -DHAVE_CONFIG_H -I.. -O2 -o rose3d rose3d.c \
img.o useful.o cmdline.o message.o filename.o osdepend.o z_getopt.o getopt1.o  -lm 

gcc -Wall -W -Wunused -Wshadow -Wpointer-arith -Wmissing-prototypes \
-Wwrite-strings -Wredundant-decls -Wnested-externs -Wcast-align -I.. -g -O2 \
-o analyse3d analyse3d.c date.o img_hosted.o useful.o cmdline.o message.o \
filename.o osdepend.o z_getopt.o getopt1.o  -lm
*/

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cmdline.h"
#include "debug.h"
#include "filelist.h"
#include "img_hosted.h"

void lprint(const char *, int, double, const char *);
double fancyround(double);
double fancyfloor(double);

static const struct option long_opts[] = {
  /* const char *name; int has_arg (0 no_argument, 1 required_*, 2 optional_*); */
  /* int *flag; int val; */
  {"survey", required_argument, 0, 's'},
  {"plot", required_argument, 0, 'p'},
  {"nsectors", required_argument, 0, 'n'},
  {"title", required_argument, 0, 't'},
  {"help", no_argument, 0, HLP_HELP},
  {"version", no_argument, 0, HLP_VERSION},
  {0, 0, 0, 0}
};

#define short_opts "s:p:n:t:"

static struct help_msg help[] = {
/*				<-- */
   {HLP_ENCODELONG(0),	      /*only load the sub-survey with this prefix*/199, 0},
   {0, 0, 0}
};

double fancyround(double x) {
  double b, c;
  b = pow(10.0, floor(log10(x)));
  c = floor(x / b + 0.5);
  /* printf("x, b, c, x/b = %g %g %g %g\n", x, b, c, x/b); */
  if ((c == 1.0 && x/b > 1.25) || (c == 2.0 && x/b < 1.75)) c = 1.5;
  if ((c == 2.0 && x/b > 2.25) || (c == 3.0 && x/b < 2.75)) c = 2.5;
  return c * b;
}

double fancyfloor(double x) {
  double b, c;
  b = pow(10.0, floor(log10(x)));
  c = floor(x / b);
  if ((c == 1.0 && x/b > 1.25) || (c == 2.0 && x/b < 1.75)) c = 1.5;
  if ((c == 2.0 && x/b > 2.25) || (c == 3.0 && x/b < 2.75)) c = 2.5;
  return c * b;
}

void lprint(const char *s, int n, double t, const char *f) {
  printf("%10s  %6i  %10.2f  %8.2f%s\n", s, n, t, n>0 ? t/(double)n : 0.0, f);
}

int main(int argc, char **argv) {
  FILE *fp = NULL;
  char *fnm, *svgfile = NULL, *title = NULL;
  img *pimg;
  img_point pt, oldpt;
  int ns = 16;
  int i, code, ival;
  int nlegs, nsurf, ndup, nsplay, nhoriz;
  int *nsector = NULL;
  double *tsector = NULL;
  double dx, dy, dz, ang, len, lhoriz, pi, rns;
  double tlegs, tsurf, tdup, tsplay, thoriz;
  double alo, ahi, tmax, rad, tmean;
  double xx0, yy0, xx1, yy1;
  const char *survey = NULL;
  bool fRewind = fFalse;

  msg_init(argv);

  cmdline_init(argc, argv, short_opts, long_opts, NULL, help, 1, 1);

  while (1) {
    int opt = cmdline_getopt();
    if (opt == EOF) break;
    if (opt == 's') survey = strdup(optarg);
    if (opt == 'p') svgfile = strdup(optarg);
    if (opt == 't') title = strdup(optarg);
    if (opt == 'n') {
      if (sscanf(optarg, "%i", &ival) == 1 && ival > 0) ns = ival;
    }
  }

  fnm = argv[optind];

  if (title == NULL) title = fnm;
   
  pi = atan2(1.0, 1.0) * 4.0;
  rns = 1.0 / (double)ns;

  nsector = (int *)malloc(ns*sizeof(int));
  tsector = (double *)malloc(ns*sizeof(double));
  if (nsector == NULL || tsector == NULL) {
    fprintf(stderr, "Ran out of space for malloc in rose3d.c\n"); 
    exit(1); /* Should handle this error better */
  }
  for (i=0; i<ns; i++) { nsector[i] = 0; tsector[i] = 0.0; }
  nlegs = nsurf = ndup = nsplay = nhoriz = 0;
  tlegs = tsurf = tdup = tsplay = thoriz = 0.0;

  pimg = img_open_survey(fnm, survey);
  if (!pimg) fatalerror(img_error(), fnm);

  code = img_BAD;
  do {
    if (code == img_STOP) {
      fprintf(stderr, "<<< REWIND <<<\n");
      fRewind = fFalse;
      if (!img_rewind(pimg)) fatalerror(img_error(), fnm);
    }
    oldpt.x = oldpt.y = oldpt.z = 0.0;
    do {
      code = img_read_item(pimg, &pt);
      switch (code) {
      case img_MOVE:
	oldpt.x = pt.x; oldpt.y = pt.y; oldpt.z = pt.z;
	break;
      case img_LINE:
	dx = pt.x - oldpt.x; dy = pt.y - oldpt.y; dz = pt.z - oldpt.z;
	oldpt.x = pt.x; oldpt.y = pt.y; oldpt.z = pt.z;
	lhoriz = sqrt(dx*dx + dy*dy);
	len = sqrt(dx*dx + dy*dy + dz*dz);
	if (pimg->flags & img_FLAG_SURFACE) { tsurf += len; nsurf++; }
	else if (pimg->flags & img_FLAG_DUPLICATE) { tdup += len; ndup++; }
	else if (pimg->flags & img_FLAG_SPLAY) { tsplay += len; nsplay++; }
	else {
	  tlegs += len; nlegs++;
	  if (lhoriz > 0.0) {
	    thoriz += lhoriz; nhoriz++;

	    /* Calculate the bearing in ang, given that x is
	       east and y is north.  Add 1/2 len to the sectors
	       containing ang and ang + 180. */

	    ang = atan2(dx, dy) * 180.0 / pi;

	    if (ang < 0.0) ang += 360.0;
	    i = (int)(ns * ang / 360.0 + 0.5);
	    if (i == ns) i = 0;
	    nsector[i]++; tsector[i] += lhoriz;

	    ang += 180.0;

	    if (ang >= 360.0) ang -= 360.0;
	    i = (int)(ns * ang / 360.0 + 0.5);
	    if (i == ns) i = 0;
	    nsector[i]++; tsector[i] += lhoriz;

	  }
	}
	break;
      case img_LABEL:
	break;
      case img_XSECT:
	break;
      case img_XSECT_END:
	break;
      case img_BAD:
	img_close(pimg);
	fatalerror(img_error(), fnm);
      }
    } while (code != img_STOP);
  } while (fRewind);

  img_close(pimg);

  printf("\nSUMMARY STATISTICS (lengths in metres)\n\n");
  printf("%10s  %6s  %10s  %8s\n", "legs", "number", "total", "mean");
  lprint("normal", nlegs, tlegs, "");
  lprint("[ horiz", nhoriz, thoriz, " ]");
  lprint("surface", nsurf, tsurf, "");
  lprint("duplicate", ndup, tdup, "");
  lprint("splay", nsplay, tsplay, "");
  lprint("TOTAL", nsurf+ndup+nsplay+nlegs, tsurf+tdup+tsplay+tlegs, "");

  if (svgfile != NULL && (fp = fopen(svgfile, "w")) != NULL) {

    fprintf(fp, "<?xml version=\"1.0\" standalone=\"no\"?>\n");
    fprintf(fp, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n");
    fprintf(fp, "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
    fprintf(fp, "<svg width=\"12cm\" height=\"12cm\" viewBox=\"0 0 1200 1200\"\n");
    fprintf(fp, "     xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n");
    fprintf(fp, "  <g transform=\"translate(600,610)\">\n");

    for (i=0, tmax=0.0; i<ns; i++) if (tmax < tsector[i]) tmax = tsector[i];
    tmean = 2.0 * thoriz / (double)ns;

    /* for (i=0; i<ns; i++) printf("sector %3i total %g\n", i, tsector[i]); */

    /* printf("tmean = %f  tmax = %f\n", tmean, tmax); */

    for (i=0; i<ns; i++) {
      alo = 2.0 * pi * (((double)i - 0.5)*rns); 
      ahi = 2.0 * pi * (((double)i + 0.5)*rns);
      rad = 420.0 * tsector[i] / tmax;
      xx0 = rad * sin(alo); yy0 = - rad * cos(alo);
      xx1 = rad * sin(ahi); yy1 = - rad * cos(ahi);
      if (i == 0) fprintf(fp, "    <path d=\"M %7.2f, %7.2f ", xx0, yy0);
      else fprintf(fp, "             L %7.2f, %7.2f ", xx0, yy0);
      fprintf(fp, "A %6.2f, %6.2f 0 0,1 %7.2f, %7.2f", rad, rad, xx1, yy1);
      if (i < ns-1) fprintf(fp, "\n");
      else fprintf(fp, " z\"\n");      
    }
    fprintf(fp, "          fill=\"yellow\" stroke=\"black\" stroke-width=\"2\" />\n");

    rad = fancyround(tmean);
    if (rad > tmax) rad = fancyfloor(tmean);
    /* printf("rad = %g\n", rad); */

    fprintf(fp, "    <circle r=\"%0.2f\" fill=\"none\" stroke=\"blue\" stroke-width=\"2\" />\n", 420.0 * rad / tmax);

    fprintf(fp, "    <path d=\"M-450,0 h900\" fill=\"none\" stroke=\"black\" stroke-width=\"2\" />\n");
    fprintf(fp, "    <path d=\"M0,-450 v900\" fill=\"none\" stroke=\"black\" stroke-width=\"2\" />\n");
    fprintf(fp, "    <path d=\"M-318.2,-318.2 L318.2,318.2\" fill=\"none\" stroke=\"black\" stroke-width=\"2\" />\n");
    fprintf(fp, "    <path d=\"M318.2,-318.2 L-318.2,318.2\" fill=\"none\" stroke=\"black\" stroke-width=\"2\" />\n");

    fprintf(fp, "    <text x=\"-480\" y=\"0\" font-family=\"Verdana\" font-size=\"50\" \n");
    fprintf(fp, "          fill=\"black\" dy=\"0.35em\" text-anchor=\"middle\">W</text>\n");
    fprintf(fp, "    <text x=\"480\" y=\"0\" font-family=\"Verdana\" font-size=\"50\" \n");
    fprintf(fp, "          fill=\"black\" dy=\"0.35em\" text-anchor=\"middle\">E</text>\n");
    fprintf(fp, "    <text x=\"0\" y=\"-480\" dy=\"0.35em\" font-family=\"Verdana\" font-size=\"50\"\n");
    fprintf(fp, "          fill=\"black\" text-anchor=\"middle\">N</text>\n");
    fprintf(fp, "    <text x=\"0\" y=\"480\" dy=\"0.35em\" font-family=\"Verdana\" font-size=\"50\"\n"); 
    fprintf(fp, "          fill=\"black\" text-anchor=\"middle\">S</text>\n");

    fprintf(fp, "    <text x=\"344.4\" y=\"344.4\" dy=\"0.35em\" font-family=\"Verdana\" font-size=\"50\"\n"); 
    fprintf(fp, "          fill=\"black\" text-anchor=\"middle\">SE</text>\n");

    fprintf(fp, "    <text x=\"-344.4\" y=\"344.4\" dy=\"0.35em\" font-family=\"Verdana\" font-size=\"50\"\n"); 
    fprintf(fp, "          fill=\"black\" text-anchor=\"middle\">SW</text>\n");

    fprintf(fp, "    <text x=\"344.4\" y=\"-344.4\" dy=\"0.3em\" font-family=\"Verdana\" font-size=\"50\"\n"); 
    fprintf(fp, "          fill=\"black\" text-anchor=\"middle\">NE</text>\n");

    fprintf(fp, "    <text x=\"-344.4\" y=\"-344.4\" dy=\"0.3em\" font-family=\"Verdana\" font-size=\"50\"\n"); 
    fprintf(fp, "          fill=\"black\" text-anchor=\"middle\">NW</text>\n");

    /* fprintf(fp, "    <path d=\"M-550,0 h1100\" fill=\"none\" stroke=\"red\" stroke-width=\"1\" />\n"); */
    /* fprintf(fp, "    <path d=\"M0,-550 v1100\" fill=\"none\" stroke=\"red\" stroke-width=\"1\" />\n"); */
    /* fprintf(fp, "    <path d=\"M-550,-550 L550,550\" fill=\"none\" stroke=\"red\" stroke-width=\"1\" />\n"); */
    /* fprintf(fp, "    <path d=\"M-550,550 L550,-550\" fill=\"none\" stroke=\"red\" stroke-width=\"1\" />\n"); */
    /* fprintf(fp, "    <circle r=\"420\" fill=\"none\" stroke=\"red\" stroke-width=\"1\" />\n"); */
    /* fprintf(fp, "    <circle r=\"450\" fill=\"none\" stroke=\"red\" stroke-width=\"1\" />\n"); */
    /* fprintf(fp, "    <circle r=\"480\" fill=\"none\" stroke=\"red\" stroke-width=\"1\" />\n"); */

    fprintf(fp, "  </g>\n");

    fprintf(fp, "    <text x=\"600\" y=\"70\" font-family=\"Verdana\" font-size=\"60\"\n"); 
    fprintf(fp, "          fill=\"black\" text-anchor=\"middle\">%s</text>\n", title);

    fprintf(fp, "    <text x=\"600\" y=\"1190\" font-family=\"Verdana\" font-size=\"50\"\n"); 
    fprintf(fp, "          fill=\"black\" text-anchor=\"middle\">length %g %s, circle radius is %g %s</text>\n",
	    thoriz > 10000 ? round(thoriz/100)/10 : round(thoriz), thoriz > 10000 ? "km" : "m",  
	    rad > 1000 ? rad/1000 : rad, rad > 1000 ? "km" : "m");
    fprintf(fp, "</svg>\n");

    fclose(fp);
  }

  free(nsector); free(tsector);

  return 0;
}

/**************************************
# download plotrix_*.tar.gz
# sudo R CMD INSTALL plotrix_*.tar.gz

require(plotrix)

fancy.round <- function (x, round.down=FALSE) {
  if (round.down) { offset <- 0.0; }
  else { offset <- 0.5; }
  b <- 10^as.integer(log10(x));
  c <- as.integer(x/b + offset);
  if ((c == 1 && x/b > 1.25) || (c == 2 && x/b < 1.75)) c = 1.5;
  if ((c == 2 && x/b > 2.25) || (c == 3 && x/b < 2.75)) c = 2.5;
  c * b;
}

as.scale <- function(x, max=NA) {
  d <- fancy.round(x);
  if (!is.na(max) && d > max) d <- fancy.round(x, round.down=TRUE);
  d;
}

rose.read <- function(file) {
  read.delim(file, header=FALSE, 
    col.names=c("ang", "len", "nlegs", "flag"))
}

# Change rose3d_executable in this, if rose3d is not in the current directory.

rose3d.read <- function(file, nsectors=16, resolution=1) {
  rose3d_executable <- "./rose3d";
  s <- paste(rose3d_executable, "-q -n", nsectors, "-r", resolution, file);
  rose.read(pipe(s));
}

rose.plot<-function(data, rp.type="p", rad.plot=NA, rad.plot.sf=1.1, ...) {
  if (is.na(rad.plot)) rad.plot <- max(data$len) * rad.plot.sf;
  polar.plot(data$len, data$ang, 
    start=90, clockwise=TRUE, rp.type=rp.type,
    labels=c("N","NE","E","SE","S","SW","W","NW"), label.pos=seq(0,315,45),
    label.prop=1.1, show.grid.labels=FALSE, radial.lim=c(0.0, rad.plot),
    grid.col=1, show.grid=FALSE, mar=c(3,2,3,2), ...);
  return(rad.plot);
}

add.circle <- function(rad.circle=1.0, rad.plot, ...) {
  polar.plot(rep(rad.plot,360), seq(0,359), rp.type="p", add=TRUE,
    radial.lim=c(0.0, rad.plot), ...)
}

rose3d.plot <- function(file, nsectors=16, resolution=1,  poly.col=7, line.col=1, 
                        rad.circle=NA, rad.plot=NA, rad.plot.sf=1.1, ...) {
  data <- rose3d.read(file, nsectors, resolution);
  rad.plot <- rose.plot(data, poly.col=poly.col, line.col=line.col, rad.plot=rad.plot, ...);
  mean <- sum(data$len[data$flag == 1]) / nsectors;
  if (is.na(rad.circle)) rad.circle <- as.scale(mean, max=max(data$len));
  if (rad.circle > 0) add.circle(rad.circle, rad.plot, line.col=4);
  if (rad.circle < 1000) { scale <- rad.circle; unit <- "m"; }
  else { scale <- rad.circle/1000; unit <- "km"; }
  return(c(rad.plot, rad.circle, paste("circle radius is", scale, unit)));
}

matienzo.process.3d <- function(site, title, len=NA, ns=NA, title.show=TRUE, 
                                png.write=FALSE, png.name="_rose.png") {
  file <- paste(sep="", site, ".3d");
  if (file.exists(file)) {
    if (is.na(ns)) {
      ns = 16;
      if (!is.na(len)) {
        if (as.integer(len) > 5000) ns = ns*2;
        if (as.integer(len) > 100000) ns = ns * 2;
      }
    }
    if (png.write) {
      img <- paste(sep="", site, png.name);
      png(img);
    }
    output <- rose3d.plot(file, nsectors=as.integer(ns));
    if (title.show) {
      dist <- as.double(output[[1]]);
      if (is.na(len)) {
        subtitle <- sprintf("%s, %s", file, output[[3]]);
      } else {
        subtitle <- sprintf("%s, %i m, %s", file, as.integer(len), output[[3]]);
      }
      text(0, dist*1.25, title, cex=1.5);
      text(0, -dist*1.25, subtitle, cex=1.2);
    }
    if (png.write) {
      dev.off();
      cat(img); 
      cat(", nsectors = "); cat(as.integer(ns));
      if (title.show) {
        cat(", title = "); cat(title); 
        cat(", subtitle = "); cat(subtitle);
      }
      cat("\n");
    }
  } else {
    cat(file); cat(" does not exist\n");
  }
}

matienzo.read.all <- function(file) {
  read.table(file, header=F, col.names=c("site", "name", "len", "ns"), sep="\t",
    colClasses = c('character', 'character', 'integer', 'integer'));
}

matienzo.process.all <- function(matienzo_caves, title.show=TRUE, png.name="_rose.png") {
  apply(matienzo_caves, 1, 
    function(x) matienzo.process.3d(x[[1]], x[[2]], x[[3]], x[[4]],
      png.write=TRUE, title.show=title.show, png.name=png.name));
}

matienzo.get.3d <- function(site, max.time=5, curl.opts="", 
    matienzo.url= "http://www.geography.lancs.ac.uk/matienzo/surveys/") {
  curl_executable <- "curl";
  file <- paste(sep="", site, ".3d");
  url <- paste(sep="", matienzo.url, file);
  temp_3d = "rose3d.R.3d";
  rc <- system(paste(curl_executable, "--max-time", max.time, curl.opts, url, ">", temp_3d));
  if (rc == 0) {
    rc <- system(paste("mv -f", temp_3d, file));
    cat(paste("Downloaded", url, "successfully\n"));
  } else {
    rc <- system(paste("rm -f", temp_3d));
    cat(paste("Attempt to download", url, "failed\n"));
  }
}

matienzo.get.all <- function(matienzo_caves, max.time=5, curl.opts="", 
    matienzo.url= "http://www.geography.lancs.ac.uk/matienzo/surveys/") {
  apply(matienzo_caves, 1, function(x) matienzo.get.3d(x[[1]], 
    max.time=max.time, curl.opts=curl.opts, matienzo.url=matienzo.url));
}

# Examples
# matienzo.read.all("matienzo_caves.dat") -> caves
# matienzo.get.all(caves)
# or matienzo.get.all(caves, "--proxy proxy.server:port")
# matienzo.process.all(caves)
# matienzo.process.all(caves, title.show=FALSE, png.name="_rose_notitle.png")
*******************************************/
