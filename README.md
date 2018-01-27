# Passage orientation diagrams

_Code being migrated to an Inkscape plugin as OpenShift 2.0 no longer exists_

This is a simple web interface to generate passage orientation 'rose' diagrams from survex .3d files that are available online.  It uses Redhat OpenShift and is cloned from https://github.com/openshift-quickstart/r-quickstart.
I can recommend reading through the docs associated with the R quickstart to see just how closely I've copied the underlying ideas!

## How it works

The heart is _analyse3d.c_ which reads survex .3d files and spits out statistics and the rose diagram in SVG format.  This is wrapped in a perl script called _wrap.pl_ which handles downloading the .3d file from a URL, running _analyse3d_, capturing the output, converting the SVG into PNG, and writing out some valid HTML.  The perl script is in turn is called by _rubyserver.rb_ which is a lightweight web server based on WEBrick.  The web server is expecting to receive POST data generated by _index.html_ on the client side.  The POST data specifies the URL for the .3d file and some other parameters.  Finally, _index.html_ on the client side contains a small amount of javascript (jquery) and the necessary gubbins to generate the POST request, and then use the HTML returned by _webserver.rb_ to update the web page in an AJAX-like manner.

The deployment contains one neat trick I found on the internet.  In the case where the user makes multiple calls it is necessary to avoid the web browser using a cached image.  This is done by appending the time (in unix seconds) to the image URL, namely (for instance)
```
<img src="diagram.png?1395845892">
```
The number after the question mark is generated by a call to the unix `time()` function from with the perl script, and is different for every call.

## Installation

First create a new application using the DIY cartridge with the R quickstart repository as a template

```
rhc app create r diy --from-code=https://github.com/openshift-quickstart/r-quickstart.git
```
Next replace all the files in the _diy_ and _.openshift/action-hooks_ directories.

### Survex

Once created, we need to set up survex in the gear. This is a bit of a pain because the normal build process looks for and fails to find the graphics libraries.  Here is a complete kludge:  SSH into the OpenShift instance.  Download the survex source tarball using curl and unpack it somewhere safe, for example in `~/app-root/data`.  Optionally run autoreconf, and then try to run _./configure_.  This will fail, because it can't find the graphics libraries.  Therefore, edit _./configure_ and delete the clause of the if/then/else statement that is causing the problems.  Once _./configure_ runs all the way through, drop into the _src_ directory and check that you can build _dump3d_ from the Makefile.  Now copy _analyse3d.c_ into the same directory and compile using the output of `make dump3d` as a template.  Specifically, replace the part that reads `-o dump3d dump3d.o` with `-o analyse3d analyse3d.c -I..` (the extra include directory is needed so that the compiler can find the header files). 

Got all that?  I told you it was a complete kludge!  But since it is server side setup, no-one need ever know ;-)

Now check you can run _analyse3d_ from within the src directory:
```
./analyse3d
```
This should produce some output complaining about too few arguments.  You check further by fetching your favourite survex .3d file from the internet, for example
```
curl http://www.geography.lancs.ac.uk/matienzo/surveys/2889.3d > 2889.3d
./analyse3d 2889.3d --title="Torca La Vaca" --nsectors=32 --plot=2889.svg
```

### Making it all work together

If you move the executable _analyse3d_ it may fail because it is looking for additional files located elsewhere on the survex file tree.  This dependency would normally be handled by `make install`.  In the present case you can either hack the installation step too, or perhaps more simply just create a link pointing to the survex src directory.  For example, from within the survex src directory

```
mkdir -p /sandbox/survex
ln -s `pwd` /sandbox/survex/bin
```
check that you can run _analyse3d_ with `/sandbox/survex/bin/analysed3d`.  Then make sure that the `$analyse3d` variable in _wrap.pl_ is pointing to the right place.  

That's it --- it should all run now.

## Usage

Point a web browser at your application URL and away you go.
