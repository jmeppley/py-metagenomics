Galaxy datatypes for lastal/lastdb databases
======================================

These Galaxy datatypes are based on the blast_datatypes tool copyright 2010-2013 by Peter Cock, The James Hutton
Institute (formerly SCRI, Scottish Crop Research Institute), UK. All rights reserved.
Contributions/revisions copyright 2012 Edward Kirton. All rights reserved.
Contributions/revisions copyright 2013 Nicola Soranzo. All rights reserved.
Adaptation for lastal/lastdb copyright 2014 John Eppley. All rights reserved.

See the licence text below.

This tool is available from:
http://edminilims.mit.edu:9009/toolshed/view/jmeppley/last_datatypes

The original is available from the Galaxy Tool Shed at:
http://toolshed.g2.bx.psu.edu/view/devteam/blast_datatypes

Installation
============

Doing this automatically via the Tool Shed is probably simplest.


Manual Installation
===================

Normally you would install this via a ToolShed, which would move
the provided lastdb.py file into a suitable location and process the
datatypes_conf.xml entry to be combined with your local configuration.

However, if you really want to this should work for a manual install. Add
the following lines from dataypes_conf.xml in this package to the datatypes_conf.xml file in the Galaxy main folder::

        <datatype extension="lastdbn" type="galaxy.datatypes.lastdb:LastalNucDb" mimetype="text/html" display_in_upload="false"/>
        <datatype extension="lastdbp" type="galaxy.datatypes.lastdb:LastalProtDb" mimetype="text/html" display_in_upload="false"/>

Also create the file lib/galaxy/datatypes/lastdb.py by moving, copying or linking
the lastdb.py file provided in this tar-ball.  Finally add 'import lastdb' near
the start of file lib/galaxy/datatypes/registry.py (after the other import
lines).


Bug Reports
===========

There is no offical bug repository for this package yet...

Licence (MIT)
=============

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

NOTE: This is the licence for the lastal/lastdb Galaxy datatypes **only**. 
