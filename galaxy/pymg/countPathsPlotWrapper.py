#!/usr/bin/python
"""
This script calls countTaxa with the given parameters and writes the results to a new directory with an HTML file that will display the results as a bar chart.
"""

import re, sys, os, shutil, subprocess, logging

def main():
    from optparse import OptionParser

    ## set up CLI
    description = __doc__
    parser = OptionParser(description=description)
    parser.disable_interspersed_args()
    parser.add_option("-M", "--htmlFile", help="Location to create output HTML file")
    parser.add_option("-O", "--outputDir", help="directory to put output and support files into")
    parser.add_option("-S","--supportFileLoc", help="directory containing support files")
    #parser.add_option("-l", "--levels", help="Levels (comma-separated list of pathway levels and/or 'ko') to create counts/plots for")
    parser.add_option("-H", "--keggFile", help="Location of pathway heirachy file. This will be passed to the script only if a pathway level output file is specified")
    parser.add_option("--koFile", help="Copy ko counts to this file")
    parser.add_option("--ko3File", help="Copy level 3 counts to this file")
    parser.add_option("--ko2File", help="Copy level 2 counts to this file")
    parser.add_option("--ko1File", help="Copy level 1 counts to this file")

    (options, args) = parser.parse_args()
    if '-v' in args:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARN)

    # create destination directory
    if not os.path.exists(options.outputDir):
        os.makedirs(options.outputDir)

    # copy support files to direcoty
    for fileName in ["stackedBars.js", "d3.v3.min.js", "style.css"]:
       shutil.copy(os.path.sep.join([options.supportFileLoc,fileName]),
                   os.path.sep.join([options.outputDir,fileName]))

    # run countTaxa (args should be a complete coommand (minus the levels and pathways file)
    cmd=list(args)
    cmd.append("-o")
    cmd.append(os.path.sep.join([options.outputDir,"taxa.count.table"]))


    # check level file options
    levelMap={}
    includesPathways=False
    for level in ['ko','1','2','3']:
        levelString = "ko"+level if level != 'ko' else level
        exec "fileName = options.%sFile" % levelString
        if fileName is not None:
            levelMap[level]=fileName

        # add level and ko options only if something other than orthologs chosen
        if level not in ('ko','gene','ortholog'):
            includesPathways=True
    levels = levelMap.keys()

    if includesPathways:
        cmd.append("-H")
        cmd.append(options.keggFile)
        for level in levels:
            cmd.append("-l")
            cmd.append(level)

    logging.debug("Command: %s" % cmd)

    exitcode=subprocess.call(cmd)
    if exitcode!=0:
        logging.warn("Error (%d) running %r " % (exitcode, cmd))
        sys.exit(exitcode)

    # Generate HTML file
    #  We could use templating here, but I think it's more straightforward to just replace the block of <option>s
    source = open(os.path.sep.join([options.supportFileLoc,"StackedBars.html"]))
    destin = open(options.htmlFile,'w')
    state="looking"
    for line in source:
        if state=="looking":
            m = fileSelectRE.search(line)
            if m:
                state="inSelect"
        elif state=="inSelect":
            m = closeSelectRE.search(line)
            if (m):
                state="copyRemainder"
                # insert option per level before closing select tag
                first=True
                for level in levels:
                    if first:
                        first=False
                        selected='selected="true"'
                    else:
                        selected=''
                    datafilename="taxa.count.table"
                    if len(levels)>1:
                        datafilename+="." + level
                    if level in levelMap:
                        logging.debug("Copying %s to %s" % (datafilename, 
                                                            levelMap[level]))
                        shutil.copy2("%s%s%s" % (options.outputDir,
                                                 os.path.sep,
                                                 datafilename),
                                     levelMap[level])
                    destin.write('<option value="%s" %s>%s</option>\n' % (datafilename, selected, level.capitalize()))
            else:
                # don't copy existing options
                continue
        else:
            # if we're in copyRemainder, don't do anything, just print the line
            pass

        # copy line from source to destination
        destin.write(line)

    source.close()
    destin.close()

fileSelectRE = re.compile(r'<select id="rankSelect"')
closeSelectRE = re.compile(r'</select>')

if __name__ == '__main__':
    main()
