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
    parser.add_option("-H", "--htmlFile", help="Location to create output HTML file")
    parser.add_option("-O", "--outputDir", help="directory to put output and support files into")
    parser.add_option("-S","--supportFileLoc", help="directory containing support files")
    parser.add_option("-d", "--debug", default=False, action='store_true')
    #parser.add_option("-r", "--ranks", help="Ranks (comma-separated) to create counts/plots for")
    parser.add_option("--domainFile", help="Copy domain counts to this file")
    parser.add_option("--phylumFile", help="Copy phylum counts to this file")
    parser.add_option("--classFile", help="Copy class counts to this file")
    parser.add_option("--orderFile", help="Copy order counts to this file")
    parser.add_option("--familyFile", help="Copy family counts to this file")
    parser.add_option("--genusFile", help="Copy genus counts to this file")
    parser.add_option("--speciesFile", help="Copy species counts to this file")
    parser.add_option("--idCountsFile", help="Copy raw counts to this file")
    parser.add_option("-T", "--taxids", help="translate hitids to taxids", 
            default=False, action='store_true')

    (options, args) = parser.parse_args()

    if options.debug:
        logging.basicConfig(level=logging.DEBUG)

    # create destination directory
    if not os.path.exists(options.outputDir):
        os.makedirs(options.outputDir)

    # copy support files to direcoty
    for fileName in ["stackedBars.js", "d3.v3.min.js", "style.css"]:
       shutil.copy(os.path.sep.join([options.supportFileLoc,fileName]),
                   os.path.sep.join([options.outputDir,fileName]))

    # run countTaxa (args should be a complete coommand (minus the ranks)
    cmd=list(args)
    cmd.append("-o")
    cmd.append(os.path.sep.join([options.outputDir,"taxa.count.table"]))

    rankMap={}
    for rank in ['domain','phylum','class','order','family','genus','species']:
        try:
            exec "fileName = options.%sFile" % rank
            if fileName is not None:
                rankMap[rank]=fileName
                cmd.append("-r")
                cmd.append(rank)
        except:
            logging.debug("%s not in %s" % (rank, dir(options)))
            pass

    rankList = rankMap.keys()

    # Hack for supporting hitid counts:
    if len(rankList)==0:
        # remove unnecessary maps
        rankList.append("taxids" if options.taxids else "hitids")
        if options.idCountsFile is not None:
            rankMap[rankList[0]]=options.idCountsFile
        i=0
        while i < len(cmd):
            arg=cmd[i]
            if arg=='-n':
                # remove taxdump location
                cmd.pop(i)
                cmd.pop(i)
            elif arg=='-m' and not options.taxids:
                # remove id->taxid map file
                cmd.pop(i)
                cmd.pop(i)
            else:
                i+=1

    #for rank in options.ranks.split(","):
    #    cmd.append("-r")
    #    cmd.append(rank)

    logging.debug("Command: %s" % cmd)
    logging.debug(rankList)

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
                # insert option per rank before closing select tag
                first=True
                #rankList = options.ranks.split(",")
                for rank in rankList:
                    if first:
                        first=False
                        selected='selected="true"'
                    else:
                        selected=''
                    datafilename="taxa.count.table"
                    if len(rankList)>1:
                        if rank.lower()=='domain':
                            datafilename+=".superkingdom"
                        else:
                            datafilename+="." + rank.lower()
                    if rank in rankMap:
                        logging.debug("Copying %s to %s" % (datafilename, 
                                                            rankMap[rank]))
                        shutil.copy2("%s%s%s" % (options.outputDir,
                                                 os.path.sep,
                                                 datafilename),
                                     rankMap[rank])
                    destin.write('<option value="%s" %s>%s</option>\n' % (datafilename, selected, rank.capitalize()))
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
