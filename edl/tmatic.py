import subprocess, logging
logger = logging.getLogger(__name__)

class TMOptions():
    baseCommand=['java', '-classpath', '/slipstream/opt/scripts/jar/trimmomatic.jar']
    #baseCommand=['java', '-classpath', '/common/lib/java/trimmomatic.jar']
    def runTmatic(self):
        command=self.buildCommand()
        logger.info("Running Trimmomatic")
        logger.debug("Tmatic command:\n%s" % (command))
        p=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE, shell=True)
        (stdout,stderr)=p.communicate()
        self.stdout=stdout
        self.stderr=stderr
        self.exitcode=p.returncode

class TMOptionsPE(TMOptions):
    javaclass='org.usadellab.trimmomatic.TrimmomaticPE'
    def __init__(self,forward,reverse,primers=None,pref=None,primerSettings=None):
        self.forward=forward
        self.reverse=reverse
        self.outpref=pref
        self.primers=primers
        self.setOutfiles()
        self.primerSettings=primerSettings

    def setOutfiles(self):
        self.outfiles={}
        if self.outpref is not None:
            base="%s/tmatic" % (self.outpref)
            for suff in ['1u','1p','2u','2p']:
                self.outfiles[suff]="%s.%s" % (base,suff)
        else:
            self.outfiles['1u']="%s.tm.u" %(self.forward)
            self.outfiles['1p']="%s.tm.p" %(self.forward)
            self.outfiles['2u']="%s.tm.u" %(self.reverse)
            self.outfiles['2p']="%s.tm.p" %(self.reverse)

    def buildCommand(self):
        command=" ".join(list(self.baseCommand))
        command+=" " + self.javaclass
        if logging.getLogger().level<=logging.DEBUG:
            command+=" -trimlog %s.log" % self.outfiles['1p']
        command='%s "%s"' % (command,self.forward)
        command='%s "%s"' % (command,self.reverse)
        command='%s "%s"' % (command, self.outfiles['1p'])
        command='%s "%s"' % (command, self.outfiles['1u'])
        command='%s "%s"' % (command, self.outfiles['2p'])
        command='%s "%s"' % (command, self.outfiles['2u'])
        if self.primers is not None:
            if self.primerSettings is None:
                clippingVals="2:40:15"
            else:
                clippingVals=self.primerSettings

            command += " ILLUMINACLIP:%s:%s" % (self.primers,clippingVals)
        return command

class TMOptionsSE(TMOptions):
    javaclass='org.usadellab.trimmomatic.Trimmomatic'
    def __init__(self,input,output,endQuality=5,minLength=45):
        self.input=input
        self.output=output
        self.endQuality=endQuality
        self.minLength=minLength

    def buildCommand(self):
        command=join(list(self.baseCommand))
        command+=" " + self.javaclass
        command = '%s "%s"' %  self.input
        command='%s "%s"' %  self.output
        if self.endQuality>0:
            command+= " " + "LEADING:%d" % (self.endQuality)
            command+= " " + "TRAILING:%d" % (self.endQuality)
        if self.minLength>0:
            command+= " " + "MINLEN:%d" % (self.minLength)
        return command
