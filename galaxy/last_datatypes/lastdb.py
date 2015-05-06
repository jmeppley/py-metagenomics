"""
Lastdb classes
"""

from galaxy.datatypes.data import Data

class _LastalDb(object):
    """Base class for LAST database datatype."""

    def set_peek( self, dataset, is_multi_byte=False ):
        """Set the peek and blurb text."""
        if not dataset.dataset.purged:
            dataset.peek  = "LASTal database (multiple files)"
            dataset.blurb = "LASTal database (multiple files)"
        else:
            dataset.peek = 'file does not exist'
            dataset.blurb = 'file purged from disk'

    def display_peek( self, dataset ):
        """Create HTML content, used for displaying peek."""
        try:
            return dataset.peek
        except:
            return "LASTal database (multiple files)"

    def display_data(self, trans, data, preview=False, filename=None,
                     to_ext=None, size=None, offset=None, **kwd):
        """Documented as an old display method, but still gets called via tests etc

        This allows us to format the data shown in the central pane via the "eye" icon.
        """
        if filename is not None and filename != "index":
            #Change nothing - important for the unit tests to access child files:
            return Data.display_data(self, trans, data, preview, filename,
                                     to_ext, size, offset, **kwd)
        if self.file_ext == "lastdbn":
            title = "This is a nucleotide LASTal database"
        elif self.file_ext =="lastdbp":
            title = "This is a protein LASTal database"
        else:
            #Error?                                                                                                                                                                     
            title = "This is a LASTal database."
        msg = ""
        try:
            #Try to use any text recorded in the dummy index file:
            handle = open(data.file_name, "rU")
            msg = handle.read().strip()
            handle.close()
        except Exception, err:
            #msg = str(err)
            pass
        if not msg:
            msg = title
        #Galaxy assumes HTML for the display of composite datatypes,
        return "<html><head><title>%s</title></head><body><pre>%s</pre></body></html>" % (title, msg)


class LastalNucDb( _LastalDb, Data ):
    """Class for nucleotide LASTal database files."""
    file_ext = 'lastdbn'
    allow_datatype_change = False
    composite_type = 'basic'

    def __init__(self, **kwd):
        Data.__init__(self, **kwd)
        self.add_composite_file('lastdb.ssp', is_binary=True)
        self.add_composite_file('lastdb.tis', is_binary=True)
        self.add_composite_file('lastdb.sds', is_binary=True)
        self.add_composite_file('lastdb.suf', is_binary=True)
        self.add_composite_file('lastdb.des', is_binary=True)
        self.add_composite_file('lastdb.bck', is_binary=True)
        self.add_composite_file('lastdb.prj', is_binary=False)
		## These next two are my own concotions
		# a map from id or accession to taxid
        self.add_composite_file('lastdb.tax', is_binary=False, optional=True)
		# a map form id or accession to description
        self.add_composite_file('lastdb.ids', is_binary=False, optional=True)
		## These are for multi-volumen DBs. I'm not sure they are necessary...
        self.add_composite_file('lastdb0.ssp', is_binary=True, optional=True)
        self.add_composite_file('lastdb0.tis', is_binary=True, optional=True)
        self.add_composite_file('lastdb0.sds', is_binary=True, optional=True)
        self.add_composite_file('lastdb0.suf', is_binary=True, optional=True)
        self.add_composite_file('lastdb0.des', is_binary=True, optional=True)
        self.add_composite_file('lastdb0.bck', is_binary=True, optional=True)
        self.add_composite_file('lastdb0.prj', is_binary=False, optional=True)
        self.add_composite_file('lastdb1.ssp', is_binary=True, optional=True)
        self.add_composite_file('lastdb1.tis', is_binary=True, optional=True)
        self.add_composite_file('lastdb1.sds', is_binary=True, optional=True)
        self.add_composite_file('lastdb1.suf', is_binary=True, optional=True)
        self.add_composite_file('lastdb1.des', is_binary=True, optional=True)
        self.add_composite_file('lastdb1.bck', is_binary=True, optional=True)
        self.add_composite_file('lastdb1.prj', is_binary=False, optional=True)
        self.add_composite_file('lastdb2.ssp', is_binary=True, optional=True)
        self.add_composite_file('lastdb2.tis', is_binary=True, optional=True)
        self.add_composite_file('lastdb2.sds', is_binary=True, optional=True)
        self.add_composite_file('lastdb2.suf', is_binary=True, optional=True)
        self.add_composite_file('lastdb2.des', is_binary=True, optional=True)
        self.add_composite_file('lastdb2.bck', is_binary=True, optional=True)
        self.add_composite_file('lastdb2.prj', is_binary=False, optional=True)
        self.add_composite_file('lastdb3.ssp', is_binary=True, optional=True)
        self.add_composite_file('lastdb3.tis', is_binary=True, optional=True)
        self.add_composite_file('lastdb3.sds', is_binary=True, optional=True)
        self.add_composite_file('lastdb3.suf', is_binary=True, optional=True)
        self.add_composite_file('lastdb3.des', is_binary=True, optional=True)
        self.add_composite_file('lastdb3.bck', is_binary=True, optional=True)
        self.add_composite_file('lastdb3.prj', is_binary=False, optional=True)

class LastalProtDb( _LastalDb, Data ):
    """Class for protein LAST database files."""
    file_ext = 'lastdbp'
    allow_datatype_change = False
    composite_type = 'basic'

    def __init__(self, **kwd):
        Data.__init__(self, **kwd)
        self.add_composite_file('lastdb.ssp', is_binary=True)
        self.add_composite_file('lastdb.tis', is_binary=True)
        self.add_composite_file('lastdb.sds', is_binary=True)
        self.add_composite_file('lastdb.suf', is_binary=True)
        self.add_composite_file('lastdb.des', is_binary=True)
        self.add_composite_file('lastdb.bck', is_binary=True)
        self.add_composite_file('lastdb.prj', is_binary=False)
        self.add_composite_file('lastdb0.ssp', is_binary=True, optional=True)
        self.add_composite_file('lastdb0.tis', is_binary=True, optional=True)
        self.add_composite_file('lastdb0.sds', is_binary=True, optional=True)
        self.add_composite_file('lastdb0.suf', is_binary=True, optional=True)
        self.add_composite_file('lastdb0.des', is_binary=True, optional=True)
        self.add_composite_file('lastdb0.bck', is_binary=True, optional=True)
        self.add_composite_file('lastdb0.prj', is_binary=False, optional=True)
        self.add_composite_file('lastdb1.ssp', is_binary=True, optional=True)
        self.add_composite_file('lastdb1.tis', is_binary=True, optional=True)
        self.add_composite_file('lastdb1.sds', is_binary=True, optional=True)
        self.add_composite_file('lastdb1.suf', is_binary=True, optional=True)
        self.add_composite_file('lastdb1.des', is_binary=True, optional=True)
        self.add_composite_file('lastdb1.bck', is_binary=True, optional=True)
        self.add_composite_file('lastdb1.prj', is_binary=False, optional=True)
        self.add_composite_file('lastdb2.ssp', is_binary=True, optional=True)
        self.add_composite_file('lastdb2.tis', is_binary=True, optional=True)
        self.add_composite_file('lastdb2.sds', is_binary=True, optional=True)
        self.add_composite_file('lastdb2.suf', is_binary=True, optional=True)
        self.add_composite_file('lastdb2.des', is_binary=True, optional=True)
        self.add_composite_file('lastdb2.bck', is_binary=True, optional=True)
        self.add_composite_file('lastdb2.prj', is_binary=False, optional=True)
        self.add_composite_file('lastdb3.ssp', is_binary=True, optional=True)
        self.add_composite_file('lastdb3.tis', is_binary=True, optional=True)
        self.add_composite_file('lastdb3.sds', is_binary=True, optional=True)
        self.add_composite_file('lastdb3.suf', is_binary=True, optional=True)
        self.add_composite_file('lastdb3.des', is_binary=True, optional=True)
        self.add_composite_file('lastdb3.bck', is_binary=True, optional=True)
        self.add_composite_file('lastdb3.prj', is_binary=False, optional=True)

