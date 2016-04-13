import edl.kegg, logging, sys

def test_kegg():
    loglevel=logging.WARN
    logging.basicConfig(stream=sys.stderr, level=loglevel)
    edl.kegg.logger.setLevel(loglevel)

    edl.kegg.kegg_nosetest('test/data/ko.map.partial',
            'test/data/kobrite/ko00001.keg')
