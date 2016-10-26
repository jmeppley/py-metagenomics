import edl.util

def test_reservoir_sampling_simple():
    ten_numbers, one_hundred = edl.util.reservoir_sample(range(100), 
                                                        N=10, 
                                                        return_count=True)
    assert len(ten_numbers) == 10
    assert len(set(ten_numbers)) == 10
    assert one_hundred == 100

    five_numbers = edl.util.reservoir_sample(range(1000),
                                            N=5)

    assert len(five_numbers) == 5

    short_test_list = [1,2,3]
    sampled_list, list_len = edl.util.reservoir_sample(short_test_list, 
                                                      N=10,
                                                      return_count=True)
    assert short_test_list == sampled_list
    assert list_len == len(short_test_list)
    
    sampled_list = edl.util.reservoir_sample(short_test_list, N=10)
    assert short_test_list == sampled_list
    

def test_parse_map_file():
    test_file='test/data/acc.to.ko.protein.plus.filtered'
    file_map=edl.util.parseMapFile(test_file)
    assert len(file_map) == 53
    assert file_map['YP_396597'] == 'K01092'

    test_file='test/data/acc.to.taxid.proetin.plus.filtered'
    file_map=edl.util.parseMapFile(test_file, valueType=int)
    assert len(file_map) == 75
    assert file_map['YP_397454'] == 74546
    
    test_file='test/data/HOT186_25m_454DNA.vs.RefSeqMicrobesProteins.250records.blastx.b50.m8.NCyanos.topHitByAcc'
    file_map=edl.util.parseMapFile(test_file, valueCol=2, keyCol=0)
    assert len(file_map) == 42
    assert file_map['000930_1547_0856'] == 'ref|YP_396597.1|'


def test_parse_file_to_set():
    test_file='test/data/acc.to.taxid.proetin.plus.filtered'
    value_set = edl.util.parse_list_to_set(test_file, keyType=int, col=1)
    assert len(value_set) == 7
    assert sorted(list(value_set)) == [59919, 74546, 93058, 93060, 146891, 167542, 167546]

