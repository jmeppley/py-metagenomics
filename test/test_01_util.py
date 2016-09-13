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
    
