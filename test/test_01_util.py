import edl.util

def test_reservoir_sampling_simple():
    ten_numbers, one_hundred = edl.util.reservoirSample(xrange(100), 
                                                        N=10, 
                                                        returnCount=True)
    assert len(ten_numbers) == 10
    assert one_hundred == 100

    five_numbers = edl.util.reservoirSample(xrange(1000),
                                            N=5)

    assert len(five_numbers) == 5

    short_test_list = [1,2,3]
    sampled_list, list_len = edl.util.reservoirSample(short_test_list, 
                                                      N=10,
                                                      returnCount=True)
    assert short_test_list == sampled_list
    assert list_len == len(short_test_list)
    
    sampled_list = edl.util.reservoirSample(short_test_list, N=10)
    assert short_test_list == sampled_list
    
